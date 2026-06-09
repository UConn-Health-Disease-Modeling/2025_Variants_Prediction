#' Select ~K informative features for regression on a continuous target
#'
#' Runs a robust feature-selection pipeline for a continuous target (e.g., \code{maximum}):
#' drop all-NA/zero-var/NZV, impute missing, high-correlation pruning, LASSO (Gaussian),
#' then top-up to \code{target_k} using |corr(x, y)|, with variance as fallback.
#' Keeps identifier columns (e.g., country, classified_label, cluster) in the output if present.
#'
#' @param model_data A data.frame/tibble containing the target column, optional ID columns
#'   (e.g., \code{country}, \code{classified_label}, \code{cluster}), and numeric features.
#' @param target_col Character. Name of the continuous target column (default \code{"maximum"}).
#' @param id_cols Character vector of identifier columns to keep in the output
#'   (default \code{c("country","classified_label","cluster")}).
#' @param target_k Integer. Target number of features to keep (default 10).
#' @param cor_cutoff Numeric in (0,1). Threshold for high-correlation pruning (default 0.95).
#' @param lasso_s Character. Which lambda from cv.glmnet: \code{"lambda.min"} or \code{"lambda.1se"} (default \code{"lambda.min"}).
#' @param nfolds Integer. Number of folds for \code{cv.glmnet} (default 5).
#' @param impute_method Character. Missing imputation for features: \code{"zero"} (default), \code{"mean"}, \code{"median"}.
#' @param seed Integer. Random seed (default 123).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{final_data}: tibble with \emph{id columns} + selected features + target column.
#'   \item \code{selected_features}: character vector of kept feature names.
#'   \item \code{cvfit}: \code{glmnet} CV object (Gaussian).
#' }
#' @export
select_features_for_regression <- function(model_data,
                                           target_col   = "maximum",
                                           id_cols      = c("country","classified_label","cluster"),
                                           target_k     = 10,
                                           cor_cutoff   = 0.95,
                                           lasso_s      = "lambda.min",
                                           nfolds       = 5,
                                           impute_method = "zero",
                                           seed         = 123) {
  # --- checks ---
  if (!(target_col %in% names(model_data))) {
    stop("`model_data` must contain the continuous target column: ", target_col)
  }
  
  # keep id columns if present (in original order)
  ids_df <- model_data |>
    dplyr::select(dplyr::any_of(id_cols))
  
  # --- prepare: drop obvious non-features from feature matrix, keep numeric only ---
  df <- model_data |>
    dplyr::select(-dplyr::any_of(setdiff(id_cols, character(0))))  # remove id cols from feature pool
  
  y <- df[[target_col]]
  if (!is.numeric(y)) {
    stop("Target column `", target_col, "` must be numeric (continuous).")
  }
  
  X <- df |>
    dplyr::select(-dplyr::all_of(target_col)) |>
    dplyr::select(where(is.numeric)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ {x <- .; x[!is.finite(x)] <- NA_real_; x}))
  
  # --- drop all-NA / zero-var / NZV ---
  non_all_na   <- sapply(X, function(col) !all(is.na(col)))
  X <- X[, non_all_na, drop = FALSE]
  non_zero_var <- sapply(X, function(col) stats::var(col, na.rm = TRUE) > 0)
  X <- X[, non_zero_var, drop = FALSE]
  nzv_idx <- caret::nearZeroVar(X)
  if (length(nzv_idx) > 0) X <- X[, -nzv_idx, drop = FALSE]
  if (ncol(X) == 0) stop("No usable numeric features remain after variance/NA filtering.")
  
  # --- impute missing ---
  X <- dplyr::as_tibble(X)
  if (impute_method == "zero") {
    X <- X |> dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.), 0, .)))
  } else if (impute_method == "mean") {
    X <- X |> dplyr::mutate(dplyr::across(dplyr::everything(),
                                          ~ {m <- mean(., na.rm = TRUE); ifelse(is.na(.), m, .)}))
  } else if (impute_method == "median") {
    X <- X |> dplyr::mutate(dplyr::across(dplyr::everything(),
                                          ~ {m <- stats::median(., na.rm = TRUE); ifelse(is.na(.), m, .)}))
  } else {
    stop("Unsupported `impute_method`. Use 'zero', 'mean', or 'median'.")
  }
  
  # --- high-correlation pruning ---
  cor_mat <- stats::cor(X, use = "pairwise.complete.obs")
  cor_mat[!is.finite(cor_mat)] <- 0; cor_mat[is.na(cor_mat)] <- 0
  hi <- caret::findCorrelation(cor_mat, cutoff = cor_cutoff)
  if (length(hi) > 0) X <- X[, -hi, drop = FALSE]
  if (ncol(X) == 0) stop("All features removed by correlation pruning; consider raising `cor_cutoff`.")
  
  # --- LASSO (Gaussian) ---
  set.seed(seed)
  X_mat <- as.matrix(X)
  cvfit <- glmnet::cv.glmnet(
    X_mat, y,
    family = "gaussian",
    alpha  = 1,
    standardize = TRUE,
    nfolds = nfolds
  )
  if (!lasso_s %in% c("lambda.min","lambda.1se")) {
    stop("`lasso_s` must be 'lambda.min' or 'lambda.1se'.")
  }
  coef_lasso <- stats::coef(cvfit, s = lasso_s)
  beta_all   <- as.numeric(coef_lasso); names(beta_all) <- rownames(coef_lasso)
  sel <- setdiff(names(beta_all)[beta_all != 0], "(Intercept)")
  
  # --- ensure ~K features: top-up with |corr(x, y)| then variance ---
  max_possible <- min(target_k, ncol(X))
  if (length(sel) >= max_possible) {
    beta_sel <- abs(beta_all[sel])
    sel_final <- names(sort(beta_sel, decreasing = TRUE))[1:max_possible]
  } else {
    remaining <- setdiff(colnames(X), sel)
    cor_score <- sapply(remaining, function(col) {
      suppressWarnings(stats::cor(X[[col]], y, use = "pairwise.complete.obs"))
    })
    cor_score[!is.finite(cor_score)] <- 0
    need_more <- max_possible - length(sel)
    add_feats <- names(sort(abs(cor_score), decreasing = TRUE))[seq_len(min(need_more, length(cor_score)))]
    sel_final <- unique(c(sel, add_feats))
    
    if (length(sel_final) < max_possible) {
      leftover <- setdiff(colnames(X), sel_final)
      if (length(leftover) > 0) {
        var_rank <- sort(sapply(X[, leftover, drop = FALSE], stats::var), decreasing = TRUE)
        fill <- names(var_rank)[seq_len(min(max_possible - length(sel_final), length(var_rank)))]
        sel_final <- c(sel_final, fill)
      }
    }
  }
  if (length(sel_final) > max_possible) sel_final <- sel_final[1:max_possible]
  
  final_features <- X[, sel_final, drop = FALSE]
  
  # --- bind IDs + features + target (keep original row order) ---
  final_data <- dplyr::bind_cols(
    ids_df,                           # keep id columns if they exist
    final_features,
    tibble::tibble(!!target_col := y) # target at the end
  )
  
  message(paste0("Selected features (", ncol(final_features), "): ",
                 paste(colnames(final_features), collapse = ", ")))
  
  list(
    final_data = tibble::as_tibble(final_data),
    selected_features = colnames(final_features),
    cvfit = cvfit
  )
}



#' Fit early logistic curves per variant and plot TTD panels by window size
#'
#' For a given \code{country}, this function fits a simple logistic curve
#' \deqn{share(day) = 1 / (1 + exp((xmid - day)/scal))} to each
#' \code{classified_label} using observations within early windows \code{ns}
#' (e.g., 14/21/28 days since first detection). It extracts the predicted
#' time-to-50\% (TTD; \code{xmid}) and draws a panel plot: rows are variants,
#' columns are different \code{n} values in \code{ns}. The red vertical line
#' marks the first-detection day from \code{refer_df}; the blue dashed line marks
#' the predicted TTD from the logistic fit.
#'
#' @param main_df A data.frame/tibble containing at least the columns:
#'   \code{country}, \code{classified_label}, \code{day} (integer or numeric),
#'   and \code{share} (0–1).
#' @param country_name Character scalar. Country to plot/fit.
#' @param ns Integer vector of early windows in days (e.g., \code{c(14, 21, 28)}).
#' @param refer_df A data.frame/tibble with columns \code{country},
#'   \code{classified_label}, and \code{first_day} (integer/numeric), used to
#'   draw the reference vertical line.
#' @param mar Numeric vector passed to \code{graphics::par(mar = ...)} for panel
#'   margins. Default \code{c(4, 4, 2, 1)}.
#'
#' @return (Invisibly) a list with:
#' \itemize{
#'   \item \code{results}: tibble \code{(classified_label, t_D_pred, n)}
#'   \item \code{curves}:  tibble \code{(classified_label, day, pred_share, n)}
#' }
#' Side effect: draws a multi-panel base R plot (rows = variants, cols = |ns|).
#'
#' @details
#' Each variant is fit only if it has at least 5 observations in the window
#' \code{0..n}. Fits are done via \code{minpack.lm::nlsLM} with starting values
#' \code{xmid = median(day)}, \code{scal = 1} and bounds \code{xmid >= 0},
#' \code{scal >= 0.01}. Predicted curves are truncated at the first time
#' \code{share >= 0.5} for clarity.
#'
#' @importFrom dplyr filter arrange group_by group_split bind_rows mutate pull
#' @importFrom tibble tibble
#' @importFrom graphics par plot lines abline mtext
#' @export
fit_logistic_TTD_and_plot <- function(
    main_df,
    country_name,
    ns = 30,
    refer_df,
    mar = c(4, 4, 2, 1)
) {
  # ---- checks ----
  req_cols <- c("country","classified_label","day","share")
  miss_main <- setdiff(req_cols, names(main_df))
  if (length(miss_main)) stop("`main_df` missing columns: ", paste(miss_main, collapse = ", "))
  if (!all(c("country","classified_label","first_day") %in% names(refer_df))) {
    stop("`refer_df` must contain columns: country, classified_label, first_day")
  }
  if (!is.character(country_name) || length(country_name) != 1) {
    stop("`country_name` must be a single character string")
  }
  if (!is.numeric(ns) || length(ns) < 1) stop("`ns` must be a numeric vector of window sizes")
  
  ns <- as.integer(ns)
  
  # ---- helpers: per-n fit ----
  fit_logistic_TTD_and_plot_single <- function(main_df, country_name, n, refer_df) {
    df <- main_df %>%
      dplyr::filter(.data$country == country_name, .data$day >= 0, .data$day <= n) %>%
      dplyr::arrange(.data$day)
    
    res_list <- list()
    cur_list <- list()
    
    grp_df <- df %>%
      dplyr::group_by(.data$classified_label) %>%
      dplyr::group_split()
    
    for (grp in grp_df) {
      lbl <- grp$classified_label[1]
      tD  <- NA_real_
      
      if (nrow(grp) >= 5) {
        fit <- try(
          minpack.lm::nlsLM(
            share ~ 1/(1 + exp((xmid - day)/scal)),
            data    = grp,
            start   = list(xmid = stats::median(grp$day), scal = 1),
            lower   = c(xmid = 0, scal = 0.01),
            upper   = c(xmid = Inf, scal = Inf),
            control = minpack.lm::nls.lm.control(maxiter = 100)
          ),
          silent = TRUE
        )
        
        if (!inherits(fit, "try-error")) {
          cf <- stats::coef(fit)
          tD <- as.numeric(cf["xmid"])
          
          days_f <- seq(0, ceiling(tD) + 1, by = 1)
          preds  <- stats::predict(fit, newdata = data.frame(day = days_f))
          idx50  <- which(preds >= 0.5)
          last_i <- if (length(idx50)) idx50[1] else length(days_f)
          
          cur_list[[lbl]] <- tibble::tibble(
            classified_label = lbl,
            day              = days_f[1:last_i],
            pred_share       = preds[1:last_i]
          )
        }
      }
      
      res_list[[lbl]] <- tibble::tibble(
        classified_label = lbl,
        t_D_pred         = tD
      )
    }
    
    list(
      results = dplyr::bind_rows(res_list),
      curves  = dplyr::bind_rows(cur_list)
    )
  }
  
  # ---- fit for each n in ns ----
  all_results <- list()
  all_curves  <- list()
  
  for (n_val in ns) {
    out <- fit_logistic_TTD_and_plot_single(main_df, country_name, n_val, refer_df)
    all_results[[as.character(n_val)]] <- out$results %>% dplyr::mutate(n = n_val)
    all_curves[[ as.character(n_val)]] <- out$curves  %>% dplyr::mutate(n = n_val)
  }
  
  results_df <- dplyr::bind_rows(all_results)
  curves_df  <- dplyr::bind_rows(all_curves)
  
  # if nothing fit, exit early (but still reset par)
  if (nrow(curves_df) == 0L) {
    warning("No curves were fit (insufficient points per variant/window?).")
    return(invisible(list(results = results_df, curves = curves_df)))
  }
  
  variants <- unique(curves_df$classified_label)
  
  # ---- plot grid (rows = variants, cols = |ns|) ----
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::par(mfcol = c(length(variants), length(ns)), mar = mar)
  
  for (n_val in ns) {
    for (lbl in variants) {
      obs <- main_df %>%
        dplyr::filter(.data$country == country_name, .data$classified_label == lbl) %>%
        dplyr::arrange(.data$day)
      
      fit <- curves_df %>%
        dplyr::filter(.data$n == n_val, .data$classified_label == lbl) %>%
        dplyr::arrange(.data$day)
      
      tD <- results_df %>%
        dplyr::filter(.data$n == n_val, .data$classified_label == lbl) %>%
        dplyr::pull(.data$t_D_pred)
      
      first_day <- refer_df %>%
        dplyr::filter(.data$country == country_name, .data$classified_label == lbl) %>%
        dplyr::pull(.data$first_day)
      
      y_min <- suppressWarnings(min(c(obs$share, fit$pred_share), na.rm = TRUE))
      y_max <- suppressWarnings(max(c(obs$share, fit$pred_share), na.rm = TRUE))
      if (!is.finite(y_min) || !is.finite(y_max)) { y_min <- 0; y_max <- 1 }
      
      graphics::plot(
        x = obs$day, y = obs$share, type = "l", lwd = 1, col = "black",
        ylim = c(y_min, y_max),
        xlab = if (lbl == variants[length(variants)]) "Day" else "",
        ylab = if (n_val == ns[1]) "Share" else "",
        main = if (lbl == variants[1]) paste0("n = ", n_val) else ""
      )
      
      if (nrow(fit)) graphics::lines(fit$day, fit$pred_share, lwd = 2, col = "blue")
      if (length(tD) == 1L && is.finite(tD)) graphics::abline(v = tD, col = "blue", lty = 2)
      if (length(first_day) == 1L && is.finite(first_day)) graphics::abline(v = first_day, col = "red", lty = 1)
      
      if (n_val == ns[1]) {
        mid <- if (is.finite(y_min + y_max)) (y_min + y_max) / 2 else 0.5
        graphics::mtext(lbl, side = 2, line = 2, at = mid, cex = 0.8)
      }
    }
  }
  
  invisible(list(results = results_df, curves = curves_df))
}


#' Fit a per-variant early logistic curve and extract TTD (time-to-50%)
#'
#' For a given \code{country}, fit a 2-parameter logistic curve
#' \deqn{share(day) = 1 / (1 + exp((xmid - day)/scal))} to each
#' \code{classified_label} over the window \code{day ∈ [0, n]}. Returns the
#' predicted time-to-50\% (\code{xmid}, denoted \code{t_D_pred}) and the
#' fitted curve up to the first crossing of 50\%.
#'
#' @param main_df A data.frame/tibble with columns:
#'   \code{country}, \code{classified_label}, \code{day} (numeric/integer),
#'   and \code{share} (in [0,1]).
#' @param country_name Character scalar; the country to fit.
#' @param n Integer; max day (inclusive) of the early window. Default \code{30}.
#'
#' @return A list with two tibbles:
#' \itemize{
#'   \item \code{results}: columns \code{country}, \code{classified_label}, \code{t_D_pred}
#'   \item \code{curves}:  columns \code{country}, \code{classified_label}, \code{day}, \code{pred_share}
#' }
#'
#' @details
#' A variant is fit only if it has at least 5 observations within \code{0..n}.
#' Fitting uses \code{minpack.lm::nlsLM} with starts \code{xmid = median(day)},
#' \code{scal = 1}, bounds \code{xmid >= 0}, \code{scal >= 0.01}. The returned
#' curve is truncated at the first day where predicted \code{share >= 0.5}.
#'
#' @importFrom dplyr filter arrange group_by group_split bind_rows
#' @importFrom tibble tibble
#' @export
fit_logistic_TTD_and_curve <- function(main_df, country_name, n = 30) {
  # ---- checks ----
  req <- c("country","classified_label","day","share")
  miss <- setdiff(req, names(main_df))
  if (length(miss)) stop("`main_df` missing columns: ", paste(miss, collapse = ", "))
  if (!is.character(country_name) || length(country_name) != 1)
    stop("`country_name` must be a single character string.")
  if (!is.numeric(n) || length(n) != 1 || n < 0)
    stop("`n` must be a single non-negative numeric value.")
  n <- as.integer(n)
  
  # ---- filter to country & window ----
  df <- main_df %>%
    dplyr::filter(.data$country == country_name,
                  .data$day >= 0, .data$day <= n) %>%
    dplyr::arrange(.data$day)
  
  results_list <- list()
  curves_list  <- list()
  
  # ---- split by variant ----
  grouped <- df %>%
    dplyr::group_by(.data$classified_label) %>%
    dplyr::group_split()
  
  for (grp in grouped) {
    lbl <- grp$classified_label[1]
    tD  <- NA_real_
    
    if (nrow(grp) >= 5) {
      fit <- try(
        minpack.lm::nlsLM(
          share ~ 1 / (1 + exp((xmid - day) / scal)),
          data    = grp,
          start   = list(xmid = stats::median(grp$day), scal = 1),
          lower   = c(xmid = 0, scal = 0.01),
          upper   = c(xmid = Inf, scal = Inf),
          control = minpack.lm::nls.lm.control(maxiter = 100)
        ),
        silent = TRUE
      )
      
      if (!inherits(fit, "try-error")) {
        cf  <- stats::coef(fit)
        tD  <- as.numeric(cf["xmid"])
        
        max_day    <- ceiling(tD) + 1
        days_full  <- seq(0, max_day, by = 1)
        preds_full <- stats::predict(fit, newdata = data.frame(day = days_full))
        
        idx50    <- which(preds_full >= 0.5)
        last_idx <- if (length(idx50)) idx50[1] else length(days_full)
        
        curves_list[[lbl]] <- tibble::tibble(
          country          = country_name,
          classified_label = lbl,
          day              = days_full[1:last_idx],
          pred_share       = preds_full[1:last_idx]
        )
      }
    }
    
    results_list[[lbl]] <- tibble::tibble(
      country          = country_name,
      classified_label = lbl,
      t_D_pred         = tD
    )
  }
  
  results <- dplyr::bind_rows(results_list)
  curves  <- dplyr::bind_rows(curves_list)
  
  list(results = results, curves = curves)
}

#' Predict TTD (time-to-50%) from an all-<50% share sequence
#'
#' Given a time-ordered share sequence \code{y_t} where **all values are < 0.5**,
#' fit a 2-parameter logistic curve
#' \deqn{share(t) = 1 / (1 + exp((xmid - t)/scal))} via \code{minpack.lm::nlsLM}
#' and return the **predicted day** when the curve reaches 50% (i.e., \code{ceiling(xmid)}).
#'
#' If the input violates the precondition (any share >= 0.5), the function will
#' return the first observed crossing day instead of fitting.
#'
#' @param share_seq Numeric vector of shares in [0,1), strictly < 0.5 by design.
#' @param time_index Optional numeric vector of the same length giving the day index
#'   for each observation. Defaults to a 0-based integer index \code{0:(n-1)}.
#' @param min_points Integer, minimum points required to attempt a fit when there is
#'   no observed crossing (default 5).
#'
#' @return Integer day of first predicted crossing to 50%:
#' \itemize{
#'   \item If any observed share >= 0.5: the first corresponding day (from \code{time_index}).
#'   \item Else if logistic fit succeeds: \code{ceiling(xmid)}.
#'   \item Else: \code{NA_integer_}.
#' }
#' @examples
#' set.seed(1)
#' y <- plogis((0:12 - 10)/2) * 0.49          # all < 0.5
#' ttd <- ttd_from_sub50_sequence(y)          # integer day when predicted crosses 0.5
#' @export
ttd_from_sub50_sequence <- function(share_seq,
                                    time_index = NULL,
                                    min_points = 5) {
  # --- checks & prep ---
  if (!is.numeric(share_seq)) stop("`share_seq` must be numeric.")
  n <- length(share_seq)
  if (n == 0) return(NA_integer_)
  if (is.null(time_index)) {
    time_index <- seq.int(0L, n - 1L)           # 0-based by default
  } else {
    if (!is.numeric(time_index) || length(time_index) != n) {
      stop("`time_index` must be numeric and the same length as `share_seq`.")
    }
  }
  
  # keep finite pairs only
  ok <- is.finite(share_seq) & is.finite(time_index)
  y  <- as.numeric(share_seq[ok])
  t  <- as.numeric(time_index[ok])
  if (length(y) == 0) return(NA_integer_)
  
  # If precondition violated (some >= 0.5), return first observed crossing day
  idx_obs <- which(y >= 0.5)
  if (length(idx_obs) > 0) {
    return(as.integer(ceiling(t[min(idx_obs)])))
  }
  
  # Need enough variability and points to fit
  if (length(y) < min_points || stats::var(y, na.rm = TRUE) == 0) {
    return(NA_integer_)
  }
  
  # --- logistic fit: share ~ 1/(1 + exp((xmid - t)/scal)) ---
  df <- data.frame(day = t, share = y)
  start_xmid <- stats::median(t, na.rm = TRUE)
  start_scal <- 1
  
  fit <- try(
    minpack.lm::nlsLM(
      share ~ 1 / (1 + exp((xmid - day) / scal)),
      data    = df,
      start   = list(xmid = start_xmid, scal = start_scal),
      lower   = c(xmid = min(t, na.rm = TRUE), scal = 0.01),
      upper   = c(xmid = max(t, na.rm = TRUE) + 1e6, scal = Inf),
      control = minpack.lm::nls.lm.control(maxiter = 100)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) return(NA_integer_)
  
  xmid <- as.numeric(stats::coef(fit)["xmid"])
  if (!is.finite(xmid)) return(NA_integer_)
  
  as.integer(ceiling(xmid))
}





# ------------------------------------------------------------------------------
# Fit several regressors with simple anti-overfitting safeguards
# - train-only imputation
# - NZV removal
# - train-only correlation filter
# - conservative hyperparameter grids
# Returns train & test metrics, fitted models, and preprocessing artifacts.
# ------------------------------------------------------------------------------
fit_eval_models <- function(X, y, seed = 123, corr_cutoff = 0.90) {
  stopifnot(is.numeric(y), length(y) == nrow(X))

  # --- numeric-only design matrix --------------------------------------------
  Xnum <- X %>% dplyr::select(where(is.numeric))
  Xnum[] <- lapply(Xnum, function(col) { col[!is.finite(col)] <- NA_real_; col })
  df_all <- data.frame(y = y, Xnum, check.names = FALSE)

  # --- split 80/20 ------------------------------------------------------------
  set.seed(seed)
  idx_tr <- caret::createDataPartition(df_all$y, p = 0.80, list = FALSE)
  tr <- df_all[idx_tr, , drop = FALSE]
  te <- df_all[-idx_tr, , drop = FALSE]
  if (stats::sd(tr$y, na.rm = TRUE) == 0) stop("Training response has zero variance.")

  # --- train-only imputation --------------------------------------------------
  pp_imp <- caret::preProcess(tr, method = "medianImpute")
  tr_imp <- stats::predict(pp_imp, tr)
  te_imp <- stats::predict(pp_imp, te)

  # --- NZV removal (train-driven) --------------------------------------------
  nzv <- caret::nearZeroVar(tr_imp, saveMetrics = TRUE)
  if (any(nzv$nzv)) {
    keep_nzv <- !nzv$nzv
    tr_imp <- tr_imp[, keep_nzv, drop = FALSE]
    te_imp <- te_imp[,  keep_nzv, drop = FALSE]
  }

  # --- correlation filter (train-only; y excluded) ---------------------------
  drop_corr <- character(0)
  if (ncol(tr_imp) > 2) {
    cor_mat <- stats::cor(tr_imp[, -1, drop = FALSE], use = "pairwise.complete.obs")
    drop_corr <- tryCatch(
      caret::findCorrelation(cor_mat, cutoff = corr_cutoff, names = TRUE, exact = TRUE),
      error = function(e) character(0)
    )
    if (length(drop_corr)) {
      keep_corr <- setdiff(colnames(tr_imp), drop_corr)
      keep_corr <- union("y", keep_corr) # always keep response
      tr_imp <- tr_imp[, keep_corr, drop = FALSE]
      te_imp <- te_imp[,  keep_corr, drop = FALSE]
    }
  }

  # --- common CV --------------------------------------------------------------
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verboseIter = FALSE,
    allowParallel = TRUE
  )

  results <- list()
  metrics_train <- list()
  metrics_test  <- list()

  add_metrics <- function(name, pred_tr, y_tr, pred_te, y_te) {
    mt_tr <- caret::postResample(pred_tr, y_tr)
    mt_te <- caret::postResample(pred_te, y_te)
    list(
      train = tibble::tibble(
        model = name,
        RMSE  = unname(mt_tr["RMSE"]),
        MSE   = (unname(mt_tr["RMSE"]))^2,
        R2    = unname(mt_tr["Rsquared"])
      ),
      test  = tibble::tibble(
        model = name,
        RMSE  = unname(mt_te["RMSE"]),
        MSE   = (unname(mt_te["RMSE"]))^2,
        R2    = unname(mt_te["Rsquared"])
      )
    )
  }

  # =================== MODELS ===================

  ## -------- Elastic Net (scaled) --------------------------------------------
  pp_cs <- caret::preProcess(tr_imp[, -1, drop = FALSE], method = c("center", "scale"))
  Xtr_glm <- stats::predict(pp_cs, tr_imp[, -1, drop = FALSE])
  Xte_glm <- stats::predict(pp_cs, te_imp[, -1, drop = FALSE])
  df_tr_glm <- data.frame(y = tr_imp$y, Xtr_glm)
  df_te_glm <- data.frame(y = te_imp$y, Xte_glm)

  glmn_grid <- expand.grid(
    alpha  = seq(0, 1, by = 0.25),
    lambda = 10^seq(-4, 0.5, length.out = 40)
  )
  set.seed(seed)
  fit_glmn <- caret::train(
    y ~ ., data = df_tr_glm, method = "glmnet",
    trControl = ctrl, tuneGrid = glmn_grid, metric = "RMSE"
  )
  pred_tr <- stats::predict(fit_glmn, newdata = df_tr_glm)
  pred_te <- stats::predict(fit_glmn, newdata = df_te_glm)
  mm <- add_metrics("ElasticNet", pred_tr, df_tr_glm$y, pred_te, df_te_glm$y)
  results$elasticnet <- fit_glmn
  metrics_train$elasticnet <- mm$train
  metrics_test$elasticnet  <- mm$test

  ## -------- Random Forest (ranger) ------------------------------------------
  p <- ncol(tr_imp) - 1
  rf_grid <- expand.grid(
    mtry = unique(pmax(1, round(c(sqrt(p), p/3)))),
    splitrule = "variance",
    min.node.size = c(5, 10)
  )
  set.seed(seed)
  fit_rf <- caret::train(
    y ~ ., data = tr_imp, method = "ranger",
    trControl = ctrl, tuneGrid = rf_grid,
    num.trees = 500,
    sample.fraction = 0.60,
    importance = "impurity"
  )
  pred_tr <- stats::predict(fit_rf, tr_imp)
  pred_te <- stats::predict(fit_rf, te_imp)
  mm <- add_metrics("RandomForest", pred_tr, tr_imp$y, pred_te, te_imp$y)
  results$random_forest <- fit_rf
  metrics_train$random_forest <- mm$train
  metrics_test$random_forest  <- mm$test

  ## -------- XGBoost (conservative) ------------------------------------------
  xgb_grid <- expand.grid(
    nrounds = c(150, 250, 350),
    max_depth = c(2, 3),
    eta = c(0.05, 0.10),
    gamma = c(0, 1),
    colsample_bytree = 0.60,
    min_child_weight = c(3, 5),
    subsample = 0.60
  )
  set.seed(seed)
  fit_xgb <- caret::train(
    y ~ ., data = tr_imp, method = "xgbTree",
    trControl = ctrl, tuneGrid = xgb_grid, metric = "RMSE",
    verbose = FALSE
  )
  # suppress ntree_limit warnings from xgboost's predict wrapper
  pred_tr <- suppressWarnings(stats::predict(fit_xgb, tr_imp))
  pred_te <- suppressWarnings(stats::predict(fit_xgb, te_imp))
  mm <- add_metrics("XGBoost", pred_tr, tr_imp$y, pred_te, te_imp$y)
  results$xgboost <- fit_xgb
  metrics_train$xgboost <- mm$train
  metrics_test$xgboost  <- mm$test

  ## -------- GBM (feasible grid for smallest fold) ---------------------------
  # Smallest training size across 5-fold CV
  min_fold_train <- floor(nrow(tr_imp) * (1 - 1/5))

  # Candidate grids (kept modest to reduce variance)
  depth_candidates   <- c(1, 2)
  trees_candidates   <- c(600, 900, 1200)
  shrink_candidates  <- c(0.03, 0.05)
  nmin_candidates    <- c(3, 5, 10)

  gbm_grid <- expand.grid(
    n.trees = trees_candidates,
    interaction.depth = depth_candidates,
    shrinkage = shrink_candidates,
    n.minobsinnode = nmin_candidates
  )

  # Choose a single bag.fraction that is feasible for the largest n.minobsinnode
  bf_candidates <- c(0.95, 0.90, 0.80, 0.70, 0.60)
  max_nmin <- max(nmin_candidates)
  feasible_bf <- bf_candidates[
    (min_fold_train * bf_candidates) > (2 * max_nmin + 1)
  ]
  safe_bf <- if (length(feasible_bf)) max(feasible_bf) else {
    # fallback: compute the largest feasible bf within (0.5, 0.95]
    bf_calc <- ( (2 * max_nmin + 2) / max(1, min_fold_train) ) + 1e-3
    bf <- max(0.50, min(0.95, 1 - bf_calc))  # conservative clamp
    if (!is.finite(bf) || bf <= 0) 0.60 else bf
  }

  set.seed(seed)
  fit_gbm <- caret::train(
    y ~ ., data = tr_imp, method = "gbm",
    trControl = ctrl, tuneGrid = gbm_grid,
    metric = "RMSE", verbose = FALSE,
    distribution = "gaussian",
    bag.fraction = safe_bf
  )
  pred_tr <- stats::predict(fit_gbm, tr_imp)
  pred_te <- stats::predict(fit_gbm, te_imp)
  mm <- add_metrics("GBM", pred_tr, tr_imp$y, pred_te, te_imp$y)
  results$gbm <- fit_gbm
  metrics_train$gbm <- mm$train
  metrics_test$gbm  <- mm$test

  ## -------- SVM Radial (scaled) ---------------------------------------------
  pp_cs_svm <- caret::preProcess(tr_imp[, -1, drop = FALSE], method = c("center","scale"))
  Xtr_svm <- stats::predict(pp_cs_svm, tr_imp[, -1, drop = FALSE])
  Xte_svm <- stats::predict(pp_cs_svm, te_imp[, -1, drop = FALSE])
  df_tr_svm <- data.frame(y = tr_imp$y, Xtr_svm)
  df_te_svm <- data.frame(y = te_imp$y, Xte_svm)

  sigma_guess <- tryCatch(kernlab::sigest(as.matrix(Xtr_svm))[[1]], error = function(e) 0.1)
  svm_grid <- expand.grid(
    sigma = sigma_guess,
    C = 2^seq(-2, 3, by = 1)
  )
  set.seed(seed)
  fit_svm <- caret::train(
    y ~ ., data = df_tr_svm, method = "svmRadial",
    trControl = ctrl, tuneGrid = svm_grid, metric = "RMSE"
  )
  pred_tr <- stats::predict(fit_svm, df_tr_svm)
  pred_te <- stats::predict(fit_svm, df_te_svm)
  mm <- add_metrics("SVM_Radial", pred_tr, df_tr_svm$y, pred_te, df_te_svm$y)
  results$svm_radial <- fit_svm
  metrics_train$svm_radial <- mm$train
  metrics_test$svm_radial  <- mm$test

  ## -------- MARS (earth) -----------------------------------------------------
  mars_grid <- expand.grid(
    nprune = seq(5, min(40, ncol(tr_imp) - 1), by = 5),
    degree = 1
  )
  set.seed(seed)
  fit_mars <- caret::train(
    y ~ ., data = tr_imp, method = "earth",
    trControl = ctrl, tuneGrid = mars_grid, metric = "RMSE"
  )
  pred_tr <- stats::predict(fit_mars, tr_imp)
  pred_te <- stats::predict(fit_mars, te_imp)
  mm <- add_metrics("MARS", pred_tr, tr_imp$y, pred_te, te_imp$y)
  results$mars <- fit_mars
  metrics_train$mars <- mm$train
  metrics_test$mars  <- mm$test

  ## -------- kNN (scaled) -----------------------------------------------------
  pp_cs_knn <- caret::preProcess(tr_imp[, -1, drop = FALSE], method = c("center","scale"))
  Xtr_knn <- stats::predict(pp_cs_knn, tr_imp[, -1, drop = FALSE])
  Xte_knn <- stats::predict(pp_cs_knn, te_imp[, -1, drop = FALSE])
  df_tr_knn <- data.frame(y = tr_imp$y, Xtr_knn)
  df_te_knn <- data.frame(y = te_imp$y, Xte_knn)

  knn_grid <- data.frame(k = seq(10, 75, by = 5))
  set.seed(seed)
  fit_knn <- caret::train(
    y ~ ., data = df_tr_knn, method = "knn",
    trControl = ctrl, tuneGrid = knn_grid, metric = "RMSE"
  )
  pred_tr <- stats::predict(fit_knn, df_tr_knn)
  pred_te <- stats::predict(fit_knn, df_te_knn)
  mm <- add_metrics("kNN", pred_tr, df_tr_knn$y, pred_te, df_te_knn$y)
  results$knn <- fit_knn
  metrics_train$knn <- mm$train
  metrics_test$knn  <- mm$test

  ## -------- Cubist -----------------------------------------------------------
  cubist_grid <- expand.grid(
    committees = c(1, 10, 25),
    neighbors  = c(0, 3, 5)
  )
  set.seed(seed)
  fit_cub <- caret::train(
    y ~ ., data = tr_imp, method = "cubist",
    trControl = ctrl, tuneGrid = cubist_grid, metric = "RMSE"
  )
  pred_tr <- stats::predict(fit_cub, tr_imp)
  pred_te <- stats::predict(fit_cub, te_imp)
  mm <- add_metrics("Cubist", pred_tr, tr_imp$y, pred_te, te_imp$y)
  results$cubist <- fit_cub
  metrics_train$cubist <- mm$train
  metrics_test$cubist  <- mm$test

  # --- assemble outputs -------------------------------------------------------
  list(
    train_metrics  = dplyr::bind_rows(metrics_train),
    test_metrics   = dplyr::bind_rows(metrics_test),
    models         = results,
    splits         = list(train_idx = idx_tr),
    preprocessors  = list(
      impute    = pp_imp,
      cs_glmnet = pp_cs,
      cs_svm    = pp_cs_svm,
      cs_knn    = pp_cs_knn
    ),
    feature_filters = list(
      removed_nzv  = if (exists("nzv")) rownames(nzv)[nzv$nzv] else character(0),
      removed_corr = drop_corr,
      corr_cutoff  = corr_cutoff
    ),
    gbm_bag_fraction = safe_bf
  )
}

