#!/usr/bin/env python3
"""Generate duration or growth manuscript CSV tables and figures.

Duration is the default outcome. Use ``--outcome growth`` for growth outputs.
Run without section arguments to regenerate every output, or pass one or more
sections: table2, candidate-windows, country, shap. The selected downstream
candidate is Random Forest for duration and Super Learner for growth.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from time import perf_counter

def generate_table2(outcome: str) -> None:
    """Create the concise manuscript-style model comparison table.

    The table reuses the locked-test-set estimates and lineage-cluster bootstrap
    intervals from the completed duration run.
    """


    from pathlib import Path

    import pandas as pd


    RESULTS_DIR = Path(__file__).resolve().parent / "modeling_results" / outcome
    OUTPUT_DIR = RESULTS_DIR / "tab_figs"

    MODEL_ORDER = {
        "logistic_regression": 0,
        "elastic_net": 1,
        "random_forest": 2,
        "xgboost": 3,
        "svm_rbf": 4,
        "gam": 5,
        "ebm": 6,
        "super_learner": 7,
    }

    MODEL_LABELS = {
        "logistic_regression": "Plain Logistic Regression",
        "elastic_net": "Elastic-net Logistic Regression",
        "random_forest": "Random Forest",
        "xgboost": "XGBoost",
        "svm_rbf": "RBF-SVM",
        "gam": "GAM",
        "ebm": "EBM",
        "super_learner": "Super Learner",
    }

    METRICS = [
        ("AUPRC", "average_precision"),
        ("AUROC", "roc_auc"),
        ("Brier score", "brier_score"),
        ("Sensitivity", "sensitivity"),
        ("PPV", "ppv"),
    ]


    def format_interval(estimate: float, lower: float, upper: float) -> str:
        if pd.isna(estimate):
            return "NA"
        if pd.isna(lower) or pd.isna(upper):
            return f"{estimate:.3f} (NA, NA)"
        return f"{estimate:.3f} ({lower:.3f}, {upper:.3f})"


    def main() -> None:
        metrics = pd.read_csv(RESULTS_DIR / "metrics.csv")
        metrics = metrics.loc[metrics["split"].eq("test")].copy()

        indexed = metrics.set_index(["window_days", "model", "metric"])
        estimate_rows: list[dict[str, object]] = []
        interval_rows: list[dict[str, object]] = []
        combinations = metrics.loc[:, ["window_days", "model"]].drop_duplicates()
        combinations["model_order"] = combinations["model"].map(MODEL_ORDER)
        combinations = combinations.sort_values(["window_days", "model_order"])

        for item in combinations.itertuples(index=False):
            window = int(item.window_days)
            model = str(item.model)
            estimate_row: dict[str, object] = {
                "Input window (days)": window,
                "Model": MODEL_LABELS[model],
            }
            interval_row = estimate_row.copy()
            for output_label, metric_name in METRICS:
                value = indexed.loc[(window, model, metric_name)]
                estimate = float(value["estimate"])
                lower = float(value["ci_95_lower"])
                upper = float(value["ci_95_upper"])
                estimate_row[output_label] = round(estimate, 3)
                interval_row[output_label] = format_interval(estimate, lower, upper)
            estimate_rows.append(estimate_row)
            interval_rows.append(interval_row)

        estimates = pd.DataFrame(estimate_rows)
        intervals = pd.DataFrame(interval_rows)
        estimates["Best metric(s)"] = ""
        intervals["Best metric(s)"] = ""
        for window, group in estimates.groupby("Input window (days)", sort=True):
            winning_metrics: dict[int, list[str]] = {index: [] for index in group.index}
            for output_label, _ in METRICS:
                best = group[output_label].min() if output_label == "Brier score" else group[output_label].max()
                for index in group.index[group[output_label].eq(best)]:
                    winning_metrics[index].append(output_label)
            for index, labels in winning_metrics.items():
                flag = "; ".join(labels)
                estimates.loc[index, "Best metric(s)"] = flag
                intervals.loc[index, "Best metric(s)"] = flag

        metric_columns = [label for label, _ in METRICS]
        if intervals[metric_columns].apply(lambda column: column.str.contains("NA").any()).any():
            raise RuntimeError("At least one requested estimate or confidence interval is missing")
        if len(estimates) != 24:
            raise RuntimeError(f"Expected 24 model/window rows, found {len(estimates)}")

        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        estimate_path = OUTPUT_DIR / "table2_model_performance.csv"
        interval_path = OUTPUT_DIR / "table2_model_performance_with_ci.csv"
        estimates.to_csv(estimate_path, index=False, float_format="%.3f")
        intervals.to_csv(interval_path, index=False)
        print(f"Saved {len(estimates)} rows to {estimate_path}")
        print(f"Saved {len(intervals)} rows to {interval_path}")
    main()


def generate_candidate_windows(outcome: str) -> None:
    """Plot locked-test candidate-model performance across input windows."""


    from pathlib import Path

    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import numpy as np
    import pandas as pd
    from sklearn.calibration import calibration_curve
    from sklearn.metrics import (
        average_precision_score,
        precision_recall_curve,
        roc_auc_score,
        roc_curve,
    )


    RESULTS_DIR = Path(__file__).resolve().parent / "modeling_results" / outcome
    OUTPUT_DIR = RESULTS_DIR / "tab_figs"
    POSITIVE_DISPLAY = "long duration" if outcome == "duration" else "large growth"
    POSITIVE_CLASS = "long" if outcome == "duration" else "large"
    NEGATIVE_CLASS = "short/medium" if outcome == "duration" else "minimal/moderate"
    OUTCOME_TITLE = "Duration" if outcome == "duration" else "Growth"
    MODEL = "random_forest" if outcome == "duration" else "super_learner"
    MODEL_LABEL = "Random Forest" if outcome == "duration" else "Super Learner"
    FILE_STEM = "rf" if outcome == "duration" else "super_learner"
    WINDOWS = (14, 21, 28)
    WINDOW_COLORS = {14: "#0072B2", 21: "#D55E00", 28: "#009E73"}
    WINDOW_MARKERS = {14: "o", 21: "s", 28: "^"}
    CLASS_COLORS = {0: "#56B4E9", 1: "#E69F00"}


    def load_predictions() -> pd.DataFrame:
        predictions = pd.read_csv(RESULTS_DIR / "predictions.csv")
        predictions = predictions.loc[
            predictions["split"].eq("test") & predictions["model"].eq(MODEL)
        ].copy()
        observed_windows = tuple(sorted(predictions["window_days"].unique()))
        if observed_windows != WINDOWS:
            raise RuntimeError(f"Expected windows {WINDOWS}, found {observed_windows}")
        for window, group in predictions.groupby("window_days"):
            if group["observed"].nunique() != 2:
                raise RuntimeError(f"Window {window} does not contain both outcome classes")
            if not group["predicted_probability"].between(0, 1).all():
                raise RuntimeError(f"Window {window} has invalid predicted probabilities")
        return predictions


    def style_axis(axis: plt.Axes) -> None:
        axis.grid(color="#D9D9D9", linewidth=0.7, alpha=0.65)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_axisbelow(True)


    def add_panel_label(axis: plt.Axes, label: str) -> None:
        axis.text(
            -0.12,
            1.06,
            label,
            transform=axis.transAxes,
            fontsize=22,
            fontweight="bold",
            va="top",
        )


    def main() -> None:
        predictions = load_predictions()
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        plt.rcParams.update(
            {
                "font.family": "DejaVu Sans",
                "font.size": 16,
                "axes.titlesize": 19,
                "axes.labelsize": 17,
                "legend.fontsize": 14.5,
                "xtick.labelsize": 15,
                "ytick.labelsize": 15,
            }
        )
        figure, axes = plt.subplots(2, 2, figsize=(15, 11), constrained_layout=True)
        roc_axis, pr_axis, calibration_axis, distribution_axis = axes.ravel()

        summaries: dict[int, dict[str, float]] = {}
        for window in WINDOWS:
            group = predictions.loc[predictions["window_days"].eq(window)]
            y = group["observed"].to_numpy(dtype=int)
            probability = group["predicted_probability"].to_numpy(dtype=float)
            color = WINDOW_COLORS[window]
            marker = WINDOW_MARKERS[window]
            prevalence = float(y.mean())
            auroc = float(roc_auc_score(y, probability))
            auprc = float(average_precision_score(y, probability))
            brier = float(np.mean((probability - y) ** 2))
            summaries[window] = {
                "n": float(len(y)),
                "prevalence": prevalence,
                "auroc": auroc,
                "auprc": auprc,
                "brier": brier,
            }

            fpr, tpr, _ = roc_curve(y, probability)
            roc_axis.plot(
                fpr,
                tpr,
                color=color,
                linewidth=2.3,
                label=f"{window} days (AUROC {auroc:.3f})",
            )

            precision, recall, _ = precision_recall_curve(y, probability)
            pr_axis.plot(
                recall,
                precision,
                color=color,
                linewidth=2.3,
                label=f"{window} days (AUPRC {auprc:.3f})",
            )
            pr_axis.axhline(prevalence, color=color, linestyle="--", linewidth=1.3, alpha=0.65)

            observed_fraction, mean_probability = calibration_curve(
                y,
                probability,
                n_bins=8,
                strategy="quantile",
            )
            calibration_axis.plot(
                mean_probability,
                observed_fraction,
                color=color,
                marker=marker,
                markersize=6,
                linewidth=2.0,
                label=f"{window} days (Brier {brier:.3f})",
            )

        roc_axis.plot([0, 1], [0, 1], color="#4D4D4D", linestyle="--", linewidth=1.2, label="No discrimination")
        roc_axis.set(
            title="ROC curves",
            xlabel="False-positive rate",
            ylabel="True-positive rate (sensitivity)",
            xlim=(0, 1),
            ylim=(0, 1),
        )
        roc_axis.legend(loc="lower right", frameon=False)
        style_axis(roc_axis)
        add_panel_label(roc_axis, "A")

        prevalence_text = "; ".join(
            f"{window}d={summaries[window]['prevalence']:.3f}" for window in WINDOWS
        )
        pr_axis.text(
            0.03,
            0.04,
            f"Dashed prevalence baselines: {prevalence_text}",
            transform=pr_axis.transAxes,
            fontsize=13.5,
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "#BFBFBF", "alpha": 0.88, "pad": 4},
        )
        pr_axis.set(
            title="Precision–recall curves",
            xlabel="Recall (sensitivity)",
            ylabel="Precision (PPV)",
            xlim=(0, 1),
            ylim=(0, 1),
        )
        pr_axis.legend(loc="upper right", frameon=False)
        style_axis(pr_axis)
        add_panel_label(pr_axis, "B")

        calibration_axis.plot([0, 1], [0, 1], color="#4D4D4D", linestyle="--", linewidth=1.2, label="Perfect calibration")
        calibration_axis.set(
            title="Calibration curves (8 quantile bins)",
            xlabel="Mean predicted probability",
            ylabel=f"Observed {POSITIVE_DISPLAY} proportion",
            xlim=(0, 1),
            ylim=(0, 1),
        )
        calibration_axis.legend(loc="upper left", frameon=False)
        style_axis(calibration_axis)
        add_panel_label(calibration_axis, "C")

        positions = np.arange(len(WINDOWS), dtype=float)
        offsets = {0: -0.17, 1: 0.17}
        for class_id in (0, 1):
            values = []
            class_positions = []
            for position, window in zip(positions, WINDOWS):
                group = predictions.loc[
                    predictions["window_days"].eq(window) & predictions["observed"].eq(class_id),
                    "predicted_probability",
                ]
                values.append(group.to_numpy(dtype=float))
                class_positions.append(position + offsets[class_id])
            violin = distribution_axis.violinplot(
                values,
                positions=class_positions,
                widths=0.30,
                showmeans=False,
                showmedians=True,
                showextrema=False,
            )
            for body in violin["bodies"]:
                body.set_facecolor(CLASS_COLORS[class_id])
                body.set_edgecolor("#333333")
                body.set_alpha(0.72)
                body.set_linewidth(0.8)
            violin["cmedians"].set_color("#1A1A1A")
            violin["cmedians"].set_linewidth(1.5)

        distribution_axis.set_xticks(
            positions,
            [f"{window} days\n(n={int(summaries[window]['n'])})" for window in WINDOWS],
        )
        distribution_axis.set(
            title="Predicted-probability distributions by observed class",
            xlabel="Input window",
            ylabel=f"Predicted probability of {POSITIVE_DISPLAY}",
            xlim=(-0.55, len(WINDOWS) - 0.45),
            ylim=(0, 1),
        )
        distribution_axis.legend(
            handles=[
                    Patch(facecolor=CLASS_COLORS[0], edgecolor="#333333", alpha=0.72, label=f"Observed {NEGATIVE_CLASS}"),
                    Patch(facecolor=CLASS_COLORS[1], edgecolor="#333333", alpha=0.72, label=f"Observed {POSITIVE_CLASS}"),
            ],
            loc="upper left",
            frameon=False,
        )
        style_axis(distribution_axis)
        add_panel_label(distribution_axis, "D")

        png_path = OUTPUT_DIR / f"figure_{FILE_STEM}_windows.png"
        pdf_path = OUTPUT_DIR / f"figure_{FILE_STEM}_windows.pdf"
        figure.savefig(png_path, dpi=320, bbox_inches="tight", facecolor="white")
        figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
        plt.close(figure)

        print(f"Saved {png_path}")
        print(f"Saved {pdf_path}")
        for window in WINDOWS:
            values = summaries[window]
            print(
                f"{window}d: n={int(values['n'])}, prevalence={values['prevalence']:.3f}, "
                f"AUROC={values['auroc']:.3f}, AUPRC={values['auprc']:.3f}, Brier={values['brier']:.3f}"
            )
    main()


def generate_country(outcome: str) -> None:
    """Country-specific 28-day model performance on the locked test set."""


    from pathlib import Path

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    import numpy as np
    import pandas as pd
    from sklearn.metrics import average_precision_score, roc_auc_score


    RESULTS_DIR = Path(__file__).resolve().parent / "modeling_results" / outcome
    OUTPUT_DIR = RESULTS_DIR / "tab_figs"
    PREDICTIONS_PATH = RESULTS_DIR / f"{outcome}_28d" / "predictions.csv"
    POSITIVE_CLASS = "long" if outcome == "duration" else "large"
    NEGATIVE_CLASS = "non-long" if outcome == "duration" else "non-large"
    OUTCOME_TITLE = "Duration" if outcome == "duration" else "Growth"
    MODEL = "random_forest" if outcome == "duration" else "super_learner"
    MODEL_LABEL = "Random Forest" if outcome == "duration" else "Super Learner"
    MODEL_ORDER = (
        "logistic_regression",
        "elastic_net",
        "random_forest",
        "xgboost",
        "svm_rbf",
        "gam",
        "ebm",
        "super_learner",
    )
    HEATMAP_MODEL_ORDER = tuple(model for model in MODEL_ORDER if model != "logistic_regression")
    MODEL_LABELS = {
        "logistic_regression": "Plain LR",
        "elastic_net": "Elastic-net LR",
        "random_forest": "RF",
        "xgboost": "XGBoost",
        "svm_rbf": "RBF-SVM",
        "gam": "GAM",
        "ebm": "EBM",
        "super_learner": "Super Learner",
    }
    WINDOW_DAYS = 28
    MIN_CLASS_N = 10
    N_BOOTSTRAP = 2_000
    RANDOM_SEED = 20260820

    PRIMARY_COLOR = "#0072B2"
    LIMITED_COLOR = "#9B9B9B"
    BASELINE_COLOR = "#5A5A5A"


    def load_all_predictions() -> pd.DataFrame:
        predictions = pd.read_csv(PREDICTIONS_PATH)
        predictions = predictions.loc[
            predictions["split"].eq("test")
            & predictions["window_days"].eq(WINDOW_DAYS)
        ].copy()
        if predictions.empty:
            raise RuntimeError("No 28-day locked-test predictions found")
        observed_models = tuple(predictions["model"].drop_duplicates())
        if set(observed_models) != set(MODEL_ORDER):
            raise RuntimeError(f"Expected models {MODEL_ORDER}, found {observed_models}")
        if not predictions["predicted_probability"].between(0, 1).all():
            raise RuntimeError("Predicted probabilities must be between 0 and 1")
        if not predictions.groupby("model")["threshold"].nunique().eq(1).all():
            raise RuntimeError("Expected one locked threshold per model")
        if predictions.duplicated(["model", "country", "lineage"]).any():
            raise RuntimeError("Expected one observation per model-country-lineage tuple")
        return predictions


    def candidate_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
        selected = predictions.loc[predictions["model"].eq(MODEL)].copy()
        if selected.empty:
            raise RuntimeError(f"No 28-day {MODEL_LABEL} locked-test predictions found")
        return selected


    def point_metrics(y: np.ndarray, probability: np.ndarray, threshold: float) -> dict[str, float]:
        predicted = probability >= threshold
        positive_n = int(y.sum())
        negative_n = int((1 - y).sum())
        tp = int(((y == 1) & predicted).sum())
        fp = int(((y == 0) & predicted).sum())
        metrics = {
            "auprc": np.nan,
            "auroc": np.nan,
            "brier": float(np.mean((probability - y) ** 2)),
            "sensitivity": float(tp / positive_n) if positive_n else np.nan,
            "ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        }
        if positive_n and negative_n:
            metrics["auprc"] = float(average_precision_score(y, probability))
            metrics["auroc"] = float(roc_auc_score(y, probability))
        return metrics


    def stratified_bootstrap(
        y: np.ndarray,
        probability: np.ndarray,
        threshold: float,
        rng: np.random.Generator,
    ) -> dict[str, tuple[float, float]]:
        """Percentile CIs from outcome-stratified resampling within a country."""
        class_indices = {value: np.flatnonzero(y == value) for value in (0, 1)}
        values: dict[str, list[float]] = {
            "auprc": [],
            "auroc": [],
            "brier": [],
            "sensitivity": [],
            "ppv": [],
        }
        for _ in range(N_BOOTSTRAP):
            pieces = [
                rng.choice(indices, size=len(indices), replace=True)
                for indices in class_indices.values()
                if len(indices)
            ]
            sampled = np.concatenate(pieces)
            bootstrap_metrics = point_metrics(y[sampled], probability[sampled], threshold)
            for metric, value in bootstrap_metrics.items():
                if np.isfinite(value):
                    values[metric].append(value)

        intervals: dict[str, tuple[float, float]] = {}
        for metric, metric_values in values.items():
            if metric_values:
                lower, upper = np.quantile(metric_values, [0.025, 0.975])
                intervals[metric] = (float(lower), float(upper))
            else:
                intervals[metric] = (np.nan, np.nan)
        return intervals


    def summarize_by_country(predictions: pd.DataFrame) -> pd.DataFrame:
        rng = np.random.default_rng(RANDOM_SEED)
        threshold = float(predictions["threshold"].iloc[0])
        records: list[dict[str, float | int | str | bool]] = []
        for country, group in predictions.groupby("country", sort=True):
            y = group["observed"].to_numpy(dtype=int)
            probability = group["predicted_probability"].to_numpy(dtype=float)
            positive_n = int(y.sum())
            negative_n = int((1 - y).sum())
            estimates = point_metrics(y, probability, threshold)
            intervals = stratified_bootstrap(y, probability, threshold, rng)
            record: dict[str, float | int | str | bool] = {
                "model": MODEL,
                "model_label": MODEL_LABEL,
                "country": country,
                "n": len(group),
                "positive_n": positive_n,
                "negative_n": negative_n,
                "prevalence": float(y.mean()),
                "primary_interpretation": positive_n >= MIN_CLASS_N and negative_n >= MIN_CLASS_N,
                "threshold": threshold,
            }
            for metric, estimate in estimates.items():
                record[metric] = estimate
                record[f"{metric}_ci_95_lower"] = intervals[metric][0]
                record[f"{metric}_ci_95_upper"] = intervals[metric][1]
            records.append(record)
        return pd.DataFrame.from_records(records)


    def country_order(summary: pd.DataFrame) -> list[str]:
        primary = summary.loc[summary["primary_interpretation"], "country"].sort_values().tolist()
        limited = summary.loc[~summary["primary_interpretation"], "country"].sort_values().tolist()
        return primary + limited


    def style_axis(axis: plt.Axes) -> None:
        axis.grid(axis="x", color="#D9D9D9", linewidth=0.7, alpha=0.75)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
        axis.set_axisbelow(True)


    def plot_country_forest(summary: pd.DataFrame) -> tuple[Path, Path]:
        ordered_countries = country_order(summary)
        plotting = summary.set_index("country").loc[ordered_countries].reset_index()
        y_positions = np.arange(len(plotting))

        plt.rcParams.update(
            {
                "font.family": "DejaVu Sans",
                "font.size": 10,
                "axes.titlesize": 12,
                "axes.labelsize": 10.5,
                "xtick.labelsize": 9.5,
                "ytick.labelsize": 9.5,
            }
        )
        figure, axes = plt.subplots(1, 3, figsize=(16.5, 9.5), sharey=True)
        panels = (
            ("auprc", "AUPRC", "AUPRC (higher is better)"),
            ("auroc", "AUROC", "AUROC (higher is better)"),
            ("brier", "Brier score", "Brier score (lower is better)"),
        )

        for panel_index, (axis, (metric, title, xlabel)) in enumerate(zip(axes, panels)):
            for row_index, row in plotting.iterrows():
                estimate = float(row[metric])
                if not np.isfinite(estimate):
                    axis.text(
                        0.02,
                        y_positions[row_index],
                        "NE",
                        color=LIMITED_COLOR,
                        va="center",
                        fontsize=9,
                        fontstyle="italic",
                    )
                    continue
                lower = float(row[f"{metric}_ci_95_lower"])
                upper = float(row[f"{metric}_ci_95_upper"])
                primary = bool(row["primary_interpretation"])
                color = PRIMARY_COLOR if primary else LIMITED_COLOR
                axis.errorbar(
                    estimate,
                    y_positions[row_index],
                    xerr=np.array([[estimate - lower], [upper - estimate]]),
                    fmt="o",
                    markersize=6.2,
                    markerfacecolor=color if primary else "white",
                    markeredgecolor=color,
                    markeredgewidth=1.3,
                    ecolor=color,
                    elinewidth=1.35,
                    capsize=2.5,
                    zorder=3,
                )

            if metric == "auprc":
                baseline_mask = plotting["prevalence"].between(0, 1, inclusive="neither")
                axis.scatter(
                    plotting.loc[baseline_mask, "prevalence"],
                    y_positions[baseline_mask.to_numpy()],
                    marker="x",
                    s=34,
                    linewidths=1.2,
                    color=BASELINE_COLOR,
                    alpha=0.85,
                    zorder=2,
                )
                axis.set_xlim(0, 1)
            elif metric == "auroc":
                axis.axvline(0.5, color=BASELINE_COLOR, linestyle="--", linewidth=1.1)
                axis.set_xlim(0, 1)
            else:
                null_brier = plotting["prevalence"] * (1 - plotting["prevalence"])
                baseline_mask = plotting["prevalence"].between(0, 1, inclusive="neither")
                axis.scatter(
                    null_brier.loc[baseline_mask],
                    y_positions[baseline_mask.to_numpy()],
                    marker="x",
                    s=34,
                    linewidths=1.2,
                    color=BASELINE_COLOR,
                    alpha=0.85,
                    zorder=2,
                )
                upper_limit = max(
                    0.35,
                    float(np.nanmax(plotting["brier_ci_95_upper"])) + 0.035,
                    float(null_brier.max()) + 0.035,
                )
                axis.set_xlim(0, min(0.55, upper_limit))

            axis.set_title(title, fontweight="bold", pad=10)
            axis.set_xlabel(xlabel)
            axis.set_ylim(len(plotting) - 0.5, -0.5)
            style_axis(axis)
            axis.text(
                -0.11,
                1.035,
                chr(ord("A") + panel_index),
                transform=axis.transAxes,
                fontsize=13,
                fontweight="bold",
                va="top",
            )

        labels = [
            f"{row.country}  ({int(row.positive_n)}/{int(row.negative_n)})"
            for row in plotting.itertuples()
        ]
        axes[0].set_yticks(y_positions, labels)
        for label, primary in zip(axes[0].get_yticklabels(), plotting["primary_interpretation"]):
            label.set_fontweight("bold" if primary else "normal")
            label.set_color("#1F1F1F" if primary else LIMITED_COLOR)

        primary_count = int(plotting["primary_interpretation"].sum())
        if 0 < primary_count < len(plotting):
            separator_y = primary_count - 0.5
            for axis in axes:
                axis.axhline(separator_y, color="#BFBFBF", linewidth=0.9)

        legend_handles = [
            Line2D(
                [0], [0], marker="o", color=PRIMARY_COLOR, markerfacecolor=PRIMARY_COLOR,
                linestyle="none", markersize=6.5,
                label=f"≥10 {POSITIVE_CLASS} and ≥10 {NEGATIVE_CLASS}",
            ),
            Line2D(
                [0], [0], marker="o", color=LIMITED_COLOR, markerfacecolor="white",
                linestyle="none", markersize=6.5, label="Limited class counts",
            ),
            Line2D(
                [0], [0], marker="x", color=BASELINE_COLOR, linestyle="none",
                markersize=7, label="Prevalence/null Brier baseline",
            ),
            Line2D(
                [0], [0], color=BASELINE_COLOR, linestyle="--", linewidth=1.1,
                label="AUROC no-skill value (0.5)",
            ),
        ]
        figure.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=4,
            frameon=False,
            bbox_to_anchor=(0.5, 0.02),
        )
        figure.subplots_adjust(left=0.19, right=0.985, top=0.96, bottom=0.12, wspace=0.16)

        png_path = OUTPUT_DIR / "figure_country_performance_28d.png"
        pdf_path = OUTPUT_DIR / "figure_country_performance_28d.pdf"
        figure.savefig(png_path, dpi=320, bbox_inches="tight", facecolor="white")
        figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
        plt.close(figure)
        return png_path, pdf_path


    def summarize_country_auroc(
        predictions: pd.DataFrame,
        eligible_summary: pd.DataFrame,
    ) -> pd.DataFrame:
        eligible = eligible_summary.loc[
            eligible_summary["primary_interpretation"],
            ["country", "positive_n", "negative_n"],
        ]
        eligible_countries = set(eligible["country"])
        records: list[dict[str, float | int | str]] = []
        for (model, country), group in predictions.loc[
            predictions["country"].isin(eligible_countries)
            & predictions["model"].isin(HEATMAP_MODEL_ORDER)
        ].groupby(["model", "country"], sort=False):
            y = group["observed"].to_numpy(dtype=int)
            probability = group["predicted_probability"].to_numpy(dtype=float)
            counts = eligible.loc[eligible["country"].eq(country)].iloc[0]
            records.append(
                {
                    "country": country,
                    "positive_n": int(counts["positive_n"]),
                    "negative_n": int(counts["negative_n"]),
                    "model": model,
                    "model_label": MODEL_LABELS[model].replace("\n", " "),
                    "auroc": float(roc_auc_score(y, probability)),
                }
            )
        summary = pd.DataFrame.from_records(records)
        expected_rows = len(eligible_countries) * len(HEATMAP_MODEL_ORDER)
        if len(summary) != expected_rows:
            raise RuntimeError(f"Expected {expected_rows} country-model rows, found {len(summary)}")
        return summary


    def plot_auroc_heatmap(auroc_summary: pd.DataFrame) -> tuple[Path, Path]:
        countries = sorted(auroc_summary["country"].unique())
        counts = auroc_summary.drop_duplicates("country").set_index("country")
        country_labels = [
            f"{country}  ({int(counts.loc[country, 'positive_n'])}/{int(counts.loc[country, 'negative_n'])})"
            for country in countries
        ]

        values = (
            auroc_summary.pivot(index="country", columns="model", values="auroc")
            .reindex(index=countries, columns=HEATMAP_MODEL_ORDER)
            .to_numpy(dtype=float)
        )
        figure, axis = plt.subplots(figsize=(13.5, 7.4))
        image = axis.imshow(values, cmap="YlGn", vmin=0.5, vmax=1.0, aspect="auto")
        axis.set_xticks(
            np.arange(len(HEATMAP_MODEL_ORDER)),
            [MODEL_LABELS[model] for model in HEATMAP_MODEL_ORDER],
        )
        axis.set_yticks(np.arange(len(countries)), country_labels)
        axis.tick_params(axis="both", length=0)
        axis.xaxis.tick_top()
        for tick, model in zip(axis.get_xticklabels(), HEATMAP_MODEL_ORDER):
            tick.set_fontweight("bold" if model == MODEL else "normal")
        for row_index in range(values.shape[0]):
            row_maximum = np.nanmax(values[row_index])
            for column_index in range(values.shape[1]):
                value = values[row_index, column_index]
                red, green, blue, _ = image.cmap(image.norm(value))
                luminance = 0.299 * red + 0.587 * green + 0.114 * blue
                is_row_maximum = np.isclose(value, row_maximum)
                axis.text(
                    column_index,
                    row_index,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    color="white" if luminance < 0.52 else "#1A1A1A",
                    fontsize=10,
                    fontweight="bold" if is_row_maximum else "normal",
                )
                if is_row_maximum:
                    axis.add_patch(
                        Rectangle(
                            (column_index - 0.5, row_index - 0.5),
                            1,
                            1,
                            fill=False,
                            edgecolor="#1A1A1A",
                            linewidth=2.1,
                        )
                    )
        for spine in axis.spines.values():
            spine.set_visible(False)

        colorbar_axis = figure.add_axes([0.91, 0.15, 0.014, 0.70])
        colorbar = figure.colorbar(image, cax=colorbar_axis)
        colorbar.set_label("AUROC (higher is better)", rotation=270, labelpad=18)
        colorbar.outline.set_visible(False)
        figure.subplots_adjust(left=0.18, right=0.88, top=0.88, bottom=0.08)

        png_path = OUTPUT_DIR / "figure_country_auroc_28d.png"
        pdf_path = OUTPUT_DIR / "figure_country_auroc_28d.pdf"
        figure.savefig(png_path, dpi=320, bbox_inches="tight", facecolor="white")
        figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
        plt.close(figure)
        return png_path, pdf_path


    def main() -> None:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        all_predictions = load_all_predictions()
        predictions = candidate_predictions(all_predictions)
        summary = summarize_by_country(predictions)
        csv_path = OUTPUT_DIR / "country_performance_28d.csv"
        summary.to_csv(csv_path, index=False, float_format="%.6f")
        forest_paths = plot_country_forest(summary)
        auroc_summary = summarize_country_auroc(all_predictions, summary)
        auroc_csv_path = OUTPUT_DIR / "country_auroc_all_models_28d.csv"
        auroc_summary.to_csv(auroc_csv_path, index=False, float_format="%.6f")
        heatmap_paths = plot_auroc_heatmap(auroc_summary)

        qualifying = summary.loc[summary["primary_interpretation"], "country"].tolist()
        print(f"Saved {csv_path}")
        print(f"Saved {auroc_csv_path}")
        for path in (*forest_paths, *heatmap_paths):
            print(f"Saved {path}")
        print(f"Primary-interpretation countries ({len(qualifying)}): {', '.join(qualifying)}")
    main()


def generate_shap(outcome: str) -> None:
    """SHAP summary plots for the fitted 28-day downstream candidate."""


    import json
    from pathlib import Path
    import warnings

    import joblib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import rdata
    import shap


    CODE_DIR = Path(__file__).resolve().parent
    RESULTS_DIR = CODE_DIR / "modeling_results" / outcome
    MODEL_DIR = RESULTS_DIR / f"{outcome}_28d"
    OUTPUT_DIR = RESULTS_DIR / "tab_figs"
    POSITIVE_DISPLAY = "long duration" if outcome == "duration" else "large growth"
    OUTCOME_TITLE = "Duration" if outcome == "duration" else "Growth"
    MODEL = "random_forest" if outcome == "duration" else "super_learner"
    MODEL_LABEL = "Random Forest" if outcome == "duration" else "Super Learner"
    FILE_STEM = "rf" if outcome == "duration" else "super_learner"
    TEST_DATA_PATH = CODE_DIR / "test.rds"
    TRAIN_DATA_PATH = CODE_DIR / "train.rds"
    PREDICTIONS_PATH = MODEL_DIR / "predictions.csv"
    SELECTION_PATH = MODEL_DIR / "selection.json"
    FEATURE_MAP_PATH = CODE_DIR / "feature_variable_map.txt"
    FEATURES = tuple(f"feature_{index:02d}" for index in range(1, 29))
    MAX_DISPLAY = 15
    N_PERMUTATIONS = 5
    BACKGROUND_SIZE = 80
    SHAP_METHOD = (
        "TreeSHAP"
        if MODEL == "random_forest"
        else f"Permutation SHAP ({N_PERMUTATIONS} permutations; {BACKGROUND_SIZE}-row training background)"
    )
    SHAP_DISPLAY = "TreeSHAP" if MODEL == "random_forest" else "Permutation SHAP"
    DISPLAY_NAME_BY_ORIGINAL = {
        "slope_log": "log_slope",
        "slope_y_1_7": "early_slope_1_7d",
        "peak_val": "peak_share",
        "auc_norm": "mean_share",
        "longest_inc_run": "longest_increase",
        "zero_count": "zero_days",
        "c22_DN_HistogramMode_5": "hist_mode_5",
        "c22_DN_HistogramMode_10": "hist_mode_10",
        "c22_CO_f1ecac": "acf_1e_lag",
        "c22_CO_FirstMin_ac": "acf_first_min",
        "c22_CO_HistogramAMI_even_2_5": "auto_MI_lag2",
        "c22_CO_trev_1_num": "time_reversibility",
        "c22_MD_hrv_classic_pnn40": "pNN40",
        "c22_SB_BinaryStats_mean_longstretch1": "above_mean_run",
        "c22_SB_TransitionMatrix_3ac_sumdiagcov": "transition_matrix",
        "c22_PD_PeriodicityWang_th0_01": "periodicity",
        "c22_CO_Embed2_Dist_tau_d_expfit_meandiff": "embedding_distance",
        "c22_IN_AutoMutualInfoStats_40_gaussian_fmmi": "auto_MI_first_min",
        "c22_FC_LocalSimple_mean1_tauresrat": "forecast_tau_change",
        "c22_DN_OutlierInclude_p_001_mdrmd": "outlier_timing_pos",
        "c22_DN_OutlierInclude_n_001_mdrmd": "outlier_timing_neg",
        "c22_SP_Summaries_welch_rect_area_5_1": "low_freq_power",
        "c22_SB_BinaryStats_diff_longstretch0": "nonpositive_diff_run",
        "c22_SB_MotifThree_quantile_hh": "motif_entropy",
        "c22_SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1": "RS_crossover",
        "c22_SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1": "DFA_crossover",
        "c22_SP_Summaries_welch_rect_centroid": "spectral_centroid",
        "c22_FC_LocalSimple_mean3_stderr": "forecast_error",
    }


    def load_feature_mapping() -> pd.DataFrame:
        records: list[dict[str, str]] = []
        for line in FEATURE_MAP_PATH.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) >= 4 and fields[0] in FEATURES:
                feature, original_variable, source, meaning = fields[:4]
                records.append(
                    {
                        "feature": feature,
                        "original_variable": original_variable,
                        "display_name": DISPLAY_NAME_BY_ORIGINAL.get(original_variable, original_variable),
                        "source": source,
                        "meaning": meaning,
                    }
                )
        mapping = pd.DataFrame.from_records(records)
        if set(mapping.get("feature", [])) != set(FEATURES) or len(mapping) != len(FEATURES):
            raise RuntimeError("feature_variable_map.txt must map each of feature_01–feature_28 exactly once")
        if mapping["original_variable"].duplicated().any() or mapping["display_name"].duplicated().any():
            raise RuntimeError("Feature mapping contains duplicate original or display names")
        missing_short_names = set(mapping["original_variable"]) - set(DISPLAY_NAME_BY_ORIGINAL)
        if missing_short_names:
            raise RuntimeError(f"Missing simplified display names for: {sorted(missing_short_names)}")
        return mapping.set_index("feature").loc[list(FEATURES)].reset_index()


    def load_features(path: Path, key: str) -> pd.DataFrame:
        warnings.filterwarnings("ignore", message="Missing constructor")
        warnings.filterwarnings("ignore", message="Unknown constructor")
        datasets = rdata.read_rds(path)
        data = datasets[key]
        missing = [feature for feature in FEATURES if feature not in data.columns]
        if missing:
            raise RuntimeError(f"Missing features in {key}: {missing}")
        features = data.loc[:, FEATURES].astype(float)
        if not np.isfinite(features.to_numpy()).all():
            raise RuntimeError(f"{key} features contain missing or non-finite values")
        return features


    def load_and_validate_predictor(features: pd.DataFrame):
        array = features.to_numpy()
        model_for_tree_shap = None
        weight_description = ""
        if MODEL == "random_forest":
            model = joblib.load(MODEL_DIR / "random_forest.joblib")
            if getattr(model, "n_features_in_", None) != len(FEATURES):
                raise RuntimeError("Saved Random Forest does not match the 28 expected features")

            def predict_probability(values: np.ndarray) -> np.ndarray:
                return model.predict_proba(np.asarray(values, dtype=float))[:, 1]

            model_for_tree_shap = model
        else:
            selection = json.loads(SELECTION_PATH.read_text(encoding="utf-8"))["super_learner"]
            weights = {name: float(weight) for name, weight in selection["weights"].items()}
            if any(weight < 0 for weight in weights.values()) or not np.isclose(sum(weights.values()), 1.0):
                raise RuntimeError(f"Invalid Super Learner weights: {weights}")
            weighted_models = []
            for name, weight in weights.items():
                if weight <= 1e-12:
                    continue
                fitted = joblib.load(MODEL_DIR / f"{name}.joblib")
                if getattr(fitted, "n_features_in_", len(FEATURES)) != len(FEATURES):
                    raise RuntimeError(f"Saved {name} model does not match the 28 expected features")
                weighted_models.append((name, weight, fitted))

            def predict_probability(values: np.ndarray) -> np.ndarray:
                values = np.asarray(values, dtype=float)
                probability = np.zeros(len(values), dtype=float)
                for _, weight, fitted in weighted_models:
                    probability += weight * fitted.predict_proba(values)[:, 1]
                return probability

            weight_description = ", ".join(
                f"{name}={weight:.3f}" for name, weight, _ in weighted_models
            )

        probability = predict_probability(array)
        saved = pd.read_csv(PREDICTIONS_PATH)
        saved = saved.loc[
            saved["split"].eq("test") & saved["model"].eq(MODEL)
        ].sort_values("row_index")
        if len(saved) != len(features):
            raise RuntimeError("Saved test predictions do not match the test feature rows")
        maximum_difference = float(
            np.max(np.abs(probability - saved["predicted_probability"].to_numpy(dtype=float)))
        )
        tolerance = 1e-10 if MODEL == "random_forest" else 1e-6
        if maximum_difference > tolerance:
            raise RuntimeError(
                f"Reconstructed {MODEL_LABEL} predictions differ from the locked results: "
                f"max absolute difference={maximum_difference:.3g}"
            )
        return predict_probability, model_for_tree_shap, maximum_difference, weight_description


    def positive_class_explanation(
        predict_probability,
        model_for_tree_shap,
        features: pd.DataFrame,
    ) -> tuple[shap.Explanation, float]:
        if MODEL == "random_forest":
            explainer = shap.TreeExplainer(
                model_for_tree_shap,
                feature_perturbation="tree_path_dependent",
                model_output="raw",
                feature_names=list(FEATURES),
            )
            explanation = explainer(features, check_additivity=False)
        else:
            training = load_features(TRAIN_DATA_PATH, "train_28")
            background = training.sample(n=min(BACKGROUND_SIZE, len(training)), random_state=20260820)
            masker = shap.maskers.Independent(background.to_numpy(), max_samples=len(background))
            explainer = shap.Explainer(
                predict_probability,
                masker=masker,
                algorithm="permutation",
                feature_names=list(FEATURES),
                seed=20260820,
            )
            explanation = explainer(
                features.to_numpy(),
                max_evals=N_PERMUTATIONS * (2 * len(FEATURES) + 1),
                batch_size=64,
                silent=False,
            )
        if explanation.values.ndim == 3:
            explanation = explanation[:, :, 1]
        if explanation.values.shape != features.shape:
            raise RuntimeError(
                f"Unexpected SHAP shape {explanation.values.shape}; expected {features.shape}"
            )

        reconstructed = (
            np.asarray(explanation.base_values, dtype=float).reshape(-1)
            + np.asarray(explanation.values, dtype=float).sum(axis=1)
        )
        probability = predict_probability(features.to_numpy())
        additivity_error = float(np.max(np.abs(reconstructed - probability)))
        additivity_tolerance = 1e-6 if MODEL == "random_forest" else 1e-5
        if additivity_error > additivity_tolerance:
            raise RuntimeError(f"SHAP additivity check failed: max error={additivity_error:.3g}")
        return explanation, additivity_error


    def with_display_names(
        explanation: shap.Explanation,
        mapping: pd.DataFrame,
    ) -> shap.Explanation:
        return shap.Explanation(
            values=np.asarray(explanation.values),
            base_values=np.asarray(explanation.base_values),
            data=np.asarray(explanation.data),
            feature_names=mapping["display_name"].tolist(),
        )


    def save_importance(explanation: shap.Explanation, mapping: pd.DataFrame) -> Path:
        importance = np.mean(np.abs(explanation.values), axis=0)
        summary = mapping.copy()
        summary["mean_abs_shap"] = importance
        summary = summary.sort_values("mean_abs_shap", ascending=False, ignore_index=True)
        summary.insert(0, "shap_method", SHAP_METHOD)
        summary.insert(0, "model_label", MODEL_LABEL)
        summary.insert(0, "model", MODEL)
        summary.insert(0, "rank", np.arange(1, len(summary) + 1))
        path = OUTPUT_DIR / f"{FILE_STEM}_shap_importance_28d.csv"
        summary.to_csv(path, index=False, float_format="%.8f")
        return path


    def plot_summary(explanation: shap.Explanation) -> tuple[Path, Path]:
        plt.rcParams.update(
            {
                "font.family": "DejaVu Sans",
                "font.size": 10,
                "axes.titlesize": 12,
                "axes.labelsize": 10.5,
                "xtick.labelsize": 9.5,
                "ytick.labelsize": 9.5,
            }
        )
        figure, (bar_axis, beeswarm_axis) = plt.subplots(
            1,
            2,
            figsize=(17.5, 9.6),
            gridspec_kw={"width_ratios": [1.0, 1.22]},
        )

        shap.plots.bar(
            explanation,
            max_display=MAX_DISPLAY,
            ax=bar_axis,
            show=False,
        )
        shap.plots.beeswarm(
            explanation,
            max_display=MAX_DISPLAY,
            ax=beeswarm_axis,
            show=False,
            plot_size=None,
            color_bar=True,
            s=18,
        )

        bar_axis.set_title("Global feature importance", fontweight="bold", pad=12)
        beeswarm_axis.set_title("Direction and distribution of effects", fontweight="bold", pad=12)
        bar_axis.text(
            -0.12,
            1.04,
            "A",
            transform=bar_axis.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
        )
        beeswarm_axis.text(
            -0.12,
            1.04,
            "B",
            transform=beeswarm_axis.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
        )
        figure.subplots_adjust(left=0.10, right=0.95, top=0.96, bottom=0.06, wspace=0.54)

        png_path = OUTPUT_DIR / f"figure_{FILE_STEM}_shap_28d.png"
        pdf_path = OUTPUT_DIR / f"figure_{FILE_STEM}_shap_28d.pdf"
        figure.savefig(png_path, dpi=320, bbox_inches="tight", facecolor="white")
        figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
        plt.close(figure)
        return png_path, pdf_path


    def main() -> None:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        features = load_features(TEST_DATA_PATH, "test_28")
        mapping = load_feature_mapping()
        predictor, tree_model, prediction_difference, weight_description = load_and_validate_predictor(features)
        explanation, additivity_error = positive_class_explanation(predictor, tree_model, features)
        importance_path = save_importance(explanation, mapping)
        figure_paths = plot_summary(with_display_names(explanation, mapping))
        print(f"Explained {len(features)} locked-test observations with {len(FEATURES)} features")
        if weight_description:
            print(f"Super Learner nonzero weights: {weight_description}")
        print(f"Prediction verification max absolute difference: {prediction_difference:.3g}")
        print(f"SHAP additivity max absolute error: {additivity_error:.3g}")
        print(f"Saved {importance_path}")
        for path in figure_paths:
            print(f"Saved {path}")
    main()


SECTION_RUNNERS: dict[str, Callable[[str], None]] = {
    "table2": generate_table2,
    "candidate-windows": generate_candidate_windows,
    "country": generate_country,
    "shap": generate_shap,
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate duration or growth manuscript CSV tables and figures."
    )
    parser.add_argument(
        "--outcome",
        choices=("duration", "growth"),
        default="duration",
        help="Response to generate; default: duration.",
    )
    parser.add_argument(
        "sections",
        nargs="*",
        choices=tuple(SECTION_RUNNERS),
        help="Sections to generate; default: all sections.",
    )
    args = parser.parse_args()
    sections = args.sections or list(SECTION_RUNNERS)
    total_start = perf_counter()
    for section in sections:
        start = perf_counter()
        print(f"\n[{section}] started")
        SECTION_RUNNERS[section](args.outcome)
        print(f"[{section}] completed in {perf_counter() - start:.1f}s")
    print(f"\nCompleted {len(sections)} section(s) in {perf_counter() - total_start:.1f}s")


if __name__ == "__main__":
    main()
