#!/usr/bin/env python3
"""Tune and compare binary classifiers for variant duration and growth.

Outcomes
--------
* duration: Long = 1; Short/Medium = 0
* growth:   Large = 1; Minimal/Moderate = 0

For each outcome and each input window (14, 21, and 28 days), models are tuned
with the same lineage-grouped, stratified five-fold cross-validation splits.
The held-out test set is used only after model/parameter/threshold selection.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
import time
import threading
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdata
from interpret.glassbox import ExplainableBoostingClassifier
from pygam import LogisticGAM, s
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import loguniform, randint, uniform
from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.calibration import calibration_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import FitFailedWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_recall_fscore_support,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    GridSearchCV,
    RandomizedSearchCV,
    StratifiedGroupKFold,
    cross_val_predict,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

warnings.filterwarnings("ignore", message="Missing constructor")
warnings.filterwarnings("ignore", message="Unknown constructor")
warnings.filterwarnings("ignore", category=FitFailedWarning)

PLOT_LOCK = threading.Lock()


SEED = 20260820
WINDOWS = (14, 21, 28)
FEATURES = tuple(f"feature_{i:02d}" for i in range(1, 29))
OUTCOMES = {
    "duration": {
        "label": "label_1",
        "positive": "long",
        "class_names": {0: "short_or_medium", 1: "long"},
        "definition": "Long = 1; Short/Medium = 0",
    },
    "growth": {
        "label": "label_2",
        "positive": "large",
        "class_names": {0: "minimal_or_moderate", 1: "large"},
        "definition": "Large = 1; Minimal/Moderate = 0",
    },
}
BASE_MODELS = (
    "logistic_regression",
    "elastic_net",
    "random_forest",
    "xgboost",
    "svm_rbf",
    "gam",
    "ebm",
)
ALL_MODELS = BASE_MODELS + ("super_learner",)


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def stable_seed(*parts: Any) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()
    return (SEED + int(digest[:8], 16)) % (2**31 - 1)


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


class LogisticGAMClassifier(ClassifierMixin, BaseEstimator):
    """Small sklearn-compatible wrapper around pyGAM's LogisticGAM."""

    def __init__(self, n_splines: int = 5, lam: float = 1.0, max_iter: int = 100):
        self.n_splines = n_splines
        self.lam = lam
        self.max_iter = max_iter

    def fit(self, X: np.ndarray, y: np.ndarray) -> "LogisticGAMClassifier":
        X = np.asarray(X, dtype=float)
        terms = s(0, n_splines=self.n_splines)
        for column in range(1, X.shape[1]):
            terms += s(column, n_splines=self.n_splines)
        self.model_ = LogisticGAM(
            terms=terms,
            lam=self.lam,
            max_iter=self.max_iter,
            verbose=False,
        ).fit(X, np.asarray(y, dtype=int))
        self.classes_ = np.array([0, 1], dtype=int)
        self.n_features_in_ = X.shape[1]
        return self

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        positive = np.asarray(self.model_.predict_proba(np.asarray(X, dtype=float)))
        return np.column_stack([1.0 - positive, positive])

    def predict(self, X: np.ndarray) -> np.ndarray:
        return (self.predict_proba(X)[:, 1] >= 0.5).astype(int)


@dataclass
class ModelSpec:
    estimator: BaseEstimator
    search_kind: str
    search_space: Any
    n_iter: int | None = None
    # Search candidates run serially. This avoids a reproducible macOS crash
    # caused by mixing joblib's process workers with XGBoost/OpenMP. Individual
    # estimators may still use their own safe thread-level parallelism.
    n_jobs: int = 1


def model_spec(name: str, seed: int, positive_rate: float, smoke: bool, model_jobs: int = 1) -> ModelSpec:
    if name == "logistic_regression":
        estimator = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        penalty=None,
                        solver="lbfgs",
                        max_iter=10_000,
                        random_state=seed,
                    ),
                ),
            ]
        )
        return ModelSpec(estimator, "none", {})

    if name == "elastic_net":
        estimator = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        penalty="elasticnet",
                        solver="saga",
                        max_iter=10_000,
                        random_state=seed,
                    ),
                ),
            ]
        )
        space = {
            "model__C": [0.03, 0.1, 0.3, 1.0, 3.0] if not smoke else [0.1, 1.0],
            "model__l1_ratio": [0.0, 0.25, 0.5, 0.75, 1.0] if not smoke else [0.0, 0.5],
            "model__class_weight": [None, "balanced"] if not smoke else [None],
        }
        return ModelSpec(estimator, "grid", space)

    if name == "random_forest":
        estimator = RandomForestClassifier(
            random_state=seed,
            n_jobs=model_jobs,
        )
        space = {
            "n_estimators": randint(300, 1001),
            "max_depth": [None, 4, 6, 8, 12, 16],
            "min_samples_split": randint(2, 21),
            "min_samples_leaf": randint(1, 16),
            "max_features": ["sqrt", "log2", 0.5, 0.75, 1.0],
            "class_weight": [None, "balanced", "balanced_subsample"],
        }
        return ModelSpec(estimator, "random", space, n_iter=24 if not smoke else 2)

    if name == "xgboost":
        estimator = XGBClassifier(
            objective="binary:logistic",
            eval_metric="logloss",
            tree_method="hist",
            random_state=seed,
            n_jobs=model_jobs,
        )
        ratio = (1.0 - positive_rate) / max(positive_rate, 1e-9)
        space = {
            "n_estimators": randint(150, 901),
            "learning_rate": loguniform(0.01, 0.2),
            "max_depth": randint(2, 9),
            "min_child_weight": loguniform(0.5, 20.0),
            "subsample": uniform(0.65, 0.35),
            "colsample_bytree": uniform(0.6, 0.4),
            "gamma": uniform(0.0, 3.0),
            "reg_alpha": loguniform(1e-4, 10.0),
            "reg_lambda": loguniform(0.1, 30.0),
            "scale_pos_weight": [1.0, ratio],
        }
        return ModelSpec(estimator, "random", space, n_iter=30 if not smoke else 2)

    if name == "svm_rbf":
        estimator = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                ("model", SVC(kernel="rbf", probability=True, random_state=seed)),
            ]
        )
        space = {
            "model__C": [0.03, 0.1, 0.3, 1, 3, 10, 30] if not smoke else [0.3, 3],
            "model__gamma": ["scale", 0.003, 0.01, 0.03, 0.1, 0.3] if not smoke else ["scale"],
            "model__class_weight": [None, "balanced"] if not smoke else [None],
        }
        return ModelSpec(estimator, "grid", space)

    if name == "gam":
        estimator = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                ("model", LogisticGAMClassifier()),
            ]
        )
        space = {
            "model__n_splines": [4, 5, 6, 8] if not smoke else [4],
            "model__lam": [0.03, 0.1, 0.3, 1.0, 3.0, 10.0] if not smoke else [1.0],
        }
        return ModelSpec(estimator, "grid", space, n_jobs=1)

    if name == "ebm":
        estimator = ExplainableBoostingClassifier(
            max_rounds=5_000,
            early_stopping_rounds=50,
            validation_size=0.15,
            n_jobs=model_jobs,
            random_state=seed,
        )
        space = {
            "learning_rate": loguniform(0.005, 0.05),
            "max_bins": [64, 128, 256],
            "max_leaves": randint(2, 9),
            "min_samples_leaf": randint(2, 21),
            "interactions": [0, 3, 5],
            "outer_bags": [8, 12, 14],
        }
        return ModelSpec(estimator, "random", space, n_iter=20 if not smoke else 2, n_jobs=1)

    raise ValueError(f"Unknown model: {name}")


def validate_data(train: pd.DataFrame, test: pd.DataFrame, label: str) -> None:
    required = {"id_1", "id_2", label, *FEATURES}
    for split_name, frame in (("train", train), ("test", test)):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{split_name} is missing columns: {sorted(missing)}")
        X = frame.loc[:, FEATURES].to_numpy(dtype=float)
        if not np.isfinite(X).all():
            raise ValueError(f"{split_name} contains missing or non-finite feature values")
    overlap = set(train["id_2"].astype(str)).intersection(test["id_2"].astype(str))
    if overlap:
        raise ValueError(f"Train/test lineage leakage detected ({len(overlap)} overlapping lineages)")


def make_binary(frame: pd.DataFrame, label: str, positive: str) -> np.ndarray:
    raw = frame[label].astype("string").str.lower().str.strip()
    if raw.isna().any():
        raise ValueError(f"Missing values found in {label}")
    observed = set(raw.unique())
    if positive not in observed:
        raise ValueError(f"Positive class {positive!r} absent from {label}: {sorted(observed)}")
    return (raw == positive).to_numpy(dtype=int)


def folds_frame(
    cv: StratifiedGroupKFold,
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    row_ids: np.ndarray,
) -> tuple[list[tuple[np.ndarray, np.ndarray]], pd.DataFrame]:
    splits = list(cv.split(X, y, groups))
    fold = np.full(len(y), -1, dtype=int)
    for fold_id, (_, validation_index) in enumerate(splits, start=1):
        fold[validation_index] = fold_id
    if (fold < 0).any():
        raise RuntimeError("Some training observations were not assigned to a fold")
    details = pd.DataFrame({"row_id": row_ids, "lineage": groups, "y": y, "fold": fold})
    return splits, details


def threshold_youden(y: np.ndarray, probability: np.ndarray) -> float:
    false_positive_rate, true_positive_rate, thresholds = roc_curve(y, probability)
    finite = np.isfinite(thresholds)
    score = np.where(finite, true_positive_rate - false_positive_rate, -np.inf)
    best = np.flatnonzero(np.isclose(score, np.nanmax(score)))
    if len(best) == 0:
        return 0.5
    candidates = thresholds[best]
    return float(candidates[np.argmin(np.abs(candidates - 0.5))])


def calibration_stats(y: np.ndarray, probability: np.ndarray) -> tuple[float, float]:
    p = np.clip(np.asarray(probability, dtype=float), 1e-6, 1 - 1e-6)
    x = logit(p).reshape(-1, 1)
    try:
        model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000).fit(x, y)
        return float(model.intercept_[0]), float(model.coef_[0, 0])
    except Exception:
        return np.nan, np.nan


def metric_bundle(
    y: np.ndarray,
    probability: np.ndarray,
    threshold: float,
    include_calibration: bool = True,
) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    p = np.clip(np.asarray(probability, dtype=float), 1e-12, 1 - 1e-12)
    predicted = (p >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, predicted, labels=[0, 1]).ravel()
    class_precision, class_recall, class_f1, _ = precision_recall_fscore_support(
        y,
        predicted,
        labels=[0, 1],
        average=None,
        zero_division=0,
    )
    macro_precision, macro_recall, macro_f1, _ = precision_recall_fscore_support(
        y,
        predicted,
        labels=[0, 1],
        average="macro",
        zero_division=0,
    )
    weighted_precision, weighted_recall, weighted_f1, _ = precision_recall_fscore_support(
        y,
        predicted,
        labels=[0, 1],
        average="weighted",
        zero_division=0,
    )

    def safe_ratio(numerator: float, denominator: float) -> float:
        return float(numerator / denominator) if denominator else np.nan

    values = {
        "n": int(len(y)),
        "positive_n": int(y.sum()),
        "prevalence": float(y.mean()),
        "threshold": float(threshold),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "brier_score": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
        "accuracy": float(accuracy_score(y, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
        "sensitivity": safe_ratio(tp, tp + fn),
        "specificity": safe_ratio(tn, tn + fp),
        "false_negative_rate": safe_ratio(fn, fn + tp),
        "false_positive_rate": safe_ratio(fp, fp + tn),
        "ppv": safe_ratio(tp, tp + fp),
        "npv": safe_ratio(tn, tn + fn),
        "f1": float(f1_score(y, predicted, zero_division=0)),
        "class_0_precision": float(class_precision[0]),
        "class_0_recall": float(class_recall[0]),
        "class_0_f1": float(class_f1[0]),
        "class_1_precision": float(class_precision[1]),
        "class_1_recall": float(class_recall[1]),
        "class_1_f1": float(class_f1[1]),
        "macro_precision": float(macro_precision),
        "macro_recall": float(macro_recall),
        "macro_f1": float(macro_f1),
        "weighted_precision": float(weighted_precision),
        "weighted_recall": float(weighted_recall),
        "weighted_f1": float(weighted_f1),
        "mcc": float(matthews_corrcoef(y, predicted)),
    }
    if include_calibration:
        calibration_intercept, calibration_slope = calibration_stats(y, p)
        values["calibration_intercept"] = calibration_intercept
        values["calibration_slope"] = calibration_slope
    return values


def cluster_bootstrap_ci(
    y: np.ndarray,
    probability: np.ndarray,
    groups: np.ndarray,
    threshold: float,
    iterations: int,
    seed: int,
) -> dict[str, tuple[float, float]]:
    if iterations <= 0:
        return {}
    y = np.asarray(y)
    probability = np.asarray(probability)
    groups = np.asarray(groups).astype(str)
    unique_groups = np.unique(groups)
    indices = {group: np.flatnonzero(groups == group) for group in unique_groups}
    rng = np.random.default_rng(seed)
    samples: dict[str, list[float]] = {}
    attempts = 0
    valid_iterations = 0
    max_attempts = max(iterations * 10, iterations)
    while valid_iterations < iterations and attempts < max_attempts:
        attempts += 1
        selected = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        sampled_index = np.concatenate([indices[group] for group in selected])
        sampled_y = y[sampled_index]
        if np.unique(sampled_y).size < 2:
            continue
        valid_iterations += 1
        # Calibration intercept/slope retain full-sample point estimates. Their
        # repeated logistic refits dominate runtime and are not part of the
        # requested discrimination/classification CI set.
        values = metric_bundle(
            sampled_y,
            probability[sampled_index],
            threshold,
            include_calibration=False,
        )
        for metric, value in values.items():
            if metric in {"n", "positive_n", "prevalence", "threshold", "tn", "fp", "fn", "tp"}:
                continue
            samples.setdefault(metric, []).append(value)
    if valid_iterations < iterations:
        warnings.warn(
            f"Only {valid_iterations} of {iterations} requested bootstrap replicates contained both classes",
            RuntimeWarning,
        )
    intervals = {}
    for metric, values in samples.items():
        finite_values = np.asarray(values, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        if finite_values.size == 0:
            intervals[metric] = (np.nan, np.nan)
        else:
            intervals[metric] = (
                float(np.percentile(finite_values, 2.5)),
                float(np.percentile(finite_values, 97.5)),
            )
    return intervals


def fit_one_model(
    name: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    splits: list[tuple[np.ndarray, np.ndarray]],
    seed: int,
    smoke: bool,
    model_jobs: int,
) -> tuple[np.ndarray, np.ndarray, BaseEstimator, dict[str, Any], pd.DataFrame]:
    spec = model_spec(name, seed, float(y_train.mean()), smoke, model_jobs)
    start = time.time()
    log(f"    {name}: tuning started ({spec.search_kind})")
    if spec.search_kind == "none":
        best_params: dict[str, Any] = {}
        best_estimator = clone(spec.estimator)
        search_results = pd.DataFrame(
            [{"model": name, "params": "{}", "mean_test_score": np.nan, "std_test_score": np.nan, "rank_test_score": 1}]
        )
        cv_score = np.nan
    else:
        common = dict(
            estimator=spec.estimator,
            scoring="roc_auc",
            cv=splits,
            refit=True,
            n_jobs=spec.n_jobs,
            return_train_score=False,
            # Some GAM configurations can legitimately diverge when smoothing
            # is too weak. Record those candidates as failed and select only
            # among stable configurations; all other estimator errors likewise
            # remain visible as NaN rows in the exported CV search table.
            error_score=np.nan,
            verbose=0,
        )
        if spec.search_kind == "grid":
            search = GridSearchCV(param_grid=spec.search_space, **common)
        elif spec.search_kind == "random":
            search = RandomizedSearchCV(
                param_distributions=spec.search_space,
                n_iter=spec.n_iter,
                random_state=seed,
                **common,
            )
        else:
            raise ValueError(spec.search_kind)
        search.fit(X_train, y_train)
        best_params = search.best_params_
        best_estimator = search.best_estimator_
        cv_score = float(search.best_score_)
        columns = [
            column
            for column in search.cv_results_
            if column == "params"
            or column.startswith("param_")
            or column in {"mean_test_score", "std_test_score", "rank_test_score", "mean_fit_time"}
        ]
        search_results = pd.DataFrame(search.cv_results_)[columns]
        search_results.insert(0, "model", name)
        search_results["params"] = search_results["params"].map(lambda value: json.dumps(jsonable(value), sort_keys=True))

    oof_probability = cross_val_predict(
        clone(best_estimator),
        X_train,
        y_train,
        cv=splits,
        method="predict_proba",
        n_jobs=spec.n_jobs,
    )[:, 1]
    fitted = clone(best_estimator).fit(X_train, y_train)
    test_probability = fitted.predict_proba(X_test)[:, 1]
    elapsed = time.time() - start
    actual_auc = roc_auc_score(y_train, oof_probability)
    log(f"    {name}: completed in {elapsed / 60:.1f} min; selected CV AUROC={cv_score:.3f}, OOF AUROC={actual_auc:.3f}")
    metadata = {
        "model": name,
        "best_params": jsonable(best_params),
        "selected_cv_roc_auc": cv_score,
        "oof_roc_auc": float(actual_auc),
        "elapsed_seconds": elapsed,
    }
    return oof_probability, test_probability, fitted, metadata, search_results


def fit_convex_super_learner(y: np.ndarray, predictions: np.ndarray) -> np.ndarray:
    n_models = predictions.shape[1]
    initial = np.repeat(1.0 / n_models, n_models)

    def objective(weights: np.ndarray) -> float:
        combined = np.clip(predictions @ weights, 1e-8, 1 - 1e-8)
        return log_loss(y, combined, labels=[0, 1])

    result = minimize(
        objective,
        initial,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * n_models,
        constraints={"type": "eq", "fun": lambda weights: weights.sum() - 1.0},
        options={"maxiter": 2000, "ftol": 1e-12},
    )
    if not result.success:
        raise RuntimeError(f"Super Learner weight optimization failed: {result.message}")
    weights = np.clip(result.x, 0.0, 1.0)
    return weights / weights.sum()


def super_learner_predictions(
    y_train: np.ndarray,
    base_oof: np.ndarray,
    base_test: np.ndarray,
    splits: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    oof = np.full(len(y_train), np.nan)
    for training_index, validation_index in splits:
        weights = fit_convex_super_learner(y_train[training_index], base_oof[training_index])
        oof[validation_index] = base_oof[validation_index] @ weights
    weights = fit_convex_super_learner(y_train, base_oof)
    return oof, base_test @ weights, weights


def metric_rows(
    outcome: str,
    window: int,
    model: str,
    split: str,
    values: dict[str, float],
    intervals: dict[str, tuple[float, float]],
) -> list[dict[str, Any]]:
    rows = []
    for metric, estimate in values.items():
        lower, upper = intervals.get(metric, (np.nan, np.nan))
        rows.append(
            {
                "outcome": outcome,
                "window_days": window,
                "model": model,
                "split": split,
                "metric": metric,
                "estimate": estimate,
                "ci_95_lower": lower,
                "ci_95_upper": upper,
            }
        )
    return rows


def class_distribution_rows(
    outcome: str,
    window: int,
    split: str,
    y: np.ndarray,
) -> list[dict[str, Any]]:
    y = np.asarray(y, dtype=int)
    counts = np.bincount(y, minlength=2)
    majority_class = int(np.argmax(counts))
    majority_accuracy = float(counts[majority_class] / len(y))
    class_names = OUTCOMES[outcome]["class_names"]
    return [
        {
            "outcome": outcome,
            "window_days": window,
            "split": split,
            "class_id": class_id,
            "class_label": class_names[class_id],
            "n": int(counts[class_id]),
            "total_n": int(len(y)),
            "proportion": float(counts[class_id] / len(y)),
            "is_majority_class": class_id == majority_class,
            "majority_class_id": majority_class,
            "majority_class_label": class_names[majority_class],
            "majority_baseline_accuracy": majority_accuracy,
        }
        for class_id in (0, 1)
    ]


def per_class_metric_rows(
    outcome: str,
    window: int,
    model: str,
    split: str,
    y: np.ndarray,
    values: dict[str, float],
    intervals: dict[str, tuple[float, float]],
) -> list[dict[str, Any]]:
    y = np.asarray(y, dtype=int)
    class_names = OUTCOMES[outcome]["class_names"]
    rows = []
    for class_id in (0, 1):
        row: dict[str, Any] = {
            "outcome": outcome,
            "window_days": window,
            "model": model,
            "split": split,
            "class_id": class_id,
            "class_label": class_names[class_id],
            "support": int((y == class_id).sum()),
        }
        for metric in ("precision", "recall", "f1"):
            key = f"class_{class_id}_{metric}"
            lower, upper = intervals.get(key, (np.nan, np.nan))
            row[metric] = values[key]
            row[f"{metric}_ci_95_lower"] = lower
            row[f"{metric}_ci_95_upper"] = upper
        rows.append(row)
    return rows


def confusion_matrix_row(
    outcome: str,
    window: int,
    model: str,
    split: str,
    values: dict[str, float],
) -> dict[str, Any]:
    return {
        "outcome": outcome,
        "window_days": window,
        "model": model,
        "split": split,
        "negative_class": OUTCOMES[outcome]["class_names"][0],
        "positive_class": OUTCOMES[outcome]["class_names"][1],
        "tn": int(values["tn"]),
        "fp": int(values["fp"]),
        "fn": int(values["fn"]),
        "tp": int(values["tp"]),
        "false_negative_rate": values["false_negative_rate"],
        "false_positive_rate": values["false_positive_rate"],
    }


def plot_task_curves(
    task_dir: Path,
    outcome: str,
    window: int,
    y_test: np.ndarray,
    probabilities: dict[str, np.ndarray],
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(17, 5))
    for model, probability in probabilities.items():
        fpr, tpr, _ = roc_curve(y_test, probability)
        axes[0].plot(fpr, tpr, label=f"{model} ({roc_auc_score(y_test, probability):.3f})")
        precision, recall, _ = precision_recall_curve(y_test, probability)
        axes[1].plot(recall, precision, label=f"{model} ({average_precision_score(y_test, probability):.3f})")
        observed, predicted = calibration_curve(y_test, probability, n_bins=8, strategy="quantile")
        axes[2].plot(predicted, observed, marker="o", label=model)
    axes[0].plot([0, 1], [0, 1], "k--", linewidth=1)
    axes[0].set(xlabel="False-positive rate", ylabel="True-positive rate", title="ROC curves")
    axes[1].axhline(y_test.mean(), color="black", linestyle="--", linewidth=1)
    axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision–recall curves")
    axes[2].plot([0, 1], [0, 1], "k--", linewidth=1)
    axes[2].set(xlabel="Mean predicted probability", ylabel="Observed proportion", title="Calibration")
    for axis in axes:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.grid(alpha=0.2)
    axes[0].legend(fontsize=7, loc="lower right")
    axes[1].legend(fontsize=7, loc="best")
    axes[2].legend(fontsize=7, loc="best")
    figure.suptitle(f"{outcome.title()} outcome — {window}-day inputs (locked test set)")
    figure.tight_layout()
    figure.savefig(task_dir / "curves.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def run_task(
    outcome: str,
    window: int,
    train: pd.DataFrame,
    test: pd.DataFrame,
    output_root: Path,
    bootstrap_iterations: int,
    smoke: bool,
    selected_models: tuple[str, ...],
    model_jobs: int = 1,
) -> None:
    task_name = f"{outcome}_{window}d"
    task_dir = output_root / task_name
    task_dir.mkdir(parents=True, exist_ok=True)
    complete_marker = task_dir / "done.json"
    if complete_marker.exists() and not smoke:
        expected_models = list(selected_models) + (["super_learner"] if len(selected_models) >= 2 else [])
        required_outputs = {
            "metrics.csv",
            "predictions.csv",
            "cv.csv",
            "class_counts.csv",
            "per_class.csv",
            "confusion.csv",
        }
        try:
            marker = json.loads(complete_marker.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            marker = {}
        marker_matches = (
            marker.get("models") == expected_models
            and marker.get("bootstrap_iterations") == bootstrap_iterations
            and all((task_dir / filename).exists() for filename in required_outputs)
        )
        if marker_matches:
            log(f"{task_name}: already complete with matching outputs; skipping (use --overwrite to rerun)")
            return
        log(f"{task_name}: prior completion is incompatible with the requested models/outputs; rerunning")

    spec = OUTCOMES[outcome]
    validate_data(train, test, spec["label"])
    X_train = train.loc[:, FEATURES].to_numpy(dtype=float)
    X_test = test.loc[:, FEATURES].to_numpy(dtype=float)
    y_train = make_binary(train, spec["label"], spec["positive"])
    y_test = make_binary(test, spec["label"], spec["positive"])
    train_groups = train["id_2"].astype(str).to_numpy()
    test_groups = test["id_2"].astype(str).to_numpy()
    train_row_ids = np.arange(len(train))
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=stable_seed(task_name, "folds"))
    splits, fold_details = folds_frame(cv, X_train, y_train, train_groups, train_row_ids)
    fold_details.to_csv(task_dir / "folds.csv", index=False)
    distribution_rows = class_distribution_rows(outcome, window, "train", y_train)
    distribution_rows.extend(class_distribution_rows(outcome, window, "test", y_test))

    log(
        f"{task_name}: started; train n={len(y_train)} ({y_train.mean():.1%} positive), "
        f"test n={len(y_test)} ({y_test.mean():.1%} positive)"
    )
    oof_probabilities: dict[str, np.ndarray] = {}
    test_probabilities: dict[str, np.ndarray] = {}
    metadata: dict[str, Any] = {}
    search_frames: list[pd.DataFrame] = []

    for model_name in selected_models:
        model_seed = stable_seed(task_name, model_name)
        oof, test_probability, fitted, details, search_results = fit_one_model(
            model_name,
            X_train,
            y_train,
            X_test,
            splits,
            model_seed,
            smoke,
            model_jobs,
        )
        oof_probabilities[model_name] = oof
        test_probabilities[model_name] = test_probability
        metadata[model_name] = details
        search_frames.append(search_results)
        joblib.dump(fitted, task_dir / f"{model_name}.joblib", compress=3)

    if len(selected_models) >= 2:
        log("    super_learner: fitting nonnegative log-loss-optimal stacking weights")
        base_oof = np.column_stack([oof_probabilities[name] for name in selected_models])
        base_test = np.column_stack([test_probabilities[name] for name in selected_models])
        sl_oof, sl_test, weights = super_learner_predictions(y_train, base_oof, base_test, splits)
        oof_probabilities["super_learner"] = sl_oof
        test_probabilities["super_learner"] = sl_test
        metadata["super_learner"] = {
            "model": "super_learner",
            "base_models": list(selected_models),
            "weights": {name: float(weight) for name, weight in zip(selected_models, weights)},
            "oof_roc_auc": float(roc_auc_score(y_train, sl_oof)),
        }
        log(
            "    super_learner: completed; OOF AUROC="
            f"{roc_auc_score(y_train, sl_oof):.3f}; weights="
            + ", ".join(f"{name}={weight:.2f}" for name, weight in zip(selected_models, weights))
        )

    prediction_rows = []
    performance_rows = []
    per_class_rows = []
    confusion_rows = []
    for model_name, oof_probability in oof_probabilities.items():
        threshold = threshold_youden(y_train, oof_probability)
        test_probability = test_probabilities[model_name]
        oof_metrics = metric_bundle(y_train, oof_probability, threshold)
        test_metrics = metric_bundle(y_test, test_probability, threshold)
        test_intervals = cluster_bootstrap_ci(
            y_test,
            test_probability,
            test_groups,
            threshold,
            bootstrap_iterations,
            stable_seed(task_name, model_name, "bootstrap"),
        )
        performance_rows.extend(metric_rows(outcome, window, model_name, "cv_oof", oof_metrics, {}))
        performance_rows.extend(metric_rows(outcome, window, model_name, "test", test_metrics, test_intervals))
        for split_name, y, values, intervals in (
            ("cv_oof", y_train, oof_metrics, {}),
            ("test", y_test, test_metrics, test_intervals),
        ):
            per_class_rows.extend(
                per_class_metric_rows(outcome, window, model_name, split_name, y, values, intervals)
            )
            confusion_rows.append(confusion_matrix_row(outcome, window, model_name, split_name, values))
        for split_name, frame, y, groups, probability in (
            ("cv_oof", train, y_train, train_groups, oof_probability),
            ("test", test, y_test, test_groups, test_probability),
        ):
            for index in range(len(y)):
                prediction_rows.append(
                    {
                        "outcome": outcome,
                        "window_days": window,
                        "model": model_name,
                        "split": split_name,
                        "row_index": index,
                        "country": str(frame.iloc[index]["id_1"]),
                        "lineage": groups[index],
                        "observed": int(y[index]),
                        "predicted_probability": float(probability[index]),
                        "threshold": threshold,
                        "predicted_class": int(probability[index] >= threshold),
                    }
                )

    pd.DataFrame(prediction_rows).to_csv(task_dir / "predictions.csv", index=False)
    pd.DataFrame(performance_rows).to_csv(task_dir / "metrics.csv", index=False)
    pd.DataFrame(distribution_rows).to_csv(task_dir / "class_counts.csv", index=False)
    pd.DataFrame(per_class_rows).to_csv(task_dir / "per_class.csv", index=False)
    pd.DataFrame(confusion_rows).to_csv(task_dir / "confusion.csv", index=False)
    pd.concat(search_frames, ignore_index=True).to_csv(task_dir / "cv.csv", index=False)
    with (task_dir / "selection.json").open("w", encoding="utf-8") as handle:
        json.dump(jsonable(metadata), handle, indent=2, sort_keys=True)
    with PLOT_LOCK:
        plot_task_curves(task_dir, outcome, window, y_test, test_probabilities)
    marker = {
        "task": task_name,
        "completed_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "models": list(oof_probabilities),
        "bootstrap_iterations": bootstrap_iterations,
    }
    with complete_marker.open("w", encoding="utf-8") as handle:
        json.dump(marker, handle, indent=2)
    log(f"{task_name}: COMPLETE; task outputs saved to {task_dir}")


def aggregate_outputs(output_root: Path) -> None:
    performance_files = sorted(output_root.glob("*_*d/metrics.csv"))
    prediction_files = sorted(output_root.glob("*_*d/predictions.csv"))
    search_files = sorted(output_root.glob("*_*d/cv.csv"))
    distribution_files = sorted(output_root.glob("*_*d/class_counts.csv"))
    per_class_files = sorted(output_root.glob("*_*d/per_class.csv"))
    confusion_files = sorted(output_root.glob("*_*d/confusion.csv"))
    if not performance_files:
        return
    companion_sets = {
        "predictions": prediction_files,
        "CV searches": search_files,
        "class distributions": distribution_files,
        "per-class metrics": per_class_files,
        "confusion matrices": confusion_files,
    }
    incomplete = [name for name, files in companion_sets.items() if len(files) != len(performance_files)]
    if incomplete:
        raise RuntimeError(
            "Aggregate output is incomplete for: "
            + ", ".join(incomplete)
            + ". Rerun the affected tasks with --overwrite or use a clean output directory."
        )
    performance = pd.concat((pd.read_csv(path) for path in performance_files), ignore_index=True)
    predictions = pd.concat((pd.read_csv(path) for path in prediction_files), ignore_index=True)
    searches = pd.concat((pd.read_csv(path) for path in search_files), ignore_index=True)
    distributions = pd.concat((pd.read_csv(path) for path in distribution_files), ignore_index=True)
    per_class = pd.concat((pd.read_csv(path) for path in per_class_files), ignore_index=True)
    confusions = pd.concat((pd.read_csv(path) for path in confusion_files), ignore_index=True)
    # Calibration intercept and slope are reported as point estimates. Keep CI
    # reporting uniform even when older resumable task files contain bootstrap
    # calibration intervals from a prior run of the script.
    calibration = performance["metric"].isin(["calibration_intercept", "calibration_slope"])
    performance.loc[calibration, ["ci_95_lower", "ci_95_upper"]] = np.nan
    performance.to_csv(output_root / "metrics.csv", index=False)
    predictions.to_csv(output_root / "predictions.csv", index=False)
    searches.to_csv(output_root / "cv.csv", index=False)
    distributions[distributions["split"] == "test"].to_csv(output_root / "class_counts.csv", index=False)
    per_class[per_class["split"] == "test"].to_csv(output_root / "per_class.csv", index=False)
    confusions[confusions["split"] == "test"].to_csv(output_root / "confusion.csv", index=False)

    point = performance.pivot_table(
        index=["outcome", "window_days", "model", "split"],
        columns="metric",
        values="estimate",
        aggfunc="first",
    ).reset_index()
    point.to_csv(output_root / "metrics_wide.csv", index=False)

    test = performance[performance["split"] == "test"].copy()
    test["estimate_ci95"] = test.apply(
        lambda row: (
            f"{row['estimate']:.3f} ({row['ci_95_lower']:.3f}, {row['ci_95_upper']:.3f})"
            if pd.notna(row["ci_95_lower"])
            else f"{row['estimate']:.3f}"
        ),
        axis=1,
    )
    standard_metrics = [
        "roc_auc",
        "average_precision",
        "brier_score",
        "accuracy",
        "balanced_accuracy",
        "macro_precision",
        "macro_recall",
        "macro_f1",
        "weighted_precision",
        "weighted_recall",
        "weighted_f1",
        "class_0_precision",
        "class_0_recall",
        "class_0_f1",
        "class_1_precision",
        "class_1_recall",
        "class_1_f1",
        "sensitivity",
        "specificity",
        "false_negative_rate",
        "false_positive_rate",
        "ppv",
        "npv",
        "f1",
        "mcc",
        "calibration_intercept",
        "calibration_slope",
    ]
    standard = (
        test[test["metric"].isin(standard_metrics)]
        .pivot_table(
            index=["outcome", "window_days", "model"],
            columns="metric",
            values="estimate_ci95",
            aggfunc="first",
        )
        .reset_index()
    )
    standard.to_csv(output_root / "comparison.csv", index=False)

    with pd.ExcelWriter(output_root / "results.xlsx", engine="openpyxl") as writer:
        standard.to_excel(writer, sheet_name="Comparison", index=False)
        distributions[distributions["split"] == "test"].to_excel(
            writer, sheet_name="Test class distribution", index=False
        )
        per_class[per_class["split"] == "test"].to_excel(writer, sheet_name="Per-class test metrics", index=False)
        confusions[confusions["split"] == "test"].to_excel(writer, sheet_name="Test confusion matrices", index=False)
        point.to_excel(writer, sheet_name="Point estimates", index=False)
        performance.to_excel(writer, sheet_name="All metrics long", index=False)
        predictions.to_excel(writer, sheet_name="Predictions", index=False)
        searches.to_excel(writer, sheet_name="CV searches", index=False)

    test_auc = point[point["split"] == "test"].sort_values(
        ["outcome", "window_days", "roc_auc"], ascending=[True, True, False]
    )
    winners = test_auc.groupby(["outcome", "window_days"], as_index=False).first()
    test_distributions = distributions[distributions["split"] == "test"].copy()
    lines = [
        "# Binary model comparison",
        "",
        "## Design",
        "",
        "- Outcomes: Duration (Long=1) and Growth (Large=1).",
        "- Inputs: all 28 features, evaluated separately at 14, 21, and 28 days.",
        "- Selection: lineage-grouped stratified 5-fold CV on the training set; AUROC was the tuning objective.",
        "- Classification threshold: Youden J selected from out-of-fold training predictions only.",
        "- Final evaluation: locked temporal test set; 95% percentile CIs use lineage-cluster bootstrap. "
        "Replicates containing only one outcome class are discarded and redrawn; metrics that remain undefined "
        "are reported as NA.",
        "- Benchmarks: majority-class accuracy is reported from the locked test-set distribution, and an "
        "unpenalized binary logistic regression is evaluated head to head with all fitted models.",
        "- Multinomial logistic regression is not applicable to these prespecified binary outcomes. It would be "
        "required if the original three-level responses were retained.",
        "- Super Learner: nonnegative weights constrained to sum to one and optimized for OOF log loss.",
        "",
        "## Locked test-set class distribution",
        "",
        "| Outcome | Window | Class 0, n (%) | Class 1, n (%) | Majority class | Majority baseline accuracy |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for (outcome, window), group in test_distributions.groupby(["outcome", "window_days"], sort=True):
        negative = group[group["class_id"] == 0].iloc[0]
        positive = group[group["class_id"] == 1].iloc[0]
        lines.append(
            f"| {outcome} | {int(window)} | {negative['class_label']}: {int(negative['n'])} "
            f"({negative['proportion']:.1%}) | {positive['class_label']}: {int(positive['n'])} "
            f"({positive['proportion']:.1%}) | {positive['majority_class_label']} | "
            f"{positive['majority_baseline_accuracy']:.3f} |"
        )
    lines.extend(
        [
        "",
        "## Highest test AUROC in each task (descriptive, not used for model selection)",
        "",
        "| Outcome | Window | Model | AUROC | AUPRC | Balanced accuracy | Brier |",
        "|---|---:|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in winners.iterrows():
        lines.append(
            f"| {row['outcome']} | {int(row['window_days'])} | {row['model']} | "
            f"{row['roc_auc']:.3f} | {row['average_precision']:.3f} | "
            f"{row['balanced_accuracy']:.3f} | {row['brier_score']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Head-to-head locked test performance",
            "",
            "All entries are estimates with lineage-cluster bootstrap 95% CIs.",
            "",
            "| Outcome | Window | Model | Accuracy | Balanced accuracy | Macro-F1 | PPV | Recall | NPV | FNR |",
            "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in standard.sort_values(["outcome", "window_days", "model"]).iterrows():
        lines.append(
            f"| {row['outcome']} | {int(row['window_days'])} | {row['model']} | {row['accuracy']} | "
            f"{row['balanced_accuracy']} | {row['macro_f1']} | {row['ppv']} | {row['sensitivity']} | "
            f"{row['npv']} | {row['false_negative_rate']} |"
        )
    lines.extend(
        [
            "",
            "The full per-class precision, recall, F1, and their 95% CIs are in `per_class.csv`. "
            "TN, FP, FN, TP, FNR, and FPR are in `confusion.csv`. The complete head-to-head table "
            "is in `comparison.csv` and `results.xlsx`.",
        ]
    )
    (output_root / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_run_manifest(output_root: Path, args: argparse.Namespace) -> None:
    packages = [
        "interpret-core",
        "joblib",
        "matplotlib",
        "numpy",
        "openpyxl",
        "pandas",
        "pygam",
        "rdata",
        "scikit-learn",
        "scipy",
        "xgboost",
    ]
    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "python": sys.version,
        "platform": platform.platform(),
        "seed": SEED,
        "arguments": vars(args),
        "outcomes": OUTCOMES,
        "windows": list(WINDOWS),
        "features": list(FEATURES),
        "models": list(ALL_MODELS),
        "bootstrap": {
            "unit": "id_2 lineage cluster",
            "method": "percentile",
            "confidence_level": 0.95,
            "requested_iterations": args.bootstrap_iterations,
            "effective_iterations_per_task": 0 if args.smoke else args.bootstrap_iterations,
            "single_class_replicates": "discarded and redrawn, up to 10 times the requested iterations",
            "undefined_metrics": "reported as NA",
        },
        "packages": {package: importlib.metadata.version(package) for package in packages},
    }
    with (output_root / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(jsonable(manifest), handle, indent=2, sort_keys=True)


def parse_args(
    fixed_outcome: str | None = None,
    default_output_dir: Path | None = None,
) -> argparse.Namespace:
    if fixed_outcome is not None and fixed_outcome not in OUTCOMES:
        raise ValueError(f"Unknown fixed outcome: {fixed_outcome}")
    description = __doc__
    if fixed_outcome is not None:
        outcome_spec = OUTCOMES[fixed_outcome]
        description = (
            f"Tune and compare binary classifiers for the {fixed_outcome} response only. "
            f"{outcome_spec['definition']}."
        )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--train-rds", type=Path, default=Path("code/train.rds"))
    parser.add_argument("--test-rds", type=Path, default=Path("code/test.rds"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir or Path("code/modeling_results/all"),
    )
    if fixed_outcome is None:
        parser.add_argument("--outcomes", nargs="+", choices=tuple(OUTCOMES), default=list(OUTCOMES))
    else:
        parser.set_defaults(outcomes=[fixed_outcome])
    parser.add_argument("--windows", nargs="+", type=int, choices=WINDOWS, default=list(WINDOWS))
    parser.add_argument("--models", nargs="+", choices=BASE_MODELS, default=list(BASE_MODELS))
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--task-jobs", type=int, default=1, help="Number of outcome/window tasks to run concurrently")
    parser.add_argument("--model-jobs", type=int, default=1, help="Threads used by each internally parallel model")
    parser.add_argument("--smoke", action="store_true", help="Use a minimal tuning grid and no bootstrap CIs")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.task_jobs < 1 or args.model_jobs < 1:
        parser.error("--task-jobs and --model-jobs must both be positive integers")
    return args


def main(
    fixed_outcome: str | None = None,
    default_output_dir: Path | None = None,
) -> None:
    args = parse_args(fixed_outcome=fixed_outcome, default_output_dir=default_output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        for outcome in args.outcomes:
            for window in args.windows:
                marker = args.output_dir / f"{outcome}_{window}d" / "done.json"
                if marker.exists():
                    marker.unlink()
    write_run_manifest(args.output_dir, args)
    log(f"Reading {args.train_rds} and {args.test_rds}")
    train_sets = rdata.read_rds(args.train_rds)
    test_sets = rdata.read_rds(args.test_rds)
    tasks = []
    for outcome in args.outcomes:
        for window in args.windows:
            train_key = f"train_{window}"
            test_key = f"test_{window}"
            if train_key not in train_sets or test_key not in test_sets:
                raise KeyError(f"Missing RDS element: {train_key} or {test_key}")
            tasks.append(
                (
                outcome,
                window,
                train_sets[train_key],
                test_sets[test_key],
                args.output_dir,
                0 if args.smoke else args.bootstrap_iterations,
                args.smoke,
                tuple(args.models),
                args.model_jobs,
                )
            )
    task_workers = min(args.task_jobs, len(tasks))
    log(
        f"Scheduling {len(tasks)} task(s) with {task_workers} concurrent task worker(s) "
        f"and up to {args.model_jobs} model thread(s) per task"
    )
    if task_workers == 1:
        for task in tasks:
            run_task(*task)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=task_workers) as executor:
            futures = [executor.submit(run_task, *task) for task in tasks]
            for future in concurrent.futures.as_completed(futures):
                future.result()
    aggregate_outputs(args.output_dir)
    log(f"ALL REQUESTED TASKS COMPLETE; aggregate outputs saved to {args.output_dir}")


if __name__ == "__main__":
    main()
