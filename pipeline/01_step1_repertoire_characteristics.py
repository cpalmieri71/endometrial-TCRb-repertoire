#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 1, v2.0: non-duplicated repertoire characterization.

One input TSV = one biological repertoire. PBMC and biopsy are NOT technical
replicates. No overlap analysis, rarefaction, hypothesis tests between tissues,
productivity filtering, or biological QC exclusions are performed here.
Input contract: prefiltered clonotype table; default clone=aaSeqCDR3,
count=readCount. Equal amino-acid CDR3 strings are collapsed within each file.

Canonical defaults:
  sampling: depth, observed_richness, singleton_count, doubleton_count,
            goods_coverage, chao1_richness;
  architecture: hill_D1, hill_D2, gini_frequency, top10_frequency,
                top100_frequency, fraction_clones_to_50pct;
  optional descriptive PL / truncated-PL / lognormal tail fits (on by default).

Paper table: depth, richness, coverage, Hill D1, Hill D2, top-10 frequency.
No separate Shannon / Simpson / D0 / evenness / clonality aliases are exported
in canonical tables. A compatibility export is opt-in and physically isolated
in legacy_compat/, never added to the canonical analysis tables.

Usage:
  python3 1-build_repertoire_characteristics_v2.py \
      --data_dir repSeq --out_dir results/01_repertoire_characterization_v2
  # Without the optional powerlaw dependency:
  python3 1-build_repertoire_characteristics_v2.py \
      --data_dir repSeq --out_dir results/01_repertoire_characterization_v2 \
      --skip_tail_fitting

Dependencies: numpy, pandas; powerlaw for tail fitting. Python >=3.9.
See README_step1_v2.md and Methods_repertoire_characterization_v2.md.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import re
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

SCRIPT_VERSION = "2.0.0"
SCHEMA_VERSION = "endometrium_step1_nonredundant_v2"
DEFAULT_PATTERN = r"^(?P<subject>\d+)_(?P<compartment>PBMC|biopsia|biopsy)\.tsv$"
COMPARTMENT_MAP = {"pbmc": "PBMC", "biopsia": "BIOPSY", "biopsy": "BIOPSY"}
DEFAULT_HILL_ORDERS = (1.0, 2.0)
DEFAULT_TOP_N = (10, 100)
DEFAULT_CUMULATIVE_TARGETS = (0.50,)
IDENTITY_COLUMNS = ["subject", "compartment", "file", "pair_complete"]

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path

def safe_div(a: float, b: float) -> float:
    if not np.isfinite(b) or b == 0:
        return np.nan
    return float(a / b)

def parse_float_list(value: str) -> List[float]:
    return [float(x.strip()) for x in str(value).split(",") if x.strip()]

def parse_int_list(value: str) -> List[int]:
    return [int(x.strip()) for x in str(value).split(",") if x.strip()]

def natural_subject_key(value: str):
    s = str(value)
    return (0, int(s)) if s.isdigit() else (1, s)

def sanitize_q(q: float) -> str:
    if abs(q - round(q)) < 1e-12:
        return str(int(round(q)))
    return str(q).replace(".", "p")

def format_number(value: float, digits: int = 3, integer: bool = False) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    if integer:
        return f"{int(round(value)):,}"
    return f"{float(value):.{digits}f}"

def format_median_iqr(
    median: float,
    q1: float,
    q3: float,
    digits: int = 3,
    integer: bool = False,
) -> str:
    if not (np.isfinite(median) and np.isfinite(q1) and np.isfinite(q3)):
        return "NA"
    return (
        f"{format_number(median, digits, integer)} "
        f"[{format_number(q1, digits, integer)}–"
        f"{format_number(q3, digits, integer)}]"
    )

def index_repertoire_files(
    data_dir: Path,
    pattern_str: str,
) -> pd.DataFrame:
    pattern = re.compile(pattern_str, flags=re.IGNORECASE)
    rows: List[Dict] = []

    for fp in sorted(data_dir.iterdir()):
        if not fp.is_file():
            continue
        match = pattern.match(fp.name)
        if not match:
            continue

        subject = str(match.group("subject"))
        raw_compartment = str(match.group("compartment")).lower()
        if raw_compartment not in COMPARTMENT_MAP:
            raise ValueError(
                f"Unrecognized compartment '{raw_compartment}' in {fp.name}"
            )

        rows.append(
            {
                "subject": subject,
                "compartment": COMPARTMENT_MAP[raw_compartment],
                "file": fp.name,
                "file_path": str(fp.resolve()),
            }
        )

    manifest = pd.DataFrame(rows)
    if manifest.empty:
        raise ValueError(
            "No input repertoire files matched the filename pattern.\n"
            f"Directory: {data_dir}\nPattern: {pattern_str}"
        )

    duplicated = manifest.duplicated(["subject", "compartment"], keep=False)
    if duplicated.any():
        bad = manifest.loc[duplicated, ["subject", "compartment", "file"]]
        raise ValueError(
            "More than one file was found for the same subject/compartment:\n"
            + bad.to_string(index=False)
        )

    completeness = (
        manifest.groupby("subject", observed=True)["compartment"]
        .agg(lambda x: set(x))
        .to_dict()
    )
    manifest["pair_complete"] = manifest["subject"].map(
        lambda s: {"PBMC", "BIOPSY"}.issubset(completeness.get(s, set()))
    )

    subject_order = {
        s: i for i, s in enumerate(sorted(manifest["subject"].unique(), key=natural_subject_key))
    }
    compartment_order = {"PBMC": 0, "BIOPSY": 1}
    manifest["_subject_order"] = manifest["subject"].map(subject_order)
    manifest["_compartment_order"] = manifest["compartment"].map(compartment_order)
    manifest = (
        manifest.sort_values(["_subject_order", "_compartment_order", "file"])
        .drop(columns=["_subject_order", "_compartment_order"])
        .reset_index(drop=True)
    )
    return manifest

def read_repertoire(
    fp: Path,
    sep: str,
    clone_col: str,
    count_col: str,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    raw = pd.read_csv(fp, sep=sep, dtype={clone_col: "string"})
    missing = {clone_col, count_col} - set(raw.columns)
    if missing:
        raise ValueError(f"{fp.name}: missing required columns {sorted(missing)}")

    n_raw_rows = int(raw.shape[0])

    clone = raw[clone_col]
    count = pd.to_numeric(raw[count_col], errors="coerce")

    clone_valid = clone.notna() & clone.astype(str).str.strip().ne("")
    count_valid = np.isfinite(count) & (count > 0)
    valid = clone_valid & count_valid

    n_invalid_rows = int((~valid).sum())
    n_valid_rows_before_collapse = int(valid.sum())

    data = pd.DataFrame(
        {
            "clonotype": clone.loc[valid].astype(str).str.strip(),
            "count": count.loc[valid].astype(float),
        }
    )

    collapsed = (
        data.groupby("clonotype", as_index=False, observed=True)["count"]
        .sum()
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )

    stats = {
        "n_raw_rows": n_raw_rows,
        "n_valid_rows_before_collapse": n_valid_rows_before_collapse,
        "n_invalid_rows": n_invalid_rows,
        "n_duplicate_rows_collapsed": int(
            max(0, n_valid_rows_before_collapse - collapsed.shape[0])
        ),
    }
    return collapsed, stats

def shannon_entropy(freqs: Sequence[float]) -> float:
    p = np.asarray(freqs, dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    if p.size == 0:
        return np.nan
    return float(-np.sum(p * np.log(p)))

def hill_number(freqs: Sequence[float], q: float) -> float:
    p = np.asarray(freqs, dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    if p.size == 0:
        return np.nan

    q = float(q)
    if abs(q) < 1e-12:
        return float(p.size)
    if abs(q - 1.0) < 1e-12:
        return float(np.exp(shannon_entropy(p)))

    moment = np.sum(np.power(p, q))
    return float(moment ** (1.0 / (1.0 - q))) if moment > 0 else np.nan

def gini_coefficient(values: Sequence[float]) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (x >= 0)]
    if x.size == 0 or np.all(x == 0):
        return np.nan
    x = np.sort(x)
    n = x.size
    return float(
        (2.0 * np.sum(np.arange(1, n + 1) * x) / (n * np.sum(x)))
        - (n + 1.0) / n
    )

def chao1_from_counts(counts: Sequence[float]) -> Dict[str, float]:
    """
    Chao1 on observed sequencing counts.

    This is retained as a descriptive richness estimator. In RNA-based RepSeq,
    singleton/doubleton read counts should not be interpreted as literal
    one-cell/two-cell clone sizes.
    """
    c = np.asarray(counts, dtype=float)
    c = c[np.isfinite(c) & (c > 0)]
    if c.size == 0:
        return {
            "chao1_richness": np.nan,
            "singleton_count": 0,
            "doubleton_count": 0,
        }

    # readCount is expected to represent integer counts
    integer_like = np.allclose(c, np.rint(c), rtol=0, atol=1e-8)
    if not integer_like:
        return {
            "chao1_richness": np.nan,
            "singleton_count": np.nan,
            "doubleton_count": np.nan,
        }

    ci = np.rint(c).astype(np.int64)
    S = int(ci.size)
    f1 = int(np.sum(ci == 1))
    f2 = int(np.sum(ci == 2))

    # Classic Chao1 form, matching the previous pipeline.
    if f2 > 0:
        chao = S + (f1 * f1) / (2.0 * f2)
    else:
        chao = S + f1 * (f1 - 1) / 2.0

    return {
        "chao1_richness": float(chao),
        "singleton_count": f1,
        "doubleton_count": f2,
    }

def clones_needed_for_fraction(
    sorted_desc_freqs: np.ndarray,
    target: float,
) -> float:
    if sorted_desc_freqs.size == 0:
        return np.nan
    cumulative = np.cumsum(sorted_desc_freqs)
    idx = int(np.searchsorted(cumulative, target, side="left"))
    return float(min(idx + 1, sorted_desc_freqs.size))

def target_label(target: float) -> str:
    """Do not round distinct targets to the same output-column name."""
    return format(float(target) * 100.0, ".12g").replace(".", "p")


def compute_repertoire_metrics(
    counts: Sequence[float],
    hill_orders: Sequence[float] = DEFAULT_HILL_ORDERS,
    top_n_values: Sequence[int] = DEFAULT_TOP_N,
    cumulative_targets: Sequence[float] = DEFAULT_CUMULATIVE_TARGETS,
) -> Dict[str, object]:
    """Canonical metrics only; all formulas match v1 for retained quantities."""
    c = np.asarray(counts, dtype=float)
    c = c[np.isfinite(c) & (c > 0)]
    if not c.size or not np.isfinite(c.sum()):
        raise ValueError("A repertoire must contain positive finite total abundance.")
    total = float(c.sum())
    integer_like = bool(abs(total - round(total)) <= 1e-8
                        and np.allclose(c, np.rint(c), rtol=0, atol=1e-8))
    p = c / total
    p_desc = np.sort(p)[::-1]
    S = int(c.size)
    chao = chao1_from_counts(c)
    f1 = chao["singleton_count"]
    out = {
        "depth": int(round(total)) if integer_like else total,
        "counts_integer_like": integer_like,
        "observed_richness": S,
        "singleton_count": f1,
        "doubleton_count": chao["doubleton_count"],
        "goods_coverage": 1.0 - float(f1) / total if np.isfinite(f1) else np.nan,
        "chao1_richness": chao["chao1_richness"],
    }
    # D1 and D2 are always available for the paper table. D0 is represented by S.
    for q in sorted(set([1.0, 2.0] + [float(q) for q in hill_orders])):
        if abs(q) < 1e-12:
            continue
        out[f"hill_D{sanitize_q(q)}"] = hill_number(p, q)
    out["gini_frequency"] = gini_coefficient(p)
    # top-10 is always present; additional cutoffs are configurable.
    for n in sorted(set([10] + list(top_n_values))):
        out[f"top{n}_frequency"] = float(np.sum(p_desc[:n]))
    for target in sorted(set(cumulative_targets)):
        k = clones_needed_for_fraction(p_desc, target)
        out[f"fraction_clones_to_{target_label(target)}pct"] = k / S
    return out


def legacy_metrics(counts: Sequence[float], canonical: Dict[str, object]) -> Dict[str, object]:
    """v1 default metric vocabulary, ONLY for explicitly requested legacy files.

    This does not add extra analysis variables to canonical tables. Extra Hill
    orders and rank cutoffs are evaluated only on this opt-in code path.
    """
    c = np.asarray(counts, dtype=float)
    p = c / c.sum()
    p_desc = np.sort(p)[::-1]
    S = int(c.size)
    H = shannon_entropy(p)
    concentration = float(np.sum(p * p))
    J = H / math.log(S) if S > 1 else np.nan
    f1, f2 = canonical["singleton_count"], canonical["doubleton_count"]
    chao = canonical["chao1_richness"]
    out = dict(canonical)
    out.update({
        "shannon_entropy": H,
        "simpson_concentration": concentration,
        "inverse_simpson_D2": canonical["hill_D2"],
        "gini_simpson": 1.0 - concentration,
        "pielou_evenness": J,
        "shannon_clonality": 1.0 - J if np.isfinite(J) else np.nan,
        "berger_parker_dominance": float(p_desc[0]),
        "singleton_fraction_clonotypes": safe_div(float(f1), S),
        "doubleton_fraction_clonotypes": safe_div(float(f2), S),
        "chao1_excess_over_observed": chao - S,
        "chao1_ratio_to_observed": safe_div(chao, S),
    })
    for q in [0.0, 0.5, 1.0, 2.0, 3.0]:
        out[f"hill_D{sanitize_q(q)}"] = hill_number(p, q)
    for n in [1, 5, 10, 20, 50, 100]:
        out[f"top{n}_frequency"] = float(p_desc[:n].sum())
    for t in [0.25, 0.50, 0.75]:
        k = clones_needed_for_fraction(p_desc, t)
        out[f"clones_to_{int(round(100*t))}pct"] = k
        out[f"fraction_clones_to_{int(round(100*t))}pct"] = k / S
    return out


TAIL_EMPTY = {
    "tail_fit_success": False,
    "tail_fit_message": "",
    "tail_fit_library": "",
    "tail_fit_data": "clone_read_counts",
    "tail_fit_discrete": True,
    "tail_xmin_count": np.nan,
    "tail_xmin_frequency": np.nan,
    "tail_n": np.nan,
    "tail_fraction_clonotypes": np.nan,
    "tail_fraction_reads": np.nan,
    "tail_pl_alpha": np.nan,
    "tail_pl_sigma": np.nan,
    "tail_pl_ks_D": np.nan,
    "tail_pl_loglik": np.nan,
    "tail_tpl_alpha": np.nan,
    "tail_tpl_lambda": np.nan,
    "tail_tpl_loglik": np.nan,
    "tail_lognormal_mu": np.nan,
    "tail_lognormal_sigma": np.nan,
    "tail_lognormal_loglik": np.nan,
    "tail_pl_aicc": np.nan,
    "tail_tpl_aicc": np.nan,
    "tail_lognormal_aicc": np.nan,
    "tail_best_model_aicc": "",
    "tail_pl_vs_tpl_R": np.nan,
    "tail_pl_vs_tpl_p": np.nan,
    "tail_pl_vs_lognormal_R": np.nan,
    "tail_pl_vs_lognormal_p": np.nan,
    "tail_min_size_requirement": np.nan,
    "tail_meets_min_size": False,
}

TAIL_EMPTY.update({"tail_fit_warnings": "", "tail_optimizer_warning": False,
                   "tail_aicc_comparison_complete": False})

def aicc(loglik: float, k: int, n: int) -> float:
    if not np.isfinite(loglik) or n <= 0:
        return np.nan
    base = 2.0 * k - 2.0 * loglik
    denom = n - k - 1
    if denom <= 0:
        return np.nan
    return float(base + (2.0 * k * (k + 1)) / denom)

def _extract_attr(obj, names: Sequence[str], default=np.nan):
    for name in names:
        if hasattr(obj, name):
            value = getattr(obj, name)
            try:
                return float(value)
            except Exception:
                return value
    return default

def _distribution_loglik(dist, tail_data: np.ndarray) -> float:
    """Never silently discard invalid per-clonotype log likelihoods.

    All candidate models must be evaluated on exactly the same observations.
    """
    try:
        vals = np.asarray(dist.loglikelihoods(tail_data), dtype=float).reshape(-1)
        if vals.size != tail_data.size or not np.all(np.isfinite(vals)):
            return np.nan
        return float(vals.sum())
    except Exception:
        return np.nan


def fit_heavy_tail(
    counts: Sequence[float],
    min_tail_size: int,
    fast_discrete: bool,
) -> Dict[str, object]:
    result = dict(TAIL_EMPTY)
    result["tail_min_size_requirement"] = int(min_tail_size)

    c = np.asarray(counts, dtype=float)
    c = c[np.isfinite(c) & (c > 0)]
    if c.size < min_tail_size:
        result["tail_fit_message"] = (
            f"Not enough clonotypes for tail fitting: {c.size} < {min_tail_size}"
        )
        return result

    if not np.allclose(c, np.rint(c), rtol=0, atol=1e-8):
        result["tail_fit_message"] = (
            "Tail fitting skipped: count column is not integer-like; "
            "discrete clone-size fitting requires integer counts."
        )
        return result

    ci = np.rint(c).astype(np.int64)

    try:
        import powerlaw  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Heavy-tail fitting requires the optional package 'powerlaw'.\n"
            "Install it with:\n"
            "    python3 -m pip install \"powerlaw==2.0.0\"\n"
            "or rerun with --skip_tail_fitting."
        ) from exc

    version = getattr(powerlaw, "__version__", "unknown")
    result["tail_fit_library"] = f"powerlaw {version}"

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = powerlaw.Fit(
                ci,
                discrete=True,
                verbose=False,
                estimate_discrete=bool(fast_discrete),
                discrete_normalization="round",
            )

            pl = fit.power_law
            xmin = float(fit.xmin)
            tail = ci[ci >= xmin]

            result["tail_xmin_count"] = xmin
            result["tail_xmin_frequency"] = safe_div(xmin, float(ci.sum()))
            result["tail_n"] = int(tail.size)
            result["tail_fraction_clonotypes"] = safe_div(tail.size, ci.size)
            result["tail_fraction_reads"] = safe_div(float(tail.sum()), float(ci.sum()))
            result["tail_meets_min_size"] = bool(tail.size >= min_tail_size)

            result["tail_pl_alpha"] = _extract_attr(pl, ["alpha"])
            result["tail_pl_sigma"] = _extract_attr(pl, ["standard_err", "sigma"])
            result["tail_pl_ks_D"] = _extract_attr(pl, ["D"])
            if not np.isfinite(result["tail_pl_ks_D"]):
                try:
                    pl.compute_distance_metrics(tail)
                    result["tail_pl_ks_D"] = _extract_attr(pl, ["D"])
                except Exception:
                    pass

            result["tail_pl_loglik"] = _distribution_loglik(pl, tail)

            # Alternative models are fitted on the same x >= xmin support
            tpl = fit.truncated_power_law
            logn = fit.lognormal

            result["tail_tpl_alpha"] = _extract_attr(tpl, ["alpha"])
            result["tail_tpl_lambda"] = _extract_attr(
                tpl, ["Lambda", "lambda_", "lambda"]
            )
            result["tail_tpl_loglik"] = _distribution_loglik(tpl, tail)

            result["tail_lognormal_mu"] = _extract_attr(logn, ["mu"])
            result["tail_lognormal_sigma"] = _extract_attr(logn, ["sigma"])
            result["tail_lognormal_loglik"] = _distribution_loglik(logn, tail)

            n_tail = int(tail.size)
            result["tail_pl_aicc"] = aicc(result["tail_pl_loglik"], 1, n_tail)
            result["tail_tpl_aicc"] = aicc(result["tail_tpl_loglik"], 2, n_tail)
            result["tail_lognormal_aicc"] = aicc(
                result["tail_lognormal_loglik"], 2, n_tail
            )

            model_aicc = {
                "power_law": result["tail_pl_aicc"],
                "truncated_power_law": result["tail_tpl_aicc"],
                "lognormal": result["tail_lognormal_aicc"],
            }
            finite_models = {
                k: v for k, v in model_aicc.items() if np.isfinite(v)
            }
            complete = len(finite_models) == 3
            result["tail_aicc_comparison_complete"] = complete
            result["tail_best_model_aicc"] = (
                min(finite_models, key=finite_models.get) if complete else ""
            )
            result["tail_optimizer_warning"] = any(
                bool(getattr(obj, "noise_flag", False)) for obj in (fit, pl, tpl, logn)
            )

            try:
                R, p = fit.distribution_compare(
                    "power_law",
                    "truncated_power_law",
                    nested=True,
                )
                result["tail_pl_vs_tpl_R"] = float(R)
                result["tail_pl_vs_tpl_p"] = float(p)
            except Exception:
                pass

            try:
                R, p = fit.distribution_compare(
                    "power_law",
                    "lognormal",
                    nested=False,
                )
                result["tail_pl_vs_lognormal_R"] = float(R)
                result["tail_pl_vs_lognormal_p"] = float(p)
            except Exception:
                pass

        result["tail_fit_warnings"] = " | ".join(dict.fromkeys(str(w.message) for w in caught))
        required = ["tail_pl_alpha", "tail_tpl_alpha", "tail_tpl_lambda",
                    "tail_lognormal_mu", "tail_lognormal_sigma",
                    "tail_pl_loglik", "tail_tpl_loglik", "tail_lognormal_loglik"]
        result["tail_fit_success"] = all(np.isfinite(result[k]) for k in required)
        notes = []
        if not result["tail_fit_success"]:
            notes.append("one or more model parameters or full-tail log likelihoods are non-finite")
        if not result["tail_meets_min_size"]:
            notes.append(f"tail smaller than min_tail_size={min_tail_size}; descriptive only")
        if result["tail_optimizer_warning"]:
            notes.append("package optimizer/noise warning; inspect tail_fit_warnings")
        if not result["tail_aicc_comparison_complete"]:
            notes.append("three-model AICc comparison unavailable")
        result["tail_fit_message"] = "; ".join(notes) if notes else "ok"
        return result

    except Exception as exc:
        result["tail_fit_success"] = False
        result["tail_fit_message"] = f"{type(exc).__name__}: {exc}"
        return result

TABLE1_METRICS = [
    ("depth", "Sequencing depth (reads)", 0, True),
    ("observed_richness", "Observed clonotypes", 0, True),
    ("goods_coverage", "Good's coverage", 4, False),
    ("hill_D1", "Hill diversity D1", 1, False),
    ("hill_D2", "Hill diversity D2", 1, False),
    ("top10_frequency", "Top-10 cumulative frequency", 3, False),
]
PAPER_PER_REPERTOIRE_COLUMNS = IDENTITY_COLUMNS + [r[0] for r in TABLE1_METRICS]
LEGACY_TABLE1_METRICS = [
    ("depth", "Sequencing depth (reads)", 0, True),
    ("observed_richness", "Observed clonotypes", 0, True),
    ("goods_coverage", "Good's coverage", 4, False),
    ("shannon_entropy", "Shannon entropy", 3, False),
    ("hill_D1", "Hill diversity D1", 1, False),
    ("inverse_simpson_D2", "Hill diversity D2 (inverse Simpson)", 1, False),
    ("shannon_clonality", "Normalized Shannon clonality", 3, False),
    ("gini_frequency", "Gini coefficient", 3, False),
    ("top10_frequency", "Top-10 cumulative frequency", 3, False),
    ("fraction_clones_to_50pct", "Fraction of clonotypes accounting for 50% of reads", 4, False),
]
LEGACY_PAPER_COLUMNS = IDENTITY_COLUMNS + [r[0] for r in LEGACY_TABLE1_METRICS] + [
    "top100_frequency", "clones_to_50pct"]

def numeric_summary_by_compartment(
    full_df: pd.DataFrame,
    metric_specs=TABLE1_METRICS,
) -> pd.DataFrame:
    rows: List[Dict] = []

    for metric, label, digits, integer in metric_specs:
        if metric not in full_df.columns:
            continue
        for compartment, group in full_df.groupby("compartment", observed=True):
            vals = pd.to_numeric(group[metric], errors="coerce")
            vals = vals[np.isfinite(vals)]
            if vals.empty:
                summary = {
                    "n": 0,
                    "median": np.nan,
                    "q1": np.nan,
                    "q3": np.nan,
                    "mean": np.nan,
                    "sd": np.nan,
                    "min": np.nan,
                    "max": np.nan,
                }
            else:
                summary = {
                    "n": int(vals.size),
                    "median": float(vals.median()),
                    "q1": float(vals.quantile(0.25)),
                    "q3": float(vals.quantile(0.75)),
                    "mean": float(vals.mean()),
                    "sd": float(vals.std(ddof=1)) if vals.size > 1 else np.nan,
                    "min": float(vals.min()),
                    "max": float(vals.max()),
                }

            rows.append(
                {
                    "metric": metric,
                    "label": label,
                    "compartment": compartment,
                    "display_digits": digits,
                    "display_integer": integer,
                    **summary,
                }
            )

    return pd.DataFrame(rows)

def build_paper_ready_table(
    numeric_summary: pd.DataFrame,
    full_df: pd.DataFrame,
    metric_specs=TABLE1_METRICS,
) -> pd.DataFrame:
    compartments = [x for x in ["PBMC", "BIOPSY"] if x in set(full_df["compartment"])]
    rows: List[Dict[str, str]] = []

    n_row: Dict[str, str] = {"Characteristic": "Repertoires, n"}
    for comp in compartments:
        n_row[comp] = str(int((full_df["compartment"] == comp).sum()))
    rows.append(n_row)

    subj_row: Dict[str, str] = {"Characteristic": "Subjects represented, n"}
    for comp in compartments:
        subj_row[comp] = str(
            int(full_df.loc[full_df["compartment"] == comp, "subject"].nunique())
        )
    rows.append(subj_row)

    for metric, label, digits, integer in metric_specs:
        row: Dict[str, str] = {"Characteristic": label}
        for comp in compartments:
            hit = numeric_summary[
                (numeric_summary["metric"] == metric)
                & (numeric_summary["compartment"] == comp)
            ]
            if hit.empty:
                row[comp] = "NA"
            else:
                r = hit.iloc[0]
                row[comp] = format_median_iqr(
                    float(r["median"]),
                    float(r["q1"]),
                    float(r["q3"]),
                    digits=digits,
                    integer=integer,
                )
        rows.append(row)

    out = pd.DataFrame(rows)
    out["Summary"] = "median [IQR]"
    out.loc[out["Characteristic"].str.endswith(", n"), "Summary"] = "n"
    return out

def selection_policy() -> pd.DataFrame:
    rows = [
        ("shannon_entropy", "hill_D1", "exact_transform", "H = ln(D1)", "D1 retained; entropy not independently exported."),
        ("simpson_concentration", "hill_D2", "exact_transform", "concentration = 1/D2", "Effective diversity D2 retained."),
        ("gini_simpson", "hill_D2", "exact_transform", "Gini-Simpson = 1 - 1/D2", "Not the Gini coefficient; redundant with D2."),
        ("inverse_simpson_D2", "hill_D2", "exact_alias", "inverse_Simpson = D2", "Single canonical column hill_D2."),
        ("hill_D0", "observed_richness", "exact_alias", "D0 = S", "Richness retained once."),
        ("pielou_evenness", "hill_D1;observed_richness", "jointly_derived", "J = ln(D1)/ln(S), S>1", "Not an extra independent architecture feature."),
        ("shannon_clonality", "hill_D1;observed_richness", "jointly_derived", "C = 1 - ln(D1)/ln(S), S>1", "Recoverable in isolated legacy export."),
        ("berger_parker_dominance", "top1_frequency (only when requested)", "exact_alias", "Berger-Parker = p_max = top1", "Do not export duplicate maximum-frequency names."),
        ("clones_to_Xpct", "fraction_clones_to_Xpct;observed_richness", "jointly_derived", "K_t = S * F_t", "Fraction retained; count omitted from canonical tables."),
        ("singleton_fraction_clonotypes", "singleton_count;observed_richness", "jointly_derived", "f1/S", "Counts retained as sampling diagnostics."),
        ("doubleton_fraction_clonotypes", "doubleton_count;observed_richness", "jointly_derived", "f2/S", "Counts retained as sampling diagnostics."),
        ("chao1_excess_over_observed", "chao1_richness;observed_richness", "jointly_derived", "Chao1-S", "Keep one Chao1 representation."),
        ("chao1_ratio_to_observed", "chao1_richness;observed_richness", "jointly_derived", "Chao1/S", "Keep one Chao1 representation."),
        ("hill_D0p5;hill_D3", "hill_D1;hill_D2", "parsimony_not_identity", "No exact identity with D1 and D2 in general", "Not defaults; may be explicitly requested."),
        ("top1;top5;top20;top50", "top10_frequency;top100_frequency", "parsimony_not_identity", "Different cutoffs of the same rank-abundance curve", "Not defaults; may be explicitly requested."),
        ("fraction_clones_to_25pct;fraction_clones_to_75pct", "fraction_clones_to_50pct", "parsimony_not_identity", "Different cumulative targets", "Not defaults; may be explicitly requested."),
        ("gini_frequency;top100_frequency;fraction_clones_to_50pct", "retained in full table", "supplementary_descriptor", "Not algebraic aliases of D1 or D2", "Not additional independent discoveries; not in default paper Table 1."),
        ("goods_coverage;chao1_richness", "retained with sampling inputs", "sampling_diagnostic", "Coverage=1-f1/N; Chao1=g(S,f1,f2)", "Intentional diagnostic formulas, not independent architecture features."),
    ]
    return pd.DataFrame(rows, columns=["original_metric", "canonical_representation", "relationship", "identity_or_note", "decision"])


def metric_definitions(hill_orders, top_n_values, cumulative_targets) -> pd.DataFrame:
    rows = [
        ("depth", "sum_i c_i", "Total positive retained read counts after within-file collapsing.", "sampling", "paper_table"),
        ("counts_integer_like", "all abs(c_i-round(c_i)) <= 1e-8; integer-like total", "Count-integrity diagnostic; not biological QC.", "input", "supplementary"),
        ("observed_richness", "S", "Observed unique aaCDR3 sequences; equals Hill D0.", "richness", "paper_table"),
        ("singleton_count", "sum_i I(c_i=1)", "Sequencing singleton counts; not one-cell clones.", "sampling", "supplementary"),
        ("doubleton_count", "sum_i I(c_i=2)", "Sequencing doubleton counts; not two-cell clones.", "sampling", "supplementary"),
        ("goods_coverage", "1-f1/N", "Read-sampling coverage descriptor, not fraction of biological clonotypes recovered.", "sampling", "paper_table"),
        ("chao1_richness", "S+f1^2/(2*f2); if f2=0: S+f1*(f1-1)/2", "Sampling-sensitive descriptive estimator on read counts.", "sampling", "supplementary"),
        ("gini_frequency", "2*sum_i i*p_(i)/(S*sum_i p_(i))-(S+1)/S; ascending order", "Uncorrected inequality over observed clonotypes; distinct from Gini-Simpson.", "inequality", "supplementary"),
    ]
    for q in sorted(set([1.0, 2.0] + list(hill_orders))):
        if abs(q) < 1e-12:
            continue
        formula = "exp(-sum_i p_i*ln(p_i))" if abs(q-1)<1e-12 else "(sum_i p_i^q)^(1/(1-q))"
        rows.append((f"hill_D{sanitize_q(q)}", formula, f"Effective number of clonotypes at q={q}; not depth-corrected.", "diversity", "paper_table" if q in (1,2) else "supplementary"))
    for n in sorted(set([10]+list(top_n_values))):
        rows.append((f"top{n}_frequency", f"sum of largest min({n},S) frequencies", "Cumulative observed read abundance; fraction, not percentage.", "dominance", "paper_table" if n==10 else "supplementary"))
    for t in sorted(set(cumulative_targets)):
        rows.append((f"fraction_clones_to_{target_label(t)}pct", f"min{{k: sum_(i=1)^k p_[i]>={t}}}/S", "Fraction of observed clonotypes needed to reach the target, descending order.", "dominance", "supplementary"))
    tail_desc = {
        "tail_fit_success": ("all required parameters and complete-tail log likelihoods finite", "Numerical completion, not a biological QC or goodness-of-fit test."),
        "tail_fit_message": ("text", "Numerical result and descriptive-size notes."),
        "tail_fit_library": ("package version", "powerlaw runtime version."),
        "tail_fit_data": ("clone_read_counts", "One observation per clonotype, not one per read."),
        "tail_fit_discrete": ("True", "Fit to integer-valued counts."),
        "tail_xmin_count": ("xmin selected by power-law KS minimization", "Common lower support used for all three models."),
        "tail_xmin_frequency": ("xmin/N", "Descriptive normalization of the fitted cutoff."),
        "tail_n": ("sum_i I(c_i>=xmin)", "Number of tail clonotypes; sample size for AICc."),
        "tail_fraction_clonotypes": ("n_tail/S", "Fraction of observed clonotypes in tail."),
        "tail_fraction_reads": ("sum_(c_i>=xmin) c_i/N", "Read abundance represented by the tail."),
        "tail_pl_alpha": ("alpha_PL", "Power-law exponent."),
        "tail_pl_sigma": ("package sigma", "Package-reported uncertainty of the power-law exponent."),
        "tail_pl_ks_D": ("sup_x |F_emp(x)-F_PL(x)|", "KS distance; not a calibrated goodness-of-fit p-value."),
        "tail_tpl_alpha": ("alpha_TPL", "Truncated power-law exponent."),
        "tail_tpl_lambda": ("lambda_TPL", "Exponential truncation rate."),
        "tail_lognormal_mu": ("mu", "Lognormal log-scale location."),
        "tail_lognormal_sigma": ("sigma", "Lognormal log-scale dispersion."),
        "tail_best_model_aicc": ("argmin AICc, only if all three finite", "Descriptive preference, not evidence of a biological mechanism."),
        "tail_min_size_requirement": ("CLI setting; default 50", "Descriptive size requirement; not repertoire exclusion."),
        "tail_meets_min_size": ("n_tail>=min_tail_size", "Descriptive tail-size flag."),
        "tail_fit_warnings": ("captured package warnings", "Inspect numerical warnings before interpreting fits."),
        "tail_optimizer_warning": ("any package noise_flag", "Numerical diagnostic, not sample QC."),
        "tail_aicc_comparison_complete": ("all three AICc values finite", "No model winner is assigned to an incomplete comparison."),
    }
    for stem,k in [("pl",1),("tpl",2),("lognormal",2)]:
        tail_desc[f"tail_{stem}_loglik"]=("sum_i log P_model(c_i|c_i>=xmin)", "All tail observations must have finite log likelihoods.")
        tail_desc[f"tail_{stem}_aicc"]=(f"2*{k}-2*loglik+2*{k}*({k}+1)/(n_tail-{k}-1)", "Conditional on selected xmin; selection uncertainty not penalized.")
    for alt in ["tpl","lognormal"]:
        tail_desc[f"tail_pl_vs_{alt}_R"]=("LL_PL-LL_alternative", "Positive favors PL in likelihood, negative alternative; not AICc.")
        tail_desc[f"tail_pl_vs_{alt}_p"]=("package likelihood-ratio comparison", "Nominal descriptive comparison; not a tissue-comparison or absolute GOF test.")
    for name in TAIL_EMPTY:
        formula, note = tail_desc[name]
        rows.append((name,formula,note,"tail","supplementary"))
    return pd.DataFrame(rows, columns=["metric","formula","interpretation","category","reporting_role"])


def package_version(name: str):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def save_tables(full_df: pd.DataFrame, out_dir: Path, legacy: bool = False) -> None:
    specs = LEGACY_TABLE1_METRICS if legacy else TABLE1_METRICS
    cols = LEGACY_PAPER_COLUMNS if legacy else PAPER_PER_REPERTOIRE_COLUMNS
    ensure_dir(out_dir)
    full_df.to_csv(out_dir / "repertoire_characteristics_full.csv", index=False)
    full_df[[c for c in cols if c in full_df.columns]].to_csv(
        out_dir / "repertoire_characteristics_table1_by_repertoire.csv", index=False)
    summary = numeric_summary_by_compartment(full_df, metric_specs=specs)
    summary.to_csv(out_dir / "table1_summary_by_compartment_numeric.csv", index=False)
    build_paper_ready_table(summary, full_df, metric_specs=specs).to_csv(
        out_dir / "table1_paper_ready.csv", index=False)


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=SCRIPT_VERSION)
    ap.add_argument("--data_dir", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--pattern", default=DEFAULT_PATTERN, help="Filename regex; named subject and compartment groups are required.")
    ap.add_argument("--file_sep", default="\t")
    ap.add_argument("--clone_col", default="aaSeqCDR3")
    ap.add_argument("--count_col", default="readCount")
    ap.add_argument("--hill_orders", default="1,2", help="Additional/configured orders; D1 and D2 always included; D0 represented only by richness.")
    ap.add_argument("--top_n", default="10,100", help="Top-N frequencies; top-10 always included.")
    ap.add_argument("--cumulative_targets", default="0.50", help="Targets in (0,1); export fractions only, not duplicate raw counts.")
    ap.add_argument("--skip_tail_fitting", action="store_true")
    ap.add_argument("--min_tail_size", type=int, default=50)
    ap.add_argument("--exact_discrete_tail_fit", action="store_true", help="Disable the package's fast discrete power-law exponent estimate; NOT a guarantee of exact discretization for every model.")
    ap.add_argument("--write_legacy_compat", action="store_true", help="Write isolated legacy_compat/ tables with original default metric names. Do not use these as a nonredundant feature matrix.")
    ap.add_argument("--overwrite", action="store_true", help="Allow rewriting an existing non-empty output directory.")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    if not args.data_dir.is_dir():
        raise NotADirectoryError(args.data_dir)
    if args.data_dir.resolve() == args.out_dir.resolve():
        raise ValueError("Input and output directories must differ.")
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"Output directory is not empty: {args.out_dir}. Use a new directory or --overwrite.")
    hill_orders = parse_float_list(args.hill_orders)
    top_n_values = parse_int_list(args.top_n)
    cumulative_targets = parse_float_list(args.cumulative_targets)
    if not hill_orders or any(not np.isfinite(q) or q<0 for q in hill_orders):
        raise ValueError("--hill_orders must contain finite nonnegative values.")
    if not top_n_values or any(n<=0 for n in top_n_values):
        raise ValueError("--top_n must contain positive integers.")
    if not cumulative_targets or any(not np.isfinite(t) or not 0<t<1 for t in cumulative_targets):
        raise ValueError("--cumulative_targets must be finite and strictly between 0 and 1.")
    if args.min_tail_size<2:
        raise ValueError("--min_tail_size must be >=2.")
    labels = [target_label(t) for t in sorted(set(cumulative_targets))]
    if len(labels)!=len(set(labels)):
        raise ValueError("Cumulative targets are too close for distinct output column names.")
    if 0.0 in hill_orders:
        print("[INFO] Hill D0 is represented by observed_richness; no duplicate hill_D0 column.")
    hill_orders = sorted(set([1.0,2.0]+[q for q in hill_orders if q>=1e-12]))
    if len({sanitize_q(q) for q in hill_orders}) != len(hill_orders):
        raise ValueError("Hill orders would create duplicate column names.")
    top_n_values = sorted(set([10]+top_n_values))
    cumulative_targets = sorted(set(cumulative_targets))
    # Fail before partially writing results when a requested dependency is absent.
    if not args.skip_tail_fitting:
        try:
            import powerlaw  # noqa: F401
        except ImportError as exc:
            raise RuntimeError('Install powerlaw (python3 -m pip install "powerlaw>=2.0.0") or use --skip_tail_fitting.') from exc
    manifest = index_repertoire_files(args.data_dir, args.pattern)
    n_complete = int(manifest.loc[manifest["pair_complete"],"subject"].nunique())
    matched = set(manifest["file"])
    ignored = sorted(p.name for p in args.data_dir.iterdir() if p.is_file() and p.suffix.lower()==".tsv" and p.name not in matched)
    print(f"[INFO] {len(manifest)} repertoires; {manifest.subject.nunique()} subjects; {n_complete} complete pairs.")
    if ignored:
        print(f"[WARNING] {len(ignored)} TSV files do not match --pattern; listed in run_config.json.")
    print("[INFO] Separate biological repertoires; no technical-replicate or sample-exclusion analysis.")
    rows, legacy_rows = [], []
    for i, rec in manifest.iterrows():
        fp = Path(rec["file_path"])
        print(f"[INFO] [{i+1}/{len(manifest)}] {fp.name}", flush=True)
        rep, stats = read_repertoire(fp, args.file_sep, args.clone_col, args.count_col)
        if rep.empty:
            raise ValueError(f"{fp.name}: no valid positive-count clonotype remains. Input error; no sample is silently dropped.")
        counts = rep["count"].to_numpy(dtype=float)
        metrics = compute_repertoire_metrics(counts, hill_orders, top_n_values, cumulative_targets)
        if args.skip_tail_fitting:
            tail = dict(TAIL_EMPTY)
            tail.update(tail_fit_message="skipped by --skip_tail_fitting", tail_min_size_requirement=args.min_tail_size)
        else:
            tail = fit_heavy_tail(counts, args.min_tail_size, not args.exact_discrete_tail_fit)
        base = {"subject":str(rec["subject"]), "compartment":str(rec["compartment"]),
                "file":str(rec["file"]), "pair_complete":bool(rec["pair_complete"]),
                "clone_definition":args.clone_col, "count_definition":args.count_col, **stats}
        rows.append({**base, **metrics, **tail})
        if args.write_legacy_compat:
            legacy_rows.append({**base, **legacy_metrics(counts,metrics), **tail})
    full_df = pd.DataFrame(rows)  # manifest already sorted by subject then PBMC/BIOPSY
    out_dir = ensure_dir(args.out_dir)
    # Invalidate an old completion record before replacing any output.
    (out_dir / "run_config.json").unlink(missing_ok=True)
    manifest.to_csv(out_dir / "input_manifest.csv", index=False)
    save_tables(full_df,out_dir)
    defs = metric_definitions(hill_orders,top_n_values,cumulative_targets)
    defs.to_csv(out_dir / "repertoire_metric_definitions.csv", index=False)
    selection_policy().to_csv(out_dir / "metric_selection_and_redundancy.csv", index=False)
    config = {
        "script":Path(__file__).name, "script_version":SCRIPT_VERSION, "schema_version":SCHEMA_VERSION,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "completed_utc":datetime.now(timezone.utc).isoformat(),
        "python_version":platform.python_version(), "platform":platform.platform(),
        "numpy_version":np.__version__, "pandas_version":pd.__version__,
        "scipy_version":package_version("scipy"), "powerlaw_version":package_version("powerlaw"),
        "data_dir":str(args.data_dir.resolve()), "out_dir":str(out_dir.resolve()),
        "pattern":args.pattern, "file_sep":args.file_sep,
        "clone_col":args.clone_col, "count_col":args.count_col,
        "hill_orders":hill_orders, "top_n":top_n_values, "cumulative_targets":cumulative_targets,
        "tail_fitting":not args.skip_tail_fitting, "tail_data":"discrete clone read counts",
        "min_tail_size":args.min_tail_size, "exact_discrete_tail_fit":bool(args.exact_discrete_tail_fit),
        "estimate_discrete":not args.exact_discrete_tail_fit, "discrete_normalization":"round",
        "n_input_repertoires":len(full_df), "n_subjects":int(full_df.subject.nunique()),
        "n_complete_pbmc_biopsy_pairs":n_complete, "ignored_tsv_files":ignored,
        "paper_table_metrics":[r[0] for r in TABLE1_METRICS],
        "primary_architecture_descriptors":["hill_D1","hill_D2","top10_frequency"],
        "legacy_compatibility_export":bool(args.write_legacy_compat),
        "productivity_filter_applied":False, "rarefaction_applied":False,
        "analysis_note":"Descriptive within-repertoire analysis. No biological QC exclusions, overlap tests, or inferential PBMC-biopsy comparisons. Sampling diagnostics and tail summaries are not independent architecture features.",
    }
    if args.write_legacy_compat:
        legacy_dir=ensure_dir(out_dir/"legacy_compat")
        save_tables(pd.DataFrame(legacy_rows),legacy_dir,legacy=True)
        manifest.to_csv(legacy_dir/"input_manifest.csv",index=False)
        legacy_config={**config,"schema_version":"legacy_metric_compatibility_only", "out_dir":str(legacy_dir.resolve()),
                       "analysis_note":"REDUNDANT COMPATIBILITY TABLES. Not the canonical nonredundant feature set. Original default metric vocabulary restored; tail diagnostics reflect v2 numerical checks."}
        (legacy_dir/"run_config.json").write_text(json.dumps(legacy_config,indent=2),encoding="utf-8")
        (legacy_dir/"README.txt").write_text("Compatibility export only. Do not use as a nonredundant feature matrix. See ../metric_selection_and_redundancy.csv.\n",encoding="utf-8")
    (out_dir/"run_config.json").write_text(json.dumps(config,indent=2,ensure_ascii=False),encoding="utf-8")
    print(f"[OK] Step 1 v{SCRIPT_VERSION} completed: {out_dir}")
    print("[OK] Canonical outputs: 7 CSV files and run_config.json; no plots.")
    if args.write_legacy_compat:
        print("[WARNING] Legacy tables contain duplicates by design; use only to migrate old downstream scripts.")


if __name__ == "__main__":
    main()
