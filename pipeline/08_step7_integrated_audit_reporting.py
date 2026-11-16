#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STEP 7 v2 — Integrated audit + publication reporting for Steps 3–6
==================================================================

Purpose
-------
Consume the FINAL frozen outputs from:
  * Step 3 v3: audited subject-level master dataset;
  * Step 4 v3: outcome-blind four-domain immune-panel audit/freeze;
  * Step 5 v1.1: outcome-blind OCA genomic-panel audit/freeze;
  * Step 6 v2: prespecified immune–genomic association analysis;

and generate publication-ready figures, a main Table 2, supplementary tables,
Supplementary Data 3, captions, numerical summaries, and provenance records.

This is a REPORTING layer only. It does not rerun feature selection, matched-depth
resampling, Spearman tests, permutations, bootstrap, LOSO, or any other inferential
analysis. All displayed statistics are read from the frozen upstream outputs.

Expected inputs
---------------
--step3  Step-3 v3 directory or ZIP.
--step4  Step-4 v3 directory or ZIP.
--step5  Step-5 v1.1 directory or ZIP.
--step6  Step-6 v2 directory or ZIP.
--subject_map (optional)  CSV/TSV/XLSX with columns subject, public_id.

Main outputs
------------
Table 2      Outcome-blind frozen immune and genomic panels, including candidate
             domains, selected representatives, LOSO stability, and audit flags.
Figure 2A    Primary 4x4 immune–genomic Spearman heatmap.
Figure 2B    Hill D1 versus TMB.
Figure 2C    Top-10 frequency versus TMB.
Figure 2D    Hill D1–TMB full-cohort and LOSO stability.

Supplementary figures
---------------------
Figure S3A   Step-4 within-domain candidate representativeness and selected features.
Figure S3B   Step-4 LOSO selection stability of the four core immune domains.
Figure S3C   Correlation structure of the four frozen immune representatives.
Figure S3D   Clonotype denominator support for Top-1%, Top-5%, and Top-10% recovery.
Figure S4A   Correlation structure of the four frozen continuous genomic estimands.
Figure S5A   Forest plot of all 16 primary immune–genomic associations.
Figure S5B   Secondary 3x4 recovery-sensitivity association heatmap.
Figure S5C   Exploratory POLE rank-biserial effects.

Tables
------
Table 2      Compact main-text frozen-panel audit.
Table S4     Step-3 cohort integration and cross-step audit summary.
Table S5     Step-4 core immune-domain candidate audit and selection.
Table S6     Step-5 genomic-panel technical audit.
Table S7     All 16 primary immune–genomic associations.
Table S8     Secondary recovery-sensitivity family (12 tests).
Table S9     Exploratory POLE family (4 tests).

Supplementary Data 3
--------------------
Machine-readable reporting copies of relevant Step-3 audit outputs, the complete
Step-4 audit/freeze outputs, Step-5 genomic audit outputs, and Step-6 association
outputs. Source identifiers are never silently relabelled; publication-facing subject
files are marked REVIEW_SOURCE_IDS unless a complete --subject_map is supplied.

Aesthetics
----------
Routine appearance is controlled by PLOT_STYLE, FIGURE_SIZE, PANEL_STYLE and LABELS
near the top of this script. Plotting functions should not need editing for ordinary
font, color, marker, size, line-width, or label changes.

Dependencies: numpy, pandas, matplotlib; openpyxl only if --subject_map is XLSX.
Python >=3.9.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import shutil
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec


SCRIPT_VERSION = "v2.0.0-integrated-step3v3-step4v3-step5v1_1-step6v2-2026-10-08"

# =============================================================================
# USER-EDITABLE AESTHETIC DICTIONARIES
# =============================================================================

PLOT_STYLE = {
    "show_panel_labels": False,
    "show_titles": False,
    "font": {
        "family": "Arial",
        "base_size": 8.0,
        "axis_label_size": 8.5,
        "tick_size": 7.5,
        "legend_size": 7.5,
        "annotation_size": 7.4,
        "panel_label_size": 11.0,
        "weight": "normal",
    },
    "colors": {
        "axis": "#2F2F2F",
        "grid": "#D9D9D9",
        "zero": "#777777",
        "ci": "#777777",
        "full_effect": "#111111",
        "nonsig_edge": "#707070",
        "audit_selected": "#1F4E79",
        "audit_unselected": "#B9C2CC",
        "audit_review": "#C17D11",
        "heatmap_cmap": "RdBu_r",
        "immune": {
            "biopsy_hill_D1": "#3B6FB6",
            "biopsy_top10_frequency": "#C34F4F",
            "within_top_5pct_pbmc_detection_fraction_std_median": "#4C956C",
            "within_morisita_horn_std_median": "#8064A2",
        },
        "genomic": {
            "oca_tmb": "#D07A00",
            "oca_msi_score": "#4E79A7",
            "oca_gim": "#59A14F",
            "oca_loh": "#9C6ADE",
        },
    },
    "markers": {
        "scatter": "o",
        "scatter_size": 34,
        "forest": "o",
        "forest_size": 34,
        "loso": "o",
        "loso_size": 30,
        "full": "D",
        "full_size": 46,
    },
    "lines": {
        "axis_width": 0.8,
        "grid_width": 0.55,
        "zero_width": 0.8,
        "ci_width": 1.2,
        "full_width": 1.2,
        "median_width": 1.0,
    },
    "axes": {
        "remove_top_spine": True,
        "remove_right_spine": True,
        "show_grid": False,
        "grid_alpha": 0.35,
    },
    "heatmap": {
        "vmin": -1.0,
        "vmax": 1.0,
        "cell_text_digits": 2,
        "bh_star_threshold": 0.05,
        "maxT_border_threshold": 0.05,
        "cell_edgecolor": "white",
        "cell_linewidth": 1.0,
    },
    "annotations": {
        "show_stats_on_scatter": True,
        "show_subject_ids": False,
        "stats_box_alpha": 0.88,
        "stats_box_edge": "#C7C7C7",
    },
    "save": {
        "dpi": 400,
        "pdf_fonttype": 42,
        "ps_fonttype": 42,
        "bbox_inches": "tight",
        "transparent": False,
        "pad_inches": 0.04,
    },
}

FIGURE_SIZE = {
    "Fig2A": (4.6, 3.7),
    "Fig2B": (3.6, 3.0),
    "Fig2C": (3.6, 3.0),
    "Fig2D": (4.1, 3.4),
    "Figure2_composite": (9.0, 7.2),
    "FigS5A": (7.2, 7.8),
    "FigS5B": (4.8, 3.4),
    "FigS5C": (5.8, 3.5),
    "FigureS5_composite": (10.2, 7.7),
    "FigS3A": (6.8, 5.8),
    "FigS3B": (5.3, 3.5),
    "FigS3C": (4.8, 4.0),
    "FigS3D": (5.5, 3.7),
    "FigureS3_composite": (10.4, 8.0),
    "FigS4A": (4.5, 3.8),
}

PANEL_STYLE = {
    "Fig2A": {}, "Fig2B": {}, "Fig2C": {}, "Fig2D": {},
    "FigS3A": {}, "FigS3B": {}, "FigS3C": {}, "FigS3D": {},
    "FigS4A": {}, "FigS5A": {}, "FigS5B": {}, "FigS5C": {},
}

LABELS = {
    "immune": {
        "biopsy_hill_D1": "Hill D1",
        "biopsy_top10_frequency": "Top-10 frequency",
        "within_top_5pct_pbmc_detection_fraction_std_median": "Top-5% PBMC recovery",
        "within_morisita_horn_std_median": "Morisita–Horn",
        "within_biopsy_clone_fraction_shared_std_median": "Global clone recovery",
        "within_top_10pct_pbmc_detection_fraction_std_median": "Top-10% PBMC recovery",
        "within_biopsy_abundance_fraction_shared_std_median": "Shared biopsy read abundance",
    },
    "genomic": {
        "oca_tmb": "TMB",
        "oca_msi_score": "MSI score",
        "oca_gim": "GIM",
        "oca_loh": "LOH",
        "oca_pole_positive": "POLE",
    },
    "axes": {
        "oca_tmb": "Tumor mutational burden (mut/Mb)",
        "biopsy_hill_D1": "Tumor Hill D1",
        "biopsy_top10_frequency": "Top-10 read frequency (%)",
    },
}

IMMUNE_ORDER = [
    "biopsy_hill_D1",
    "biopsy_top10_frequency",
    "within_top_5pct_pbmc_detection_fraction_std_median",
    "within_morisita_horn_std_median",
]
GENOMIC_ORDER = ["oca_tmb", "oca_msi_score", "oca_gim", "oca_loh"]
SECONDARY_ORDER = [
    "within_biopsy_clone_fraction_shared_std_median",
    "within_top_10pct_pbmc_detection_fraction_std_median",
    "within_biopsy_abundance_fraction_shared_std_median",
]


# =============================================================================
# Generic helpers
# =============================================================================

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def configure_font(strict: bool) -> bool:
    family = PLOT_STYLE["font"]["family"]
    try:
        font_manager.findfont(family, fallback_to_default=False)
        available = True
    except Exception:
        available = False
    if strict and not available:
        raise RuntimeError(
            f"Requested font '{family}' not found. Install/enable it or use --allow_font_fallback."
        )
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [family, "Liberation Sans", "DejaVu Sans"],
        "font.size": PLOT_STYLE["font"]["base_size"],
        "font.weight": PLOT_STYLE["font"]["weight"],
        "axes.labelsize": PLOT_STYLE["font"]["axis_label_size"],
        "xtick.labelsize": PLOT_STYLE["font"]["tick_size"],
        "ytick.labelsize": PLOT_STYLE["font"]["tick_size"],
        "legend.fontsize": PLOT_STYLE["font"]["legend_size"],
        "pdf.fonttype": PLOT_STYLE["save"]["pdf_fonttype"],
        "ps.fonttype": PLOT_STYLE["save"]["ps_fonttype"],
    })
    return available


def apply_axes_style(ax, grid: Optional[bool] = None):
    cfg = PLOT_STYLE["axes"]
    colors = PLOT_STYLE["colors"]
    if cfg["remove_top_spine"]:
        ax.spines["top"].set_visible(False)
    if cfg["remove_right_spine"]:
        ax.spines["right"].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_linewidth(PLOT_STYLE["lines"]["axis_width"])
        ax.spines[side].set_color(colors["axis"])
    ax.tick_params(
        axis="both",
        width=PLOT_STYLE["lines"]["axis_width"],
        colors=colors["axis"],
        labelsize=PLOT_STYLE["font"]["tick_size"],
    )
    show_grid = cfg["show_grid"] if grid is None else bool(grid)
    if show_grid:
        ax.grid(True, color=colors["grid"], linewidth=PLOT_STYLE["lines"]["grid_width"],
                alpha=cfg["grid_alpha"], zorder=0)
    else:
        ax.grid(False)


def add_panel_label(ax, panel_id: str):
    if not PLOT_STYLE.get("show_panel_labels", False):
        return
    m = re.search(r"([A-Z])$", panel_id)
    if m:
        ax.text(-0.08, 1.06, m.group(1), transform=ax.transAxes,
                ha="left", va="bottom", fontweight="bold",
                fontsize=PLOT_STYLE["font"]["panel_label_size"])


def save_figure(fig, stem: Path):
    cfg = PLOT_STYLE["save"]
    common = dict(
        bbox_inches=cfg["bbox_inches"], transparent=cfg["transparent"],
        pad_inches=cfg["pad_inches"]
    )
    fig.savefig(str(stem) + ".pdf", **common)
    fig.savefig(str(stem) + ".png", dpi=cfg["dpi"], **common)
    plt.close(fig)


def source_members(source: Path) -> List[str]:
    if source.is_dir():
        return [str(p.relative_to(source)) for p in source.rglob("*") if p.is_file()]
    if source.is_file() and source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source, "r") as zf:
            return zf.namelist()
    raise ValueError(f"{source} must be a directory or ZIP")


def _find_member(source: Path, filename: str) -> str:
    if source.is_dir():
        hits = [p for p in source.rglob(filename)
                if "__MACOSX" not in p.parts and not p.name.startswith("._")]
        if len(hits) != 1:
            raise FileNotFoundError(f"Expected exactly one {filename} under {source}; found {len(hits)}")
        return str(hits[0])
    if source.is_file() and source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source, "r") as zf:
            hits = [n for n in zf.namelist()
                    if Path(n).name == filename and "__MACOSX" not in Path(n).parts
                    and not Path(n).name.startswith("._")]
        if len(hits) != 1:
            raise FileNotFoundError(f"Expected exactly one {filename} in {source}; found {len(hits)}")
        return hits[0]
    raise ValueError(f"{source} must be a directory or ZIP")


def read_csv_from_source(source: Path, filename: str, subject_cols: Sequence[str] = ("subject", "left_out_subject")) -> pd.DataFrame:
    member = _find_member(source, filename)
    if source.is_dir():
        return pd.read_csv(member, dtype={c: str for c in subject_cols})
    with zipfile.ZipFile(source, "r") as zf:
        data = zf.read(member)
    return pd.read_csv(io.BytesIO(data), dtype={c: str for c in subject_cols})


def read_json_from_source(source: Path, filename: str = "00_run_config.json") -> Dict:
    member = _find_member(source, filename)
    if source.is_dir():
        return json.loads(Path(member).read_text(encoding="utf-8"))
    with zipfile.ZipFile(source, "r") as zf:
        return json.loads(zf.read(member).decode("utf-8"))


def read_bytes_from_source(source: Path, filename: str) -> bytes:
    member = _find_member(source, filename)
    if source.is_dir():
        return Path(member).read_bytes()
    with zipfile.ZipFile(source, "r") as zf:
        return zf.read(member)


def sha256_source_member(source: Path, filename: str) -> str:
    return hashlib.sha256(read_bytes_from_source(source, filename)).hexdigest()


def require_columns(df: pd.DataFrame, cols: Sequence[str], label: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing columns {missing}")


def read_subject_map(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    if path.suffix.lower() in {".xlsx", ".xls"}:
        d = pd.read_excel(path, dtype=str)
    elif path.suffix.lower() == ".csv":
        d = pd.read_csv(path, dtype=str)
    else:
        d = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(d, ["subject", "public_id"], "subject map")
    d = d[["subject", "public_id"]].copy()
    d["subject"] = d["subject"].astype(str).str.strip()
    d["public_id"] = d["public_id"].astype(str).str.strip()
    if d["subject"].duplicated().any() or d["public_id"].duplicated().any():
        raise ValueError("subject map must be one-to-one")
    return d


def map_subject_column(df: pd.DataFrame, mapping: Optional[pd.DataFrame], column: str) -> pd.DataFrame:
    d = df.copy()
    if column not in d.columns or mapping is None:
        return d
    mp = dict(zip(mapping["subject"], mapping["public_id"]))
    missing = sorted(set(d[column].dropna().astype(str)) - set(mp))
    if missing:
        raise ValueError(f"subject map lacks IDs present in {column}: {missing}")
    d[column] = d[column].astype(str).map(mp)
    return d


def fmt_p(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    if x < 0.001:
        return f"{x:.2e}"
    if x < 0.01:
        return f"{x:.4f}"
    return f"{x:.3f}"


def fmt_num(x: float, digits: int = 3) -> str:
    return "NA" if not np.isfinite(x) else f"{x:.{digits}f}"


def human_immune(feature: str) -> str:
    return LABELS["immune"].get(feature, feature)


def human_genomic(feature: str) -> str:
    return LABELS["genomic"].get(feature, feature)


# =============================================================================
# Validation / preparation
# =============================================================================

def load_inputs(step3: Path, step4: Path, step5: Path, step6: Path) -> Dict[str, object]:
    c3 = read_json_from_source(step3)
    c4 = read_json_from_source(step4)
    c5 = read_json_from_source(step5)
    c6 = read_json_from_source(step6)

    if not str(c3.get("script_version", "")).startswith("v3.0"):
        raise ValueError(f"Step 3 must be v3.0; found {c3.get('script_version')}")
    if not str(c4.get("script_version", "")).startswith("v3.0"):
        raise ValueError(f"Step 4 must be v3.0; found {c4.get('script_version')}")
    if not str(c5.get("script_version", "")).startswith("v1.1"):
        raise ValueError(f"Step 5 must be v1.1; found {c5.get('script_version')}")
    if not str(c6.get("script_version", "")).startswith("v2.0"):
        raise ValueError(f"Step 6 must be v2.0; found {c6.get('script_version')}")

    # Step 3: subject-master and integration audits.
    step3_merge_audit = read_csv_from_source(step3, "04_dataset_merge_audit.csv")
    step3_identity_audit = read_csv_from_source(step3, "05b_repseq_step1_identity_audit.csv")
    step3_denominator_audit = read_csv_from_source(step3, "05c_step1_step2_denominator_audit.csv")
    step3_oca_audit = read_csv_from_source(step3, "05_oca_import_audit.csv")
    step3_dictionary = read_csv_from_source(step3, "06_column_dictionary.csv", subject_cols=())

    # Step 4: complete immune audit/freeze.
    immune_registry = read_csv_from_source(step4, "01_domain_registry.csv", subject_cols=())
    immune_candidate_audit = read_csv_from_source(step4, "02_domain_candidate_audit.csv", subject_cols=())
    immune_depth = read_csv_from_source(step4, "03_depth_associations.csv", subject_cols=())
    immune_within_corr = read_csv_from_source(step4, "04_within_domain_correlations_long.csv", subject_cols=())
    immune_selection = read_csv_from_source(step4, "05_domain_selection_full.csv", subject_cols=())
    immune_loso = read_csv_from_source(step4, "06_domain_selection_loso.csv")
    immune_stability = read_csv_from_source(step4, "07_domain_selection_stability.csv", subject_cols=())
    immune_selected_corr = read_csv_from_source(step4, "08_selected_representative_correlations.csv", subject_cols=())
    immune_redundancy = read_csv_from_source(step4, "09_cross_domain_redundancy_pairs.csv", subject_cols=())
    immune_rec = read_csv_from_source(step4, "10_immune_panel_recommendation.csv", subject_cols=())
    immune_values = read_csv_from_source(step4, "11_immune_panel_values.csv")
    immune_secondary = read_csv_from_source(step4, "12_secondary_domain_audit.csv", subject_cols=())
    immune_support = read_csv_from_source(step4, "14_recovery_denominator_support_subject.csv")

    # Step 5: genomic audit/freeze.
    genomic_registry = read_csv_from_source(step5, "01_genomic_variable_registry.csv", subject_cols=())
    genomic_numeric = read_csv_from_source(step5, "02_genomic_numeric_audit.csv", subject_cols=())
    genomic_categorical = read_csv_from_source(step5, "03_genomic_categorical_audit.csv", subject_cols=())
    genomic_corr = read_csv_from_source(step5, "04_genomic_correlation_spearman.csv", subject_cols=())
    genomic_redundancy = read_csv_from_source(step5, "05_genomic_redundancy_pairs.csv", subject_cols=())
    genomic_rec = read_csv_from_source(step5, "06_genomic_panel_recommendation.csv", subject_cols=())
    genomic_primary = read_csv_from_source(step5, "07_genomic_primary_values.csv")
    genomic_secondary = read_csv_from_source(step5, "08_genomic_secondary_values.csv")

    # Step 6: fixed association families.
    primary = read_csv_from_source(step6, "01_primary_associations.csv")
    rho_matrix = read_csv_from_source(step6, "02_primary_spearman_rho_matrix.csv", subject_cols=())
    perm_matrix = read_csv_from_source(step6, "03_primary_permutation_p_matrix.csv", subject_cols=())
    q_matrix = read_csv_from_source(step6, "04_primary_bh_q_matrix.csv", subject_cols=())
    maxT_matrix = read_csv_from_source(step6, "05_primary_maxT_fwer_p_matrix.csv", subject_cols=())
    primary_loso_long = read_csv_from_source(step6, "06_primary_loso_rho_long.csv")
    primary_loso_summary = read_csv_from_source(step6, "07_primary_loso_summary.csv")
    secondary = read_csv_from_source(step6, "08_secondary_recovery_sensitivity_associations.csv")
    secondary_loso_summary = read_csv_from_source(step6, "09_secondary_recovery_loso_summary.csv")
    secondary_loso_long = read_csv_from_source(step6, "09a_secondary_recovery_loso_rho_long.csv")
    pole = read_csv_from_source(step6, "10_exploratory_pole_associations.csv")
    contract = read_csv_from_source(step6, "11_frozen_panel_contract.csv", subject_cols=())
    matrix = read_csv_from_source(step6, "12_subject_analysis_matrix.csv")
    provenance = read_csv_from_source(step6, "14_upstream_provenance_contract.csv", subject_cols=())

    require_columns(primary, [
        "immune_feature", "genomic_feature", "spearman_rho", "permutation_p_two_sided",
        "bh_q_within_family", "maxT_fwer_p", "bootstrap_ci_low", "bootstrap_ci_high",
        "loso_min_rho", "loso_max_rho", "loso_sign_consistency_fraction"
    ], "primary associations")
    require_columns(matrix, ["subject"] + IMMUNE_ORDER + GENOMIC_ORDER + ["oca_pole_positive"], "subject analysis matrix")
    require_columns(immune_rec, ["domain_id", "selected_feature", "most_selected_fraction_all_runs", "freeze_status"], "Step-4 recommendation")

    if len(primary) != 16:
        raise ValueError(f"Expected 16 primary associations; found {len(primary)}")
    if len(secondary) != 12:
        raise ValueError(f"Expected 12 secondary recovery associations; found {len(secondary)}")
    if len(pole) != 4:
        raise ValueError(f"Expected 4 POLE associations; found {len(pole)}")
    if matrix["subject"].duplicated().any():
        raise ValueError("Subject analysis matrix contains duplicate subjects")
    if matrix[IMMUNE_ORDER + GENOMIC_ORDER + ["oca_pole_positive"]].apply(pd.to_numeric, errors="coerce").isna().any().any():
        raise ValueError("Subject analysis matrix has missing/non-numeric frozen variables")

    expected_pairs = {(i, g) for i in IMMUNE_ORDER for g in GENOMIC_ORDER}
    observed_pairs = set(zip(primary["immune_feature"], primary["genomic_feature"]))
    if observed_pairs != expected_pairs:
        raise ValueError("Primary 4x4 feature family does not match the frozen reporting contract")

    # Cross-step subject sets and Step-3 analysis-dataset provenance.
    subject_sets = {
        "step3": set(step3_merge_audit["subject"].astype(str)),
        "step4": set(immune_values["subject"].astype(str)),
        "step5": set(genomic_primary["subject"].astype(str)),
        "step6": set(matrix["subject"].astype(str)),
    }
    if len({frozenset(v) for v in subject_sets.values()}) != 1:
        raise ValueError(f"Step 3/4/5/6 subject sets differ: {subject_sets}")

    step3_analysis_sha = sha256_source_member(step3, "03_analysis_dataset_subject_level.csv")
    prov_sha_rows = provenance.loc[
        provenance.astype(str).apply(lambda col: col.str.contains(step3_analysis_sha, regex=False, na=False)).any(axis=1)
    ]
    if prov_sha_rows.empty:
        raise ValueError("Step-6 provenance contract does not contain the Step-3 analysis-dataset SHA256")

    # Step-3 audits must be fully successful for the integrated 13-subject cohort.
    if not step3_merge_audit["complete_integrated"].astype(bool).all():
        raise ValueError("Step-3 merge audit is not complete for all subjects")
    if not step3_identity_audit["filename_match"].astype(bool).all():
        raise ValueError("Step-3 RepSeq/Step-1 identity audit contains a mismatch")
    if not step3_denominator_audit["all_checks_pass"].astype(bool).all():
        raise ValueError("Step-3 Step1/Step2 denominator audit contains a failure")

    return locals()


# =============================================================================
# Plotting primitives
# =============================================================================

def association_matrix(df: pd.DataFrame, value: str,
                       row_order: Sequence[str], col_order: Sequence[str]) -> np.ndarray:
    p = df.pivot(index="immune_feature", columns="genomic_feature", values=value)
    return p.reindex(index=row_order, columns=col_order).to_numpy(float)


def draw_assoc_heatmap(ax, results: pd.DataFrame, immune_order: Sequence[str], genomic_order: Sequence[str],
                       panel_id: str, title: str = ""):
    rho = association_matrix(results, "spearman_rho", immune_order, genomic_order)
    q = association_matrix(results, "bh_q_within_family", immune_order, genomic_order)
    maxT = association_matrix(results, "maxT_fwer_p", immune_order, genomic_order)
    cfg = PLOT_STYLE["heatmap"]
    norm = TwoSlopeNorm(vmin=cfg["vmin"], vcenter=0.0, vmax=cfg["vmax"])
    im = ax.imshow(rho, cmap=PLOT_STYLE["colors"]["heatmap_cmap"], norm=norm, aspect="auto")

    ax.set_xticks(np.arange(len(genomic_order)))
    ax.set_yticks(np.arange(len(immune_order)))
    ax.set_xticklabels([human_genomic(x) for x in genomic_order])
    ax.set_yticklabels([human_immune(x) for x in immune_order])
    ax.tick_params(axis="x", rotation=0)

    # White cell boundaries.
    ax.set_xticks(np.arange(-.5, len(genomic_order), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(immune_order), 1), minor=True)
    ax.grid(which="minor", color=cfg["cell_edgecolor"], linestyle="-", linewidth=cfg["cell_linewidth"])
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(rho.shape[0]):
        for j in range(rho.shape[1]):
            val = rho[i, j]
            star = "*" if q[i, j] < cfg["bh_star_threshold"] else ""
            text_color = "white" if abs(val) >= 0.55 else PLOT_STYLE["colors"]["axis"]
            ax.text(j, i, f"{val:.{cfg['cell_text_digits']}f}{star}",
                    ha="center", va="center", color=text_color,
                    fontsize=PLOT_STYLE["font"]["annotation_size"], fontweight="bold" if star else "normal")
            if maxT[i, j] < cfg["maxT_border_threshold"]:
                rect = plt.Rectangle((j - 0.47, i - 0.47), 0.94, 0.94,
                                     fill=False, edgecolor=PLOT_STYLE["colors"]["axis"], linewidth=1.6)
                ax.add_patch(rect)

    if PLOT_STYLE.get("show_titles", False) and title:
        ax.set_title(title)
    add_panel_label(ax, panel_id)
    return im


def make_primary_heatmap(primary: pd.DataFrame, outstem: Path):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["Fig2A"])
    im = draw_assoc_heatmap(ax, primary, IMMUNE_ORDER, GENOMIC_ORDER, "Fig2A", "Primary immune–genomic associations")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Spearman ρ")
    fig.tight_layout()
    save_figure(fig, outstem)


def association_row(primary: pd.DataFrame, immune_feature: str, genomic_feature: str) -> pd.Series:
    d = primary[(primary["immune_feature"] == immune_feature) & (primary["genomic_feature"] == genomic_feature)]
    if len(d) != 1:
        raise ValueError(f"Expected one association for {immune_feature} x {genomic_feature}")
    return d.iloc[0]


def scatter_panel(matrix: pd.DataFrame, primary: pd.DataFrame, immune_feature: str,
                  genomic_feature: str, panel_id: str, outstem: Path,
                  ylabel: str, percent_y: bool = False, log_y: bool = False,
                  public_map: Optional[pd.DataFrame] = None):
    d = matrix[["subject", genomic_feature, immune_feature]].copy()
    x = pd.to_numeric(d[genomic_feature], errors="raise").to_numpy(float)
    y = pd.to_numeric(d[immune_feature], errors="raise").to_numpy(float)
    if percent_y:
        y = 100.0 * y

    fig, ax = plt.subplots(figsize=FIGURE_SIZE[panel_id])
    color = PLOT_STYLE["colors"]["immune"].get(immune_feature, "#4C78A8")
    ax.scatter(x, y, s=PLOT_STYLE["markers"]["scatter_size"], marker=PLOT_STYLE["markers"]["scatter"],
               facecolor=color, edgecolor="white", linewidth=0.7, zorder=3)
    ax.set_xlabel(LABELS["axes"].get(genomic_feature, human_genomic(genomic_feature)))
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")

    if PLOT_STYLE["annotations"]["show_subject_ids"]:
        labels = d["subject"].astype(str)
        if public_map is not None:
            mp = dict(zip(public_map["subject"], public_map["public_id"]))
            labels = labels.map(mp)
        for xi, yi, lab in zip(x, y, labels):
            ax.annotate(str(lab), (xi, yi), xytext=(4, 4), textcoords="offset points",
                        fontsize=PLOT_STYLE["font"]["annotation_size"] - 0.6)

    r = association_row(primary, immune_feature, genomic_feature)
    if PLOT_STYLE["annotations"]["show_stats_on_scatter"]:
        stats = (
            f"ρ = {r['spearman_rho']:.3f}\n"
            f"Permutation P = {fmt_p(float(r['permutation_p_two_sided']))}\n"
            f"BH q = {fmt_p(float(r['bh_q_within_family']))}\n"
            f"maxT P = {fmt_p(float(r['maxT_fwer_p']))}"
        )
        ax.text(0.93, 0.97, stats, transform=ax.transAxes, ha="left", va="top",
                fontsize=PLOT_STYLE["font"]["annotation_size"],
                bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                          edgecolor=PLOT_STYLE["annotations"]["stats_box_edge"],
                          alpha=PLOT_STYLE["annotations"]["stats_box_alpha"], linewidth=0.6))

    apply_axes_style(ax)
    add_panel_label(ax, panel_id)
    fig.tight_layout()
    save_figure(fig, outstem)


def loso_stability_panel(primary: pd.DataFrame, loso_long: pd.DataFrame, outstem: Path,
                         immune_feature: str = "biopsy_hill_D1", genomic_feature: str = "oca_tmb"):
    """Compact robustness strip: full-cohort bootstrap CI plus all LOSO point estimates."""
    full = association_row(primary, immune_feature, genomic_feature)
    d = loso_long[(loso_long["immune_feature"] == immune_feature) &
                  (loso_long["genomic_feature"] == genomic_feature)].copy()
    if d.empty:
        raise ValueError("LOSO long table lacks the primary D1–TMB result")
    d = d.reset_index(drop=True)
    vals = pd.to_numeric(d["rho"], errors="raise").to_numpy(float)
    color = PLOT_STYLE["colors"]["immune"][immune_feature]

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["Fig2D"])
    # Full cohort with bootstrap 95% CI.
    full_rho = float(full["spearman_rho"])
    lo, hi = float(full["bootstrap_ci_low"]), float(full["bootstrap_ci_high"])
    ax.plot([lo, hi], [2.0, 2.0], color=color, linewidth=2.2, alpha=0.85, zorder=2)
    ax.scatter([full_rho], [2.0], s=PLOT_STYLE["markers"]["full_size"],
               marker=PLOT_STYLE["markers"]["full"], color=PLOT_STYLE["colors"]["full_effect"],
               edgecolor="white", linewidth=0.7, zorder=4)

    # LOSO cloud and min-max segment. Deterministic vertical jitter is presentation only.
    jitter = np.linspace(-0.20, 0.20, len(vals))
    ax.plot([float(np.min(vals)), float(np.max(vals))], [1.0, 1.0],
            color=color, linewidth=2.2, alpha=0.55, zorder=1)
    ax.scatter(vals, 1.0 + jitter, s=PLOT_STYLE["markers"]["loso_size"],
               color=color, edgecolor="white", linewidth=0.6, zorder=3)

    ax.axvline(0, color=PLOT_STYLE["colors"]["zero"],
               linewidth=PLOT_STYLE["lines"]["zero_width"], linestyle="--", zorder=0)
    ax.set_xlim(-1.02, 0.08)
    ax.set_ylim(0.55, 2.45)
    ax.set_xlabel("Spearman ρ")
    ax.set_yticks([1.0, 2.0])
    ax.set_yticklabels(["LOSO (13 runs)", "Full cohort"])
    apply_axes_style(ax)
    add_panel_label(ax, "Fig2D")
    fig.tight_layout()
    save_figure(fig, outstem)

def forest_primary(primary: pd.DataFrame, outstem: Path):
    d = primary.copy()
    d["_io"] = d["immune_feature"].map({f: i for i, f in enumerate(IMMUNE_ORDER)})
    d["_go"] = d["genomic_feature"].map({f: i for i, f in enumerate(GENOMIC_ORDER)})
    d = d.sort_values(["_io", "_go"]).reset_index(drop=True)
    labels = [f"{human_immune(i)} — {human_genomic(g)}" for i, g in zip(d.immune_feature, d.genomic_feature)]
    y = np.arange(len(d))[::-1]

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS5A"])
    ax.axvline(0, color=PLOT_STYLE["colors"]["zero"], linestyle="--", linewidth=0.8, zorder=0)
    for k, row in d.iterrows():
        yy = y[k]
        rho = float(row["spearman_rho"])
        lo = float(row["bootstrap_ci_low"])
        hi = float(row["bootstrap_ci_high"])
        color = PLOT_STYLE["colors"]["immune"].get(row["immune_feature"], "#4C78A8")
        ax.plot([lo, hi], [yy, yy], color=color, linewidth=PLOT_STYLE["lines"]["ci_width"], alpha=0.82)
        face = color if float(row["bh_q_within_family"]) < 0.05 else "white"
        ax.scatter([rho], [yy], s=PLOT_STYLE["markers"]["forest_size"], marker=PLOT_STYLE["markers"]["forest"],
                   facecolor=face, edgecolor=color, linewidth=1.1, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(-1.02, 1.02)
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    apply_axes_style(ax)
    add_panel_label(ax, "FigS5A")
    fig.tight_layout()
    save_figure(fig, outstem)


def make_secondary_heatmap(secondary: pd.DataFrame, outstem: Path):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS5B"])
    im = draw_assoc_heatmap(ax, secondary, SECONDARY_ORDER, GENOMIC_ORDER, "FigS5B", "Secondary recovery sensitivity")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Spearman ρ")
    fig.tight_layout()
    save_figure(fig, outstem)


def pole_forest(pole: pd.DataFrame, outstem: Path):
    d = pole.copy()
    d["_o"] = d["immune_feature"].map({f: i for i, f in enumerate(IMMUNE_ORDER)})
    d = d.sort_values("_o").reset_index(drop=True)
    y = np.arange(len(d))[::-1]

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS5C"])
    ax.axvline(0, color=PLOT_STYLE["colors"]["zero"], linestyle="--", linewidth=0.8)
    for k, row in d.iterrows():
        yy = y[k]
        eff = float(row["rank_biserial_effect"])
        lo = float(row["bootstrap_ci_low"])
        hi = float(row["bootstrap_ci_high"])
        color = PLOT_STYLE["colors"]["immune"].get(row["immune_feature"], "#4C78A8")
        ax.plot([lo, hi], [yy, yy], color=color, linewidth=PLOT_STYLE["lines"]["ci_width"])
        ax.scatter([eff], [yy], s=PLOT_STYLE["markers"]["forest_size"], color=color,
                   edgecolor="white", linewidth=0.7, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([human_immune(x) for x in d["immune_feature"]])
    ax.set_xlim(-1.02, 1.02)
    ax.set_xlabel("Rank-biserial effect (bootstrap 95% CI)")
    apply_axes_style(ax)
    add_panel_label(ax, "FigS5C")
    fig.tight_layout()
    save_figure(fig, outstem)


def genomic_corr_heatmap(genomic_corr: pd.DataFrame, outstem: Path):
    d = genomic_corr.copy()
    # Pandas CSV of a correlation matrix stores row names in the first unnamed column.
    first = d.columns[0]
    if first not in GENOMIC_ORDER:
        d = d.set_index(first)
    else:
        d.index = d.columns
    d.index = d.index.astype(str)
    d = d.reindex(index=GENOMIC_ORDER, columns=GENOMIC_ORDER)
    arr = d.to_numpy(float)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS4A"])
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax.imshow(arr, cmap=PLOT_STYLE["colors"]["heatmap_cmap"], norm=norm)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels([human_genomic(x) for x in GENOMIC_ORDER])
    ax.set_yticklabels([human_genomic(x) for x in GENOMIC_ORDER])
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(4):
        for j in range(4):
            val = arr[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    color="white" if abs(val) > 0.55 else PLOT_STYLE["colors"]["axis"],
                    fontsize=PLOT_STYLE["font"]["annotation_size"])
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Spearman ρ")
    add_panel_label(ax, "FigS4A")
    fig.tight_layout()
    save_figure(fig, outstem)


# =============================================================================
# Composite figures
# =============================================================================

def composite_figure2(primary, matrix, loso_long, outstem, public_map=None):
    fig = plt.figure(figsize=FIGURE_SIZE["Figure2_composite"])
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.15, 1.0], height_ratios=[1.0, 1.0],
                  wspace=0.38, hspace=0.38)

    axA = fig.add_subplot(gs[0, 0])
    im = draw_assoc_heatmap(axA, primary, IMMUNE_ORDER, GENOMIC_ORDER, "Fig2A")
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04)
    cb.set_label("Spearman ρ")

    def scatter_into(ax, immune_feature, panel_id, ylabel, percent=False, logy=False):
        d = matrix[["subject", "oca_tmb", immune_feature]].copy()
        x = pd.to_numeric(d["oca_tmb"]).to_numpy(float)
        y = pd.to_numeric(d[immune_feature]).to_numpy(float)
        if percent:
            y *= 100
        color = PLOT_STYLE["colors"]["immune"][immune_feature]
        ax.scatter(x, y, s=PLOT_STYLE["markers"]["scatter_size"], facecolor=color,
                   edgecolor="white", linewidth=0.7, zorder=3)
        ax.set_xlabel(LABELS["axes"]["oca_tmb"])
        ax.set_ylabel(ylabel)
        if logy:
            ax.set_yscale("log")
        r = association_row(primary, immune_feature, "oca_tmb")
        if PLOT_STYLE["annotations"]["show_stats_on_scatter"]:
            stats = (f"ρ = {r['spearman_rho']:.3f}\n"
                     f"Permutation P = {fmt_p(float(r['permutation_p_two_sided']))}\n"
                     f"BH q = {fmt_p(float(r['bh_q_within_family']))}\n"
                     f"maxT P = {fmt_p(float(r['maxT_fwer_p']))}")
            ax.text(0.03, 0.97, stats, transform=ax.transAxes, ha="left", va="top",
                    fontsize=PLOT_STYLE["font"]["annotation_size"],
                    bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                              edgecolor=PLOT_STYLE["annotations"]["stats_box_edge"],
                              alpha=PLOT_STYLE["annotations"]["stats_box_alpha"], linewidth=0.6))
        apply_axes_style(ax)
        add_panel_label(ax, panel_id)

    axB = fig.add_subplot(gs[0, 1])
    scatter_into(axB, "biopsy_hill_D1", "Fig2B", LABELS["axes"]["biopsy_hill_D1"], logy=True)

    axC = fig.add_subplot(gs[1, 0])
    scatter_into(axC, "biopsy_top10_frequency", "Fig2C", LABELS["axes"]["biopsy_top10_frequency"], percent=True)

    axD = fig.add_subplot(gs[1, 1])
    full = association_row(primary, "biopsy_hill_D1", "oca_tmb")
    d = loso_long[(loso_long["immune_feature"] == "biopsy_hill_D1") &
                  (loso_long["genomic_feature"] == "oca_tmb")].reset_index(drop=True)
    vals = pd.to_numeric(d["rho"]).to_numpy(float)
    color = PLOT_STYLE["colors"]["immune"]["biopsy_hill_D1"]
    full_rho=float(full["spearman_rho"]); lo=float(full["bootstrap_ci_low"]); hi=float(full["bootstrap_ci_high"])
    axD.plot([lo,hi],[2,2],color=color,linewidth=2.2,alpha=.85,zorder=2)
    axD.scatter([full_rho],[2],s=PLOT_STYLE["markers"]["full_size"],marker=PLOT_STYLE["markers"]["full"],
                color=PLOT_STYLE["colors"]["full_effect"],edgecolor="white",linewidth=.7,zorder=4)
    jitter=np.linspace(-.20,.20,len(vals))
    axD.plot([float(np.min(vals)),float(np.max(vals))],[1,1],color=color,linewidth=2.2,alpha=.55,zorder=1)
    axD.scatter(vals,1+jitter,s=PLOT_STYLE["markers"]["loso_size"],color=color,edgecolor="white",linewidth=.6,zorder=3)
    axD.axvline(0,color=PLOT_STYLE["colors"]["zero"],linewidth=.8,linestyle="--",zorder=0)
    axD.set_xlim(-1.02,.08); axD.set_ylim(.55,2.45); axD.set_xlabel("Spearman ρ")
    axD.set_yticks([1,2]); axD.set_yticklabels(["LOSO (13 runs)","Full cohort"])
    apply_axes_style(axD); add_panel_label(axD,"Fig2D")

    save_figure(fig, outstem)


def composite_figureS3(primary, secondary, pole, outstem):
    fig = plt.figure(figsize=FIGURE_SIZE["FigureS5_composite"])
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.45, 1.0], height_ratios=[1.0, 1.0],
                  wspace=0.44, hspace=0.42)

    # A forest
    ax = fig.add_subplot(gs[:, 0])
    d = primary.copy()
    d["_io"] = d["immune_feature"].map({f:i for i,f in enumerate(IMMUNE_ORDER)})
    d["_go"] = d["genomic_feature"].map({f:i for i,f in enumerate(GENOMIC_ORDER)})
    d = d.sort_values(["_io","_go"]).reset_index(drop=True)
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color=PLOT_STYLE["colors"]["zero"], linestyle="--", linewidth=0.8)
    for k,row in d.iterrows():
        yy=y[k]; color=PLOT_STYLE["colors"]["immune"][row["immune_feature"]]
        ax.plot([row["bootstrap_ci_low"],row["bootstrap_ci_high"]],[yy,yy],color=color,linewidth=1.2,alpha=.82)
        face=color if row["bh_q_within_family"]<.05 else "white"
        ax.scatter([row["spearman_rho"]],[yy],s=34,facecolor=face,edgecolor=color,linewidth=1.1,zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{human_immune(i)} — {human_genomic(g)}" for i,g in zip(d.immune_feature,d.genomic_feature)])
    ax.set_xlim(-1.02,1.02); ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    apply_axes_style(ax); add_panel_label(ax,"FigS5A")

    axB = fig.add_subplot(gs[0,1])
    im = draw_assoc_heatmap(axB, secondary, SECONDARY_ORDER, GENOMIC_ORDER, "FigS5B")
    cb=fig.colorbar(im,ax=axB,fraction=.046,pad=.04); cb.set_label("Spearman ρ")

    axC = fig.add_subplot(gs[1,1])
    d2=pole.copy(); d2["_o"]=d2["immune_feature"].map({f:i for i,f in enumerate(IMMUNE_ORDER)}); d2=d2.sort_values("_o").reset_index(drop=True)
    yy=np.arange(len(d2))[::-1]; axC.axvline(0,color=PLOT_STYLE["colors"]["zero"],linestyle="--",linewidth=.8)
    for k,row in d2.iterrows():
        c=PLOT_STYLE["colors"]["immune"][row["immune_feature"]]
        axC.plot([row["bootstrap_ci_low"],row["bootstrap_ci_high"]],[yy[k],yy[k]],color=c,linewidth=1.2)
        axC.scatter([row["rank_biserial_effect"]],[yy[k]],s=34,color=c,edgecolor="white",linewidth=.7,zorder=3)
    axC.set_yticks(yy); axC.set_yticklabels([human_immune(x) for x in d2["immune_feature"]])
    axC.set_xlim(-1.02,1.02); axC.set_xlabel("Rank-biserial effect (bootstrap 95% CI)")
    apply_axes_style(axC); add_panel_label(axC,"FigS5C")

    save_figure(fig,outstem)


# =============================================================================
# Step-3 / Step-4 audit figures and tables
# =============================================================================

def human_candidate(feature: str) -> str:
    if feature in LABELS["immune"]:
        return human_immune(feature)
    replacements = {
        "biopsy_hill_D2": "Hill D2",
        "biopsy_top100_frequency": "Top-100 frequency",
        "biopsy_gini_frequency": "Gini frequency",
        "biopsy_fraction_clones_to_50pct": "Fraction of clonotypes to 50% reads",
        "within_top_1pct_pbmc_detection_fraction_std_median": "Top-1% PBMC recovery",
        "within_top_10pct_pbmc_detection_fraction_std_median": "Top-10% PBMC recovery",
        "within_biopsy_clone_fraction_shared_std_median": "Global clone recovery",
        "within_biopsy_abundance_fraction_shared_std_median": "Shared biopsy read abundance",
        "pbmc_hill_D1": "PBMC Hill D1",
        "biopsy_public_raw_public_clone_fraction_ge2subjects": "Public clonotype fraction",
        "biopsy_public_raw_public_abundance_fraction_ge2subjects": "Public read-abundance fraction",
    }
    return replacements.get(feature, feature.replace("_", " "))


def immune_candidate_representativeness_plot(audit: pd.DataFrame, outstem: Path):
    d = audit[audit["panel_type"].astype(str).eq("core")].copy()
    d = d.sort_values(["domain_order", "selection_rank", "feature"], na_position="last").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS3A"])
    y = np.arange(len(d))[::-1]
    selected = d["selected_full"].astype(bool).to_numpy()
    x = pd.to_numeric(d["representativeness"], errors="coerce").fillna(0).to_numpy(float)
    colors = [PLOT_STYLE["colors"]["audit_selected"] if z else PLOT_STYLE["colors"]["audit_unselected"] for z in selected]
    sizes = [54 if z else 28 for z in selected]
    ax.scatter(x, y, s=sizes, c=colors, edgecolor="white", linewidth=0.6, zorder=3)
    # Domain separators.
    domains = d["domain_order"].to_numpy()
    for i in range(1, len(d)):
        if domains[i] != domains[i-1]:
            ax.axhline(len(d)-i-0.5, color=PLOT_STYLE["colors"]["grid"], linewidth=0.7, zorder=0)
    labels = []
    for r in d.itertuples(index=False):
        star = "● " if bool(r.selected_full) else "   "
        labels.append(star + human_candidate(str(r.feature)))
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(-0.03, 1.04)
    ax.set_xlabel("Within-domain representativeness (median |Spearman ρ|)")
    apply_axes_style(ax, grid=False)
    add_panel_label(ax, "FigS3A")
    fig.tight_layout()
    save_figure(fig, outstem)


def immune_loso_stability_plot(stability: pd.DataFrame, rec: pd.DataFrame, outstem: Path):
    d = stability[stability["panel_type"].astype(str).eq("core")].copy()
    order_map = dict(zip(rec["domain_id"].astype(str), rec["domain_order"]))
    label_map = dict(zip(rec["domain_id"].astype(str), rec["domain_label"].astype(str)))
    d["_order"] = d["domain_id"].astype(str).map(order_map)
    d = d.sort_values("_order", ascending=False).reset_index(drop=True)
    y = np.arange(len(d))
    vals = 100 * pd.to_numeric(d["most_selected_fraction_all_runs"], errors="coerce").to_numpy(float)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS3B"])
    ax.barh(y, vals, color=PLOT_STYLE["colors"]["audit_selected"], alpha=0.88, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([label_map.get(str(x), str(x)) for x in d["domain_id"]])
    ax.set_xlim(0, 104)
    ax.set_xlabel("LOSO runs selecting full-cohort representative (%)")
    for yy, val, cnt, nrun in zip(y, vals, d["most_selected_count"], d["n_loso_runs"]):
        ax.text(min(val + 1.2, 101), yy, f"{int(cnt)}/{int(nrun)}", va="center", ha="left", fontsize=PLOT_STYLE["font"]["annotation_size"])
    apply_axes_style(ax, grid=False)
    add_panel_label(ax, "FigS3B")
    fig.tight_layout()
    save_figure(fig, outstem)


def immune_selected_corr_heatmap(corr: pd.DataFrame, outstem: Path):
    d = corr.copy()
    first = d.columns[0]
    if str(first).startswith("Unnamed") or first not in d.columns[1:]:
        d = d.set_index(first)
    else:
        d = d.set_index(first)
    d.index = d.index.astype(str)
    d.columns = d.columns.astype(str)
    d = d.reindex(index=IMMUNE_ORDER, columns=IMMUNE_ORDER)
    vals = d.apply(pd.to_numeric, errors="coerce").to_numpy(float)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS3C"])
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax.imshow(vals, cmap=PLOT_STYLE["colors"]["heatmap_cmap"], norm=norm)
    ax.set_xticks(range(len(IMMUNE_ORDER))); ax.set_yticks(range(len(IMMUNE_ORDER)))
    ax.set_xticklabels([human_immune(x) for x in IMMUNE_ORDER], rotation=35, ha="right")
    ax.set_yticklabels([human_immune(x) for x in IMMUNE_ORDER])
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v=vals[i,j]
            ax.text(j,i,f"{v:.2f}",ha="center",va="center",color="white" if abs(v)>=.55 else PLOT_STYLE["colors"]["axis"],fontsize=PLOT_STYLE["font"]["annotation_size"])
    cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04); cb.set_label("Spearman ρ")
    add_panel_label(ax,"FigS3C")
    fig.tight_layout()
    save_figure(fig,outstem)


def recovery_denominator_support_plot(support: pd.DataFrame, outstem: Path):
    specs=[("Top 1%","within_raw_top_1pct_n_clonotypes",False),("Top 5%","within_raw_top_5pct_n_clonotypes",True),("Top 10%","within_raw_top_10pct_n_clonotypes",False)]
    rows=[]
    for label,col,selected in specs:
        x=pd.to_numeric(support[col],errors="coerce").dropna().to_numpy(float)
        rows.append((label,float(np.min(x)),float(np.median(x)),float(np.max(x)),selected))
    fig,ax=plt.subplots(figsize=FIGURE_SIZE["FigS3D"])
    y=np.arange(len(rows))[::-1]
    for yy,(label,mn,med,mx,selected) in zip(y,rows):
        c=PLOT_STYLE["colors"]["audit_selected"] if selected else PLOT_STYLE["colors"]["audit_unselected"]
        ax.plot([mn,mx],[yy,yy],color=c,linewidth=3 if selected else 2,alpha=.9)
        ax.scatter([med],[yy],s=56 if selected else 38,color=c,edgecolor="white",linewidth=.7,zorder=3)
        ax.scatter([mn,mx],[yy,yy],s=18,color=c,zorder=3)
        ax.text(mx*1.06,yy,f"min {mn:g}; median {med:g}; max {mx:g}",va="center",fontsize=PLOT_STYLE["font"]["annotation_size"])
    ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel("Biopsy clonotypes supporting abundance-ranked recovery estimand")
    apply_axes_style(ax,grid=False); add_panel_label(ax,"FigS3D")
    fig.tight_layout(); save_figure(fig,outstem)


def composite_figureS3_immune_audit(inputs: Dict[str, object], outstem: Path):
    # Recreate the four audit panels in a single canvas using compact internal drawing.
    fig=plt.figure(figsize=FIGURE_SIZE["FigureS3_composite"])
    gs=GridSpec(2,2,figure=fig,wspace=.55,hspace=.45)
    # A
    ax=fig.add_subplot(gs[0,0]); d=inputs["immune_candidate_audit"]; d=d[d["panel_type"].astype(str).eq("core")].copy(); d=d.sort_values(["domain_order","selection_rank","feature"],na_position="last").reset_index(drop=True); yy=np.arange(len(d))[::-1]; sel=d["selected_full"].astype(bool).to_numpy(); xx=pd.to_numeric(d["representativeness"],errors="coerce").fillna(0).to_numpy(float); ax.scatter(xx,yy,s=[50 if z else 25 for z in sel],c=[PLOT_STYLE["colors"]["audit_selected"] if z else PLOT_STYLE["colors"]["audit_unselected"] for z in sel],edgecolor="white",linewidth=.5); ax.set_yticks(yy); ax.set_yticklabels([("● " if bool(r.selected_full) else "   ")+human_candidate(str(r.feature)) for r in d.itertuples(index=False)],fontsize=6.7); ax.set_xlim(-.03,1.04); ax.set_xlabel("Representativeness (median |ρ|)"); apply_axes_style(ax); add_panel_label(ax,"FigS3A")
    # B
    ax=fig.add_subplot(gs[0,1]); st=inputs["immune_stability"]; rec=inputs["immune_rec"]; st=st[st["panel_type"].astype(str).eq("core")].copy(); om=dict(zip(rec.domain_id.astype(str),rec.domain_order)); lm=dict(zip(rec.domain_id.astype(str),rec.domain_label.astype(str))); st["_o"]=st.domain_id.astype(str).map(om); st=st.sort_values("_o",ascending=False); yy=np.arange(len(st)); vv=100*pd.to_numeric(st.most_selected_fraction_all_runs).to_numpy(float); ax.barh(yy,vv,color=PLOT_STYLE["colors"]["audit_selected"],height=.62); ax.set_yticks(yy); ax.set_yticklabels([lm.get(str(x),str(x)) for x in st.domain_id],fontsize=6.7); ax.set_xlim(0,104); ax.set_xlabel("LOSO selection stability (%)"); [ax.text(min(v+1,101),y,f"{int(c)}/{int(n)}",va="center",fontsize=6.7) for y,v,c,n in zip(yy,vv,st.most_selected_count,st.n_loso_runs)]; apply_axes_style(ax); add_panel_label(ax,"FigS3B")
    # C
    ax=fig.add_subplot(gs[1,0]); cc=inputs["immune_selected_corr"].copy(); cc=cc.set_index(cc.columns[0]); cc=cc.reindex(index=IMMUNE_ORDER,columns=IMMUNE_ORDER); vals=cc.apply(pd.to_numeric,errors="coerce").to_numpy(float); im=ax.imshow(vals,cmap=PLOT_STYLE["colors"]["heatmap_cmap"],norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)); ax.set_xticks(range(4)); ax.set_yticks(range(4)); ax.set_xticklabels([human_immune(x) for x in IMMUNE_ORDER],rotation=35,ha="right",fontsize=6.5); ax.set_yticklabels([human_immune(x) for x in IMMUNE_ORDER],fontsize=6.5); [[ax.text(j,i,f"{vals[i,j]:.2f}",ha="center",va="center",color="white" if abs(vals[i,j])>=.55 else PLOT_STYLE["colors"]["axis"],fontsize=6.8) for j in range(4)] for i in range(4)]; cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04); cb.set_label("Spearman ρ"); add_panel_label(ax,"FigS3C")
    # D
    ax=fig.add_subplot(gs[1,1]); sp=inputs["immune_support"]; specs=[("Top 1%","within_raw_top_1pct_n_clonotypes",False),("Top 5%","within_raw_top_5pct_n_clonotypes",True),("Top 10%","within_raw_top_10pct_n_clonotypes",False)]; rr=[];
    for lab,col,selected in specs:
        x=pd.to_numeric(sp[col],errors="coerce").dropna().to_numpy(float); rr.append((lab,float(x.min()),float(np.median(x)),float(x.max()),selected))
    yy=np.arange(3)[::-1]
    for y,(lab,mn,med,mx,selected) in zip(yy,rr):
        c=PLOT_STYLE["colors"]["audit_selected"] if selected else PLOT_STYLE["colors"]["audit_unselected"]; ax.plot([mn,mx],[y,y],color=c,linewidth=3 if selected else 2); ax.scatter([med],[y],s=50 if selected else 34,color=c,edgecolor="white",linewidth=.6,zorder=3); ax.scatter([mn,mx],[y,y],s=15,color=c,zorder=3)
    ax.set_xscale("log"); ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rr]); ax.set_xlabel("Supporting biopsy clonotypes (min–median–max)"); apply_axes_style(ax); add_panel_label(ax,"FigS3D")
    save_figure(fig,outstem)


def step3_audit_summary(inputs: Dict[str, object]) -> pd.DataFrame:
    c3=inputs["c3"]; merge=inputs["step3_merge_audit"]; ident=inputs["step3_identity_audit"]; den=inputs["step3_denominator_audit"]; oca=inputs["step3_oca_audit"]
    rows=[
        ("Curated subjects", int(c3.get("n_curated_subjects", len(merge))), "Clinical/pathology table defines subject universe"),
        ("RepSeq repertoires", int(c3.get("n_repseq_rows", len(ident))), "Paired PBMC and tumor-biopsy repertoires"),
        ("Subjects with complete TCR core", int(merge["complete_tcr_core"].astype(bool).sum()), f"of {len(merge)}"),
        ("Subjects with complete integrated TCR + OCA data", int(merge["complete_integrated"].astype(bool).sum()), f"of {len(merge)}"),
        ("RepSeq–Step 1 filename identity checks passed", int(ident["filename_match"].astype(bool).sum()), f"of {len(ident)} repertoires"),
        ("Step 1–Step 2 denominator audits passed", int(den["all_checks_pass"].astype(bool).sum()), f"of {len(den)} subjects"),
        ("OCA rows mapped to curated cohort", int((oca["action"].astype(str)=="MERGE").sum()), f"of {len(oca)} OCA rows"),
        ("Step-2 matched-depth resamples", int(c3.get("frozen_upstream_contract",{}).get("step2_n_resamples", np.nan)), "per PBMC–biopsy pair"),
    ]
    return pd.DataFrame(rows,columns=["Audit item","Result","Detail"])


def immune_core_audit_table(inputs: Dict[str, object]) -> pd.DataFrame:
    rec=inputs["immune_rec"].copy()
    cand=inputs["immune_candidate_audit"].copy()
    groups=cand[cand["panel_type"].astype(str).eq("core")].groupby("domain_id",sort=False)["feature"].apply(list).to_dict()
    rec["Candidate measures"]=["; ".join(human_candidate(x) for x in groups.get(str(d),[])) for d in rec["domain_id"]]
    rec["Selected representative"]=rec["selected_feature"].map(human_candidate)
    rec["LOSO selection stability"] = rec.apply(lambda r:f"{int(r.most_selected_count)}/{int(r.n_loso_runs)}",axis=1)
    rec["Max |rho| with relevant depth"] = rec["max_abs_depth_rho"]
    rec["Median resampling interval width"] = rec["resampling_median_interval_width"]
    rec["Minimum clonotype denominator"] = rec["raw_denominator_min"]
    rec["Review flags"] = rec["review_flags"].fillna("none").replace("","none")
    return rec[["domain_label","Candidate measures","Selected representative","selection_mode","LOSO selection stability","Max |rho| with relevant depth","Median resampling interval width","Minimum clonotype denominator","Review flags","freeze_status"]].rename(columns={"domain_label":"Immune domain","selection_mode":"Selection mode","freeze_status":"Freeze status"})


def main_table2_frozen_panels(inputs: Dict[str, object]) -> pd.DataFrame:
    immune=immune_core_audit_table(inputs)
    rows=[]
    for r in immune.itertuples(index=False):
        rows.append({"Panel":"Immune","Biological domain":r[0],"Candidate / source-defined measures":r[1],"Frozen representative":r[2],"Selection / audit":r[4]+" LOSO; "+("no review flag" if r[8]=="none" else r[8]),"Status":r[9]})
    gnum=inputs["genomic_numeric"].copy(); grec=inputs["genomic_rec"].copy(); m=gnum.merge(grec[["feature","recommendation","audit_flags"]],on="feature",how="left",suffixes=("","_rec"))
    domain_map=dict(zip(inputs["genomic_registry"]["feature"].astype(str),inputs["genomic_registry"]["label"].astype(str)))
    for _,r in m.iterrows():
        f=str(r["feature"]); flags=str(r.get("audit_flags_rec",r.get("audit_flags","none"))); flags="none" if flags in {"", "nan", "None"} else flags
        rows.append({"Panel":"Genomic","Biological domain":domain_map.get(f,human_genomic(f)),"Candidate / source-defined measures":"Source-defined continuous OCA estimand","Frozen representative":human_genomic(f),"Selection / audit":f"n={int(r['n_nonmissing'])}; unique={int(r['n_unique'])}; audit flag: {flags}","Status":str(r["recommendation"])})
    return pd.DataFrame(rows)


# =============================================================================
# Publication tables / data
# =============================================================================

def compact_primary_table(primary: pd.DataFrame) -> pd.DataFrame:
    d = primary.copy()
    d["Immune descriptor"] = d["immune_feature"].map(human_immune)
    d["Genomic estimand"] = d["genomic_feature"].map(human_genomic)
    d["Spearman rho"] = d["spearman_rho"]
    d["Permutation P"] = d["permutation_p_two_sided"]
    d["BH q"] = d["bh_q_within_family"]
    d["maxT FWER P"] = d["maxT_fwer_p"]
    d["Bootstrap 95% CI"] = d.apply(lambda r: f"{r.bootstrap_ci_low:.3f} to {r.bootstrap_ci_high:.3f}", axis=1)
    d["LOSO rho range"] = d.apply(lambda r: f"{r.loso_min_rho:.3f} to {r.loso_max_rho:.3f}", axis=1)
    d["LOSO sign consistency"] = d["loso_sign_consistency_fraction"].map(lambda x: f"{int(round(13*x))}/13" if np.isfinite(x) else "NA")
    return d[["Immune descriptor","Genomic estimand","Spearman rho","Permutation P","BH q","maxT FWER P","Bootstrap 95% CI","LOSO rho range","LOSO sign consistency"]]


def compact_secondary_table(d: pd.DataFrame) -> pd.DataFrame:
    x=d.copy(); x["Immune descriptor"]=x["immune_feature"].map(human_immune); x["Genomic estimand"]=x["genomic_feature"].map(human_genomic)
    x["Bootstrap 95% CI"]=x.apply(lambda r:f"{r.bootstrap_ci_low:.3f} to {r.bootstrap_ci_high:.3f}",axis=1)
    return x[["Immune descriptor","Genomic estimand","spearman_rho","permutation_p_two_sided","bh_q_within_family","maxT_fwer_p","Bootstrap 95% CI","loso_sign_consistency_fraction"]].rename(columns={
        "spearman_rho":"Spearman rho","permutation_p_two_sided":"Permutation P","bh_q_within_family":"BH q","maxT_fwer_p":"maxT FWER P","loso_sign_consistency_fraction":"LOSO sign consistency fraction"})


def compact_pole_table(pole: pd.DataFrame) -> pd.DataFrame:
    d=pole.copy(); d["Immune descriptor"]=d["immune_feature"].map(human_immune)
    d["Bootstrap 95% CI"]=d.apply(lambda r:f"{r.bootstrap_ci_low:.3f} to {r.bootstrap_ci_high:.3f}",axis=1)
    return d[["Immune descriptor","n_pole_positive","n_pole_negative","rank_biserial_effect","exact_label_permutation_p_two_sided","bh_q_within_exploratory_pole_family","Bootstrap 95% CI"]].rename(columns={
        "n_pole_positive":"POLE positive, n","n_pole_negative":"POLE negative, n","rank_biserial_effect":"Rank-biserial effect","exact_label_permutation_p_two_sided":"Exact permutation P","bh_q_within_exploratory_pole_family":"BH q"})


def genomic_audit_table(numeric: pd.DataFrame, rec: pd.DataFrame) -> pd.DataFrame:
    d=numeric.merge(rec[["feature","recommendation","audit_flags"]],on="feature",how="left",suffixes=("","_rec"))
    d["Genomic estimand"]=d["feature"].map(human_genomic)
    d["Median [IQR]"]=d.apply(lambda r:f"{r['median']:.2f} [{r['q1']:.2f}–{r['q3']:.2f}]",axis=1)
    d["Range"]=d.apply(lambda r:f"{r['min']:.2f}–{r['max']:.2f}",axis=1)
    d["Zero fraction (%)"]=100*d["zero_fraction"]
    flags=d["audit_flags_rec"].fillna(d.get("audit_flags",pd.Series(index=d.index,dtype=str))).fillna("none")
    d["Audit flags"]=flags
    return d[["Genomic estimand","n_nonmissing","n_unique","Median [IQR]","Range","Zero fraction (%)","recommendation","Audit flags"]].rename(columns={"n_nonmissing":"n","n_unique":"Unique values","recommendation":"Recommendation"})


def table_to_markdown(df: pd.DataFrame) -> str:
    cols=list(df.columns); lines=["| "+" | ".join(cols)+" |","| "+" | ".join(["---"]*len(cols))+" |"]
    for _,r in df.iterrows():
        vals=[]
        for c in cols:
            v=r[c]
            if isinstance(v,float):
                vals.append("NA" if not np.isfinite(v) else f"{v:.4g}")
            else:
                vals.append(str(v).replace("|","\\|"))
        lines.append("| "+" | ".join(vals)+" |")
    return "\n".join(lines)+"\n"


def write_supplementary_data(data_dir: Path, inputs: Dict[str, object], subject_map: Optional[pd.DataFrame]) -> Path:
    ensure_dir(data_dir)
    files = {
        # Step 3: publication-relevant cohort/cross-step audit material.
        "01_step3_cohort_integration_audit.csv": inputs["step3_merge_audit"],
        "02_step3_step1_step2_denominator_audit.csv": inputs["step3_denominator_audit"],
        "03_step3_column_dictionary.csv": inputs["step3_dictionary"],
        # Step 4: complete outcome-blind immune audit/freeze.
        "04_step4_domain_registry.csv": inputs["immune_registry"],
        "05_step4_candidate_audit.csv": inputs["immune_candidate_audit"],
        "06_step4_depth_associations.csv": inputs["immune_depth"],
        "07_step4_within_domain_correlations.csv": inputs["immune_within_corr"],
        "08_step4_domain_selection_full.csv": inputs["immune_selection"],
        "09_step4_domain_selection_loso.csv": inputs["immune_loso"],
        "10_step4_domain_selection_stability.csv": inputs["immune_stability"],
        "11_step4_selected_representative_correlations.csv": inputs["immune_selected_corr"],
        "12_step4_cross_domain_redundancy.csv": inputs["immune_redundancy"],
        "13_step4_immune_panel_recommendation.csv": inputs["immune_rec"],
        "14_step4_immune_panel_values.csv": inputs["immune_values"],
        "15_step4_secondary_domain_audit.csv": inputs["immune_secondary"],
        "16_step4_recovery_denominator_support.csv": inputs["immune_support"],
        # Step 5 genomic audit/freeze.
        "17_step5_genomic_numeric_audit.csv": inputs["genomic_numeric"],
        "18_step5_genomic_categorical_audit.csv": inputs["genomic_categorical"],
        "19_step5_genomic_redundancy_pairs.csv": inputs["genomic_redundancy"],
        "20_step5_genomic_panel_recommendation.csv": inputs["genomic_rec"],
        "21_step5_genomic_primary_values.csv": inputs["genomic_primary"],
        "22_step5_genomic_secondary_values.csv": inputs["genomic_secondary"],
        # Step 6 fixed association families.
        "23_step6_subject_analysis_matrix.csv": inputs["matrix"],
        "24_step6_primary_associations.csv": inputs["primary"],
        "25_step6_primary_loso_summary.csv": inputs["primary_loso_summary"],
        "26_step6_primary_loso_long.csv": inputs["primary_loso_long"],
        "27_step6_secondary_recovery_associations.csv": inputs["secondary"],
        "28_step6_secondary_recovery_loso_summary.csv": inputs["secondary_loso_summary"],
        "29_step6_secondary_recovery_loso_long.csv": inputs["secondary_loso_long"],
        "30_step6_exploratory_pole_associations.csv": inputs["pole"],
        "31_step6_frozen_panel_contract.csv": inputs["contract"],
        "32_step6_upstream_provenance_contract.csv": inputs["provenance"],
    }
    source_id_present = subject_map is None
    for name, df in files.items():
        d=df.copy()
        d=map_subject_column(d,subject_map,"subject")
        d=map_subject_column(d,subject_map,"left_out_subject")
        if source_id_present and ("subject" in d.columns or "left_out_subject" in d.columns):
            name=name.replace(".csv","_REVIEW_SOURCE_IDS.csv")
        d.to_csv(data_dir/name,index=False)

    # Aggregated no-ID summaries are always publication safe.
    step3_audit_summary(inputs).to_csv(data_dir/"00_step3_audit_summary.csv",index=False)
    main_table2_frozen_panels(inputs).to_csv(data_dir/"00_frozen_panel_summary.csv",index=False)

    readme=(
        "Supplementary Data 3 — integrated Steps 3–6 audit and immune–genomic analysis\n"
        "===========================================================================\n\n"
        "This directory contains reporting copies of publication-relevant Step-3 audit outputs, "
        "the complete Step-4 immune audit/freeze, Step-5 genomic audit/freeze, and Step-6 association outputs.\n"
        "The reporting script does not recompute feature selection or inferential statistics.\n\n"
    )
    if subject_map is None:
        readme += "WARNING: subject-level files retain source identifiers and are marked REVIEW_SOURCE_IDS. Supply --subject_map before publication.\n"
    else:
        readme += "Subject identifiers were replaced using the user-supplied subject/public_id mapping.\n"
    readme += "Detailed Step-3 filename/OCA source-identifier audits are intentionally not copied into the publication Supplementary Data archive; they remain internal provenance records in Step 3.\n"
    (data_dir/"README.txt").write_text(readme,encoding="utf-8")

    zip_path=data_dir.parent/"Supplementary_Data_3.zip"
    with zipfile.ZipFile(zip_path,"w",compression=zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(data_dir.iterdir()):
            zf.write(p,arcname=p.name)
    return zip_path


def write_captions(path: Path, n: int, primary: pd.DataFrame):
    d1=association_row(primary,"biopsy_hill_D1","oca_tmb")
    t10=association_row(primary,"biopsy_top10_frequency","oca_tmb")
    text=f"""# Draft figure legends

## Figure 2. Frozen immune repertoire descriptors and OCA genomic features in endometrial cancer.
**(A)** Spearman correlation matrix for the prespecified primary family comprising four frozen immune descriptors and four frozen continuous OCA genomic estimands ({n} subjects; 16 tests). Cell values are Spearman rho. An asterisk indicates BH-FDR q < 0.05 within the 16-test primary family; a black cell border indicates max-|rho| permutation FWER P < 0.05. **(B)** Intratumoral Hill D1 versus tumor mutational burden (TMB). Each point represents one subject; no regression model is fitted. The association was rho={d1['spearman_rho']:.3f}, permutation P={fmt_p(float(d1['permutation_p_two_sided']))}, BH q={fmt_p(float(d1['bh_q_within_family']))}, and maxT FWER P={fmt_p(float(d1['maxT_fwer_p']))}. **(C)** Tumor Top-10 cumulative read frequency versus TMB (rho={t10['spearman_rho']:.3f}, permutation P={fmt_p(float(t10['permutation_p_two_sided']))}, BH q={fmt_p(float(t10['bh_q_within_family']))}, maxT FWER P={fmt_p(float(t10['maxT_fwer_p']))}). **(D)** Full-cohort and leave-one-subject-out stability of the Hill D1–TMB Spearman association. The full-cohort marker is accompanied by its subject-level bootstrap 95% interval; individual points show the 13 rho values obtained after omitting one subject. Primary P values are two-sided subject-label permutation P values. BH-FDR and block-preserving max-|rho| permutation FWER were calculated across the 16-test primary family.

## Supplementary Figure S3. Outcome-blind audit and freeze of the primary immune panel.
**(A)** Within-domain representativeness of eligible candidate measurements in the four prespecified core immune domains. Representativeness is the median absolute Spearman correlation with the other eligible candidates in the same domain; filled/dark points identify the selected representative. Single-estimand domains have representativeness fixed at 1 by definition. **(B)** Leave-one-subject-out re-selection stability of the four frozen immune representatives; labels report the number of LOSO runs selecting the full-cohort representative. **(C)** Pairwise Spearman correlation matrix among the four frozen immune representatives. Correlation was audited for redundancy but was not used to select features according to genomic or clinical associations. **(D)** Subject-level clonotype-count support for the Top-1%, Top-5%, and Top-10% abundance-conditioned tumor-to-blood recovery definitions. Lines span the observed minimum to maximum denominator across subjects and points denote medians; Top-5% is the outcome-blind selected primary recovery definition.

## Supplementary Figure S4. Correlation structure of the frozen continuous genomic panel.
**(A)** Pairwise Spearman correlations among TMB, MSI score, GIM, and LOH in the {n}-subject cohort. The genomic panel was frozen outcome-blind; correlation was audited but was not used to select variables based on immune associations.

## Supplementary Figure S5. Complete effect-size and sensitivity summaries for the frozen immune–genomic analysis.
**(A)** All 16 prespecified primary Spearman associations with subject-level bootstrap 95% confidence intervals. Filled points denote BH-FDR q < 0.05. **(B)** Spearman correlation matrix for the separate 12-test recovery-sensitivity family comprising global tumor-to-blood clone recovery, Top-10% PBMC recovery, and shared biopsy read abundance against TMB, MSI score, GIM, and LOH. **(C)** Exploratory associations between POLE source status and the four primary immune descriptors. Points show rank-biserial effects and horizontal lines show stratified-bootstrap 95% intervals. Exact P values enumerate all binary-label assignments with the observed POLE group sizes and BH-FDR is applied within the four-test exploratory family.
"""
    path.write_text(text,encoding="utf-8")


# =============================================================================
# Main
# =============================================================================

def main():
    ap=argparse.ArgumentParser(description="Integrated publication reporting for frozen Steps 3 v3, 4 v3, 5 v1.1 and 6 v2; no new inference.")
    ap.add_argument("--step3",type=Path,required=True,help="Step-3 v3 directory or ZIP")
    ap.add_argument("--step4",type=Path,required=True,help="Step-4 v3 directory or ZIP")
    ap.add_argument("--step5",type=Path,required=True,help="Step-5 v1.1 directory or ZIP")
    ap.add_argument("--step6",type=Path,required=True,help="Step-6 v2 directory or ZIP")
    ap.add_argument("--out_dir",type=Path,required=True)
    ap.add_argument("--subject_map",type=Path,help="Optional mapping table with columns subject, public_id")
    ap.add_argument("--show_subject_ids",action="store_true",help="Show subject/public IDs on scatter and LOSO panels")
    ap.add_argument("--allow_font_fallback",action="store_true")
    ap.add_argument("--overwrite",action="store_true")
    ap.add_argument("--dpi",type=int,default=PLOT_STYLE["save"]["dpi"])
    args=ap.parse_args()

    for p in [args.step3,args.step4,args.step5,args.step6,args.subject_map]:
        if p is not None and not p.exists():
            raise FileNotFoundError(p)
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite")
    if args.dpi<72:
        raise ValueError("--dpi must be >=72")
    PLOT_STYLE["save"]["dpi"]=args.dpi
    PLOT_STYLE["annotations"]["show_subject_ids"]=bool(args.show_subject_ids)
    font_available=configure_font(strict=not args.allow_font_fallback)

    inputs=load_inputs(args.step3,args.step4,args.step5,args.step6)
    subject_map=read_subject_map(args.subject_map)
    if subject_map is not None:
        source_subjects=set(inputs["matrix"]["subject"].astype(str))
        if set(subject_map["subject"]) != source_subjects:
            raise ValueError("subject_map must contain exactly the Step-6 subject set")

    out=ensure_dir(args.out_dir)
    main_fig=ensure_dir(out/"main_figures")
    supp_fig=ensure_dir(out/"supplementary_figures")
    main_tab=ensure_dir(out/"main_tables")
    supp_tab=ensure_dir(out/"supplementary_tables")
    supp_data=ensure_dir(out/"Supplementary_Data_3")
    support=ensure_dir(out/"manuscript_support")

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    make_primary_heatmap(inputs["primary"],main_fig/"Figure2A_primary_4x4_heatmap")
    scatter_panel(inputs["matrix"],inputs["primary"],"biopsy_hill_D1","oca_tmb","Fig2B",
                  main_fig/"Figure2B_HillD1_vs_TMB",LABELS["axes"]["biopsy_hill_D1"],log_y=True,public_map=subject_map)
    scatter_panel(inputs["matrix"],inputs["primary"],"biopsy_top10_frequency","oca_tmb","Fig2C",
                  main_fig/"Figure2C_Top10_vs_TMB",LABELS["axes"]["biopsy_top10_frequency"],percent_y=True,public_map=subject_map)
    loso_stability_panel(inputs["primary"],inputs["primary_loso_long"],main_fig/"Figure2D_HillD1_TMB_LOSO")
    composite_figure2(inputs["primary"],inputs["matrix"],inputs["primary_loso_long"],main_fig/"Figure2_composite",public_map=subject_map)

    immune_candidate_representativeness_plot(inputs["immune_candidate_audit"],supp_fig/"FigureS3A_immune_candidate_representativeness")
    immune_loso_stability_plot(inputs["immune_stability"],inputs["immune_rec"],supp_fig/"FigureS3B_immune_LOSO_selection_stability")
    immune_selected_corr_heatmap(inputs["immune_selected_corr"],supp_fig/"FigureS3C_frozen_immune_correlation_heatmap")
    recovery_denominator_support_plot(inputs["immune_support"],supp_fig/"FigureS3D_recovery_denominator_support")
    composite_figureS3_immune_audit(inputs,supp_fig/"FigureS3_composite")

    genomic_corr_heatmap(inputs["genomic_corr"],supp_fig/"FigureS4A_genomic_correlation_heatmap")

    forest_primary(inputs["primary"],supp_fig/"FigureS5A_primary_forest")
    make_secondary_heatmap(inputs["secondary"],supp_fig/"FigureS5B_secondary_recovery_heatmap")
    pole_forest(inputs["pole"],supp_fig/"FigureS5C_POLE_effects")
    composite_figureS3(inputs["primary"],inputs["secondary"],inputs["pole"],supp_fig/"FigureS5_composite")

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------
    table2=main_table2_frozen_panels(inputs)
    table2.to_csv(main_tab/"Table2_frozen_immune_genomic_panels.csv",index=False)
    (main_tab/"Table2_frozen_immune_genomic_panels.md").write_text(table_to_markdown(table2),encoding="utf-8")

    s4=step3_audit_summary(inputs)
    s4.to_csv(supp_tab/"TableS4_step3_cohort_integration_audit.csv",index=False)
    (supp_tab/"TableS4_step3_cohort_integration_audit.md").write_text(table_to_markdown(s4),encoding="utf-8")

    s5=immune_core_audit_table(inputs)
    s5.to_csv(supp_tab/"TableS5_immune_feature_panel_audit.csv",index=False)
    (supp_tab/"TableS5_immune_feature_panel_audit.md").write_text(table_to_markdown(s5),encoding="utf-8")

    s6=genomic_audit_table(inputs["genomic_numeric"],inputs["genomic_rec"])
    s6.to_csv(supp_tab/"TableS6_genomic_panel_audit.csv",index=False)

    s7=compact_primary_table(inputs["primary"])
    s7.to_csv(supp_tab/"TableS7_primary_immune_genomic_associations.csv",index=False)
    (supp_tab/"TableS7_primary_immune_genomic_associations.md").write_text(table_to_markdown(s7),encoding="utf-8")

    s8=compact_secondary_table(inputs["secondary"])
    s8.to_csv(supp_tab/"TableS8_secondary_recovery_sensitivity.csv",index=False)

    s9=compact_pole_table(inputs["pole"])
    s9.to_csv(supp_tab/"TableS9_exploratory_POLE.csv",index=False)

    disp=inputs["primary"].copy()
    disp["cell"]=disp.apply(lambda r:f"{r.spearman_rho:.2f} (q={fmt_p(float(r.bh_q_within_family))})",axis=1)
    dm=disp.pivot(index="immune_feature",columns="genomic_feature",values="cell").reindex(index=IMMUNE_ORDER,columns=GENOMIC_ORDER)
    dm.index=[human_immune(x) for x in dm.index]; dm.columns=[human_genomic(x) for x in dm.columns]
    dm.to_csv(support/"Primary_4x4_display_matrix_rho_BHq.csv")

    # ------------------------------------------------------------------
    # Supplementary Data 3
    # ------------------------------------------------------------------
    zip_path=write_supplementary_data(supp_data,inputs,subject_map)

    # ------------------------------------------------------------------
    # Manuscript support
    # ------------------------------------------------------------------
    d1=association_row(inputs["primary"],"biopsy_hill_D1","oca_tmb")
    t10=association_row(inputs["primary"],"biopsy_top10_frequency","oca_tmb")
    sig_bh=inputs["primary"][inputs["primary"]["bh_q_within_family"]<0.05]
    sig_maxt=inputs["primary"][inputs["primary"]["maxT_fwer_p"]<0.05]
    numbers={
        "n_subjects":int(inputs["matrix"]["subject"].nunique()),
        "step3_complete_integrated":int(inputs["step3_merge_audit"]["complete_integrated"].astype(bool).sum()),
        "step3_denominator_audit_pass":int(inputs["step3_denominator_audit"]["all_checks_pass"].astype(bool).sum()),
        "immune_core_loso_stability":{str(r.domain_id):float(r.most_selected_fraction_all_runs) for r in inputs["immune_stability"].query("panel_type == 'core'").itertuples(index=False)},
        "n_primary_tests":16,
        "n_primary_bh_q_lt_0_05":int(len(sig_bh)),
        "n_primary_maxT_fwer_p_lt_0_05":int(len(sig_maxt)),
        "hill_D1_vs_TMB":{k:(float(d1[k]) if pd.notna(d1[k]) and isinstance(d1[k],(int,float,np.integer,np.floating)) else d1[k]) for k in ["spearman_rho","permutation_p_two_sided","bh_q_within_family","maxT_fwer_p","bootstrap_ci_low","bootstrap_ci_high","loso_min_rho","loso_max_rho","loso_sign_consistency_fraction"]},
        "top10_vs_TMB":{k:(float(t10[k]) if pd.notna(t10[k]) and isinstance(t10[k],(int,float,np.integer,np.floating)) else t10[k]) for k in ["spearman_rho","permutation_p_two_sided","bh_q_within_family","maxT_fwer_p","bootstrap_ci_low","bootstrap_ci_high","loso_min_rho","loso_max_rho","loso_sign_consistency_fraction"]},
    }
    (support/"manuscript_numbers.json").write_text(json.dumps(numbers,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    write_captions(support/"Figure_legends_draft.md",numbers["n_subjects"],inputs["primary"])

    table_legends=(
        "# Draft table legends\n\n"
        "**Table 2. Outcome-blind audit and freeze of the primary immune and genomic panels.** The four immune domains were prespecified before genomic or histopathological association testing. Competitive immune domains were resolved using TCR-only technical properties and leave-one-subject-out re-selection; source-defined continuous OCA estimands were audited independently of immune associations. Review flags are reported without outcome-driven deletion.\n\n"
        "**Supplementary Table S4. Step-3 cohort integration and cross-step audit.** Aggregate completeness and consistency checks linking the curated cohort, paired RepSeq inputs, Steps 1 and 2, and OCA Plus data.\n\n"
        "**Supplementary Table S5. Outcome-blind immune feature-panel audit.** Core biological domains, eligible candidate measurements, selected representatives, LOSO selection stability, depth-association audit, matched-depth resampling precision where applicable, denominator support, review flags, and final freeze status.\n\n"
        "**Supplementary Table S6. Technical audit of the frozen continuous genomic panel.** Availability, value resolution, distribution summaries, zero mass, technical recommendation, and audit flags for TMB, MSI score, GIM, and LOH.\n\n"
        "**Supplementary Table S7. Primary immune–genomic associations.** Spearman effect sizes for the 16 prespecified primary tests, two-sided subject-label permutation P values, BH-FDR q values across the 16-test family, block-preserving max-|rho| permutation FWER P values, subject-level bootstrap 95% intervals, and leave-one-subject-out effect-size stability.\n\n"
        "**Supplementary Table S8. Secondary recovery-sensitivity associations.** Separate 12-test family comprising three prespecified tumor–blood recovery sensitivity descriptors against the four frozen continuous genomic estimands. Multiplicity corrections are calculated within this secondary family only.\n\n"
        "**Supplementary Table S9. Exploratory POLE associations.** Rank-biserial effects between POLE source status and the four primary immune descriptors. Exact P values enumerate all label assignments with the observed group sizes; BH-FDR is applied across the four-test exploratory family.\n"
    )
    (support/"Table_legends_draft.md").write_text(table_legends,encoding="utf-8")

    checks=[
        {"check":"Step3 version v3.0","status":"PASS"},
        {"check":"Step4 version v3.0","status":"PASS"},
        {"check":"Step5 version v1.1","status":"PASS"},
        {"check":"Step6 version v2.0","status":"PASS"},
        {"check":"Step3 complete integrated cohort","status":"PASS"},
        {"check":"Step3 RepSeq/Step1 identity audit","status":"PASS"},
        {"check":"Step3 Step1/Step2 denominator audit","status":"PASS"},
        {"check":"Step3 analysis-dataset SHA present in Step6 provenance contract","status":"PASS"},
        {"check":"Step4 four core domains frozen","status":"PASS"},
        {"check":"Primary family exactly 4x4 / 16 tests","status":"PASS"},
        {"check":"Secondary recovery family exactly 3x4 / 12 tests","status":"PASS"},
        {"check":"Exploratory POLE family exactly 4 tests","status":"PASS"},
        {"check":"Subject analysis matrix complete for frozen features","status":"PASS"},
        {"check":"No inferential statistic recomputed by reporting layer","status":"PASS"},
    ]
    pd.DataFrame(checks).to_csv(support/"reporting_validation.csv",index=False)

    source_rows=[]
    for role,p in [("step3",args.step3),("step4",args.step4),("step5",args.step5),("step6",args.step6),("subject_map",args.subject_map)]:
        if p is None: continue
        if p.is_file():
            source_rows.append({"source_role":role,"path":str(p.resolve()),"size_bytes":p.stat().st_size,"sha256":sha256_file(p)})
        else:
            source_rows.append({"source_role":role,"path":str(p.resolve()),"size_bytes":np.nan,"sha256":"DIRECTORY_INPUT"})
    pd.DataFrame(source_rows).to_csv(support/"source_sha256_manifest.csv",index=False)

    fig_manifest=[]
    for folder,placement in [(main_fig,"MAIN"),(supp_fig,"SUPPLEMENTARY")]:
        for pp in sorted(folder.glob("*.pdf")):
            fig_manifest.append({"placement":placement,"file":pp.name,"path":str(pp)})
    pd.DataFrame(fig_manifest).to_csv(support/"figure_manifest.csv",index=False)
    table_manifest=[]
    for folder,placement in [(main_tab,"MAIN"),(supp_tab,"SUPPLEMENTARY")]:
        for pp in sorted(folder.glob("*.csv")):
            table_manifest.append({"placement":placement,"file":pp.name,"path":str(pp)})
    pd.DataFrame(table_manifest).to_csv(support/"table_manifest.csv",index=False)

    config={
        "script":Path(__file__).name,
        "script_version":SCRIPT_VERSION,
        "step3":str(args.step3.resolve()),"step4":str(args.step4.resolve()),"step5":str(args.step5.resolve()),"step6":str(args.step6.resolve()),"out_dir":str(out.resolve()),
        "step3_script_version":inputs["c3"].get("script_version"),"step4_script_version":inputs["c4"].get("script_version"),"step5_script_version":inputs["c5"].get("script_version"),"step6_script_version":inputs["c6"].get("script_version"),
        "n_subjects":numbers["n_subjects"],"subject_map_supplied":subject_map is not None,
        "source_ids_shown_on_figures":bool(args.show_subject_ids),"requested_font":PLOT_STYLE["font"]["family"],"font_available":font_available,
        "plot_style":PLOT_STYLE,"figure_size":FIGURE_SIZE,"panel_style":PANEL_STYLE,
        "reporting_contract":"Consumes frozen Steps 3 v3, 4 v3, 5 v1.1 and 6 v2. No upstream audit, feature-selection, permutation, bootstrap, LOSO, or association statistic is recomputed."
    }
    (support/"00_run_config.json").write_text(json.dumps(config,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

    print(f"[OK] Reporting {SCRIPT_VERSION}")
    print(f"[OK] Subjects: {numbers['n_subjects']}; Step-3 integrated complete: {numbers['step3_complete_integrated']}/{numbers['n_subjects']}")
    print(f"[OK] Step-4 core domains: 4; primary tests: 16; secondary tests: 12; POLE tests: 4")
    print(f"[OK] BH q<0.05 primary associations: {numbers['n_primary_bh_q_lt_0_05']}; maxT FWER P<0.05: {numbers['n_primary_maxT_fwer_p_lt_0_05']}")
    print(f"[OK] Supplementary Data ZIP: {zip_path}")
    if subject_map is None:
        print("[WARNING] No subject map supplied: subject-level Supplementary Data files retain source IDs and are marked REVIEW_SOURCE_IDS.")
    print(f"[OUTPUT] {out.resolve()}")


if __name__ == "__main__":
    main()
