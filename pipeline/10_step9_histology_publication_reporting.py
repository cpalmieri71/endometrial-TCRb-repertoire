#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STEP 9 — Publication reporting for frozen immune/genomic panels vs histology
=========================================================================== 

Reporting only. NO NEW INFERENTIAL ANALYSIS is performed.

Input
-----
--step8 : Step-8 v3 output directory or ZIP archive.

Optional
--------
--subject_map : CSV with columns subject,public_id. If supplied, a subject-level
                publication matrix is exported using public_id. Without a map,
                source subject IDs are never included in publication-ready tables.

Outputs
-------
main_tables/
    Supplementary_Table_S10_histology_endpoint_multiplicity_summary.csv/.md
    Supplementary_Table_S11_all_histology_associations.csv/.md
    Supplementary_Table_S12_histotype_p53_descriptive.csv/.md

supplementary_figures/
    FigureS6_composite.pdf/.png
    FigureS6A_effect_heatmap.pdf/.png
    FigureS6B_grade_morisita.pdf/.png
    FigureS6C_invasion_loh.pdf/.png
    FigureS6D_multiplicity_summary.pdf/.png
    FigureS7_forest_all_associations.pdf/.png

publication_data/
    machine-readable publication-ready CSVs

Supplementary_Data_4.zip
00_reporting_config.json
01_reporting_manifest.csv
02_reporting_validation.csv

Scientific contract
-------------------
* Reads effect estimates, bootstrap intervals, permutation P values, BH q values,
  max-|effect| FWER P values and LOSO summaries from Step 8.
* Does not recompute any inferential statistic.
* Figure S6B/C highlight the smallest nominal-P association in grade and invasion
  only as descriptive visualization; inferential interpretation remains based on
  the prespecified corrected results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

SCRIPT_VERSION = "v1.0.0-step8-v3-publication-reporting-2026-10-08"

# =====================================================================
# USER-EDITABLE AESTHETICS
# =====================================================================
PLOT_STYLE = {
    "font_family": "Arial",
    "font_fallback": "DejaVu Sans",
    "font_size": 9.0,
    "axis_label_size": 10.0,
    "tick_label_size": 8.5,
    "title_size": 10.5,
    "annotation_size": 8.0,
    "panel_letter_size": 12.0,
    "line_width": 1.1,
    "axis_line_width": 0.8,
    "marker_size": 46,
    "marker_edge_width": 0.8,
    "ci_line_width": 1.25,
    "ci_capsize": 2.5,
    "grid_alpha": 0.18,
    "dpi": 400,
    # Colors can be edited freely.
    "immune_color": "#2F6B8A",
    "genomic_color": "#B55D45",
    "neutral_color": "#5C5C5C",
    "light_neutral": "#D9D9D9",
    "reference_color": "#222222",
    "heatmap_cmap": "RdBu_r",
    "highlight_grade_color": "#3B7F8C",
    "highlight_invasion_color": "#B66A3C",
    "group0_color": "#A8B3BC",
    "group1_color": "#536B78",
    "box_facecolor": "white",
    "box_edgecolor": "#B8B8B8",
    "box_alpha": 0.92,
}

FIGURE_SIZE = {
    "FigureS6_composite": (13.2, 8.6),
    "FigureS6A_effect_heatmap": (10.0, 4.5),
    "FigureS6B_grade_morisita": (5.4, 4.8),
    "FigureS6C_invasion_loh": (5.4, 4.8),
    "FigureS6D_multiplicity_summary": (6.3, 4.8),
    "FigureS7_forest_all_associations": (10.2, 11.8),
}

# Panel-specific settings. Positions are axes fractions, so (0,0) is bottom-left
# and (1,1) is top-right. Change these values only if you want to move a box.
PANEL_STYLE = {
    "show_panel_letters": True,
    "show_titles": False,
    "show_stats_boxes": True,
    "FigureS6B": {
        "stats_box_x": 0.03,
        "stats_box_y": 0.97,
        "stats_box_ha": "left",
        "stats_box_va": "top",
    },
    "FigureS6C": {
        "stats_box_x": 0.03,
        "stats_box_y": 0.97,
        "stats_box_ha": "left",
        "stats_box_va": "top",
    },
    "heatmap_show_effect": True,
    "heatmap_show_corrected_marker": True,
    "forest_show_global_q": False,
}

LABELS = {
    "biopsy_hill_D1": "Hill D1",
    "biopsy_top10_frequency": "Top-10",
    "within_top_5pct_pbmc_detection_fraction_std_median": "Top-5% recovery",
    "within_morisita_horn_std_median": "Morisita–Horn",
    "oca_tmb": "TMB",
    "oca_msi_score": "MSI score",
    "oca_gim": "GIM",
    "oca_loh": "LOH",
    "grade": "Endometrioid grade",
    "myometrial_invasion": "Myometrial invasion",
    "lvsi": "LVSI",
    "IMMUNE": "Immune",
    "GENOMIC": "Genomic",
}

DESCRIPTOR_ORDER = [
    "biopsy_hill_D1",
    "biopsy_top10_frequency",
    "within_top_5pct_pbmc_detection_fraction_std_median",
    "within_morisita_horn_std_median",
    "oca_tmb",
    "oca_msi_score",
    "oca_gim",
    "oca_loh",
]
ENDPOINT_ORDER = ["grade", "myometrial_invasion", "lvsi"]

# =====================================================================
# Helpers
# =====================================================================

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_font_family() -> str:
    try:
        from matplotlib import font_manager
        wanted = PLOT_STYLE["font_family"]
        names = {f.name for f in font_manager.fontManager.ttflist}
        return wanted if wanted in names else PLOT_STYLE["font_fallback"]
    except Exception:
        return PLOT_STYLE["font_fallback"]


def apply_plot_style() -> None:
    family = resolve_font_family()
    matplotlib.rcParams.update({
        "font.family": family,
        "font.size": PLOT_STYLE["font_size"],
        "axes.labelsize": PLOT_STYLE["axis_label_size"],
        "axes.titlesize": PLOT_STYLE["title_size"],
        "xtick.labelsize": PLOT_STYLE["tick_label_size"],
        "ytick.labelsize": PLOT_STYLE["tick_label_size"],
        "axes.linewidth": PLOT_STYLE["axis_line_width"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.transparent": False,
    })


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=PLOT_STYLE["axis_line_width"])


def panel_letter(ax, letter: str):
    if PANEL_STYLE["show_panel_letters"]:
        ax.text(
            -0.12, 1.05, letter,
            transform=ax.transAxes,
            fontsize=PLOT_STYLE["panel_letter_size"],
            fontweight="bold",
            va="top", ha="left",
        )


def save_figure(fig, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=PLOT_STYLE["dpi"], bbox_inches="tight")
    plt.close(fig)


def as_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if pd.isna(x):
        return False
    return str(x).strip().lower() in {"true", "1", "yes", "y"}


def source_root(path: Path) -> Tuple[Path, Optional[tempfile.TemporaryDirectory]]:
    """Return a directory containing Step-8 files; extract ZIPs to a temp dir."""
    if path.is_dir():
        return path, None
    if path.is_file() and path.suffix.lower() == ".zip":
        td = tempfile.TemporaryDirectory(prefix="step8_reporting_")
        with zipfile.ZipFile(path, "r") as zf:
            zf.extractall(td.name)
        root = Path(td.name)
        # Find the unique directory containing 00_run_config.json.
        hits = [p.parent for p in root.rglob("00_run_config.json") if "__MACOSX" not in p.parts]
        if len(hits) != 1:
            td.cleanup()
            raise FileNotFoundError(f"Expected one Step-8 root in ZIP; found {len(hits)}")
        return hits[0], td
    raise FileNotFoundError(path)


def require_files(root: Path, names: Sequence[str]):
    missing = [n for n in names if not (root / n).exists()]
    if missing:
        raise FileNotFoundError("Missing Step-8 files: " + ", ".join(missing))


def read_json(p: Path) -> Dict:
    return json.loads(p.read_text(encoding="utf-8"))


def fmt_p(x) -> str:
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"


def fmt_num(x, digits=3) -> str:
    return "NA" if pd.isna(x) else f"{float(x):.{digits}f}"


def write_markdown_table(df: pd.DataFrame, path: Path):
    try:
        text = df.to_markdown(index=False)
    except Exception:
        # Fallback if tabulate is unavailable.
        cols = list(df.columns)
        rows = ["| " + " | ".join(map(str, cols)) + " |",
                "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, r in df.iterrows():
            rows.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        text = "\n".join(rows)
    path.write_text(text + "\n", encoding="utf-8")


def jitter_positions(n: int, width: float = 0.08) -> np.ndarray:
    if n <= 1:
        return np.zeros(n)
    return np.linspace(-width, width, n)


def add_stats_box(ax, text: str, panel_key: str):
    if not PANEL_STYLE["show_stats_boxes"]:
        return
    cfg = PANEL_STYLE[panel_key]
    ax.text(
        cfg["stats_box_x"], cfg["stats_box_y"], text,
        transform=ax.transAxes,
        ha=cfg["stats_box_ha"], va=cfg["stats_box_va"],
        fontsize=PLOT_STYLE["annotation_size"],
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor=PLOT_STYLE["box_facecolor"],
            edgecolor=PLOT_STYLE["box_edgecolor"],
            alpha=PLOT_STYLE["box_alpha"],
            linewidth=0.7,
        ),
    )

# =====================================================================
# Validation and publication tables
# =====================================================================

def validate_step8(root: Path) -> Dict:
    cfg = read_json(root / "00_run_config.json")
    assoc = pd.read_csv(root / "17_all_histology_associations.csv")
    family = pd.read_csv(root / "20_family_multiplicity_summary.csv")
    contract = pd.read_csv(root / "04_frozen_descriptor_panel_contract.csv")
    matrix = pd.read_csv(root / "05_histology_descriptor_analysis_matrix.csv", dtype={"subject": str})

    checks = {
        "step8_v3": str(cfg.get("script_version", "")).startswith("v3.0.0"),
        "n_subjects_13": int(cfg.get("n_subjects", -1)) == 13,
        "n_descriptors_per_endpoint_8": int(cfg.get("n_descriptors_per_endpoint", -1)) == 8,
        "n_histology_tests_total_24": int(cfg.get("n_histology_tests_total", -1)) == 24,
        "association_rows_24": len(assoc) == 24,
        "contract_rows_8": len(contract) == 8,
        "matrix_subjects_13": matrix["subject"].nunique() == 13,
        "three_endpoint_families": set(assoc["analysis_family"].astype(str)) == set(ENDPOINT_ORDER),
        "four_immune_four_genomic": (
            (contract["panel"].eq("IMMUNE").sum() == 4)
            and (contract["panel"].eq("GENOMIC").sum() == 4)
        ),
        "no_family_bh_positive": int((pd.to_numeric(assoc["family_bh_q"], errors="coerce") < 0.05).sum()) == 0,
        "no_family_maxT_positive": int((pd.to_numeric(assoc["family_maxT_fwer_p"], errors="coerce") < 0.05).sum()) == 0,
    }
    if not all(checks.values()):
        bad = [k for k, v in checks.items() if not v]
        raise ValueError("Step-8 reporting validation failed: " + ", ".join(bad))
    return {"config": cfg, "checks": checks, "assoc": assoc, "family": family,
            "contract": contract, "matrix": matrix}


def endpoint_distribution(endpoint_audit: pd.DataFrame, endpoint: str) -> Tuple[str, int]:
    if endpoint == "grade":
        s = endpoint_audit.loc[endpoint_audit["hist_grade_analysis_status"].eq("evaluable"), "hist_grade"].astype(str)
        vc = s.value_counts()
        parts = [f"{g}={int(vc.get(g, 0))}" for g in ["G1", "G2", "G3"]]
        return ", ".join(parts), int(len(s))
    if endpoint == "myometrial_invasion":
        s = pd.to_numeric(endpoint_audit["hist_invasion_gt50"], errors="coerce").dropna()
        return f">50%={int((s==1).sum())}, ≤50%={int((s==0).sum())}", int(len(s))
    if endpoint == "lvsi":
        s = pd.to_numeric(endpoint_audit["hist_lvsi_positive"], errors="coerce").dropna()
        return f"positive={int((s==1).sum())}, negative={int((s==0).sum())}", int(len(s))
    return "", 0


def make_table_s10(root: Path, assoc: pd.DataFrame, family: pd.DataFrame) -> pd.DataFrame:
    audit = pd.read_csv(root / "02_histology_endpoint_audit.csv")
    rows = []
    for endpoint in ENDPOINT_ORDER:
        fam = family[family["analysis_family"].astype(str).eq(endpoint)].iloc[0]
        a = assoc[assoc["analysis_family"].astype(str).eq(endpoint)].copy()
        a["permutation_p"] = pd.to_numeric(a["permutation_p"], errors="coerce")
        top = a.sort_values("permutation_p").iloc[0]
        dist, n = endpoint_distribution(audit, endpoint)
        rows.append({
            "Endpoint": LABELS[endpoint],
            "n": n,
            "Distribution": dist,
            "Tests": int(fam["n_tests"]),
            "Smallest nominal-P descriptor": LABELS.get(str(top["descriptor_feature"]), str(top["descriptor_label"])),
            "Effect": round(float(top["effect_size"]), 3),
            "Nominal permutation P": round(float(top["permutation_p"]), 6),
            "Minimum family BH q": round(float(fam["min_family_bh_q"]), 6),
            "Minimum family maxT P": round(float(fam["min_family_maxT_fwer_p"]), 6),
            "Minimum global BH q": round(float(fam["min_global_bh_q"]), 6),
        })
    return pd.DataFrame(rows)


def make_table_s11(assoc: pd.DataFrame) -> pd.DataFrame:
    d = assoc.copy()
    d["endpoint_order"] = d["analysis_family"].map({k: i for i, k in enumerate(ENDPOINT_ORDER)})
    d["descriptor_order"] = d["descriptor_feature"].map({k: i for i, k in enumerate(DESCRIPTOR_ORDER)})
    d = d.sort_values(["endpoint_order", "descriptor_order"])
    out = pd.DataFrame({
        "Endpoint": d["analysis_family"].map(LABELS),
        "Panel": d["panel"].map(LABELS),
        "Descriptor": d["descriptor_feature"].map(lambda x: LABELS.get(x, x)),
        "n": d["n"].astype(int),
        "Effect measure": d["effect_type"].map({"spearman_rho": "Spearman ρ", "rank_biserial": "Rank-biserial"}),
        "Effect": d["effect_size"].round(3),
        "Bootstrap 95% CI low": pd.to_numeric(d["bootstrap_ci_low"], errors="coerce").round(3),
        "Bootstrap 95% CI high": pd.to_numeric(d["bootstrap_ci_high"], errors="coerce").round(3),
        "Permutation P": pd.to_numeric(d["permutation_p"], errors="coerce").round(6),
        "Family BH q": pd.to_numeric(d["family_bh_q"], errors="coerce").round(6),
        "Family maxT P": pd.to_numeric(d["family_maxT_fwer_p"], errors="coerce").round(6),
        "Global BH q": pd.to_numeric(d["global_bh_q"], errors="coerce").round(6),
        "LOSO min": pd.to_numeric(d["loso_min_effect"], errors="coerce").round(3),
        "LOSO max": pd.to_numeric(d["loso_max_effect"], errors="coerce").round(3),
        "LOSO sign consistency": pd.to_numeric(d["loso_sign_consistency_fraction"], errors="coerce").round(3),
    })
    return out.reset_index(drop=True)


def make_table_s12(root: Path) -> pd.DataFrame:
    parts = []
    for filename, label in [
        ("15_histotype_descriptive_summary.csv", "Histotype"),
        ("16_p53_descriptive_summary.csv", "p53 pattern"),
    ]:
        d = pd.read_csv(root / filename)
        if d.empty:
            continue
        x = pd.DataFrame({
            "Descriptive variable": label,
            "Level": d["level"].astype(str),
            "Panel": d["panel"].map(LABELS),
            "Descriptor": d["descriptor_feature"].map(lambda z: LABELS.get(z, z)),
            "n": d["n"].astype(int),
            "Median": pd.to_numeric(d["median"], errors="coerce").round(4),
            "Q1": pd.to_numeric(d["q1"], errors="coerce").round(4),
            "Q3": pd.to_numeric(d["q3"], errors="coerce").round(4),
        })
        parts.append(x)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

# =====================================================================
# Plot data and panels
# =====================================================================

def ordered_assoc(assoc: pd.DataFrame) -> pd.DataFrame:
    d = assoc.copy()
    d["endpoint_order"] = d["analysis_family"].map({k: i for i, k in enumerate(ENDPOINT_ORDER)})
    d["descriptor_order"] = d["descriptor_feature"].map({k: i for i, k in enumerate(DESCRIPTOR_ORDER)})
    return d.sort_values(["endpoint_order", "descriptor_order"]).reset_index(drop=True)


def plot_effect_heatmap(ax, assoc: pd.DataFrame, add_letter=True):
    d = ordered_assoc(assoc)
    mat = np.full((3, 8), np.nan)
    qmat = np.full((3, 8), np.nan)
    pmax = np.full((3, 8), np.nan)
    for i, ep in enumerate(ENDPOINT_ORDER):
        for j, feat in enumerate(DESCRIPTOR_ORDER):
            hit = d[(d["analysis_family"] == ep) & (d["descriptor_feature"] == feat)]
            if len(hit) == 1:
                mat[i, j] = float(hit.iloc[0]["effect_size"])
                qmat[i, j] = float(hit.iloc[0]["family_bh_q"])
                pmax[i, j] = float(hit.iloc[0]["family_maxT_fwer_p"])
    im = ax.imshow(mat, aspect="auto", cmap=PLOT_STYLE["heatmap_cmap"], vmin=-1, vmax=1)
    ax.set_xticks(np.arange(8))
    ax.set_xticklabels([LABELS[x] for x in DESCRIPTOR_ORDER], rotation=35, ha="right")
    ax.set_yticks(np.arange(3))
    ax.set_yticklabels([LABELS[x] for x in ENDPOINT_ORDER])
    # Group divider.
    ax.axvline(3.5, color="white", lw=2.2)
    ax.text(1.5, -0.72, "IMMUNE", ha="center", va="center", fontsize=8, fontweight="bold")
    ax.text(5.5, -0.72, "GENOMIC", ha="center", va="center", fontsize=8, fontweight="bold")
    if PANEL_STYLE["heatmap_show_effect"]:
        for i in range(3):
            for j in range(8):
                val = mat[i, j]
                if not np.isfinite(val):
                    continue
                txt = f"{val:.2f}"
                # Corrected markers only. In this dataset none are expected.
                if PANEL_STYLE["heatmap_show_corrected_marker"] and (qmat[i,j] < .05 or pmax[i,j] < .05):
                    txt += "*"
                color = "white" if abs(val) >= 0.48 else "black"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7.3, color=color)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if PANEL_STYLE["show_titles"]:
        ax.set_title("Frozen immune/genomic descriptors versus histopathology")
    if add_letter:
        panel_letter(ax, "A")
    return im


def plot_grade_morisita(ax, matrix: pd.DataFrame, grade_assoc: pd.DataFrame, add_letter=True):
    feat = "within_morisita_horn_std_median"
    d = matrix[np.isfinite(pd.to_numeric(matrix["hist_grade_numeric"], errors="coerce"))].copy()
    d["grade"] = pd.to_numeric(d["hist_grade_numeric"], errors="coerce")
    for grade in [1,2,3]:
        y = pd.to_numeric(d.loc[d["grade"] == grade, feat], errors="coerce").dropna().sort_values()
        jit = jitter_positions(len(y), 0.07)
        ax.scatter(np.full(len(y), grade) + jit, y, s=PLOT_STYLE["marker_size"],
                   color=PLOT_STYLE["highlight_grade_color"], edgecolor="white",
                   linewidth=PLOT_STYLE["marker_edge_width"], zorder=3)
        if len(y):
            med = float(np.median(y))
            ax.plot([grade-0.16, grade+0.16], [med, med], color=PLOT_STYLE["reference_color"], lw=1.5)
    r = grade_assoc[grade_assoc["descriptor_feature"].eq(feat)].iloc[0]
    stats = (
        f"Spearman ρ = {float(r['spearman_rho']):.3f}\n"
        f"exact P = {fmt_p(r['permutation_p'])}\n"
        f"BH q = {fmt_p(r['family_bh_q'])}\n"
        f"maxT P = {fmt_p(r['family_maxT_fwer_p'])}"
    )
    add_stats_box(ax, stats, "FigureS6B")
    ax.set_xticks([1,2,3])
    ax.set_xticklabels(["G1", "G2", "G3"])
    ax.set_xlabel("Endometrioid grade")
    ax.set_ylabel("Morisita–Horn similarity")
    ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
    clean_axes(ax)
    if PANEL_STYLE["show_titles"]:
        ax.set_title("Morisita–Horn similarity across endometrioid grade")
    if add_letter:
        panel_letter(ax, "B")


def plot_invasion_loh(ax, matrix: pd.DataFrame, invasion_assoc: pd.DataFrame, add_letter=True):
    feat = "oca_loh"
    d = matrix[np.isfinite(pd.to_numeric(matrix["hist_invasion_gt50"], errors="coerce"))].copy()
    g = pd.to_numeric(d["hist_invasion_gt50"], errors="coerce")
    for group, x, color in [(0,0,PLOT_STYLE["group0_color"]),(1,1,PLOT_STYLE["highlight_invasion_color"])]:
        y = pd.to_numeric(d.loc[g == group, feat], errors="coerce").dropna().sort_values()
        jit = jitter_positions(len(y), 0.07)
        ax.scatter(np.full(len(y), x) + jit, y, s=PLOT_STYLE["marker_size"],
                   color=color, edgecolor="white", linewidth=PLOT_STYLE["marker_edge_width"], zorder=3)
        if len(y):
            q1, med, q3 = np.quantile(y, [0.25,0.5,0.75])
            ax.plot([x-0.17,x+0.17],[med,med], color=PLOT_STYLE["reference_color"], lw=1.5)
            ax.plot([x,x],[q1,q3], color=PLOT_STYLE["reference_color"], lw=1.0)
    r = invasion_assoc[invasion_assoc["descriptor_feature"].eq(feat)].iloc[0]
    stats = (
        f"Rank-biserial = {float(r['rank_biserial']):.3f}\n"
        f"exact P = {fmt_p(r['permutation_p'])}\n"
        f"BH q = {fmt_p(r['family_bh_q'])}\n"
        f"maxT P = {fmt_p(r['family_maxT_fwer_p'])}"
    )
    add_stats_box(ax, stats, "FigureS6C")
    ax.set_xticks([0,1])
    ax.set_xticklabels(["≤50%", ">50%"])
    ax.set_xlabel("Myometrial invasion")
    ax.set_ylabel("LOH")
    ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
    clean_axes(ax)
    if PANEL_STYLE["show_titles"]:
        ax.set_title("LOH by myometrial invasion")
    if add_letter:
        panel_letter(ax, "C")


def plot_multiplicity_summary(ax, family: pd.DataFrame, add_letter=True):
    d = family[family["analysis_family"].isin(ENDPOINT_ORDER)].copy()
    d["order"] = d["analysis_family"].map({x:i for i,x in enumerate(ENDPOINT_ORDER)})
    d = d.sort_values("order")
    y = np.arange(len(d))[::-1]
    metrics = [
        ("min_permutation_p", "Nominal P", "o", PLOT_STYLE["neutral_color"]),
        ("min_family_bh_q", "BH q", "s", PLOT_STYLE["immune_color"]),
        ("min_family_maxT_fwer_p", "maxT P", "D", PLOT_STYLE["genomic_color"]),
    ]
    offsets = [-0.16, 0.0, 0.16]
    for (col, label, marker, color), off in zip(metrics, offsets):
        vals = pd.to_numeric(d[col], errors="coerce").to_numpy(float)
        ax.scatter(vals, y+off, marker=marker, s=46, color=color, label=label, zorder=3)
    ax.axvline(0.05, color=PLOT_STYLE["reference_color"], lw=0.9, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(0.03, 1.1)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS[x] for x in d["analysis_family"]])
    ax.set_xlabel("Minimum P or q within each 8-test family")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=7.5)
    ax.grid(axis="x", alpha=PLOT_STYLE["grid_alpha"], which="both")
    clean_axes(ax)
    if PANEL_STYLE["show_titles"]:
        ax.set_title("Multiplicity summary")
    if add_letter:
        panel_letter(ax, "D")


def make_figure_s6(root: Path, outdir: Path, assoc: pd.DataFrame, family: pd.DataFrame, matrix: pd.DataFrame):
    grade = pd.read_csv(root / "06_grade_associations.csv")
    inv = pd.read_csv(root / "09_myometrial_invasion_associations.csv")

    # Individual panels.
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigureS6A_effect_heatmap"])
    im = plot_effect_heatmap(ax, assoc, add_letter=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Effect size")
    fig.tight_layout()
    save_figure(fig, outdir / "FigureS6A_effect_heatmap")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigureS6B_grade_morisita"])
    plot_grade_morisita(ax, matrix, grade, add_letter=False)
    fig.tight_layout()
    save_figure(fig, outdir / "FigureS6B_grade_morisita")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigureS6C_invasion_loh"])
    plot_invasion_loh(ax, matrix, inv, add_letter=False)
    fig.tight_layout()
    save_figure(fig, outdir / "FigureS6C_invasion_loh")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigureS6D_multiplicity_summary"])
    plot_multiplicity_summary(ax, family, add_letter=False)
    fig.tight_layout()
    save_figure(fig, outdir / "FigureS6D_multiplicity_summary")

    # Composite: full-width heatmap + three equally sized lower panels.
    fig = plt.figure(figsize=FIGURE_SIZE["FigureS6_composite"])
    gs = fig.add_gridspec(2, 3, height_ratios=[1.02, 1.0], hspace=0.48, wspace=0.40)
    axA = fig.add_subplot(gs[0, :])
    im = plot_effect_heatmap(axA, assoc, add_letter=True)
    cbar = fig.colorbar(im, ax=axA, fraction=0.022, pad=0.02)
    cbar.set_label("Effect size")
    axB = fig.add_subplot(gs[1, 0])
    plot_grade_morisita(axB, matrix, grade, add_letter=True)
    axC = fig.add_subplot(gs[1, 1])
    plot_invasion_loh(axC, matrix, inv, add_letter=True)
    axD = fig.add_subplot(gs[1, 2])
    plot_multiplicity_summary(axD, family, add_letter=True)
    save_figure(fig, outdir / "FigureS6_composite")


def plot_forest_all(ax, endpoint_df: pd.DataFrame, endpoint: str):
    d = endpoint_df.copy()
    d["order"] = d["descriptor_feature"].map({k:i for i,k in enumerate(DESCRIPTOR_ORDER)})
    d = d.sort_values("order", ascending=False).reset_index(drop=True)
    y = np.arange(len(d))
    for i, r in d.iterrows():
        est = float(r["effect_size"])
        lo = float(r["bootstrap_ci_low"])
        hi = float(r["bootstrap_ci_high"])
        panel = str(r["panel"])
        color = PLOT_STYLE["immune_color"] if panel == "IMMUNE" else PLOT_STYLE["genomic_color"]
        marker = "o" if panel == "IMMUNE" else "s"
        ax.plot([lo, hi], [i, i], color=color, lw=PLOT_STYLE["ci_line_width"], zorder=2)
        ax.scatter([est], [i], marker=marker, s=42, color=color, edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(0, color=PLOT_STYLE["reference_color"], lw=0.9)
    ax.set_xlim(-1.05,1.05)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(x,x) for x in d["descriptor_feature"]])
    ax.set_xlabel("Spearman ρ" if endpoint == "grade" else "Rank-biserial effect")
    ax.set_title(LABELS[endpoint], loc="left", fontweight="bold", fontsize=PLOT_STYLE["title_size"])
    ax.grid(axis="x", alpha=PLOT_STYLE["grid_alpha"])
    clean_axes(ax)
    # Optional q labels at right.
    if PANEL_STYLE["forest_show_global_q"]:
        for i, r in d.iterrows():
            ax.text(1.02, i, f"q={float(r['global_bh_q']):.2f}", va="center", fontsize=7)


def make_figure_s7(outdir: Path, assoc: pd.DataFrame):
    fig, axes = plt.subplots(3, 1, figsize=FIGURE_SIZE["FigureS7_forest_all_associations"], sharex=False)
    for ax, endpoint, letter in zip(axes, ENDPOINT_ORDER, ["A","B","C"]):
        plot_forest_all(ax, assoc[assoc["analysis_family"].eq(endpoint)], endpoint)
        panel_letter(ax, letter)
    handles = [
        Line2D([0],[0], marker='o', linestyle='none', markerfacecolor=PLOT_STYLE["immune_color"], markeredgecolor='white', label='Immune'),
        Line2D([0],[0], marker='s', linestyle='none', markerfacecolor=PLOT_STYLE["genomic_color"], markeredgecolor='white', label='Genomic'),
    ]
    axes[0].legend(handles=handles, frameon=False, loc="lower right")
    fig.subplots_adjust(hspace=0.46)
    save_figure(fig, outdir / "FigureS7_forest_all_associations")

# =====================================================================
# Publication data package
# =====================================================================

def read_subject_map(path: Optional[Path], source_subjects: Sequence[str]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    m = pd.read_csv(path, dtype={"subject": str, "public_id": str})
    if not {"subject", "public_id"}.issubset(m.columns):
        raise ValueError("--subject_map must contain subject,public_id")
    if m["subject"].duplicated().any() or m["public_id"].duplicated().any():
        raise ValueError("subject_map must be one-to-one")
    if set(m["subject"]) != set(map(str, source_subjects)):
        raise ValueError("subject_map subjects do not exactly match Step-8 subjects")
    return m[["subject","public_id"]].copy()


def build_publication_data(root: Path, outdir: Path, subject_map: Optional[pd.DataFrame]):
    ensure_dir(outdir)
    copies = {
        "01_frozen_descriptor_panel_contract.csv": "04_frozen_descriptor_panel_contract.csv",
        "02_all_histology_associations.csv": "17_all_histology_associations.csv",
        "03_grade_loso_summary.csv": "08_grade_loso_summary.csv",
        "04_myometrial_invasion_loso_summary.csv": "11_myometrial_invasion_loso_summary.csv",
        "05_lvsi_loso_summary.csv": "14_lvsi_loso_summary.csv",
        "06_histotype_descriptive_summary.csv": "15_histotype_descriptive_summary.csv",
        "07_p53_descriptive_summary.csv": "16_p53_descriptive_summary.csv",
        "08_family_multiplicity_summary.csv": "20_family_multiplicity_summary.csv",
        "09_upstream_provenance_contract.csv": "01_upstream_provenance_contract.csv",
    }
    for dst, src in copies.items():
        shutil.copy2(root / src, outdir / dst)

    matrix = pd.read_csv(root / "05_histology_descriptor_analysis_matrix.csv", dtype={"subject": str})
    if subject_map is not None:
        m = matrix.merge(subject_map, on="subject", validate="one_to_one")
        m = m.drop(columns=["subject"]).rename(columns={"public_id":"subject"})
        cols = ["subject"] + [c for c in m.columns if c != "subject"]
        m[cols].to_csv(outdir / "10_subject_level_analysis_matrix_public_ids.csv", index=False)
    else:
        (outdir / "10_SUBJECT_LEVEL_MATRIX_NOT_EXPORTED.txt").write_text(
            "The Step-8 subject-level matrix was intentionally not copied into the publication data package because no subject-to-public-ID map was supplied.\n",
            encoding="utf-8",
        )


def zip_directory(src: Path, zip_path: Path):
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(src.rglob("*")):
            if p.is_file():
                zf.write(p, p.relative_to(src.parent))

# =====================================================================
# Main
# =====================================================================

def main():
    ap = argparse.ArgumentParser(description="Publication reporting for Step-8 v3 frozen panels vs histology; no new inference.")
    ap.add_argument("--step8", type=Path, required=True, help="Step-8 output directory or ZIP")
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--subject_map", type=Path, default=None, help="Optional CSV: subject,public_id")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite.")
    if args.overwrite and args.out_dir.exists():
        shutil.rmtree(args.out_dir)
    ensure_dir(args.out_dir)

    apply_plot_style()
    root, temp_handle = source_root(args.step8)
    try:
        required = [
            "00_run_config.json", "01_upstream_provenance_contract.csv",
            "02_histology_endpoint_audit.csv", "04_frozen_descriptor_panel_contract.csv",
            "05_histology_descriptor_analysis_matrix.csv", "06_grade_associations.csv",
            "08_grade_loso_summary.csv", "09_myometrial_invasion_associations.csv",
            "11_myometrial_invasion_loso_summary.csv", "12_lvsi_associations.csv",
            "14_lvsi_loso_summary.csv", "15_histotype_descriptive_summary.csv",
            "16_p53_descriptive_summary.csv", "17_all_histology_associations.csv",
            "18_source_sha256_manifest.csv", "20_family_multiplicity_summary.csv",
        ]
        require_files(root, required)
        v = validate_step8(root)
        assoc = v["assoc"]
        family = v["family"]
        matrix = v["matrix"]

        tables_dir = ensure_dir(args.out_dir / "main_tables")
        figures_dir = ensure_dir(args.out_dir / "supplementary_figures")
        data_dir = ensure_dir(args.out_dir / "publication_data")

        # Tables.
        s10 = make_table_s10(root, assoc, family)
        s10.to_csv(tables_dir / "Supplementary_Table_S10_histology_endpoint_multiplicity_summary.csv", index=False)
        write_markdown_table(s10, tables_dir / "Supplementary_Table_S10_histology_endpoint_multiplicity_summary.md")

        s11 = make_table_s11(assoc)
        s11.to_csv(tables_dir / "Supplementary_Table_S11_all_histology_associations.csv", index=False)
        write_markdown_table(s11, tables_dir / "Supplementary_Table_S11_all_histology_associations.md")

        s12 = make_table_s12(root)
        s12.to_csv(tables_dir / "Supplementary_Table_S12_histotype_p53_descriptive.csv", index=False)
        write_markdown_table(s12, tables_dir / "Supplementary_Table_S12_histotype_p53_descriptive.md")

        # Figures.
        make_figure_s6(root, figures_dir, assoc, family, matrix)
        make_figure_s7(figures_dir, assoc)

        # Publication data.
        subject_map = read_subject_map(args.subject_map, matrix["subject"].astype(str).tolist())
        build_publication_data(root, data_dir, subject_map)
        zip_directory(data_dir, args.out_dir / "Supplementary_Data_4.zip")

        # Validation.
        validation_rows = [{"check": k, "pass": bool(val)} for k, val in v["checks"].items()]
        validation_rows += [
            {"check": "Table_S10_rows_3", "pass": len(s10) == 3},
            {"check": "Table_S11_rows_24", "pass": len(s11) == 24},
            {"check": "No_family_corrected_associations", "pass": int((assoc["family_bh_q"] < .05).sum()) == 0 and int((assoc["family_maxT_fwer_p"] < .05).sum()) == 0},
            {"check": "FigureS6_created", "pass": (figures_dir / "FigureS6_composite.pdf").exists()},
            {"check": "FigureS7_created", "pass": (figures_dir / "FigureS7_forest_all_associations.pdf").exists()},
        ]
        validation = pd.DataFrame(validation_rows)
        validation.to_csv(args.out_dir / "02_reporting_validation.csv", index=False)
        if not validation["pass"].all():
            raise ValueError("Reporting validation failed.")

        # Manifest.
        manifest_rows = []
        for p in sorted(args.out_dir.rglob("*")):
            if p.is_file() and p.name not in {"01_reporting_manifest.csv", "00_reporting_config.json"}:
                manifest_rows.append({
                    "relative_path": str(p.relative_to(args.out_dir)),
                    "size_bytes": int(p.stat().st_size),
                    "sha256": sha256_file(p),
                })
        pd.DataFrame(manifest_rows).to_csv(args.out_dir / "01_reporting_manifest.csv", index=False)

        cfg = {
            "script": Path(__file__).name,
            "script_version": SCRIPT_VERSION,
            "matplotlib_version": matplotlib.__version__,
            "step8_input": str(args.step8.resolve()),
            "step8_source_script_version": v["config"].get("script_version"),
            "n_subjects": int(v["config"].get("n_subjects", 0)),
            "n_histology_tests": int(len(assoc)),
            "n_family_bh_q_lt_0_05": int((pd.to_numeric(assoc["family_bh_q"], errors="coerce") < .05).sum()),
            "n_family_maxT_p_lt_0_05": int((pd.to_numeric(assoc["family_maxT_fwer_p"], errors="coerce") < .05).sum()),
            "reporting_contract": "Reporting only; no inferential statistic is recomputed. Effect estimates and inferential quantities are read from frozen Step-8 outputs.",
            "subject_level_publication_matrix_exported": subject_map is not None,
            "font_requested": PLOT_STYLE["font_family"],
            "font_used_in_this_environment": resolve_font_family(),
            "plot_style": PLOT_STYLE,
            "figure_size": FIGURE_SIZE,
            "panel_style": PANEL_STYLE,
        }
        (args.out_dir / "00_reporting_config.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")

        print(f"[OK] {SCRIPT_VERSION}")
        print(f"[OK] Step-8 subjects: {v['config'].get('n_subjects')}")
        print(f"[OK] Histology tests: {len(assoc)}")
        print(f"[OK] Family BH q < 0.05: {(assoc['family_bh_q'] < .05).sum()}")
        print(f"[OK] Family maxT P < 0.05: {(assoc['family_maxT_fwer_p'] < .05).sum()}")
        print(f"[OK] Output: {args.out_dir}")
    finally:
        if temp_handle is not None:
            temp_handle.cleanup()

if __name__ == "__main__":
    main()
