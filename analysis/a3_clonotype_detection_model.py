"""A3. Clonotype-level model of peripheral detection, and intratumoral enrichment of
shared clonotypes (Results, second section; Figure 2D-E).

1. Per patient: logistic regression of PBMC detection on log10 intratumoral frequency.
2. All biopsy clonotypes together: the same model with patient fixed effects, so that the
   slope is estimated within patients; cluster-robust (sandwich) standard errors with the
   13 patients as clusters and a G/(G-1) small-sample correction.
3. Among shared clonotypes: fraction at least twofold more frequent in tumor than in blood
   and vice versa, compared across patients by the Wilcoxon signed-rank test.

Input : Step 2 output (02_exact_clonotype_overlap_500).
Output: <output_dir>/analysis/clonotype_detection_model.json, per_patient_slopes.csv
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats
from _common import load_config, irls


def main():
    cfg = load_config(description=__doc__); OUT = cfg.out_dir("analysis")
    D = cfg.result("02_exact_clonotype_overlap_500")
    t = pd.read_csv(os.path.join(D, "within_subject_clonotype_table.csv.gz"), low_memory=False)
    b = t[t.in_biopsy].copy(); b["y"] = b.in_pbmc.astype(int); b["logf"] = np.log10(b.biopsy_frequency)
    b = b[np.isfinite(b.logf)]

    rows = []
    for s, g in b.groupby("subject"):
        if g.y.nunique() < 2: rows.append(dict(subject=s, n=len(g), slope=np.nan)); continue
        be, _, H = irls(np.column_stack([np.ones(len(g)), g.logf.values]), g.y.values.astype(float))
        rows.append(dict(subject=s, n=len(g), slope=be[1], se=np.sqrt(np.linalg.inv(H)[1, 1]), or_per_decade=np.exp(be[1])))
    r = pd.DataFrame(rows); sl = r.slope.dropna()
    r.drop(columns="subject").to_csv(os.path.join(OUT, "per_patient_slopes.csv"), index=False)
    print(f"1. per-patient slopes: positive in {int((sl > 0).sum())}/{len(sl)}; median OR per decade {np.exp(sl.median()):.2f}; "
          f"Wilcoxon vs 0 P = {stats.wilcoxon(sl).pvalue:.5f}")

    subs = sorted(b.subject.unique())
    X = np.column_stack([np.column_stack([(b.subject == s).astype(float).values for s in subs]), b.logf.values])
    y = b.y.values.astype(float); be, p, H = irls(X, y); k = X.shape[1] - 1
    Hinv = np.linalg.inv(H); meat = np.zeros((X.shape[1], X.shape[1])); res = y - p
    for s in subs:
        m = (b.subject == s).values; u = X[m].T @ res[m]; meat += np.outer(u, u)
    G = len(subs); V = Hinv @ meat @ Hinv * G / (G - 1); se = np.sqrt(V[k, k])
    OR, lo, hi = np.exp(be[k]), np.exp(be[k] - 1.96 * se), np.exp(be[k] + 1.96 * se)
    pval = 2 * stats.t.sf(abs(be[k] / se), G - 1)
    print(f"2. fixed-effects model: {len(b):,} clonotypes, {G} clusters; OR per decade {OR:.2f} (95% CI {lo:.2f}-{hi:.2f}); "
          f"z = {be[k]/se:.2f}, P (t, {G-1} df) = {pval:.2g}")

    raw = pd.read_csv(os.path.join(D, "within_subject_overlap_summary_raw.csv"))
    bi, pb = raw.shared_fraction_biopsy_enriched_ge2x, raw.shared_fraction_pbmc_enriched_ge2x
    lr = raw.shared_median_log2_biopsy_pbmc_ratio; pw = stats.wilcoxon(bi, pb).pvalue
    print(f"3. shared clonotypes >=2-fold enriched in tumor: median {bi.median():.3f} ({bi.min():.3f}-{bi.max():.3f}); "
          f"in blood: median {pb.median():.3f} ({pb.min():.3f}-{pb.max():.3f}); tumor > blood in {int((bi > pb).sum())}/{len(bi)}; "
          f"Wilcoxon P = {pw:.5f}; median log2 tumor/blood {lr.median():.2f}")

    json.dump(dict(n_clonotypes=int(len(b)), n_patients=G, or_per_decade=OR, ci95=[lo, hi], p_t=pval,
                   slopes_positive=int((sl > 0).sum()), slopes_n=int(len(sl)),
                   enriched_tumor_median=float(bi.median()), enriched_blood_median=float(pb.median()),
                   tumor_gt_blood=int((bi > pb).sum()), enrichment_wilcoxon_p=float(pw),
                   median_log2_ratio=float(lr.median())),
              open(os.path.join(OUT, "clonotype_detection_model.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
