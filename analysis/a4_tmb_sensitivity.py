"""A4. Sensitivity of the TMB-diversity association (Results, genomic features; Figure 3).

Checks whether the inverse association between tumor mutational burden (TMB) and
intratumoral Hill D1 is explained by sequencing depth or driven by a few tumors:
  - Spearman correlations of TMB and D1 with biopsy depth and Good's coverage;
  - partial Spearman correlation of TMB and D1 controlling for biopsy depth;
  - exact-permutation P values (100,000 permutations, seed 1) after excluding the two
    highest-TMB tumors, all tumors carrying the assay POLE flag, and the one or two
    shallowest biopsies;
  - D1 in the hypermutated tumors versus the rest (Mann-Whitney);
  - correlations among the four genomic measures.

Inputs : Step 6 subject matrix (06_immune_genomic_associations/12_subject_analysis_matrix.csv);
         A1 output recomputed_metrics.csv.
Output : <output_dir>/analysis/tmb_sensitivity.csv, tmb_sensitivity.json
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats
from _common import load_config


def sp(x, y):
    r = stats.spearmanr(x, y); return float(r.statistic), float(r.pvalue)


def partial(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r = stats.pearsonr(ex, ey); return float(r.statistic), float(r.pvalue)


def perm_p(x, y, n=100000, seed=1):
    x = np.asarray(x, float); y = np.asarray(y, float); obs = abs(stats.spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    rx = stats.rankdata(x); ry = stats.rankdata(y); rx = rx - rx.mean(); ry = ry - ry.mean()
    den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    P = np.argsort(rng.random((n, len(x))), axis=1)
    vals = np.abs((rx[None, :] * ry[P]).sum(axis=1) / den)
    return float((np.sum(vals >= obs - 1e-15) + 1) / (n + 1))


def main():
    cfg = load_config(description=__doc__); OUT = cfg.out_dir("analysis")
    m = pd.read_csv(cfg.result("06_immune_genomic_associations", "12_subject_analysis_matrix.csv"))
    rep = pd.read_csv(os.path.join(OUT, "recomputed_metrics.csv"))
    b = rep[rep.compartment == "BIOPSY"][["subject", "depth", "goods_coverage", "hill_D1"]]
    d = m.merge(b, on="subject", suffixes=("", "_rep"))
    assert np.allclose(d.biopsy_hill_D1, d.hill_D1, atol=1e-3), "Step 6 D1 differs from recomputed D1"
    res = {}
    print("Depth confounding")
    for lab, (x, y) in {"TMB vs biopsy depth": (d.oca_tmb, d.depth), "D1 vs biopsy depth": (d.biopsy_hill_D1, d.depth),
                        "TMB vs D1": (d.oca_tmb, d.biopsy_hill_D1), "TMB vs Good's coverage": (d.oca_tmb, d.goods_coverage)}.items():
        r, p = sp(x, y); res[lab] = dict(rho=r, p=p); print(f"  {lab:24s} rho = {r:+.3f}  P = {p:.4f}")
    r, p = partial(d.oca_tmb, d.biopsy_hill_D1, d.depth); res["partial TMB vs D1 | depth"] = dict(rho=r, p=p)
    print(f"  partial TMB vs D1 | depth rho = {r:+.3f}  P = {p:.4f}")

    top_tmb = d.nlargest(2, "oca_tmb").subject.tolist()
    shallow = d.nsmallest(2, "depth").subject.tolist()
    cases = {"all tumors": [], "without the two highest-TMB tumors": top_tmb,
             "without all tumors with the assay POLE flag": d[d.oca_pole_positive == 1].subject.tolist(),
             "without the shallowest biopsy": shallow[:1], "without the two shallowest biopsies": shallow}
    rows = []
    print("Exclusions (exact permutation P, 100,000 permutations)")
    for lab, drop in cases.items():
        s = d[~d.subject.isin(drop)]; r, p = sp(s.oca_tmb, s.biopsy_hill_D1); pp = perm_p(s.oca_tmb, s.biopsy_hill_D1)
        rows.append(dict(analysis=lab, n=len(s), rho=r, p_asymptotic=p, p_permutation=pp))
        print(f"  {lab:44s} n = {len(s):2d}  rho = {r:+.3f}  permutation P = {pp:.4f}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "tmb_sensitivity.csv"), index=False)

    hyper = d.subject.isin(top_tmb)
    mw = stats.mannwhitneyu(d[hyper].biopsy_hill_D1, d[~hyper].biopsy_hill_D1).pvalue
    res["hypermutated vs rest, D1 Mann-Whitney P"] = float(mw)
    print(f"Hypermutated (n=2) vs rest D1: median {d[hyper].biopsy_hill_D1.median():.0f} vs {d[~hyper].biopsy_hill_D1.median():.0f}; Mann-Whitney P = {mw:.3f}")
    G = ["oca_tmb", "oca_msi_score", "oca_gim", "oca_loh"]; C = d[G].corr(method="spearman")
    print("Spearman correlations among genomic measures\n" + C.round(2).to_string())
    res["genomic_measures_max_abs_rho"] = float(np.abs(C.values[~np.eye(4, dtype=bool)]).max())
    json.dump(res, open(os.path.join(OUT, "tmb_sensitivity.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
