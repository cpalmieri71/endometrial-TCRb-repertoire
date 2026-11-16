"""A2. Sensitivity of the compartmental comparisons to sequencing-depth and coverage
standardization (Supplementary Data 3; Supplementary Figure S1D-E).

Three standardizations, each compared with the observed values by a two-sided Wilcoxon
signed-rank test across pairs:
  1. per-pair rarefaction of the deeper repertoire to the depth of the shallower one
     (100 subsamples without replacement, metrics averaged);
  2. cohort-wide rarefaction to one common depth, excluding the pairs that contain the
     shallowest repertoires (0, 1 or 2 excluded);
  3. standardization to a common Good's coverage, the highest value attainable in every
     repertoire (Chao & Jost).
Rarefaction curves are exact hypergeometric expectations.

Inputs : MiXCR tables (repseq_dir); A1 output recomputed_metrics.csv.
Outputs: <output_dir>/analysis/ rarefaction_curves.csv, paired_rarefied_perpair.csv,
         paired_rarefied_cohort_<n>.csv, paired_coverage_standardized.csv,
         depth_standardisation.csv, depth_standardisation.xlsx, *_results.json
Seeds  : 20261110 (rarefaction), 7 (coverage standardization). The order of random draws
         is part of the result; do not reorder the steps.
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats
from scipy.special import gammaln
from _common import load_config, patient_codes, repseq_files

REPS = 100


def exact_rich(c, n):
    """E[S_n] = sum_i 1 - C(N-c_i, n)/C(N, n): exact hypergeometric rarefaction."""
    N = c.sum()
    if n >= N: return float((c > 0).sum())
    lgN = gammaln(N + 1) - gammaln(n + 1) - gammaln(N - n + 1)
    M = N - c; ok = M >= n
    lg = np.full(len(c), -np.inf)
    lg[ok] = gammaln(M[ok] + 1) - gammaln(n + 1) - gammaln(M[ok] - n + 1)
    return float(np.sum(1 - np.exp(lg - lgN)))


def hill(c):
    c = c[c > 0]; N = c.sum(); p = c / N; S = len(c); ps = np.sort(p)[::-1]
    return S, np.exp(-(p * np.log(p)).sum()), 1.0 / (p ** 2).sum(), ps[:min(10, S)].sum()


def mc(rng, c, n, reps=REPS):
    if n >= c.sum():
        S, D1, D2, t10 = hill(c)
        return dict(S=S, D1=D1, D2=D2, T10=t10, S_sd=0, D1_sd=0, D2_sd=0, T10_sd=0)
    a = np.array([hill(rng.multivariate_hypergeometric(c, n, method="count")) for _ in range(reps)])
    return dict(S=a[:, 0].mean(), D1=a[:, 1].mean(), D2=a[:, 2].mean(), T10=a[:, 3].mean(),
                S_sd=a[:, 0].std(ddof=1), D1_sd=a[:, 1].std(ddof=1), D2_sd=a[:, 2].std(ddof=1), T10_sd=a[:, 3].std(ddof=1))


def report(df, label):
    print(f"\n--- {label} (n={len(df)}) ---"); res = {}
    for m, nm in [("S", "observed richness"), ("D1", "Hill D1"), ("D2", "Hill D2"), ("T10", "Top-10")]:
        b, p = df[f"{m}_B"].values, df[f"{m}_P"].values; w = stats.wilcoxon(b, p)
        print(f"  {nm:18s} lower in biopsy {int((b < p).sum()):2d}/{len(df)}  median B {np.median(b):10.4g}  P {np.median(p):10.4g}  Wilcoxon P = {w.pvalue:.5f}")
        res[m] = dict(lower=int((b < p).sum()), n=len(df), med_B=float(np.median(b)), med_P=float(np.median(p)), p=float(w.pvalue))
    return res


def stats_of(c):
    c = c[c > 0]; N = c.sum(); p = c / N; S = len(c); ps = np.sort(p)[::-1]
    return dict(S=S, D1=np.exp(-(p * np.log(p)).sum()), D2=1.0 / (p ** 2).sum(), T10=ps[:min(10, S)].sum(), cov=1 - (c == 1).sum() / N)


def main():
    cfg = load_config(description=__doc__); OUT = cfg.out_dir("analysis")
    dat = {k: pd.read_csv(f, sep="\t").readCount.values.astype(np.int64) for k, f in repseq_files(cfg).items()}
    subjects = sorted({k[0] for k in dat}); depths = {k: int(v.sum()) for k, v in dat.items()}

    # ---------------- rarefaction (seed 20261110) ----------------
    rng = np.random.default_rng(20261110)
    curve = []
    for k, c in dat.items():
        N = depths[k]
        grid = np.unique(np.concatenate([np.array([1, 100, 500]), np.round(np.logspace(np.log10(1000), np.log10(N), 28)).astype(int)]))
        for n in grid[grid <= N]:
            curve.append(dict(subject=k[0], compartment=k[1], n_reads=int(n), richness=exact_rich(c, int(n))))
    pd.DataFrame(curve).to_csv(f"{OUT}/rarefaction_curves.csv", index=False)
    # Agreement of the exact expectation with Monte Carlo rarefaction, checked on one PBMC
    # repertoire (patient code P5) at 63,199 reads. This check draws from the same random
    # stream as the analyses below and is kept so that their results are reproduced exactly.
    code2id = {v: k for k, v in patient_codes(cfg).items()}
    c = dat[(code2id["P5"], "PBMC")]
    print("exact vs Monte Carlo expected richness (P5 PBMC, 63,199 reads): %.1f vs %.1f" % (exact_rich(c, 63199), mc(rng, c, 63199, 30)["S"]))

    rows = []
    for s in subjects:
        nB, nP = depths[(s, "BIOPSY")], depths[(s, "PBMC")]; n = min(nB, nP)
        b = mc(rng, dat[(s, "BIOPSY")], n); p = mc(rng, dat[(s, "PBMC")], n)
        rows.append(dict(subject=s, common_n=n, shallower="BIOPSY" if nB <= nP else "PBMC",
                         S_B=b["S"], S_P=p["S"], D1_B=b["D1"], D1_P=p["D1"], D2_B=b["D2"], D2_P=p["D2"], T10_B=b["T10"], T10_P=p["T10"],
                         S_B_sd=b["S_sd"], S_P_sd=p["S_sd"], D1_B_sd=b["D1_sd"], D1_P_sd=p["D1_sd"]))
    A = pd.DataFrame(rows); A.to_csv(f"{OUT}/paired_rarefied_perpair.csv", index=False)
    resA = report(A, "per-pair rarefaction to the shallower member")

    # pairs ordered by their shallowest repertoire; exclude the 0, 1 or 2 shallowest
    order = sorted(subjects, key=lambda s: min(depths[(s, "BIOPSY")], depths[(s, "PBMC")]))
    resB = {}
    for k in (0, 1, 2):
        subs = [s for s in subjects if s not in order[:k]]
        n = min(depths[(s, c)] for s in subs for c in ("BIOPSY", "PBMC")); rr = []
        for s in subs:
            b = mc(rng, dat[(s, "BIOPSY")], n); p = mc(rng, dat[(s, "PBMC")], n)
            rr.append(dict(subject=s, S_B=b["S"], S_P=p["S"], D1_B=b["D1"], D1_P=p["D1"], D2_B=b["D2"], D2_P=p["D2"], T10_B=b["T10"], T10_P=p["T10"]))
        df = pd.DataFrame(rr); df.to_csv(f"{OUT}/paired_rarefied_cohort_{len(subs)}.csv", index=False)
        resB[f"{len(subs)} pairs"] = dict(common_n=int(n), res=report(df, f"cohort-wide common depth, {n:,} reads"))
    json.dump(dict(perpair=resA, cohort=resB), open(f"{OUT}/rarefaction_results.json", "w"), indent=2)

    # ---------------- coverage standardization (seed 7) ----------------
    rng = np.random.default_rng(7)

    def sub_stats(c, n, reps=60):
        if n >= c.sum(): return stats_of(c)
        acc = [stats_of(rng.multivariate_hypergeometric(c, int(n), method="count")) for _ in range(reps)]
        return {k: float(np.mean([a[k] for a in acc])) for k in acc[0]}

    def depth_for_coverage(c, target, lo=50):
        hi = int(c.sum())
        if stats_of(c)["cov"] < target: return hi
        for _ in range(26):
            mid = (lo + hi) // 2
            if mid <= lo: break
            if sub_stats(c, mid, reps=15)["cov"] < target: lo = mid
            else: hi = mid
        return hi

    TARGET = float(min(stats_of(v)["cov"] for v in dat.values()))
    print(f"\ncommon attainable Good's coverage = {TARGET:.4f}")
    rows = []
    for s in subjects:
        r = {}
        for comp in ("BIOPSY", "PBMC"):
            c = dat[(s, comp)]; n = depth_for_coverage(c, TARGET); st = sub_stats(c, n, reps=60); t = comp[0]
            r.update({f"n_{t}": n, f"S_{t}": st["S"], f"D1_{t}": st["D1"], f"D2_{t}": st["D2"], f"T10_{t}": st["T10"], f"cov_{t}": st["cov"]})
        r["subject"] = s; rows.append(r)
    C = pd.DataFrame(rows); C.to_csv(f"{OUT}/paired_coverage_standardized.csv", index=False)
    resC = report(C, f"coverage standardization (Good's coverage = {TARGET:.3f})")
    json.dump(dict(target_coverage=TARGET, results=resC), open(f"{OUT}/coverage_results.json", "w"), indent=2)

    # depth no longer drives biopsy richness after standardization
    ob = pd.read_csv(f"{OUT}/recomputed_metrics.csv"); ob = ob[ob.compartment == "BIOPSY"]
    p11 = pd.read_csv(f"{OUT}/paired_rarefied_cohort_11.csv").merge(ob[["subject", "depth"]], on="subject")
    print("\nSpearman, biopsy depth vs richness: observed rho=%.3f P=%.4f" % stats.spearmanr(ob.depth, ob.observed_richness))
    print("                                   rarefied (11 pairs) rho=%.3f P=%.4f" % stats.spearmanr(p11.depth, p11.S_B))
    mcv = C.merge(ob[["subject", "depth"]], on="subject")
    print("                                   coverage-standardized rho=%.3f P=%.4f" % stats.spearmanr(mcv.depth, mcv.S_B))

    # ---------------- summary table ----------------
    w = ob.pipe(lambda _: pd.read_csv(f"{OUT}/recomputed_metrics.csv")).pivot(index="subject", columns="compartment",
                values=["observed_richness", "hill_D1", "hill_D2", "top10_frequency"])
    M = [("observed_richness", "S", "Observed clonotypes"), ("hill_D1", "D1", "Hill D1"), ("hill_D2", "D2", "Hill D2"), ("top10_frequency", "T10", "Top-10 frequency")]
    rows = []
    def add(label, n, common, b, p, metric):
        rows.append(dict(analysis=label, standardisation=common, metric=metric, n_pairs=n, n_lower_in_biopsy=int((b < p).sum()),
                         median_biopsy=float(np.median(b)), median_pbmc=float(np.median(p)), wilcoxon_p=float(stats.wilcoxon(b, p).pvalue)))
    for col, sh, nm in M: add("Observed (as reported)", 13, "none (raw depth)", w[(col, "BIOPSY")].values, w[(col, "PBMC")].values, nm)
    for col, sh, nm in M: add("Per-pair rarefaction", 13, "equal depth within pair", A[f"{sh}_B"].values, A[f"{sh}_P"].values, nm)
    for nsub in (12, 11):
        Cc = pd.read_csv(f"{OUT}/paired_rarefied_cohort_{nsub}.csv"); dep = f"{resB[f'{nsub} pairs']['common_n']:,} reads (cohort-wide)"
        for col, sh, nm in M: add(f"Cohort common depth (n={nsub})", nsub, dep, Cc[f"{sh}_B"].values, Cc[f"{sh}_P"].values, nm)
    for col, sh, nm in M: add("Coverage standardisation", 13, f"Good's coverage = {TARGET:.3f}", C[f"{sh}_B"].values, C[f"{sh}_P"].values, nm)
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/depth_standardisation.csv", index=False)
    write_xlsx(t, f"{OUT}/depth_standardisation.xlsx")
    print("\nwrote", f"{OUT}/depth_standardisation.xlsx")


def write_xlsx(t, path):
    import openpyxl
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    order = ["Observed clonotypes", "Hill D1", "Hill D2", "Top-10 frequency"]
    an = ["Observed (as reported)", "Per-pair rarefaction", "Cohort common depth (n=12)", "Cohort common depth (n=11)", "Coverage standardisation"]
    t = t.copy(); t["metric"] = pd.Categorical(t.metric, order, ordered=True); t["analysis"] = pd.Categorical(t.analysis, an, ordered=True)
    t = t.sort_values(["analysis", "metric"])
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Depth standardisation"
    H = Font(name="Arial", size=10, bold=True, color="FFFFFF"); N = Font(name="Arial", size=10)
    HF = PatternFill("solid", fgColor="2F4858"); thin = Side(style="thin", color="D6D6D2")
    ws.append(["Analysis", "Standardisation", "Metric", "Pairs (n)", "Lower in biopsy (n)", "Median biopsy", "Median PBMC", "Wilcoxon P"])
    for c in ws[1]: c.font = H; c.fill = HF; c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for _, r in t.iterrows():
        ws.append([str(r.analysis), r.standardisation, str(r.metric), int(r.n_pairs), int(r.n_lower_in_biopsy), float(r.median_biopsy), float(r.median_pbmc), float(r.wilcoxon_p)])
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = N; c.border = Border(bottom=thin)
            if c.column in (4, 5): c.number_format = "0"; c.alignment = Alignment(horizontal="center")
            elif c.column in (6, 7): c.number_format = "#,##0.000"; c.alignment = Alignment(horizontal="right")
            elif c.column == 8: c.number_format = "0.0000"; c.alignment = Alignment(horizontal="right")
    for col, wd in zip("ABCDEFGH", [26, 30, 20, 9, 12, 14, 14, 12]): ws.column_dimensions[col].width = wd
    ws.row_dimensions[1].height = 34; ws.freeze_panes = "D2"
    n = wb.create_sheet("Notes"); n.append(["Item", "Value"])
    for c in n[1]: c.font = H; c.fill = HF
    notes = [("Purpose", "Sensitivity of the paired tumour-biopsy / PBMC comparisons to sequencing-depth and coverage standardisation."),
             ("Observed (as reported)", "Measures at the achieved sequencing depth, as reported in the main text and Supplementary Data 1."),
             ("Per-pair rarefaction", "Within each pair, the deeper repertoire was randomly subsampled without replacement to the depth of the shallower member (common depth 6,194-657,278 reads; 100 subsamples, metrics averaged)."),
             ("Cohort common depth", "All repertoires rarefied to one depth across the cohort: 12,792 reads (12 pairs) or 63,199 reads (11 pairs), retaining pairs in which both members reached that depth."),
             ("Coverage standardisation", "All 26 repertoires standardised to Good's coverage = 0.838, the highest value attainable in every repertoire (Chao & Jost)."),
             ("Test", "Two-sided Wilcoxon signed-rank test on paired values. P values are descriptive; no correction for multiple comparisons."),
             ("Rarefaction curves", "Computed as exact hypergeometric expectations, not by simulation. Shown in Supplementary Figure S1D."),
             ("Software", "Python 3.13; NumPy; pandas; SciPy. Subsampling seed recorded in the analysis script.")]
    for k, v in notes:
        n.append([k, v]); r = n.max_row
        n.cell(r, 1).font = Font(name="Arial", size=10, bold=True); n.cell(r, 1).alignment = Alignment(vertical="top")
        n.cell(r, 2).font = N; n.cell(r, 2).alignment = Alignment(wrap_text=True, vertical="top")
    n.column_dimensions["A"].width = 26; n.column_dimensions["B"].width = 104; n.freeze_panes = "A2"
    wb.save(path)


if __name__ == "__main__":
    main()
