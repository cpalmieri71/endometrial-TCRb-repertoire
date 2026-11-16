---
title-meta: "Step 2: Exact TCRβ clonotype overlap"
author: []
date: ""
lang: en-GB
papersize: a4
fontsize: 10pt
geometry:
  - top=18mm
  - bottom=18mm
  - left=18mm
  - right=18mm
mainfont: DejaVu Serif
sansfont: DejaVu Sans
monofont: DejaVu Sans Mono
monofontoptions:
  - Scale=0.76
colorlinks: true
header-includes:
  - |
    \usepackage{titlesec}
    \usepackage{fancyhdr}
    \usepackage{microtype}
    \usepackage{needspace}
    \titleformat{\section}{\large\sffamily\bfseries}{\thesection}{0.5em}{}
    \titleformat{\subsection}{\normalsize\sffamily\bfseries}{\thesubsection}{0.5em}{}
    \titlespacing*{\section}{0pt}{9pt}{5pt}
    \titlespacing*{\subsection}{0pt}{7pt}{4pt}
    \setlength{\parskip}{4pt}
    \setlength{\emergencystretch}{3em}
    \AtBeginDocument{\setlength{\abovedisplayskip}{5pt}\setlength{\belowdisplayskip}{5pt}}
    \pagestyle{fancy}
    \fancyhf{}
    \fancyhead[L]{\small\sffamily Step 2 | Exact clonotype overlap}
    \fancyhead[R]{\small\sffamily Script reference}
    \fancyfoot[L]{\footnotesize\sffamily Read-based analysis}
    \fancyfoot[R]{\small\sffamily\thepage}
    \renewcommand{\headrulewidth}{0.3pt}
    \setlength{\headheight}{14pt}
---

# Step 2: Exact TCRβ clonotype overlap

**Script:** `2-analyze_exact_clonotype_overlap(2).py`  
**Source label:** frozen candidate v1.2, as stated in the script header.

# 1. What it does

The script quantifies exact clonotype sharing between paired tumor biopsies and peripheral blood mononuclear cells (PBMC), between biopsies from different subjects, and across the biopsy cohort (*publicness*). It retains observed metrics and repeats selected calculations after read downsampling without replacement. It uses **read counts, not cell counts**. Exact sharing means identical CDR3 amino acid strings, not necessarily an identical complete TCR.

It calculates descriptive statistics and conditional resampling distributions, but performs **no hypothesis tests, p-value calculations, regressions or overlap-based sample exclusions**. It does not calculate Hill diversity or generate figures or Excel workbooks.

# 2. How it works

## 2.1 Input preparation and notation

The script scans one directory, without recursion, for `<subject>_PBMC.tsv`, `<subject>_biopsia.tsv` or `<subject>_biopsy.tsv`; the subject identifier is numeric by default. Filenames are matched case-insensitively. Duplicate subject–compartment files stop execution. Complete pairs enter within-subject analyses; all available biopsies enter between-subject and publicness analyses. At least one complete pair and two biopsies are required.

By default, `aaSeqCDR3` identifies clonotypes and `readCount` supplies abundance. Rows with missing/blank identifiers or nonnumeric, nonfinite, zero or negative counts are removed. Identifiers are trimmed; identical strings are collapsed within each file by summing counts. Aggregated counts must be integer-like within $10^{-8}$ and are rounded to integers. The default input check requires the 20 canonical uppercase amino acid symbols; `--allow_noncanonical_aaseq` bypasses this check. Productive-only input is assumed upstream: this is an input-contract check, not a biological productivity filter.

For an aligned pair, let $b_i$ and $p_i$ be biopsy and PBMC read counts over the union of identifiers, with zero for an undetected sequence. Define:

$$
N_B=\sum_i b_i,\quad N_P=\sum_i p_i,\qquad
f_i=b_i/N_B,\quad g_i=p_i/N_P.
$$

Let $B=\{i:b_i>0\}$, $P=\{i:p_i>0\}$, $S_B=|B|$, $S_P=|P|$, $I=B\cap P$, $n_I=|I|$ and $n_U=|B\cup P|$. Fractions are exported on a **0–1 scale**. Undefined divisions return `NaN`.

## 2.2 Global overlap measures

| Output | Formula |
|:--|:--|
| `n_shared`; `n_union` | $n_I$; $n_U=S_B+S_P-n_I$ |
| `jaccard` | $n_I/n_U$ |
| `overlap_coefficient` | $n_I/\min(S_B,S_P)$ |
| `biopsy_clone_fraction_shared`; `pbmc_clone_fraction_shared` | $n_I/S_B$; $n_I/S_P$ |
| `biopsy_abundance_fraction_shared`; `pbmc_abundance_fraction_shared` | $F_B=\sum_{i\in I}f_i$; $F_P=\sum_{i\in I}g_i$ |
| `morisita_horn` | $\displaystyle 2\sum_i f_i g_i\big/\left(\sum_i f_i^2+\sum_i g_i^2\right)$ |

The raw pair summary also exports $N_B,N_P,S_B,S_P$, the compartment-only counts $S_B-n_I$ and $S_P-n_I$, clone-only fractions $1-n_I/S_B$ and $1-n_I/S_P$, and read-abundance-only fractions $1-F_B$ and $1-F_P$. Directional sharing refers to the denominator, not a direction of cellular migration.

\newpage

## 2.3 Abundance ranks and PBMC detection

Biopsy clonotypes with positive counts are ordered by decreasing read abundance, with ties broken lexicographically by sequence. Let $i_{(r)}$ identify the clonotype at rank $r$, for $r=1,\ldots,S_B$. The clonotype table stores ranks and rank fractions $r/S_B$; PBMC ranks are calculated analogously. Missing-compartment ranks are `NaN`.

**Disjoint rank strata.** Default cutoffs are $a=(0.01,0.05,0.10,0.25,0.50,1.00)$. For $S_B>0$, set $h_0=0$, $h_j=\max(1,\lceil a_j S_B\rceil)$ and force the last boundary to $S_B$. The implementation also clamps boundaries to preserve order and stay within the repertoire. Stratum $H_j$ contains ranks $h_{j-1}+1,\ldots,h_j$:

**top 1%, 1–5%, 5–10%, 10–25%, 25–50%, and 50–100%.**

These are percentiles of **clonotype number**, not cumulative read mass. For any stratum $H$, define $d_i=\mathbf{1}(p_i>0)$ and:

$$
n_H=|H|,\qquad k_H=\sum_{i\in H}d_i,\qquad D_H=k_H/n_H,
$$

$$
M_H=\sum_{i\in H}f_i,\qquad
M_H^{\mathrm{shared}}=\sum_{i\in H}f_i d_i,\qquad
W_H=M_H^{\mathrm{shared}}/M_H.
$$

| Stratum output | Calculation |
|:--|:--|
| `n_biopsy_clonotypes`; `n_detected_in_pbmc` | $n_H$; $k_H$ |
| `pbmc_detection_fraction` | $D_H$: fraction of stratum clonotypes detected in PBMC |
| `biopsy_abundance_mass`; `shared_biopsy_abundance_mass` | $M_H$; $M_H^{\mathrm{shared}}$ |
| `fraction_of_stratum_biopsy_mass_detected_in_pbmc` | $W_H$: shared fraction of that stratum's biopsy read mass |
| `median_biopsy_frequency` | $\operatorname{median}_{i\in H}(f_i)$ |

Empty strata have zero counts and masses, but missing detection fractions and median frequency. Small repertoires can produce empty strata after boundary rounding.

**Cumulative top fractions.** For each $t\in\{0.01,0.05,0.10,0.25,0.50\}$, the script selects the first $K_t=\max(1,\lceil tS_B\rceil)$ clonotypes. It exports $K_t$, $D_{T_t}$, $M_{T_t}$ and $W_{T_t}$ for this set $T_t$. These groups overlap: the top 5% includes the top 1%, unlike the disjoint 1–5% stratum. Their column suffixes are `n_clonotypes`, `pbmc_detection_fraction`, `biopsy_abundance_mass` and `fraction_biopsy_mass_detected_in_pbmc`, prefixed by `top_<percentage>pct_`.

## 2.4 Relative abundance among shared clonotypes

For shared clonotypes only, the script computes:

$$
L_i=\log_2(f_i/g_i),\qquad i\in I.
$$

No pseudocount is used. Nonshared clonotypes are labelled `BIOPSY_ONLY` or `PBMC_ONLY`. Shared clonotypes are labelled `SHARED_BIOPSY_ENRICHED_GE2X` when $L_i\ge1$, `SHARED_PBMC_ENRICHED_GE2X` when $L_i\le-1$, and `SHARED_BALANCED_LT2X` when $-1<L_i<1$. Separate biopsy-enrichment flags use $L_i\ge2$ and $L_i\ge3$ (fourfold and eightfold).

The pair summary reports the median finite $L_i$ and the fractions of finite shared-clonotype ratios satisfying $L_i\ge1$, $L_i\le-1$, $L_i\ge2$ and $L_i\ge3$. All are missing when no finite shared ratio is available. **These are raw, descriptive quantities; they are not resampled or significance-tested.**

\newpage

## 2.5 Matched-depth downsampling

For each biopsy–PBMC pair, the target is $T=\min(N_B,N_P)$. For an original count vector $\mathbf c=(c_1,\ldots,c_K)$ with $N=\sum_i c_i$, `rarefy()` calls NumPy's multivariate hypergeometric generator. Its sampling law is:

$$
\mathbf X\sim\operatorname{MultiHypergeometric}(\mathbf c,T),
\qquad
\Pr(\mathbf X=\mathbf x)=
\frac{\prod_{i=1}^{K}\binom{c_i}{x_i}}{\binom{N}{T}},
\quad \sum_i x_i=T,\quad 0\le x_i\le c_i.
$$

This formula describes the delegated random draw; the script does not separately calculate or export the probability mass. Sampling is **without replacement**, and each iteration starts from the original counts. A repertoire already at the target depth is unchanged. Targets are pair-specific, not common to all patients.

Each iteration recomputes the eight selected pair metrics: shared-clonotype count, Jaccard, overlap coefficient, both directional clone-sharing fractions, both directional shared-abundance fractions and Morisita–Horn. It also recomputes $D_H,W_H$ for each disjoint stratum and $D_{T_t},W_{T_t}$ for each cumulative top fraction. **Biopsy presence, ranks and stratum membership are recalculated from the downsampled biopsy.** Resampling summaries are not exported for every raw descriptor, for example stratum size or enrichment ratios.

### Resampling summary formulas

For a metric, let $x_1,\ldots,x_m$ be the **finite** values from the requested $R$ iterations, and let $x_0$ be its observed value. The function `summarize()` returns:

$$
\bar x=\frac{1}{m}\sum_{r=1}^{m}x_r,\qquad
\widetilde x=\operatorname{median}(x_1,\ldots,x_m),
$$

$$
s_x=\sqrt{\frac{\sum_{r=1}^{m}(x_r-\bar x)^2}{m-1}},
\qquad \mathrm{CV}=s_x/|\bar x|,
\qquad \Delta=x_0-\widetilde x.
$$

| Output field | Definition |
|:--|:--|
| `raw_value`; `n_resamples_valid` | $x_0$; $m$ |
| `std_mean`; `std_median` | $\bar x$; $\widetilde x$ |
| `std_sd`; `std_cv` | Sample standard deviation $s_x$ (`ddof=1`); $s_x/|\bar x|$ |
| `std_q025`; `std_q25`; `std_q75`; `std_q975` | NumPy empirical quantiles at 0.025, 0.25, 0.75 and 0.975 |
| `raw_minus_std_median` | $\Delta=x_0-\widetilde x$ |

`std_` means **depth-standardized**, not z-scored. The standard deviation is missing for $m<2$; CV is missing if the standard deviation is unavailable or the mean is zero. With no finite resamples, all standardized summaries are missing. Requested iteration count and target depth accompany the summaries. Individual iteration values are not written to disk.

The default is **200 iterations**, with seed **20260927**; both are configurable. A single seeded random generator is used sequentially for within-subject, between-biopsy and publicness resampling. The actual iteration count and software versions are recorded in `run_config.json`.

The 2.5–97.5% range describes conditional variation under read downsampling, **not a confidence interval across biological subjects**. Standardization addresses sequencing depth only, not RNA input, T-cell content, extraction, PCR or library preparation. This script does not calculate cohort medians of patient-specific overlap estimates; that is a downstream reporting operation.

\newpage

## 2.6 Between-subject biopsy overlap

All $M(M-1)/2$ unordered pairs of the $M$ available biopsies use the Section 2.2 formulas, with directional output suffixes `_a` and `_b`. Raw depth, richness, shared count and union size are retained. The same eight overlap metrics are resampled at $T_{ab}=\min(N_a,N_b)$ and summarized as in Section 2.5.

Symmetric matrices contain raw values or resampling medians for shared count, Jaccard and Morisita–Horn. Their diagonals are missing, 1 and 1, respectively: display conventions, not additional self-comparison calculations.

## 2.7 Biopsy publicness

Let $B_s$ be the observed clonotype set in biopsy $s$, $f_{is}$ its read-frequency distribution, and $U=\bigcup_{s=1}^{M}B_s$. For each distinct sequence:

$$
v_i=\sum_{s=1}^{M}\mathbf{1}(i\in B_s),\qquad
\pi_i=v_i/M.
$$

The catalogue exports `n_subjects` ($v_i$), `prevalence` ($\pi_i$), the subject list and `total_count` ($\sum_s c_{is}$). Mean, median and maximum frequencies use **only biopsies where the sequence is present**:

$$
\overline f_{i,\mathrm{present}}
=\frac{1}{v_i}\sum_{s:i\in B_s}f_{is}.
$$

For each biopsy $s$, and each recurrence threshold $k\in\{2,3\}$, the script exports:

$$
C_{s,k}=\sum_{i\in B_s}\mathbf{1}(v_i\ge k),\qquad
Q_{s,k}=C_{s,k}/|B_s|,\qquad
A_{s,k}=\sum_{i\in B_s}f_{is}\mathbf{1}(v_i\ge k).
$$

These correspond to `n_public_clonotypes_ge<k>subjects`, `public_clone_fraction_ge<k>subjects` and `public_abundance_fraction_ge<k>subjects`. It also exports biopsy depth, richness and $\max_{i\in B_s}v_i$: despite the word “prevalence” in `max_cross_subject_prevalence_for_any_clone`, this last quantity is a **subject count**, not a fraction.

The cohort recurrence distribution is:

$$
H_j=\sum_{i\in U}\mathbf{1}(v_i=j),\qquad
h_j=H_j/|U|.
$$

`n_clonotypes` is $H_j$; `fraction_of_unique_biopsy_clonotypes` is $h_j$. Each distinct CDR3 contributes once to the denominator. Only observed recurrence classes are exported. “Private” means $v_i=1$; “public” at threshold two means $v_i\ge2$ in this cohort.

### Global-minimum-depth sensitivity analysis

All biopsies are downsampled to $T_{\mathrm{global}}=\min_s N_s$. Presence and $v_i$ are recalculated at each iteration. Per-subject summaries cover $Q_{s,2},A_{s,2},Q_{s,3},A_{s,3}$, but not their underlying counts. Cohort-level iteration quantities are:

$$
|U|,\qquad H_1,\qquad \sum_{j\ge2}H_j,\qquad
\sum_{j\ge3}H_j,\qquad \frac{\sum_{j\ge2}H_j}{|U|}.
$$

These are summarized using the same formulas as other resampled metrics. The global summary passes `NaN` as its raw reference, so its `raw_value` and `raw_minus_std_median` fields remain missing. Per-subject summaries use the corresponding raw publicness value.

The script explicitly retains **raw publicness as the primary descriptive estimate** and labels global-minimum-depth resampling as a sensitivity analysis, not the primary biological publicness feature.

\newpage

# 3. What it produces

A successful run writes **27 files: 24 plain CSV files, two gzip-compressed CSV files and one JSON configuration file**. The inventories below list the exact filenames. No PDF, figure or Excel workbook is produced by the Python script itself.

## 3.1 Input records and within-subject results

**`input_manifest.csv`**  
Matched input filenames and paths, subject and compartment identifiers, and pair-completeness flags.

**`input_repertoire_audit.csv`**  
One row per input repertoire: retained depth, clonotype count, noncanonical-sequence count and input-contract status. It does not export the numbers of removed input rows or collapsed duplicates.

**`within_subject_overlap_summary_raw.csv`**  
One row per complete pair: raw global overlap, compartment-only measures, shared-clonotype enrichment summaries and cumulative top-fraction measures.

**`within_subject_clonotype_table.csv.gz`**  
One row per subject–clonotype in the paired union: counts, frequencies, presence flags, sharing status, compartment ranks and rank fractions, biopsy stratum, shared log-ratio, overlap category and enrichment flags.

**`within_subject_biopsy_abundance_strata_raw.csv`**  
One row per subject and disjoint stratum; all stratum descriptors from Section 2.3.

**`within_subject_biopsy_top_fraction_detection_raw.csv`**  
One row per subject and cumulative top fraction: clonotype count, PBMC detection fraction, biopsy read mass and fraction of that mass detected in PBMC.

**`within_subject_resampling_summary.csv`**  
Long-format summaries of the eight global pair metrics, including valid iteration counts and all resampling statistics in Section 2.5.

**`within_subject_biopsy_abundance_strata_resampling_summary.csv`**  
Long-format summaries of $D_H$ and $W_H$ for each subject and disjoint stratum.

**`within_subject_top_fraction_resampling_summary.csv`**  
Long-format summaries of $D_{T_t}$ and $W_{T_t}$ for each subject and cumulative top fraction.

**`within_subject_overlap_features.csv`**  
One row per complete pair: target depth, requested iteration count, and raw value, standardized median, 2.5th/97.5th percentiles and CV for each global pair and cumulative-top detection metric. Disjoint-stratum metrics are stored separately.

**`within_subject_biopsy_abundance_strata_features.csv`**  
One row per complete pair: all seven raw stratum descriptors, plus standardized median, 2.5th/97.5th percentiles and CV for $D_H,W_H$. Column names start with `within_stratum_<stratum>_`.

## 3.2 Between-biopsy pair tables

**`between_subject_biopsy_overlap_raw.csv`**  
One row per unordered biopsy pair: raw counts, depths, richness, directional sharing and similarity metrics.

**`between_subject_biopsy_resampling_summary.csv`**  
One row per biopsy pair and selected metric: pair-specific target depth and full resampling summaries for the eight metrics.

\newpage

## 3.3 Between-biopsy matrices and publicness files

The six symmetric matrices are:

```text
biopsy_shared_clonotypes_matrix_raw.csv
biopsy_jaccard_matrix_raw.csv
biopsy_morisita_horn_matrix_raw.csv
biopsy_shared_clonotypes_matrix_std_median.csv
biopsy_jaccard_matrix_std_median.csv
biopsy_morisita_horn_matrix_std_median.csv
```

Their rows and columns are subject identifiers; `_std_median` denotes resampling medians.

**`biopsy_public_clonotypes_raw.csv.gz`**  
One row per distinct biopsy CDR3: subject recurrence, prevalence, subject list, total reads and frequency summaries among positive biopsies.

**`biopsy_publicness_by_subject_raw.csv`**  
Raw per-biopsy publicness counts and fractions at thresholds of two and three subjects.

**`biopsy_publicness_distribution_raw.csv`**  
The cohort recurrence distribution $H_j,h_j$.

**`biopsy_publicness_global_min_depth_sensitivity_by_subject.csv`**  
Per-subject global-minimum-depth summaries of the four publicness fractions.

**`biopsy_publicness_global_min_depth_sensitivity_summary.csv`**  
Cohort-level summaries of unique, private and recurrent sequence counts, and the fraction of unique sequences observed in at least two subjects.

## 3.4 Integration and metadata

**`subject_overlap_features.csv`**  
An outer merge of raw paired results (`within_raw_`), selected standardized pair/top-fraction results (`within_`), disjoint-stratum features (`within_stratum_`), raw publicness (`biopsy_public_raw_`) and standardized publicness sensitivity (`biopsy_public_sensitivity_globalmin_`). Biopsy-only subjects can appear with missing paired-analysis fields; PBMC-only subjects do not enter this integration table.

**`overlap_metric_definitions.csv`**  
A selected metric dictionary with formulas, interpretations and categories; it is not a complete dictionary of every exported column.

**`run_config.json`**  
Script filename, Python/NumPy/pandas versions, paths, column definitions, rank cutoffs, top fractions, resampling settings, seed, standardization rules, subject counts and input-contract flags. This script does not record a source hash or completion timestamp.

# 4. Run and export

The script imports `numpy` and `pandas`. From the directory containing the script:

```bash
python3 "2-analyze_exact_clonotype_overlap(2).py" \
  --data_dir "repSeq" \
  --out_dir "results/02_exact_clonotype_overlap" \
  --n_resamples 200 --seed 20260927
```

Use `--skip_resampling` for raw-only analysis. Resampling-summary files are then written empty; standardized matrices retain missing off-diagonal entries and their fixed diagonal conventions. Filename pattern, delimiter, input columns, stratum cutoffs and cumulative top fractions are configurable. **Existing output filenames are overwritten without an `--overwrite` safeguard; use a fresh output directory.**

To regenerate this PDF from the Markdown with Pandoc and XeLaTeX:

```bash
pandoc Step2_exact_clonotype_overlap_guide_EN.md \
  --pdf-engine=xelatex -o Step2_exact_clonotype_overlap_guide_EN.pdf
```

**Source map:** input: lines 151–277; overlap: 284–366; ranks/enrichment: 373–578; within-pair resampling: 612–707; between-biopsy analysis: 714–786; publicness: 793–918; integration: 925–1035; CLI and file writing: 1083–1289. This guide documents the uploaded script, not a specific patient-data execution.
