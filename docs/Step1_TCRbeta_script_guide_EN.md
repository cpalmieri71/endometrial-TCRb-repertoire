---
title-meta: "Step 1: TCRβ repertoire characterization"
author: []
date: ""
lang: en-GB
papersize: a4
fontsize: 10pt
geometry:
  - top=17mm
  - bottom=18mm
  - left=19mm
  - right=19mm
mainfont: DejaVu Serif
sansfont: DejaVu Sans
monofont: DejaVu Sans Mono
monofontoptions:
  - Scale=0.80
colorlinks: true
header-includes:
  - |
    \usepackage{titlesec}
    \usepackage{fancyhdr}
    \usepackage{microtype}
    \titleformat{\section}{\large\sffamily\bfseries}{\thesection}{0.5em}{}
    \titleformat{\subsection}{\normalsize\sffamily\bfseries}{\thesubsection}{0.5em}{}
    \titlespacing*{\section}{0pt}{10pt}{5pt}
    \titlespacing*{\subsection}{0pt}{8pt}{4pt}
    \setlength{\parskip}{4pt}
    \setlength{\emergencystretch}{3em}
    \AtBeginDocument{\setlength{\abovedisplayskip}{5pt}\setlength{\belowdisplayskip}{5pt}}
    \pagestyle{fancy}
    \fancyhf{}
    \fancyhead[L]{\small\sffamily Step 1 | Script reference}
    \fancyhead[R]{\small\sffamily v2.0.0}
    \fancyfoot[R]{\small\sffamily\thepage}
    \renewcommand{\headrulewidth}{0.3pt}
    \setlength{\headheight}{14pt}
---

# Step 1: TCRβ repertoire characterization

**Script reference · v2.0.0**  
Source: `1-build_repertoire_characteristics_v2(4).py`.

# 1. What it does

The script characterizes each input TCRβ repertoire independently, calculates read-based sampling and repertoire-architecture measures, optionally fits abundance-tail models, and summarizes selected measures by PBMC or BIOPSY compartment.

It does **not** perform productivity filtering, rarefaction, biological sample exclusions, clonotype-overlap analysis, or inferential PBMC–biopsy comparisons. PBMC and biopsy are separate biological repertoires, not technical replicates.

# 2. How it works

## 2.1 Input processing

The script scans the input directory, without recursion, for names such as `1_PBMC.tsv`, `1_biopsia.tsv`, or `1_biopsy.tsv` (case-insensitive). The numeric prefix identifies the subject; both biopsy labels become `BIOPSY`. Duplicate subject–compartment inputs cause an error. Incomplete pairs are flagged in `pair_complete`, but retained.

By default, `aaSeqCDR3` identifies clonotypes and `readCount` supplies abundances. Rows with missing or blank identifiers, or non-finite, non-numeric, zero or negative counts, are removed. Identifiers are stripped of surrounding whitespace; identical strings within each file are collapsed by summing counts. No V/J matching or productivity check is added. An empty cleaned repertoire stops execution.

For the retained clonotypes, let $c_i$ be the collapsed count, $S$ the number of distinct identifiers, $N=\sum_{i=1}^{S}c_i$, and $p_i=c_i/N$. All logarithms are natural. Frequencies are recalculated from counts; an input `readFraction` column is not used.

## 2.2 Sampling measures

| Output column | Calculation |
|:--|:--|
| `depth` | $N=\sum_i c_i$ |
| `observed_richness` | $S$, the number of retained distinct clonotypes |
| `singleton_count` | $f_1=\sum_i\mathbf{1}(c_i=1)$ |
| `doubleton_count` | $f_2=\sum_i\mathbf{1}(c_i=2)$ |
| `goods_coverage` | $\widehat C=1-f_1/N$ |
| `chao1_richness` | $\widehat S_{\mathrm{Chao1}}$, defined below |

$$
\widehat S_{\mathrm{Chao1}}=
\begin{cases}
S+\dfrac{f_1^2}{2f_2}, & f_2>0,\\[4pt]
S+\dfrac{f_1(f_1-1)}{2}, & f_2=0.
\end{cases}
$$

Singletons and doubletons are defined after collapsing. If any collapsed count differs from its nearest integer by more than $10^{-8}$, these counts, Chao1 and Good's coverage are missing, and discrete tail fitting is skipped. Other repertoire measures are still calculated. The `counts_integer_like` flag additionally checks the total count against the same tolerance.

Good's coverage and Chao1 are read-sampling descriptors. Coverage is not the fraction of biological clonotypes recovered; singleton/doubleton read counts do not identify one-cell/two-cell clones.

\newpage

## 2.3 Repertoire-architecture measures

**Hill diversity.** The default outputs are `hill_D1` and `hill_D2`:

$$
D_1=\exp\!\left(-\sum_{i=1}^{S}p_i\ln p_i\right),
\qquad
D_2=\frac{1}{\sum_{i=1}^{S}p_i^2}.
$$

Additional nonnegative orders can be requested using `--hill_orders`. For $q\ne1$, $D_q=(\sum_i p_i^q)^{1/(1-q)}$. $D_1$ and $D_2$ are always included; $D_0$ is represented only by `observed_richness`. These are effective numbers of clonotypes, calculated without depth correction.

**Gini coefficient** (`gini_frequency`). With frequencies sorted in ascending order, $p_{(1)}\le\cdots\le p_{(S)}$:

$$
G=\frac{2\sum_{i=1}^{S}i\,p_{(i)}}{S\sum_{i=1}^{S}p_{(i)}}-\frac{S+1}{S}.
$$

This is the uncorrected Gini inequality coefficient over observed clonotypes, not the Gini–Simpson index.

**Top-N abundance.** With descending frequencies $p_{[1]}\ge\cdots\ge p_{[S]}$:

$$
T_n=\sum_{i=1}^{\min(n,S)}p_{[i]}.
$$

Defaults are `top10_frequency` and `top100_frequency`. `--top_n` changes the requested cutoffs, but Top-10 is always included. Values are read fractions, not percentages or cell fractions.

**Fraction of clonotypes reaching a cumulative abundance target.**

$$
K_t=\min\!\left\{k:\sum_{i=1}^{k}p_{[i]}\ge t\right\},
\qquad F_t=\frac{K_t}{S}.
$$

The default is $t=0.50$, exported as `fraction_clones_to_50pct`. Other targets in $(0,1)$ can be set with `--cumulative_targets`. Canonical tables export $F_t$, not the redundant raw count $K_t$.

## 2.4 Optional abundance-tail fitting

Enabled by default; disabled with `--skip_tail_fitting`. The script delegates fitting to `powerlaw`, using one integer-valued read count per clonotype, `discrete=True`, and `discrete_normalization="round"`. The package selects $x_{\min}$ by power-law Kolmogorov–Smirnov distance minimization; power-law, truncated power-law and lognormal models are evaluated on the same tail $T=\{i:c_i\ge x_{\min}\}$.

The exported tail descriptors are $x_{\min}$, $x_{\min}/N$, $n_T=|T|$, $n_T/S$, and $\sum_{i\in T}c_i/N$. Model outputs include fitted parameters, power-law exponent uncertainty, the power-law KS distance, and full-tail log-likelihoods:

$$
\ell_M=\sum_{i\in T}\log P_M(c_i\mid c_i\ge x_{\min}),
\qquad
\mathrm{AICc}_M=2k_M-2\ell_M+\frac{2k_M(k_M+1)}{n_T-k_M-1}.
$$

Here $k_M=1$ for power law and $2$ for each alternative. AICc is missing if $n_T\le k_M+1$ or the log-likelihood is invalid. A lowest-AICc model is named only when all three AICc values are finite. The calculation is conditional on the selected cutoff; cutoff-selection uncertainty is not penalized. Package likelihood-ratio comparisons return $R=\ell_{\mathrm{PL}}-\ell_{\mathrm{alternative}}$ and nominal $p$-values, using nested comparison for truncated power law and non-nested comparison for lognormal.

With `--min_tail_size 50` (default), fewer than 50 total clonotypes prevents fitting; fewer than 50 selected tail clonotypes only triggers a descriptive flag. `tail_fit_success` requires finite model parameters and full-tail log-likelihoods; it is not a goodness-of-fit test. Warnings and incomplete comparisons are recorded. `--exact_discrete_tail_fit` disables the fast exponent estimate, but does not guarantee exact discretization for every model.

\newpage

## 2.5 Compartment summaries and compatibility export

The six reporting measures are **depth, observed richness, Good's coverage, D₁, D₂ and Top-10**. For each measure and compartment, finite values are summarized by valid $n$, median, first and third quartiles, mean, sample standard deviation (`ddof=1`), minimum and maximum. The formatted table reports median [Q1–Q3], plus repertoire and subject counts. All available repertoires contribute, including incomplete pairs; this is not a paired comparison.

Canonical tables avoid duplicate or derived aliases such as Shannon entropy, inverse Simpson, D₀, evenness and clonality. `--write_legacy_compat` additionally writes separate tables under `legacy_compat/`, restoring those fields, additional Hill orders and Top-N cutoffs, and cumulative counts/fractions at 25%, 50% and 75%. These compatibility tables do not replace the canonical outputs.

# 3. What it produces

A successful default run writes **seven CSV files and one JSON file** to `--out_dir`.

| File | Contents |
|:--|:--|
| `input_manifest.csv` | Subject, compartment, filename, absolute input path and pair-completeness flag. |
| `repertoire_characteristics_full.csv` | One row per repertoire: identifiers, cleaning diagnostics, all configured metrics and tail-fit results/status. |
| `repertoire_characteristics_table1_by_repertoire.csv` | Identifiers and the six reporting measures for each repertoire. |
| `table1_summary_by_compartment_numeric.csv` | Numeric compartment summaries of the six reporting measures. |
| `table1_paper_ready.csv` | Formatted PBMC/BIOPSY comparison of descriptive summaries; no hypothesis tests. |
| `repertoire_metric_definitions.csv` | Metric formulas, interpretations, categories and reporting roles. |
| `metric_selection_and_redundancy.csv` | Retained representations, mathematical relationships and selection decisions. |
| `run_config.json` | Settings, software versions, script SHA-256, completion time, input counts, ignored TSV filenames and analysis flags. |

No plots or Excel workbook are generated. Cleaning diagnostics include the numbers of raw, invalid, valid pre-collapse and duplicate-collapsed rows. When tail fitting is skipped, tail columns remain present with missing/default values and an explanatory status message.

## Run command

The script specifies Python ≥3.9, `numpy` and `pandas`; tail fitting additionally requires `powerlaw`.

```bash
python3 "1-build_repertoire_characteristics_v2(4).py" \
  --data_dir "repSeq" \
  --out_dir "results/01_repertoire_characterization_v2"
```

Append `--skip_tail_fitting` to omit tail models. Use a new output directory, or add `--overwrite` to allow rewriting an existing non-empty directory. Input and output directories must differ. Filename matching, delimiter and column names are configurable through `--pattern`, `--file_sep`, `--clone_col` and `--count_col`.
