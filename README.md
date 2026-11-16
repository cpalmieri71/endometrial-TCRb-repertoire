# Paired tumor-biopsy and PBMC TCRβ repertoires in endometrial cancer

Code for the manuscript *(title, authors and DOI to be added on acceptance)*.
It reproduces every number, figure and Supplementary Data file from the processed
clonotype tables, the clinical table and the Oncomine Comprehensive Assay Plus (OCA Plus) results.

## Repository layout

| folder | content |
|---|---|
| `pipeline/` | Steps 1–9: the frozen analysis chain (ten scripts, run in the order of their prefix). Each step writes a `run_config.json` with inputs, parameters, seeds and library versions. `SHA256SUMS` lists their checksums. |
| `analysis/` | Analyses added during manuscript revision: A1 recomputation of the repertoire metrics, A2 depth/coverage standardization (Supplementary Data 3), A3 clonotype-level model of peripheral detection, A4 TMB sensitivity analyses. |
| `figures/` | Main Figures 1–3 and Supplementary Figures S1–S6. |
| `supplementary/` | S1 patient-level table; S2 assembly and de-identification of Supplementary Data 1–15. |
| `molecular/` | Standard-library BAM tools used to verify POLE/TP53 variants on the OCA Plus DNA BAMs, and the WHO 2020 classifier. |
| `data/` | `sra_accessions.csv` (BioProject PRJNA1474700). |
| `config/` | `config.example.json` and templates for the two private files. |
| `docs/` | Detailed guides for Steps 1 and 2. |

## Data

* Processed MiXCR clonotype tables (26 files, `<id>_PBMC.tsv` / `<id>_biopsia.tsv`):
  Zenodo, https://doi.org/10.5281/zenodo.18246999
* Raw reads: NCBI SRA, BioProject PRJNA1474700 (`data/sra_accessions.csv`).
* Clinical table (`tabella_endo.xlsx`) and OCA Plus summary (`oca_plus_new.txt`): patient-level
  data, available from the corresponding author under a data-use agreement.

## Private files (never committed)

Published outputs identify patients only as P1–P13. The link to internal identifiers is kept in
`private/`, which is excluded by `.gitignore`:

* `private/patient_key.csv` — columns `internal_id,patient_code`
  (template: `config/templates/patient_key.example.csv`);
* `private/curated_annotations.csv` — per-patient curated fields not produced by the pipeline:
  MSI status (assay), POLE exonuclease-domain variant, TP53 variant, WHO 2020 class, regional nodal
  status and nodes examined (template: `config/templates/curated_annotations.example.csv`).
  The molecular fields come from the Ion Reporter exports and the BAM verification in `molecular/`.

Pipeline outputs under `results/` contain internal identifiers in several audit tables; only the
de-identified Supplementary Data produced by `supplementary/s2_build_supplementary_data.py` are
intended for publication. Do not use the `--show_subject_ids` option of Step R2 for public outputs.

## Running

```bash
cp config/config.example.json config/config.json     # edit project_root and paths
pip install -r requirements.txt
./run_all.sh                 # Steps 1–9, analyses, figures, Supplementary Data
./run_all.sh --no-pipeline   # only the downstream part, when results/ already exists
```

All scripts read the configuration from `config/config.json`, or from the file given by
`--config` or the `ENDO_CONFIG` environment variable. Outputs of `analysis/`, `figures/` and
`supplementary/` go to `<output_dir>/analysis`, `/figures`, `/supplementary`, `/supplementary_data`.

Run order and dependencies:

```
repSeq ─┬─ Step 1 ─┬──────────────── Step R2 (repertoire figures/tables, Supplementary Data 1–2 sources)
        └─ Step 2 ─┘
Steps 1, 2 + clinical table + OCA Plus ── Step 3 ─┬─ Step 4 ─┐
                                                  └─ Step 5 ─┼─ Step 6 ── Step 7 (reporting)
                                                             └─ Step 8 ── Step 9 (reporting)
A1 (repSeq, Step 1) → A2 → figure S1;   A3 (Step 2) → Figure 2;   A4 (Step 6, A1)
S1 (Step 3 + private files) → S2 (Steps R2, 7, 9, A2, S1) → Supplementary Data 1–15
```

Random seeds: Step 2 resampling 20260927; A2 rarefaction 20261110 and coverage
standardization 7; A4 permutations 1. In A2 the order of random draws is part of the result.

### Molecular verification (optional)

```bash
molecular/run_pileups.sh sample_DNA.bam out_dir/            # POLE, TP53, MMR genes, POLD1 (hg19)
python molecular/crosscons.py chr17 7577000 7579000 a.bam b.bam ...   # variants with VAF > 50%
python molecular/classify_who2020.py private/curated_annotations.csv --check
```
Tumor-only data: germline and somatic variants cannot be distinguished. Indels in homopolymers
are unreliable on Ion Torrent data and were not used for calls.

## Environment

Pipeline Steps 1–9: Python 3.9.22, NumPy 2.0.2, pandas 2.3.3, SciPy 1.13.1, matplotlib 3.9.4,
powerlaw 2.0.0 (MiXCR 4.1.2 upstream). Downstream scripts: Python 3.13, openpyxl 3.1.
Figures use the FreeSans font when installed (falls back to the matplotlib default).

## License

*(to be chosen by the authors before the repository is made public)*
