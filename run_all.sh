#!/usr/bin/env bash
# Run the full chain after configuring config/config.json (see README).
# Usage: ./run_all.sh            (pipeline + analyses + figures + supplementary data)
#        ./run_all.sh --no-pipeline   (skip Steps 1-9 when results/ already exists)
set -euo pipefail
cd "$(dirname "$0")"
CFG="${ENDO_CONFIG:-config/config.json}"; export ENDO_CONFIG="$CFG"
get() { python3 -c "import json,os;c=json.load(open('$CFG'));v=c['$1'];print(v if os.path.isabs(v) else os.path.join(c['project_root'],v))"; }
D=$(get repseq_dir); R=$(get results_dir); CLIN=$(get clinical_table); OCA="$(get project_root)/oca_plus_new.txt"
if [[ "${1:-}" != "--no-pipeline" ]]; then
  P=pipeline
  python3 $P/01_step1_repertoire_characteristics.py --data_dir "$D" --out_dir "$R/01_repertoire_characterization"
  python3 $P/02_step2_exact_clonotype_overlap.py --data_dir "$D" --out_dir "$R/02_exact_clonotype_overlap_500" --n_resamples 500 --seed 20260927
  python3 $P/03_stepR2_repertoire_figures_tables.py --step1 "$R/01_repertoire_characterization" --step2 "$R/02_exact_clonotype_overlap_500" --out_dir "$R/publication_reporting_v1" --histology_xlsx "$CLIN"
  python3 $P/04_step3_subject_master_dataset.py --histology_file "$CLIN" --repseq_dir "$D" --step1 "$R/01_repertoire_characterization" --step2 "$R/02_exact_clonotype_overlap_500" --oca_file "$OCA" --out_dir "$R/03_subject_master"
  A="$R/03_subject_master/03_analysis_dataset_subject_level.csv"
  python3 $P/05_step4_freeze_immune_panel.py --analysis_dataset "$A" --out_dir "$R/04_immune_feature_panel"
  python3 $P/06_step5_freeze_genomic_panel.py --analysis_dataset "$A" --out_dir "$R/05_genomic_feature_panel"
  python3 $P/07_step6_immune_genomic_associations.py --analysis_dataset "$A" --immune_panel_dir "$R/04_immune_feature_panel" --genomic_panel_dir "$R/05_genomic_feature_panel" --out_dir "$R/06_immune_genomic_associations"
  python3 $P/08_step7_integrated_audit_reporting.py --step3 "$R/03_subject_master" --step4 "$R/04_immune_feature_panel" --step5 "$R/05_genomic_feature_panel" --step6 "$R/06_immune_genomic_associations" --out_dir "$R/2_figures_tables_steps3_6"
  python3 $P/09_step8_histology_associations.py --step3 "$R/03_subject_master" --step4 "$R/04_immune_feature_panel" --step5 "$R/05_genomic_feature_panel" --out_dir "$R/08_histology_associations"
  python3 $P/10_step9_histology_publication_reporting.py --step8 "$R/08_histology_associations" --out_dir "$R/09_histology_reporting"
fi
for s in a1_recompute_repertoire_metrics a2_depth_standardisation a3_clonotype_detection_model a4_tmb_sensitivity; do python3 analysis/$s.py; done
for s in fig1_fig2 fig3_figS6 figS1 figS2_montage figS3_S5; do python3 figures/$s.py; done
python3 supplementary/s1_build_patient_table.py
python3 supplementary/s2_build_supplementary_data.py
