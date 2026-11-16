#!/usr/bin/env bash
# Pile up the classifier genes (gene_regions.tsv, hg19) in one DNA BAM.
# Usage: run_pileups.sh <sample_DNA.bam> <output_dir>
# Writes <output_dir>/<sample>_<GENE>.json. Probe positions: POLE exonuclease-domain hotspots
# (c.857 P286, c.833 V278, c.1231 V411, c.1346 S459 region, c.1373 A456 region; LOVD, NC_000012.11).
set -euo pipefail
bam="$1"; out="$2"; here="$(cd "$(dirname "$0")" && pwd)"; mkdir -p "$out"
s="$(basename "$bam" .bam)"
tail -n +2 "$here/gene_regions.tsv" | while IFS=$'\t' read -r gene chrom start end; do
  python3 "$here/pileup.py" "$bam" "$chrom" "$start" "$end" 133253184 133253208 133250289 133250174 133249850 > "$out/${s}_${gene}.json"
done
