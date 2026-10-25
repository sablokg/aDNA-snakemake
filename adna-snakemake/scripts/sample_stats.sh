#!/usr/bin/env bash
# Per-sample summary line.
# Args: ID LIB FQ1 MERGED_FQ FILT_BAM DEDUP_BAM FINAL_BAM MT X Y OUT
set -euo pipefail
ID=$1; LIB=$2; FQ1=$3; MERGED=$4; FILT=$5; DEDUP=$6; FINAL=$7; MT=$8; X=$9; Y=${10}; OUT=${11}

count_fq() { gzip -dc "$1" | awk 'END{printf "%d", NR/4}'; }
RAW=$(count_fq "$FQ1")
MRG=$(count_fq "$MERGED")
Q=$(samtools view -c "$FILT")
U=$(samtools view -c "$DEDUP")

ENDO=$(awk -v a="$Q" -v b="$MRG" 'BEGIN{printf "%.2f", b ? 100*a/b : 0}')
DUPR=$(awk -v a="$Q" -v b="$U" 'BEGIN{printf "%.3f", a ? 1-b/a : 0}')
MEANLEN=$(samtools view "$DEDUP" | awk '{s+=length($10); n++} END{printf "%.1f", n ? s/n : 0}')

# genome-wide mean depth and 1x breadth over all reference contigs
read -r DEPTH COV1X < <(samtools coverage "$FINAL" | awk 'NR>1{L=$3-$2+1; cb+=$5; dp+=$7*L; len+=L}
                         END{printf "%.4f %.2f\n", len ? dp/len : 0, len ? 100*cb/len : 0}')

MTD=$(samtools coverage -r "$MT" "$FINAL" | awk 'NR>1{printf "%.2f", $7}')
MTD=${MTD:-NA}

# Ry (Skoglund et al. 2013): XY > 0.075, XX < 0.016
IDX=$(samtools idxstats "$FINAL")
NX=$(awk -v c="$X" '$1==c || $1=="chr"c {s+=$3} END{print s+0}' <<<"$IDX")
NY=$(awk -v c="$Y" '$1==c || $1=="chr"c {s+=$3} END{print s+0}' <<<"$IDX")
RY=$(awk -v x="$NX" -v y="$NY" 'BEGIN{printf "%.4f", (x+y) ? y/(x+y) : 0}')
SEX=$(awk -v r="$RY" -v n="$((NX+NY))" 'BEGIN{ if (n<100) print "NA(low)";
        else if (r>0.075) print "XY"; else if (r<0.016) print "XX"; else print "ambiguous" }')

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$ID" "$LIB" "$RAW" "$MRG" "$Q" "$ENDO" "$U" "$DUPR" "$MEANLEN" "$DEPTH" "$COV1X" "$RY" "$SEX" "$MTD" > "$OUT"
