#!/bin/bash
set -euo pipefail

REP="$1"

MODEL="lobster_patch_2.slim"
SCENARIO_TSV="$(pwd)/scenario_table.tsv"

OUTDIR="$(pwd)/results"
mkdir -p "$OUTDIR"
OUTFILE="$(printf "%s/rep%04d.tsv" "$OUTDIR" "$REP")"

BASESEED=100000
MPAS=(block_up block_down alt_A alt_B)

cond_to_id() {
  case "$1" in
    Low)    echo 1 ;;
    Medium) echo 2 ;;
    High)   echo 3 ;;
    *)      echo 99 ;;
  esac
}

mpa_to_id() {
  case "$1" in
    block_up)   echo 1 ;;
    block_down) echo 2 ;;
    alt_A)      echo 3 ;;
    alt_B)      echo 4 ;;
    *)          echo 9 ;;
  esac
}




# cond の読み込み（既存の読み込み処理をそのまま使ってOK）
conds=()
{
  IFS=$'\t' read -r -a header
  cond_col=-1
  for i in "${!header[@]}"; do
    [[ "${header[$i]}" == "cond" ]] && cond_col=$i && break
  done
  [[ $cond_col -lt 0 ]] && exit 1

  while IFS=$'\t' read -r -a cols; do
    [[ ${#cols[@]} -eq 0 ]] && continue
    [[ -z "${cols[0]:-}" ]] && continue
    [[ "${cols[0]}" =~ ^# ]] && continue
    [[ ${#cols[@]} -lt ${#header[@]} ]] && continue

    cond_val="${cols[$cond_col]}"
    [[ -z "$cond_val" ]] && continue
    [[ "$cond_val" =~ ^# ]] && continue
    conds+=("$cond_val")
  done
} < "$SCENARIO_TSV"

# rep ファイルはジョブ開始時に初期化（ヘッダが必要ならここで1回だけ書く）
: > "$OUTFILE"

for cond in "${conds[@]}"; do
  cond_id="$(cond_to_id "$cond")"

  for mpa in "${MPAS[@]}"; do
    mpa_id="$(mpa_to_id "$mpa")"
    seed=$((BASESEED + cond_id * 10000 + mpa_id * 1000 + REP))

    slim \
      -s "$seed" \
      -d COND="\"$cond\"" \
      -d MPA="\"$mpa\"" \
      -d REP="$REP" \
      -d SCENARIO_TSV="\"$SCENARIO_TSV\"" \
      -d OUTFILE="\"$OUTFILE\"" \
      "$MODEL"
  done
done



if [[ ${#conds[@]} -eq 0 ]]; then
  echo "ERROR: scenario_table.tsv から cond を1つも読めませんでした: $SCENARIO_TSV" >&2
  exit 1
fi