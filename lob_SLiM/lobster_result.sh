#!/bin/bash
set -euo pipefail

DIR="${1:-results}"

SUFFIXES=(
  ".tsv"
  ".tsv.natal.tsv"
  ".tsv.popN.tsv"
  ".tsv.sizebin.tsv"
  ".tsv.trap.tsv"
)

mkdir -p "${DIR}/aggregate"

for suf in "${SUFFIXES[@]}"; do
  out="${DIR}/aggregate/all${suf}"
  first="${DIR}/rep0001${suf}"

  if [[ ! -f "$first" ]]; then
    echo "ERROR: missing file: $first" >&2
    exit 1
  fi

  head -n 1 "$first" > "$out"

  for rep in $(seq 1 200); do
    rep4=$(printf "%04d" "$rep")
    f="${DIR}/rep${rep4}${suf}"
    if [[ ! -f "$f" ]]; then
      echo "WARN: missing $f (skip)" >&2
      continue
    fi
    tail -n +2 "$f" >> "$out"
  done

  echo "Wrote: $out"
done
