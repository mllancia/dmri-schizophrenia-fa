#!/bin/bash
set -u
BUCKET="s3://openneuro.org/ds000030"
mkdir -p raw logs
tail -n +2 cohort.tsv | cut -f1 | while read -r sub; do
  dest="raw/$sub/dwi"
  n=$(ls "$dest" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$n" -eq 4 ]; then echo "skip $sub"; continue; fi
  mkdir -p "$dest"
  if aws s3 cp --no-sign-request --recursive "$BUCKET/$sub/dwi/" "$dest/" \
       >> logs/download.log 2>&1; then echo "ok   $sub"; else echo "FAIL $sub"; fi
done
