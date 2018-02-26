#! /bin/bash

mkdir -p /app/report
Rscript /app/scripts/make-graphical-summary.r
Rscript /app/scripts/make-naive-studies-figs.r
Rscript /app/scripts/make-naive-indel-figs.r
python3 /app/scripts/generate-report.py > /app/report/README.md
