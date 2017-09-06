#! /bin/bash

mkdir -p /app/report
Rscript /app/scripts/make-graphical-summary.r
python3 /app/scripts/generate-report.py > /app/report/README.md
