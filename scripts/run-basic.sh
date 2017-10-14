#! /bin/bash

set -e

APPDIR="/app"
cd $(dirname $0)
python3 analysis-naive-studies.py
python3 run-basic.py
mkdir -p $HYPHYOUT
mkdir -p "/app/report"
