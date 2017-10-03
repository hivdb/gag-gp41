#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
cd $(dirname $0)
python3 make-naive-tables.py
python3 run-basic.py
mkdir -p $HYPHYOUT
mkdir -p "/app/report"
