#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
cd $(dirname $0)
python3 run.py
mkdir -p $HYPHYOUT
