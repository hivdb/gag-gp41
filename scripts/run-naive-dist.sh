#! /bin/bash

set -e

APPDIR="/app"
cd $(dirname $0)
python3 calc-naive-distances.py
Rscript analyze-naive-diversify.r
