#! /bin/bash

set -e

APPDIR="/app"
cd $(dirname $0)
python3 run-basic.py
mkdir -p "/app/report"
