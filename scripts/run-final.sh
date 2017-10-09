#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
CLEANOUT="$APPDIR/result_data/hyphy_cleaned_output"
cd $(dirname $0)
mkdir -p $CLEANOUT

for RX in NNRTIs PIs; do
    for GENE in gag gp41; do
        python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/${GENE}${RX}.pairwise.output.txt > $CLEANOUT/${GENE}${RX}.pairwise.tsv
        echo "$CLEANOUT/${GENE}${RX}.pairwise.tsv created"

        python3 $APPDIR/scripts/clean-fel.py $HYPHYOUT/${GENE}${RX}.fel.output.txt > $CLEANOUT/${GENE}${RX}.fel.tsv
        echo "$CLEANOUT/${GENE}${RX}.fel.tsv created"
    done
    python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/gag${RX}-MA.pairwise.output.txt > $CLEANOUT/gag${RX}-MA.pairwise.tsv
    echo "$CLEANOUT/gag${RX}-MA.pairwise.tsv created"
    python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/gag${RX}-CTerminal.pairwise.output.txt > $CLEANOUT/gag${RX}-CTerminal.pairwise.tsv
    echo "$CLEANOUT/gag${RX}-CTerminal.pairwise.tsv created"
    python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/gp41${RX}-CD.pairwise.output.txt > $CLEANOUT/gp41${RX}-CD.pairwise.tsv
    echo "$CLEANOUT/gp41${RX}-CD.pairwise.tsv created"
done
