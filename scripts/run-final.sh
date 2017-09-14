#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
CLEANOUT="$APPDIR/result_data/hyphy_cleaned_output"
cd $(dirname $0)
mkdir -p $CLEANOUT

for RX in NNRTIs PIs; do
    for GENE in Gag Gp41; do
        python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/${GENE}${RX}.pairwise.output.txt > $CLEANOUT/${GENE}${RX}.pairwise.tsv
        echo "$CLEANOUT/${GENE}${RX}.pairwise.tsv created"

        python3 $APPDIR/scripts/clean-fel.py $HYPHYOUT/${GENE}${RX}.fel.output.txt > $CLEANOUT/${GENE}${RX}.fel.tsv
        echo "$CLEANOUT/${GENE}${RX}.fel.tsv created"
    done
    python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/Gag${RX}-MA.pairwise.output.txt > $CLEANOUT/Gag${RX}-MA.pairwise.tsv
    echo "$CLEANOUT/Gag${RX}-MA.pairwise.tsv created"
    python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/Gag${RX}-CTerminal.pairwise.output.txt > $CLEANOUT/Gag${RX}-CTerminal.pairwise.tsv
    echo "$CLEANOUT/Gag${RX}-CTerminal.pairwise.tsv created"
    python3 $APPDIR/scripts/clean-pairwise.py $HYPHYOUT/Gp41${RX}-CD.pairwise.output.txt > $CLEANOUT/Gp41${RX}-CD.pairwise.tsv
    echo "$CLEANOUT/Gp41${RX}-CD.pairwise.tsv created"
done
