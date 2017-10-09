#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
cd $(dirname $0)

cd /usr/local/lib/hyphy
for GENE in gag gp41; do
    for RX in NNRTIs PIs; do
        (
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "$APPDIR/result_data/${GENE}${RX}.meds.tre"; # Select tree file
          echo "$HYPHYOUT/${GENE}${RX}.meds.csv"; # Specify output (.csv) file
          echo "1"; # Select [Universal] code mode
        ) |
        HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/MEDS.bf > $HYPHYOUT/${GENE}${RX}.meds.output.txt 2> /dev/null
        python3 $APPDIR/scripts/medsproc.py ${HYPHYOUT}/${GENE}${RX}.meds.csv 0.05 > ${HYPHYOUT}/${GENE}${RX}.meds.result.csv
        mv messages.log $HYPHYOUT/${GENE}${RX}.meds.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.meds.output.txt created"
        echo "$HYPHYOUT/${GENE}${RX}.meds.csv created"
        echo "$HYPHYOUT/${GENE}${RX}.meds.result.csv created"
    done
done
