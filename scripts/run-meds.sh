#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/internalFiles/hyphyOutput"
INTERNAL_PHYLO="$APPDIR/internalFiles/phylo"
HYPHY="HYPHYMP"
cd $(dirname $0)
mkdir -p $INTERNAL_PHYLO

cd /usr/local/lib/hyphy
for GENE in gag gp41; do
    for RX in NNRTIs PIs; do
        (
          echo "$APPDIR/internalFiles/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "$INTERNAL_PHYLO/${GENE}${RX}.meds.tre"; # Select tree file
          echo "$HYPHYOUT/${GENE}${RX}.meds.csv"; # Specify output (.csv) file
          echo "1"; # Select [Universal] code mode
        ) |
        $HYPHY LIBPATH=`pwd` TemplateBatchFiles/MEDS.bf > $HYPHYOUT/${GENE}${RX}.meds.output.txt 2> /dev/null
        python3 $APPDIR/scripts/medsproc.py ${HYPHYOUT}/${GENE}${RX}.meds.csv 0.05 > ${HYPHYOUT}/${GENE}${RX}.meds.result.csv
        echo "$HYPHYOUT/${GENE}${RX}.meds.output.txt created"
        echo "$HYPHYOUT/${GENE}${RX}.meds.csv created"
        echo "$HYPHYOUT/${GENE}${RX}.meds.result.csv created"
    done
done
