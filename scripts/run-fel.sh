#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
cd $(dirname $0)

cd /usr/local/lib/hyphy
for GENE in gag gp41; do
    for RX in NNRTIs PIs; do
        (
          echo 1; # Select [Universal] code mode
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "$APPDIR/result_data/${GENE}${RX}.tre"; # Select tree file
          echo 5; # Choose set of [Post] to test for selection
          echo 1; # Use synonymous rate variation
          echo 0.1;
        ) |
        HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/SelectionAnalyses/FEL.bf > $HYPHYOUT/${GENE}${RX}.fel.output.txt 2> /dev/null
        mv messages.log $HYPHYOUT/${GENE}${RX}.fel.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.fel.output.txt created"
    done
done
