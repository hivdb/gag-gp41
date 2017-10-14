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
          echo 1; # Select [Universal] code mode
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "$INTERNAL_PHYLO/${GENE}${RX}.tre"; # Select tree file
          echo 5; # Choose set of [Post] to test for selection
          echo 1; # Use synonymous rate variation
          echo 0.1;
        ) |
        $HYPHY LIBPATH=`pwd` TemplateBatchFiles/SelectionAnalyses/FEL.bf > $HYPHYOUT/${GENE}${RX}.fel.output.txt 2> /dev/null
        echo "$HYPHYOUT/${GENE}${RX}.fel.output.txt created"
    done
done
