#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
cd $(dirname $0)

cd /usr/local/lib/hyphy
for GENE in Gag Gp41; do
    for RX in NNRTIs PIs; do
        (
          echo 1; # Select [Universal] code mode
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "all"; # Include all codons
        ) |
        HYPHYMP LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/${GENE}${RX}.pairwise.output.txt 2> /dev/null
        mv messages.log $HYPHYOUT/${GENE}${RX}.pairwise.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.pairwise.output.txt created"
    done
done
