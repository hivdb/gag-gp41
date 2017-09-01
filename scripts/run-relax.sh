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
          echo "$APPDIR/result_data/${GENE}${RX}.tre"; # Select tree file
          echo 2; # Choose set of [Post] to test for relaxed selection
          echo 1; # Choose set of [Pre] as the reference branches
          echo 2; # Chose [Minimal] analysis type
        ) |
        HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/SelectionAnalyses/RELAX.bf > $HYPHYOUT/${GENE}${RX}.relax.output.txt 2> /dev/null
        mv messages.log $HYPHYOUT/${GENE}${RX}.relax.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.relax.output.txt created"
    done
done
