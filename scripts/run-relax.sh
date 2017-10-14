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
          echo 2; # Choose set of [Post] to test for relaxed selection
          echo 1; # Choose set of [Pre] as the reference branches
          echo 2; # Chose [Minimal] analysis type
        ) |
        $HYPHY LIBPATH=`pwd` TemplateBatchFiles/SelectionAnalyses/RELAX.bf > $HYPHYOUT/${GENE}${RX}.relax.output.txt 2> /dev/null
        mv messages.log $HYPHYOUT/${GENE}${RX}.relax.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.relax.output.txt created"
    done
done
