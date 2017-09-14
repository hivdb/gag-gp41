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
for RX in NNRTIs PIs; do
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/data/fasta/Gag${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "0-395"; # MA domain: 1 - 132; (1 - 1) * 3 = 0; 132 * 3 - 1 = 395
    ) |
    HYPHYMP LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/Gag${RX}-MA.pairwise.output.txt 2> /dev/null
    mv messages.log $HYPHYOUT/Gag${RX}-MA.pairwise.messages.log
    echo "$HYPHYOUT/Gag${RX}-MA.pairwise.output.txt created"
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/data/fasta/Gag${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "1089-1499"; # C-terminal domain: 364 - 500; (364 - 1) * 3 = 1089; 500 * 3 - 1 = 1499
    ) |
    HYPHYMP LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/Gag${RX}-CTerminal.pairwise.output.txt 2> /dev/null
    mv messages.log $HYPHYOUT/Gag${RX}-CTerminal.pairwise.messages.log
    echo "$HYPHYOUT/Gag${RX}-CTerminal.pairwise.output.txt created"
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/data/fasta/Gp41${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "582-1034"; # CD domain: 195 - 345; (195 - 1) * 3 = 582; 345 * 3 - 1 = 1034
    ) |
    HYPHYMP LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/Gp41${RX}-CD.pairwise.output.txt 2> /dev/null
    mv messages.log $HYPHYOUT/Gp41${RX}-CD.pairwise.messages.log
    echo "$HYPHYOUT/Gp41${RX}-CD.pairwise.output.txt created"
done
