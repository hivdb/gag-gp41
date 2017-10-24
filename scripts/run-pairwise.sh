#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/internalFiles/hyphyOutput"
HYPHY="HYPHYMP"
cd $(dirname $0)

cd /usr/local/lib/hyphy
for GENE in gag gp41; do
    for RX in NNRTIs PIs; do
        (
          echo 1; # Select [Universal] code mode
          echo "$APPDIR/internalFiles/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "all"; # Include all codons
        ) |
        $HYPHY LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/${GENE}${RX}.pairwise.output.txt 2> /dev/null
        echo "$HYPHYOUT/${GENE}${RX}.pairwise.output.txt created"
    done
done
for RX in NNRTIs PIs; do
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/internalFiles/fasta/gag${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "0-395"; # MA domain: 1 - 132; (1 - 1) * 3 = 0; 132 * 3 - 1 = 395
    ) |
    $HYPHY LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/gag${RX}-MA.pairwise.output.txt 2> /dev/null
    echo "$HYPHYOUT/gag${RX}-MA.pairwise.output.txt created"
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/internalFiles/fasta/gag${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "1089-1499"; # C-terminal domain: 364 - 500; (364 - 1) * 3 = 1089; 500 * 3 - 1 = 1499
    ) |
    $HYPHY LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/gag${RX}-CTerminal.pairwise.output.txt 2> /dev/null
    echo "$HYPHYOUT/gag${RX}-CTerminal.pairwise.output.txt created"
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/internalFiles/fasta/gp41${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "582-1034"; # CD domain: 195 - 345; (195 - 1) * 3 = 582; 345 * 3 - 1 = 1034
    ) |
    $HYPHY LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/gp41${RX}-CD.pairwise.output.txt 2> /dev/null
    echo "$HYPHYOUT/gp41${RX}-CD.pairwise.output.txt created"
done
