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
      echo "396-1088"; # CA domain: 133 - 363; (133 - 1) * 3 = 396; 363 * 3 - 1 = 1088
    ) |
    HYPHYMP LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/Gag${RX}-CA.pairwise.output.txt 2> /dev/null
    mv messages.log $HYPHYOUT/Gag${RX}-CA.pairwise.messages.log
    echo "$HYPHYOUT/Gag${RX}-CA.pairwise.output.txt created"
    (
      echo 1; # Select [Universal] code mode
      echo "$APPDIR/data/fasta/Gp41${RX}.aln.fasta.txt"; # Select sequence alignment file
      echo "582-1034"; # CD domain: 195 - 345; (195 - 1) * 3 = 582; 345 * 3 - 1 = 1034
    ) |
    HYPHYMP LIBPATH=`pwd` $APPDIR/scripts/pairwise-estimator-dnds.bf > $HYPHYOUT/Gp41${RX}-CD.pairwise.output.txt 2> /dev/null
    mv messages.log $HYPHYOUT/Gp41${RX}-CD.pairwise.messages.log
    echo "$HYPHYOUT/Gp41${RX}-CD.pairwise.output.txt created"
done
