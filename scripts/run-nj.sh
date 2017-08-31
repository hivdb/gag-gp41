#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/result_data/hyphy_output"
cd $(dirname $0)

cd /usr/local/lib/hyphy
for GENE in Gag Gp41; do
    for RX in NNRTIs PIs; do
        (
          echo 1; # Select [Distance formulae] computation
          echo 1; # Select [Nucleotide/Protein] data type
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          # echo "$APPDIR/local/dotfasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo 1; # Select [Keep Negative] branch lengths
          echo 7; # Select [TN93] distance formula
          echo "y"; # Save the tree to a file
          echo "$APPDIR/result_data/${GENE}${RX}.tre"; # The filename

        ) |
        HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/NeighborJoining.bf > $HYPHYOUT/${GENE}${RX}.nj.output.txt 2> /dev/null
        mv messages.log $HYPHYOUT/${GENE}${RX}.nj.messages.log
        cp $APPDIR/result_data/${GENE}${RX}.tre $APPDIR/result_data/${GENE}${RX}.meds.tre
        sed -i 's/Pre:/Pre{Pre}:/g' $APPDIR/result_data/${GENE}${RX}.tre
        sed -i 's/Post:/Post{Post}:/g' $APPDIR/result_data/${GENE}${RX}.tre
        sed -i 's/Post:/Post{FG}:/g' $APPDIR/result_data/${GENE}${RX}.meds.tre
        echo "$APPDIR/result_data/${GENE}${RX}.tre created"
        echo "$APPDIR/result_data/${GENE}${RX}.meds.tre created"
    done
done
