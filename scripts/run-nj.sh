#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/internalFiles/hyphyOutput"
INTERNAL_PHYLO="$APPDIR/internalFiles/phylo"
HYPHY="HYPHYMP"
cd $(dirname $0)
mkdir -p $INTERNAL_PHYLO
mkdir -p $HYPHYOUT
mkdir -p $APPDIR/data/phylo

cd /usr/local/lib/hyphy
for GENE in gag gp41; do
    for RX in NNRTIs PIs; do
        (
          echo 1; # Select [Distance formulae] computation
          echo 1; # Select [Nucleotide/Protein] data type
          echo "$APPDIR/internalFiles/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo 1; # Select [Keep Negative] branch lengths
          echo 7; # Select [TN93] distance formula
          echo "y"; # Save the tree to a file
          echo "$INTERNAL_PHYLO/${GENE}${RX}.tre"; # The filename

        ) |
        $HYPHY LIBPATH=`pwd` TemplateBatchFiles/NeighborJoining.bf > $HYPHYOUT/${GENE}${RX}.nj.output.txt 2> /dev/null
        cp $INTERNAL_PHYLO/${GENE}${RX}.tre $INTERNAL_PHYLO/${GENE}${RX}.meds.tre
        sed -i 's/Pre:/Pre{Pre}:/g' $INTERNAL_PHYLO/${GENE}${RX}.tre
        sed -i 's/Post:/Post{Post}:/g' $INTERNAL_PHYLO/${GENE}${RX}.tre
        sed -i 's/Post:/Post{FG}:/g' $INTERNAL_PHYLO/${GENE}${RX}.meds.tre
        echo "$INTERNAL_PHYLO/${GENE}${RX}.tre created"
        echo "$INTERNAL_PHYLO/${GENE}${RX}.meds.tre created"
    done
    (
      echo 1; # Select [Distance formulae] computation
      echo 1; # Select [Nucleotide/Protein] data type
      echo "$APPDIR/data/fasta/${GENE}Aligned.fas"; # Select sequence alignment file
      echo 1; # Select [Keep Negative] branch lengths
      echo 7; # Select [TN93] distance formula
      echo "y"; # Save the tree to a file
      echo "$APPDIR/data/phylo/${GENE}Aligned.tre"; # The filename

    ) |
    $HYPHY LIBPATH=`pwd` TemplateBatchFiles/NeighborJoining.bf > $HYPHYOUT/${GENE}All.nj.output.txt 2> /dev/null
    echo "$APPDIR/data/phylo/${GENE}Aligned.tre created"
done
