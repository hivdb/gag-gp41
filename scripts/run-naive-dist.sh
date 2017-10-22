#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/local/hyphyOutput"
HYPHY="mpirun -np $((`nproc --all`-2)) HYPHYMPI"
cd $(dirname $0)
mkdir -p $HYPHYOUT
mkdir -p $APPDIR/data/phylo

cd /usr/local/lib/hyphy
cp /app/scripts/DistanceMatrix.bf TemplateBatchFiles/DistanceMatrix.bf
for GENE in gag gp41; do
    (
      echo 1; # Select [Distance formulae] computation
      echo 1; # Select [Nucleotide/Protein] data type
      echo "$APPDIR/data/naiveStudies/${GENE}NaiveAligned.fas"; # Select sequence alignment file
      # echo "$APPDIR/data/fasta/${GENE}Aligned.fas"; # Select sequence alignment file
      echo 4; # Select [TAB] matrix format
      echo "$HYPHYOUT/${GENE}NaiveDistanceMatrix.txt"; # save file path
      echo 7; # Select [TN93] distance formula
    ) |
    $HYPHY LIBPATH=`pwd` TemplateBatchFiles/DistanceMatrix.bf > $HYPHYOUT/${GENE}.naivedist.output.txt 2> /dev/null
done
