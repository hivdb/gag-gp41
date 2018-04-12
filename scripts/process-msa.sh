#! /bin/bash

set -e

APPDIR="/app"

for GENE in gag gp41; do
    python $APPDIR/scripts/stripins.py $GENE CONSENSUS_B $APPDIR/data/fasta/${GENE}MSA.fas $APPDIR/data/fasta/${GENE}Aligned.fas $APPDIR/local/${GENE}Insertions.csv
    sed -i 's/Header,Gene,/Accession,Gene,PID,TimePoint,/' $APPDIR/local/${GENE}Insertions.csv
    sed -Ei 's/^(.+?)\|([0-9]+)_.+?_(Pre|Post),([a-z]+),/\1,\4,\2,\3,/g' $APPDIR/local/${GENE}Insertions.csv
done

cat $APPDIR/local/gagInsertions.csv > $APPDIR/data/insertions.csv
tail -n +2 $APPDIR/local/gp41Insertions.csv >> $APPDIR/data/insertions.csv

