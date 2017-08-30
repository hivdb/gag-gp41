#! /bin/bash

set -e

APPDIR="/app"
HYPHYOUT="$APPDIR/local/hyphy-output"
cd $(dirname $0)
python3 run.py
mkdir -p $HYPHYOUT

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
        # (
        #   echo 1; # Select [Universal] code mode
        #   echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
        #   echo "$APPDIR/result_data/${GENE}${RX}.tre"; # Select tree file
        #   echo 2; # Choose set of [Post] to test for relaxed selection
        #   echo 2; # Choose set of [Pre] as the reference branches
        #   echo 2; # Chose [Minimal] analysis type
        # ) |
        # HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/SelectionAnalyses/RELAX.bf > $HYPHYOUT/${GENE}${RX}.relax.output.txt 2> /dev/null
        # mv messages.log $HYPHYOUT/${GENE}${RX}.relax.messages.log
        # echo "$HYPHYOUT/${GENE}${RX}.relax.output.txt created"
        (
          echo 1; # Select [Universal] code mode
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "$APPDIR/result_data/${GENE}${RX}.tre"; # Select tree file
          echo 5; # Choose set of [Post] to test for selection
          echo 0.1;
        ) |
        HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/SelectionAnalyses/FEL.bf > $HYPHYOUT/${GENE}${RX}.fel.output.txt 2> /dev/null
        mv messages.log $HYPHYOUT/${GENE}${RX}.fel.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.fel.output.txt created"
        (
          echo "$APPDIR/data/fasta/${GENE}${RX}.aln.fasta.txt"; # Select sequence alignment file
          echo "$APPDIR/result_data/${GENE}${RX}.meds.tre"; # Select tree file
          echo "$HYPHYOUT/${GENE}${RX}.meds.csv"; # Specify output (.csv) file
          echo "1"; # Select [Universal] code mode
        ) |
        HYPHYMP LIBPATH=`pwd` TemplateBatchFiles/MEDS.bf > $HYPHYOUT/${GENE}${RX}.meds.output.txt 2> /dev/null
        python3 $APPDIR/scripts/medsproc.py ${HYPHYOUT}/${GENE}${RX}.meds.csv 0.05 > ${HYPHYOUT}/${GENE}${RX}.meds.result.csv
        mv messages.log $HYPHYOUT/${GENE}${RX}.meds.messages.log
        echo "$HYPHYOUT/${GENE}${RX}.meds.output.txt created"
        echo "$HYPHYOUT/${GENE}${RX}.meds.csv created"
        echo "$HYPHYOUT/${GENE}${RX}.meds.result.csv created"
    done
done
