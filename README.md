# The Data and their analysis pipeline of "Evolution of *gag* and *gp41* in Patients Receiving Ritonavir-Boosted Protease Inhibitors"

## Data Description

```
.
+-- data
|   +-- fasta
|   |   +-- *.aln.fasta.txt
|   +-- consensus.csv
|   +-- *NaiveSequences.csv
|   +-- patients.csv
|   +-- samples.csv
|   +-- sequences.csv
+-- result_data
|   +-- hyphy_output
|   |   +-- *.fel.output.txt
|   |   +-- *.meds.csv
|   |   +-- *.meds.output.txt
|   |   +-- *.meds.result.csv
|   |   +-- *.nj.output.txt
|   |   +-- *.pairwise.output.txt
|   |   +-- *.relax.output.txt
|   +-- hyphy_cleaned_output
|   |   +-- *.fel.tsv
|   |   +-- *.pairwise.tsv
|   +-- *AAChangesByPosWPrev.csv
|   +-- *CodonChangesByPt.csv
|   +-- MutPrevalence.csv
|   +-- *.meds.tre
|   +-- *.tre
+-- report
|   +-- README.md
|   +-- gag-mutations.png
|   +-- gp41-mutations.png
+-- scripts
    +-- analysis_functions.py
    +-- codonutils.py
    +-- data_reader.py
    +-- run-basic.py
    +-- run-basic.sh
    +-- run-nj.sh
    +-- run-fel.sh
    +-- run-meds.sh
    +-- medsproc.py
    +-- run-pairwise.sh
    +-- pairwise-estimator-dnds.bf
    +-- run-relax.sh
    +-- run-final.sh
    +-- clean-fel.py
    +-- clean-pairwise.py
    +-- generate-report.sh
    +-- generate-report.py
    +-- make-graphical-summary.r
```

There are three folders in this repository which store the data. The `data`
folder stores unprocessed raw data with one exception that the `fasta` folder
stores aligned sequence results from `sequences.csv`:

- `data/fasta/*.aln.fasta.txt` files: aligned sequences which were grouped by
  their gene and Rx.
- `data/consensus.csv`: consensus/wild type Gag, gp41, PR and RT genes which
  were used by the scripts.
- `data/*NaiveSequences.csv`: naive Gag/Gp41 sequences which were manually
  selected from over ??? published manuscripts.
- `data/patients.csv`: information of anonymous patients whoes sequence samples
  were used in this research.
- `data/samples.csv`: information of samples which were used in this research,
  included information of before and after treatments.
- `data/sequences.csv`: unaligned sequences which were used in this research.

The `result_data` folder stores processed data:

- `result_data/hyphy_output` folder: data which were created by HyPhy program:
  - `*.fel.output.txt`: the raw outputs from FEL analysis by HyPhy. See
    `scripts/run-fel.sh` for the parameters being used for calculation.
    The results were grouped by the gene and Rx.
  - `*.meds.output.txt`, `*.meds.csv`: the raw outputs from MEDS analysis. See
    `scripts/run-meds.sh` for the parameters.
  - `*.meds.result.csv`: the post-processed results by `scripts/medsproc.py`.
  - `*.nj.output.txt`: the neighbour joining outputs by HyPhy. See
    `scripts/run-nj.sh` for the parameters.
  - `*.pairwise.output.txt`: the outputs of pairwise dN/dS estimator. See
    `scripts/run-pairwise.sh` for the parameters.
  - `*.relax.output.txt`: the outputs of RELAX analysis. See
    `scripts/run-relax.sh` for the parameters.
- `*AAChangesByPosWPrev.csv`: aggregated synonymous and non-synonymous mutations
  by position and each amino acid with the prevalence data from
  `MutPrevalence.csv`. The results were grouped by gene.
- `*CodonChangesByPt.csv`: synonymous and non-synonymous mutations by position,
  each amino acid and each patient. The results were grouped by gene.
- `MutPrevalence.csv`: amino acid prevalence of each position. The prevalence
  data were calculated based on the manual selected naive sequences.
- `*.meds.tre`, `*.tre`: the Newick format tree files which were generated by
  neighbour joining process. The `*.meds.tre` were only used for MEDS analysis.
  
The `report` folder stores automatical report result:

- `README.md`: a markdown file contains readable tables and figures.
- `gag-mutations.png`, `gp41-mutations.png`: the graphical summary of gag/gp41
  sites at which amino acid mutations developed during therapy. There two
  figures were lately modified manually to highlight the MEDS, FEL and cleavage
  sites to be represented as the Figure 3 and 4 of the manuscript.

The `scripts folder stores script files:
- `analysis_functions.py`, `codonutils.py`, `data_reader.py`: internal Python
  libraries.
- `run-basic.sh`, `run-basic.py`: basic analysis. This generated the
  `MutPrevalence.csv`, `*AAChangesByPosWPrev.csv`, `*CodonChangesByPt.csv`.
- `run-nj.sh`: neighbour joining process. This called HyPhy to generate the
  newick trees for the aligned sequences from `data/fasta/` folder.
- `run-fel.sh`: run Fixed Effects Likelihood (FEL) method analysis.
- `run-meds.sh`: run Model of Episodic Directional Selection (MEDS) method
  analysis.
- `medsproc.py`: run post process for MEDS result. The script was from [HyPhy
  Wiki](http://www.hyphy.org/w/index.php/MEDS). The script was slightly changed
  to be worked with Python 3.5.
- `run-pairwise.sh`, `pairwise-estimator-dnds.bf`: run the pairwise dN/dS ratio
  estimator.
- `run-relax.sh`: run RELAX analysis.
- `run-final.sh`, `clean-fel.py`, `clean-pairwise.py`: clean the HyPhy outputs
  and store the cleaned result into `hyphy_cleaned_output` folder.
- `generate-report.sh`, `generate-report.py`: generate the report markdown file.
- `make-graphical-summary.r`: generate the two figures in `report` folder.


## Pipeline

We created an automatical pipeline to analysis the raw data, generate and store
the intermediate and final result. We believe that researchers can reuse this
pipeline easily by updating/replacing the raw files in the `data` folder.

The computer need to have [Docker](https://www.docker.com/) installed to use the
pipeline. All required softwares by this pipeline were carefully packed in a
[single Docker image](https://hub.docker.com/r/hivdb/gaggp41-runtime/). To
process the data, a high performance computer is highly recommended since HyPhy
is a compute-intensive program and it usually takes hours to finish the analysis
for ~100 sequences.

To start the automatical pipeline, just type `make` at the root folder of this
repository:

```sh
make
```

If everything goes well, you will see something like this:

```
/app/result_data/MutPrevalence.csv created
/app/result_data/GagAAChangesByPosWPrev.csv created
/app/result_data/Gp41AAChangesByPosWPrev.csv created
/app/result_data/GagCodonChangesByPt.csv created
...
```

### Pipeline steps

#### Step 1: basic

This step generated three types of basic files for further analysis:

- `result_data/MutPrevalence.csv`
- `result_data/*AAChangesByPosWPrev.csv`
- `result_data/*CodonChangesByPt.csv`

To run this step solely, type this command:

```sh
make basic
```

### Step 2: pairwise

This step calculated the pairwise dN/dS ratio with HyPhy. The output files:

- `result_data/hyphy_output/*.pairwise.output.txt`

We analyzed different regions of Gag and Gp41 sequences by using the different
parameter of codon range. For the Gag MA domain, the parameter was 0-395 (AA
position 1 - 132). For the Gag CA domain, the parameter was 396-1088 (AA
position 133 - 363). For the Gp41 CD domain, the parameter was 582-1034 (AA
position 195 - 345).

To run this step solely:

```sh
make pairwise
```

### Step 3: nj

This step calculated the newick tree from the aligned sequences with HyPhy.
The output files:

- `result_data/*.tre`
- `result_data/*.meds.tre`

The `*.tre` files were lately used by steps "relax" and "fel". The `*.meds.tre`
were lately used by "meds".

The parameters we used:

- branch lengths handling: \[Keep Negative\]
- distance formula: \[TN93\]

To run this step solely:

```sh
make nj
```

Noted the difference between our manuscript and the data in this repository was
mainly caused by this step. In the manuscript we used an old version of HyPhy
which has a slightly different approach of handling ambiguous DNA codes than
the newer version of HyPhy (see https://github.com/veg/hyphy/issues/618).

### Step 4a: relax

This step ran the RELAX method analysis with HyPhy. The output files:

- `result_data/hyphy_output/*.relax.output.txt`

This step usually takes an hour or more. To run this step solely:

```sh
make relax
```

### Step 4b: fel

This step ran the FEL method analysis with HyPhy. The output files:

- `result_data/hyphy_output/*.fel.output.txt`

This step usually takes an hour or more. To run this step solely:

```sh
make fel
```

### Step 4c: meds

This step ran the MEDS method analysis with HyPhy. The output files:

- `result_data/hyphy_output/*.meds.output.txt`
- `result_data/hyphy_output/*.meds.csv`
- `result_data/hyphy_output/*.meds.result.csv`

This step usually takes an hour or more. To run this step solely:

```sh
make meds
```

### Step 5: final

This step cleaned the output of step 4. The output files:

- `result_data/hyphy_cleaned_output/*`

To run this step solely:

```sh
make final
```

### Step 6: report

This step generated the automatical report. The output files:

- `report/*`

To run this step solely:

```sh
make report
```

### For Windows users

You can run the pipeline in Windows since there is a Windows version of Docker.
However, the `Makefile` is written for Linux/Mac computers and is not compatible
with Windows. If you are familar with `Makefile` and Docker I believe you can
figure this out easily without any help. If not, please [submit an
issue](https://github.com/hivdb/gag-gp41/issues/new) and I will come to to help
you.
