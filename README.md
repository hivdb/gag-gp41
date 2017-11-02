# Selection Analyses of Paired HIV-1 gag and gp41 Sequences Obtained Before and After Antiretroviral Therapy

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

### Step 1: basic

This step generated three types of basic files for further analysis:

- `internalFiles/aaChangesByPosWPrev/*.csv`
- `internalFiles/codonChangesByPt/*.csv`

To run this step solely, type this command:

```sh
make basic
```

### Step 2: pairwise

This step calculated the pairwise dN/dS ratio with HyPhy. The output files:

- `internalFiles/hyphyOutput/*.pairwise.output.txt`

To run this step solely:

```sh
make pairwise
```

### Step 3: nj

This step calculated the newick tree from the aligned sequences with HyPhy.
The output files:

- `data/phylo/*.tre`
- `internalFiles/phylo/*.tre`

The tree files in `internalFiles` were lately used by step "meds".

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

### Step 4a: fel

This step ran the FEL method analysis with HyPhy. The output files:

- `internalFiles/hyphyOutput/*.fel.output.txt`

This step usually takes an hour or more. To run this step solely:

```sh
make fel
```

### Step 4b: meds

This step ran the MEDS method analysis with HyPhy. The output files:

- `internalFiles/hyphyOutput/*.meds.output.txt`
- `internalFiles/hyphyOutput/*.meds.csv`
- `internalFiles/hyphyOutput/*.meds.result.csv`

This step usually takes an hour or more. To run this step solely:

```sh
make meds
```

### Step 5: final

This step cleaned the output of step 4. The output files:

- `internalFiles/hyphyCleanedOutput/*`

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

## Pipeline for Processing PI-naïve Studies

We reviewed and identified studies in which all of the individuals in the study
were PI-naïve.

### Step 1: (manually) download sequences from LANL

We first access LANL HIV Sequence Database for retriving full-size gag and gp41
sequences. Following described the way to download these two files:

1.  Open https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html.
2.  In "define start and end", type 790-2289 for gag, 7758-8792 for gp41 (both
    without the ending stop codons).
3.  Check "include problematic sequences", we will reevaluate every sequence
    later by our scripts and remove disqualified ones.
4.  Click "search", and wait until the list loaded.
5.  Click "download sequences".
6.  In "gap handling", select "sequeeze".
7.  Check "include HXB2 Reference Sequence (K03455)".
8.  Click "compose a label".
9.  Delete the default numbers 1, 2, 3, 4 from the input boxes first. Then add
    numbers like this:
    - Accession: 1
    - PAT id(SSAM): 2
    - Subtype: 3
10. DO NOT delete/change "field separator" or "missing charater".
11. Click OK.
12. Input your contact information.
13. Check your email inbox in 1-3 hours for FASTA-format files.
14. rename the FASTA files to `hiv-db_gag_squeeze.fasta` (for gag) and
    `hiv-db_gp41_squeeze.fasta` (for gp41). Put both files under the
    `gag-gp41/local` folder.

### Step 2: Run the main analysis

This step runs an automatic script to clean and filter the sequences retrieved
by the first step. Then it exports the FASTA files, AA prevalence files and
other statistics files.

Every file in `data/naiveStudies` were created by this step.

```sh
make naive
```

### Step 3 (optional): Diversify analysis

This step calculates the pairwise distances of all naive sequences found by
step 2. Then it generates the two `reports/*-naive-diversify.png` figures.

We take this step out of step 2 since the distance calculation is very CPU
intensive. It may take hours if you don't have a powerful computer.

```sh
make naive-dist
```

## For Windows users

You can run the pipeline in Windows since there is a Windows version of Docker.
However, the `Makefile` is written for Linux/Mac computers and is not compatible
with Windows. If you are familar with `Makefile` and Docker I believe you can
figure this out easily without any help. If not, please [submit an
issue](https://github.com/hivdb/gag-gp41/issues/new) and I will come to to help
you.
