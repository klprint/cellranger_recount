# Rationale
- recount cellranger output bam files using different gtf files
- in cellranger this feature is missing and the alignment needs to be redone when recount is wanted
- should only need this folder, the bam file and a gtf file
  - obviously this gtf file needs to be compatible to the aligned genome (version)

# Installation
`cellranger-recount` was programmed using the conda framework, so I recommend to install `conda` on your computer (Mac or Linux) and create a new conda-environment using the `environment.yml` file:

```bash
conda env create --file=environment.yml
```

This will create a conda-environment called "cellranger-recount".
Activate the environment using:

```bash
source activate cellranger_recount
```



# How to use
To get a help message for usage, enter:

```bash
python ./cellranger_recount.py --help
```

`cellranger_recount` relies on a BAM file, produced by cellranger, and a aligned genome compatible GTF file. 
By default, `cellranger_recount`, counts all genewise features flagged as `transcript` in the GTF file.
This means, that if multiple transcripts are known, all reads/UMI of these transcripts will be aggregated in the associated gene (GTF file's `geneID` metadata).

A run with 140e6 alginments in the BAM file took from start to end seven hours, 16GB of ram and a single core (2.7Ghz).

Read the **Installation** section to activate the conda-environment, to be sure that you use the right dependencies on your machine.

To count UMIs for all transcripts of all annotated genes in the genome use:

```bash
python cellranger_recount.py \
    path/to/bam/file.bam \
    path/to/GTF.gtf.gz \
    path/to/output/folder
```

To count UMIs for exons only:

```bash
python cellranger_recount.py \
    -f exon \
    path/to/bam/file.bam \
    path/to/GTF.gtf.gz \
    path/to/output/folder
```

The output will be, as cellranger does by default, in the Matrix Market format.
All tools which are able to read cellrangers output files (`genes.tsv`, `barcodes.tsv` and `matrix.mtx`) should be able to open `cellranger_recount`'s output.

# Todo
## TODO handle multiple occuring IDs
Here, it is important to provide an easy way to allow users to specify which ID should be used for counting.
There might be a problem, if features overlap (i.e. exons) so each UMI might be counted multiple times.
Maybe using the geneID in any case might be the best solution and only accept unique UMI - cell barcodes pairs.
