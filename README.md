# Motif discovery pipeline

## Introduction
With ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) the genome-wide chromatin accessibility can be detected. This is done by a hyperactive Tn5 transposase, which cuts fragments only in open regions of the DNA. After amplification, sequencing, and following analysis of the ATAC-seq data, open chromatin regions are identified as an accumulation of reads. Further focus on open chromatin regions reveals footprints as small spaces of less read coverage, where transcription factors were bound. Applied to a motif database, footprints can be linked to a motif, and assumptions about transcription factor roles, in the regulatory network can be inferred. But not all of the footprints are traceable to motifs, implying the existence of unknown motifs. To unfold function and binding of de novo motifs, we present a pipeline that integrates with the [TOBIAS framework](https://github.com/loosolab/TOBIAS/) and enhances it for motif generation.

## How to run
To run this pipeline make sure Snakemake and Conda are installed and working on your machine.
The pipeline can be installed by cloning the repository.
```
git clone https://gitlab.gwdg.de/loosolab/software/motif-discovery-pipeline.git
```
Next the configuration file `config.yml` needs to be adjusted for the data. After this is done the pipeline can be started with the following command:
```
snakemake --configfile [path_to/]config.yml --cores [number of cores] --use-conda 
```
If Mamba is installed the building time of environments can be greatly reduced running this command instead:
```
snakemake --configfile [path_to/]config.yml --cores [number of cores] --use-conda --conda-frontend mamba
```

## Contact
Hendrik Schultheis (hendrik.schultheis@mpi-bn.mpg.de)
