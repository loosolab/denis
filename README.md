# DENIS - DE Novo motIf diScovery
<img width="12%" src="/figures/Denis_V3.png">

# Introduction
With ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) the genome-wide chromatin accessibility can be detected. This is done by a hyperactive Tn5 transposase, which cuts fragments only in open regions of the DNA. After amplification, sequencing, and following analysis of the ATAC-seq data, open chromatin regions are identified as an accumulation of reads. Further focus on open chromatin regions reveals footprints as small spaces of less read coverage, where transcription factors were bound. Applied to a motif database, footprints can be linked to a motif, and assumptions about transcription factor roles, in the regulatory network can be inferred. But not all of the footprints are traceable to motifs, implying the existence of unknown motifs. To unfold function and binding of de novo motifs, we present a pipeline that integrates with the [TOBIAS framework](https://github.com/loosolab/TOBIAS/) and enhances it for motif generation.


# Overview
<img align="right" width="50%" src="/figures/figure1.png">
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<ol type="a">
    <li>The framework receives continuous binding scores as minimal input. In the first step, the binding scores are scanned and FPs are extracted in a genome wide scale. If binding scores for two conditions are provided as input, the framework calculates differential binding scores.</li>
    <li>Sequences at detected FPs are extracted and used to perform data driven de novo motif prediction. In an iterative process, most prominent motifs are extracted and corresponding sequences are removed from the motif discovery pool. FPs not (yet) used are kept for the following discovery run.</li>
    <li>Resulting motifs are subjected to downstream analysis. DENIS provides i) an individual motif site annotation module (e.g. annotation to genes), ii) a motif database comparison module, iii) a module to compute global binding site enrichment in open chromatin, and iv) based on the annotation module a motif based gene set enrichment analysis (GSEA) module.</li>
</ol>

<br clear="right"/>

# How to run
To run this pipeline make sure Snakemake and Conda are installed and working on your machine.
The pipeline can be installed by cloning the repository.
```
git clone https://github.com/loosolab/denis.git
```
Next the configuration file `config.yml` needs to be adjusted for the data. After this is done the pipeline can be started with the following command:
```
snakemake --configfile [path_to/]config.yml --cores [number of cores] --use-conda --conda-frontend conda
```
If [Mamba](https://github.com/mamba-org/mamba) is installed the building time of environments can be greatly reduced running this command instead:
```
snakemake --configfile [path_to/]config.yml --cores [number of cores] --use-conda --conda-frontend mamba
```

# Example run
Data in the `example` folder can be used to do a small example run. Note that GO-enrichment analysis will be skipped unless an email is added in `example/example_config.yml`. Start the example run with:
```
snakemake --configfile example/example_config.yml --use-conda --cores [number of cores]
```
Pipeline output will be added to `example_output` folder located in the main directory.

# Output
After the pipeline is done, an output folder in the following format will be produced:

- <path_to_output>/
  - **1_footprints/**
    - 1_extraction/
    - 2_fasta/
  - **2_discovery/**
    - 1_meme/
      - 1_motifs/
      - iteration_0/, iteration_1/, ...
    - 2_processed_motifs/
      - motifs/
      - dendrogram.png, distance.tsv
    - 3_control_motifs/
  - **3_evaluation/**
    - motif_evaluation/
    - rescan/
    - venn/
    - motif_ranks.tsv, motifs.meme, ranks_barplot.pdf
  - **4_annotation/**
    - open_binding_sites/
    - uropa/
    - feature_enrichment_plot.pdf, feature_enrichment_table.tsv
  - **5_gene_ontology/**
  - logs/

## 1_footprints
This folder contains the results of the footprint extraction step. It holds statistics diagrams as well as bed files with all footprints, footprints assigned to motifs, and footprints after motif filtering (`1_extraction`). The footprints without motifs are also available in fasta format, which is used in the next step (`2_fasta`).

## 2_discovery
This folder is divided into three parts. `1_meme` contains the results of each motif discovery iteration, as well as a summary directory with all motifs and their binding sites. `2_processed_motifs` takes the motifs from the previous step and builds a consensus for similar motifs. The processed motifs are stored in `motifs/` with a dendrogram and distances next to the folder. The third directory `3_control_motifs/` contains motifs selected from the database (if available), which are used to provide the context in the following steps.

## 3_evaluation
Anything related to the evaluation of motifs is provided in this folder. The `motif_evaluation` folder contains plots showing the distance to known motifs in multiple ways. The `rescan`, `venn` folder, and `ranks_barplot.pdf` give insights on genome-wide binding, binding in open chromatin, and enrichment between the two, which is also accompanied by a table.

## 4_annotation
Each motif is annotated through a [UROPA](https://uropa-manual.readthedocs.io/) run. To do this the binding sites for each motif within open chromatin are collected (`open_binding_sites`). The binding sites are then annotated using the provided UROPA template and the results for each run are stored in `uropa`. All annotations are then summarized into a feature enrichment plot and a table.

## 5_gene_ontology
A gene ontology enrichment analysis is conducted based on the annotated genes of each motif found in the previous step. The analysis creates multiple trees of significant ontologies (celluar components, biological process & molecular function) as well as a table containing all information needed to create said trees.

# Contact
Hendrik Schultheis (hendrik.schultheis@mpi-bn.mpg.de)
