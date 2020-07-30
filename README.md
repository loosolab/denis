# Motif discovery pipeline

## Introduction
With ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) the genome-wide chromatin accessibility can be detected. This is done by a hyperactive Tn5 transposase, which cuts fragments only in open regions of the DNA. After amplification, sequencing, and following analysis of the ATAC-seq data, open chromatin regions are identified as an accumulation of reads. Further focus on open chromatin regions reveals footprints as small spaces of less read coverage, where transcription factors were bound. Applied to a motif database, footprints can be linked to a motif, and assumptions about transcription factor roles, in the regulatory network can be inferred. But not all of the footprints are traceable to motifs, implying the existence of unknown motifs. To unfold function and binding of de novo motifs, we present a pipeline that integrates with the [TOBIAS framework](https://github.com/loosolab/TOBIAS/) and enhances it for motif generation.

## Snakemake setup
TODO

## Contact
Hendrik Schultheis (hendrik.schultheis@mpi-bn.mpg.de)
