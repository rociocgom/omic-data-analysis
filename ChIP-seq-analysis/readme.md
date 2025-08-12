LEC1 ChIP-seq Analysis in Arabidopsis thaliana
================
Rocío Cañete Gómez, Diana Andreea Baicea, Alberto Barrero Gonzalez,
Modesto Berraquero Perez,Pablo Cano Jimenez

# Overview

This project analyzes ChIP-seq data to identify genome-wide binding
sites of the transcription factor **LEAFY COTYLEDON1 (LEC1)** in
*Arabidopsis thaliana* seeds. The goal is to elucidate LEC1 target genes
and regulatory networks involved in seed development and stress
responses.

# Background

LEC1 is a key regulator of seed maturation controlling storage
accumulation and desiccation tolerance. Beyond maturation, LEC1
influences hormonal and light-response genes, but its genome-wide
targets and precise regulatory roles remain incompletely understood.
This study compares wild-type and lec1 mutant seeds, integrating RNA-seq
and ChIP-seq data to characterize LEC1-dependent gene regulation.

# Experimental Design

- Organism: *Arabidopsis thaliana*  
- GEO accession: **GSE99587**  
- Samples:
  - 2 biological replicates of lec1 mutants expressing LEC1-GFP (ChIP
    samples)  
  - 2 biological replicates input controls with nonspecific antibodies  
- Technique: Chromatin immunoprecipitation followed by sequencing
  (ChIP-seq) using anti-GFP antibody  
- Objective: Identify LEC1 binding sites genome-wide in developing seeds

# Workflow

1.  Download raw fastq files for ChIP and input samples.  
2.  Build Bowtie2 index and align reads to the reference genome.  
3.  Process alignments with samtools for sorting and indexing.  
4.  Generate normalized coverage (.bw) files with bamCoverage.  
5.  Perform peak calling with macs3 using input as control.  
6.  Intersect peaks across replicates to obtain consensus binding
    sites.  
7.  Annotate peaks using ChIPseeker and analyze genomic distribution
    relative to TSS.  
8.  Identify LEC1 target genes and perform GO functional enrichment
    analysis with clusterProfiler.  
9.  Visualize motif enrichment and peak distribution using R and IGV.

