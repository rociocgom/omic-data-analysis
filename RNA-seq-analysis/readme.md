RNA-seq Analysis of LEC1: Transcriptome Profiling and Differential
Expression
================
Rocío Cañete Gómez, Diana Andreea Baicea, Alberto Barrero Gonzalez,
Modesto Berraquero Perez,Pablo Cano Jimenez

# Overview

This project analyzes gene expression changes in *Arabidopsis thaliana*
following induction of the transcription factor HB21. Data from GEO
accession **GSE249766** were processed to identify differentially
expressed genes associated with floral arrest.

# Background

HB21 is a key transcription factor regulating floral arrest at the end
of flowering. Using an inducible system (HB21ind), HB21 expression can
be activated by dexamethasone (Dex), enabling study of its specific
effects without interference from related genes HB40 and HB53.

# Experimental Design

- Organism: *Arabidopsis thaliana* (genome version TAIR10.1)  
- Treatments:
  - Mock control (0.015% Silwet L-77)  
  - Dex treatment (10 μM dexamethasone + 0.015% Silwet L-77)  
- Sample collection: Inflorescence apices collected 6-8 hours
  post-treatment, excluding open flowers  
- Replicates: 3 biological replicates per condition  
- Sequencing: Illumina NovaSeq 6000, paired-end reads

# Workflow

1.  Download raw reads, reference genome, and annotations.  
2.  Generate STAR genome index and align reads.  
3.  Perform quality control with FASTQC.  
4.  Assemble transcripts and quantify expression with StringTie.  
5.  Create count matrix using prepDE.py.  
6.  Normalize data with NormalyzerDE and perform PCA with FactoMineR.  
7.  Identify differentially expressed genes using DESeq2.  
8.  Conduct functional enrichment analyses (GO, KEGG) with
    clusterProfiler.  
9.  Visualize results with plots such as volcano plots.
