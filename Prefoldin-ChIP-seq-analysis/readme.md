Constructing an Integrated Model of the Transcriptional Role of Human
Prefoldin in Chromatin
================
Rocío Cañete Gómez Master’s Thesis – Omics Data Analysis and Systems
Biology, University of Seville & International University of Andalusia

# Overview

This repository contains four R scripts used for the analysis of
**ChIP-seq** data for three proteins:

- **RPB1**
- **SPT16**
- **macroH2A1**

The study focuses on the effect of **PFDN5 subunit knockout** on the
transcriptional removal of the histone variant **macroH2A1** by the
**FACT complex**.  
Input samples were included as controls.

------------------------------------------------------------------------

## Scripts

1.  **alignment**  
    Aligns raw FASTQ reads to the human reference genome and the
    *Saccharomyces cerevisiae* spike-in normalization.

2.  **reads count**  
    Counts aligned reads over genomic regions (genes).

3.  **differential analysis**  
    Performs differential binding analysis between KO and WT conditions
    using the **limma** tool.

4.  **overlapping and hypergeometric analysis**  
    Computes overlaps between differentially occupied gene lists and
    performs a hypergeometric analysis to assess the significance of
    overlaps.

------------------------------------------------------------------------

## Notes

- **R version:** ≥ 4.4.2  
- **Tested on:** Windows 10  
- **Language:** R  
- Dependencies and parameters are detailed within each script.ç
