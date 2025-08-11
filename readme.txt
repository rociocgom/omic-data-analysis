Master's thesis: "Constructing an integrated model of the transcriptional role of human prefoldin in chromatin: the FACT histone chaperone and macroH2A1" - Omics Data Analysis and Systems Biology, University of Seville and International University of Andalusia. 

This folder contains four scripts used in the analysis of ChIP-seq data for three proteins (RPB1, SPT16 and macroH2A1) for my Master's Thesis. The study focuses on the effect of PFDN5 subunit knockout on transcriptional removal of the histone variant macroH2A1 by the FACT complex. Input samples were included as controls.

Scripts included:

1. alignment: Aligns raw FASTQ reads to the human reference genome and the "Saccharomyces cerevisiae" spike-in normalization.

2. reads count: Counts aligned reads over genomic regions (genes).

3. differential analysis: Performs differential binding analysis between KO and WT conditions using the limma tool. 

4. overlapping and hypergeometric analysis: Computes overlaps between differential occupied gene lists and performs hypergeometric analysis to assess significance to overlaps. 

Notes:
- R version: ≥ 4.4.2, tested on Windows 10. 
- Scripts are written in R.
- Dependencies and parameters are detailed within each script.