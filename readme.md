Characterization of the multi-omic signature of Minimal Hepatic
Encephalopathy in the hippocampus of the *Rattus norvegicus* model
induced by diet.
================
Patricia Checa Gómez, Ana Clavijo Pizzano, José Luis del Río Vázquez,
Rocío Cañete Gómez

# Multi-Omic Integration Study in Hyperammonemic Rats (HA Rats 1 Dataset)

This repository contains the integrative multi-omic analysis performed
on a rat model of minimal hepatic encephalopathy (MHE), specifically
under a hyperammonemic diet (HA rats 1 dataset). The main goal is to
characterize molecular changes in the hippocampus through the
integration of transcriptomics, regulomics, and metabolomics.

------------------------------------------------------------------------

## Project Overview

### 1. Data Description, Objectives, and Biological Question

- Multi-omic data from *Rattus norvegicus* hippocampus of rats under
  hyperammonemic diet and controls.
- Omics types: RNA-Seq, miRNA-Seq, and metabolomics.
- General objective: To characterize integrated changes in
  transcriptome, regulome, and metabolome to discriminate between
  control and disease groups.
- Key question: Is there a characteristic multi-omic signature in MHE
  that can distinguish the groups?

### 2. Integration Methods and Tools Used

- Statistical integration with R, using **MORE** and **MixOmics**
  packages (MultiBlock sPLS-DA - DIABLO).
- Conceptual analysis with the web tool **PaintOmics4**.
- Statistical power calculation with **MultiPower**.
- Regulatory network visualization using **Cytoscape**.

### 3. Initial data exploration

- Boxplots for distribution assessment and transformations.
- RNA-Seq and miRNA-Seq data centered.
- Metabolomics data centered and scaled.
- Good normalization and data consistency observed.

### 4. Sample size and statistical power calculation

- Desired power: 60% per omic, 80% overall.
- Sample size estimated with MultiPower, controlling FDR at 5%.
- Assumed unit experimental cost due to lack of specific cost data.

### 5. MultiBlock sPLS-DA model for integration and discrimination

- DIABLO model to identify key variables discriminating between groups.
- Validation with Leave-One-Out technique due to small sample size.
- RNA-Seq and miRNA-Seq showed better discriminative ability than
  metabolomics.

### 6. Regulatory analysis with MORE

- Identification of significant regulators (miRNAs) via multiple linear
  regression with elastic net.
- Selection of relevant regulators modulating multiple genes.
- Differential regulatory network visualization in Cytoscape.

### 7. Functional and conceptual integration with PaintOmics

- Integration of RNA, miRNA, and metabolites to identify metabolic and
  signaling pathways.

### 8. Discussion on alternative methods
