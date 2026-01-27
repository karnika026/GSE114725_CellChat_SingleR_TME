# GSE114725_CellChat_SingleR_TME
Single-cell RNA-seq analysis of breast cancer tumor microenvironment (GSE114725) using SingleR and CellChat.

# Cell–Cell Communication Analysis in Breast Cancer (GSE114725)

This repository contains a reproducible workflow for analyzing **cell–cell communication** in the breast cancer tumor microenvironment (TME) using **3′ single-cell RNA-seq data (GSE114725)**.  
The analysis integrates **automated cell type annotation using SingleR** with **ligand–receptor interaction inference using CellChat**, enabling unbiased characterization of intercellular signaling in tumor and normal tissues.

## Dataset Information
- **Accession:** GSE114725  
- **Data type:** 3′ single-cell RNA-seq  
- **Biological context:** Tumor and normal breast tissue samples  
- **Source:** NCBI Gene Expression Omnibus (GEO)  

Raw data are available from GEO. No raw data files are stored in this repository.

## Research Objective
> To identify and characterize cell–cell communication networks within the breast cancer tumor microenvironment and to compare signaling pathways between tumor and normal tissues.

## Tools & Packages
- **R (≥ 4.5)**
- **Seurat** – preprocessing, clustering, and UMAP visualization  
- **SingleCellExperiment** – data structure interoperability  
- **SingleR** – automated, reference-based cell type annotation  
- **celldex** – reference datasets for SingleR  
- **CellChat** – inference of ligand–receptor–mediated cell–cell communication  
- **ggplot2, patchwork** – visualization  


## Analysis Workflow

1. **Input Data**
   - A preprocessed Seurat object generated in a previous project
   - Data normalized, clustered, and embedded using UMAP

2. **Conversion to SingleCellExperiment**
   - Seurat object converted to SingleCellExperiment for compatibility with SingleR and CellChat

3. **Automated Cell Type Annotation (SingleR)**
   - Reference dataset: Human Primary Cell Atlas (celldex)
   - Cell identities assigned in an automated and unbiased manner
   - Generates a `SingleR_label` column in the metadata
   - No manual marker-based annotation performed

4. **CellChat Object Construction**
   - Cells grouped using SingleR-derived cell type labels
   - Human ligand–receptor interaction database applied

5. **Cell–Cell Communication Analysis**
   - Inference of communication probabilities and signaling pathways
   - Comparison of interaction networks between tumor and normal samples
   - Visualization of global and pathway-specific communication patterns

## Repository Structure
GSE114725_CellChat_TME/
├── scripts/
│ ├── 01_Seurat_to_SCE.R # Convert Seurat object to SingleCellExperiment
│ ├── 02_SingleR_Annotation.R # Automated cell type annotation using SingleR
│ ├── 03_Create_CellChat_Object.R # CellChat object construction
│ ├── 04_CellChat_Analysis.R # Communication probability and pathway analysis
│ └── 05_Visualization.R # Network and UMAP visualizations
├── figures/ # Output figures
└── README.md

## Key Outputs
- Automated cell type annotations (`SingleR_label`)
- CellChat object built using SingleR-derived cell identities
- Ligand–receptor interaction networks across major cell populations
- Visualizations including:
  - UMAP colored by SingleR annotations
  - Cell–cell communication network plots
  - Signaling pathway heatmaps

## Limitations
- Automated annotation may misclassify rare cell populations or tumor subtypes  
- CellChat focuses on coarse cell type interactions rather than fine-grained subclustering  
- Ligand–receptor interactions are computational predictions and require experimental validation  

## Future Directions
- Subtype-specific analysis of tumor epithelial cells  
- Integration with spatial transcriptomics or other multi-omics data  
- Validation of key signaling pathways in independent breast cancer datasets  
- Dynamic communication analysis using trajectory or pseudotime approaches  

---

## Data Availability
Raw data are available from the Gene Expression Omnibus under accession **GSE114725**.
