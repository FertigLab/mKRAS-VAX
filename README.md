# J1994.scRNA-TCR
This repository contains processed data and analysis scripts used in the manuscript by Huff et. al, "An off-the-shelf vaccine activates a broad repertoire of mutant KRAS-specific T cells in patients with resected pancreatic cancer". Code written by Alexander Girgis and Ludmila Danilova. 

## Directory Structure 
R code files are located in parent directory. These files should be executed in their current location relative to all subfolders in directory. 

/code_outputs/: Empty folder which will become populated with resluts files upon executing R code files. 
/_rfunctions/: Code folder contaning miscellaneous helper functions in processing TRBV nomenclature. 
/raw_data/: Contains scRNA and scVDJ output from CellRanger pipeline. Currently not available due to data sharing permissions, but will be uploaded to dbGaP and appropriate repositories. 
/code_inputs/: Processed data objects loaded and used in at least one of the R code files. Currently not available due to data sharing permissions, but will be uploaded to dbGaP and appropriate repositories. 
## Scope of Analysis Scripts
The contents of this repo include: 

I. QC and annotation of single-cell RNA data for 5 subjects receiving mKRAS-Vax. Isolation of annotated T cells from single-cell data. 
		A.2023_05_29_Read10x.r: Process CellRanger output into aggregate Seurat object.
		B.2023_06_06_sc-by-batch.r: Split aggregate object by sequencing batch. Perform separate QC and filtering. 
		C.2023_11_08_sc-cluster-assignment-batch1.r: Annotate putative T cell clusters based on cluster markers and Azimuth annotation for sequencing batch 1. 
		D.2023_11_08_sc-cluster-assignment-batch2.r: Annotate putative T cell clusters based on cluster markers and Azimuth annotation for sequencing batch 2. 
II. Analysis of T-cell expansion assays to identify clonotypes expanded to mKRAS peptides. 
		E.2023_05_05_DiffExFESTData.r: Analyze expansion assay data to identify T cell receptor clonotypes expanded to mKRAS peptides. 
		F.2024_04_18_MakeFilteredDiffExTable.r: Generate compiled table of all mKRAS expanded TCRs 
III. Integration of single-cell and TCR expansion data to identify anti-mKRAS T cells in single-cell data. 
		G.2024_04_17_scKRASAnnot.r: Integrate mKRAS expansion data and single-cell data to create a new Seurat object with annotations for anti-mKRAS T cell reactivity. 
		H.2024_05_08_pullAllmKRAS_withExpansionInfo-all.r: Generate composite table containing paired chain information, T cell phenotype, and mKRAS reactivity for anti-mKRAS TCRs in scRNA data.
		I.2024_04_27_MinorEditingObjects.r: Combine annotated T cell objects for batch 1 and batch2 back into one seurat object for convenience. 
IV. Analysis of gene expression profile of mKRAS-specific T-cells. 
		J.2024-04-23_Seurat_viz.Rmd
V. Comparison of mKRAS TCRs with TCGA data. 
		K.024-04-29_TCGA.Rmd
## Other Files
		LICENSE: GNU v3.0 license governing the use and sharing of code in this repository. 
		SESSIONINFO: R Package versions for libraries used in analysis scripts for I-III described above. 
