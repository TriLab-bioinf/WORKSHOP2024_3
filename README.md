## WORKSHOP2024_3

## Differential gene expression analysis with R

### Workshop outline

**A) Basic differential gene expression analysis:**

1. How to load the input data into R
2. Clean gene expression data
3. Read counts normalization (TPM, FPKM/RPKM, CPM) -> use featureCounts?
4. Generate DESeq object
5. Basic experimental design formula (1 categorical variable (e.g. Treatment) with 2 levels (e.g. Case, Control)).
6. Vst/rlog normalization for exploratory analysis
7. Exploratory analysis: PCA, correlation or dist analyses.
8. Run DE analysis with `results()` function
9. Interpret meaning of DE-table columns
10. MA-plot to see dependency of Log2FC on gene expression
11. Shrinkage
12. Save results table into a file

**B) Advanced DE analysis**

1. Experimental design with 1 categorical variable with three levels
2. Paired experimental design
3. Paired experimental design with 2 categorical variables
4. Longitudinal design (LRT) with one categorical variable

**C) Convert gene accession numbers to gene symbols**

1. Function biomart
2. 


