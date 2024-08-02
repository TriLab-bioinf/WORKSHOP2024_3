# Guide 1: Basic differential gene expression analysis
## Comparison between two groups (e.g. Case vs Control ; Female vs Male ; Mutant vs WT)

For any differential expression analysis we will need to prepare and load into R the following two files:
a. A **read counts per gene** table, containing `genes` as rows and `samples` as columns:

![image](https://github.com/user-attachments/assets/94d0390b-a2cc-4ecc-8b3d-cd47a6f276d3)

b. A **metadata table** containing `samples` as rows and `variables` as columns:

![image](https://github.com/user-attachments/assets/887a2415-fcdb-43ae-a902-d68c6fb66b1d)

**Note:** `sample-columns` in the `read_counts` table have to be in the same order as `sample-rows` in the `metadata` table. The order can be fixed using R.

Let's start the analysis by opening a new **R notebook** in Rstudio:

In Rstudio go to `File > New File > R Notebook`. A new notebook template will open. Delete all lines from line 6 to the end and save it as `guide_1.Rmd`.

Go to `Code > Insert Chunck`

Within the new area that appeared, **named chunk**, you can enter your R code. Above or below the chunks you can add regular text to comment your code or conclussions.

![](https://github.com/TriLab-bioinf/WORKSHOP2024_3/blob/main/Figures/chuhk.png)

Now let's start processing the data.

# A. Load required R packages

```{r}
pacman::p_load(BiocManager, DESeq2, readxl, tidyverse, pheatmap, ashr)
```

# B. Load read counts table

```{r}
# Load raw read counts
counts.raw <- read.delim(file = "data/counts_M4.txt", 
                         header = T, sep = "\t",row.names = 1)

counts.raw
```

# C. Load metadata

```{r}

# Load metadata and clean sample_ids
metadata <- read.delim(file = "data/metadata_M4.txt", 
                         header = T, sep = "\t",row.names = 1)

# Check metadata row names == counts column names
same_order <- all(rownames(metadata) == colnames(counts.raw))

# Sort counts columns to match order in metadata rows
if(!same_order){
  print("Sorting columns in counts.all")
  counts.raw <- select(counts.raw, rownames(metadata))
}

# Include total read counts in metadata
metadata$read_counts <- colSums(counts.raw, na.rm = TRUE)

metadata

# Convert categorical variables to factors
metadata$Treatment <- factor(metadata$Treatment)
```

# D. Clean gene expression data

## 1. Filter genes based on read counts per sample

```{r}
# Print out dimension of the counts.raw dataframe 
dim(counts.raw)

# For each gene, set the minimum number of reads per sample and the minimum group size
min_num_reads_per_sample <- 10
min_group_size <- 6

# Filter genes
keep <- rowCounts(counts.raw > min_num_reads_per_sample) > min_group_size
counts.fil <- counts.raw[keep,]

# Print out dimension of the counts.fil dataframe
dim(counts.fil)
```

## 2. Optional: Filter genes based on their standard deviation across samples

```{r}
# Calculate Standard Deviations for each gene across all samples
gene.sd <- rowSds(as.matrix(counts.fil))

# Plot histogram of Std Dev:
hist(gene.sd, breaks = 100000, xlim = c(0,50))

summary(gene.sd > 10)

# Filter genes
keep <- gene.sd > 10
counts.fil <- counts.fil[keep,]

# Print out new dimension of counts.fil dataframe
dim(counts.fil)

```

# E. Generate DESeq object

```{r message=FALSE, warning=FALSE}
dds <- DESeqDataSetFromMatrix(countData = counts.fil,
                              colData = metadata,
                              design = ~1) 

dds <- DESeq(dds)
```

# F. Exploratory analysis

## 1. Principal component analysis

```{r}
# 1. variance stabilizing transformation of counts (optinally you can use log transformation with rlog)
dds.vst <- vst(dds, blind=TRUE)

# 2. Generate PCA plot
pca.p <- plotPCA(dds.vst, 
        intgroup = c("Treatment"),
        ntop = 500,
        returnData = FALSE,
        pcsToUse = 1:2) 

# 3. Save PCA plot
ggsave(filename = "pca.pdf", plot = pca.p)

# 4. Print out PCA plot
pca.p

```

## 2. Heatmap

```{r}
# 1. Normalize counts using Variance Stabilizing Transformation
dds.vst <- vst(dds, blind=TRUE)

# 3. Calculate distances between samples
sampleDists <- dist(t(assay(dds.vst)))

# 4. Calculate inter-sample distances
sampleDist.mat <- as.matrix(sampleDists)

# 5. Add Treatment levels to sample names
rownames(sampleDist.mat) <- paste(rownames(sampleDist.mat), dds.vst$Treatment)

# 6. Generate heatmap plot
hm.all.p <- pheatmap(mat = sampleDist.mat,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists)

# 7. Save heatmap plot
ggsave(filename = "heatmap.pdf", plot = hm.all.p)
```

# G. Run differential expression

We will contrast *`Drug`* vs *`Control`* factors of the *`Treatment`* variable.

```{r message=FALSE, warning=FALSE}
# 1. Add design formula
design(dds) <- ~Treatment

# 2. Recompute dispersions
dds <- DESeq(dds)

# 3. Print out coeficients
resultsNames(dds)

# 4. Get DE results ------------------------------------------------

# 4.1 Method 1: Factor pairwise comparison 
res.1 <- results(dds, contrast = c("Treatment", "Drug", "Control"))

# 4.2 Method 2: Using coeficients
res.2 <- results(dds, name=c("Treatment_Drug_vs_Control"))
res.3 <- results(dds, contrast=list("Treatment_Drug_vs_Control"))


```

# H. Log2FC Shrinkage

Distribution of Log2FC values usually is directly correlated with gene expression levels. Therefore, it is advisable to correct Log2FC using the DESeq2 function `lfcShrink()`.

```{r message=FALSE, warning=FALSE}
# 1. Plot Log2FC as a function of mean normalized gene expression
plotMA(res.2, ylim=c(-3,3))

# 2. Calculate DE results, adjusting Log2FC based on gene expression
res.shk <- lfcShrink(dds, contrast = c("Treatment", "Drug", "Control"), type = "ashr")

# 3. Plot adjusted Log2FC as a function of mean normalized gene expression
plotMA(res.shk, ylim=c(-3,3))
```

# I. Clean and save results table

```{r}
# 1. Sort by sdj.p
res.shk.sorted <- res.shk[order(res.shk$padj), ]

# 2. Select genes with  -1 <= Log2FC >= 1
res.shk.sorted.filt <- subset(res.shk.sorted, abs(log2FoldChange)>=1)

# 3. Save results to a file
write.table(x = as.data.frame(res.shk.sorted.filt), file = "DE_res_shk_filt.txt", sep = "\t", col.names = NA)
```
