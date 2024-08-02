# Guide 1: Differential gene expression analysis with R

# A. Load libraries

Load required R packages

```{r Step_1, message=FALSE, warning=FALSE}

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
metadata$Sbj_id <- factor(metadata$Sbj_id)
```

# D. Clean gene expression data

## 1. Filter genes based on read counts per sample

```{r}

dim(counts.raw)
keep <- rowCounts(counts.raw > 10) > 4
counts.fil <- counts.raw[keep,]
dim(counts.fil)
```

## 2. Optional: Filter genes based on their standard deviation across samples

```{r}

hist(rowSds(as.matrix(counts.fil)), breaks = 100000, xlim = c(0,50))

summary(rowSds(as.matrix(counts.fil)) > 10)

keep <- rowSds(as.matrix(counts.fil)) > 10

counts.fil <- counts.fil[keep,]
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

# Method 1: ---------------------------------
# easy

DESeq2::plotPCA(dds.vst, 
              intgroup = c("Treatment"),
              ntop = 500,
              returnData = FALSE,
              pcsToUse = 1:2) 

# Method 2: : ---------------------------------
# more customizable

# 1. Compute PCA 
pca <- prcomp(t(assay(dds.vst)))

# 2. Compute contribution of each component to the total variance
percentVar <- pca$sdev^2 / sum( pca$sdev^2)
pc1 <- round(percentVar*100, digits = 1)[1]
pc2 <- round(percentVar*100, digits = 1)[2]

# 3. assembly the data for the plot
pca.df <- data.frame(PC1=pca$x[,"PC1"], 
                   PC2=pca$x[,"PC2"], 
                   Treatment=dds.vst@colData$Treatment,
                   Time_id=dds.vst@colData$Time_id
                   )

# 4. Generate plot
pca.p <- ggplot(data=pca.df, aes(x=PC1, y=PC2, color=Treatment, shape=Time_id)) + 
    geom_point(size=3) + 
    xlab(paste0("PC1:",pc1,"% variance")) + 
    ylab(paste0("PC2:",pc2,"% variance")) 

# 5. Save PCA plot
ggsave(filename = "pca.pdf", plot = pca.p)

```

## 2. Heatmap

```{r}
# 1. Normalize counts
dds.vst <- vst(dds, blind=TRUE)

# 2. Add new Sample_id column to metadata in vst object
dds.vst$Sample_id <- as.factor(paste(dds.vst$Treatment,
                                     dds.vst$Time_id,
                                     dds.vst$Sbj_id, 
                                     sep = "_"))


# 3. Calculate distances between samples
sampleDists <- dist(t(assay(dds.vst)))

# 4. Plot inter-sample distances
sampleDist.mat <- as.matrix(sampleDists)
rownames(sampleDist.mat) <- paste(rownames(sampleDist.mat), dds.vst$Treatment)

hm.all.p <- pheatmap(mat = sampleDist.mat,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                labels_row = dds.vst@colData$Sample_id,
                labels_col = dds.vst@colData$Time_id,
                fontsize_row = 8,
                fontsize_col = 8
                )
# 6. Save heatmap plot
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

To control for dependency of Log2FC with gene expression level.

```{r message=FALSE, warning=FALSE}
plotMA(res.2, ylim=c(-3,3))

res.shk <- lfcShrink(dds, contrast = c("Treatment", "Drug", "Control"), type = "ashr")

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
