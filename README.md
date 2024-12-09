SRSTCM: Single-cell Recognition of Senescent Tumor Cells based on Multinomial-distribution Model and Maximum Likelihood Estimation 
===
Senescence has been listed as one of the cancer hallmarks, and numerous anti-tumor treatments could cause tumor cell senescence, including radiotherapy, chemotherapy, and targeted therapies. Meanwhile, senescent tumor cells (STCs) have been found in diverse cancer stages. STCs could induce immunosuppression of cancer, and then promote the development of tumor. It is therefore critical to understand how STCs contribute to cancer progression by identifying STCs in cancer. To solve this problem, we proposed SRSTCM, a supervised STCs identifier that could accurately identify STCs in cancer single-cell RNA-seq data based on Multinomial-distribution Model and Maximum Likelihood Estimation. In this vignette, we will predict senescent tumor cells by calculating 10X single-cell RNA data using the {SRSTCM} R package.<br>

![image text](https://github.com/yjybio/SRSTCM/blob/master/workflows/workflows.png)

Step 1: installation
--
Installing SRSTCM from GitHub <br>
```R
library(devtools) 
install_github("yjybio/SRSTCM")
```
We start by loading the packages needed for the analyses. Please install them if you haven't.<br>
```R
library("copykat")
library("dplyr")
library("ggplot2")
library("harmony")
library("Hmisc")
library("patchwork")
library("pheatmap")
library("Seurat")
library("SeuratObject")
library("survival")
library("survminer")
```
Step 2: single-cell RNA-seq transcriptome data preprocessing
--
Single-cell RNA-seq transcriptome data should be preprocessed before identifying senescent tumor cells in single-cell RNA-seq transcriptome data. The single-cell RNA-seq transcriptome data read by default is UMI count matrix generated from 10X output.<br>

An example input 10x raw data matrix
```R
data <- subset(data, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
data <- NormalizeData(data) %>% FindVariableFeatures(nfeatures = ceiling(nrow(data@assays$RNA) * 0.3)) %>% ScaleData()
data <- RunPCA(data, features = VariableFeatures(object = data), reduction.name = "pca")
data <- RunUMAP(data, reduction = "pca", dims = 1:10, reduction.name = "umap")
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, res = 0.5)
```
Step 3: identify malignant cells in single-cell RNA-seq transcriptome data
--
```R
tumor.cell <- as.matrix(data@assays$RNA@counts)
copykat.test <- copykat(rawmat = tumor.cell, sam.name = cancer)
sample.name <- paste(cancer, "_copykat_", sep = "")
res <- paste(sample.name, "prediction.txt", sep = "")
copykat.prediction <- read.table(res)
copykat.prediction <- copykat.prediction[which(copykat.prediction[, 2] == "aneuploid"), ]
copykat.prediction <- data@assays$RNA@counts[ ,which(colnames(data@assays$RNA@counts) %in% copykat.prediction[, 1])]
```
Step 4: running SRSTCM
--
Taking Ovarian cancer as an example, the most important input data is the malignant cells obtained in the previous step, and the malignant cells mentioned above are extracted from the preprocessed single-cell RNA-seq transcriptome data and recorded as copykat-prediction.<br>
Then, the corresponding senescence-related genes in Ovarian cancer were extracted from the data, and senescent tumor cells in malignant cells were identified according to the maximum likelihood function principle according to the multinomial-distribution model of Ovarian cancer constructed by us.<br>

Now run the code:
```R
copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% OVCA.sen.gene), ]
dir <- which(OVCA.sen.gene %in% rownames(copykat.prediction))
if (length(dir) == 0 || length(dir) == 1)
     stop("too few senescence-related genes;cannot be calculated")
if (length(dir) > 1) {
     copykat.prediction <- copykat.prediction[match(OVCA.sen.gene[dir], rownames(copykat.prediction)), ]
     OVCA.sen.average <- OVCA.sen.average[dir]
     OVCA.no.sen.average <- OVCA.no.sen.average[dir]
     p.no.sen <- ( OVCA.no.sen.average + 1 ) / ( sum(OVCA.no.sen.average) + length(OVCA.no.sen.average) )
     p.sen <- ( OVCA.sen.average + 1 ) / ( sum(OVCA.sen.average) + length(OVCA.sen.average) )
     q.OVCA <- apply(copykat.prediction, 2, q)
     Senescence <- c()
     for (i in 1 : ncol(q.OVCA)) {
          if (argmax(p.sen, q.OVCA[, i]) > argmax(p.no.sen, q.OVCA[, i])) {
                 Senescence <- c(Senescence, colnames(q.OVCA)[i])
           }
     }
}
```
I can save identified senescent tumor cells.
```R
 write.table(Senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
```     
Annotate the position of senescent tumor cells in the single-cell RNA-seq transcriptome data cluster map.<br>
```R
data@meta.data$cell.type <- data@meta.data$seurat_clusters
data@meta.data$cell.type <- as.character(data@meta.data$cell.type)
data@meta.data$cell.type[which(rownames(data@meta.data) %in% aging)] <- "STCs"
data@meta.data$cell.type <- as.factor(data@meta.data$cell.type)
DimPlot(data, reduction = "umap", label = TRUE, group.by = "cell.type", label.size = 3)
```
Identify differentially expressed genes between senescent tumor cells and non-senescent tumor cells.
```R
data@meta.data$cell.type <- as.character(data@meta.data$cell.type)
no.sen <- colnames(copykat.prediction)[-which(colnames(copykat.prediction) %in% sen)]
data@meta.data$cell.type[which(rownames(data@meta.data) %in% no.sen)] <- "NSTCs"
data@meta.data$cell.type <- as.factor(data@meta.data$cell.type)
Idents(data) <- data@meta.data$cell.type
mydeg <- FindMarkers(data, ident.1 = 'STCs', ident.2 = 'NSTCs', verbose = FALSE, test.use = 'wilcox', min.pct = 0.1)
write.csv(mydeg, paste(cancer, "_STCs_NSTCs_deg.csv", sep = ""), row.names = TRUE)
```
Plot the expression levels of the top ten differentially expressed genes in senescent tumor cells.
```R
top10 <- mydeg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
VlnPlot(data, features = top10, split.by = 'cell.type', idents = 'STCs')  
```
Plot the expression levels of the top ten differentially expressed genes in non-senescent tumor cells.
```R
top10 <- mydeg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
VlnPlot(data, features = top10, split.by = 'cell.type', idents = 'NSTCs')  
```
Step 5: if the data type is bulk, the senescence score is calculated
--
For bulk data, we can calculate the senescence score of each sample according to the senescence multinomial-distribution model. The likelihood function of the senescence multinomial distribution model corresponding to each sample is used as the score to evaluate the senescence degree of the sample.<br>

Take skin cutaneous melanoma, for example:<br>
```R
data <- data[which(rownames(data) %in% SKCM.sen.gene), ]
dir <- which(SKCM.sen.gene %in% rownames(data))
if (length(dir) == 0 || length(dir) == 1)
	stop("too few senescence-related genes;cannot be calculated")
if (length(dir) > 1) {
        data <- data[match(SKCM.sen.gene[dir], rownames(data)), ]
        SKCM.sen.average <- SKCM.sen.average[dir]
        p.sen <- ( SKCM.sen.average + 1 ) / ( sum(SKCM.sen.average) + length(SKCM.sen.average) )
        q.SKCM <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.SKCM), ncol = 2)
	for (i in 1 : ncol(q.SKCM)) {
             Senescence.score[i, 1] <- colnames(q.SKCM)[i]
             Senescence.score[i, 2] <- argmax(p.sen, q.SKCM[, i])
         }
        colnames(Senescence.score) <- c("sample", "senescence_score")
	Senescence.score <- as.data.frame(Senescence.score)
}
```
I can save the aging score.
```R
write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
```
If the survival data is not empty, we can calculate the prognostic survival curves of the two groups with high and low senescence scores.<br>
```R
Senescence.score$senescence_score <- as.numeric(Senescence.score$senescence_score)
Senescence.score <- Senescence.score[match(sort(Senescence.score$senescence_score), Senescence.score$senescence_score), ]
survival <- survival[survival$sample %in% Senescence.score$sample, ]
survival <- survival[match(Senescence.score$sample, survival$sample), ]
survival$senescence_score <- Senescence.score$senescence_score
count <- ceiling(nrow(survival) / 2)
class <- c(rep("low_score", count), rep("high_score", (nrow(survival) - count)))
survival$class <- class
fit <- survfit(Surv(OS.time, OS) ~ class, data = survival)
ggsurvplot(fit,
           data = survival,        
           palette = c("#E7B800", "#2E9FDF"),   
           conf.int = FALSE,  
           pval = T,         
           pval.method = T,   
           risk.table = T,    
           cumevents = F,    
           cumcensor = F,    
)
```
