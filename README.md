ScrSTCs: Inference of human senescent tumor cells from high-throughput single cell RNA-seq data based on senescence marker genes
===
More and more evidence shows that senescent tumor cells (STCs) can induce immunosuppression of cancer and promote the development of tumor. It is therefore critical to understand how senescent tumor cells contribute to disease progression by identifying senescent tumor cells in a disease or condition. To address this problem, we propose ScrSTCs, a supervised senescent tumor cell identifier that can accurately identify senescent tumor cells in newly sequenced cancer single-cell data. ScrSTCs annotates senescent tumor cells in cancer single-cell data. We use a multinomial distribution model and maximum likelihood function estimation to develop ScrSTCs for accurate, rapid, and robust identification of senescent tumor cells. Due to the limited senescence induced cancer RNA-seq dataset, we highlighted senescent tumor cells in the single cell transcriptomic data for six cancers: lung cancer, colorectal cancer, breast cancer, ovarian cancer, hepatocellular carcinoma, and cutaneous melanoma. In addition, for bulk data, we can calculate the aging score of each sample according to the multinomial distribution model. The likelihood function of the aging multinomial distribution model corresponding to each sample is used as the score to evaluate the aging degree of the sample. In this vignette, we will predict senescent tumor cells by calculating 10X single-cell RNA data using the {ScrSTCs} R package.<br>

Step 1: installation
--
Installing ScrSTCs from GitHub <br>
```R
library(devtools) 
install_github("yjybio/ScrSTCs")
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
Step 2: single-cell transcriptome data preprocessing
--
Single-cell sequencing data should be preprocessed before identifying senescent tumor cells in single-cell data. The single-cell sequencing data read by default is UMI count matrix generated from 10X output.<br>

An example input 10x raw data matrix
```R
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
data <- subset(data, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
data <- NormalizeData(data) %>% FindVariableFeatures(nfeatures = ceiling(nrow(data@assays$RNA) * 0.3)) %>% ScaleData()
data <- RunPCA(data, features = VariableFeatures(object = data), reduction.name = "pca")
data <- RunUMAP(data, reduction = "pca", dims = 1:10, reduction.name = "umap")
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, res = 0.5)
```
Step 3: identify malignant cells in single-cell sequencing data
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
Step 4: running ScrSTCs
--
Taking colorectal cancer as an example, the most important input data is the malignant cells obtained in the previous step, and the malignant cells mentioned above are extracted from the preprocessed single-cell data and recorded as copykat-prediction.<br>
Then, the corresponding aging marker genes in colorectal cancer were extracted from the data, and senescent tumor cells in malignant cells were identified according to the maximum likelihood function principle according to the multinomial distribution model of colorectal cancer constructed by us.<br>

Now run the code:
```R
copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% CRC.aging.gene), ]
dir <- which(CRC.aging.gene %in% rownames(copykat.prediction))
if (length(dir) == 0 || length(dir) == 1)
     stop("too few aging marker genes;cannot be calculated")
if (length(dir) > 1) {
     copykat.prediction <- copykat.prediction[match(CRC.aging.gene[dir], rownames(copykat.prediction)), ]
     CRC.aging.average <- CRC.aging.average[dir]
     CRC.no.aging.average <- CRC.no.aging.average[dir]
     p.no.aging <- ( CRC.no.aging.average + 1 ) / ( sum(CRC.no.aging.average) + length(CRC.no.aging.average) )
     p.aging <- ( CRC.aging.average + 1 ) / ( sum(CRC.aging.average) + length(CRC.aging.average) )
     q.CRC <- apply(copykat.prediction, 2, q)
     aging <- c()
     for (i in 1 : ncol(q.CRC)) {
          if (argmax(p.aging, q.CRC[, i]) > argmax(p.no.aging, q.CRC[, i])) {
                 aging <- c(aging, colnames(q.CRC)[i])
           }
     }
}
```
I can save identified senescent tumor cells.
```R
 write.table(aging, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
```     
Annotate the position of senescent tumor cells in the single-cell cluster map.<br>
```R
data@meta.data$cell.type <- data@meta.data$seurat_clusters
data@meta.data$cell.type <- as.character(data@meta.data$cell.type)
data@meta.data$cell.type[which(rownames(data@meta.data) %in% aging)] <- "STCs"
data@meta.data$cell.type <- as.factor(data@meta.data$cell.type)
DimPlot(data, reduction = "umap", label = TRUE, group.by = "cell.type", label.size = 3)
```
Identify differentially expressed genes between senescent tumor cells and non-aging tumor cells.
```R
data@meta.data$cell.type <- as.character(data@meta.data$cell.type)
no.aging <- colnames(copykat.prediction)[-which(colnames(copykat.prediction) %in% aging)]
data@meta.data$cell.type[which(rownames(data@meta.data) %in% no.aging)] <- "NSTCs"
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
Plot the expression levels of the top ten differentially expressed genes in non-aging tumor cells.
```R
top10 <- mydeg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
VlnPlot(data, features = top10, split.by = 'cell.type', idents = 'NSTCs')  
```
Step 5: if the data type is bulk, the aging score is calculated
--
For bulk data, we can calculate the aging score of each sample according to the multinomial distribution model. The likelihood function of the aging multinomial distribution model corresponding to each sample is used as the score to evaluate the aging degree of the sample.<br>

Take colorectal cancer, for example:<br>
```R
data <- data[which(rownames(data) %in% CRC.aging.gene), ]
dir <- which(CRC.aging.gene %in% rownames(data))
if (length(dir) == 0 || length(dir) == 1)
	stop("too few aging marker genes;cannot be calculated")
if (length(dir) > 1) {
        data <- data[match(CRC.aging.gene[dir], rownames(data)), ]
        CRC.aging.average <- CRC.aging.average[dir]
        p.aging <- ( CRC.aging.average + 1 ) / ( sum(CRC.aging.average) + length(CRC.aging.average) )
        q.CRC <- apply(data, 2, q)
        Aging.score <- matrix(nrow = ncol(q.CRC), ncol = 2)
	for (i in 1 : ncol(q.CRC)) {
             Aging.score[i, 1] <- colnames(q.CRC)[i]
             Aging.score[i, 2] <- argmax(p.aging, q.CRC[, i])
         }
        colnames(Aging.score) <- c("sample", "aging_score")
	Aging.score <- as.data.frame(Aging.score)
}
```
I can save the aging score.
```R
write.table(Aging.score, paste(cancer, "_aging_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
```
If the survival data is not empty, we can calculate the prognostic survival curves of the two groups with high and low aging scores.<br>
```R
Aging.score$aging_score <- as.numeric(Aging.score$aging_score)
Aging.score <- Aging.score[match(sort(Aging.score$aging_score), Aging.score$aging_score), ]
survival <- survival[survival$sample %in% Aging.score$sample, ]
survival <- survival[match(Aging.score$sample, survival$sample), ]
survival$aging_score <- Aging.score$aging_score
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
