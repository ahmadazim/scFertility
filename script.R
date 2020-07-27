
library(Seurat)
library(reticulate)
list.files("Data/GSM2928377_HumanSpermatogonia_17-1")
list.files("Data/GSM2928377_HumanSpermatogonia_17-2")
list.files("Data/GSM2928377_HumanSpermatogonia_17-5")
list.files("Data/GSM2928380_AdultHuman_17-3")
list.files("Data/GSM2928380_AdultHuman_17-4")
list.files("Data/GSM2928380_AdultHuman_17-5")
list.files("Data/GSM2928383_AdultHuman-Spermatocytes_17-6")
list.files("Data/GSM2928384_AdultHuman-Spermatids_17-6")
list.files("Data/GSM2928385_Human-Mouse-multiplet-hg19")
list.files("Data/GSM2928385_Human-Mouse-multiplet-mm10")


gonia1 <- Read10X(data.dir = "Data/GSM2928377_HumanSpermatogonia_17-1")
gonia2 <- Read10X(data.dir = "Data/GSM2928377_HumanSpermatogonia_17-2")
gonia5 <- Read10X(data.dir = "Data/GSM2928377_HumanSpermatogonia_17-5")
ad3 <- Read10X(data.dir = "Data/GSM2928380_AdultHuman_17-3")
ad4 <- Read10X(data.dir = "Data/GSM2928380_AdultHuman_17-4")
ad5 <- Read10X(data.dir = "Data/GSM2928380_AdultHuman_17-5")
adSpermatocytes <- Read10X(data.dir = "Data/GSM2928383_AdultHuman-Spermatocytes_17-6")
adSpermatids <- Read10X(data.dir = "Data/GSM2928384_AdultHuman-Spermatids_17-6")
mpHG19 <- Read10X(data.dir = "Data/GSM2928385_Human-Mouse-multiplet-hg19")
mpMM10 <- Read10X(data.dir = "Data/GSM2928385_Human-Mouse-multiplet-mm10")


#Creating the Seurat object
ad3.obj <- CreateSeuratObject(counts = ad3, min.cells = 3, min.features = 200)
ad4.obj <- CreateSeuratObject(counts = ad4, min.cells = 3, min.features = 200)
ad5.obj <- CreateSeuratObject(counts = ad5, min.cells = 3, min.features = 200)

adult.unnorm <- merge(ad3.obj, y = c(ad4.obj, ad5.obj), add.cell.ids = c("R3", "R4", "R5"), project = "sperm")

#Normalizing the data
adult <- NormalizeData(adult.unnorm, normalization.method = "LogNormalize", scale.factor = 10000)

#Examining nFeature_RNA, nCount_RNA, and percent.mt for quality control
adult[["percent.mt"]] <- PercentageFeatureSet(adult, pattern = "^MT-")
VlnPlot(adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
adult <- subset(adult, subset = nFeature_RNA > 200 & percent.mt < 50)

library(sctransform)
adult <- SCTransform(adult, verbose = TRUE)

adult <- FindVariableFeatures(adult, selection.method = "vst", nfeatures = 5000)

adult <- RunPCA(adult, verbose = TRUE)
adult <- RunTSNE(adult, verbose = TRUE)

adult <- JackStraw(adult, num.replicate = 100)
adult <- ScoreJackStraw(adult, dims = 1:20)
JackStrawPlot(adult, dims = 1:20)

adult <- FindNeighbors(adult, dims = 1:18, verbose = TRUE)
adult <- FindClusters(adult, verbose = TRUE, resolution = 0.6)

#Find markers for each cluster
adult.markers <- FindAllMarkers(adult, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
library(dplyr)
top10 <- adult.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top50 <- adult.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

#Vizulizing clusters
png("./fig1.png", units = "in", width = 10, height = 7, res = 500)
DimPlot(adult, reduction = "tsne", label = TRUE, label.size = 5)
dev.off()

# Upload clusterID.csv file
clustMark <- read.csv("Spreadsheets/clustMark.csv", header = T)
clustMark50 <- clustMark %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
x <- table(rep(0:20))
for(i in 1:16){
  x <- rbind(x, round(table(adult.markers$cluster[adult.markers$gene %in% as.character(clustMark$gene[clustMark$cluster == i])])/table(adult.markers$cluster), digits = 2))
}
x <- x[2:17,]
table(top50$cluster[top50$gene %in% as.character(clustMark50$gene[clustMark50$cluster == i])])   # i from 1 to 16

## Cluster ID's Found...
#Cluster0 = 2/13
#Cluster1 = 8
#Cluster2 = 10
#Cluster3 = 6
#Cluster4 = 4
#Cluster5 = 12
#Cluster6 = 5
#Cluster7 = 5
#Cluster8 = 4
#Cluster9 = 3                           # Interested in clusters that matched to clusters
#Cluster10 = 7                            15 and 16, which are my clusters 17-20
#Cluster11 = 11
#Cluster12 = 9
#Cluster13 = 14
#Cluster14 = 7
#Cluster15 = 1
#Cluster16 = 2
#Cluster17 = 16
#Cluster18 = 15
#Cluster19 = 16
#Cluster20 = 15

cells <- names(adult$seurat_clusters[adult$seurat_clusters == 17|adult$seurat_clusters == 18|adult$seurat_clusters == 19|adult$seurat_clusters == 20])
finalClust <- subset(adult.unnorm, cells = cells)

finalClust <- NormalizeData(finalClust, normalization.method = "LogNormalize", scale.factor = 10000)
finalClust[["percent.mt"]] <- PercentageFeatureSet(finalClust, pattern = "^MT-")
VlnPlot(finalClust, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
finalClust <- subset(finalClust, subset = nFeature_RNA < 7500 & percent.mt < 30)
library(sctransform)
finalClust <- SCTransform(finalClust, verbose = TRUE)
finalClust <- FindVariableFeatures(finalClust, selection.method = "vst", nfeatures = 5000)

finalClust <- RunPCA(finalClust, verbose = TRUE)
finalClust <- RunTSNE(finalClust, verbose = TRUE)
finalClust <- JackStraw(finalClust, num.replicate = 100)
finalClust <- ScoreJackStraw(finalClust, dims = 1:20)
JackStrawPlot(finalClust, dims = 1:20)
finalClust <- FindNeighbors(finalClust, dims = 1:8, verbose = TRUE)
finalClust <- FindClusters(finalClust, verbose = TRUE, resolution = 0.9)

finalClust.markers <- FindAllMarkers(finalClust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
library(dplyr)
top10 <- finalClust.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top50 <- finalClust.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

png("./fig2.png", units = "in", width = 10, height = 7, res = 500)
DimPlot(finalClust, reduction = "tsne", label = TRUE, label.size = 5)
dev.off()

is.element("CA2", row.names(finalClust))
FeaturePlot(finalClust, feature= "CA2")

is.element("TNP1", row.names(finalClust))
is.element("TNP2", row.names(finalClust))
FeaturePlot(finalClust, feature= c("TNP1", "TNP2"))

is.element("PRM1", row.names(finalClust))
is.element("PRM2", row.names(finalClust))
FeaturePlot(finalClust, feature= c("PRM1", "PRM2"))

top10 <- as.data.frame(top10)
top10 <- top10[top10$pct.2 < 0.1,]
