library(data.table)
library(tidyverse)
library(Seurat)



ff <- list.files()
scRNAlist <- list()

for(i in 1:4){
  f <- ff[i]
  b <- str_match(f, "HPV")[1]
  a <- data.table::fread(input = f, sep = "\t", skip = 7)
  a <- as.data.frame(a)
  a <- remove_rownames(a) %>% column_to_rownames("V1")
  colnames(a) <- paste0(b, "_", i, "_", colnames(a))

  scRNAlist[[i]] <- CreateSeuratObject(counts = a, 
                                       project = paste0(b, "_", i))
}

for(i in 1:length(scRNAlist)){
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]],
                                                         pattern = "^MT-")
  scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]],
                                                         pattern = "^RP[SL]")
}
scRNA$sample <- str_match(colnames(scRNA), "[A-Z].*?_[1-9]") %>% unlist

scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
plot.featrures <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
plots <- list()
for(i in seq_along(plot.featrures)){
  plots[[i]] <- VlnPlot(scRNA, group.by = "sample", 
                        pt.size = 0,
                        cols = c("#2c69b0", "#b5c8e2", "#f02720", "#ffb6b0", 
                                 "#e9c39b", "#6ba3d6", "#b5dffd", "#ac8763", 
                                 "#ac8763", "#ddc9b4", "#bd0a36", "#f4737a"), 
                        features = plot.featrures[i]) + NoLegend()
}
scRNA <- subset(scRNA, percent.mt < 10)

# 校正前

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- NormalizeData(object = scRNA, verbose = FALSE)
scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes)
scRNA <- SCTransform(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 5000)
sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(scRNA), value = TRUE, invert = TRUE)
scRNA <- RunPCA(scRNA, features = sele)
ElbowPlot(scRNA,ndims = 20)
show_tree <- function(data, dims = 1:15){
  library(clustree)
  test <- FindNeighbors(data, dims = dims)
  sce <- FindClusters(
    object = test, resolution = c(seq(.1,1.2,.2))
  )
  clustree(sce@meta.data, prefix = "SCT_snn_res.")
}
show_tree(scRNA, dims = 1:10)
scRNA <- FindNeighbors(scRNA, dims = 1:10)
scRNA <- FindClusters(object = scRNA, resolution = 0.5)
scRNA <- scRNA %>% RunTSNE(dims = 1:10) %>% RunUMAP(dims = 1:10)
DimPlot(scRNA, reduction = "umap", group.by = "sample")
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters")

# 矫正后

library(harmony)
scRNA <- RunHarmony(object = scRNA, group.by.vars = "sample", plot_convergence = TRUE)
ElbowPlot(scRNA,ndims = 25, reduction = "harmony")
show_tree <- function(data, dims = 1:15){
  test <- FindNeighbors(data, dims = dims, reduction = "harmony")
  sce <- FindClusters(
    object = test, resolution = c(seq(.1,1.2,.2))
  )
  clustree(sce@meta.data, prefix = "SCT_snn_res.")
}
show_tree(scRNA)
sc <- FindNeighbors(scRNA, dims = 1:15, reduction = "harmony")
sc <- FindClusters(object = sc, resolution = 0.7)
sc <- sc %>% RunTSNE(dims = 1:15, reduction = "harmony") %>% RunUMAP(dims = 1:15, reduction = "harmony")
DimPlot(sc, reduction = "umap", group.by = "sample")
DimPlot(sc, reduction = "umap", group.by = "seurat_clusters")

FeaturePlot(sc, features = c("CD4", "CD8A", "CD3D"))

if(T){
  Epithelial <- c("KRT5", "KRT6A","IL1RN", "EMP1", "KRT13")
  `T cells` <- c("PTPRC", "CD3E", "CD3D", "TRBC1", "CD4", "CD8A")
  `Myeloid cells` <- c("LYZ", "CD86", "CD68", "FCGR3A")
  `B cells` <- c("CD79A", "CD79B", "JCHAIN", "IGKC","IGHG3") 
  Fibroblasts <- c("DCN", "C1R", "COL1A1","FGF7")
  genes_to_check <- list(Epithelial, `T cells`, `Myeloid cells`,
                         `B cells`, Fibroblasts)
  names(genes_to_check) <- c("Epithelial", "T cells", "Myeloid cells",
                             "B cells", "Fibroblasts")
  DotPlot(sc, features = genes_to_check,group.by = "seurat_clusters", assay = "RNA")+ RotatedAxis() +
    scale_x_discrete("") + 
    scale_y_discrete("") + 
    scale_color_gradientn(colours = viridis(256))
}
library(viridis)


sc <- subset(sc, seurat_clusters != 11)
sc$seurat_clusters <- factor(sc$seurat_clusters)
sc$HPV <- str_split(sc$sample, "_", simplify = T)[,1]

sc$cell_type <- fct_collapse(sc$seurat_clusters, 
                             Epithelial = c("1", "3", "6", "9", "10", "14"),
                             T_cells = c("0", "2", "4", "5"),
                             B_cells = c("7", "12", "13"),
                             Myeloid = c("8"),
                             Fibroblasts = c("15")
                             )
DotPlot(sc, features = genes_to_check,group.by = "cell_type", assay = "RNA")+ RotatedAxis() +
  scale_x_discrete("")+scale_y_discrete("")
saveRDS(sc, file = "step_1.rds")
FeaturePlot(sc, features = delt)
DimPlot(sc, group.by = "cell_type")

sc <- readRDS("step_1.rds")

Idents(sc) <- sc$cell_type
sc_marker <- COSG::cosg(sc, n_genes_user = 15)
Fib <- subset(sc, cell_type == "Fibroblasts")
Fib$cell_type <- factor(Fib$cell_type)

# T细胞分群 -------------------------------------------------------------------

Tcells <- subset(sc, cell_type == "T_cells")
Tcells <- CreateSeuratObject(counts = Tcells@assays$RNA@counts, meta.data = Tcells@meta.data)
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
Tcells <- NormalizeData(object = Tcells, verbose = FALSE)
Tcells <- CellCycleScoring(Tcells, s.features = s.genes, g2m.features = g2m.genes)
Tcells <- SCTransform(Tcells, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(Tcells), value = TRUE, invert = TRUE)
Tcells <- RunPCA(Tcells, features = sele)
Tcells <- RunHarmony(object = Tcells, group.by.vars = "sample", plot_convergence = TRUE)
ElbowPlot(Tcells,ndims = 25, reduction = "harmony")
show_tree <- function(data, dims = 1:14){
  test <- FindNeighbors(data, dims = dims, reduction = "harmony")
  Tcellse <- FindClusters(
    object = test, resolution = c(seq(.1,0.8,.2))
  )
  clustree(Tcellse@meta.data, prefix = "SCT_snn_res.")
}
show_tree(Tcells)
Tcells <- FindNeighbors(Tcells, dims = 1:14, reduction = "harmony")
Tcells <- FindClusters(object = Tcells, resolution = 0.5)
Tcells <- Tcells %>% RunTSNE(dims = 1:14, reduction = "harmony") %>% RunUMAP(dims = 1:14, reduction = "harmony")
DimPlot(Tcells, reduction = "umap", group.by = "sample")
DimPlot(Tcells, reduction = "tsne", group.by = "seurat_clusters")
DimPlot(Tcells, reduction = "umap", group.by = "seurat_clusters")


FeaturePlot(Tcells, features = c("FOXP3", "CD4", "CD8A", "IL7R"))
Idents(Tcells) <- Tcells$seurat_clusters
T_marker <- COSG::cosg(Tcells, n_genes_user = 15)
DotPlot(Tcells, features = unique(unlist(T_marker$names)), group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
DotPlot(Tcells, features = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA", 
                             "FOXP3", "FCGR3A", "CD28", "TRDC", "TRGV9", "TRDV2",
                             "NKG7", "GNLY", "XCL1", "XCL2"), group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
DotPlot(Tcells, features = c(c("ENTPD1", "LAYN", "CD3D", "CD3E",
                               "HAVCR2", "LAG3", "CTLA4", "PDCD1", "GZMK", "IFNG", "FOXP3", "GNLY", "NKG7", "CD4", "CD8A")), group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
levels(Tcells$seurat_clusters) <- c("Naive_T", "cyto1_T", "Treg", "gammaT", "Tem",
                                    "Tem", "cyto2_T", "Th17", "NKT", "ISG1")
Tcells$cell_type <- Tcells$seurat_clusters
saveRDS(Tcells, file = "Tcells.rds")
Tcells <- readRDS("Tcells.rds")



# Epi细胞 -------------------------------------------------------------------

Epi <- subset(sc, cell_type == "Epithelial")
Epi <- CreateSeuratObject(counts = Epi@assays$RNA@counts, meta.data = Epi@meta.data)
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
Epi <- NormalizeData(object = Epi, verbose = FALSE)
Epi <- CellCycleScoring(Epi, s.features = s.genes, g2m.features = g2m.genes)
Epi <- SCTransform(Epi, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000)
sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(Epi), value = TRUE, invert = TRUE)
Epi <- RunPCA(Epi, features = sele)
Epi <- RunHarmony(object = Epi, group.by.vars = "sample", plot_convergence = TRUE)
ElbowPlot(Epi,ndims = 25, reduction = "harmony")
show_tree <- function(data, dims = 1:16){
  test <- FindNeighbors(data, dims = dims, reduction = "harmony")
  Epie <- FindClusters(
    object = test, resolution = c(seq(.1,0.8,.2))
  )
  clustree(Epie@meta.data, prefix = "SCT_snn_res.")
}
show_tree(Epi)
Epi <- FindNeighbors(Epi, dims = 1:16, reduction = "harmony")
Epi <- FindClusters(object = Epi, resolution = 0.3)
Epi <- Epi %>% RunTSNE(dims = 1:16, reduction = "harmony") %>% RunUMAP(dims = 1:16, reduction = "harmony")
DimPlot(Epi, reduction = "umap", group.by = "sample")
DimPlot(Epi, reduction = "tsne", group.by = "seurat_clusters")
DimPlot(Epi, reduction = "umap", group.by = "seurat_clusters")


FeaturePlot(Epi, features = c("SAA1", "SAA2"))
Idents(Epi) <- Epi$seurat_clusters
marker <- COSG::cosg(Epi, n_genes_user = 10, expressed_pct = 0.2)
DotPlot(Epi, features = unique(unlist(marker$names)), group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
Epi <- subset(Epi, seurat_clusters != "5")
Epi$seurat_clusters <- factor(Epi$seurat_clusters)
levels(Epi$seurat_clusters) <- paste0("Epi_", 1:5)
Epi$cell_type <- Epi$seurat_clusters
Idents(Epi) <- Epi$cell_type
Epi_marker <- COSG::cosg(Epi, n_genes_user = 20)
Epi_marker$names$Epi_5
saveRDS(Epi, file = "Epi.rds")
Epi <- readRDS("Epi.rds")


# My细胞 --------------------------------------------------------------------

My <- subset(sc, cell_type == "Myeloid")
My <- CreateSeuratObject(counts = My@assays$RNA@counts, meta.data = My@meta.data)
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
My <- NormalizeData(object = My, verbose = FALSE)
My <- CellCycleScoring(My, s.features = s.genes, g2m.features = g2m.genes)
My <- SCTransform(My, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
My <- FindVariableFeatures(My, selection.method = "vst", nfeatures = 2000)
sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(My), value = TRUE, invert = TRUE)
My <- RunPCA(My, features = sele)
My <- RunHarmony(object = My, group.by.vars = "sample", plot_convergence = TRUE)
ElbowPlot(My,ndims = 25, reduction = "harmony")
show_tree <- function(data, dims = 1:16){
  test <- FindNeighbors(data, dims = dims, reduction = "harmony")
  Mye <- FindClusters(
    object = test, resolution = c(seq(.1,0.8,.2))
  )
  clustree(Mye@meta.data, prefix = "SCT_snn_res.")
}
show_tree(My)
My <- FindNeighbors(My, dims = 1:16, reduction = "harmony")
My <- FindClusters(object = My, resolution = 0.3)
My <- My %>% RunTSNE(dims = 1:16, reduction = "harmony") %>% RunUMAP(dims = 1:16, reduction = "harmony")
DimPlot(My, reduction = "umap", group.by = "sample")
DimPlot(My, reduction = "tsne", group.by = "seurat_clusters")
DimPlot(My, reduction = "umap", group.by = "seurat_clusters")

FeaturePlot(My, features = c("SAA1", "SAA2"))
Idents(My) <- My$seurat_clusters
marker <- COSG::cosg(My, n_genes_user = 20, expressed_pct = 0.2)
DotPlot(My, features = unique(unlist(marker$names)), group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
My <- subset(My, seurat_clusters %in% c(0, 1, 2))
My$seurat_clusters <- factor(My$seurat_clusters)
levels(My$seurat_clusters) <- c("pDC", "TAM", "monocyte")
My$cell_type <- My$seurat_clusters
Idents(My) <- My$cell_type
My_marker <- COSG::cosg(My, n_genes_user = 15)
saveRDS(My, file = "My.rds")
My <- readRDS("My.rds")




# B细胞 ---------------------------------------------------------------------

Bcells <- subset(sc, cell_type == "B_cells")
Bcells <- CreateSeuratObject(counts = Bcells@assays$RNA@counts, meta.data = Bcells@meta.data)
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
Bcells <- NormalizeData(object = Bcells, verbose = FALSE)
Bcells <- CellCycleScoring(Bcells, s.features = s.genes, g2m.features = g2m.genes)
Bcells <- SCTransform(Bcells, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
Bcells <- FindVariableFeatures(Bcells, selection.method = "vst", nfeatures = 2000)
sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(Bcells), value = TRUE, invert = TRUE)
Bcells <- RunPCA(Bcells, features = sele)
Bcells <- RunHarmony(object = Bcells, group.by.vars = "sample", plot_convergence = TRUE)
ElbowPlot(Bcells,ndims = 25, reduction = "harmony")
show_tree <- function(data, dims = 1:15){
  test <- FindNeighbors(data, dims = dims, reduction = "harmony")
  Bcellse <- FindClusters(
    object = test, resolution = c(seq(.1,0.8,.2))
  )
  clustree(Bcellse@meta.data, prefix = "SCT_snn_res.")
}
show_tree(Bcells)
Bcells <- FindNeighbors(Bcells, dims = 1:15, reduction = "harmony")
Bcells <- FindClusters(object = Bcells, resolution = 0.3)
Bcells <- Bcells %>% RunTSNE(dims = 1:15, reduction = "harmony") %>% RunUMAP(dims = 1:15, reduction = "harmony")
DimPlot(Bcells, reduction = "umap", group.by = "sample")
DimPlot(Bcells, reduction = "tsne", group.by = "seurat_clusters")
DimPlot(Bcells, reduction = "umap", group.by = "seurat_clusters")

Idents(Bcells) <- Bcells$seurat_clusters
marker <- COSG::cosg(Bcells, n_genes_user = 60, expressed_pct = 0.2)
DotPlot(Bcells, features = marker$names$`2`, group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()

DotPlot(Bcells, features = unique(unlist(marker$names)), group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
Bcells <- subset(Bcells, seurat_clusters != "4")
Bcells$seurat_clusters <- factor(Bcells$seurat_clusters)
Bcells$cell_type <- factor(Bcells$seurat_clusters, levels = paste0("Bcells_", 1:4))
saveRDS(Bcells, file = "Bcells.rds")

mar <- c("TCL1A", "CD79B", "CD79A", "MS4A1", "MKI67")
DotPlot(Bcells, features = mar, group.by = "seurat_clusters", assay = "RNA") + RotatedAxis()
# B cells (general) # 0, 1
mar <- c("MS4A1", "CD79A", "CD79B", "VPREB3", "CD19")
# B_activated # 0
mar <- c("CD69", "CD83", "JUN", "IRF4", "MYC")
# Plasma B cells # 1
mar <- c("JCHAIN", "IGHGP", "IGKC", "IGHM", "IGHG3", "IGKV1D-39",
         "IGHG2", "IGLC3", "IGHG1", "IGHG4", 
         "IGHV4-4", "XBP1", "CD79A", "CD38", "MZB1", 
         "SSR4", "DERL3", "FKBP11", "PRDX4")
Bcells <- subset(Bcells, seurat_clusters %in% c(0, 1))
Bcells$seurat_clusters <- factor(Bcells$seurat_clusters)
levels(Bcells$seurat_clusters) <- c("B_activated", "B_Plasma")
Bcells$cell_type <- Bcells$seurat_clusters
saveRDS(Bcells, file = "B_cells.rds")
Bcells <- readRDS("B_cells.rds")
Idents(Bcells) <- Bcells$cell_type
B_marker <- COSG::cosg(Bcells, n_genes_user = 15)


mer <- merge(Tcells, list(Epi, My, Bcells, Fib))

mer <- NormalizeData(object = mer, verbose = FALSE)
mer <- CellCycleScoring(mer, s.features = s.genes, g2m.features = g2m.genes)
mer <- SCTransform(mer, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
mer <- FindVariableFeatures(mer, selection.method = "vst", nfeatures = 2000)
sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(mer), value = TRUE, invert = TRUE)
mer <- RunPCA(mer, features = sele)
mer <- RunHarmony(object = mer, group.by.vars = "sample", plot_convergence = TRUE)
mer <- mer %>% RunTSNE(dims = 1:10, reduction = "harmony") %>% RunUMAP(dims = 1:10, reduction = "harmony")
mer$cells <- fct_collapse(factor(mer$cell_type),
                          T_cells = c("cyto1_T", "cyto2_T", "gammaT", "ISG1", "Naive_T", "NKT",  "Tem", "Th17", "Treg"),
                          Myeloid = c("pDC", "TAM", "monocyte"),
                          B_cells = c("B_activated", "B_Plasma"),
                          Epithelial = c("Epi_1", "Epi_2", "Epi_3", "Epi_4", "Epi_5"),
                          Fibroblasts = c("Fibroblasts"))
saveRDS(mer, file = "merge.sc.rds")
dim(mer)
table(mer$cell_type)

Idents(mer) <- mer$cell_type 
expr <- AverageExpression(mer, assays = "RNA", slot = "count")[[1]]
col_lis <- c("Naive_T", "cyto1_T", "Treg", "gammaT", 
         "Tem", "cyto2_T", "Th17", "NKT", "ISG1", paste0("Epi_", 1:5), c("pDC", "TAM", "monocyte"),
         c("B_activated", "B_Plasma"), "Fibroblasts")
row_list <- c(unlist(T_marker$names), unlist(Epi_marker$names), unlist(My_marker$names),
              unlist(B_marker$names), sc_marker$names$Fibroblasts)

expr <- expr[row_list, col_lis]

library(ComplexHeatmap)
library(RColorBrewer)
m <- t(scale(t(expr), center = T))
mm <- m
mm[mm < -1] <- -1
Heatmap(mm, cluster_rows = F, cluster_columns = F, show_row_names = F,
        col = viridis::magma(10)
)



mer$cell_type
dat <- dplyr::select(mer@meta.data, cell_type) %>% group_by(cell_type) %>% summarise(n = n())
dat$cell_type <- factor(dat$cell_type, levels = col_lis)
dat <- dat[order(dat$cell_type),]
dat$type <- c(rep("Tcells", 9), rep("Epithelith", 5), rep("Bcells", 2), rep("My", 3), rep("Fib", 1))
ggplot(dat, aes(cell_type, n, fill = type)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  RotatedAxis() 

# 总的细胞比例

library(tidydr)
library(scales)
library(Seurat)
library(ggplot2)
library(ggalluvial)

sc <- readRDS("merge.sc.rds")
table(sc$cell_type)
color_cluster=c("#AFD4A0FF", "#1A7832FF")
sc$group <- case_when(sc$sample %in% c("HPV_1", "HPV_2") ~ "HPV_pos",
                      sc$sample %in% c("HPV_3", "HPV_4") ~ "HPV_neg")
scRNA <- sc
Idents(scRNA) <- scRNA$cell_type
table(Idents(scRNA), scRNA$group)

Cellratio <- prop.table(table(scRNA$group, Idents(scRNA)), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("sample","celltype","ratio")

Cellratio <- prop.table(table(Idents(scRNA), scRNA$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","sample","ratio")
col_lis <- c("Naive_T", "cyto1_T", "Treg", "gammaT", 
             "Tem", "cyto2_T", "Th17", "NKT", "ISG1", paste0("Epi_", 1:5), c("pDC", "TAM", "monocyte"),
             c("B_activated", "B_Plasma"), "Fibroblasts")
Cellratio$celltype <- factor(Cellratio$celltype, levels = col_lis)
p <- ggplot(Cellratio,aes(x=celltype,y=ratio,fill=sample,stratum=sample,alluvium=sample))+
  geom_col(width = 0.6,color=NA)+
  geom_flow(width=0.6,alpha=0.2,knot.pos=0) +
  scale_fill_manual(values=color_cluster) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 45, hjust = 1)); p
theme(legend.position = 'none',); p  

# gammaT 细胞比例

dat <- data.frame(
  stringsAsFactors = FALSE,
  patient = c("HPV+.1", "HPV+.2", "HPV-.1", "HPV-.2"),
  Tcells = c(1391L, 2061L, 1963L, 1613L),
  gammaT = c(77L, 94L, 225L, 315L)
)
ggplot(data = dat) +
  geom_col(aes(Tcells, patient), width = 0.8) +
  geom_col(aes(gammaT, patient, fill = patient), width = 0.3) + 
  coord_flip()



options(
  ggplot2.discrete.colour = ggsci::scale_colour_cosmic,
  ggplot2.discrete.fill = ggsci::scale_fill_lancet
)


DimPlot(Tcells, group.by = "seurat_clusters")
DimPlot(Bcells, group.by = "seurat_clusters")
DimPlot(Epi, group.by = "seurat_clusters")
DimPlot(My, group.by = "seurat_clusters")
DimPlot(mer, group.by = "cells")
DimPlot(mer, group.by = "sample")

# 拟时序分析 -------------------------------------------------------------------


library(monocle)
library(tidyverse)
test <- Tcells
data <- as(as.matrix(test@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = test@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
test <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
test <- estimateSizeFactors(test) 
test <- estimateDispersions(test)

library(Seurat)
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 3000)
sele <- grep("(^MT[0-9])|(^AL[0-9])|(^AC[0-9])", VariableFeatures(Tcells), value = TRUE, invert = TRUE)
marker <- COSG::cosg(Tcells, n_genes_user = 200)
cds <- setOrderingFilter(test, ordering_genes = unique(unlist(marker$names)))
cds <- reduceDimension(cds,reduction_method = "DDRTree",
                    max_components = 2)

cds <- orderCells(cds)
plot_cell_trajectory(cds,color_by = "cell_type", size = 1, show_branch_points = T)
plot_cell_trajectory(cds,color_by = "Pseudotime", size = 1, show_branch_points = T)
plot_cell_trajectory(cds,color_by = "State", size = 2, show_branch_points = T)
table(cds$State, cds$cell_type)


# gamma细胞 ---------------------------------------------------------------

library(monocle)
gammaT <- subset(Tcells, cell_type == "gammaT")
test <- gammaT
data <- as(as.matrix(test@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = test@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
test <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
test <- estimateSizeFactors(test) 
test <- estimateDispersions(test)

gammaT <- FindVariableFeatures(gammaT, selection.method = "vst", nfeatures = 2000)
sele <- grep("(^MT[0-9])|(^AL[0-9])|(^AC[0-9])", VariableFeatures(gammaT), value = TRUE, invert = TRUE)

cds <- setOrderingFilter(test, ordering_genes = sele)
cds <- reduceDimension(cds, reduction_method = "DDRTree",
                       max_components = 2)
cds <- orderCells(cds)
saveRDS(cds, file = "cds.rds")
cds <- readRDS(file = "cds.rds")
plot_cell_trajectory(cds,color_by = "cell_type", size = 1, show_branch_points = T)
plot_cell_trajectory(cds,color_by = "Pseudotime", size = 1, show_branch_points = T) +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62))

plot_cell_trajectory(cds,color_by = "State", size = 2, show_branch_points = T)
plot_cell_trajectory(cds,color_by = "sample", size = 2, show_branch_points = T)
table(cds$State, cds$sample)

plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "sample")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "State")

pData(cds)$NKG7 <- exprs(cds)["NKG7",]
pData(cds)$IL7R <- exprs(cds)["IL7R",]
pData(cds)$CD8A <- exprs(cds)["CD8A",]
plot_cell_trajectory(cds,color_by = "CD8A", size = 2, show_branch_points = T)

BEAM_res2 <- BEAM(cds,branch_point = 1,cores = 10)
BEAM_res3 <- BEAM_res2[,c("gene_short_name","pval","qval")]
save(BEAM_res2, BEAM_res3, file = "BEAM_1.Rdata")

plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res3,qval < 0.05)),], 
                            branch_point = 2, num_clusters = 2, 
                            cores = 1,
                            branch_labels = c("Cell fate 1", "Cell fate 2"),
                            hmcols = colorRampPalette(brewer.pal(9, "Spectral"))(50),
                            #hmcols = paletteer_c("scico::berlin", n = 50),
                            use_gene_short_name = T, show_rownames = T)

tmp <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res3,qval<0.05)),], 
                            branch_point = 2, num_clusters = 2, 
                            cores = 1,
                            branch_labels = c("Cell fate 1", "Cell fate 2"),
                            hmcols = colorRampPalette(brewer.pal(9, "Spectral"))(50),
                            #hmcols = paletteer_c("scico::berlin", n = 50),
                            use_gene_short_name = T,
                            return_heatmap = T)
tmp$annotation_row

keygenes <- c("CD8A", "XCL1", "XCL2", "NR4A1", "EGR1", "NKG7", "GNLY")
cds_subset <- cds[keygenes,]
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
p1|p2
p1 <- plot_genes_jitter(cds_subset, grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds_subset, grouping = "State", color_by = "State") + theme_minimal()
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62))
p1|p2|p3

# 分支1 ---------------------------------------------------------------------

BEAM_res4 <- BEAM(cds,branch_point = 2,cores = 10)
BEAM_res5 <- BEAM_res4[,c("gene_short_name","pval","qval")]
save(BEAM_res4, BEAM_res5, file = "BEAM_2.Rdata")

plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res5,qval < 0.05)),], 
                            branch_point = 1, num_clusters = 2, 
                            cores = 1,
                            branch_labels = c("Cell fate 1", "Cell fate 2"),
                            hmcols = colorRampPalette(brewer.pal(9, "Spectral"))(50) %>% rev,
                            #hmcols = paletteer_c("scico::berlin", n = 50),
                            use_gene_short_name = T, show_rownames = T)

tmp <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res3,qval<0.05)),], 
                                   branch_point = 2, num_clusters = 2, 
                                   cores = 1,
                                   branch_labels = c("Cell fate 1", "Cell fate 2"),
                                   hmcols = colorRampPalette(brewer.pal(9, "Spectral"))(50),
                                   #hmcols = paletteer_c("scico::berlin", n = 50),
                                   use_gene_short_name = T,
                                   return_heatmap = T)
tmp$annotation_row

# gammaT细胞聚类 --------------------------------------------------------------

library(harmony)
gammaT <- SCTransform(gammaT, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
# gammaT <- FindVariableFeatures(gammaT, selection.method = "vst", nfeatures = 2000)
# sele <- grep("(^MT-)|(^RP[SL]))", VariableFeatures(gammaT), value = TRUE, invert = TRUE)
sele <- unique(unlist(gamma_marker$names))
gammaT <- RunPCA(gammaT, features = sele)
# gammaT <- RunHarmony(object = gammaT, group.by.vars = "sample", plot_convergence = TRUE)

gammaT <- FindNeighbors(gammaT, dims = 1:10) # , reduction = "harmony"
gammaT <- gammaT %>% RunTSNE(dims = 1:10) %>% RunUMAP(dims = 1:10)
gammaT$State <- paste0("S_", cds$State)
gammaT$time <- cds$Pseudotime
saveRDS(gammaT, file = "gammaT.rds")

DimPlot(gammaT, reduction = "umap", group.by = "State")
DimPlot(gammaT, reduction = "tsne", group.by = "State")

dat <- data.frame(gammaT@reductions$umap@cell.embeddings, time = gammaT$time)
Idents(gammaT) <- paste0("S_", gammaT$State)
gamma_marker <- COSG::cosg(gammaT, n_genes_user = 15)
DotPlot(gammaT, features = unique(unlist(gamma_marker$names)), assay = "RNA") + RotatedAxis()

library(RColorBrewer)
ggplot(dat, aes(UMAP_1, UMAP_2, color = time)) + 
  geom_point() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62)) +
  theme_minimal()

p1 <- FeaturePlot(gammaT, features = c("GZMA")) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62))
p2 <- FeaturePlot(gammaT, features = c("CCL5")) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62))
p3 <- FeaturePlot(gammaT, features = c("CSF1")) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62))
p1|p2|p3



# WGCNA -------------------------------------------------------------------

library(COSG)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
sc <- readRDS("Tcells.rds")
sc$HPV <- case_when(sc$sample %in% c("HPV_1", "HPV_2") ~ "HPV_pos",
          sc$sample %in% c("HPV_3", "HPV_4") ~ "HPV_neg")

theme_set(theme_cowplot())

set.seed(12345)
enableWGCNAThreads(nThreads = 16)

seurat_obj <- sc
options(repr.plot.height = 6, repr.plot.width = 8)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "tutorial" 
)

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj)
Idents(seurat_obj) <- seurat_obj$cell_type
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type"),
  reduction = "pca", 
  k = 25, 
  max_shared = 10, 
  ident.group = "cell_type"  
)

seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj@misc$tutorial$wgcna_metacell_obj

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("gammaT"), # 挑选感兴趣的细胞类型
  group.by = "cell_type", 
  assay = "RNA",
  slot = "data"
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed" 
)

options(repr.plot.height = 10, repr.plot.width = 16)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power = 4, 
  setDatExpr = FALSE,
  tom_name = 'gamma' 
)
options(repr.plot.height = 4, repr.plot.width = 10)

PlotDendrogram(seurat_obj, main='gammaT hdWGCNA Dendrogram')

seurat_obj <- ModuleEigengenes(
  seurat_obj
  #group.by.vars = "patient"  # 该数据不需要进行批次校正，因此注释掉这一步
)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'gammaT'  
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "gamma-M"
)

saveRDS(seurat_obj, file = "final.seurat_obj.rds") 
seurat_obj <- readRDS("final.seurat_obj.rds")


#### 可视化部分 
# 可视化每个模块中的基因的kME
options(repr.plot.height = 10, repr.plot.width = 16)
PlotKMEs(seurat_obj, ncol = 6)

modules <- GetModules(seurat_obj)
hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 30)
write.csv(hub_df, file = "top30_hub_genes.csv", row.names = F)

hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 40)
write.csv(hub_df, file = "top40_hub_genes.csv", row.names = F)

hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 60)
write.csv(hub_df, file = "top60_hub_genes.csv", row.names = F)

# 根据每个模块的topN个基因对细胞进行打分
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25, 
  method = "Seurat" 
)

# 根据hMEs进行FeaturePlot
options(repr.plot.height = 10, repr.plot.width = 16)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs', 
  order = TRUE 
)
wrap_plots(plot_list, ncol = 6)

# 模块相关性
options(repr.plot.height = 6, repr.plot.width = 6)
col2 <- colorRampPalette(c("#fde7cf","white" ,"#6aa7ef"),alpha = TRUE)

ModuleCorrelogram(seurat_obj, col = col2(10))

# 利用Seurat自带的可视化方法
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
options(repr.plot.height = 4, repr.plot.width = 10)
DotPlot(seurat_obj, features = mods, group.by = "seurat_clusters") + 
  RotatedAxis() + 
  scale_color_gradient(low = "white", high = "#6ba7ef")
FeaturePlot(seurat_obj, features = mods, cols = c("white", "white", "red"))

options(repr.plot.height = 4, repr.plot.width = 10)

dat <- FetchData(seurat_obj, vars = c("gamma-M4", "cell_type"))
dat <- dat %>% group_by("cell_type") %>% 
  mutate(`gamma-M4` = remove_outliers(`gamma-M4`))
dat <- na.omit(dat)
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x <= (qnt[1] - H)] <- NA
  y[x >= (qnt[2] + H)] <- NA
  y
}
library(ggnormalviolin)
d <- group_by(dat, cell_type) %>% 
  summarise(dist_mean = mean(`gamma-M4`, na.rm = T), 
                                            dist_sd = sd(`gamma-M4`, na.rm = T))

p <- ggplot(data = d, 
            aes(x = cell_type, mu = dist_mean,
                sigma = dist_sd,
                fill = cell_type)) +
  theme(legend.position = "none") +
  theme_minimal()

p + 
  geom_normalviolin(
    p_tail = 0.05, alpha = 0.7,
    tail_fill = "white", 
    tail_alpha = 0.8,
    color = "gray20",
    size = 0.1
  ) + scale_fill_viridis_d(option = "D") + 
  geom_boxplot(data = dat, 
               width = 0.3, notch = T, alpha = 0.5, 
               mapping = aes(x = `cell_type`, y = `gamma-M4`), 
               inherit.aes = F, outlier.shape = NA) + 
  theme(legend.position = "NULL")


ggplot(dat) +
geom_boxplot(data = dat, mapping = aes(x = `cell_type`, y = `gamma-M6`))


theme_set(theme_cowplot())
set.seed(12345)
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "gamma-M6")

top <- read.csv("top40_hub_genes.csv", header = T)
gene_list <- split(top$gene_name, top$module)
M <- lapply(gene_list, function(x){
  GO_enrich(gene = x)
})

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

M1 <- as.data.frame(M$`gamma-M6`) %>% 
  separate(GeneRatio, into = c("number", "all")) %>% 
  mutate(number = as.numeric(number)) %>% 
  arrange(-number) %>% 
  top_n(n = 5, wt = number) %>% 
  head() %>% 
  mutate(Description = fct_reorder(Description, number), 
         nc = nchar(.$Description)) %>% 
  filter(nc < 60)

ggplot(M1, aes(y = number, x = Description)) + geom_bar(stat = "identity", fill = "grey", alpha = 0.3) +
  coord_flip() +
  scale_y_continuous(breaks = integer_breaks()) +
  theme_bw() +
  theme(panel.grid =element_blank())




options(repr.plot.height = 10, repr.plot.width = 10)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, 
  n_other = 5, 
  edge_prop = 0.75, 
  mods = "all"
)

seurat_obj <- RunModuleUMAP(
  seurat_obj, 
  n_hubs = 10, 
  n_neighbors = 15, 
  min_dist = 0.1 
)

umap_df <- GetModuleUMAP(seurat_obj)
options(repr.plot.height = 6, repr.plot.width = 8)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(
    color = umap_df$color, 
    size = umap_df$kME * 2 
  ) + umap_theme()

options(repr.plot.height = 10, repr.plot.width = 10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.25,
  sample_edges = T,
  edge_prop = 0.1, 
  label_hubs = 2, 
  keep_grey_edges = F
)

cur_traits <- c("HPV")
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits, 
  group.by = 'cell_type', subset_groups = "node" 
)


# 提取相关性矩阵
mt_cor <- GetModuleTraitCorrelation(seurat_obj)

# 模块-特征相关性热图
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',  # 可选pval、fdr作为显著性分数
  label_symbol = 'stars',  # 以*号作为显著性标记，numeric则显示数值
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = '#e98095',
  mid_color = '#ffd9c4',
  low_color = '#427db8',
  plot_max = 0.3,
  combine=TRUE
)

# 模块相关性 -------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(AUCell)
modul <- read.csv("top30_hub_genes.csv", header = T)
genelist <- split(x = modul$gene_name, f = modul$module)

gammaT <- readRDS("gammaT.rds")

cells_rankings <- AUCell_buildRankings(gammaT@assays$RNA@data) 

AUC <- AUCell_calcAUC(genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- data.frame(t(getAUC(AUC)))

gammaT@meta.data <- cbind(gammaT@meta.data, aucs)
aucs <- data.frame(t(getAUC(AUC))) %>% 
  cbind(clusters = gammaT$State, .) %>% 
  tidyr::gather(value = "value", key = "group", -clusters) %>% 
  group_by(clusters, group) %>% 
  summarise(median = median(value)) %>% 
  pivot_wider(values_from = median, names_from = "group") %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("clusters") %>% t() %>% 
  as.data.frame()
aucs <- aucs[c(1,2,4,6),]

aucs <- t(scale(t(aucs)))
library(ComplexHeatmap)
Heatmap(aucs, col = c("#21114EFF", "#2F1163FF", "#3F0F72FF", "#4D117BFF", "#5B167EFF", "#681C81FF", "#FDDD9FFF", "#FCECAEFF", "#FCFDBFFF"))

if(T){
  aucs <- data.frame(t(getAUC(AUC)))
  aucs <- cbind(gammaT$State, aucs)
  colnames(aucs)[1] <- "clusters"
  aucs <- tidyr::gather(aucs, value = "value", key = "group", -clusters)
  lis <- c("S_1", "S_2", "S_3", "S_4", "S_5")
  aucs$clusters <- factor(aucs$clusters, levels = lis)
  color_ct <- c("#e2a000", "#cc79a7", "#d55e01", "#039d71", "black")
  names(color_ct) <- lis
  aucs %>% ggplot(aes(clusters,value))+
    geom_violin(aes(fill=clusters, alpha = 0.3),scale = "width")+
    geom_boxplot(width = 0.15, outlier.shape = NA,mapping = aes(alpha = 0.3))+
    facet_grid(group~.,scales = "free_y")+
    scale_fill_manual(values = color_ct)+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = NULL,color = "black",size = 14),
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank(),
      legend.position = "none",
      panel.spacing.y = unit(0, "cm"),
      strip.text.y = element_text(angle=0,size = 14,hjust = 0),
      strip.background.y = element_blank()
    )
}

# gammaT 亚群富集分析 -----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(COSG)
gammaT <- readRDS("gammaT.rds")
gammaT$State
Idents(gammaT) <- gammaT$State
marker <- cosg(gammaT)
s_1 <- marker$names$S_1
s1 <- GO_enrich(s_1)
s1 <- arrange(s1, -Count) %>% remove_rownames() %>% filter(ONTOLOGY == "BP") %>% mutate(id = "S1")
write.csv(s5, file = "s_5.enrich.csv", row.names = F)
s_2 <- marker$names$S_2
s2 <- GO_enrich(s_2)
s2 <- arrange(s2, -Count) %>% remove_rownames() %>% filter(ONTOLOGY == "BP") %>% mutate(id = "S2")

s_3 <- marker$names$S_3
s3 <- GO_enrich(s_3)
s3 <- arrange(s3, -Count) %>% remove_rownames() %>% filter(ONTOLOGY == "BP") %>% mutate(id = "S3")

s_4 <- marker$names$S_4
s4 <- GO_enrich(s_4)
s4 <- arrange(s4, -Count) %>% remove_rownames() %>% mutate(id = "S4")

s_5 <- marker$names$S_5
s5 <- GO_enrich(s_5)
s5 <- arrange(s5, -Count) %>% remove_rownames() %>% filter(ONTOLOGY == "BP") %>% mutate(id = "S5")

dat <- rbind(s1[1:5,], s2[c(9, 15, 17, 18, 19),], s3[1:5,], s4[1:5,], s5[1:5,])
dat$Description <- fct_inorder(dat$Description)
dat$Description <- paste0(dat$Description, "_", dat$id)
dat$p <- -log10(dat$pvalue)
ggplot(dat, aes(Description, Count)) + 
  geom_bar(stat = "identity", width = 0.8, mapping = aes(fill = id))+ 
  theme_minimal() + RotatedAxis()
options(
  ggplot2.discrete.colour = ggsci::scale_colour_d3,
  ggplot2.discrete.fill = ggsci::scale_fill_d3
)
ggplot(dat, aes(Description, p, group = id, color = id)) + 
  geom_point(stat = "identity", width = 0.8, mapping = aes(fill = id))+ 
  geom_line(color = "grey", cex = 0.5) +
  theme_minimal() + RotatedAxis()


