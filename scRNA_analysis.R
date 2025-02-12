#package
library(Seurat)
library(tidyverse)
library(patchwork) 
library(harmony)
library(DoubletFinder)
library(glmGamPoi)
#Find_doublet
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  doublets.percentage = 0.08  
  nExp_poi <- round(as.numeric(doublets.percentage)*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE) # 需修改PC， SCT
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data
}
#profile location
dir = '/datapool/bioinfo/wangjj/work/05ShouyiLiujie/ResultIncludeIntron/'
info <- read.table("/datapool/bioinfo/wangjj/work/05ShouyiLiujie/ResultIncludeIntron/03TotalMerge/sample14/samInfo.txt", header = T)
sample <- info$orig.ident
#sample <- list.files(dir, '^cDNA.*')
Obj.list = lapply(sample, function(folder){
	statRaw <- as.data.frame(matrix(nrow=0,ncol=0))
	statFilter <- as.data.frame(matrix(nrow=0,ncol=0))
	statDelDoublet <- as.data.frame(matrix(nrow=0,ncol=0))
	# Create seurat object
	seurat <- CreateSeuratObject(counts = Read10X(paste0(dir,folder,"/output/filter_matrix/", sep = "")), project = folder, min.cells = 3)	
	# stat raw cell number, calculate the percentage of ribo, mt
	statRaw <- rbind(statRaw, as.data.frame(table(seurat$orig.ident)))
	# mitochondiral and ribosomal gene percent
	seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
	ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(seurat), value = TRUE)
	seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, features = ribo.genes)
    seurat[["percent.hbb"]] <- PercentageFeatureSet(seurat, pattern = "^HBA|^HBB")
	VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "ercent.pribo", "percent.hbb"), ncol = 3, pt.size = 0) 
	ggplot2::ggsave(filename = paste(folder,".raw.QC.png", sep = "") ,width = 10, height = 10)
	# filter low quality cell and stat the cell number after filtering
	seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 & percent.hbb < 5) 
	VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.ribo", "percent.hbb"), ncol = 3, pt.size = 0) 
	ggplot2::ggsave(filename = paste(folder, ".filtered.QC.png", sep = ""),width = 10, height = 10)
#	
	statFilter <- rbind(statFilter, as.data.frame(table(seurat$orig.ident)))
	seurat <- SCTransform(seurat, method = "glmGamPoi", verbose = T, vars.to.regress = c("nCount_RNA", "percent.mt"), conserve.memory = T)  
	seurat <- RunPCA(seurat, verbose = F)
	pc.num = 1:30 
	seurat <- seurat %>% 
	  RunTSNE(dims = pc.num) %>% 
	  RunUMAP(dims = pc.num)
	# remove the doublet cell
	seurat <- Find_doublet(seurat)
	write.table(seurat@meta.data,paste(folder, ".doublets_info.txt", sep = ""),sep="\t",quote=FALSE)
	DimPlot(seurat, reduction = "umap", group.by = "doublet_info")
	ggplot2::ggsave(filename = paste(folder, ".doublets.umap.png", sep = ""),width = 10, height = 10)
	seurat <- subset(seurat,subset=doublet_info=="Singlet")
	statDelDoublet <- rbind(statDelDoublet, as.data.frame(table(seurat$orig.ident)))

	stat <- left_join(x = statRaw, y = statFilter, by = 'Var1')
	stat <- left_join(x = stat, y = statDelDoublet, by = 'Var1')
	colnames(stat) <- c("SampleId", "RawCellNum", "AfterQC", "AfterDelDoublet")
	write.table(stat,paste(folder, ".statInfo.txt", sep = ""),sep="\t",quote=FALSE, row.names = F)		
	saveRDS(seurat, paste(folder,".filtered.rds", sep = ""))
	seurat
})

	# 将seurat对象进行merge
mergeRaw <- merge(Obj.list[[1]], y = c(Obj.list[2:length(Obj.list)]), 
                  add.cell.ids = sample, project = "Lung")
	## add sample group info
metatable <- read.table("samInfo.txt",header=T)
metadata <- FetchData(mergeRaw, 'orig.ident')
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x = metadata, y = metatable, by = 'orig.ident')
rownames(metadata) <- metadata$cell_id
mergeRaw <- AddMetaData(mergeRaw, metadata = metadata)
mergeRaw$cell_id <- NULL # delete cell_id column
#
table(mergeRaw$orig.ident)
#
mergeRaw[["percent.mt"]] <- PercentageFeatureSet(mergeRaw, pattern = "^MT-")
ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(mergeRaw), value = TRUE)
mergeRaw[["percent.ribo"]] <- PercentageFeatureSet(mergeRaw, features = ribo.genes)
mergeRaw[["percent.hb"]] <- PercentageFeatureSet(mergeRaw, pattern = "^HBA|^HBB") 
# order
Idents(mergeRaw) <- "orig.ident"
my_levels <- c("cDNA43", "cDNA40", "cDNA41", "cDNA44", "cDNA46", "cDNA47") # 
Idents(mergeRaw) <- factor(Idents(mergeRaw), levels= my_levels)
VlnPlot(mergeRaw, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.ribo", "percent.hbb"), ncol = 3, pt.size = 0) 
ggsave("Lung.filtered.QC.png")

saveRDS(mergeRaw, file = "./Lung.filtered.rds")

scRNA <- mergeRaw

scRNA@active.assay

DefaultAssay(scRNA) <- "RNA"
#SCTransform	
scRNA <- SCTransform(scRNA, verbose = T, vars.to.regress = c("nCount_RNA", "percent.mt"), conserve.memory = T)
#RunPCA
scRNA<- RunPCA(scRNA, npcs = 50, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num = 1:30
#Harmony
library(harmony)
scRNA <- RunHarmony(scRNA, dims = pc.num, group.by.vars = "orig.ident", verbose = FALSE)
#RunUMAP
scRNA <- RunUMAP(scRNA,reduction = "harmony", dims = pc.num)
#FindNeighbors
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution =c (0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0))
#
c <- grep("pANN_",colnames(scRNA@meta.data))
scRNA@meta.data <- scRNA@meta.data[,-c]
#another order
ord = c("WTNS" ,"WTIL33",'WTSP',"WTIL33SP")
scRNA$orig.ident = factor(scRNA$orig.ident ,levels = ord)
table(WT_ILC2$orig.ident)

scRNA$RNA_snn_res.1 = factor(scRNA$RNA_snn_res.1,levels = c(0:25))

#命名
Idents(scRNA)<-scRNA$SCT_snn_res.1 

DimPlot(scRNA, reduction = "umap", label = T, pt.size = 1, label.size = 4, raster=FALSE)

scRNA <- RenameIdents(scRNA,
 `0` = "NK",
 `1` = "CD4 T",
 `2` = "Alveolar Macrophage",
 `3` = "CD8 T",
 `4` = "Neutrophil",
 `5` = "CD4 T",
 `6` = "MThigh T",
 `7` = "Capillary",
 `8` = "CD8 T",
 `9` = "cMonocyte",
 `10` = "AT2",
 `11` = "Basophil",
 `12` = "Fibroblast",
 `13` = "Venous",
 `14` = "AT1",
 `15` = "Neutrophil",
 `16` = "Aerocyte",
 `17` = "ncMonocyte",
 `18` = "Interstitial Macrophage",
 `19` = "Alveolar Macrophage",
 `20` = "Club Cell",
 `21` = "mDC2",
    `22` = "CD8 T",
    `23` = "AT2",
    `24` = "Artery",
    `25` = "B",
    `26` = "AT1",
    `27` = "Lymphatic Cell",
    `28` = "SMC",
    `29` = "Plasma B",
    `30` = "Mesothelial cell",
    `31` = "Pericyte",
    `32` = "Ciliated Cell",
    `33` = "Neutrophil",
    `34` = "pDC",
    `35` = "Neutrophil",
    `36` = "Neutrophil",
    `37` = "AT2/AT1"
)
#Name the file.
Epithelial@meta.data$EpiType <- Epithelial@active.ident

#DotPlot
DotPlot(scRNA, features =unique(marker),dot.scale=3.5) +
  coord_flip() + RotatedAxis()+theme(axis.text.y = element_text(size = 10))+theme(axis.text.x = element_text(size = 10))

DotPlot(scRNA, features=unique(marker), dot.scale = 3, cols = c("lightgray", "red")) + RotatedAxis() + theme(panel.grid = element_line(color = 'grey'))
  
#ratio

COPD$cellType
#The number of cells in different cell populations for each sample.
table(COPD$Groups,COPD$cellType)
#Calculate the proportion of different cell populations in each group of samples.
Cellratio <- prop.table(table(COPD$cellType,COPD$Groups), margin = 2)
#Convert to a data frame format.
Cellratio <- as.data.frame(Cellratio)
#Group by cell type.
Cellratio$cluster <- Cellratio$celltype
ids <- which(Cellratio$cluster %in% c('AT1','AT2','Ciliated Cell','Club Cell'))
Cellratio$cluster <- as.character(Cellratio$cluster)
Cellratio$cluster[ids] <- 'Epithelial'

ids <- which(Cellratio$cluster %in% c('Aerocyte','Artery','Capillary','Lymphatic Cell','Venous'))

Cellratio$cluster[ids] <- 'Endothelial'

table(Cellratio$cluster)
ids <- which(Cellratio$cluster %in% c('Fibroblast','Mesothelial cell','Pericyte','SMC'))

Cellratio$cluster[ids] <- 'Stromal'

ids <- which(Cellratio$cluster %in% c('B','Plasma B','Basophil','cMonocyte','Interstitial Macrophage','Alveolar Macrophage','mDC2','Neutrophil','NK','ncMonocyte','pDC','T'))

Cellratio$cluster[ids] <- 'Immune'

Cellratio$sample <- factor(Cellratio$sample,levels = c('Ctl','A','S'))  

ED.fig2 <- Cellratio %>% group_by(cluster) %>% arrange(cluster, -ratio) %>% ungroup() 
ED.fig2$celltype <- factor(ED.fig2$celltype, levels = unique(ED.fig2$celltype))
ED.fig2 

# fig
ggplot(ED.fig2, aes(ratio, celltype, fill = cluster)) + 
  geom_col() + 
  geom_text(aes(label = paste(round(ratio*100,2),'%'),vjust=-0.2), hjust = -0.5, size = 3)  + 
  facet_wrap(~sample, scales = "free") + 
  xlim(0,0.5)+
  theme_bw() + 
  scale_fill_manual(values = c("#9370DB","#98FB98","#F08080","#1E90FF","#0000CD")) + 
  scale_x_continuous(expand = c(0.02, 0), limits = c(0,0.5)) + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.ticks = element_line(size = 2),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(angle = -30, vjust  = 1,colour = "black"), 
        axis.text = element_text(size = 12, colour = "black")) + 
  labs(x = "Percentage(%)",y='')
#ratio1
plot_group<-ggplot(ED.fig2,aes(x=sample,fill=cluster,weight=ratio))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")
plot_group
#ratio2
plot_group<-ggplot(ED.fig2,aes(x=sample,fill=celltype,weight=ratio))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")
plot_group
#VLNplot
VlnPlot(COPD,features=c("NKG7","CD3D","MARCO","FCGR3B","CA4","VCAN","SFTPD","MS4A2","PDGFRA"),pt.size=0,ncol = 3,raster=FALSE)
VlnPlot(COPD,features=c("VWF","AGER","EDNRB","LILRB2","MAF","SCGB3A1","HLA-DRB5","IGFBP3","MS4A1"),pt.size=0,ncol = 3,raster=FALSE)
VlnPlot(COPD,features=c("TFF3","CNN1","IGHG1","WT1","TRPC6","FOXJ1","LILRB4"),pt.size=0,ncol = 3,raster=FALSE)










