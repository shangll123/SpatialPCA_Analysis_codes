---
layout: default
title: Slide-Seq Analysis
nav_order: 3
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/Slideseq
---

Below are the analysis codes for the Slide-seq cerebellum data.

##### Load package
```R
library(SpatialPCA)
```
##### Load data
Slide-seq data is available at [Broad Institute’s single-cell repository](https://singlecell.broadinstitute.org/single_cell/) with ID SCP354. We also saved the raw data that we used in our examples in RData format, which can be downloaded from [here](https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT?usp=sharing).

```R
load(paste0("slideseq.rds"))
```

##### Run SpatialPCA
```R
mem_sparse1 <- peakRAM({
start_time <- Sys.time()
slideseq = CreateSpatialPCAObject(counts=sp_count, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5, customGenelist=NULL,min.loctions = 20, min.features=20)
end_time <- Sys.time()
T_sparse1 = end_time - start_time
})
mem_sparse1
T_sparse1

# > mem_sparse1
#   Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
# 1          147.783              141.9            9963.9
# > T_sparse1
# Time difference of 2.463042 mins

mem_sparse1 <- peakRAM({
start_time <- Sys.time()
slideseq = SpatialPCA_buildKernel(slideseq, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
slideseq = SpatialPCA_EstimateLoading(slideseq,fast=TRUE,SpatialPCnum=20)
slideseq = SpatialPCA_SpatialPCs(slideseq, fast=TRUE)
end_time <- Sys.time()
T_sparse1 = end_time - start_time
})
mem_sparse1
T_sparse1

## Selected kernel type is:  gaussian  
## The bandwidth is:  0.00870514399377566  
## Calculating sparse kernel matrix
#   Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
# 1         1329.051               4539           11439.8
# > T_sparse1
# Time difference of 22.15085 mins

# save results
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = as.matrix(slideseq@SpatialPCs)
SpatialPCA_result$normalized_expr  = slideseq@normalized_expr
SpatialPCA_result$location = slideseq@location
save(SpatialPCA_result, file = "slideseq_SpatialPCA_result.RData")

```


##### Obtain clustering result
```R
clusterlabel= louvain_clustering(clusternum=8,latent_dat=slideseq@SpatialPCs,knearest=round(sqrt(dim(slideseq@SpatialPCs)[2])) )

# Visualize spatial domains detected by SpatialPCA.
cbp_spatialpca <- c("lightyellow2", "coral", "lightcyan2" ,"#66C2A5", "cornflowerblue" ,"#FFD92F" ,"#E78AC3", "skyblue1")
pdf("slideseq_SpatialPCA_cluster8.pdf",width=5,height=5)
clusterlabel = SpatialPCA_result$clusterlabel
p=plot_cluster(legend="right",location=SpatialPCA_result$location,clusterlabel,pointsize=1,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=cbp_spatialpca)
p
dev.off()

```


##### cell type deconvolution
```R
# devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
library(RCTD)
library(Matrix)
# make reference data
# try DropViz
# install.packages('DropSeq.util_2.0.tar.gz', repos=NULL)
library(DropSeq.util)
 dge.path <- "~/F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt"
 dge <- loadSparseDge(dge.path) 
cell_cluster_outcomes = readRDS("~/F_GRCm38.81.P60Cerebellum_ALT.cell_cluster_outcomes.RDS")
subcluster = readRDS("~/F_GRCm38.81.P60Cerebellum_ALT.subcluster.assign.RDS")
clusterassign = readRDS("~/F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS")
annotationBraincellAtlas = readRDS("~/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
reference_clusters = as.character(paste0(1:11))
reference_clusters_common_names = c("GranularNeuron","PurkinjeNeuron","Interneurons","Interneurons_and_Other_Nnat","Microglia_Macrophage",
  "Oligodendrocyte_Polydendrocyte","BergmannGlia","Astrocyte","Choroid_Plexus","Endothelial","Fibroblast")
datt = data.frame(reference_clusters,reference_clusters_common_names )
# meta_data$subcluster_names = as.character(meta_data$subcluster)
# meta_data$subcluster_names[which(meta_data$subcluster == "1")] = "Granule_cells"
# meta_data$subcluster_names[which(meta_data$subcluster == "1-1")] = "Granule_cells"
# meta_data$subcluster_names[which(meta_data$subcluster == "2")] = "Purkinje_Neurons"
# meta_data$subcluster_names[which(meta_data$subcluster == "3")] = "Interneurons"
# meta_data$subcluster_names[which(meta_data$subcluster == "3-1")] = "Interneurons"
# meta_data$subcluster_names[which(meta_data$subcluster == "3-2")] = "Interneurons"
# meta_data$subcluster_names[which(meta_data$subcluster == "3-3")] = "Interneurons"
# meta_data$subcluster_names[which(meta_data$subcluster == "3-4")] = "Interneurons"
# meta_data$subcluster_names[which(meta_data$subcluster == "4")] = "Interneurons_and_Other_Nnat"
# meta_data$subcluster_names[which(meta_data$subcluster == "4-1")] = "Interneurons_and_Other_Nnat"
# meta_data$subcluster_names[which(meta_data$subcluster == "4-2")] = "Interneurons_and_Other_Nnat"
# meta_data$subcluster_names[which(meta_data$subcluster == "4-3")] = "Interneurons_and_Other_Nnat"
# meta_data$subcluster_names[which(meta_data$subcluster == "5")] = "Microglia_Macrophage"
# meta_data$subcluster_names[which(meta_data$subcluster == "6")] = "Oligodendrocyte_Polydendrocyte"
# meta_data$subcluster_names[which(meta_data$subcluster == "6-1")] = "Oligodendrocyte"
# meta_data$subcluster_names[which(meta_data$subcluster == "6-2")] = "Oligodendrocyte"
# meta_data$subcluster_names[which(meta_data$subcluster == "6-3")] = "Polydendrocyte"
# meta_data$subcluster_names[which(meta_data$subcluster == "7")] = "Bergmann_Glia"
# meta_data$subcluster_names[which(meta_data$subcluster == "8")] = "Astrocytes"
# meta_data$subcluster_names[which(meta_data$subcluster == "8-1")] = "Astrocytes"
# meta_data$subcluster_names[which(meta_data$subcluster == "8-2")] = "Astrocytes"
# meta_data$subcluster_names[which(meta_data$subcluster == "9")] = "Choroid_Plexus"
# meta_data$subcluster_names[which(meta_data$subcluster == "10")] = "Endothelial_Stalk"
# meta_data$subcluster_names[which(meta_data$subcluster == "10-1")] = "Endothelial_Stalk"
# meta_data$subcluster_names[which(meta_data$subcluster == "10-2")] = "Endothelial_Stalk"
# meta_data$subcluster_names[which(meta_data$subcluster == "10-3")] = "Endothelial_Stalk"
# meta_data$subcluster_names[which(meta_data$subcluster == "11")] = "Fibroblast_Like"
# meta_data$subcluster_names[which(meta_data$subcluster == "11-1")] = "Endothelial_Tip"
# meta_data$subcluster_names[which(meta_data$subcluster == "11-2")] = "Endothelial_Tip"
# meta_data$subcluster_names[which(meta_data$subcluster == "11-3")] = "Mural"
# meta_data$subcluster_names[which(meta_data$subcluster == "11-4")] = "Mural"
# meta_data$subcluster_names[which(meta_data$subcluster == "11-5")] = "Mural"
# meta_data$subcluster_names = as.factor(meta_data$subcluster_names)
# cell_cluster_outcomes = meta_data
cell_cluster_outcomes$reason = NULL
cell_cluster_outcomes = na.omit(cell_cluster_outcomes)
raw.data = dge[,match(rownames(cell_cluster_outcomes),colnames(dge))]
cell_cluster_outcomes$nUMI = colSums(raw.data)
cell_cluster_outcomes$liger_ident_coarse = cell_cluster_outcomes$cluster
reference = Seurat::CreateSeuratObject(raw.data, meta.data = cell_cluster_outcomes)
saveRDS(reference, paste(getwd(),"Ref_cerebellum.RDS",sep="/"))
ref_celltype = cell_cluster_outcomes$cluster
names(ref_celltype) = rownames(cell_cluster_outcomes)
reference <- Reference(raw.data, cell_types=ref_celltype, nUMI=colSums(raw.data))
SCTcount = sp_count[,match(colnames(expr),colnames(sp_count))]
location = location[match(colnames(SCTcount),rownames(location)),]
colnames(location) = c("x","y")
barcodes = rownames(location)
coords = data.frame(barcodes, "x" = location[,1], "y" = location[,2])
rownames(coords) = NULL
coords[,2:3] = scale(coords[,2:3])
write.csv(location, "BeadLocationsForR.csv",quote=F)
colnames(coords)[2] = "x"
colnames(coords)[3] = "y"
SCTcount_dataframe = as.data.frame(SCTcount)
rownames(SCTcount_dataframe) = NULL
coords = tibble::column_to_rownames(coords, var = "barcodes")
coords$barcodes <- NULL
nUMI <- colSums(SCTcount)
puck = SpatialRNA(coords, SCTcount,nUMI)
library(doParallel)
myRCTD <- create.RCTD(puck, reference, max_cores = 5)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results

metadata_RCTD = myRCTD@results$results_df
meta_data$cellid = rownames(meta_data)
metadata_RCTD$cellid = rownames(meta_data_RCTD)
metadata_RCTD = merge(meta_data,metadata_RCTD, by="cellid" )
metadata_RCTD$celltype = datt$reference_clusters_common_names[match(metadata_RCTD$first_type, datt$reference_clusters)]
metadata_RCTD$celltype = as.factor(metadata_RCTD$celltype)
metadataSlideseq$clusterlabel=SpatialPCA_result$clusterlabel
# save(metadataSlideseq, file="metadataSlideseq.RData")
metadata_RCTD$clusterlabel=SpatialPCA_result$clusterlabel[match(rownames(metadata_RCTD), rownames(SpatialPCA_result$location))]

#--------- Stack_barplot SpatialPCA
method="SpatialPCA"
metadata_RCTD$celltype = metadata_RCTD$first_type
metadata_RCTD$celltype = as.character(metadata_RCTD$celltype)
metadata_RCTD$second_type=NULL
metadata_RCTD$first_class=NULL
metadata_RCTD$second_class=NULL
metadata_RCTD$min_score=NULL
metadata_RCTD$singlet_score=NULL
metadata_RCTD$conv_all=NULL
metadata_RCTD$conv_doublet=NULL
metadata_RCTD$cellid = rownames(metadata_RCTD)
metadataSlideseq$cellid = rownames(metadataSlideseq)
metadata_RCTD$celltype = as.factor(metadata_RCTD$celltype)

pdf(paste0("Stack_barplot_",method,"_add.pdf"),width=8,height=5)
percentage = matrix(0,8,16)
for(k in 1:8){
metadata_sub = metadata_RCTD[which(metadata_RCTD$clusterlabel==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}

datt=preparedata(percentage)
makefigure(datt)+ggtitle(paste0(method))
percentage = matrix(0,8,16)
for(k in 1:8){
metadata_sub = metadata_RCTD[which(metadata_RCTD$clusterlabel==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadata_RCTD)[1]*100,2)
}
datt=preparedata(percentage)
datt$cluster_vec = as.character(datt$cluster_vec)
datt$cluster_vec = factor(datt$cluster_vec, levels=c(paste0("Cluster",c(1,5,2,4,3,7,8,6))),order=T)
makefigure(datt)+ggtitle(paste0(method))+coord_flip()
dev.off()
```

##### cell type proportion in spatial domains
```R
metadataRCTDSlideseq = merge(metadata_RCTD, metadataSlideseq,by="cellid")

metadataRCTDSlideseq$clusterlabel.x=NULL
# metadataRCTDSlideseq$clusterlabel
metatable = table(metadataRCTDSlideseq$clusterlabel,metadataRCTDSlideseq$celltype)
metatable_proportion = matrix(0 ,dim(metatable)[1],dim(metatable)[2])
for(celltypeid in 1:dim(metatable)[2]){
  metatable_proportion[,celltypeid] = metatable[,celltypeid]/sum(metatable[,celltypeid])
}
colnames(metatable_proportion) = colnames(metatable)
rownames(metatable_proportion) = c("Choroid plexus","White matter","Granule middle","Granule inner","Cerebellar nucleus","Molecular layer","Granule outer","Purkinje layer")
regionname=c("Choroid plexus","White matter","Granule middle","Granule inner","Cerebellar nucleus","Molecular layer","Granule outer","Purkinje layer")
cbp_spatialpca <- c("#FFD92F", "skyblue1", "#E78AC3" ,"lightcyan2", "#66C2A5" ,"coral" ,"cornflowerblue", "lightyellow2")
pdf(paste0("Slideseq cell type proportion in spatial domains SpatialPCA.pdf"),width=6,height=6)
for(celltypeid in 1:dim(metatable)[2]){
  dat = data.frame("Regions"=regionname,"Proportion"=metatable_proportion[,celltypeid])
  dat$Regions = factor(dat$Regions, levels=c("Molecular layer","Purkinje layer","Granule outer","Granule middle","Granule inner","White matter","Cerebellar nucleus","Choroid plexus"),order=T)
  p2=ggplot(dat, aes(x=Regions, y=Proportion,fill=Proportion)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
      #scale_fill_brewer(palette="Paired")+
      #scale_fill_manual(values = cbp_spatialpca)+
      scale_fill_distiller(palette = "Greens")+
      labs(title=paste0("SpatialPCA ",colnames(metatable)[celltypeid]),x="Clusters", y = "Proprotions")+
      theme(legend.position="right") +
      theme_classic()+
      #scale_fill_manual(values = method_color)+
      theme(axis.text.x = element_text(angle = 60,  hjust=1))+
      theme(plot.title = element_text(size = 22),
                  text = element_text(size = 22),
                  #axis.title = element_text(face="bold"),
                  #axis.text.x=element_text(size = 20) ,
                  legend.position = "right")

  print(p2)
  }
dev.off()
```

##### Marker gene mean expression in each spatial domain (line plot)
```R
cbp_spatialpca=c("#66C2A5", "lightyellow2", "cornflowerblue" ,"#E78AC3", "skyblue1" ,"#FFD92F" ,"lightcyan2", "coral")
countmat = sp_count[,match(colnames(SpatialPCA_result$normalized_expr), colnames(sp_count))]
dat_in = dat_in[which(dat_in$genename %in% c("Olig2","Cbln3","Eomes","Otx2","Calb1","Ntn1")),]

pdf("Markergenes_lineplot.pdf",width=3,height=10)
for(k in 1:dim(dat_in)[1]){
#for(k in 1:2){
	print(k)
	counts = countmat[which(rownames(countmat) %in% paste0(dat_in$genename[k])),match(rownames(SpatialPCA_result$location), colnames(countmat))]
	cluster = as.character(metadataSlideseq$clusterlabel)
	data=data.frame(counts, cluster)
	dat = data
	dat=dat[order(dat$counts),]
	dat$cellid = factor(paste0(1:dim(dat)[1]),levels=c(paste0(1:dim(dat)[1])),order=T)
	dat$cluster = factor(dat$cluster,levels=c(paste0(c(1,5,2,4,3,7,8,6))),order=T)
	ca <- dat %>%group_by(cluster) %>%summarise(
    	mean = mean(counts),
    	sd = sd(counts),
    	n = n(),
    	se = sd / sqrt(n))
dattt =as.data.frame(ca)
p<- ggplot(dattt, aes(x=cluster, y=mean, color=cluster,group=cluster)) + 
  #geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = cbp_spatialpca)+
  geom_vline(data=dattt, aes(xintercept=cluster, color=cluster),
             linetype="dashed",alpha=0.5)+
   geom_rect(
    aes(xmin = 0.5, xmax = 1.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 1.5, xmax = 2.5, fill = cbp_spatialpca[3],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 2.5, xmax = 3.5, fill = cbp_spatialpca[8],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 3.5, xmax = 4.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 4.5, xmax = 5.5, fill = cbp_spatialpca[7],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 5.5, xmax = 6.5, fill = cbp_spatialpca[4],alpha = 0.05), colour ="white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 6.5, xmax = 7.5, fill = cbp_spatialpca[5],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 7.5, xmax = 8.5, fill = cbp_spatialpca[6],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   coord_flip()+
   theme_bw(base_size = 22)+
  theme(plot.title = element_text(size = 22),
              legend.position = "none")+
   geom_line(aes(group=1),color="black", size=1) + ### !!!!!!!!!! aes(group=1) is important
  geom_point( size=3, color="#20854E")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
                 position=position_dodge(.9),color="#20854E") +
  labs(title=paste0(dat_in$genename[k]),x="", y = "Expression")
print(p)
}
dev.off()
```


##### GSEA analysis of region specific genes
```R
Seu <- CreateSeuratObject(counts = sp_count, project = "slideseq", min.cells = 20, min.features = 20)
Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)

Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
Seu <- FindNeighbors(Seu, dims = 1:10)
Seu <- FindClusters(Seu, resolution = 0.5)

Idents(Seu) = paste0("cluster",SpatialPCA_result$clusterlabel)
DE_gene = list()
DEgene_spatialPCA=c()
for(cluster in 1:8){
	print(cluster)
	DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = paste0("cluster",cluster), ident.2 =NULL, test.use = "MAST")
	DEgene_spatialPCA = c(DEgene_spatialPCA, rownames(DE_gene[[cluster]]))
}

> unique(DEgene_spatialPCA)
 [1] "Cst3"   "Mt1"    "Pcp2"   "Fth1"   "Mbp"    "Plp1"   "Trf"    "Apod"  
 [9] "Snap25" "Calm2"  "Pcp4"   "Dbi"    "Aldoc"  "Nefl"   "Malat1" "S100b" 
[17] "Gng13"  "Ptgds"  "Calb1"  "Car8"   "Ywhah"  "Pvalb"  "Nsg1"  


# Use gprofiler to do GSEA
library(MAST)
library(fdrtool)
library(qvalue)
library(gprofiler2)
library(ggplot2)
library(tidyverse)
library(forcats)
strsplit_func = function(input){
  strsplit(input, split = "_")[[1]][1]
}

upload_GMT_file(gmtfile = "~/h.c255all.v7.4.symbols.filter.gmt")

# make figures
datt = list()
for(num in 1:8){
  print(num)
topGOnum = 10
bg = rownames(sp_count)
gostres <- gost(query = rownames(DE_gene[[num]]), 
                organism = "gp__odZm_mQw4_yzA", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "", 
                sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                as_short_link = FALSE)
gostres$result = gostres$result[order(gostres$result$p_value),]
gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
datt[[num]] = gostres$result[1:10,]
datt[[num]]$log10p = -log10(datt[[num]]$p_value)
}



spatial_domain = c("Choroid plexus","White matter","Granule middle sublayer",
	"Granule inner sublayer","Cerebellar nucleus",
	"Molecular layer",
  "Granule outer sublayer","Purkinje layer")


p=list()
for(num in 1:8){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0(spatial_domain[num]),x="Biological terms", y = "-log10(p value)")+
  coord_flip()+
  theme_classic()+
  #geom_text(vjust = 0.5) +
  geom_text(vjust = 1, nudge_y = 0.5)+
  #ylim(0,1)+
    theme(plot.title = element_text(size = 30,color="black",face="bold"),
              text = element_text(size = 30,color="black",face="bold"),
              #axis.title = element_text(size = 25,color="black",face="bold"),
              axis.text.x=element_text(size = 30,color="black",face="bold") ,
              legend.position = "right")# +
}

pdf(paste0("GSEA_region_specific_Slideseq.pdf"),width=20,height=5)
for(num in 1:8){
  print(p[[num]])
}
dev.off()

```

##### RGB plots
```R
Z_SpatialPCA=as.matrix(slideseq@SpatialPCs)
Z_pca=get_PCA(SpatialPCA_result$normalized_expr,20)
count_use=sp_count[na.omit(match(rownames(SpatialPCA_result$normalized_expr), rownames(sp_count))),na.omit(match(colnames(SpatialPCA_result$normalized_expr), colnames(sp_count)))]
Z_NMF=get_NMF(as.matrix(count_use),20)

p1 = plot_RGB_tSNE(slideseq@location,Z_SpatialPCA,"SpatialPCA",pointsize=0.8,textsize=15)
p2 = plot_RGB_tSNE(slideseq@location,Z_pca,"PCA",pointsize=0.8,textsize=15)
p3 = plot_RGB_tSNE(slideseq@location,Z_NMF,"NMF",pointsize=0.8,textsize=15)
pdf("slideseq_RGB_tSNE.pdf",width=15, height=5)
ggarrange(p1, p2, p3, 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

p1 = plot_RGB_UMAP(slideseq@location,Z_SpatialPCA,"SpatialPCA",pointsize=0.8,textsize=15)
p2 = plot_RGB_UMAP(slideseq@location,Z_pca,"PCA",pointsize=0.8,textsize=15)
p3 = plot_RGB_UMAP(slideseq@location,Z_NMF,"NMF",pointsize=0.8,textsize=15)
pdf("slideseq_RGB_UMAP.pdf",width=15, height=5)
ggarrange(p1, p2, p3, 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

```

##### RGB variation between spatial domains
```R
 fx_dist = function(i,location,knearest=6){
                        line_i = rep(0,dim(location)[1])
                        line_i[i] = 1
                        ED=pdist(location[i,],location[-i,])@dist
                        #line_i[-i] = 1 - ED2/(ED2 + as.numeric(bandwidth))
                        ind_i=order(ED)[1:knearest]
                        return(list("ind_i"=ind_i))
}


get_RGB_var = function(p1,p2,p3){


rgb_r = p1$RGB$r
rgb_g = p1$RGB$g
rgb_b = p1$RGB$b
p1$RGB$combine = ( rgb_r*var(rgb_r)+rgb_g*var(rgb_g)+rgb_b*var(rgb_b) )/(var(rgb_r)+var(rgb_g)+var(rgb_b))
p1$RGB$combine = (p1$RGB$combine-mean(p1$RGB$combine))/sd(p1$RGB$combine)
rgb_r = p2$RGB$r
rgb_g = p2$RGB$g
rgb_b = p2$RGB$b
p2$RGB$combine = ( rgb_r*var(rgb_r)+rgb_g*var(rgb_g)+rgb_b*var(rgb_b) )/(var(rgb_r)+var(rgb_g)+var(rgb_b))
p2$RGB$combine = (p2$RGB$combine-mean(p2$RGB$combine))/sd(p2$RGB$combine)
rgb_r = p3$RGB$r
rgb_g = p3$RGB$g
rgb_b = p3$RGB$b
p3$RGB$combine = ( rgb_r*var(rgb_r)+rgb_g*var(rgb_g)+rgb_b*var(rgb_b) )/(var(rgb_r)+var(rgb_g)+var(rgb_b))
p3$RGB$combine = (p3$RGB$combine-mean(p3$RGB$combine))/sd(p3$RGB$combine)


        # SpatialPCA
        SpatialPCA_RGB_var = c()
        for(cell in 1:dim(p1$RGB)[1]){
          nearest = fx_dist(cell,p1$RGB[,1:2])$ind_i
          SpatialPCA_RGB_var[cell] = var(p1$RGB[c(nearest,cell),]$combine)
        }

        # PCA
        PCA_RGB_var = c()
        for(cell in 1:dim(p2$RGB)[1]){
          nearest = fx_dist(cell,p2$RGB[,1:2])$ind_i
          PCA_RGB_var[cell] = var(p2$RGB[c(nearest,cell),]$combine)
        }

        # NMF
        NMF_RGB_var = c()
        for(cell in 1:dim(p3$RGB)[1]){
          nearest = fx_dist(cell,p3$RGB[,1:2])$ind_i
          NMF_RGB_var[cell] = var(p3$RGB[c(nearest,cell),]$combine)
        }

  
        var_sum = c(SpatialPCA_RGB_var,PCA_RGB_var,NMF_RGB_var)
        method = c(rep("SpatialPCA",length(SpatialPCA_RGB_var)),rep("PCA",length(PCA_RGB_var)),rep("NMF",length(NMF_RGB_var)))
        dat = data.frame(method,var_sum)
        return(dat)
}


method_color=c("#1F77B4", "#FF7F0E", "#2CA02C" ,"#D62728", "#9467BD" ,"#8C564B", "#E377C2","#7F7F7F", "#BCBD22", "#17BECF")


dat_tSNE = get_RGB_var(p1,p2,p3)
dat_tSNE$method = factor(dat_tSNE$method,levels=c("SpatialPCA","PCA","NMF"),order=T)
pdf("~/slideseq_tSNE_RGB_varsum.pdf",width=7,height=7)
p=ggplot(dat_tSNE, aes(x=method, y=var_sum)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6,fill = method_color[1:3])+
  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  geom_hline(yintercept=median(dat_tSNE$var_sum[which(dat_tSNE$method=="SpatialPCA")]), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=25)+
  #ylim(0,0.5)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("Slide-seq: RGB in tSNE"),
       x="", y = "Variance in RGB")
 p
dev.off()



dat_umap = get_RGB_var(p1,p2,p3)
dat_umap$method = factor(dat_umap$method,levels=c("SpatialPCA","PCA","NMF"),order=T)
pdf("~/slideseq_UMAP_RGB_varsum.pdf",width=7,height=7)
p=ggplot(dat_umap, aes(x=method, y=var_sum)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6,fill = method_color[1:3])+
  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  geom_hline(yintercept=median(dat_umap$var_sum[which(dat_umap$method=="SpatialPCA")]), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=25)+
  #ylim(0,1)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("Slide-seq: RGB in UMAP"),
       x="", y = "Variance in RGB")
 p
dev.off()


```


