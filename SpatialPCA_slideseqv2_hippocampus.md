---
layout: default
title: Slide-Seq V2 Analysis
nav_order: 4
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/SlideseqV2
---

Below are the analysis codes for the Slide-seq V2 hippocampus data.

##### Load package
```R
library(SpatialPCA)
```
##### Load data
Slide-seq V2 hippocampus data is available at [Broad Institute’s single-cell repository](https://singlecell.broadinstitute.org/single_cell/) with ID SCP815. We also saved the raw data that we used in our examples in RData format, which can be downloaded from [here](https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT?usp=sharing).
```R
library(SpatialPCA)
library(peakRAM)

#Puck_200115_08 (mouse hippocampus), process downloaded raw data from Broad Institute’s single-cell repository:
dataset="Puck_200115_08"
Puck_200115_08_expr = read.table("Puck_200115_08.digital_expression.txt",header=T)
Puck_200115_08_loc = read.csv("Puck_200115_08_bead_locations.csv",sep=",")
countmat = Puck_200115_08_expr[,-1]
countmat = as.matrix(countmat)
rownames(countmat) = as.character(Puck_200115_08_expr[,1])
location = Puck_200115_08_loc[,2:3]
rownames(location) = colnames(countmat)

# > dim(countmat)
# [1] 23264 53208

location = as.matrix(location)
library(Matrix)
countmat=Matrix(countmat, sparse=TRUE)

save(countmat, location, file = "Puck_200115_08_count_location.RData") # this data can be downloaded from above google drive link.

```

##### Run SpatialPCA
```R
slideseqv2 = CreateSpatialPCAObject(counts=countmat, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=10, gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)

dim(slideseqv2@normalized_expr)
# > dim(slideseqv2@normalized_expr)
# [1]  3000 51398

mem_sparse1 <- peakRAM({
start_time <- Sys.time()
	slideseqv2 = SpatialPCA_buildKernel(slideseqv2, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
end_time <- Sys.time()
T_sparse1 = end_time - start_time
	slideseqv2 = SpatialPCA_EstimateLoading(slideseqv2,fast=TRUE,SpatialPCnum=20)
end_time <- Sys.time()
T_sparse2 = end_time - start_time
	slideseqv2 = SpatialPCA_SpatialPCs(slideseqv2, fast=TRUE)
end_time <- Sys.time()
T_sparse3 = end_time - start_time
})

# mem_sparse1
#   Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
# 1         22211.96            28787.2           71993.6
# # [1] 70.30625
# T_sparse1
# T_sparse2
# T_sparse3

# > T_sparse1
# Time difference of 2.821149 mins
# > T_sparse2
# Time difference of 44.05035 mins
# > T_sparse3
# Time difference of 6.16999 hours


SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = as.matrix(slideseqv2@SpatialPCs)
SpatialPCA_result$normalized_expr  = slideseqv2@normalized_expr
SpatialPCA_result$location = slideseqv2@location
save(SpatialPCA_result, file = "SlideseqV2_SpatialPCA_result.RData")


```


##### Obtain clustering result
```R
SpatialPCA_louvain = louvain_clustering(14,latent_dat=as.matrix(slideseqv2@SpatialPCs), 310) 
# spatial domain cluster label for each location

cbp_spatialpca = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
 "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
 "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
 "#6B6ECF", "#7B4173" )

  pdf(paste0("Figure_slideseqV2_SpatialPCA_cluster.pdf"))
      loc1 = unlist(slideseqv2@location[,1])
      loc2 = unlist(slideseqv2@location[,2])
      cluster = as.character(SpatialPCA_louvain)
      datt = data.frame(cluster, loc1, loc2)
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
        geom_point( alpha = 1,size=0.5) +
        scale_color_manual(values = cbp_spatialpca)+
        theme_void()+
        theme(plot.title = element_text(size = 20,  face = "bold"),
              text = element_text(size = 20),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 15) ,
              legend.position = "bottom")
    print(p)
  dev.off() 

metadataSlideseqV2 = data.frame(SpatialPCA_louvain)

```


##### Region specific gene detection
```R
Seu <- CreateSeuratObject(counts = countmat, project = "xx", min.cells = 20, min.features = 20)
Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
Seu <- FindNeighbors(Seu, dims = 1:10)
Seu <- FindClusters(Seu, resolution = 0.5)

Idents(Seu) = paste0("cluster",metadataSlideseqV2$SpatialPCA_louvain)
DE_gene = list()
DEgene_spatialPCA=c()
each_num = c()
for(cluster in 1:14){
  print(cluster)

  DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = paste0("cluster",cluster), ident.2 =NULL, test.use = "MAST")
  each_num[cluster] = dim(DE_gene[[cluster]])[1]
  DEgene_spatialPCA = c(DEgene_spatialPCA, rownames(DE_gene[[cluster]]))
}
DEgene_spatialPCA=unique(DEgene_spatialPCA)
length(DEgene_spatialPCA)
each_num

> length(DEgene_spatialPCA)
[1] 210
> each_num
 [1]  72  61  34 114  20  66  39  24  32  38  32  16  14  17

```

##### Marker gene mean expression in each spatial domain (line plot)
```R

make_lineplot = function(genename,clusterlabel,cluster_order){

  counts = countmat[which(rownames(countmat) %in% paste0(genename)),match(rownames(SpatialPCA_result$location), colnames(countmat))]
  cluster = as.character(clusterlabel)
  data=data.frame(counts, cluster)
  dat = data
  dat=dat[order(dat$counts),]
  dat$cellid = factor(paste0(1:dim(dat)[1]),levels=c(paste0(1:dim(dat)[1])),order=T)
  dat$cluster = factor(dat$cluster,levels=c(paste0(clusterorder)),order=T)
  ca <- dat %>%group_by(cluster) %>%summarise(
      mean = mean(counts),
      sd = sd(counts),
      n = n(),
      se = sd / sqrt(n))
  dattt =as.data.frame(ca)
    cbp_spatialpca = c("#FFDC91","#DADAEB")
  p<- ggplot(dattt, aes(x=cluster, y=mean, color=cluster,group=cluster)) + 
  #geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = cbp_spatialpca)+
  geom_vline(data=dattt, aes(xintercept=cluster, color=cluster),
           linetype="dashed",alpha=0.5)+
   geom_rect(
    aes(xmin = 0.5, xmax = 1.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 1.5, xmax = 2.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 2.5, xmax = 3.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 3.5, xmax = 4.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 4.5, xmax = 5.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 5.5, xmax = 6.5, fill = cbp_spatialpca[2],alpha = 0.05), colour ="white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 6.5, xmax = 7.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 7.5, xmax = 8.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 8.5, xmax = 9.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 9.5, xmax = 10.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 10.5, xmax = 11.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 11.5, xmax = 12.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 12.5, xmax = 13.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
   geom_rect(
    aes(xmin = 13.5, xmax = 14.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
    ymin = -Inf, ymax = Inf, alpha = 0.05) +
  coord_flip()+
  theme_bw(base_size = 22)+
  theme(plot.title = element_text(size = 22),
            legend.position = "none")+
  geom_line(aes(group=1),color="black", size=1) + ### !!!!!!!!!! aes(group=1) is important
  geom_point( size=3, color="#20854E")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
               position=position_dodge(.9),color="#20854E") +
  labs(title=paste0(genename),x="", y = "Expression")
  return(p)
}

clusterorder=c(9,13,5,6,1,4,2,11,7,14,8,12,10,3)
pdf("Markergenes_lineplot_SpatialPCA_slideseqV2.pdf",width=3,height=10)
make_lineplot("Itpka",metadataSlideseqV2$SpatialPCA_louvain,clusterorder)
make_lineplot("Igfbp4",metadataSlideseqV2$SpatialPCA_louvain,clusterorder)
make_lineplot("Map4",metadataSlideseqV2$SpatialPCA_louvain,clusterorder)
dev.off()


```


##### GSEA analysis of region specific genes
```R
library(MAST)
library(fdrtool)
library(qvalue)

library(gprofiler2)
library(ggplot2)
library(tidyverse)

library(gprofiler2)
library(forcats)
strsplit_func = function(input){
  strsplit(input, split = "_")[[1]][1]
}

upload_GMT_file(gmtfile = "h.c255all.v7.4.symbols.filter.gmt")

# make figures
datt = list()
for(num in 1:14){
  print(num)
topGOnum = 10
bg = rownames(countmat)
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


p=list()
for(num in 1:length(datt)){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0("cluster",num),x="Biological terms", y = "-log10(p value)")+
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

pdf(paste0("GSEA_SpatialPCA_region_specific_gene_slideseqV2.pdf"),width=20,height=5)
for(num in 1:length(datt)){
  print(p[[num]])
}
dev.off()


```

##### cell type deconvolution
```R
# following https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html

library(RCTD)
library(Matrix)

# make reference single cell data
library(DropSeq.util)
# download reference data from http://dropviz.org
#dge.path <- "/net/mulan/disk2/shanglu/Projects/SpatialPCA/sparse_kernel/Puck_200115_08/RCTD/F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz"
# dge <- loadSparseDge(dge.path) 
refdir = "~/RCTD"
dge.path <- file.path(refdir,"F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz")
class_path <- file.path(refdir,"F_GRCm38.81.P60Hippocampus.subcluster.assign.RDS")
dge <- loadSparseDge(dge.path)
genes = read.table("gene_name_Puck_200115_08.txt")
genesname = as.vector(t(genes))
genesname[which(genesname=="")]=NA
genesname = na.omit(genesname)
cells = read.table("cell_name_Puck_200115_08.txt")
cellsname = as.vector(t(cells))
cellsname[which(cellsname=="")]=NA
cellsname = na.omit(cellsname)
rownames(dge) <- genesname
colnames(dge) <- cellsname
cluster <- readRDS(class_path)
common_barcodes = intersect(names(cluster), colnames(dge))
raw.data = dge[,common_barcodes]
cluster = cluster[common_barcodes]
meta_data = as.data.frame(cluster)
meta_data$nUMI = colSums(raw.data)
meta_data$barcode = rownames(meta_data)
cell_types <- meta_data$cluster
names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list
reference <- Reference(raw.data, cell_types, nUMI)
print(dim(reference@counts))
table(reference@cell_types)

# > print(dim(reference@counts))
# [1] 27953 53204
# > table(reference@cell_types)

#   1-1  1-10  1-11  1-12  1-13  1-14  1-15  1-16  1-17  1-18  1-19   1-2  1-20 
#   108    83   134    88    85   304   171    63   179   101    42    40    84 
#  1-21  1-22  1-23  1-24  1-25  1-26  1-27   1-3   1-4   1-5   1-6   1-7   1-8 
#    73    86   153   148    52   101    49    31   160   101   244   111   105 
#   1-9  10-1  10-2  11-1  13-1  13-2  13-3  13-4  13-5  13-6  15-1  15-2  16-1 
#   139   467    78   337    75   110    41    21   130    81  1816    87   332 
#  16-2  16-3  16-4  17-1  17-2  17-3  17-4   2-1   2-2   2-3   2-4   2-5   2-6 
#   222   118    33   148    19    69    60   196   187   150   105   854   864 
#   2-7   3-1  3-10  3-11   3-2   3-3   3-4   3-5   3-6   3-7   3-8   3-9   4-1 
#   731   107   398    57   269   240    15    27   549   379   131    51 11727 
#   4-2   5-1  5-10  5-11  5-12  5-13  5-14  5-15   5-2   5-3   5-4   5-5   5-6 
#  1538  3088    56   787    80    43   689   327   500   621   588    74  1556 
#   5-7   5-8   5-9   6-1   6-2   6-3   6-4   6-5   6-6   6-7   6-8   7-1   7-2 
#   406   144   156  5063   243   649   330   114   416   246    71  6767   491 
#   7-3   8-1   8-2   8-3   8-4   8-5   9-1   9-2   9-3   9-4  12-1  14-1 
#   245   556   227  1100   292    21   208   155   387   296    22   336 


# spatial data
nUMI <- colSums(countmat) # In this case, total counts per pixel is nUMI
### Create SpatialRNA object
location = as.data.frame(location)

# split entire data into four parts, run RCTD on each part. otherwise the sample size is too large to run RCTD.
puck1of4 <- SpatialRNA(location[1:10000,], countmat[,1:10000], nUMI[1:10000])
myRCTD1of4 <- create.RCTD(puck1of4, reference, CELL_MIN_INSTANCE = 15,max_cores = 30)
myRCTD1of4 <- run.RCTD(myRCTD1of4, doublet_mode = 'doublet')
result1of4 <- myRCTD1of4@results

puck2of4 <- SpatialRNA(location[10001:20000,], countmat[,10001:20000], nUMI[10001:20000])
myRCTD2of4 <- create.RCTD(puck2of4, reference, CELL_MIN_INSTANCE = 15,max_cores = 30)
myRCTD2of4 <- run.RCTD(myRCTD2of4, doublet_mode = 'doublet')
result2of4 <- myRCTD2of4@results


puck3of4 <- SpatialRNA(location[20001:30000,], countmat[,20001:30000], nUMI[20001:30000])
myRCTD3of4 <- create.RCTD(puck3of4, reference, CELL_MIN_INSTANCE = 15,max_cores = 30)
myRCTD3of4 <- run.RCTD(myRCTD3of4, doublet_mode = 'doublet')
result3of4 <- myRCTD3of4@results


puck4of4 <- SpatialRNA(location[30001:40000,], countmat[,30001:40000], nUMI[30001:40000])
myRCTD4of4 <- create.RCTD(puck4of4, reference, CELL_MIN_INSTANCE = 15,max_cores = 30)
myRCTD4of4 <- run.RCTD(myRCTD4of4, doublet_mode = 'doublet')
result4of4 <- myRCTD4of4@results


puck5of4 <- SpatialRNA(location[40001:dim(location)[1],], countmat[,40001:dim(location)[1]], nUMI[40001:dim(location)[1]])
myRCTD5of4 <- create.RCTD(puck5of4, reference, CELL_MIN_INSTANCE = 15,max_cores = 30)
myRCTD5of4 <- run.RCTD(myRCTD5of4, doublet_mode = 'doublet')
result5of4 <- myRCTD5of4@results



Puck_200115_08_RCTD_result = rbind(result1of4$results_df,result2of4$results_df,result3of4$results_df,result4of4$results_df,result5of4$results_df)



metaRCTD = data.frame("CellID"=rownames(Puck_200115_08_RCTD_result), 
						first_type=Puck_200115_08_RCTD_result$first_type, 
						second_type=Puck_200115_08_RCTD_result$second_type)


RCTDcell_iddex = match(rownames(Puck_200115_08_RCTD_result), colnames(SpatialPCA_result$normalized_expr))
#> length(RCTDcell_iddex)
#[1] 31798

metaRCTD$x_loc = SpatialPCA_result$location[RCTDcell_iddex,1]
metaRCTD$y_loc = SpatialPCA_result$location[RCTDcell_iddex,2]

celltypeanno = read.csv("../RCTD/Dropviz_Hippocampus.csv")
celltype_index = match(metaRCTD$first_type, celltypeanno$subcluster)
metaRCTD$celltype = celltypeanno$common_name.1[celltype_index]
metaRCTD$celltype = as.factor(metaRCTD$celltype)



metadataSlideseqV2$CellID = rownames(metadataSlideseqV2)

metaRCTD = merge(metadataSlideseqV2,metaRCTD, by="CellID")
save(metaRCTD, file = "metaRCTD_slideseqv2.RData")
```

##### Cell type proportion in each spatial domain
```R

colors = c("#F39B7F","#91D1C2",
  "#E64B35","#BC3C29",
  "#E18727","#EFC000",
  "#9C9EDE","#8F7700",
  "#C7E9C0","#FF9896",
  "#E377C2","#71D0F5",
  "#709AE1","#17BECF",
  "#0099B4","#79AF97",
  "#31A354","#FED439",
  "#00A087","#74C476",
  "#C6DBEF","#C49C94",
  "#D9D9D9","#8CA252",
  "#FFDC91","#6B6ECF",
  "#9632B8","#DBDB8D",
  "#E7CB94","#9ECAE1")


preparedata = function(percentage){
  rownames(percentage) = c("CA1","Dentate gyrus","Third ventricle","CA3","Layer6","Corpus callosum",
  "Hippocampus(slm)","Thalamus subregion1","Layer4","Thalamus subregion3","Hippocampus(so/sr)","Thalamus subregion2","Layer5","Hippocampus(so)")
  celltype = names(table(metaRCTD$celltype))
  colnames(percentage) = celltype
  rownames(percentage) = c("CA1","Dentate gyrus","Third ventricle","CA3","Layer6","Corpus callosum",
  "Hippocampus(slm)","Thalamus subregion1","Layer4","Thalamus subregion3","Hippocampus(so/sr)","Thalamus subregion2","Layer5","Hippocampus(so)")
  percentage_vec = c(percentage)
  cluster_vec = c(rep(c("CA1","Dentate gyrus","Third ventricle","CA3","Layer6","Corpus callosum",
  "Hippocampus(slm)","Thalamus subregion1","Layer4","Thalamus subregion3","Hippocampus(so/sr)","Thalamus subregion2","Layer5","Hippocampus(so)"),length(celltype)))
  CellType = c(rep(celltype,each=14))
  datt = data.frame(cluster_vec, percentage_vec,CellType)
  datt$cluster_vec = factor(cluster_vec, level=c("Third ventricle","Thalamus subregion3","Thalamus subregion2","Thalamus subregion1","Hippocampus(so)","Hippocampus(slm)","Hippocampus(so/sr)","Dentate gyrus","CA3","CA1","Corpus callosum","Layer6","Layer5","Layer4"),order=T)
  return(datt)
}
makefigure = function(datt){
p=ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = CellType)) +        ## global aes
  scale_fill_manual(values=colors)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black",alpha=0.7) +
  theme_bw()+xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90,  hjust=1))+
  theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
return(p)
}
clusternum=14
celltypes=30
# SpatialPCA
method="SpatialPCA"
pdf(paste0("Stack_barplot_",method,"_Jan10.pdf"),width=15,height=8)
      percentage = matrix(0,clusternum,celltypes)
      for(k in 1:clusternum){
      metadata_sub = metaRCTD[which(metaRCTD$SpatialPCA_louvain==k ),]
      match_type = metadata_sub$celltype
      percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
      }
      
datt=preparedata(percentage)
makefigure(datt)+ggtitle(paste0(method))

      percentage = matrix(0,clusternum,celltypes)
      for(k in 1:clusternum){
      metadata_sub = metaRCTD[which(metaRCTD$SpatialPCA_louvain==k ),]
      match_type = metadata_sub$celltype
      percentage[k,] = round(unlist(table(match_type))/dim(metaRCTD)[1]*100,2)
      }
      
datt=preparedata(percentage)
makefigure(datt)+ggtitle(paste0(method))

dev.off()
```

##### For each cell type, visualize its proportion distribution across spatial domains
```R
celltypeanno = read.csv("~/Dropviz_Hippocampus_v2.csv")
celltype_index = match(metaRCTD$first_type, celltypeanno$subcluster)
metaRCTD$celltype = celltypeanno$use_name_merge[celltype_index]
metaRCTD$celltype = as.factor(metaRCTD$celltype)

# SpatialPCA
#metaRCTD$SpatialPCA_louvain = factor(metaRCTD$SpatialPCA_louvain,levels=c(1:14),order=T)
metatable = table(metaRCTD$SpatialPCA_louvain,metaRCTD$celltype)
metatable_proportion = matrix(0 ,dim(metatable)[1],dim(metatable)[2])
for(celltypeid in 1:dim(metatable)[2]){
  metatable_proportion[,celltypeid] = metatable[,celltypeid]/sum(metatable[,celltypeid])
}
colnames(metatable_proportion) = colnames(metatable)
rownames(metatable_proportion) = c("CA1","Thalamus subregion3","Hippocampus(so/sr)","Thalamus subregion2","Layer5","Hippocampus(so)","Dentate gyrus","Third ventricle","CA3","Layer6","Corpus callosum",
  "Hippocampus(slm)","Thalamus subregion1","Layer4")
regionname = c("CA1","Thalamus subregion3","Hippocampus(so/sr)","Thalamus subregion2","Layer5","Hippocampus(so)","Dentate gyrus","Third ventricle","CA3","Layer6","Corpus callosum",
  "Hippocampus(slm)","Thalamus subregion1","Layer4")

cbp_spatialpca = c( "#6B6ECF" ,"#DE9ED6", "#9EDAE5","#31A354","#CE6DBD" ,
  "#91D1C2", "#FF9896" ,"#7B4173" ,"yellow", "#709AE1",
  "#DADAEB", "#C7E9C0" ,"#BCBD22", "#FD7446" )

pdf(paste0("SlideseqV2 cell type proportion in spatial domains SpatialPCA.pdf"),width=8,height=6)
for(celltypeid in 1:dim(metatable)[2]){
  dat = data.frame("Regions"=regionname,"Proportion"=metatable_proportion[,celltypeid])
  dat$Regions = factor(dat$Regions, levels=c("Third ventricle","Thalamus subregion3","Thalamus subregion2","Thalamus subregion1","Hippocampus(so)","Hippocampus(slm)","Hippocampus(so/sr)","Dentate gyrus","CA3","CA1","Corpus callosum","Layer6","Layer5","Layer4"),order=T)
  p2=ggplot(dat, aes(x=Regions, y=Proportion,fill=Proportion)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
      #scale_fill_brewer(palette="Paired")+
      #scale_fill_manual(values = "#6B6ECF")+
      scale_fill_distiller(palette = "Greens")+
      #scale_colour_brewer(palette = "Reds")+
      labs(title=paste0("SpatialPCA - ",colnames(metatable)[celltypeid]),x="", y = "Proprotions")+
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


##### Trajectory analysis focusing on cortical layers
```R
# cortical layer

library(slingshot)
library(SingleCellExperiment)


cellid = which(metadataSlideseqV2$SpatialPCA_louvain %in% c("9","13","5"))
#> length(cellid)
#[1] 13195
Z_SpatialPCA = SpatialPCA_result$SpatialPCs
sim_spatialPCA_layer = SingleCellExperiment(assays = list(counts = as.matrix(Z_SpatialPCA)[,cellid]))
reducedDims(sim_spatialPCA_layer) <- SimpleList(DRM = t(as.matrix(Z_SpatialPCA)[,cellid]))
colData(sim_spatialPCA_layer)$louvain <- factor(metadataSlideseqV2$SpatialPCA_louvain[cellid])    
sim_spatialPCA_layer  <-slingshot(sim_spatialPCA_layer, clusterLabels = 'louvain', reducedDim = 'DRM',start.clus="5" )

summary(sim_spatialPCA_layer@colData@listData)
# > summary(sim_spatialPCA_layer@colData@listData)
#                   Length Class  Mode   
# louvain           13195  factor numeric
# slingPseudotime_1 13195  -none- numeric

save(sim_spatialPCA_layer, file = "slideseqV2_slingshot_sim_SpatialPCA_layer.RData")
spatialPCA_layer_cluster = metadataSlideseqV2$SpatialPCA_louvain[cellid]

```

##### Gene enrichment analysis with pseudotime associated genes
```R

library(parallel)
library(MASS)

get_pesudotime_gene = function(sim,count_use){
  genelist=list()
  for(traj in 1:(length(sim@colData@listData)-1)){
    print(paste0("traj ",traj))
    pseudotime = sim@colData@listData[[traj+1]]

    pseudotime_gene = function(ind){
      expr = count_use[ind,]
      dat = data.frame(expr,pseudotime)
      colnames(dat) = c("expr","pseudotime")
        dat = na.omit(dat)
        res = try(summary(m1 <- glm.nb(expr ~ pseudotime, data = dat)))
        if(isTRUE(class(res)=="try-error")) { next } else { 
            p_value_traj =  res$coefficients[2,4] 
            gene = rownames(count_use)[ind]
        } 
        if(p_value_traj<0.05/length(expr)){
          return(gene)
        }else {
          return(NA)
        }
    }
  
  results = unlist(mclapply(1:dim(count_use)[1], pseudotime_gene, mc.cores = 50))
  errorid = which(results=="Error in FUN(X[[i]], ...) : no loop for break/next, jumping to top level\n")
  if(length(errorid)>0){
      results=results[-errorid]
      genelist[[traj]] = na.omit(results)
    }else{
      genelist[[traj]] = na.omit(results)
    }
  }
  
  return(genelist)
}

load("slideseqV2_slingshot_sim_SpatialPCA_layer.RData")

expr_counts = countmat[,match(colnames(result$normalized_expr),colnames(countmat))]
cellid = which(metadataSlideseqV2$SpatialPCA_louvain %in% c("9","13","5"))
count_use = expr_counts[,cellid]

slideseqV2_genelist_all = c()
for(i in 1:233){
  print(i)
  start = (i-1)*100+1
  end= i*100
  tmp=get_pesudotime_gene(sim_spatialPCA_layer,count_use[start:end,])
  slideseqV2_genelist_all = c(slideseqV2_genelist_all,unlist(tmp[[1]]))
}


length(slideseqV2_genelist_all)
> length(slideseqV2_genelist_all)
[1] 883


slideseq_genelist = list(slideseqV2_genelist_all)
datt = list()
for(num in 1){
  print(num)
topGOnum = 10
bg = rownames(expr_counts)
gostres <- gost(query = slideseq_genelist[[num]], 
                organism = "gp__800L_KrJY_4lE", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "", 
                #sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                as_short_link = FALSE)
gostres$result = gostres$result[order(gostres$result$p_value),]
gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
#gostres$result$Source[which(gostres$result$Source=="MODULE")]="CANCER MODULE"
datt[[num]] = gostres$result[1:10,]
datt[[num]]$log10p = -log10(datt[[num]]$p_value)
}


p=list()
for(num in 1){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0("Pseudotime trajectory",num),x="Biological terms", y = "-log10(p value)")+
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

pdf(paste0("GSEA_SpatialPCA_msigdb_slideseqV2_pseudotime_gene.pdf"),width=20,height=5)
print(p[[1]])
dev.off()

```







