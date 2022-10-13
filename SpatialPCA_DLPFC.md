---
layout: default
title: DLPFC Analysis
nav_order: 2
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/DLPFC
---

Below are the analysis codes for the DLPFC data.

##### Load package
```R
library(SpatialPCA)
```

We used sample 151673 in DLPFC data as a main example which contains expression measurement of 33,538 genes on 3,639 spots. All 12 DLPFC samples can be downloaded from their original study [spatialLIBD](http://research.libd.org/spatialLIBD/index.html). We also saved the raw data that we used in our examples in RData format, which can be downloaded from [here](https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT?usp=sharing).


##### LIBD data
Below are the codes I used when exploring the LIBD data and save the ready to use datasets.
```R

#----------------
# LIBD
#----------------

dataset="LIBD"

# BiocManager::install("spatialLIBD")
library('spatialLIBD')
sce <- fetch_data(type = 'sce')
metaData = SingleCellExperiment::colData(sce)
expr = SingleCellExperiment::counts(sce)
sample_names <- paste0("sample_", unique(colData(sce)$sample_name))

# first try on sample 151673
i=9 # same as sample used in Bayespace, sample_151673
sce_sub <- sce[, colData(sce)$sample_name == gsub("^sample_", "", sample_names[i])]
dim(sce_sub)
count_sub = SingleCellExperiment::counts(sce_sub)
xy_coords <- data.frame(
        x_coord = colData(sce_sub)[, c("imagecol")], 
        y_coord = -colData(sce_sub)[, c("imagerow")]
    )

## Load known marker genes (from Kristen)
# load spreadsheet of known marker genes (from Kristen)
# download from https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/KRM_Layer_Markers.xlsx
library( readxl)
KRM_Layer_Markers <- read_xlsx("/net/mulan/disk2/shanglu/Projects/SpatialPCA/data/LIBD/KRM_Layer_Markers.xlsx")
KRM_Layer_Markers
dim(KRM_Layer_Markers)
marker_genes_KRM <- KRM_Layer_Markers$Gene
# get gene IDs (note: not all available in "sce" object)
sum(toupper(marker_genes_KRM) %in% rowData(sce)$gene_name)
genes_markers <- data.frame(gene_name = toupper(marker_genes_KRM))
genes_markers$gene_id <- rowData(sce)$gene_id[match(genes_markers$gene_name, rowData(sce)$gene_name)]
dim(genes_markers)
sum(is.na(genes_markers$gene_id))
# > dim(genes_markers)
# [1] 81  2
# > sum(is.na(genes_markers$gene_id))
# [1] 4

# check ground truth layer number
true_num=c()
 metaData1 = as.data.frame(metaData)
 name=unique(colData(sce)$sample_name)
for(i in 1:12){
	true_num[i] = length(na.omit(unique(metaData1$layer_guess_reordered[which(metaData1$sample_name==name[i])])))

}
> true_num
 [1] 7 7 7 7 5 5 5 5 7 7 7 7

# > dim(count_sub)
# [1] 33538  3639
# > dim(xy_coords)
# [1] 3639    2
# > count_sub[1:4,1:4]
# 4 x 4 sparse Matrix of class "dgCMatrix"
#                 AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1 AAACACCAATAACTGC-1
# ENSG00000243485                  .                  .                  .
# ENSG00000237613                  .                  .                  .
# ENSG00000186092                  .                  .                  .
# ENSG00000238009                  .                  .                  .
#                 AAACAGAGCGACTCCT-1
# ENSG00000243485                  .
# ENSG00000237613                  .
# ENSG00000186092                  .
# ENSG00000238009                  .

# then generate ready to use data for sample i from 1 to 12
# also get ground truth labels for each sample
    metaData1 = as.data.frame(metaData)
    KRM_manual_layers_sub <- filter(metaData1, sample_name == gsub("^sample_", "", sample_names[i]))
    dim(KRM_manual_layers_sub)
     ground_truth_sub <- data.frame(
         truth = rep(NA, ncol(sce_sub)))
     rownames(ground_truth_sub) <- colnames(sce_sub)
    ground_truth_sub$truth = KRM_manual_layers_sub$layer_guess_reordered
    Layer_sub = KRM_manual_layers_sub$layer_guess_reordered
    sum(as.character(KRM_manual_layers_sub$barcode) == colnames(count_sub)) # yes all matched
    save(Layer_sub, xy_coords, KRM_manual_layers_sub, count_sub, file = paste0("LIBD_sample",i,".RData") ) 

```


#### Load data

Or you can directly download the data I processed and saved in the [Google drive](https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT?usp=sharing) mentioned above.
```R
i=9 # use sample 9 as an example
path_figure="~/SpatialPCA/LIBD/SpatialPCA"
load(paste0("~/LIBD_sample",i,".RData") ) 
dataset=paste0("LIBD_sample",i,"_SpatialPCA")
clusterNum=c(7,7,7,7,5,5,5,5,7,7,7,7)

```

##### Run SpatialPCA
```R
library(SPARK)
library(Seurat)
library(peakRAM)
library(SpatialPCA)
library(ggplot2)

xy_coords = as.matrix(xy_coords)
rownames(xy_coords) = colnames(count_sub)

LIBD = CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)
mem <- peakRAM({
start_time <- Sys.time()
LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
LIBD = SpatialPCA_EstimateLoading(LIBD,fast=FALSE,SpatialPCnum=20)
LIBD = SpatialPCA_SpatialPCs(LIBD, fast=FALSE)
end_time <- Sys.time()
T = end_time - start_time
})
T

# Time difference of 4.509904 mins

# save results
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = LIBD@SpatialPCs # extracted spatial PCs
SpatialPCA_result$location = LIBD@location
truth = KRM_manual_layers_sub$layer_guess_reordered[match(colnames(LIBD@normalized_expr),colnames(count_sub))]
SpatialPCA_result$truth = truth
pred_cluster= walktrap_clustering(clusterNum[i],SpatialPCA_result$SpatialPCs,knearest=70 ) # I set knearest=70 for all 12 samples in the DLPFC data. User can try other knearest number in other datasets, e.g. knearest=round(sqrt(dim(SpatialPCA_result$SpatialPCs)[2]))
SpatialPCA_result$clusterlabel = pred_cluster
SpatialPCA_result$clusterlabel_refine=refine_cluster_10x(pred_cluster,SpatialPCA_result$location,shape="hexagon") 
# Here we followed SpaGCN to refine the clusters.
ind_na=which(is.na(SpatialPCA_result$truth)) # we remove NA samples in the original annotation to calculate ARI.
SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine[-ind_na],SpatialPCA_result$truth[-ind_na])
# for sample i, calculate the ARI.
# PAS and CHAOS scores:
SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$location)
SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$location)
# LISI score calculate:
library(lisi)
metadata_DLPFC = data.frame("SpatialPCA" = SpatialPCA_result$clusterlabel_refine)
SpatialPCA_result$LISI <- compute_lisi(SpatialPCA_result$location, metadata_DLPFC, c('SpatialPCA'))

save(SpatialPCA_result, file = paste0("SpatialPCA_DLPFC_sample",i,"result.RData"))
```

##### visualize ground truth annotations
```R
###################
ground truth
###################
path_LIBD="~/SpatialPCA/main_figures/LIBD/"
cbp=c("#5CB85C" ,"#9C9EDE" ,"#FFDC91", "#4DBBD5" ,"#FF9896" ,"#FED439", "#E377C2", "#FED439")

p=list()
for(i in 1:12){
load(paste0("SpatialPCA_DLPFC_sample",i,"result.RData"))
clusterlabel=SpatialPCA_result$truth
p[[i]] = plot_cluster(location=SpatialPCA_result$location,clusterlabel,pointsize=4,text_size=40 ,title_in=paste0("Groundtruth Sample ",sample_names[i]),color_in=cbp)
}
pdf(paste0(path_LIBD,"Sample12_groundtruth_clustering_1to6.pdf"),width=50,height=10)
print(ggarrange(p[[1]], p[[2]], p[[3]], 
 p[[4]], p[[5]], p[[6]], 
          # labels = c("A", "B", "C"),
          ncol = 6, nrow = 1))
dev.off()

pdf(paste0(path_LIBD,"Sample12_groundtruth_clustering_7to12.pdf"),width=50,height=10)
print(ggarrange(
 p[[7]], p[[8]], p[[9]], 
 p[[10]], p[[11]], p[[12]], 
          # labels = c("A", "B", "C"),
          ncol = 6, nrow = 1))
dev.off()

```

##### Obtain SpatialPCA clustering results
```R
###################
SpatialPCA results
###################
cbp=c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91")

p=list()
ARI=c()
for(i in 1:12){
load( paste0("SpatialPCA_DLPFC_sample",i,"result.RData"))
clusterlabel=SpatialPCA_result$clusterlabel_refine
ARI[i] = SpatialPCA_result$ARI
p[[i]] = plot_cluster(location=SpatialPCA_result$location,clusterlabel,pointsize=4,text_size=40 ,title_in=paste0("SpatialPCA Sample ",sample_names[i]),color_in=cbp)
}
print("ARI:")
print(ARI)
print(median(ARI))

> print(ARI)
 [1] 0.5399951 0.5174669 0.4557821 0.4678340 0.3760976 0.5222957 0.5974674
 [8] 0.5806230 0.5770166 0.5543583 0.5446831 0.6352343
> print(median(ARI))
[1] 0.5423391
```

##### tSNE and UMAP
```R
#################################################################

			tSNE & UMAP

#################################################################
i=9 # use 9th sample
library(ggpubr)
library(umap)

xy_coords = as.matrix(xy_coords)
rownames(xy_coords) = colnames(count_sub)
load( paste0("SpatialPCA_DLPFC_sample",i,"result.RData"))
Z_pca=get_PCA(SpatialPCA_result$normalized_expr,20)
count_use=count_sub[na.omit(match(rownames(SpatialPCA_result$normalized_expr), rownames(count_sub))),na.omit(match(colnames(SpatialPCA_result$normalized_expr), colnames(count_sub)))]
Z_NMF=get_NMF(as.matrix(count_use),20)

set.seed(1234)
p1 = plot_RGB_tSNE(SpatialPCA_result$location,SpatialPCA_result$SpatialPCs,pointsize=2,textsize=15)
p2 = plot_RGB_tSNE(SpatialPCA_result$location,Z_pca,pointsize=2,textsize=15)
p3 = plot_RGB_tSNE(SpatialPCA_result$location,Z_NMF,pointsize=2,textsize=15)
save(p1,p2,p3,file=paste0("tSNE_1x.RData"))
pdf("LIBD_RGB_tSNE_1x.pdf",width=12, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

set.seed(1234)
p1 = plot_RGB_UMAP(SpatialPCA_result$location,SpatialPCA_result$SpatialPCs,pointsize=2,textsize=15)
p2 = plot_RGB_UMAP(SpatialPCA_result$location,Z_pca,pointsize=2,textsize=15)
p3 = plot_RGB_UMAP(SpatialPCA_result$location,Z_NMF,pointsize=2,textsize=15)
save(p1,p2,p3,file=paste0("UMAP_1x.RData"))
pdf("LIBD_RGB_UMAP_1x.pdf",width=12, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()
```

##### RGB variance between ground truth clusters
Compare SpatialPCA, PCA, and NMF.
```R
# average variation in groundtruth clusters
# load("tSNE_1x.RData")
dat_tSNE = get_RGB_var(p1,p2,p3)
dat_tSNE$method = factor(dat_tSNE$method,levels=c("SpatialPCA","PCA","NMF"),order=T)
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


dat_tSNE$truth = rep(SpatialPCA_result$truth,3)
dat_tSNE$RGBcombine = c(p1$RGB$combine,p2$RGB$combine,p3$RGB$combine )
dat_tSNE_use = na.omit(dat_tSNE)


pdf("~/SpatialPCA/main_figures/DLPFC_tSNE_RGB_RGBcombine_between_methods.pdf",width=15,height=7)
p=ggplot(dat_tSNE_use, aes(x=truth, y=RGBcombine,fill=method)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6)+
  facet_wrap(~method)+
  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  #geom_hline(yintercept=median(dat_tSNE$var_sum[which(dat$method=="SpatialPCA")]), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=25)+
 # ylim(0,0.5)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("DLPFC (by ground truth clusters): RGB in tSNE"),
       x="", y = "RGB combined values")
 p
dev.off()


# load(paste0("UMAP_1x.RData"))
dat_umap = get_RGB_var(p1,p2,p3)
dat_umap$method = factor(dat_umap$method,levels=c("SpatialPCA","PCA","NMF"),order=T)
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


dat_umap$truth = rep(SpatialPCA_result$truth,3)
dat_umap$RGBcombine = c(p1$RGB$combine,p2$RGB$combine,p3$RGB$combine )
dat_umap_use = na.omit(dat_umap)

pdf("~/SpatialPCA/main_figures/DLPFC_UMAP_RGB_RGBcombine_between_methods.pdf",width=15,height=7)
p=ggplot(dat_umap_use, aes(x=truth, y=RGBcombine,fill=method)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6)+
  facet_wrap(~method)+
  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  #geom_hline(yintercept=median(dat_tSNE$var_sum[which(dat$method=="SpatialPCA")]), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=25)+
 # ylim(0,0.5)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("DLPFC (by ground truth clusters): RGB in UMAP"),
       x="", y = "RGB combined values")
 p
dev.off()



```


##### Trajectory analysis
```R

library(slingshot)

# SpatialPCA
i=9
sim<- SingleCellExperiment(assays = count_sub)
reducedDims(sim) <- SimpleList(DRM = t(SpatialPCA_result$SpatialPCs))
colData(sim)$Walktrap <- factor(SpatialPCA_result$clusterlabel_refine)    
sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="3" ) 
# in this data we set white matter region as start cluster, one can change to their preferred start region 
summary(sim@colData@listData)
# > summary(sim@colData@listData)
#                   Length Class              Mode   
# Walktrap          3639   factor             numeric
# slingshot         3639   PseudotimeOrdering S4     
# slingPseudotime_1 3639   -none-             numeric

sim_SpatialPCA=sim
# save(sim_SpatialPCA, file = "sim_SpatialPCA.RData")

pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
clusterlabels = SpatialPCA_result$clusterlabel_refine
gridnum = 10
color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
pdf("Trajectory_SpatialPCA.pdf",width=7,height=5)
p_traj1 = plot_trajectory(pseudotime_traj1, SpatialPCA_result$location,clusterlabels,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj1[[1]],
          ncol = 2, nrow = 1))
dev.off()	  

```

##### Detect region differential expressed genes in Seurat
```R
Seu <- CreateSeuratObject(counts = count_sub, project = "LIBD", min.cells = 20, min.features = 20)
Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
Seu <- FindNeighbors(Seu, dims = 1:10)
Seu <- FindClusters(Seu, resolution = 0.5)

Idents(Seu) = paste0("cluster",SpatialPCA_result$clusterlabel_refine)
DE_gene = list()
DEgene_spatialPCA=c()
each_num = c()
for(cluster in 1:7){
	print(cluster)
	DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = paste0("cluster",cluster), ident.2 =NULL, test.use = "MAST")
	DEgene_spatialPCA = c(DEgene_spatialPCA, rownames(DE_gene[[cluster]]))
	each_num[cluster] = dim(DE_gene[[cluster]])[1]
}
DEgene_spatialPCA=unique(DEgene_spatialPCA)

> each_num
[1]  252  265 1327  172  395   68  196

#region specific gene GSEA
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

# match ENSG with symbol
gene_dict = read.table("~/gene_annotation_downloaded_from_ENCODE.txt")
colnames(gene_dict)=gene_dict[1,]
gene_dict=gene_dict[-1,]
gene_name = as.character(gene_dict$ENSG)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
gene_name02=unlist(gene_name02)
gene_dict$ENSG = gene_name02


diff_gene_symbol = list()
for(num in 1:7){
	diff_gene_symbol[[num]] = gene_dict$gene[match(rownames(DE_gene[[num]]), gene_dict$ENSG)]
	print(length(diff_gene_symbol[[num]]))
}


upload_GMT_file(gmtfile = "~/h.c255all.v7.4.symbols.filter.gmt")

# make figures
datt = list()
for(num in 1:7){
  print(num)
topGOnum = 10
bg = as.vector(na.omit(gene_dict$gene[match(rownames(count_sub), gene_dict$ENSG)] ))
gostres <- gost(query = na.omit(diff_gene_symbol[[num]]), 
                organism = "gp__YCER_H1Gl_DW0", ordered_query = FALSE, 
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

pdf(paste0("GSEA_SpatialPCA_region_specific_gene_LIBD.pdf"),width=20,height=5)
for(num in 1:length(datt)){
  print(p[[num]])
}
dev.off()
```

##### Gene enrichment analysis with pseudotime associated genes
```R
library(parallel)
library(MASS)
length(ST_genelist)
library(gprofiler2)
library(forcats)

strsplit_func = function(input){
  strsplit(input, split = "_")[[1]][1]
}


get_pesudotime_gene = function(sim,count_use){
	genelist=list()
		pseudotime = sim@colData@listData[[3]]

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
			genelist = na.omit(results)
		}else{
			genelist = na.omit(results)
		}
	
	
	return(genelist)
}

count_use = count_sub
LIBD_genelist_all = c()
for(i in 1:336){
	print(i)
	start = (i-1)*100+1
	end= i*100
	LIBD_genelist_all = c(LIBD_genelist_all,unlist(get_pesudotime_gene(sim_SpatialPCA,count_use[start:end,])))
}

length(LIBD_genelist_all)
# > length(LIBD_genelist_all)
# [1] 2763

# match ENSG with symbol
gene_dict = read.table("~/gene_annotation_downloaded_from_ENCODE.txt", header=T)
gene_name = as.character(gene_dict$ENSG)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
gene_name02=unlist(gene_name02)
gene_dict$ENSG = gene_name02
LIBD_gene_symbol = gene_dict$gene[match(LIBD_genelist_all, gene_dict$ENSG)]

topGOnum = 10
bg = as.vector(na.omit(gene_dict$gene[match(rownames(count_use), gene_dict$ENSG)] ))
gostres <- gost(query = na.omit(LIBD_gene_symbol), 
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
datt = gostres$result[1:10,]
datt$log10p = -log10(datt$p_value)

p=ggplot(data=datt, aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0("Pseudotime trajectory 1"),x="Biological terms", y = "-log10(p value)")+
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


pdf(paste0("GSEA_SpatialPCA_pseudotime_gene_LIBD_Nov8.pdf"),width=20,height=5)
print(p)
dev.off()

```












