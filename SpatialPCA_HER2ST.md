---
layout: default
title: Breast Tumor Analysis
nav_order: 5
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/HER2ST
---

## Table of Contents
- [Load package](#load-package)
- [Data preprocessing](#data-preprocessing)
- [Run SpatialPCA](#run-spatialpca)
	- [Run other methods](#run-other-methods)
	- [Annotate ground truth for 8 annotated samples](annotate-ground-truth-for-8-annotated-samples)
	- [Eight annotated samples](#eight-annotated-samples)
- [Analysis](#analysis)
	- [Find metagenes](#find-metagenes)
	- [Functions in metagene visualization](#functions-in-metagene-visualization)
	- [Compare ARI](#compare-ari)
	- [Visualize each method](#visualize-each-method)
	- [PseudoR2 and ARI vs PCs](#pseudor2-and-ari-vs-pcs) 
	- [Moran's I of top spatial PCs](#morans-i-of-top-spatial-pcs )
	- [TLS region marker genes](#tls-region-marker-genes)
	- [Deconvolution](#deconvolution)
	- [RGB plots](#rgb-plots)
	- [Trajectory analysis](#trajectory-analysis)
	- [High resolution spatial map reconstrucion](#high-resolution-spatial-map-reconstrucion)


### Load package
```R
library(SpatialPCA)
```

### Data preprocessing
Following data processing codes from [https://github.com/almaan/her2st](https://github.com/almaan/her2st):

ST breast tumor data are downloaded from [https://github.com/almaan/her2st](https://github.com/almaan/her2st). We also saved the raw data that we used in our examples in RData format, which can be downloaded from [here](https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT?usp=sharing).

```R
#   library(Seurat)
  # library(data.table)
  # library(ggplot2)
  # library(plotly)
  # library(STutility)
  # library(zeallot)
  # library(openxlsx)
  # meta_data <- read.csv("../res/ST-cluster/sample.csv",header=TRUE, stringsAsFactors=FALSE,sep=",")
# meta_data$patient_id = c()
# for(i in 1:dim(meta_data)[1]){
#   meta_data$patient_id[i] = paste0(meta_data$patient[i],meta_data$cluster[i] )
# }
# rownames(meta_data) = meta_data$patient_id
# samples <- list.files(pattern = ".tsv", path = "../data/ST-cnts/", full.names = T)
# names(samples) <- substr(do.call(rbind, strsplit(samples, split = "/"))[, 5], start = 1, stop = 2)
# imgs <- list.files(path = "../data/ST-imgs/", recursive = T, full.names = T, pattern = ".jpg")
# names(imgs) <- do.call(rbind, strsplit(imgs, split = "/"))[, 6]
# ids <- names(samples)
# infoTable <- data.frame(samples, imgs = imgs[ids], ids, patient_id = substr(x = ids, start = 1, stop = 1), stringsAsFactors = FALSE)

# tmp = meta_data[match(infoTable$ids, meta_data$patient_id),]

# infoTable <- cbind(infoTable, tmp)
# infoTable[, -c(11:28)]

# #Subset infoTable to include specified datasets.
# infoTable$spotfiles <- list.files(path = "../data/ST-spotfiles", full.names = T)[1:36]
# head(infoTable)

# ## Load data
# #Load all patient datasets and merge into one Seurat object per patient. Each gene has to bre present in at least 20 spots per sample and each spot has to have at least 300 unique features (genes).
# seu.list <- lapply(unique(infoTable$patient_id), function(s) {
#     InputFromTable(infotable = subset(infoTable, patient_id == s), 
#                       min.gene.spots = 20,
#                       min.spot.feature.count = 300,
#                       platform = "1k")
# }
# )

# # remove ring genes
# seu.list <- lapply(seu.list, function(seu) {
#   subset(seu, features = setdiff(rownames(seu@assays$RNA@counts), ring.genes))
# })


# #Calculate some QC metrics
# total.qc <- do.call(rbind, lapply(seu.list, function(se) {
#   data.frame(total_UMIs = sum(se@assays$RNA@counts), nSpots = ncol(se))
# }))


# # QC

# qcMat <- do.call(rbind, lapply(1:length(seu.list), function(i) {
#     seu <- seu.list[[i]]
#     do.call(rbind, lapply(unique(seu[["ids", drop = T]]), function(id) {
#         repMat <- seu@assays$RNA@counts[, seu[["ids", drop = T]] == id]
#         nUMI <- Matrix::colSums(repMat)
#         nGene <- apply(repMat, 2, function(x) sum(x > 0))
#         data.frame(sample = id, 
#                    avg.nUMI = round(mean(nUMI)),
#                    median.nUMI = median(nUMI),
#                    max.nUMI = max(nUMI),
#                    min.nUMI = min(nUMI),
#                    avg.nGene = round(mean(nGene)),
#                    median.nGene = median(nGene),
#                    min.nGene = min(nGene),
#                    max.nGene = max(nGene),
#                    nSpots = ncol(repMat))
#     }))
# }))

# qcMat

# ```

# Prepare count matrix and normalized matrix. We used H1 sample in the manuscript. We took the raw count matrix and location matrix as input in our paper. 
# ```R


# # default variable.features.rv.th is 1.3, original https://github.com/almaan/her2st used 1.1
# seu.list.1.3 = seu.list
# seu.list.1.3 <- lapply(seu.list.1.3, function(seu) {
#   SCTransform(seu, 
#               vars.to.regress = c("ids"), 
#               return.only.var.genes = FALSE, 
#               variable.features.n = NULL, 
#               variable.features.rv.th = 1.3)
# })

# for(num in 1:length(seu.list.1.3)){
#   print(num)
#   seu.list.single = seu.list.1.3[[num]]
#   save(seu.list.single, file = paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
# }

# her2stdatanum = 0
# for(num in 1:8){
#   load(paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
#   count = 0
#   for(id in unique(seu.list.single[["ids", drop = T]])){
#     count = count + 1
#     her2stdatanum = her2stdatanum + 1
#     print(her2stdatanum)
#     SCTcounts = seu.list.single@assays$SCT@counts[, seu.list.single[["ids", drop = T]] == id]
#     SCTscaled = seu.list.single@assays$SCT@scale.data[, seu.list.single[["ids", drop = T]] == id]
#     ind = which(seu.list.single@tools$Staffli@meta.data$sample == count)
#     metadata = seu.list.single@tools$Staffli@meta.data[ind,]
#     print(dim(metadata))
#     save(SCTcounts, SCTscaled, metadata, file = paste0("~/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
#   }
# }
```


### Run SpatialPCA
```R

num= c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8)
clusternum = c( 6,6,6,6,6,6, 5,5,5,5,5,5, 4,4,4,4,4,4, 4,4,4,4,4,4, 4 ,4,4,4,4,4,7,7,7, 7,7,7)
zimu = c(rep("A",6),rep("B",6),rep("C",6),rep("D",6),rep("E",3),rep("F",3),rep("G",3),rep("H",3))

her2stdatanum = c(1:36)
for(kk in 1:36){
print(kk)
load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",her2stdatanum[kk],".rv1.3.RData"))
load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/seu.list.single",num[kk],".rv1.3.RData"))
H_count = seu.list.single@assays$RNA@counts
rawcount = H_count[match(rownames(SCTcounts),rownames(H_count)),match(colnames(SCTcounts),colnames(H_count))]
location=metadata[,5:6] # in our previous version, we used image coordinates rather than array row and column numbers, will update this later. The main results didn't change.
location=as.matrix(location)
ST = CreateSpatialPCAObject(counts=rawcount, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
# ind_na = which(metadata$truth2020=="undetermined")
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = ST@SpatialPCs
SpatialPCA_result$normalized_expr  = ST@normalized_expr
SpatialPCA_result$location = ST@location
pred_cluster= walktrap_clustering(clusternum[kk], ST@SpatialPCs,round(sqrt(dim(ST@location)[1])))
SpatialPCA_result$clusterlabel = pred_cluster
SpatialPCA_result$clusterlabel_refine=refine_cluster_10x(pred_cluster,SpatialPCA_result$location,shape="square")
save(SpatialPCA_result, file = paste0("ST_SpatialPCA_sample",kk,"result.RData"))

}

```

### Run other methods
```R
args <- as.numeric(commandArgs(TRUE))
i = args[1] # data
j = args[2] # method
print(i)
print(j)

dataset="her2st"

suppressMessages(require(SPARK))
library(Seurat)
library(peakRAM)
library(ggplot2)
library(ggpubr)
library(mclust) # ARI
library(aricode)# NMI
library(SpatialPCA)# NMI
library(lisi) # lisi score
library(reticulate) 

source("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility.R")

num= c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8)
clusternum = c( 6,6,6,6,6,6, 5,5,5,5,5,5, 4,4,4,4,4,4, 4,4,4,4,4,4, 4 ,4,4,4,4,4,7,7,7, 7,7,7)
zimu = c(rep("A",6),rep("B",6),rep("C",6),rep("D",6),rep("E",3),rep("F",3),rep("G",3),rep("H",3))
her2stdatanum = c(1:36)
samplename = c(paste0("A",1:6),paste0("B",1:6),paste0("C",1:6),paste0("D",1:6),paste0("E",1:3),paste0("F",1:3),paste0("G",1:3),paste0("H",1:3))
path_anno = "/net/mulan/disk2/shanglu/Projects/spatialPCA/data/her2st/her2st-master/data/ST-pat/lbl"
name = c("A1","B1","C1","D1","E1","F1","G2","H1")

kk=match(name,samplename)[i]

if(j ==1){
    # SpatialPCA
    load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
    load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",kk,".rv1.3.RData"))
    nameid = name[i]
    anno = read.csv(paste0(path_anno,"/",nameid,"_labeled_coordinates.tsv"),sep="\t")
    colnames(anno) = c("Row.names"  , "adj_x" ,"adj_y",  "pixel_x"  ,"pixel_y" ,"label")
    anno$x=round(anno$adj_x)
    anno$y=round(anno$adj_y)
    myanno = merge(metadata, anno, by=c("x","y"))
    anno_label = myanno
    labels = anno_label$label
    table(labels)
    location = myanno[,5:6] 
    # in a previous version of manuscript, we used image col and image row as coordinates, later we switched to array row and col, the main results didn't change, as we scaled the locations in the algorithm
    location[,2]=-location[,2]
	myanno_reorder = myanno
	myanno_reorder$spotid = paste0(myanno$x,"x",myanno$y,"_",myanno$sample)
	ST = CreateSpatialPCAObject(counts=SpatialPCA_result$rawcount, location=SpatialPCA_result$location, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
	ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
	ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
	ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)

	# ind_na = which(metadata$truth2020=="undetermined")
	SpatialPCA_result$SpatialPCs  = ST@SpatialPCs
	SpatialPCA_result$normalized_expr  = ST@normalized_expr
	SpatialPCA_result$location = ST@location
	myanno_reorder = myanno_reorder[match(rownames(SpatialPCA_result$location),myanno_reorder$spotid),]
	SpatialPCA_result$annotation = myanno_reorder
	cluster_num = length(table(SpatialPCA_result$annotation$label))
	ind_na = which(SpatialPCA_result$annotation$label=="undetermined")
	SpatialPCA_result$clusterlabel_walktrap = walktrap_clustering(cluster_num, SpatialPCA_result$SpatialPCs,25)
	SpatialPCA_result$clusterlabel = refine_cluster_10x(SpatialPCA_result$clusterlabel_walktrap ,SpatialPCA_result$location,shape="square")
	SpatialPCA_result$Truth = SpatialPCA_result$annotation$label
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth[-ind_na],"clusterlabel"=SpatialPCA_result$clusterlabel[-ind_na]))

	SpatialPCA_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	SpatialPCA_result$NMI = NMI(tabb[,1],tabb[,2])
	SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
	SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
	SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
	save(SpatialPCA_result, file = paste0("ST_SpatialPCA_sample",kk,"result.RData"))

}else if(j==2){

	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	count_in = SpatialPCA_result$rawcount[,match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))]
	loc_in=SpatialPCA_result$annotation[,1:2]
	clusternum=length(table(SpatialPCA_result$annotation$label))
	count_in = as.matrix(count_in)
	rownames(loc_in) = colnames(count_in)
	res = BayesSpace_func(count_in, loc_in, clusternum=clusternum,nrep=10000,ifLIBD=FALSE,platform="ST")
	ind_na = which(SpatialPCA_result$annotation$label=="undetermined")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$annotation$label[-ind_na],"clusterlabel"=res$clusterlabel[-ind_na]))
	BayesSpace_result = list()
	BayesSpace_result$clusterlabel = as.character(res$clusterlabel)
	BayesSpace_result$location = SpatialPCA_result$location
	BayesSpace_result$Truth = SpatialPCA_result$annotation$label
	BayesSpace_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	BayesSpace_result$NMI = NMI(tabb[,1],tabb[,2])
	BayesSpace_result$CHAOS = fx_CHAOS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$PAS = fx_PAS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$LISI = fx_lisi(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$normalized_expr = SpatialPCA_result$normalized_expr
	BayesSpace_result$rawcount = SpatialPCA_result$rawcount
	BayesSpace_result$rawlocation = SpatialPCA_result$rawlocation
	BayesSpace_result$annotation = SpatialPCA_result$annotation

	save(BayesSpace_result, file = paste0("ST_BayesSpace_sample",kk,"result.RData"))

	}else if(j==3){

	library(reticulate)
	use_python("/net/mulan/home/shanglu/anaconda3/envs/SpaGCN/bin/python3.7", required=TRUE)
	source_python("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility_python.py")
	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	count_in = SpatialPCA_result$rawcount[,match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))]
	loc_in=as.data.frame(SpatialPCA_result$annotation[,1:2])
	rownames(loc_in) = colnames(count_in)
	clusternum=length(table(SpatialPCA_result$annotation$label))
	count_in = t(as.matrix(count_in))
	ind_na = which(SpatialPCA_result$annotation$label=="undetermined")
	colnames(loc_in)=c("x","y") # important
	res = run_SpaGCN_py(count_in, loc_in , clusternum) # cell by gene 
	ind_na = which(SpatialPCA_result$annotation$label=="undetermined")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$annotation$label[-ind_na],"clusterlabel"=res$refined_pred[-ind_na]))
	SpaGCN_result = list()
	SpaGCN_result$clusterlabel = res$refined_pred
	SpaGCN_result$location = res[,c("x","y")]
	SpaGCN_result$Truth = SpatialPCA_result$annotation$label
	SpaGCN_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	SpaGCN_result$NMI = NMI(tabb[,1],tabb[,2])
	SpaGCN_result$CHAOS = fx_CHAOS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$PAS = fx_PAS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$LISI = fx_lisi(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$normalized_expr = SpatialPCA_result$normalized_expr
	SpaGCN_result$rawcount = SpatialPCA_result$rawcount
	SpaGCN_result$rawlocation = SpatialPCA_result$rawlocation
	SpaGCN_result$annotation = SpatialPCA_result$annotation
	save(SpaGCN_result, file = paste0("ST_SpaGCN_sample",kk,"result.RData"))

	}else if(j==4){
	
	library(reticulate)
	library(peakRAM)
	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id = match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))
	loc_in = as.data.frame(SpatialPCA_result$rawlocation[match_id,])
	count_in = as.matrix(SpatialPCA_result$rawcount[,match_id])
	ind_na = which(SpatialPCA_result$annotation$label=="undetermined")
	path_out = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST/sample_annotated",i)
	res = HMRF_func(count_in=t(count_in), location_in=loc_in, clusternum=clusterNum,path_out=path_out,betas = c(0,2,6)) # betas = [start, step, number]
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$annotation$label[-ind_na],"clusterlabel"=res[-ind_na,6]))
	HMRF_result = list()
	HMRF_result$clusterlabel = as.character(res[,6])
	HMRF_result$location = SpatialPCA_result$location
	HMRF_result$Truth = SpatialPCA_result$Truth
	HMRF_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	HMRF_result$NMI = NMI(tabb[,1],tabb[,2])
	HMRF_result$CHAOS = fx_CHAOS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$PAS = fx_PAS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$LISI = fx_lisi(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$normalized_expr = SpatialPCA_result$normalized_expr
	HMRF_result$rawcount = SpatialPCA_result$rawcount
	HMRF_result$rawlocation = SpatialPCA_result$rawlocation
	HMRF_result$rawmeta = SpatialPCA_result$rawmeta
	HMRF_result$res = res
	HMRF_result$annotation = SpatialPCA_result$annotation
	save(HMRF_result, file = paste0("ST_HMRF_sample",kk,"result.RData"))
	
	}else if(j==5){ # NMF

	library(peakRAM)
	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id = match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))
	loc_in = as.data.frame(SpatialPCA_result$rawlocation[match_id,])
	count_in = as.matrix(SpatialPCA_result$rawcount[,match_id])
	res = NMF_cluster_func(count_in=count_in, genelist=rownames(SpatialPCA_result$normalized_expr),PCnum=20,cluster_method="walktrap",clusternum=clusterNum)
	ind_na = which(SpatialPCA_result$Truth=="undetermined")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth[-ind_na],"clusterlabel"=res$clusterlabel[-ind_na]))

	NMF_result = list()
	NMF_result$clusterlabel = res$clusterlabel
	NMF_result$location = SpatialPCA_result$location
	NMF_result$Truth = SpatialPCA_result$Truth
	NMF_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	NMF_result$NMI = NMI(tabb[,1],tabb[,2])
	NMF_result$CHAOS = fx_CHAOS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PAS = fx_PAS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$LISI = fx_lisi(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$normalized_expr = SpatialPCA_result$normalized_expr
	NMF_result$PCs = res$PC
	NMF_result$normalized_expr = SpatialPCA_result$normalized_expr
	NMF_result$rawcount = SpatialPCA_result$rawcount
	NMF_result$rawlocation = SpatialPCA_result$rawlocation
	NMF_result$annotation = SpatialPCA_result$annotation

	save(NMF_result, file = paste0("ST_NMF_sample",kk,"result.RData"))

}else if(j==6){ # PCA

	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id = match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))
	loc_in = as.data.frame(SpatialPCA_result$rawlocation[match_id,])
	count_in = as.matrix(SpatialPCA_result$rawcount[,match_id])
	res = PCA_cluster_func(expr=SpatialPCA_result$normalized_expr,PCnum=20,cluster_method="walktrap",clusternum=clusterNum)
	ind_na = which(SpatialPCA_result$Truth=="undetermined")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth[-ind_na],"clusterlabel"=res$clusterlabel[-ind_na]))
	PCA_result = list()
	PCA_result$clusterlabel = res$clusterlabel
	PCA_result$location = SpatialPCA_result$location
	PCA_result$Truth = SpatialPCA_result$Truth
	PCA_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	PCA_result$NMI = NMI(tabb[,1],tabb[,2])
	PCA_result$CHAOS = fx_CHAOS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PAS = fx_PAS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$LISI = fx_lisi(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$normalized_expr = SpatialPCA_result$normalized_expr
	PCA_result$PCs = res$PC
	PCA_result$normalized_expr = SpatialPCA_result$normalized_expr
	PCA_result$rawcount = SpatialPCA_result$rawcount
	PCA_result$rawlocation = SpatialPCA_result$rawlocation
	PCA_result$annotation = SpatialPCA_result$annotation

	save(PCA_result, file = paste0("ST_PCA_sample",kk,"result.RData"))

}
	

```

### Annotate ground truth for 8 annotated samples
```R

name = c("A1","B1","C1","D1","E1","F1","G2","H1")
samplename = c(paste0("A",1:6),paste0("B",1:6),paste0("C",1:6),paste0("D",1:6),paste0("E",1:3),paste0("F",1:3),paste0("G",1:3),paste0("H",1:3))
path_anno = "/net/mulan/disk2/shanglu/Projects/spatialPCA/data/her2st/her2st-master/data/ST-pat/lbl"

pdf("Groundtruth_8samples.pdf",width=10,height=6)
for(nameid in name){
    kk=which(samplename %in% nameid)
    load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
    load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",kk,".rv1.3.RData"))
    anno = read.csv(paste0(path_anno,"/",nameid,"_labeled_coordinates.tsv"),sep="\t")
    colnames(anno) = c("Row.names"  , "adj_x" ,"adj_y",  "pixel_x"  ,"pixel_y" ,"label")
    anno$x=round(anno$adj_x)
    anno$y=round(anno$adj_y)
    myanno = merge(metadata, anno, by=c("x","y"))
    anno_label = myanno
    labels = anno_label$label
    table(labels)
    location = myanno[,5:6]
    location[,2]=-location[,2]
	myanno_reorder = myanno
	myanno_reorder$spotid = paste0(myanno$x,"x",myanno$y,"_",myanno$sample)
	myanno_reorder = myanno_reorder[match(rownames(SpatialPCA_result$location),myanno_reorder$spotid),]
	SpatialPCA_result$annotation = myanno_reorder
	save(SpatialPCA_result, file = paste0("ST_SpatialPCA_sample",kk,"result.RData"))
    p=plot_cluster(legend="right",location=location,myanno$label,pointsize=6,text_size=40 ,title_in=paste0("Groundtruth ",nameid),color_in=cbp_spatialpca)
    print(p)
}
dev.off()


```

### Eight annotated samples
```R
name = c("A1","B1","C1","D1","E1","F1","G2","H1")
samplename = c(paste0("A",1:6),paste0("B",1:6),paste0("C",1:6),paste0("D",1:6),paste0("E",1:3),paste0("F",1:3),paste0("G",1:3),paste0("H",1:3))

i=0
for(nameid in name){
	i=i+1
        kk=which(samplename %in% nameid)
	kk=match(name,samplename)[i]
	load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",kk,".rv1.3.RData"))
	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST/ST_SpatialPCA_sample",kk,"result.RData"))
	count_in = SpatialPCA_result$rawcount[,match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))]
	loc_in=as.matrix(SpatialPCA_result$annotation[,1:2]) # we updated results using locations as array row and columns
	clusternum=length(table(SpatialPCA_result$annotation$label))
	count_in = as.matrix(count_in)
	rownames(loc_in) = colnames(count_in)
	ST = CreateSpatialPCAObject(counts=count_in, location=loc_in, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
	ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype = "SJ")
	ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
	ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
	ind_na = which(SpatialPCA_result$Truth=="undetermined")
	SpatialPCA_result$SpatialPCs  = ST@SpatialPCs
	SpatialPCA_result$normalized_expr  = ST@normalized_expr
	SpatialPCA_result$location = ST@location
	clusternum=length(table(SpatialPCA_result$annotation$label))
	knn = c(5,18,12,12,11,74,13,25)[i]
	pred_cluster= walktrap_clustering(clusternum, ST@SpatialPCs,knn)
	SpatialPCA_result$clusterlabel = pred_cluster
	SpatialPCA_result$clusterlabel_refine=refine_cluster_10x(pred_cluster,SpatialPCA_result$location,shape="square")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$annotation$label[-ind_na],"clusterlabel"=SpatialPCA_result$clusterlabel[-ind_na]))
	SpatialPCA_result$Truth = SpatialPCA_result$annotation$label
	SpatialPCA_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	print(SpatialPCA_result$ARI)
	SpatialPCA_result$NMI = NMI(tabb[,1],tabb[,2])
	SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel_refine, SpatialPCA_result$location)
	SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel_refine, SpatialPCA_result$location)
	SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel_refine, SpatialPCA_result$location)
	save(SpatialPCA_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST/ST_SpatialPCA_sample",kk,"result.RData"))

}

```

## Analysis

### Find metagenes
```R

# get annotated gene expressions
name = c("A1","B1","C1","D1","E1","F1","G2","H1")
samplename = c(paste0("A",1:6),paste0("B",1:6),paste0("C",1:6),paste0("D",1:6),paste0("E",1:3),paste0("F",1:3),paste0("G",1:3),paste0("H",1:3))
path_anno = "/net/mulan/disk2/shanglu/Projects/spatialPCA/data/her2st/her2st-master/data/ST-pat/lbl"

count_each = list()
annotate_region = list()
gene_list = c()
tmp = 0
for(kk in match(name,samplename)){
	tmp = tmp + 1
	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	count_each[[tmp]] = SpatialPCA_result$rawcount
	gene_list[[tmp]] = rownames(SpatialPCA_result$rawcount)
}

gene_common = Reduce(intersect, gene_list)
count_each = list()
annotate_region = list()
gene_list = c()
sample_id = c()
tmp = 0
for(kk in match(name,samplename)){
	tmp = tmp + 1
	load(paste0("ST_SpatialPCA_sample",kk,"result.RData"))
	count_each[[tmp]] = SpatialPCA_result$rawcount[match(gene_common, rownames(SpatialPCA_result$rawcount)),]
	gene_list[[tmp]] = rownames(SpatialPCA_result$rawcount)
	annotate_region[[tmp]] = SpatialPCA_result$annotation$label
	sample_id = c(sample_id, rep(paste0("sample",tmp),dim(count_each[[tmp]])[2]))
	print(sum(length(annotate_region[[tmp]])==dim(count_each[[tmp]])[2]))
}

count_all = do.call("cbind",count_each)
region_label = unlist(annotate_region)

# create seurat object
Seu <- CreateSeuratObject(counts = count_all, project = "ST", min.cells = 0, min.features = 0)
Seu = SCTransform(Seu, return.only.var.genes = TRUE, variable.features.n = NULL,  variable.features.rv.th = 1.3)	 
Idents(Seu) = region_label
DE_gene = list()
each_num = c()
for(cluster in 1:length(table(region_label))){
	print(cluster)
	DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = names(table(region_label))[cluster], ident.2 =NULL, test.use = "MAST")
	each_num[cluster] = dim(DE_gene[[cluster]])[1]
}

names(DE_gene) = names(table(region_label))

# create metagene as mean expression of genes with positive log2FC in each region annotation

metagene = list()
metagene_count = list()
markergenes = list()
for(cluster in 1:length(table(region_label))){
	print(cluster)
	# positive fold-change:
	index = which(DE_gene[[cluster]]$avg_log2FC>0 & DE_gene[[cluster]]$p_val_adj<0.001)
	markergenes[[cluster]] = genenames=rownames(DE_gene[[cluster]])[index]
	metagene[[cluster]] = colMeans(Seu@assays$SCT@scale.data[na.omit(match(genenames,rownames(Seu@assays$SCT@scale.data))),])
	metagene_count[[cluster]] = colSums(count_all[na.omit(match(genenames,rownames(count_all))),])

}
names(markergenes) = names(table(region_label))
names(metagene_count) = names(table(region_label))
names(metagene) = names(table(region_label))



```

### Functions in metagene visualization
```R

metagene_visualize = function(region_label,dat,name,sample,coloroption){
	figure = list()
	tmp = 0
	figure_count = 0
	# mean normalized metagene expression
	for(regionID in names(table(region_label))){
		figure_count = figure_count + 1
		tmp = tmp + 1
		datt = dat[,c(1,2,tmp+2)]
		datt = as.data.frame(datt)
		colnames(datt) = c("x","y","region")
		figure[[figure_count]] =ggplot(datt, aes(x = x, y = -y, color = region)) +
						geom_point(size=pointsize[sample], alpha = 1) +
						    scale_color_viridis(option=coloroption)+
						  #ggtitle(paste0("sample ",name[sample],": ",regionID))+
						  ggtitle(paste0(regionID))+
						  theme_void()+
						  theme(plot.title = element_text(size = textsize),
						        text = element_text(size = textsize),
						        legend.position = "bottom")
		# visualize mean expression of metagene in each region
	}

	return(figure)
}
	


method_cluster = function(location,cluster,sample,method){
	loc = location
	loc[,2]=-loc[,2]
	figure = plot_cluster(loc,
								clusterlabel = as.character(cluster),
								pointsize=pointsize[sample],
								text_size=textsize,color_in=D3,
								title_in=paste0(method),
								legend="bottom")
}


metagene_cluster = function(region_label,ground_truth,dat,cl,name,sample,method,coloroption){

		var_cluster_mean = c()
		range_cluster_mean = c()
		var_allspots=c()
		figure = list()
		ratio = c()
		diff = c()
		method_ratio = c()
		sample_ratio = c()

		table_dat = table(cl,ground_truth)
		proportion = table_dat/rowSums(table_dat)

		tmp=0
		for(regionID in names(table(ground_truth))){ # find each groundtruth corresponding tissue region
		tmp = tmp + 1
		datt = dat[,c(1,2,which(colnames(dat)==regionID))]
		datt = as.data.frame(datt)
		dattt=datt
		dattt$cluster = unlist(cl)
		dattt$clustermean = dattt[,3]
		dattt$truth = ground_truth

		if( sum(proportion[,tmp]>0.35)>0){
		idx = which.max(proportion[,tmp])
		scaled_clustermean = (dattt$clustermean-min(dattt$clustermean))/(max(dattt$clustermean)-min(dattt$clustermean))
		#scaled_clustermean = scale(dattt$clustermean)
		within_region_mean = mean(scaled_clustermean[which(dattt$cluster == idx)])
		out_region_mean = mean(scaled_clustermean[which(dattt$cluster != idx)])
		ratio_eachregion = within_region_mean/out_region_mean
		diff_eachregion = within_region_mean-out_region_mean
		}else{
		ratio_eachregion= 0 # set NA as 0
		diff_eachregion = 0
		}

		ratio = c(ratio, ratio_eachregion)
		diff = c(diff,diff_eachregion)

		dattt$cluster_mean = dattt$clustermean
		for(clustereach in names(table(dattt$cluster))){
			id_tmp = which(dattt$cluster == clustereach)
			dattt$cluster_mean[id_tmp] = mean(dattt$clustermean[id_tmp])
		}


		figure[[tmp]] =ggplot(dattt, aes(x = x, y = -y, color = cluster_mean)) +
		geom_point(size=pointsize[sample], alpha = 1) +
		    scale_color_viridis(option=coloroption)+
		  #ggtitle(paste0(method,": sample ",name[sample],": ",regionID))+
		  ggtitle(paste0(method))+
		  theme_void()+
		  theme(plot.title = element_text(size = textsize),
		        text = element_text(size = textsize),
		        legend.position = "bottom")

		}

		method_ratio = c(method_ratio, rep(method,length(ratio_eachregion)))
		sample_ratio = c(sample_ratio, rep(sample, length(ratio_eachregion)))
		ratio_dat = data.frame(ratio,method_ratio,sample_ratio,diff)

		return(list("figure"=figure,
		"ratio"=ratio_dat))
}



library(viridis)
textsize=40
coloroption="viridis"
# textsize=25
# coloroption="inferno"
pointsize=c(9,9,10,9,7,7,8,7)

pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Test_HER2ST_DEmarkerGENEs_bycluster_",coloroption,".pdf"),width=40,height=50)

  method_name = c()
  samplenum = c()
  ratios = list()

  for(sample in 1:8){
  	print(sample)
	load(paste0("ST_SpatialPCA_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_BayesSpace_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_SpaGCN_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_HMRF_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_PCA_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_NMF_sample",match(name,samplename)[sample],"result.RData"))
	SpaGCN_result$clusterlabel = as.factor(as.integer(SpaGCN_result$clusterlabel))
	SpatialPCA_result$clusterlabel = SpatialPCA_result$clusterlabel_refine
	BayesSpace_result$location = SpaGCN_result$location = HMRF_result$location = PCA_result$location = NMF_result$location = SpatialPCA_result$location = as.data.frame(SpatialPCA_result$location)

	sample_index = which(sample_id==paste0("sample",sample))
	metagene_each = matrix(0,length(sample_index),length(table(region_label)))
	for(cluster in 1:length(table(region_label))){
		metagene_each[,cluster] = metagene[[cluster]][sample_index]
	}
	colnames(metagene_each) = names(table(region_label))
	dat = cbind(SpatialPCA_result$location,metagene_each)
	clusternum = length(unique(SpatialPCA_result$annotation$label))

	figure_metagene = metagene_visualize(region_label,dat,name,sample,coloroption)
	figure_metagene_cluster = list()
	ratio_eachsample = list()
	nn = 0
	for(method in c("Hand_annotate","SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF")){
		if(method=="Hand_annotate"){
			ground_truth = SpatialPCA_result$Truth
			cl=as.integer(as.factor(SpatialPCA_result$Truth))
			metagene_cluster_res = metagene_cluster(region_label,ground_truth,dat,cl,name,sample,method,coloroption)
			figure_metagene_cluster = c(figure_metagene_cluster, metagene_cluster_res$figure)
			nn = nn + 1
			ratio_eachsample[[nn]] = metagene_cluster_res$ratio
			method_name = c(method_name, rep(method, dim( metagene_cluster_res$ratio)[1]))
			figure_method = method_cluster(SpatialPCA_result$location,cl,sample,method)
		}else{
			ground_truth = SpatialPCA_result$Truth
			cl=get(paste0(method,"_result"))$clusterlabel
			metagene_cluster_res = metagene_cluster(region_label,ground_truth,dat,cl,name,sample,method,coloroption)
			figure_metagene_cluster = c(figure_metagene_cluster, metagene_cluster_res$figure)
			nn = nn + 1
			ratio_eachsample[[nn]] = metagene_cluster_res$ratio
			method_name = c(method_name, rep(method, dim( metagene_cluster_res$ratio)[1]))
			figure_method = method_cluster(SpatialPCA_result$location,cl,sample,method)
		}

	}

	ratios[[sample]] = do.call("rbind",ratio_eachsample)
	print(ggarrange(plotlist=c(figure_metagene, figure_metagene_cluster),
					ncol = 7, nrow = 8))

}

dev.off()


# calculate metagene enrichment score in all samples
ratio_table = do.call("rbind",ratios)
# all samples
taba = ratio_table %>% group_by(method_ratio) %>%  summarise(median=median(na.omit(ratio)))
taba = taba[order(-taba$median),]
as.data.frame(taba)

# > as.data.frame(taba)
#    method_ratio    median
# 1 Hand_annotate 1.4386369
# 2    SpatialPCA 1.1405555
# 3    BayesSpace 1.0816380
# 4           NMF 0.9828833
# 5          HMRF 0.9155259
# 6        SpaGCN 0.8216376
# 7           PCA 0.0000000


# calculate metagene enrichment score in sample 34
ratio_table_sample34 = ratio_table[which(ratio_table$sample_ratio==8),]
taba = ratio_table_sample34 %>% group_by(method_ratio) %>%  summarise(median=median(na.omit(ratio)))
taba = taba[order(-taba$median),]
as.data.frame(taba)

# > as.data.frame(taba)
#    method_ratio   median
# 1    SpatialPCA 1.417109
# 2 Hand_annotate 1.380699
# 3           NMF 1.345310
# 4    BayesSpace 1.269250
# 5           PCA 1.248220
# 6          HMRF 1.202340
# 7        SpaGCN 1.095634


# order metagene enrichment score
order_dat = list()
for(sample_index in 1:8){
ratio_table_sample34 = ratio_table[which(ratio_table$sample_ratio==sample_index),]
taba = ratio_table_sample34 %>% group_by(method_ratio) %>%  summarise(median=median(na.omit(ratio)))
taba = taba[order(-taba$median),]
taba$order = c(1:7)
order_dat[[sample_index]] = taba
}
order_dat = do.call("rbind",order_dat)
order_dat$method_ratio = factor(order_dat$method_ratio, levels=c("Hand_annotate","SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF"),order=T)
taba = order_dat %>% group_by(method_ratio) %>%  summarise(median_score=median(order),
	mean_score=mean(order)
	)
taba = taba[order(taba$median_score),]
taba

# > taba
# # A tibble: 7 × 3
#   method_ratio  median_score mean_score
#   <ord>                <dbl>      <dbl>
# 1 Hand_annotate          1         1.75
# 2 SpatialPCA             3         3.75
# 3 NMF                    3.5       4   
# 4 PCA                    4         3.88
# 5 BayesSpace             4.5       4.12
# 6 HMRF                   5.5       5   
# 7 SpaGCN                 6.5       5.5 
```


### Compare ARI
```R

library(reshape2)

dats = list()
p1=list()
p2=list()

tmp = 0

for(sample in 1:8){

  	print(sample)

	load(paste0("ST_SpatialPCA_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_BayesSpace_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_SpaGCN_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_HMRF_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_PCA_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_NMF_sample",match(name,samplename)[sample],"result.RData"))
	SpaGCN_result$clusterlabel = as.factor(as.integer(SpaGCN_result$clusterlabel))

	ARI = c()
	NMI = c()
	CHAOS = c()
	PAS = c()
	LISI = c()
	methods = c()
	methods_LISI = c()

	for(method in c("SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF")){
		ARI = c(ARI, get(paste0(method,"_result"))$ARI)
		NMI = c(NMI, get(paste0(method,"_result"))$NMI)
		CHAOS = c(CHAOS, get(paste0(method,"_result"))$CHAOS)
		PAS = c(PAS, get(paste0(method,"_result"))$PAS)
		LISI = c(LISI, get(paste0(method,"_result"))$LISI)
		methods = c(methods, method)
		methods_LISI = c(methods_LISI, rep(method,length(get(paste0(method,"_result"))$LISI)))
	}


	dat_LISI = data.frame(LISI,methods_LISI)
	dat_LISI$methods_LISI = factor(dat_LISI$methods_LISI , levels=c("SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF"),order=T)

	dat = data.frame("methods"=c("SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF"),ARI, NMI, CHAOS, PAS)
	dat <- melt(dat, id.vars = c("methods"))
	dat$variable = as.character(dat$variable)
	dat$methods = factor(dat$methods , levels=c("SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF"),order=T)
	dat$sample = sample
	dats[[sample]] = dat

	save(dat, file = paste0("Sample_",sample,"_ARINMICHAOSPAS.RData"))

		p1[[sample]] =ggplot(dat, aes(x=methods, y=value,fill=methods)) +
	  	  geom_bar(stat="identity", color="black",width=0.8)+
	      #scale_fill_brewer(palette="Paired")+
	      #scale_fill_distiller(palette = "Oranges")+
	      scale_fill_manual(values = D3)+
	      labs(title=paste0("Sample ",name[sample]))+
	      theme(legend.position="right") +
	      theme_classic()+
	      #geom_hline(yintercept=median(dat[1,metrics+1]), linetype="dashed", color = "red",size=1)+
	      #scale_fill_manual(values = method_color)+
	      theme(axis.text.x = element_text(angle = 60,  hjust=1))+
	      theme(plot.title = element_text(size = 30),
	                  text = element_text(size = 30),
	                  #axis.title = element_text(face="bold"),
	                  #axis.text.x=element_text(size = 20) ,
	                  legend.position = "right")+
	      facet_wrap(~variable)

	
  

  p2[[sample]] =ggplot(dat_LISI, aes(x=methods_LISI, y=LISI,fill=methods_LISI)) +
  	  geom_boxplot()+
      #scale_fill_brewer(palette="Paired")+
      #scale_fill_distiller(palette = "Oranges")+
      scale_fill_manual(values = D3)+
      labs(title=paste0("Sample ",name[sample]," LISI"))+
      theme(legend.position="right") +
      theme_classic()+
      geom_hline(yintercept=median(dat_LISI[which(dat_LISI$methods_LISI=="SpatialPCA"),1]), linetype="dashed", color = "red",size=1)+
      #scale_fill_manual(values = method_color)+
      theme(axis.text.x = element_text(angle = 60,  hjust=1))+
      theme(plot.title = element_text(size = 30),
                  text = element_text(size = 30),
                  #axis.title = element_text(face="bold"),
                  #axis.text.x=element_text(size = 20) ,
                  legend.position = "right")
# print(p2[[sample]])

}

pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/ARINMICHAOSPAS.pdf",width=20,height=40)
print(ggarrange(plotlist=c(p1),
					ncol = 2, nrow = 4))
dev.off()


pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/LISI_8samples.pdf",width=15,height=30)
print(ggarrange(plotlist=c(p2),
					ncol = 2, nrow = 4))
dev.off()


# result in sample 8( the H1 sample we used in the main analysis)
# > sample
# [1] 8
# > dat
#       methods variable     value sample
# 1  SpatialPCA      ARI 0.4452831      8
# 2  BayesSpace      ARI 0.4180627      8
# 3      SpaGCN      ARI 0.3763917      8
# 4        HMRF      ARI 0.3500935      8
# 5         PCA      ARI 0.3420037      8
# 6         NMF      ARI 0.3197780      8
# 7  SpatialPCA      NMI 0.5471441      8
# 8  BayesSpace      NMI 0.5448721      8
# 9      SpaGCN      NMI 0.4811621      8
# 10       HMRF      NMI 0.4607784      8
# 11        PCA      NMI 0.4043444      8
# 12        NMF      NMI 0.3962226      8
# 13 SpatialPCA    CHAOS 0.1378253      8
# 14 BayesSpace    CHAOS 0.1393070      8
# 15     SpaGCN    CHAOS 0.1525325      8
# 16       HMRF    CHAOS 0.1486891      8
# 17        PCA    CHAOS 0.1483046      8
# 18        NMF    CHAOS 0.1536805      8
# 19 SpatialPCA      PAS 0.1976936      8
# 20 BayesSpace      PAS 0.1762768      8
# 21     SpaGCN      PAS 0.2652389      8
# 22       HMRF      PAS 0.2718287      8
# 23        PCA      PAS 0.3542010      8
# 24        NMF      PAS 0.4365733      8


# result in all 8 annotated samples 
order_dats = do.call("rbind",dats)
order_dats$methods = factor(order_dats$methods, levels=unique(order_dats$methods),order=T)
# order_dats$methods = factor(order_dats$methods, levels=rev(unique(order_dats$methods)),order=T)
order_dats_ari = order_dats[which(order_dats$variable=="ARI"),]

pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Supplementary_Boxplot_metagene_ARI.pdf",width=6,height=5)
p2=ggplot(order_dats_ari, aes(x=methods, y=value)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6,fill = cbp[1:6])+
  #stat_summary(fun=mean, geom="point", shape=20, size=10, color="red", fill="red") +
  # geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  #geom_hline(yintercept=median(SpatialPCA_ratio), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=20)+
  #ylim(0,1)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("ARI in 8 annotated samples"),
       x="", y = "ARI")
  #coord_flip()
print(p2)
dev.off()


taba = order_dats_ari %>% group_by(methods) %>%  summarise(median=median(value))
# > taba
# # A tibble: 6 × 2
#   methods    median
#   <ord>       <dbl>
# 1 SpatialPCA  0.300
# 2 BayesSpace  0.234
# 3 SpaGCN      0.258
# 4 HMRF        0.212
# 5 PCA         0.140
# 6 NMF         0.210

```

### Visualize each method
```R

plot_cluster = function (location, clusterlabel, pointsize = 3, text_size = 15, shape=16,title_in, color_in, legend = "none"){
    cluster = clusterlabel
    loc_x = location[, 1]
    loc_y = location[, 2]
    datt = data.frame(cluster, loc_x, loc_y)
    p = ggplot(datt, aes(x = location[, 1], y = location[, 2], 
        color = cluster)) + geom_point(alpha = 1, size = pointsize,shape=shape) + 
        scale_color_manual(values = color_in) + ggtitle(paste0(title_in)) + 
        theme_void() + theme(plot.title = element_text(size = text_size, 
        face = "bold"), text = element_text(size = text_size), 
        legend.position = legend)
    p
}


sample=8
	load(paste0("ST_SpatialPCA_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_BayesSpace_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_SpaGCN_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_HMRF_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_PCA_sample",match(name,samplename)[sample],"result.RData"))
	load(paste0("ST_NMF_sample",match(name,samplename)[sample],"result.RData"))
	SpaGCN_result$clusterlabel = as.factor(as.integer(SpaGCN_result$clusterlabel))

pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/ST_Fig5b_clusters_square.pdf"),width=30,height=5)

    legend_pos = "none"
    point_size = 6
    text_size_in = 30
    loc = SpatialPCA_result$location
    loc[,2]=-loc[,2]
    shapein=15

    cbp  <- c(  "mediumaquamarine", "chocolate1","dodgerblue",  "#F0E442","palegreen4","lightblue2","plum1","red","#CC79A7","mediumpurple","seagreen1")
 	p1=plot_cluster(legend=legend_pos,location=loc,SpatialPCA_result$clusterlabel_refine,pointsize=point_size,text_size=text_size_in ,title_in=paste0("SpatialPCA"),color_in=cbp,shape=shapein)

    cbp  <- c(  "mediumaquamarine", "lightblue2","palegreen4",  "#F0E442","plum1","dodgerblue","chocolate1","red","#CC79A7","mediumpurple","seagreen1")
    p2=plot_cluster(location=loc,clusterlabel=BayesSpace_result$clusterlabel,pointsize=point_size,text_size=text_size_in ,title_in=paste0("BayesSpace"),color_in=cbp,legend=legend_pos,shape=shapein)
    
    cbp  <- c(  "mediumaquamarine", "palegreen4","chocolate1",  "dodgerblue","lightblue2","plum1","#F0E442","red","#CC79A7","mediumpurple","seagreen1")
    p3=plot_cluster(location=loc,clusterlabel=SpaGCN_result$clusterlabel,pointsize=point_size,text_size=text_size_in ,title_in=paste0("SpaGCN"),color_in=cbp,legend=legend_pos,shape=shapein)
    
    cbp  <- c(  "palegreen4", "#F0E442","chocolate1",  "dodgerblue","mediumaquamarine","lightblue2","plum1","red","#CC79A7","mediumpurple","seagreen1")
    p4=plot_cluster(location=loc,clusterlabel=HMRF_result$clusterlabel,pointsize=point_size,text_size=text_size_in ,title_in=paste0("HMRF"),color_in=cbp,legend=legend_pos,shape=shapein)
    
    cbp  <- c(  "dodgerblue", "plum1","chocolate1",  "mediumaquamarine","#F0E442","palegreen4","lightblue2","red","#CC79A7","mediumpurple","seagreen1")
    p5=plot_cluster(location=loc,clusterlabel=PCA_result$clusterlabel,pointsize=point_size,text_size=text_size_in ,title_in=paste0("PCA"),color_in=cbp,legend=legend_pos,shape=shapein)
    
    cbp  <- c(  "palegreen4", "chocolate1","lightblue2",  "dodgerblue","mediumaquamarine","plum1","#F0E442","palegreen4","red","#CC79A7","mediumpurple","seagreen1")
    p6=plot_cluster(location=loc,clusterlabel=NMF_result$clusterlabel,pointsize=point_size,text_size=text_size_in ,title_in=paste0("NMF"),color_in=cbp,legend=legend_pos,shape=shapein)
 
    print(ggarrange(p1, p2,p3,p4,p5,p6, ncol = 6, nrow = 1))
    
dev.off()

```

### PseudoR2 and ARI vs PCs

```R

library(mclust)
library(nnet) # multinom function
library(DescTools) # PseudoR2 function

#--------------------------------
# Pseudo R2 change with PC number
#--------------------------------

SpaGCN_embed = read.csv("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/SpaGCN_embed/SpaGCN_ST_embed.csv",header=F)

SpatialPCA_R2=c()
PCA_R2=c()
NMF_R2=c()
SpaGCN_R2=c()

for(i in 1:20){
	if(i==1){
		fit <- multinom(SpatialPCA_result$Truth ~ SpatialPCA_result$SpatialPCs[1:i,], maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpatialPCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$Truth ~ PCA_result$PCs[1:i,], maxit = 1000, MaxNWts = 2000,model = TRUE)
		PCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$Truth~ SpaGCN_embed[,1:i], maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpaGCN_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$Truth~ NMF_result$PCs[1:i,], maxit = 1000, MaxNWts = 2000,model = TRUE)
		NMF_R2[i] = PseudoR2(fit,c("McFaddenAdj"))

	}else{
		fit <- multinom(SpatialPCA_result$Truth ~ as.matrix(as.data.frame(t(SpatialPCA_result$SpatialPCs[1:i,]))), maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpatialPCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$Truth ~ as.matrix(as.data.frame(t(PCA_result$PCs[1:i,]))), maxit = 1000, MaxNWts = 2000,model = TRUE)
		PCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$Truth~ as.matrix(as.data.frame(SpaGCN_embed[,1:i])), maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpaGCN_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$Truth~ as.matrix(as.data.frame(t(NMF_result$PCs[1:i,]))), maxit = 1000, MaxNWts = 2000,model = TRUE)
		NMF_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
	}
}

fit <- multinom(SpatialPCA_result$Truth ~ BayesSpace_result$clusterlabel, maxit = 1000, MaxNWts = 2000,model = TRUE)
BayesSpace_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(HMRF_result$Truth ~ HMRF_result$clusterlabel, maxit = 1000, MaxNWts = 2000,model = TRUE)
HMRF_R2 = PseudoR2(fit,c("McFaddenAdj"))



PseudoR2=c(SpatialPCA_R2, PCA_R2, NMF_R2,SpaGCN_R2,rep(BayesSpace_R2,20),rep(HMRF_R2,20))
method=c(rep("SpatialPCA",20),rep("PCA",20),rep("NMF",20),rep("SpaGCN",20),rep("BayesSpace",20),rep("HMRF",20))
PCnum=c(1:20,1:20,1:20,1:20,1:20,1:20)
method=factor(method, levels=c("SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF"))
dat=data.frame(PseudoR2, method, PCnum)

pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Fig5_PseudoR2.pdf"),width=10,height=5)
p<- ggplot(dat, aes(x=PCnum, y=PseudoR2, color=method,group=method)) + 
   #theme_bw(base_size = 22)+
   theme_bw(base_size=25)+
   geom_line(size=1.5)+
  #ylim(0,1)+
  scale_colour_manual(values=c("#264653", "#219ebc", "#ffb703" ,"#fb8500", "#ffafcc" ,"#3a86ff"))+
  geom_point( size=1.5, color="black")+
  labs(title=paste0("Pseudo R2"),x="Low dimensional components number", y = "Pseudo R2")+
   theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              legend.position = "right")# +
print(p)
dev.off()


SpatialPCA_R2
PCA_R2
NMF_R2
SpaGCN_R2


# > SpatialPCA_R2
#  [1] 0.2960768 0.4215008 0.4793948 0.5396736 0.5505426 0.5929037 0.5991844
#  [8] 0.6107530 0.6148023 0.6241743 0.6246817 0.6246543 0.6274298 0.6388320
# [15] 0.6359479 0.6399741 0.6430037 0.6456862 0.6546647 0.6660248
# > PCA_R2
#  [1] 0.2552689 0.3574143 0.4097657 0.4571835 0.4644409 0.4763442 0.4767187
#  [8] 0.4829494 0.4827772 0.4811759 0.4780582 0.4850526 0.4817216 0.4784048
# [15] 0.4759626 0.4742448 0.4787981 0.4751303 0.4727905 0.4766256
# > NMF_R2
#  [1] 0.1482301 0.2260518 0.2406384 0.2587563 0.2732112 0.2963646 0.3237220
#  [8] 0.3332762 0.3566377 0.3619794 0.3685203 0.3852103 0.4441064 0.4556819
# [15] 0.4740000 0.4695245 0.4708137 0.4727923 0.4719767 0.4718803
# > SpaGCN_R2
#  [1] 0.2360305 0.3685295 0.4185808 0.4672689 0.4910370 0.4987424 0.4983591
#  [8] 0.5067656 0.5054406 0.5035914 0.5003980 0.4971404 0.4963184 0.4911112
# [15] 0.4883992 0.4862277 0.4826128 0.4791062 0.4769042 0.4729235



#--------------------------------
# ARI change with PC number
#--------------------------------

ind_na = which((SpatialPCA_result$Truth=="undetermined")) 
BayesSpace_ari_addPCs = adjustedRandIndex(SpatialPCA_result$Truth[-ind_na],BayesSpace_result$clusterlabel[-ind_na])
HMRF_ari_addPCs = adjustedRandIndex(SpatialPCA_result$Truth[-ind_na],HMRF_result$clusterlabel[-ind_na])
SpatialPCA_ari_addPCs=c()
PCA_ari_addPCs=c()
NMF_ari_addPCs=c()
SpaGCN_ari_addPCs=c()
BayesSpace_ari_addPCs = rep(BayesSpace_ari_addPCs, 19)
HMRF_ari_addPCs = rep(HMRF_ari_addPCs, 19)


for(i in 2:20){

		print(i)

		labels=walktrap_clustering(7,SpatialPCA_result$SpatialPCs[1:i,],25)
		labels = refine_cluster_10x(labels ,SpatialPCA_result$location,shape="square")
		SpatialPCA_ari_addPCs[i] = adjustedRandIndex(SpatialPCA_result$Truth[-ind_na],labels[-ind_na])

		labels=walktrap_clustering(7,PCA_result$PCs[1:i,],25)
		PCA_ari_addPCs[i] = adjustedRandIndex(SpatialPCA_result$Truth[-ind_na],labels[-ind_na])

		labels=walktrap_clustering(7,NMF_result$PCs[1:i,],25)
		NMF_ari_addPCs[i] = adjustedRandIndex(SpatialPCA_result$Truth[-ind_na],labels[-ind_na])

		labels=walktrap_clustering(7,t(SpaGCN_embed[,1:i]),25)
		SpaGCN_ari_addPCs[i] = adjustedRandIndex(SpatialPCA_result$Truth[-ind_na],labels[-ind_na])

}


ARIs = c(SpatialPCA_ari_addPCs[-1], PCA_ari_addPCs[-1], NMF_ari_addPCs[-1],SpaGCN_ari_addPCs[-1],BayesSpace_ari_addPCs,HMRF_ari_addPCs)
method=c(rep("SpatialPCA",19),rep("PCA",19),rep("NMF",19),rep("SpaGCN",19),rep("BayesSpace",19),rep("HMRF",19))
PCnum=c(2:20,2:20,2:20,2:20,2:20,2:20)
method=factor(method, levels=c("SpatialPCA","BayesSpace","SpaGCN","HMRF","PCA","NMF"))
dat=data.frame(ARIs, method, PCnum)

pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Fig1f_ARI_vs_PCs.pdf"),width=10,height=5)
p<- ggplot(dat, aes(x=PCnum, y=ARIs, color=method,group=method)) + 
   theme_bw(base_size=25)+
   geom_line(size=1.5)+
  ylim(0,0.5)+
  scale_colour_manual(values=c("#264653", "#219ebc", "#ffb703" ,"#fb8500", "#ffafcc" ,"#3a86ff", "#43AA8B"))+
  geom_point( size=1.5, color="black")+
  labs(title=paste0("ARI"),x="Low dimensional components number", y = "ARI")+
   theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              legend.position = "right")# +
print(p)
dev.off()


```


### Moran's I of top spatial PCs 
```R

get_moranI = function(count_in,location){
    # count_in: gene by spot
    count_in = as.matrix(count_in)
    if(length(which(rowSums(count_in)==0))>0){
        count_in = count_in[-which(rowSums(count_in)==0),]
    }
    library(moranfast)
    output <- apply(count_in, 1, function(x) moranfast(x, location[,1], location[,2]))
    out =  sapply(output, '[[', 1)
    return(out)
}   



spatialpc_moranI = get_moranI(SpatialPCA_result$SpatialPCs[1:5,], SpatialPCA_result$location)
pc_moranI = get_moranI(PCA_result$PCs[1:5, ], SpatialPCA_result$location)


> spatialpc_moranI
[1] 0.18784394 0.14004474 0.12060055 0.17407981 0.09836492
> pc_moranI
[1] 0.15514189 0.10754111 0.08593695 0.13617112 0.05794850


(spatialpc_moranI-pc_moranI)/pc_moranI
[1] 0.2107880 0.3022438 0.4033609 0.2783900 0.6974543

summary((spatialpc_moranI-pc_moranI)/pc_moranI)
> summary((spatialpc_moranI-pc_moranI)/pc_moranI)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.2108  0.2784  0.3022  0.3784  0.4034  0.6975 

```

### TLS region marker genes
```R

#----------------------------------------------------------------
# TLS region identification
#----------------------------------------------------------------

mapDrugToColor<-function(annotations){
    colorsVector = ifelse(annotations["category"]=="Cluster1", 
        "red", ifelse(annotations["category"]=="Cluster2", 
        "orange",ifelse(annotations["category"]=="Cluster3", 
        "yellow",ifelse(annotations["category"]=="Cluster4", 
        "green",ifelse(annotations["category"]=="Cluster5", 
        "blue",ifelse(annotations["category"]=="Cluster6", 
        "purple",ifelse(annotations["category"]=="Cluster7", 
           "skyblue",ifelse(annotations["category"]=="Cluster8", 
        "black","grey"))))))))
    return(colorsVector)
}

testHeatmap3<-function(logCPM, annotations) {    
    sampleColors = mapDrugToColor(annotations)
    # Assign just column annotations
    heatmap3(logCPM, margins=c(10,10), 
    	ColSideColors=sampleColors,scale="none",
    	col = colorRampPalette(c( "#0072B2","#F0E442", "#D16103"))(1024),
    	Rowv=NA,
    	Colv=NA,
        xlab = "Cell ID",
    ylab = "Marker genes",
    showColDendro = F,
  showRowDendro = F) 
    #Assign column annotations and make a custom legend for them
    heatmap3(logCPM, margins=c(10,10), ColSideColors=sampleColors, 
    	scale="none",
    	col = colorRampPalette(c( "#0072B2", "#F0E442","#D16103"))(1024),
        legendfun=function()showLegend(legend=paste0("Cluster",1:7), col=c("red", "orange", "yellow","green","blue","purple","skyblue"), cex=1),
            	Rowv=NA,
    	Colv=NA,
        xlab = "Cell ID",
    ylab = "Marker genes",
    showColDendro = F,
  showRowDendro = F)
    
    #Assign column annotations as a mini-graph instead of colors,
    #and use the built-in labeling for them
    ColSideAnn<-data.frame(Cluster=annotations[["category"]])
    heatmap3(logCPM, ColSideAnn=ColSideAnn,
        #ColSideFun=function(x)showAnn(x),
        margins=c(10,10),
        ColSideWidth=0.8,
        Rowv=NA,
    	Colv=NA,
        xlab = "Cell ID",
    	ylab = "Marker genes",
    	showColDendro = F,
  		showRowDendro = F)
}


# library(BiRewire)
library(corrplot)
library(heatmap3)

# use all genes normalized expression to make marker gene heatmap
seu <- CreateSeuratObject(counts = SpatialPCA_result$rawcount, project = "forheatmap", min.cells = 0, min.features = 0)
seu=   SCTransform(seu,          
              return.only.var.genes = FALSE, 
              variable.features.n = NULL, 
              variable.features.rv.th = 1.3)     
SCTscaled = seu@assays$SCT@scale.data

load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/HER2ST/SpatialPCA/metadataST.RData")
marker_TLS = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10",
    "CXCL11","CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1") # breast cancer TLS genes
heat_expr = SCTscaled[which(rownames(SCTscaled) %in% marker_TLS),]
anno = paste0("Cluster",as.character(SpatialPCA_result$clusterlabel_refine))
myorder = order(anno)
category = anno[myorder]
gAnnotationData = data.frame("cells"=colnames(heat_expr)[myorder], category)
colnames(gAnnotationData) = c("cells","category")
gLogCpmData = heat_expr[,myorder]
genenum = dim(heat_expr)[1]
datt = matrix(0,genenum,7)
for(i in 1:7){
    for(j in 1:genenum){
        datt[j,i] = mean(gLogCpmData[j,which(category %in% paste0("Cluster",i))])
    }
}
colnames(datt) = paste0(1:7)
rownames(datt) = rownames(heat_expr)
category = paste0("Cluster",1:7)
gAnnotationData_mean = data.frame("cells"=colnames(datt), category)
colnames(gAnnotationData_mean) = c("cells","category")
na.omit(match(marker_TLS,rownames(datt)))
datt = datt[rev(na.omit(match(marker_TLS,rownames(datt)))),]

pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/TLS_marker_heatmap_breast_genes.pdf",width=6, height=6)  #heat_expr = scale(SCTscaled)
testHeatmap3(datt, gAnnotationData_mean)
dev.off()



#------------------------------------------------------------
# Obtain TLS score from original paper
# I rewrote their python code in R
#------------------------------------------------------------

H1_proportion = read.table(file = '/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/H1-proportion.tsv', sep = '\t', header = TRUE)
H1_proportion_607 = H1_proportion[match(colnames(SpatialPCA_result$rawcount),paste0(H1_proportion$X,"_1")),-1]

n_spots = 607
jprod = c()
pos_1=1
pos_2=8
for(s in 1:n_spots){
	vec = as.matrix(unlist(H1_proportion_607[s,]),8,1)
	prod = vec %*% t(vec)
	nprod = prod / sum(prod)
	N=8
	jprod[s] = nprod[pos_1,pos_2] * 2
}

plot_factor_value= function (location, feature, textmethod, pointsize = 2, textsize = 15,shape=15) 
{
    location = as.data.frame(location)
        locc1 = location[, 1]
        locc2 = location[, 2]
        datt = data.frame(feature, locc1, locc2)
        p = ggplot(datt, aes(x = locc1, y = locc2, color = feature)) + 
            geom_point(size = pointsize, alpha = 1,shape=shape) + 
            #scale_color_viridis(option="magma")+
			scale_color_gradient2( low="#05a8aa",mid="#edede9", high="#bc412b")+
            ggtitle(paste0(textmethod)) + 
            theme_void() + 
            theme(plot.title = element_text(size = textsize), 
            text = element_text(size = textsize), legend.position = "right")
  
    return(p)
}


pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/ST_Fig_TLS_score.pdf"),width=6,height=5)

    point_size = 6
    text_size_in = 30
    loc = SpatialPCA_result$location
    shapein=15
    loc[,2]=-loc[,2]
    p1=plot_factor_value(loc,jprod,"TLS score",pointsize=6,textsize = 15,shape=16)
	print(p1)
dev.off()

#------------------------------------------------------------
# Boxplot for TLS score in each spatial domain
#------------------------------------------------------------
sample=8
load(paste0("ST_SpatialPCA_sample",match(name,samplename)[sample],"result.RData"))
cbp=c(  "mediumaquamarine", "chocolate1","dodgerblue",  "#F0E442","palegreen4","lightblue2","plum1","red","#CC79A7","mediumpurple","seagreen1")
pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/BOXplot_TLSscore_in_spatialdomain.pdf",width=6,height=4)
dat_TLS = data.frame(loc,"TLSscore"=jprod,"SpatialDomain"=SpatialPCA_result$clusterlabel_refine)
	p2=ggplot(dat_TLS, aes(x=SpatialDomain, y=TLSscore,fill=SpatialDomain)) +
	  # geom_boxplot(alpha = 0.6,fill = "lightgreen")+
	  geom_boxplot(alpha = 1)+
	  #geom_bar(alpha = 1,position="dodge", stat="identity")+
	  # scale_fill_brewer(palette="Set3")+
	  scale_fill_manual(values = cbp[1:7])+
	  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
	  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
	  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
	  #ylim(0,1)+
	  #geom_hline(yintercept=dat1_dataframe$NMI[1], linetype="dashed", color = "red",size=1)+
	  theme_bw(base_size=25)+
	  theme(plot.title = element_text(size = 22),
	              text = element_text(size = 22),
	              axis.text.x = element_text(angle = 60,  hjust=1),
	              #axis.title = element_text(face="bold"),
	              #axis.text.x=element_text(size = 20) ,
	              legend.position = "none")+
	    #facet_wrap(~genetype, ncol = 1)+
	  labs(title=paste0("TLS score"),
	       x="Spatial Domains", y = paste0("TLS score"))
print(p2)
dev.off()


taba = dat_TLS %>% group_by(SpatialDomain) %>%  summarise(median=median(na.omit(TLSscore)))
as.data.frame(taba)

# > as.data.frame(taba)
#   SpatialDomain       median
# 1             1 0.0001256088
# 2             2 0.0021306888
# 3             3 0.0023813420
# 4             4 0.0480196187
# 5             5 0.0010816596
# 6             6 0.0002913864
# 7             7 0.0022554624

```

### Deconvolution
```R

#----------------------------------------------------------------
# Deconvolution
#----------------------------------------------------------------

# Wu et al. EMBO 2020
# Stromal cell diversity associated with immune evasion in human triple-negative breast cancer
# https://www.embopress.org/doi/epdf/10.15252/embj.2019104063
# found at https://singlecell.broadinstitute.org/single_cell/study/SCP1106/stromal-cell-diversity-associated-with-immune-evasion-in-human-triple-negative-breast-cancer#study-download



library(Matrix)
library(readr)

# Read in `matrix.mtx`
counts <- readMM("./counts_matrix/matrix.mtx.gz")
# Read in `genes.tsv`
genes <- read_tsv("./counts_matrix/features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1
# Read in `barcodes.tsv`
cell_ids <- read_tsv("./counts_matrix/barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

meta_ref = read.csv("Wu_EMBO_metadata.csv")
meta_ref_singlecell = meta_ref[-1,]

library(RCTD)

# reference single cell data
meta_data <- meta_ref_singlecell # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$celltype_final; names(cell_types) <- as.character(meta_data$NAME) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- as.integer(meta_data$nCount_RNA); names(nUMI) <- as.character(meta_data$NAME)  # create nUMI named list
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
#> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
#>             is intended, there is no problem.
## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type
# > print(dim(reference@counts)) #observe Digital Gene Expression matrix
# [1] 28118 24271
# > table(reference@cell_types) #number of occurences for each cell type

#                   B_Cells              CD4+ T-cells              CD8+ T-cells 
#                      1245                      2003                      3691 
#                      dPVL               Endothelial          Epithelial_Basal 
#                       214                       610                      4095 
#  Epithelial_Basal_Cycling Epithelial_Luminal_Mature                     iCAFs 
#                       614                       277                      1129 
#                     imPVL                    myCAFs                   Myeloid 
#                       106                       280                      4606 
#             Myoepithelial                  NK cells                 NKT cells 
#                       212                       358                       164 
#              Plasma_Cells        T_cells_unassigned           T-cells Cycling 
#                      1955                       938                       605 
#                    T-Regs                 Tfh cells 
#                       994                       175 
#----------------------------------------------------------------
# run RCTD 
#----------------------------------------------------------------
# spatial data
library(DropletUtils)
num=8
her2stdatanum=34
load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/seu.list.single",num,".rv1.3.RData"))
H_count = seu.list.single@assays$RNA@counts
rawcount = H_count[match(rownames(SCTcounts),rownames(H_count)),match(colnames(SCTcounts),colnames(H_count))]

location=metadata[,5:6]
location=as.matrix(location)
counts <- rawcount
coords <- as.data.frame(location)
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results
# > table(results$results_df$first_type)

#                   B_Cells              CD4+ T-cells              CD8+ T-cells 
#                        32                         4                         1 
#                      dPVL               Endothelial          Epithelial_Basal 
#                         7                         9                       112 
#  Epithelial_Basal_Cycling Epithelial_Luminal_Mature                     iCAFs 
#                       233                         0                        14 
#                     imPVL                    myCAFs                   Myeloid 
#                         0                        26                         5 
#             Myoepithelial                  NK cells                 NKT cells 
#                        17                         0                         0 
#              Plasma_Cells        T_cells_unassigned           T-cells Cycling 
#                       120                        27                         0 
#                    T-Regs                 Tfh cells 
#                         0                         0 

metadataST$celltype = as.character(results$results_df$first_type)

library(ggplot2)
library(viridis)
library(hrbrthemes)
library(ggthemes)

cbp2 = c("#FFDB6D", "#C4961A", "#F4EDCA", "tomato","#C3D7A4",  "#4E84C4","#52854C",
"#D16103", "deepskyblue1", "cadetblue3","lightblue1","plum1","chartreuse3")

preparedata = function(percentage,celltypenames){
rownames(percentage) = paste0("Cluster",1:7)
celltype = celltypenames
colnames(percentage) = celltype
rownames(percentage) = paste0("Cluster",1:7)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:7)),length(celltype)))
CellType = c(rep(celltype,each=7))
datt = data.frame(cluster_vec, percentage_vec,CellType)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:7))
return(datt)
}

makefigure = function(datt){
p=ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = CellType)) +        ## global aes
  scale_fill_manual(values=cbp2)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_bw()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
return(p)
}

#--------- Stack_barplot SpatialPCA

method="SpatialPCA"
percentage = matrix(0,7,13)
metadataST$celltype=as.factor(metadataST$celltype)
for(k in 1:7){
metadata_sub = metadataST[which(SpatialPCA_result$clusterlabel_refine ==k),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}
celltypenames =  names(table(metadataST$celltype))
datt=preparedata(percentage,celltypenames)
pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Stack_barplot_",method,".pdf"),width=8,height=5)
makefigure(datt)+ggtitle(paste0(method))+  theme(axis.text.x = element_text(angle = 60,  hjust=1))
percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = metadataST[which(SpatialPCA_result$clusterlabel_refine==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadataST)[1]*100,2)
}
datt=preparedata(percentage,celltypenames)
makefigure(datt)+ggtitle(paste0(method))+  theme(axis.text.x = element_text(angle = 60,  hjust=1))
dev.off()

method_clusters = data.frame(
	"SpatialPCA"=SpatialPCA_result$clusterlabel_refine,
	"BayesSpace"=BayesSpace_result$clusterlabel,
    "SpaGCN"=SpaGCN_result$clusterlabel,
    "HMRF"=HMRF_result$clusterlabel,
    "PCA"=PCA_result$clusterlabel,
    "NMF"=NMF_result$clusterlabel
    )

rownames(method_clusters) = rownames(SpatialPCA_result$location)
metaRCTD_ST = metadataST[,c("pixel_x","pixel_y","celltype")]
metaRCTD_ST$celltype = as.factor(metaRCTD_ST$celltype)
metaRCTD_ST$spotID = rownames(metaRCTD_ST)
method_clusters$spotID = rownames(method_clusters)
metaRCTD_ST = merge(metaRCTD_ST,method_clusters, by="spotID" )

save(metaRCTD_ST, file = "/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/metaRCTD_ST.RData")

cbp2 = c("#FFDB6D", "#C4961A", "#F4EDCA", "tomato","#C3D7A4",  "#4E84C4","#52854C",
"#D16103", "deepskyblue1", "cadetblue3","lightblue1","plum1","chartreuse3")

preparedata = function(percentage){
rownames(percentage) = paste0("Cluster",1:7)
celltype = names(table(metaRCTD_ST$celltype))
colnames(percentage) = celltype
rownames(percentage) = paste0("Cluster",1:7)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:7)),length(celltype)))
CellType = c(rep(celltype,each=7))
datt = data.frame(cluster_vec, percentage_vec,CellType)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:7))
return(datt)
}

makefigure= function(datt){
p=ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = CellType)) +        ## global aes
  scale_fill_manual(values=cbp2)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_bw()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 40),
              text = element_text(size = 40),
              #axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 30,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "none")# +
return(p)
}

#----------------------------------------------------------------
# make figures for all methods
#----------------------------------------------------------------

p1=list()
p2=list()
count = 0
for(method in 5:10){
    count = count + 1
percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = metaRCTD_ST[which(metaRCTD_ST[,method]==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}
datt=preparedata(percentage)
p1[[count]] = makefigure(datt)+ggtitle(paste0(colnames(metaRCTD_ST)[method]))
percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = metaRCTD_ST[which(metaRCTD_ST[,method]==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metaRCTD_ST)[1]*100,2)
}
datt=preparedata(percentage)
datt$cluster_vec = as.character(datt$cluster_vec)
datt$cluster_vec = factor(datt$cluster_vec, levels=c(paste0("Cluster",c(1:7))),order=T)
p2[[count]]=makefigure(datt)+ggtitle(paste0(colnames(metaRCTD_ST)[method]))#+coord_flip()
}


pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/ST_Cell_type_Stack_barplot.pdf"),width=40,height=15)
print(ggarrange(plotlist=c(p1,p2),
    ncol = 6, nrow = 2))
dev.off()


```

### RGB plots
```R

#----------------------------------------------------------------
# RGB
#----------------------------------------------------------------



library(pdist)

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


loc = SpatialPCA_result$location
loc[,2]=-loc[,2]


method_color=c("#1F77B4", "#FF7F0E", "#2CA02C" ,"#D62728", "#9467BD" ,"#8C564B", "#E377C2","#7F7F7F", "#BCBD22", "#17BECF")


library(ggpubr)

set.seed(2333)
p1 = plot_RGB_tSNE(loc,SpatialPCA_result$SpatialPCs,pointsize=5,textsize=20)
p2 = plot_RGB_tSNE(loc,PCA_result$PCs,pointsize=5,textsize=20)
p3 = plot_RGB_tSNE(loc,NMF_result$PCs,pointsize=5,textsize=20)
save(p1,p2,p3,file=paste0("tSNE_1x.RData"))
pdf("ST_RGB_tSNE_1x.pdf",width=15, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

set.seed(2333)
p1 = plot_RGB_UMAP(loc,SpatialPCA_result$SpatialPCs,pointsize=5,textsize=20)
p2 = plot_RGB_UMAP(loc,PCA_result$PCs,pointsize=5,textsize=20)
p3 = plot_RGB_UMAP(loc,NMF_result$PCs,pointsize=5,textsize=20)
save(p1,p2,p3,file=paste0("UMAP_1x.RData"))
pdf("ST_RGB_UMAP_1x.pdf",width=15, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

set.seed(2333)
p1 = plot_RGB_tSNE(loc,10*SpatialPCA_result$SpatialPCs,pointsize=5,textsize=20)
p2 = plot_RGB_tSNE(loc,10*PCA_result$PCs,pointsize=5,textsize=20)
p3 = plot_RGB_tSNE(loc,10*NMF_result$PCs,pointsize=5,textsize=20)
save(p1,p2,p3,file=paste0("tSNE_10x.RData"))
pdf("ST_RGB_tSNE_10x.pdf",width=15, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

set.seed(2333)
p1 = plot_RGB_UMAP(loc,10*SpatialPCA_result$SpatialPCs,pointsize=5,textsize=20)
p2 = plot_RGB_UMAP(loc,10*PCA_result$PCs,pointsize=5,textsize=20)
p3 = plot_RGB_UMAP(loc,10*NMF_result$PCs,pointsize=5,textsize=20)
save(p1,p2,p3,file=paste0("UMAP_10x.RData"))
pdf("ST_RGB_UMAP_10x.pdf",width=15, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

set.seed(2333)
p1 = plot_RGB_tSNE(loc,20*SpatialPCA_result$SpatialPCs,pointsize=5,textsize=20)
p2 = plot_RGB_tSNE(loc,20*PCA_result$PCs,pointsize=5,textsize=20)
p3 = plot_RGB_tSNE(loc,20*NMF_result$PCs,pointsize=5,textsize=20)
save(p1,p2,p3,file=paste0("tSNE_20x.RData"))
pdf("ST_RGB_tSNE_20x.pdf",width=15, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

set.seed(2333)
p1 = plot_RGB_UMAP(loc,20*SpatialPCA_result$SpatialPCs,pointsize=5,textsize=20)
p2 = plot_RGB_UMAP(loc,20*PCA_result$PCs,pointsize=5,textsize=20)
p3 = plot_RGB_UMAP(loc,20*NMF_result$PCs,pointsize=5,textsize=20)
save(p1,p2,p3,file=paste0("UMAP_20x.RData"))
pdf("ST_RGB_UMAP_20x.pdf",width=15, height=5)
ggarrange(p1[[2]], p2[[2]], p3[[2]], 
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()


load(paste0("tSNE_1x.RData"))
dat = get_RGB_var(p1,p2,p3)
dat$method = factor(dat$method,levels=c("SpatialPCA","PCA","NMF"),order=T)
pdf("ST_tSNE_RGB_varsum.pdf",width=7,height=7)
p=ggplot(dat, aes(x=method, y=var_sum)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6,fill = method_color[1:3])+
  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  geom_hline(yintercept=median(dat$var_sum[which(dat$method=="SpatialPCA")]), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=25)+
  #ylim(0,0.4)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("ST: RGB in tSNE"),
       x="", y = "Variance in RGB")
 p
dev.off()


load(paste0("UMAP_1x.RData"))
dat = get_RGB_var(p1,p2,p3)
dat$method = factor(dat$method,levels=c("SpatialPCA","PCA","NMF"),order=T)
pdf("ST_UMAP_RGB_varsum.pdf",width=7,height=7)
p=ggplot(dat, aes(x=method, y=var_sum)) +
  #geom_boxplot(alpha = 0.6,fill = "lightgreen")+
  geom_boxplot(alpha = 0.6,fill = method_color[1:3])+
  #geom_abline(intercept = coefs[1], slope = coefs[2],color="orange",size=2)+
  #geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  geom_hline(yintercept=median(dat$var_sum[which(dat$method=="SpatialPCA")]), linetype="dashed", color = "red",size=1)+
  #geom_jitter(position = position_jitter(w = 0.1, h = 0),fill="#D55E00",color="#D55E00", size=1)+
  theme_bw(base_size=25)+
  #ylim(0,0.6)+
  theme(axis.text.x = element_text(angle = 60,  hjust=1))+
  labs(title=paste0("ST: RGB in UMAP"),
       x="", y = "Variance in RGB")
 p
dev.off()


```

### Trajectory analysis
```R

#------------------------------------------------------------
 #     trajectory on tumor and surrounding region
#------------------------------------------------------------


library(slingshot)

# SpatialPCA
tumor_ind = which(SpatialPCA_result$clusterlabel_refine %in% c(2,3,7))
sim_tumor <- SingleCellExperiment(assays = SpatialPCA_result$rawcount[,tumor_ind])
reducedDims(sim_tumor) <- SimpleList(DRM = t(SpatialPCA_result$SpatialPCs[,tumor_ind]))
colData(sim_tumor)$Walktrap <- factor(SpatialPCA_result$clusterlabel_refine[tumor_ind])    
sim_tumor  <-slingshot(sim_tumor, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="2" ) 
summary(sim_tumor@colData@listData)


pseudotime_traj1_tumor = sim_tumor@colData@listData$slingPseudotime_1
clusterlabels_tumor = SpatialPCA_result$clusterlabel_refine[tumor_ind]
gridnum = 10
color_in=c(  "chocolate1", "dodgerblue","plum1",  "black")
pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Trajectory_SpatialPCA_tumor_region.pdf",width=8,height=5)
p_traj1 = plot_trajectory(pseudotime_traj1_tumor, loc[tumor_ind,],clusterlabels_tumor,gridnum,color_in,pointsize=7 ,arrowlength=0.3,arrowsize=1,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj1[[1]],
          ncol = 2, nrow = 1))
dev.off()   


# visualize on whole tissue
tumor_ind = which(SpatialPCA_result$clusterlabel_refine %in% c(2,3,7))
pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
pseudotime_traj1[-tumor_ind]=NA
pseudotime_traj1[tumor_ind]=pseudotime_traj1_tumor
clusterlabels = SpatialPCA_result$clusterlabel_refine
gridnum = 10
color_in=c(  "mediumaquamarine", "chocolate1","dodgerblue",  "#F0E442","palegreen4","lightblue2","plum1","black","#CC79A7","mediumpurple","seagreen1")
pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Trajectory_SpatialPCA_focus_on_tumor.pdf",width=10,height=5)
p_traj1 = plot_trajectory(pseudotime_traj1, loc,clusterlabels,gridnum,color_in,pointsize=5 ,arrowlength=0.3,arrowsize=1.3,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj1[[1]],
          ncol = 2, nrow = 1))
dev.off()   
```

### High resolution spatial map reconstrucion
```R

	kk=34
	load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",kk,".rv1.3.RData"))
	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST/ST_SpatialPCA_sample",kk,"result.RData"))
	count_in = SpatialPCA_result$rawcount[,match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$rawcount))]
	loc_in=as.matrix(SpatialPCA_result$annotation[,1:2]) # array row and columns
	clusternum=length(table(SpatialPCA_result$annotation$label))
	count_in = as.matrix(count_in)
	rownames(loc_in) = colnames(count_in)
	ST = CreateSpatialPCAObject(counts=count_in, location=loc_in, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
	ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype = "SJ")
	ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
	ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)

STsimu_high_ST = SpatialPCA_highresolution(ST, platform="ST",newlocation=NULL)
cluster_SpatialPCA_high = walktrap_clustering(7, latent_dat=STsimu_high_ST@highPCs,200)
location = STsimu_high_ST@highPos
color_in=c(  "palegreen4", "chocolate1","plum1",  "#F0E442","mediumaquamarine","dodgerblue","lightblue2")
pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/High_resolution_cluster_software_test.pdf",width=5,height=5)
pointsize=1.5
text_size=20
legend="bottom"
shape = 16
p = plot_cluster(location, as.character(cluster_SpatialPCA_high), pointsize=pointsize,text_size=text_size ,title_in=paste0("High resolution spatial domain"),color_in,shape=shape,legend=legend)
print(p)
dev.off()


#-------------------------
# High resolution trajectory
#-------------------------

library(slingshot)
sim<- SingleCellExperiment(assays = STsimu_high_ST@highPCs)
reducedDims(sim) <- SimpleList(DRM = t(STsimu_high_ST@highPCs))
colData(sim)$Walktrap <- factor(cluster_high)    
sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="2" ) 
summary(sim@colData@listData)
tumor_ind = which(cluster_high %in% c(2,3,6))
sim_tumor <- SingleCellExperiment(assays = STsimu_high_ST@highPCs[,tumor_ind])
reducedDims(sim_tumor) <- SimpleList(DRM = t(STsimu_high_ST@highPCs[,tumor_ind]))
colData(sim_tumor)$Walktrap <- factor(cluster_high[tumor_ind])    
sim_tumor  <-slingshot(sim_tumor, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="2" ) 
summary(sim_tumor@colData@listData)

pseudotime_traj1_tumor = sim_tumor@colData@listData$slingPseudotime_1
clusterlabels_tumor = cluster_high[tumor_ind]
gridnum = 10
color_in=c( "chocolate1","plum1","dodgerblue",  "black")
# visualize on whole tissue
tumor_ind = which(cluster_high %in% c(2,3,6))
pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
pseudotime_traj1[-tumor_ind]=NA
pseudotime_traj1[tumor_ind]=pseudotime_traj1_tumor
# pseudotime_traj2 = sim@colData@listData$slingPseudotime_2
# pseudotime_traj2[-tumor_ind]=NA
# pseudotime_traj2[tumor_ind]=pseudotime_traj2_tumor
clusterlabels = cluster_high
gridnum = 20
location = STsimu_high_ST@highPos
location[,2]=-location[,2]

color_in=c(  "palegreen4", "chocolate1","plum1",  "#F0E442","mediumaquamarine","dodgerblue","lightblue2","black")
pdf("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/ST_update/Trajectory_SpatialPCA_high_focus_on_tumor.pdf",width=10,height=5)
p_traj1 = plot_trajectory(pseudotime_traj1, location,clusterlabels,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj1[[1]],
          ncol = 2, nrow = 1))
dev.off()   


```








