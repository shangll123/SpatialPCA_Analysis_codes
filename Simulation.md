---
layout: default
title: Simulation
nav_order: 6
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/Simulation
---


## Table of Contents

- [Prepare simulation datasets](#prepare-simulation-datasets)
- [Generate single cell data](#generate-single-cell-data)
- [Main functions in simulation](#main-functions-in-simulation)
- [Simulate data, 10000 cells](#simulate-data-10000-cells)
- [Simulate data, spot level](#simulate-data-spot-level)
- [Simulate data, spot level with stripe pattern](#simulate-data-spot-level-with-stripe-pattern)
- [Simulate data, spot level with arbitrary spatial correlation](#simulate-data-spot-level-with-arbitrary-spatial-correlation)
- [High resolution spatial map reconstruction](#high-resolution-spatial-map-reconstruction)

## Prepare simulation datasets

The locations are generated from pixels in a image (LIBDsimu.jpg), which can be downloaded from [here](https://drive.google.com/drive/folders/18rwQjB3-g86A-M9xYPPJlHz60bfMABdE?usp=sharing).

```R
library("grid")
library("gridExtra")

library(jpeg)
mandrill = readJPEG("LIBDsimu.jpg", native = FALSE)
# > dim(mandrill)
# [1] 1100  984    3
# copy the image three times
mandrill.R = mandrill
mandrill.G = mandrill
mandrill.B = mandrill
# zero out the non-contributing channels for each image copy
mandrill.R[,,2:3] = 0
mandrill.G[,,1]=0
mandrill.G[,,3]=0
mandrill.B[,,1:2]=0

df = data.frame(
  red = matrix(mandrill[,,1], ncol=1),
  green = matrix(mandrill[,,2], ncol=1),
  blue = matrix(mandrill[,,3], ncol=1)
)

### compute the k-means clustering, to obtain regions from image by color
K = kmeans(df,5)

df$label = K$cluster

table(df$label)

### Replace the color of each pixel in the image with the mean 
### R,G, and B values of the cluster in which the pixel resides:
# get the coloring
colors = data.frame(
  label = 1:nrow(K$centers), 
  R = K$centers[,"red"],
  G = K$centers[,"green"],
  B = K$centers[,"blue"]
)
# merge color codes on to df
# IMPORTANT: we must maintain the original order of the df after the merge!
df$order = 1:nrow(df)
df = merge(df, colors)
df = df[order(df$order),]
df$order = NULL

# get mean color channel values for each row of the df.
R = matrix(df$R, nrow=dim(mandrill)[1])
G = matrix(df$G, nrow=dim(mandrill)[1])
B = matrix(df$B, nrow=dim(mandrill)[1])
  
# reconstitute the segmented image in the same shape as the input image
mandrill.segmented = array(dim=dim(mandrill))
mandrill.segmented[,,1] = R
mandrill.segmented[,,2] = G
mandrill.segmented[,,3] = B

RGB_label = R*B*G

# View the result
pdf("segmented.pdf")
grid.raster(mandrill.segmented)
dev.off()

#save(df, mandrill.segmented, R,G,B, file = "From_image_region.RData")

#> df[1:4,]
#       label red green blue         R         G         B
#200066     2   1     1    1 0.9997333 0.9994392 0.9992982
#200067     2   1     1    1 0.9997333 0.9994392 0.9992982
#200068     2   1     1    1 0.9997333 0.9994392 0.9992982
#200069     2   1     1    1 0.9997333 0.9994392 0.9992982

pixel_ind = c(1:dim(df)[1])
x_coor = rep(1:dim(mandrill.segmented)[2], each=dim(mandrill.segmented)[1])
y_coor = -rep(1:dim(mandrill.segmented)[1], times=dim(mandrill.segmented)[2])

data_groundtruth = data.frame(pixel_ind, x_coor,y_coor )
data_groundtruth$label = as.integer(as.factor(c(RGB_label)))

# remove background 
data_groundtruth_use = data_groundtruth[-which(data_groundtruth$label==5),]

set.seed(1234)
subsample = data_groundtruth_use[sample(1:dim(data_groundtruth_use)[1],10000,replace=F),]
# > min(subsample$y_coor)
# [1] -1591
subsample$y_coor = subsample$y_coor+1591 # make coordinates >=0
save(data_groundtruth, data_groundtruth_use,subsample, file = "simulate_spatial_cell_label.RData")

save(subsample, file = "LIBDsubsample.RData") 

```


## Generate single cell data

Single cell data are downloaded from [GSE104276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276).

```R
library(Matrix)
library(readr)
library(openxlsx)

PFC_ref = read.table("GSE104276_all_pfc_2394_UMI_count_NOERCC.txt")
PFC_ref = as.matrix(PFC_ref)
celltype_ref = read.xlsx("GSE104276_readme_sample_barcode.xlsx",5)
PFC_ref_use = PFC_ref[,match(as.character(celltype_ref[,1]), colnames(PFC_ref))]
raw.data = PFC_ref_use
cell_types <- celltype_ref$cell_types
names(cell_types) <- as.character(celltype_ref[,1]) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- colSums(PFC_ref_use); names(nUMI) <- as.character(celltype_ref[,1]) # create nUMI named list
library(RCTD)
reference <- Reference(raw.data, cell_types, nUMI)
print(dim(reference@counts))
table(reference@cell_types)
#> table(reference@cell_types)
#       Astrocytes GABAergic neurons         Microglia           Neurons 
#               76               701                68              1057 
#              OPC        Stem cells 
#              117               290 
# > head(subsample)
#        pixel_ind x_coor y_coor label
# 858114    858114    781   1477     3
# 723450    723450    658    841     1
# 380568    380568    346    523     2
# 371366    371366    338    925     2
# 169923    169923    155   1068     1
# 466557    466557    425   1434     3

library(splatter)
cnts=as.matrix(reference@counts) 
# [1] 24153  2309
celltype=reference@cell_types
init_params <- splatEstimate(cnts[,which(reference@cell_types=="Neurons")])
save(init_params, file = "init_params_LIBD.RData")

```

## Main functions in simulation
code_utility.R
```R
"%&%" <- function(x, y) paste0(x, y)

BayesSpace_read <- function(count_in, location_in,platform="ST") {
    spot = rownames(location_in)
    in_tissue = rep(1,length(spot))
    row = location_in[,1]
    col = location_in[,2]
    imagerow = location_in[,1]
    imagecol = location_in[,2]
    colData=data.frame(spot,in_tissue, row, col, imagerow, imagecol)
    colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
    rownames(colData) <- colData$spot
    colData <- colData[colData$in_tissue > 0, ]
    gene_id = gene_name = rownames(count_in)
    feature_type = rep("",length(gene_id))
    rowData = data.frame(gene_id, gene_name, feature_type)
    colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
    rowData <- rowData[, c("gene_id", "gene_name")]
    rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
    sce <- SingleCellExperiment(assays=list(counts=count_in),
                                rowData=rowData,
                                colData=colData)
    sce <- scater::logNormCounts(sce)
    metadata(sce)$BayesSpace.data <- list()
    print("BayesSpace platform")
    print(platform)
    metadata(sce)$BayesSpace.data$platform <- platform
    metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    
    return(sce)
}


BayesSpace_func = function(count_in, location_in, clusternum,nrep=50000,ifLIBD = FALSE,LIBD_sampleid=NULL,platform = "Visium"){
	
	if(ifLIBD==TRUE){

		library(igraph)
		library(peakRAM)
		library(SingleCellExperiment)
		library(ggplot2)
		library(BayesSpace)
		library('spatialLIBD')
		# library(BiocFileCache)

		clusternums = c(7,7,7,7,5,5,5,5,7,7,7,7)
		# sce <- fetch_data(type = 'sce')
		sce <- spatialLIBD::fetch_data(type = 'sce',destdir="/net/mulan/disk2/shanglu/Projects/SpatialPCA/LIBD_dat")
		metaData = SingleCellExperiment::colData(sce)
		sample_names <- paste0( unique(colData(sce)$sample_name))

		dlpfc <- getRDS("2020_maynard_prefrontal-cortex", sample_names[i],cache=FALSE)
		load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/LIBD/LIBD_sample",LIBD_sampleid,".RData") )
		ground_truth = KRM_manual_layers_sub$layer_guess_reordered

		# use top 2000 genes
		set.seed(101)
		dec <- scran::modelGeneVar(dlpfc)
		top <- scran::getTopHVGs(dec, n = 2000)
		set.seed(102)
		dlpfc <- scater::runPCA(dlpfc, subset_row=top)
		## Add BayesSpace metadata
		dlpfc <- spatialPreprocess(dlpfc, platform="Visium", skip.PCA=TRUE)
		#q <- length(table(ground_truth))  # Number of clusters
		q=clusternums[LIBD_sampleid]
		d <- 15  # Number of PCs
		## Run BayesSpace clustering
		set.seed(104)
		dlpfc <- spatialCluster(dlpfc, q=q, d=d, platform='Visium',
		                        nrep=50000, gamma=3, save.chain=TRUE) 

		clusterlabel = dlpfc$spatial.cluster
		loc1 = unlist(dlpfc@colData@listData$imagecol)
		loc2 = unlist(dlpfc@colData@listData$imagerow)
		location = data.frame(loc1, loc2)

		return(list("clusterlabel"=clusterlabel,"location"=location))

	}else if(ifLIBD==FALSE){

		library(peakRAM)
		library(SingleCellExperiment)
		library(ggplot2)
		library(BayesSpace)
		library(Matrix)
    print("platform")
    print(platform)
    platform_in = platform
		sce = BayesSpace_read(count_in,location_in,platform_in)
		set.seed(1234)
		top <- scran::getTopHVGs(sce, n = 2000)
		sce <- scater::runPCA(sce, subset_row=top)
		## Add BayesSpace metadata
		sce <- spatialPreprocess(sce, platform=platform_in, skip.PCA=TRUE)
		q <- clusternum # Number of clusters
		d <- 15  # Number of PCs
		## Run BayesSpace clustering
		sce <- spatialCluster(sce, q=q, d=d, platform=platform_in,
		                        nrep=nrep, gamma=3, save.chain=FALSE)  # slow
		clusterlabel = sce$spatial.cluster
		loc1 = unlist(sce@colData@listData$imagecol)
		loc2 = unlist(sce@colData@listData$imagerow)
		location = data.frame(loc1, loc2)

		return(list("clusterlabel"=clusterlabel,"location"=location))
	}

}


HMRF_func = function(count_in, location_in, clusternum,path_out,betas=c(28,2,2)){
	
	# HMRF, giotto
	library(Giotto)
	library(peakRAM)
	library(tidyverse)
	my_python_path = "/net/mulan/home/shanglu/anaconda3/envs/SpaGCN/bin/python3.7"
	instrs <- createGiottoInstructions(python_path = my_python_path)
	hmrf <- createGiottoObject(raw_exprs = as.data.frame(t(count_in)), 
	spatial_locs = location_in, instructions = instrs)
	hmrf <- filterGiotto(gobject = hmrf, expression_threshold = 0.5, 
	gene_det_in_min_cells = 20, min_det_genes_per_cell = 0)
	hmrf <- normalizeGiotto(hmrf) 
	hmrf <- addStatistics(hmrf)
	hmrf <- createSpatialNetwork(hmrf, minimum_k = 2)
	# select SE genes using
	genes <- binSpect(hmrf, bin_method = 'kmeans')$genes
	genes <- genes[1:100]
	outfile <- file.path(path_out)
	R=clusternum
	verbose <- doHMRF(
	gobject = hmrf, 
	spatial_genes = genes, 
	expression_values = "scaled", 
	k = clusternum, 
	betas=betas,
	#betas = c(0, 2, 26), # betas = [start, step, number]
	output_folder = outfile,
	python_path = my_python_path, 
	overwrite_output = T)
	files <- list.files(file.path(outfile, "result.spatial.zscore", 
	"k_" %&% R), pattern = "\\.prob\\.txt", full.names = T)
	labels <- lapply(files, function(f){
	res.i <- read.table(f, row.names = 1)
	apply(res.i, 1, which.max)
	})
	labels <- do.call(cbind, labels)
	betas <- sapply(strsplit(files, "beta\\.|\\.0\\.prob"), `[`, 2)
	colnames(labels) <- paste0("beta", betas)
	labels <- labels[, order(as.numeric(betas))]

	return(labels)

}


get_NMF = function(count, PCnum){
  suppressMessages(require(scater))
  suppressMessages(require(NMF))
  suppressMessages(require(SingleCellExperiment))
  suppressMessages(require(RcppML))
  #expr = log(count+1) # non negative
  sce <- SingleCellExperiment(list(counts=count))
  sce <- logNormCounts(sce)
  #sce <- runNMF(sce,ncomponents = PCnum)
  res <- scater::calculateNMF(sce, ncomponents = PCnum)
  Z_NMF = t(res)
  return(Z_NMF)
}



NMF_cluster_func = function(count_in, genelist,PCnum,cluster_method="walktrap",clusternum){
# install.packages("RcppML")
# devtools::install_github("zdebruine/RcppML")
	library(RcppML)
	library(SpatialPCA)
	count_use = count_in[match(genelist,rownames(count_in)),]
	#res <- RcppML::nmf(count_use,  k = PCnum, nonneg = TRUE)
	# NMF_pcs = res$h

  NMF_pcs = get_NMF(count_use,PCnum)

	if(cluster_method=="walktrap"){
		clusterlabel = walktrap_clustering(clusternum=clusternum,latent_dat=as.matrix(NMF_pcs),knearest=round(sqrt(dim(NMF_pcs)[2])) ) 
	}else if(cluster_method=="louvain"){
		clusterlabel = louvain_clustering(clusternum=clusternum,latent_dat=as.matrix(NMF_pcs),knearest=round(sqrt(dim(NMF_pcs)[2])) ) 
	}

	return(list("PC"=NMF_pcs,"clusterlabel"=clusterlabel))

}



PCA_cluster_func = function(expr, PCnum,cluster_method="walktrap",clusternum){
	output=expr
	d=PCnum
	n=dim(expr)[2]
	k = dim(output)[1]
	output_sub_mean=matrix(0,k,n)
	for(i_k in 1:k){
  		output_sub_mean[i_k,]=output[i_k,]-mean(output[i_k,])
	}
	svd_output_sub_mean=svd(output_sub_mean)
	A_ini=svd_output_sub_mean$u[,1:d] 
	Z_pca = t(A_ini) %*% output_sub_mean

	if(cluster_method=="walktrap"){
		clusterlabel = walktrap_clustering(clusternum=clusternum,latent_dat=as.matrix(Z_pca),knearest=round(sqrt(dim(Z_pca)[2])) ) 
	}else if(cluster_method=="louvain"){
		clusterlabel = louvain_clustering(clusternum=clusternum,latent_dat=as.matrix(Z_pca),knearest=round(sqrt(dim(Z_pca)[2])) ) 
	}

	return(list("PC"=Z_pca,"clusterlabel"=clusterlabel))
}


fx_lisi = function(clusterlabel_in, location){
	library(lisi)
	dat = data.frame("clusterlabel"=clusterlabel_in)
	lisi=compute_lisi(location, dat, c("clusterlabel"))
	return("lisi"=lisi$clusterlabel)
}


	library(BayesSpace)
	library(assertthat)
	library(RCurl)

getRDS <- function(dataset, sample, cache=TRUE) {


    url <- "https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SpatialTranscriptomes/%s/%s.rds"
    url <- sprintf(url, dataset, sample)
    assert_that(url.exists(url), msg="Dataset/sample not available")
    
    if (cache) {
        bfc <- BiocFileCache()
        local.path <- bfcrpath(bfc, url)
    } else {
        local.path <- tempfile(fileext=".rds")
        download.file(url, local.path, quiet=TRUE, mode="wb")
    }
    
    readRDS(local.path)
}

simu_stripe = function(
  location,
  label,
  init_params,
  scenario,
  J,
  batch_facLoc,
  de_prop,
  de_facLoc,
  de_facScale,
  sim_seed,
  debug = FALSE
  )
{ 

  N <- nrow(location)
  set.seed(sim_seed)
  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  params <- setParams(
    init_params,
    batchCells = rep(3 * N, 1),
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    nGenes = J,
    group.prob = c(0.25, 0.25, 0.25,0.25),
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)

  sim_groups <- splatSimulate(
    params = params,
    method = "groups",
    verbose = FALSE)

  # remove cells having no expressed genes
    idx_zerosize <- apply(counts(sim_groups), MARGIN = 2, sum) == 0
    sim_groups <- sim_groups[, !idx_zerosize]

   if(scenario == 1){
      prop <- c(0.6, 0.3, 0.05,0.05)
    }else if(scenario == 2){
      prop <- c(0.45, 0.45, 0.05,0.05)
    }else if(scenario == 3){
      prop <- c(0.35, 0.30, 0.30,0.05)
    }

  # 3.generate cell types
  print("generate cell types")

    ztrue <- label
    c_simu <- rep(NA, length(ztrue))
    if(scenario == 1){ # scenario 1: 4 cell types
      print("scenario == 1")
      c_simu <- ztrue # cell type assignment same as region assignment
    }else if(scenario != 1){
      print("scenario != 1")
      for(z in unique(label)){
        zi_idx <- ztrue == z # zi_idx is index for region z
        c_simu[zi_idx] <- sample(map_z2c(z), sum(zi_idx), prob = prop, replace = T)
      }
    }

    # 4.assign count data
    groups <- as.data.frame(colData(sim_groups)) %>% filter(Batch == "Batch" %&% 1)
    sim_cnt <- array(NA, c(J, N))
    for(c in 1:4){
      c_size <- sum(c_simu == c)  # true number of cells in cell type c
      c_cells <- groups$Cell[grepl(c, groups$Group)] # generated cells in group c
      cells_select <- sample(as.character(c_cells), c_size, replace = F)
      # sample same number of group c cells in real data from generated cells
      sim_cnt[, c_simu == c] <- as.matrix(counts(sim_groups)[, cells_select])
      # for positions of original cell type c cells, assign generated counts of group c
    }
    colnames(sim_cnt) <- "Cell" %&% 1:N
    rownames(sim_cnt) <- "Gene" %&% 1:J

  return(list(sim_cnt, c_simu, sim_seed))
}




map_z2c = function(z)
{
  case_when(
    z == 1 ~ c(1, 2, 3, 4),
    z == 2 ~ c(2, 3, 4, 1),
    z == 3 ~ c(3, 4, 1, 2),
    z == 4 ~ c(4, 1, 2, 3)
    )
}


louvain_clustering_new=function(clusternum, latent_dat,knearest=10){
set.seed(1234)
suppressMessages(require(FNN))
suppressMessages(require(igraph))
suppressMessages(require(bluster))
 
  PCvalues = latent_dat
info.spatial = as.data.frame(t(PCvalues))
colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))

  knn.norm = get.knn(as.matrix(t(PCvalues)), k = knearest)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index), 
  k=10), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = simplify(nw.norm)
  lc.norm = cluster_louvain(nw.norm)
  merged <- mergeCommunities(nw.norm, lc.norm$membership, number=clusternum)
  clusterlabel = as.character(as.integer(as.factor(paste0("cluster",merged))))
return("cluster_label"=clusterlabel)
}


refine_v2 = function(spotlist, pred_cluster, dist, shape="square"){
  # refined_pred = c(NA,length(pred_cluster))
  #pred = data.frame("spot"=spotlist, "cluster"=pred_cluster)
  dis_df = as.matrix(dist)
  if(shape=="square"){
    num_obs = 4
  }else if(shape == "hexagon"){
    num_nbs = 6
  }else{
    print("Shape is not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
  }
  refined_pred = pred_cluster
  for(i in 1:length(pred_cluster)){
    nearby_spot_ind = order(dist[i,])[1:(num_obs+1)]
    labels_nearby_spot_ind = refined_pred[nearby_spot_ind]  # use updated cluster
    spot_of_interest = refined_pred[i]
    labels_table = table(labels_nearby_spot_ind)
    if( labels_table[spot_of_interest]<num_obs/2 & max(labels_table)>num_obs/2){
      refined_pred[i] = names(labels_table)[which.max(labels_table)]
    }else{
      refined_pred[i] = spot_of_interest
    }

  }

  return(refined_pred)

}


library(tidyverse)
library(mclust)
library(CVXR)
permLabel <- function(labels, perm, zeroidx = FALSE)
{
  if(zeroidx){
    labels <- labels + 1
  }
  K <- length(perm)
  idx <- list()
  for(i in 1:K)
  {
    idx[[i]] <- labels == i
  }
  for(i in 1:K)
  {
    labels[idx[[i]]] <- perm[i]
  }
  labels
}

matchLabel <- function(est_labels, true_labels)
{
  C <- length(unique(true_labels))
  est_labels <- factor(est_labels, 1:C)
  A <- matrix(table(est_labels, true_labels), C, C)
  P <- Variable(C, C, boolean = T) # row permutation matrix
  Y <- Variable(C, C)
  problem <- suppressMessages(
    Problem(Maximize(matrix_trace(Y)), list(Y == P %*% A, 
      sum_entries(P, axis = 1) == 1, sum_entries(P, axis = 2) == 1))
  ) 
  result <- suppressMessages(solve(problem))
  P <- result$getValue(P)
  perm <- apply(P, MARGIN = 2, function(x) which(x==1))
  est_labels <- permLabel(est_labels, perm)
}


simu = function(
  location,
  label,
  init_params,
  scenario,
  J,
  batch_facLoc,
  de_prop,
  de_facLoc,
  de_facScale,
  sim_seed,
  debug = FALSE
  )
{ 

  N <- nrow(location)
  set.seed(sim_seed)
  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  params <- setParams(
    init_params,
    batchCells = rep(3 * N, 1),
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    nGenes = J,
    group.prob = c(0.25, 0.25, 0.25,0.25),
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)

  sim_groups <- splatSimulate(
    params = params,
    method = "groups",
    verbose = FALSE)

  # remove cells having no expressed genes
    idx_zerosize <- apply(counts(sim_groups), MARGIN = 2, sum) == 0
    sim_groups <- sim_groups[, !idx_zerosize]


   if(scenario == 1){
      prop <- c(0.7, 0.1, 0.1,0.1)
    }else if(scenario == 2){
      prop <- c(0.6, 0.3, 0.05,0.05)
    }else if(scenario == 3){
      prop <- c(0.45, 0.45, 0.05,0.05)
    }else if(scenario == 4){
      prop <- c(0.35, 0.30, 0.30,0.05)
    }

  # 3.generate cell types
  print("generate cell types")

    ztrue <- label
    c_simu <- rep(NA, length(ztrue))
    if(scenario == 1){ # scenario 1: 4 cell types
      print("scenario == 1")
      c_simu <- ztrue # cell type assignment same as region assignment
    }else if(scenario != 1){
      print("scenario != 1")
      for(z in unique(label)){
        zi_idx <- ztrue == z # zi_idx is index for region z
        c_simu[zi_idx] <- sample(map_z2c(z), sum(zi_idx), prob = prop, replace = T)
      }
    }

    # 4.assign count data
    groups <- as.data.frame(colData(sim_groups)) %>% filter(Batch == "Batch" %&% 1)
    sim_cnt <- array(NA, c(J, N))
    for(c in 1:4){
      c_size <- sum(c_simu == c)  # true number of cells in cell type c
      c_cells <- groups$Cell[grepl(c, groups$Group)] # generated cells in group c
      cells_select <- sample(as.character(c_cells), c_size, replace = F)
      # sample same number of group c cells in real data from generated cells
      sim_cnt[, c_simu == c] <- as.matrix(counts(sim_groups)[, cells_select])
      # for positions of original cell type c cells, assign generated counts of group c
    }
    colnames(sim_cnt) <- "Cell" %&% 1:N
    rownames(sim_cnt) <- "Gene" %&% 1:J

  return(list(sim_cnt, c_simu, sim_seed))
}



make_grid = function(square_size, location){
  x=location[,1]
  y=location[,2]
  max_x = max(x)
  min_x = min(x)
  max_y = max(y)
  min_y = min(y)
  grid_num_x = round((max_x-min_x)/square_size)+1
  grid_num_y = round((max_y-min_y)/square_size)+1



  grid = list()
  grid_ind = 0
  for(grid_x_id in 1:grid_num_x){
    for(grid_y_id in 1:grid_num_y){
      grid_ind = grid_ind + 1
      x_left = min_x + (grid_x_id-1)*square_size
      x_right = min_x + (grid_x_id)*square_size
      y_low = min_y + (grid_y_id-1)*square_size
      y_up = min_y + (grid_y_id)*square_size
      cell_ids = which(x>=x_left & x<=x_right & y>=y_low & y<=y_up)
      center_x = (x_left+x_right)/2
      center_y = (y_low+y_up)/2
      grid[[grid_ind]] = list()
      grid[[grid_ind]]$cell_ids = cell_ids
      grid[[grid_ind]]$cell_num = length(cell_ids)
      grid[[grid_ind]]$center_x = center_x
      grid[[grid_ind]]$center_y = center_y
      grid[[grid_ind]]$pseudo_x = grid_x_id
      grid[[grid_ind]]$pseudo_y = grid_y_id
    }
  }
  return(grid)
}


make_spot = function(grid_obj, count_mat,celltypes, celltruths){
  gene_num = dim(count_mat)[1]
  spot_num = length(grid_obj)
  count_spot = matrix(0,gene_num,spot_num)
  location_spot = matrix(NA,spot_num,2)
  pseudo_location_spot = matrix(NA,spot_num,2)

  celltype = c()
  subspottruth = c()
  for(spot in 1:spot_num){
    if(length(grid_obj[[spot]]$cell_ids)==1){
      count_spot[,spot] = count_mat[,grid_obj[[spot]]$cell_ids]
      location_spot[spot,1] = grid_obj[[spot]]$center_x
      location_spot[spot,2] = grid_obj[[spot]]$center_y
      pseudo_location_spot[spot,1] = grid_obj[[spot]]$pseudo_x
      pseudo_location_spot[spot,2] = grid_obj[[spot]]$pseudo_y
      celltype[spot] = names(table(celltypes[grid_obj[[spot]]$cell_ids]))[which.max(table(celltypes[grid_obj[[spot]]$cell_ids]))]
      subspottruth[spot] = names(table(celltruths[grid_obj[[spot]]$cell_ids]))[which.max(table(celltruths[grid_obj[[spot]]$cell_ids]))]

    }else if(length(grid_obj[[spot]]$cell_ids)>=2){
      count_spot[,spot] = rowSums(count_mat[,grid_obj[[spot]]$cell_ids])
      location_spot[spot,1] = grid_obj[[spot]]$center_x
      location_spot[spot,2] = grid_obj[[spot]]$center_y
      pseudo_location_spot[spot,1] = grid_obj[[spot]]$pseudo_x
      pseudo_location_spot[spot,2] = grid_obj[[spot]]$pseudo_y
      celltype[spot] = names(table(celltypes[grid_obj[[spot]]$cell_ids]))[which.max(table(celltypes[grid_obj[[spot]]$cell_ids]))]
      subspottruth[spot] = names(table(celltruths[grid_obj[[spot]]$cell_ids]))[which.max(table(celltruths[grid_obj[[spot]]$cell_ids]))]
    }else{
      count_spot[,spot] = 0
      location_spot[spot,1] = grid_obj[[spot]]$center_x
      location_spot[spot,2] = grid_obj[[spot]]$center_y
      pseudo_location_spot[spot,1] = grid_obj[[spot]]$pseudo_x
      pseudo_location_spot[spot,2] = grid_obj[[spot]]$pseudo_y
      celltype[spot] = "empty"
      subspottruth[spot] = "empty"

    }
  }


  return(list("count_spot"=count_spot,"location_spot"=location_spot,"pseudo_location_spot"=pseudo_location_spot,"subspottruth"=subspottruth,"celltype"=celltype))
}


make_spot_from_subspot = function(count_location_subspot){
  count_subspot = count_location_subspot[[1]]
  location_subspot = count_location_subspot[[3]]/3
  truth_subspot = count_location_subspot[[4]]

  spot_x_max = floor(max(location_subspot[,1]))
  spot_y_max = floor(max(location_subspot[,2]))

  spot_num = spot_x_max * spot_y_max
  used_subspots = c()
  truth_spot = c()
  count_spot = matrix(0, dim(count_subspot)[1], spot_num)
  location_spot = matrix(NA, spot_num, 2)
  count = 0
  for(spot_id_x in 1:spot_x_max){
    for(spot_id_y in 1:spot_y_max){
      count = count + 1
      subspot_id = which(location_subspot[,1]>=spot_id_x-0.5 & location_subspot[,1]<=spot_id_x+0.5 & location_subspot[,2]>=spot_id_y-0.5 & location_subspot[,2]<=spot_id_y+0.5)
      used_subspots = c(used_subspots, subspot_id)
      count_spot[, count] = rowSums(count_subspot[,subspot_id])
      location_spot[count,1] = mean(location_subspot[subspot_id,1])
      location_spot[count,2] = mean(location_subspot[subspot_id,2])
      truth_spot[count] = names(table(count_location_subspot$subspottruth[subspot_id]))
    }
  }

  return(list("count_spot"=count_spot, "location_spot" = location_spot, "truth_spot"=truth_spot,"used_subspots" = used_subspots))

}



```

code_utility_python.py
```python
import anndata as ad
import SpaGCN as spg
import pandas as pd
import numpy as np
import scanpy as sc
from scanpy import read_10x_h5
import random
import torch
import cv2
import os

def run_SpaGCN_py(sim_cnt, info, R, p = 0.5):
  adata = ad.AnnData(sim_cnt, obs = info)
  # calculate adjacency matrix
  x_array = adata.obs["x"].tolist()
  y_array = adata.obs["y"].tolist()
  adj = spg.calculate_adj_matrix(x = x_array, y = y_array, histology = False)
  
  # Expression data preprocessing
  adata.var_names_make_unique()
  spg.prefilter_genes(adata)
  spg.prefilter_specialgenes(adata)
  sc.pp.normalize_per_cell(adata)
  sc.pp.log1p(adata)

  # set hyper-parameters
  l = spg.search_l(p, adj, start = 0.01, end = 1000, tol = 0.01, max_run = 100)
  n_clusters = R
  r_seed = t_seed = n_seed = 0
  res = spg.search_res(adata, adj, l, n_clusters, 
    start = 0.7, step = 0.1, tol = 5e-3, lr = 0.05, 
    max_epochs = 20, r_seed = r_seed, t_seed = t_seed, 
    n_seed = n_seed)

  # run spaGCN
  clf = spg.SpaGCN()
  clf.set_l(l)
  random.seed(r_seed)
  torch.manual_seed(t_seed)
  np.random.seed(n_seed)
  clf.train(adata, adj, init_spa = True, init = "louvain", res = res, 
    tol = 5e-3, lr = 0.05, max_epochs = 200)
  y_pred, prob = clf.predict()
  adata.obs["pred"]= y_pred
  adata.obs["pred"] = adata.obs["pred"].astype('category')

  # cluster refinement
  adj_2d = spg.calculate_adj_matrix(x = x_array, y = y_array, histology = False)
  refined_pred = spg.refine(sample_id = adata.obs.index.tolist(), 
    pred = adata.obs["pred"].tolist(), dis = adj_2d, shape = "hexagon")
  adata.obs["refined_pred"] = refined_pred
  adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')

  return(adata.obs)
```


## Simulate data, 10000 cells

```R
args <- as.numeric(commandArgs(TRUE))
i = args[1] # scenario
j = args[2] # repeat
print(i)

  library(tidyverse)
  library(Giotto)
  library(scater)
  library(Seurat)
  library(mclust)
  library(SC3)
  library(BayesSpace)
  library(gtools)
  library(splatter)
  library(reticulate)
  library(mclust )
  library(igraph)
  library(assertthat)
  library(SpatialPCA)
  library(ggplot2)
  
source("code_utility.R")  
# can be found in the above sections
load("LIBDsubsample.RData")
load("init_params_LIBD.RData")
# These two R object can be downloaded from https://drive.google.com/drive/folders/18rwQjB3-g86A-M9xYPPJlHz60bfMABdE?usp=sharing.

res = simu(location=subsample[,2:3],label = subsample$label,init_params,
    scenario=i,J=5000, batch_facLoc=0, de_prop=0.5, de_facLoc=0.5, de_facScale=0.5,sim_seed=j, debug = FALSE)

save(res, file = paste0("LIBD_res_scenairo_",i,"_rep_",j,".RData"))


count_mat = res[[1]]
truth = subsample$label
location=as.matrix(subsample[,2:3])
# truth = LIBDsimu_pseudo$info$z
location=as.matrix(location)
celltypes = res[[2]]
count_spot = count_mat
location_spot = location
rownames(location_spot) = colnames(count_spot) = paste0("spot",1:dim(count_spot)[2])
rownames(count_spot) = paste0("gene",1:dim(count_spot)[1])
dim(count_spot)


###############
# Run SpatialPCA
###############

LIBDsimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
LIBDsimu = SpatialPCA_buildKernel(LIBDsimu, kerneltype="gaussian", bandwidthtype="Silverman")
LIBDsimu = SpatialPCA_EstimateLoading(LIBDsimu,fast=TRUE,SpatialPCnum=20)
LIBDsimu = SpatialPCA_SpatialPCs(LIBDsimu, fast=TRUE)


# Collect results
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = LIBDsimu@SpatialPCs
SpatialPCA_result$normalized_expr  = LIBDsimu@normalized_expr
SpatialPCA_result$location = LIBDsimu@location
pred_cluster= louvain_clustering_new(4,SpatialPCA_result$SpatialPCs,500 )
SpatialPCA_result$clusterlabel = pred_cluster
SpatialPCA_result$truth = truth[match(rownames(LIBDsimu@location),rownames(location_spot))]
SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$truth)
SpatialPCA_result$NMI = compare(as.factor(SpatialPCA_result$clusterlabel),as.factor(SpatialPCA_result$truth), method = "nmi")
SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)

save(SpatialPCA_result, file = paste0("celllevel_SpatialPCA_spatialgene_result_scenario_",i,"_rep_",j,".RData"))

###############
# BayesSpace
###############

print("run BayesSpace, ST, hvg genes")

# filter out spots with 0 counts
ind_keep=which(colSums(count_spot) > 0)
location_spot_bayesSpace=location_spot[ind_keep,]
count_spot_BayesSpace=count_spot[,ind_keep]
colnames(location_spot_bayesSpace) <- c("row", "col")

sce_LIBDsimu_ST <- SingleCellExperiment(assays = list(counts = count_spot_BayesSpace), colData = location_spot_bayesSpace)
sce_LIBDsimu_ST <- spatialPreprocess(sce_LIBDsimu_ST, platform="ST",n.PCs = 15, n.HVGs = 2000, log.normalize = T)
sce_LIBDsimu_ST <- spatialCluster(sce_LIBDsimu_ST, q=4, d=15, platform='ST',nrep=10000, gamma=3, save.chain=FALSE) 
sce_labels=sce_LIBDsimu_ST$spatial.cluster

BayesSpace_ST_result = list()
BayesSpace_ST_result$sce_LIBDsimu_ST = sce_LIBDsimu_ST
BayesSpace_ST_result$clusterlabel = sce_labels
BayesSpace_ST_result$location = location_spot_bayesSpace
BayesSpace_ST_result$truth = truth[ind_keep]
BayesSpace_ST_result$ARI = adjustedRandIndex(BayesSpace_ST_result$clusterlabel,BayesSpace_ST_result$truth)
BayesSpace_ST_result$NMI = compare(BayesSpace_ST_result$clusterlabel,BayesSpace_ST_result$truth, method = "nmi")
BayesSpace_ST_result$CHAOS = fx_CHAOS(BayesSpace_ST_result$clusterlabel, BayesSpace_ST_result$location)
BayesSpace_ST_result$PAS = fx_PAS(BayesSpace_ST_result$clusterlabel, BayesSpace_ST_result$location)

print(BayesSpace_ST_result$ARI)

save(BayesSpace_ST_result,  file = paste0("celllevel_BayesSpace_hvggene_result_scenario_",i,"_rep_",j,".RData"))

###############
# prepare SpaGCN data
# and then run in python
###############

print("SpaGCN")
# SpaGCN uses all genes

library(DropletUtils)
write10xCounts(
  paste0("celllevel_SpaGCN_scenairo_",i,"_rep_",j,"_count_allgene.h5"),
  count_spot,
  barcodes = colnames(count_spot),
  gene.id = rownames(count_spot),
  gene.symbol = rownames(count_spot),
  gene.type = "Gene Expression",
  overwrite = FALSE,
  type = c( "HDF5"),
  genome = "unknown",
  #version = c("2", "3"),
  #chemistry = "Single Cell 3' v3",
  original.gem.groups = 1L,
  library.ids = "custom"
)



```



## Simulate data, spot level

```R

args <- as.numeric(commandArgs(TRUE))
i = args[1] # scenario
j = args[2] # repeatment
k = args[3] # method
num = args[4] # spot level dataset
print(i)
print(j)
print(k)
print(num)

if(num == 1){
	spotnum=5077 # 90um
}else if(num==2){
	spotnum=3602 # 107um
}else{
	spotnum=1948 # 145um
}


setwd("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation")

library(tidyverse)
library(Giotto)
library(scater)
library(SC3)
library(BayesSpace)
library(gtools)
library(splatter)
library(igraph)
library(assertthat)
library(SPARK)
library(Seurat)
library(peakRAM)
library(ggplot2)
library(ggpubr)
library(mclust) # ARI
library(aricode)# NMI
library(SpatialPCA)# NMI
library(lisi) # lisi score
library(reticulate) 
library(FNN)

source("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility.R")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/LIBDsubsample.RData")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/init_params_LIBD.RData")

load(paste0("LIBD_res_scenairo_",i,"_rep_",j,".RData"))


D3=c("#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080",
	"#F94144", "#dda15e", "#F8961E" ,"#bc6c25", "#F9C74F" ,"#90BE6D", "#43AA8B")


# subspot level
count_mat = res[[1]]
truth = subsample$label
location=as.matrix(subsample[,2:3])
location=as.matrix(location)
celltypes = res[[2]]
# first generate subspot level data
if(spotnum == 5077){
	grid_subspot = make_grid(square_size = 4,location)
}else if(spotnum == 3602){
	grid_subspot = make_grid(square_size = 5,location)
}else if(spotnum == 1948){
	grid_subspot = make_grid(square_size = 7,location)
}
count_location_subspot = make_spot(grid_subspot,count_mat,celltypes,subsample$label)
count_subspot = count_location_subspot$count_spot
location_subspot = count_location_subspot$pseudo_location_spot/3
truth_subspot = count_location_subspot$subspottruth
spot_level_data = make_spot_from_subspot(count_location_subspot)
# truth_subspot[spot_level_data$used_subspots]
# spot level
truth=spot_level_data$truth_spot
truth_empty = which(truth=="empty")
count_spot = spot_level_data$count_spot[,-truth_empty]
location_spot = spot_level_data$location_spot[-truth_empty,]
truth=truth[-truth_empty]
rownames(location_spot) = colnames(count_spot) = paste0("spot",1:dim(count_spot)[2])
rownames(count_spot) = paste0("gene",1:dim(count_spot)[1])

save(count_spot,location_spot, truth,celltypes,file=paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/data/res_scenairo_",i,"_rep_",j,"_dataset_",num,"_count_location.RData"))



if(k ==1){

	if(num==1){

		simu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
		simu = SpatialPCA_buildKernel(simu, kerneltype="gaussian", bandwidthtype="Silverman")
		simu = SpatialPCA_EstimateLoading(simu,fast=FALSE,SpatialPCnum=20)
		simu = SpatialPCA_SpatialPCs(simu, fast=FALSE)
		SpatialPCA_result = list()
		SpatialPCA_result$SpatialPCs  = simu@SpatialPCs
		SpatialPCA_result$location = simu@location
		SpatialPCA_result$celltypes = paste0("celltype",celltypes)
		SpatialPCA_result$Truth = truth[match(rownames(simu@location),rownames(location_spot))]
		SpatialPCA_result$clusterlabel = louvain_clustering(4,SpatialPCA_result$SpatialPCs,500) # original:500
		SpatialPCA_result$clusterlabel_refine = refine_cluster_10x(SpatialPCA_result$clusterlabel, SpatialPCA_result$location,shape="square") 
		SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
		SpatialPCA_result$ARI_refine = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$Truth)
		SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
		SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$count_spot = count_spot
		SpatialPCA_result$location_spot = location_spot
		SpatialPCA_result$normalized_expr = simu@normalized_expr
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		save(SpatialPCA_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

		#load(paste0("cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/SpatialPCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
		p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
		p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
		print(ggarrange(p1, p2, ncol = 2, nrow = 1))
		dev.off()

	}else if (num==2){

		simu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
		simu = SpatialPCA_buildKernel(simu, kerneltype="gaussian", bandwidthtype="Silverman")
		simu = SpatialPCA_EstimateLoading(simu,fast=FALSE,SpatialPCnum=20)
		simu = SpatialPCA_SpatialPCs(simu, fast=FALSE)
		SpatialPCA_result = list()
		SpatialPCA_result$SpatialPCs  = simu@SpatialPCs
		SpatialPCA_result$location = simu@location
		SpatialPCA_result$celltypes = paste0("celltype",celltypes)
		SpatialPCA_result$Truth = truth[match(rownames(simu@location),rownames(location_spot))]
		SpatialPCA_result$clusterlabel = louvain_clustering(4,SpatialPCA_result$SpatialPCs,500) # original:500
		SpatialPCA_result$clusterlabel_refine = refine_cluster_10x(SpatialPCA_result$clusterlabel, SpatialPCA_result$location,shape="square") 
		SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
		SpatialPCA_result$ARI_refine = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$Truth)
		SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
		SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$count_spot = count_spot
		SpatialPCA_result$location_spot = location_spot
		SpatialPCA_result$normalized_expr = simu@normalized_expr
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		save(SpatialPCA_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	#load(paste0("cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/SpatialPCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
		p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
		p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
		print(ggarrange(p1, p2, ncol = 2, nrow = 1))
		dev.off()

	}else if (num==3){

		simu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
		simu = SpatialPCA_buildKernel(simu, kerneltype="gaussian", bandwidthtype="Silverman")
		simu = SpatialPCA_EstimateLoading(simu,fast=FALSE,SpatialPCnum=20)
		simu = SpatialPCA_SpatialPCs(simu, fast=FALSE)
		SpatialPCA_result = list()
		SpatialPCA_result$SpatialPCs  = simu@SpatialPCs
		SpatialPCA_result$location = simu@location
		SpatialPCA_result$celltypes = paste0("celltype",celltypes)
		SpatialPCA_result$Truth = truth[match(rownames(simu@location),rownames(location_spot))]
		SpatialPCA_result$clusterlabel = walktrap_clustering(4,SpatialPCA_result$SpatialPCs,100) 
		SpatialPCA_result$clusterlabel_refine = refine_cluster_10x(SpatialPCA_result$clusterlabel, SpatialPCA_result$location,shape="square") 
		SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
		SpatialPCA_result$ARI_refine = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$Truth)
		SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
		SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$count_spot = count_spot
		SpatialPCA_result$location_spot = location_spot
		SpatialPCA_result$normalized_expr = simu@normalized_expr
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		save(SpatialPCA_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	#load(paste0("cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/SpatialPCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
		p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
		p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
		print(ggarrange(p1, p2, ncol = 2, nrow = 1))
		dev.off()

	}



}else if (k==2){ # BayesSpace

	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important

	res = BayesSpace_func(count_in, loc_in, clusternum=clusterNum,nrep=50000,ifLIBD=FALSE,platform="ST")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	save(res, file = paste0("cells_BayesSpace_result_scenario_",i,"_rep_",j,".RData"))

	BayesSpace_result = list()
	BayesSpace_result$clusterlabel = as.character(res$clusterlabel)
	BayesSpace_result$location = SpatialPCA_result$location
	BayesSpace_result$Truth = SpatialPCA_result$Truth
	BayesSpace_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	BayesSpace_result$NMI = NMI(tabb[,1],tabb[,2])
	BayesSpace_result$CHAOS = fx_CHAOS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$PAS = fx_PAS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$LISI = fx_lisi(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$count_spot = SpatialPCA_result$count_spot
	BayesSpace_result$location_spot = SpatialPCA_result$location_spot
	BayesSpace_result$res = res

	save(BayesSpace_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_BayesSpace_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	clusterNum=length(unique(BayesSpace_result$clusterlabel))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/BayesSpace_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=BayesSpace_result$location,clusterlabel=BayesSpace_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=BayesSpace_result$location,clusterlabel=BayesSpace_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("BayesSpace"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==3){ # SpaGCN

	library(reticulate)
	library(tidyverse)
	use_python("/net/mulan/home/shanglu/anaconda3/envs/SpaGCN/bin/python3.7", required=TRUE)
	source_python("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility_python.py")

	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important
	res = run_SpaGCN_py(t(count_in), loc_in , clusterNum) # cell by gene 
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$refined_pred))

	SpaGCN_result = list()
	SpaGCN_result$clusterlabel = res$refined_pred
	SpaGCN_result$location = res[,c("x","y")]
	SpaGCN_result$Truth = SpatialPCA_result$Truth
	SpaGCN_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	SpaGCN_result$NMI = NMI(tabb[,1],tabb[,2])
	SpaGCN_result$CHAOS = fx_CHAOS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$PAS = fx_PAS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$LISI = fx_lisi(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$count_spot = SpatialPCA_result$count_spot
	SpaGCN_result$location_spot = SpatialPCA_result$location_spot

	save(SpaGCN_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpaGCN_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/SpaGCN_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=SpaGCN_result$location,clusterlabel=SpaGCN_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=SpaGCN_result$location,clusterlabel=SpaGCN_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpaGCN"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()


}else if(k==4){ # HMRF
	
	library(reticulate)
	library(tidyverse)
	library(peakRAM)	
	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important

	path_out = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/HMRF/scenario_",i,"_rep_",j,"_dataset_",num)
	res = HMRF_func(count_in=t(count_in), location_in=loc_in, clusternum=clusterNum,path_out=path_out,betas = c(0,2,6)) # betas = [start, step, number]
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res[,6]))

	HMRF_result = list()
	HMRF_result$clusterlabel = as.character(res[,6])
	HMRF_result$location = SpatialPCA_result$location
	HMRF_result$Truth = SpatialPCA_result$Truth
	HMRF_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	HMRF_result$NMI = NMI(tabb[,1],tabb[,2])
	HMRF_result$CHAOS = fx_CHAOS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$PAS = fx_PAS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$LISI = fx_lisi(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$res = res

	save(HMRF_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_HMRF_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/HMRF_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=HMRF_result$location,clusterlabel=HMRF_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=HMRF_result$location,clusterlabel=HMRF_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("HMRF"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==5){ # NMF

	library(peakRAM)

	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])

	res = NMF_cluster_func(count_in=count_in, genelist=rownames(count_in),PCnum=20,cluster_method="louvain",clusternum=clusterNum)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	NMF_result = list()
	NMF_result$clusterlabel = res$clusterlabel
	NMF_result$location = SpatialPCA_result$location
	NMF_result$Truth = SpatialPCA_result$Truth
	NMF_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	NMF_result$NMI = NMI(tabb[,1],tabb[,2])
	NMF_result$CHAOS = fx_CHAOS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PAS = fx_PAS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$LISI = fx_lisi(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PCs = res$PC

	save(NMF_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_NMF_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/NMF_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=NMF_result$location,clusterlabel=NMF_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=NMF_result$location,clusterlabel=NMF_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("NMF"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==6){ # PCA

	
	load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id = match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$count_spot))
	loc_in = as.data.frame(SpatialPCA_result$location_spot[match_id,])
	count_in = as.matrix(SpatialPCA_result$count_spot[,match_id])

	res = PCA_cluster_func(expr=SpatialPCA_result$normalized_expr,PCnum=20,cluster_method="louvain",clusternum=clusterNum)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	PCA_result = list()
	PCA_result$clusterlabel = res$clusterlabel
	PCA_result$location = SpatialPCA_result$location
	PCA_result$Truth = SpatialPCA_result$Truth
	PCA_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	PCA_result$NMI = NMI(tabb[,1],tabb[,2])
	PCA_result$CHAOS = fx_CHAOS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PAS = fx_PAS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$LISI = fx_lisi(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PCs = res$PC
	PCA_result$count_spot = SpatialPCA_result$count_spot
	PCA_result$location_spot = SpatialPCA_result$location_spot


	save(PCA_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/result/cells_PCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/simulation/figure/PCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=PCA_result$location,clusterlabel=PCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=PCA_result$location,clusterlabel=PCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("PCA"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()


}



```

## Simulate data, spot level with stripe pattern
```R
args <- as.numeric(commandArgs(TRUE))
i = args[1] # scenario
j = args[2] # repeatment
k = args[3] # method
print(i)
print(j)
print(k)

setwd("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2")

library(tidyverse)
library(Giotto)
library(scater)
library(SC3)
library(BayesSpace)
library(gtools)
library(splatter)
library(igraph)
library(assertthat)
library(SPARK)
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


D3=c("#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080",
	"#F94144", "#dda15e", "#F8961E" ,"#bc6c25", "#F9C74F" ,"#90BE6D", "#43AA8B",
 "#4D908E", "#577590", "#277DA1", "#AEC7E8" ,"#FFBB78" ,"#98DF8A", "#FF9896",
"#C5B0D5" ,"#C49C94", "#F7B6D2", "#C7C7C7" ,"#DBDB8D" ,"#9EDAE5", "#393B79",
"#637939", "#8C6D31", "#843C39", "#7B4173" ,"#5254A3" ,"#8CA252", "#BD9E39",
"#AD494A", "#A55194", "#6B6ECF", "#B5CF6B", "#E7BA52" ,"#D6616B", "#CE6DBD",
"#9C9EDE", "#CEDB9C" ,"#E7CB94", "#E7969C", "#DE9ED6" ,"#3182BD", "#E6550D",
"#31A354", "#756BB1" ,"#636363", "#6BAED6" ,"#FD8D3C" ,"#74C476", "#9E9AC8",
"#969696", "#9ECAE1" ,"#FDAE6B", "#A1D99B" ,"#BCBDDC" ,"#BDBDBD", "#C6DBEF",
"#FDD0A2" ,"#C7E9C0" ,"#DADAEB", "#D9D9D9")

load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/LIBDsubsample.RData")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/init_params_LIBD.RData")

# 1,2,1,2,1,2 regions
# still 10000 cells
a = t(combn(1:1200,2))
b=a
b[,1] = a[,2]
b[,2] = a[,1]
c=rbind(a,b)
# x: 1:1000
# y: 1:200
a1=c[which(c[,2]<=200),]
a1 = as.data.frame(a1)
a1$label = 1
a1$label[which(a1[,1]>200 & a1[,1]<=400)] =2
a1$label[which(a1[,1]>600 & a1[,1]<=800)] =2
a1$label[which(a1[,1]>1000 & a1[,1]<=1200)] =2
a1$label = as.character(a1$label)

set.seed(06112022)
subsample = a1[sample(c(1:dim(a1)[1]),10000),]

save(subsample, file = "subsample_stripe.RData")


load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/subsample_stripe.RData")


seed=06112022+j

res = simu_stripe(location=subsample[,1:2],label = subsample[,3],init_params,
    scenario=i,J=5000, batch_facLoc=0, de_prop=0.95, de_facLoc=0.95, de_facScale=0.05,sim_seed=seed, debug = FALSE)

save(res, file=paste0("stripe_res_scenairo_",i,"_rep_",j,".RData"))

# subspot level
count_mat = res[[1]]
truth = subsample[,3]
location=as.matrix(subsample[,1:2])
# truth = stripesimu_pseudo$info$z
location=as.matrix(location)
celltypes = res[[2]]

count_spot = count_mat
location_spot = location
rownames(location_spot) = colnames(count_spot) = paste0("spot",1:dim(count_spot)[2])
rownames(count_spot) = paste0("gene",1:dim(count_spot)[1])
dim(count_spot)

save(count_spot,location_spot, truth,celltypes,file=paste0("stripe_res_scenairo_",i,"_rep_",j,"_count_location.RData"))


print("data saved!")

print("data loading!")

load(paste0("stripe_res_scenairo_",i,"_rep_",j,"_count_location.RData"))

if(k ==1){

stripesimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
stripesimu = SpatialPCA_buildKernel(stripesimu, kerneltype="gaussian", bandwidthtype="Silverman")
stripesimu = SpatialPCA_EstimateLoading(stripesimu,fast=TRUE,SpatialPCnum=20)
stripesimu = SpatialPCA_SpatialPCs(stripesimu, fast=TRUE)

SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = stripesimu@SpatialPCs
SpatialPCA_result$location = stripesimu@location
SpatialPCA_result$clusterlabel = louvain_clustering(2,SpatialPCA_result$SpatialPCs,100) # original:500
SpatialPCA_result$celltypes = paste0("celltype",celltypes)
SpatialPCA_result$Truth = truth[match(rownames(stripesimu@location),rownames(location_spot))]
SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
SpatialPCA_result$count_spot = count_spot
SpatialPCA_result$location_spot = location_spot
SpatialPCA_result$normalized_expr = stripesimu@normalized_expr

clusterNum=length(unique(SpatialPCA_result$clusterlabel))

save(SpatialPCA_result, file = paste0("stripe_pattern_SpatialPCA_result_scenario_",i,"_rep_",j,".RData"))
	
	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/figure/SpatialPCA_scenario_",i,"_rep_",j,".pdf"),width=20,height=5)
	p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()

}else if (k==2){ # BayesSpace


	load(paste0("stripe_pattern_SpatialPCA_result_scenario_",i,"_rep_",j,".RData"))
	count_in = SpatialPCA_result$count_spot
	loc_in=SpatialPCA_result$location_spot
	clusterNum = length(unique(SpatialPCA_result$Truth))
	res = BayesSpace_func(count_in, loc_in, clusternum=clusterNum,nrep=50000,ifLIBD=FALSE)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	save(res, file = paste0("stripe_pattern_BayesSpace_result_scenario_",i,"_rep_",j,".RData"))

	BayesSpace_result = list()
	BayesSpace_result$clusterlabel = as.character(res$clusterlabel)
	BayesSpace_result$location = SpatialPCA_result$location
	BayesSpace_result$Truth = SpatialPCA_result$Truth
	BayesSpace_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	BayesSpace_result$NMI = NMI(tabb[,1],tabb[,2])
	BayesSpace_result$CHAOS = fx_CHAOS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$PAS = fx_PAS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$LISI = fx_lisi(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$count_spot = SpatialPCA_result$count_spot
	BayesSpace_result$location_spot = SpatialPCA_result$location_spot
	BayesSpace_result$res = res

	save(BayesSpace_result, file = paste0("stripe_pattern_BayesSpace_result_scenario_",i,"_rep_",j,".RData"))

	clusterNum=length(unique(BayesSpace_result$clusterlabel))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/figure/BayesSpace_scenario_",i,"_rep_",j,".pdf"),width=20,height=5)
	p1=plot_cluster(location=BayesSpace_result$location,clusterlabel=BayesSpace_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=BayesSpace_result$location,clusterlabel=BayesSpace_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("BayesSpace"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==3){ # SpaGCN

	library(reticulate)
	library(tidyverse)
	use_python("/net/mulan/home/shanglu/anaconda3/envs/SpaGCN/bin/python3.7", required=TRUE)
	source_python("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility_python.py")

	load(paste0("stripe_pattern_SpatialPCA_result_scenario_",i,"_rep_",j,".RData"))
	count_in = SpatialPCA_result$count_spot
	loc_in=as.data.frame( SpatialPCA_result$location_spot)
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important
	res = run_SpaGCN_py(t(count_in), loc_in , clusterNum) # cell by gene 
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$refined_pred))

	SpaGCN_result = list()
	SpaGCN_result$clusterlabel = res$refined_pred
	SpaGCN_result$location = res[,c("x","y")]
	SpaGCN_result$Truth = SpatialPCA_result$Truth
	SpaGCN_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	SpaGCN_result$NMI = NMI(tabb[,1],tabb[,2])
	SpaGCN_result$CHAOS = fx_CHAOS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$PAS = fx_PAS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$LISI = fx_lisi(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$count_spot = SpatialPCA_result$count_spot
	SpaGCN_result$location_spot = SpatialPCA_result$location_spot

	save(SpaGCN_result, file = paste0("stripe_pattern_SpaGCN_result_scenario_",i,"_rep_",j,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/figure/SpaGCN_scenario_",i,"_rep_",j,".pdf"),width=20,height=5)
	p1=plot_cluster(location=SpaGCN_result$location,clusterlabel=SpaGCN_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=SpaGCN_result$location,clusterlabel=SpaGCN_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpaGCN"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()




}else if(k==4){ # HMRF
	
	library(reticulate)
	library(tidyverse)
	library(peakRAM)

	load(paste0("stripe_pattern_SpatialPCA_result_scenario_",i,"_rep_",j,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	loc_in=as.data.frame( SpatialPCA_result$location_spot)
	count_in =  SpatialPCA_result$count_spot

	path_out = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/HMRF/scenario_",i,"_rep_",j)
	res = HMRF_func(count_in=t(count_in), location_in=loc_in, clusternum=clusterNum,path_out=path_out,betas = c(0,2,6)) # betas = [start, step, number]
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res[,6]))

	HMRF_result = list()
	HMRF_result$clusterlabel = as.character(res[,6])
	HMRF_result$location = SpatialPCA_result$location
	HMRF_result$Truth = SpatialPCA_result$Truth
	HMRF_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	HMRF_result$NMI = NMI(tabb[,1],tabb[,2])
	HMRF_result$CHAOS = fx_CHAOS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$PAS = fx_PAS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$LISI = fx_lisi(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$res = res

	save(HMRF_result, file = paste0("stripe_pattern_HMRF_result_scenario_",i,"_rep_",j,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/figure/HMRF_scenario_",i,"_rep_",j,".pdf"),width=20,height=5)
	p1=plot_cluster(location=HMRF_result$location,clusterlabel=HMRF_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=HMRF_result$location,clusterlabel=HMRF_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("HMRF"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==5){ # NMF

	library(peakRAM)

	load(paste0("stripe_pattern_SpatialPCA_result_scenario_",i,"_rep_",j,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	loc_in=as.data.frame( SpatialPCA_result$location_spot)
	count_in =  SpatialPCA_result$count_spot

	res = NMF_cluster_func(count_in=count_in, genelist=rownames(count_in),PCnum=20,cluster_method="louvain",clusternum=clusterNum)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	NMF_result = list()
	NMF_result$clusterlabel = res$clusterlabel
	NMF_result$location = SpatialPCA_result$location
	NMF_result$Truth = SpatialPCA_result$Truth
	NMF_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	NMF_result$NMI = NMI(tabb[,1],tabb[,2])
	NMF_result$CHAOS = fx_CHAOS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PAS = fx_PAS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$LISI = fx_lisi(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PCs = res$PC

save(NMF_result, file = paste0("stripe_pattern_NMF_result_scenario_",i,"_rep_",j,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/figure/NMF_scenario_",i,"_rep_",j,".pdf"),width=20,height=5)
	p1=plot_cluster(location=NMF_result$location,clusterlabel=NMF_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=NMF_result$location,clusterlabel=NMF_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("NMF"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==6){ # PCA

	
	load(paste0("stripe_pattern_SpatialPCA_result_scenario_",i,"_rep_",j,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id = match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$count_spot))
	loc_in = as.data.frame(SpatialPCA_result$location_spot[match_id,])
	count_in = as.matrix(SpatialPCA_result$count_spot[,match_id])

	res = PCA_cluster_func(expr=SpatialPCA_result$normalized_expr,PCnum=20,cluster_method="louvain",clusternum=clusterNum)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	PCA_result = list()
	PCA_result$clusterlabel = res$clusterlabel
	PCA_result$location = SpatialPCA_result$location
	PCA_result$Truth = SpatialPCA_result$Truth
	PCA_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	PCA_result$NMI = NMI(tabb[,1],tabb[,2])
	PCA_result$CHAOS = fx_CHAOS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PAS = fx_PAS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$LISI = fx_lisi(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PCs = res$PC
	PCA_result$count_spot = SpatialPCA_result$count_spot
	PCA_result$location_spot = SpatialPCA_result$location_spot


	save(PCA_result, file = paste0("stripe_pattern_PCA_result_scenario_",i,"_rep_",j,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q2/figure/PCA_scenario_",i,"_rep_",j,".pdf"),width=20,height=5)
	p1=plot_cluster(location=PCA_result$location,clusterlabel=PCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=PCA_result$location,clusterlabel=PCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("PCA"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()


}

```


## Simulate data, spot level with arbitrary spatial correlation
```R

args <- as.numeric(commandArgs(TRUE))
i = args[1] # scenario
j = args[2] # repeatment
k = args[3] # method
num = args[4] # spot level dataset
print(i)
print(j)
print(k)
print(num)

if(num == 1){
	spotnum=5077
}else if(num==2){
	spotnum=3602
}else{
	spotnum=1948
}

# cd /net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3

#--------------------------------------------------------------------


setwd("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29")



library(tidyverse)
library(Giotto)
library(scater)
library(SC3)
library(BayesSpace)
library(gtools)
library(splatter)
library(igraph)
library(assertthat)
library(SPARK)
library(Seurat)
library(peakRAM)
library(ggplot2)
library(ggpubr)
library(mclust) # ARI
library(aricode)# NMI
library(SpatialPCA)# NMI
library(lisi) # lisi score
library(reticulate) 
library(FNN)

source("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility.R")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/LIBDsubsample.RData")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/init_params_LIBD.RData")

load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/results_5077_spark/LIBD_res_scenairo_",i,"_rep_",j,".RData"))


D3=c("#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080",
	"#F94144", "#dda15e", "#F8961E" ,"#bc6c25", "#F9C74F" ,"#90BE6D", "#43AA8B")


# subspot level
count_mat = res[[1]]
truth = subsample$label
location=as.matrix(subsample[,2:3])
knearest = 4
indNN = FNN::get.knn(location,k=knearest)[[1]] 
# randomly select 2500 cells to fuse with its nearest 4 cells
set.seed(j+07132022)
tofuse_cells = sample(1:10000,2500)
count_mat_use = count_mat
for(cell_id in tofuse_cells){
	set.seed(cell_id+07292022)
	sample_tmp = sample(1:10,knearest)
	prop = matrix(sample_tmp/sum(sample_tmp),1,knearest) 
	# for each nearest cell, randomly split the current cell to 4 components, and fuse to the 4 nearest cells
	count_mat_use[,indNN[cell_id,]] = count_mat[,indNN[cell_id,]]+ count_mat[,cell_id] %*% prop
	count_mat_use[,cell_id]=0
}
count_mat= count_mat_use
location=as.matrix(location)
celltypes = res[[2]]
# first generate subspot level data
if(spotnum == 5077){
	grid_subspot = make_grid(square_size = 4,location)
}else if(spotnum == 3602){
	grid_subspot = make_grid(square_size = 5,location)
}else if(spotnum == 1948){
	grid_subspot = make_grid(square_size = 7,location)
}
count_location_subspot = make_spot(grid_subspot,count_mat,celltypes,subsample$label)
count_subspot = count_location_subspot$count_spot
location_subspot = count_location_subspot$pseudo_location_spot/3
truth_subspot = count_location_subspot$subspottruth
spot_level_data = make_spot_from_subspot(count_location_subspot)
# truth_subspot[spot_level_data$used_subspots]
# spot level
truth=spot_level_data$truth_spot
truth_empty = which(truth=="empty")
count_spot = spot_level_data$count_spot[,-truth_empty]
location_spot = spot_level_data$location_spot[-truth_empty,]
truth=truth[-truth_empty]
rownames(location_spot) = colnames(count_spot) = paste0("spot",1:dim(count_spot)[2])
rownames(count_spot) = paste0("gene",1:dim(count_spot)[1])

save(count_spot,location_spot, truth,celltypes,file=paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/fuse_res_scenairo_",i,"_rep_",j,"_dataset_",num,"_count_location.RData"))



load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/fuse_res_scenairo_",i,"_rep_",j,"_dataset_",num,"_count_location.RData"))


if(k ==1){

	if(num==1){

		fusesimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
		fusesimu = SpatialPCA_buildKernel(fusesimu, kerneltype="gaussian", bandwidthtype="Silverman")
		fusesimu = SpatialPCA_EstimateLoading(fusesimu,fast=FALSE,SpatialPCnum=20)
		fusesimu = SpatialPCA_SpatialPCs(fusesimu, fast=FALSE)
		SpatialPCA_result = list()
		SpatialPCA_result$SpatialPCs  = fusesimu@SpatialPCs
		SpatialPCA_result$location = fusesimu@location
		SpatialPCA_result$celltypes = paste0("celltype",celltypes)
		SpatialPCA_result$Truth = truth[match(rownames(fusesimu@location),rownames(location_spot))]
		SpatialPCA_result$clusterlabel = louvain_clustering(4,SpatialPCA_result$SpatialPCs,500) # original:500
		SpatialPCA_result$clusterlabel_refine = refine_cluster_10x(SpatialPCA_result$clusterlabel, SpatialPCA_result$location,shape="square") 
		SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
		SpatialPCA_result$ARI_refine = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$Truth)
		SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
		SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$count_spot = count_spot
		SpatialPCA_result$location_spot = location_spot
		SpatialPCA_result$normalized_expr = fusesimu@normalized_expr
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		save(SpatialPCA_result, file = paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

		#load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/SpatialPCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
		p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
		p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
		print(ggarrange(p1, p2, ncol = 2, nrow = 1))
		dev.off()

	}else if (num==2){

		fusesimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
		fusesimu = SpatialPCA_buildKernel(fusesimu, kerneltype="gaussian", bandwidthtype="Silverman")
		fusesimu = SpatialPCA_EstimateLoading(fusesimu,fast=FALSE,SpatialPCnum=20)
		fusesimu = SpatialPCA_SpatialPCs(fusesimu, fast=FALSE)
		SpatialPCA_result = list()
		SpatialPCA_result$SpatialPCs  = fusesimu@SpatialPCs
		SpatialPCA_result$location = fusesimu@location
		SpatialPCA_result$celltypes = paste0("celltype",celltypes)
		SpatialPCA_result$Truth = truth[match(rownames(fusesimu@location),rownames(location_spot))]
		SpatialPCA_result$clusterlabel = louvain_clustering(4,SpatialPCA_result$SpatialPCs,500) # original:500
		SpatialPCA_result$clusterlabel_refine = refine_cluster_10x(SpatialPCA_result$clusterlabel, SpatialPCA_result$location,shape="square") 
		SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
		SpatialPCA_result$ARI_refine = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$Truth)
		SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
		SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$count_spot = count_spot
		SpatialPCA_result$location_spot = location_spot
		SpatialPCA_result$normalized_expr = fusesimu@normalized_expr
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		save(SpatialPCA_result, file = paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	#load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/SpatialPCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
		p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
		p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
		print(ggarrange(p1, p2, ncol = 2, nrow = 1))
		dev.off()

	}else if (num==3){

		fusesimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
		fusesimu = SpatialPCA_buildKernel(fusesimu, kerneltype="gaussian", bandwidthtype="Silverman")
		fusesimu = SpatialPCA_EstimateLoading(fusesimu,fast=FALSE,SpatialPCnum=20)
		fusesimu = SpatialPCA_SpatialPCs(fusesimu, fast=FALSE)
		SpatialPCA_result = list()
		SpatialPCA_result$SpatialPCs  = fusesimu@SpatialPCs
		SpatialPCA_result$location = fusesimu@location
		SpatialPCA_result$celltypes = paste0("celltype",celltypes)
		SpatialPCA_result$Truth = truth[match(rownames(fusesimu@location),rownames(location_spot))]
		SpatialPCA_result$clusterlabel = walktrap_clustering(4,SpatialPCA_result$SpatialPCs,100) 
		SpatialPCA_result$clusterlabel_refine = refine_cluster_10x(SpatialPCA_result$clusterlabel, SpatialPCA_result$location,shape="square") 
		SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$Truth)
		SpatialPCA_result$ARI_refine = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine,SpatialPCA_result$Truth)
		SpatialPCA_result$NMI = NMI(as.character(SpatialPCA_result$clusterlabel),as.character(SpatialPCA_result$Truth))
		SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$LISI = fx_lisi(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
		SpatialPCA_result$count_spot = count_spot
		SpatialPCA_result$location_spot = location_spot
		SpatialPCA_result$normalized_expr = fusesimu@normalized_expr
		clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		save(SpatialPCA_result, file = paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	#load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum=length(unique(SpatialPCA_result$clusterlabel))
		pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/SpatialPCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
		p1=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
		p2=plot_cluster(location=SpatialPCA_result$location,clusterlabel=SpatialPCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=D3[1:clusterNum])
		print(ggarrange(p1, p2, ncol = 2, nrow = 1))
		dev.off()

	}



}else if (k==2){ # BayesSpace

	load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important

	res = BayesSpace_func(count_in, loc_in, clusternum=clusterNum,nrep=50000,ifLIBD=FALSE,platform="ST")
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	save(res, file = paste0("fused_cells_BayesSpace_result_scenario_",i,"_rep_",j,".RData"))

	BayesSpace_result = list()
	BayesSpace_result$clusterlabel = as.character(res$clusterlabel)
	BayesSpace_result$location = SpatialPCA_result$location
	BayesSpace_result$Truth = SpatialPCA_result$Truth
	BayesSpace_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	BayesSpace_result$NMI = NMI(tabb[,1],tabb[,2])
	BayesSpace_result$CHAOS = fx_CHAOS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$PAS = fx_PAS(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$LISI = fx_lisi(res$clusterlabel, SpatialPCA_result$location)
	BayesSpace_result$count_spot = SpatialPCA_result$count_spot
	BayesSpace_result$location_spot = SpatialPCA_result$location_spot
	BayesSpace_result$res = res

	save(BayesSpace_result, file = paste0("fused_cells_BayesSpace_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	clusterNum=length(unique(BayesSpace_result$clusterlabel))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/BayesSpace_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=BayesSpace_result$location,clusterlabel=BayesSpace_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=BayesSpace_result$location,clusterlabel=BayesSpace_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("BayesSpace"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==3){ # SpaGCN

	library(reticulate)
	library(tidyverse)
	use_python("/net/mulan/home/shanglu/anaconda3/envs/SpaGCN/bin/python3.7", required=TRUE)
	source_python("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility_python.py")

	load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important
	res = run_SpaGCN_py(t(count_in), loc_in , clusterNum) # cell by gene 
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$refined_pred))

	SpaGCN_result = list()
	SpaGCN_result$clusterlabel = res$refined_pred
	SpaGCN_result$location = res[,c("x","y")]
	SpaGCN_result$Truth = SpatialPCA_result$Truth
	SpaGCN_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	SpaGCN_result$NMI = NMI(tabb[,1],tabb[,2])
	SpaGCN_result$CHAOS = fx_CHAOS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$PAS = fx_PAS(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$LISI = fx_lisi(res$refined_pred, SpaGCN_result$location)
	SpaGCN_result$count_spot = SpatialPCA_result$count_spot
	SpaGCN_result$location_spot = SpatialPCA_result$location_spot

	save(SpaGCN_result, file = paste0("fused_cells_SpaGCN_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/SpaGCN_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=SpaGCN_result$location,clusterlabel=SpaGCN_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=SpaGCN_result$location,clusterlabel=SpaGCN_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("SpaGCN"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()


}else if(k==4){ # HMRF
	
	library(reticulate)
	library(tidyverse)
	library(peakRAM)	
	load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])
	clusterNum = length(unique(SpatialPCA_result$Truth))
	colnames(loc_in)=c("x","y") # important

	path_out = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/HMRF/scenario_",i,"_rep_",j,"_dataset_",num)
	res = HMRF_func(count_in=t(count_in), location_in=loc_in, clusternum=clusterNum,path_out=path_out,betas = c(0,2,6)) # betas = [start, step, number]
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res[,6]))

	HMRF_result = list()
	HMRF_result$clusterlabel = as.character(res[,6])
	HMRF_result$location = SpatialPCA_result$location
	HMRF_result$Truth = SpatialPCA_result$Truth
	HMRF_result$ARI = adjustedRandIndex(tabb[,1], tabb[,2])
	HMRF_result$NMI = NMI(tabb[,1],tabb[,2])
	HMRF_result$CHAOS = fx_CHAOS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$PAS = fx_PAS(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$LISI = fx_lisi(HMRF_result$clusterlabel, HMRF_result$location)
	HMRF_result$res = res

	save(HMRF_result, file = paste0("fused_cells_HMRF_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/HMRF_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=HMRF_result$location,clusterlabel=HMRF_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=HMRF_result$location,clusterlabel=HMRF_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("HMRF"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==5){ # NMF

	library(peakRAM)

	load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id=match(colnames(SpatialPCA_result$normalized_expr), colnames(SpatialPCA_result$count_spot))
	count_in = SpatialPCA_result$count_spot[,match_id]
	loc_in=as.data.frame( SpatialPCA_result$location_spot[match_id,])

	res = NMF_cluster_func(count_in=count_in, genelist=rownames(count_in),PCnum=20,cluster_method="louvain",clusternum=clusterNum)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	NMF_result = list()
	NMF_result$clusterlabel = res$clusterlabel
	NMF_result$location = SpatialPCA_result$location
	NMF_result$Truth = SpatialPCA_result$Truth
	NMF_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	NMF_result$NMI = NMI(tabb[,1],tabb[,2])
	NMF_result$CHAOS = fx_CHAOS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PAS = fx_PAS(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$LISI = fx_lisi(NMF_result$clusterlabel, SpatialPCA_result$location)
	NMF_result$PCs = res$PC

	save(NMF_result, file = paste0("fused_cells_NMF_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/NMF_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=NMF_result$location,clusterlabel=NMF_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=NMF_result$location,clusterlabel=NMF_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("NMF"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()



}else if(k==6){ # PCA

	
	load(paste0("fused_cells_SpatialPCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))
	clusterNum = length(unique(na.omit(SpatialPCA_result$Truth)))
	match_id = match(colnames(SpatialPCA_result$normalized_expr),colnames(SpatialPCA_result$count_spot))
	loc_in = as.data.frame(SpatialPCA_result$location_spot[match_id,])
	count_in = as.matrix(SpatialPCA_result$count_spot[,match_id])

	res = PCA_cluster_func(expr=SpatialPCA_result$normalized_expr,PCnum=20,cluster_method="louvain",clusternum=clusterNum)
	tabb = na.omit(data.frame("Truth"=SpatialPCA_result$Truth,"clusterlabel"=res$clusterlabel))

	PCA_result = list()
	PCA_result$clusterlabel = res$clusterlabel
	PCA_result$location = SpatialPCA_result$location
	PCA_result$Truth = SpatialPCA_result$Truth
	PCA_result$ARI =  adjustedRandIndex(tabb[,1], tabb[,2])
	PCA_result$NMI = NMI(tabb[,1],tabb[,2])
	PCA_result$CHAOS = fx_CHAOS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PAS = fx_PAS(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$LISI = fx_lisi(PCA_result$clusterlabel, PCA_result$location)
	PCA_result$PCs = res$PC
	PCA_result$count_spot = SpatialPCA_result$count_spot
	PCA_result$location_spot = SpatialPCA_result$location_spot


	save(PCA_result, file = paste0("fused_cells_PCA_result_scenario_",i,"_rep_",j,"_dataset_",num,".RData"))

	pdf(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer2/Q3/July29/figure/PCA_scenario_",i,"_rep_",j,"_dataset_",num,".pdf"),width=10,height=5)
	p1=plot_cluster(location=PCA_result$location,clusterlabel=PCA_result$Truth,pointsize=2,text_size=20 ,title_in=paste0("Truth"),color_in=D3[1:clusterNum])
	p2=plot_cluster(location=PCA_result$location,clusterlabel=PCA_result$clusterlabel,pointsize=2,text_size=20 ,title_in=paste0("PCA"),color_in=D3[1:clusterNum])
	print(ggarrange(p1, p2, ncol = 2, nrow = 1))
	dev.off()


}




```

## High resolution spatial map reconstruction
```R
args <- as.numeric(commandArgs(TRUE))
i = args[1] # scenario
j = args[2] # repeat
print(i)

  library(tidyverse)
  library(Giotto)
  library(scater)
  library(Seurat)
  library(mclust)
  library(SC3)
  library(BayesSpace)
  library(gtools)
  library(splatter)
  library(reticulate)
  library(mclust )
  library(igraph)
  library(assertthat)
  library(SpatialPCA)
  library(ggplot2)
  
 source("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/code_utility.R")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/LIBDsubsample.RData")
load("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/init_params_LIBD.RData")
load(paste0("LIBD_res_scenairo_",i,"_rep_",j,".RData"))


# subspot level
count_mat = res[[1]]
truth = subsample$label
location=as.matrix(subsample[,2:3])
# truth = LIBDsimu_pseudo$info$z
location=as.matrix(location)
celltypes = res[[2]]
# first generate subspot level data
grid_subspot = make_grid(square_size = 7,location)
count_location_subspot = make_spot(grid_subspot,count_mat,celltypes,subsample$label)
count_subspot = count_location_subspot$count_spot
location_subspot = count_location_subspot$pseudo_location_spot/3
truth_subspot = count_location_subspot$subspottruth
spot_level_data = make_spot_from_subspot(count_location_subspot)
truth_subspot[spot_level_data$used_subspots]
# spot level
truth=spot_level_data$truth_spot
truth_empty = which(truth=="empty")
count_spot = spot_level_data$count_spot[,-truth_empty]
location_spot = spot_level_data$location_spot[-truth_empty,]
truth=truth[-truth_empty]
rownames(location_spot) = colnames(count_spot) = paste0("spot",1:dim(count_spot)[2])
rownames(count_spot) = paste0("gene",1:dim(count_spot)[1])

dim(count_spot)

save(count_spot, location_spot, truth, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/data/Data_count_location_truth_scenario_",i,"_rep_",j,".RData"))


#----------------------------------  
# Run SpatialPCA
#----------------------------------  
print("run SpatialPCA")

LIBDsimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
LIBDsimu = SpatialPCA_buildKernel(LIBDsimu, kerneltype="gaussian", bandwidthtype="SJ")
LIBDsimu = SpatialPCA_EstimateLoading(LIBDsimu,fast=FALSE,SpatialPCnum=20)
LIBDsimu = SpatialPCA_SpatialPCs(LIBDsimu, fast=FALSE)
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = LIBDsimu@SpatialPCs
SpatialPCA_result$normalized_expr  = LIBDsimu@normalized_expr
SpatialPCA_result$location = LIBDsimu@location
pred_cluster= walktrap_clustering_new(4,SpatialPCA_result$SpatialPCs,100 )
spotlist=rownames(SpatialPCA_result$location)
dist = as.matrix(dist(SpatialPCA_result$location))
pred_refine = refine_v2(spotlist, pred_cluster, dist,shape="square") 
SpatialPCA_result$pred_cluster = pred_cluster
SpatialPCA_result$clusterlabel = pred_refine
SpatialPCA_result$truth = truth[match(rownames(LIBDsimu@location),rownames(location_spot))]
SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$truth)
SpatialPCA_result$NMI = compare(as.factor(SpatialPCA_result$clusterlabel),as.factor(SpatialPCA_result$truth), method = "nmi")
SpatialPCA_result$CHAOS = fx_chaos(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
SpatialPCA_result$PAS = fx_PAEP(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
print(SpatialPCA_result$ARI )

save(SpatialPCA_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/result/spotlevel_SpatialPCA_spatialgene_result_scenario_",i,"_rep_",j,".RData"))

# load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/reviewer/simulation/LIBDsimu/results_1948/spotlevel_SpatialPCA_spatialgene_result_scenario_",i,"_rep_",j,".RData"))

#----------------------------------  
# run BayesSpace, ST, hvg genes
#----------------------------------  

print("run BayesSpace, ST, hvg genes")
# filter out spots with 0 counts
ind_keep=which(colSums(count_spot) > 0)
location_spot_bayesSpace=location_spot[ind_keep,]
count_spot_BayesSpace=count_spot[,ind_keep]
colnames(location_spot_bayesSpace) <- c("row", "col")

sce_LIBDsimu_ST <- SingleCellExperiment(assays = list(counts = count_spot_BayesSpace), colData = location_spot_bayesSpace)
sce_LIBDsimu_ST <- spatialPreprocess(sce_LIBDsimu_ST, platform="ST",n.PCs = 15, n.HVGs = 2000, log.normalize = T)
sce_LIBDsimu_ST <- spatialCluster(sce_LIBDsimu_ST, q=4, d=15, platform='ST',nrep=10000, gamma=3, save.chain=FALSE) 
sce_labels=sce_LIBDsimu_ST$spatial.cluster

save(sce_LIBDsimu_ST,  file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/result/BayesSpace_sce_LIBDsimu_ST_scenario_",i,"_rep_",j,".RData"))


#-------------------------------------
# run BayesSpace, high resolution, ST
#-------------------------------------

print("BayesSpace, high resolution, ST")
sce_LIBDsimu_ST_enhanced <- spatialEnhance(sce_LIBDsimu_ST, q=4, platform="ST", d=15,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, 
                                    burn.in=1000,
                                    save.chain=FALSE)
sce_LIBDsimu_ST_enhanced=BayesSpace_high_ST_result$sce_LIBDsimu_ST_enhanced
sce_LIBDsimu_ST_enhanced_labels = sce_LIBDsimu_ST_enhanced$spatial.cluster
vertices_ST <- unique(.make_vertices(sce_LIBDsimu_ST_enhanced, sce_LIBDsimu_ST_enhanced_labels, "ST", is.enhanced="TRUE"))
vertices_loc = unique(vertices_ST[,c(1,2,4)])
vertices_location = vertices_loc[,1:2]
BayesSpace_high_ST_result = list()
BayesSpace_high_ST_result$sce_LIBDsimu_ST_enhanced = sce_LIBDsimu_ST_enhanced
BayesSpace_high_ST_result$sce_cluster = sce_LIBDsimu_ST_enhanced_labels
BayesSpace_high_ST_result$location = unique(data.frame("x"=vertices_location$y.pos,"y"=vertices_location$x.pos))
coordinate_allsubspot = c()
for(k in 1:dim(location_subspot)[1]){
  x = round(location_subspot[k,1],3)
  y = round(location_subspot[k,2],3)
  coordinate_allsubspot[k] = paste0(x,"_",y)
}
coordinate_BayesSpacesubspot = c()
for(k in 1:dim(BayesSpace_high_ST_result$location)[1]){
  x = round(BayesSpace_high_ST_result$location[k,1],3)
  y = round(BayesSpace_high_ST_result$location[k,2],3)
  coordinate_BayesSpacesubspot[k] = paste0(x,"_",y)
}
id_empty = which(BayesSpace_high_ST_result$truth=="empty")
ind_match_totalsubspots = match(coordinate_BayesSpacesubspot, coordinate_allsubspot)
BayesSpace_high_ST_result$truth =truth_subspot[ind_match_totalsubspots]
BayesSpace_high_ST_result$ARI = adjustedRandIndex(BayesSpace_high_ST_result$sce_cluster[-id_empty],BayesSpace_high_ST_result$truth[-id_empty])
BayesSpace_high_ST_result$NMI = compare(as.factor(BayesSpace_high_ST_result$sce_cluster[-id_empty]),as.factor(BayesSpace_high_ST_result$truth[-id_empty]), method = "nmi")
BayesSpace_high_ST_result$CHAOS = fx_chaos(BayesSpace_high_ST_result$sce_cluster[-id_empty], BayesSpace_high_ST_result$location[-id_empty,])
BayesSpace_high_ST_result$PAS = fx_PAEP(BayesSpace_high_ST_result$sce_cluster[-id_empty], BayesSpace_high_ST_result$location[-id_empty,])
print(BayesSpace_high_ST_result$ARI )

save(BayesSpace_high_ST_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/result/BayesSpace_sce_LIBDsimu_ST_enhanced_scenario_",i,"_rep_",j,".RData"))

#-------------------------------------- 
# run SpatialPCA, high resolution, ST
#-------------------------------------- 

print("run SpatialPCA, high resolution, ST")
newloc_ST = BayesSpace_high_ST_result$location[-id_empty,] # note, here x and y switched in BayesSpace vertex
LIBDsimu_high_ST = SpatialPCA_highresolution(LIBDsimu,newloc_ST)
expr_high = LIBDsimu@W %*% LIBDsimu_high_ST@highPCs
rownames(expr_high) = rownames(LIBDsimu@normalized_expr)
colnames(expr_high) = vertices_loc$spot[-id_empty]
SpatialPCA_high_ST_result = list()
SpatialPCA_high_ST_result$LIBDsimu_high_ST = LIBDsimu_high_ST
pred_cluster= louvain_clustering_new(4,LIBDsimu_high_ST@highPCs,1000 )
SpatialPCA_high_ST_result$clusterlabel = pred_cluster
SpatialPCA_high_ST_result$highPCs = LIBDsimu_high_ST@highPCs
SpatialPCA_high_ST_result$location = LIBDsimu_high_ST@highPos
SpatialPCA_high_ST_result$expr_high = expr_high
SpatialPCA_high_ST_result$truth =BayesSpace_high_ST_result$truth[-id_empty]
SpatialPCA_high_ST_result$ARI = adjustedRandIndex(SpatialPCA_high_ST_result$clusterlabel,SpatialPCA_high_ST_result$truth)
SpatialPCA_high_ST_result$NMI = compare(as.factor(SpatialPCA_high_ST_result$clusterlabel),as.factor(SpatialPCA_high_ST_result$truth), method = "nmi")
SpatialPCA_high_ST_result$CHAOS = fx_chaos(SpatialPCA_high_ST_result$clusterlabel, SpatialPCA_high_ST_result$location)
SpatialPCA_high_ST_result$PAS = fx_PAEP(SpatialPCA_high_ST_result$clusterlabel, SpatialPCA_high_ST_result$location)
print(SpatialPCA_high_ST_result$ARI)
save(SpatialPCA_high_ST_result,  file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/result/spotlevel_SpatialPCA_highST_result_scenario_",i,"_rep_",j,".RData"))

#-------------------------------------- 
# run BayesSpace, prediction gene expression
#-------------------------------------- 

print("run BayesSpace, prediction gene expression")
gene_input = rownames(SpatialPCA_result$normalized_expr)
sce_pred = function(index){
  sce_enhanced_each <- enhanceFeatures(sce_LIBDsimu_ST_enhanced, sce_LIBDsimu_ST,
                                model="xgboost",
                                feature_names=gene_input[index],
                                nrounds=0)
  logcounts = logcounts(sce_enhanced_each)
  return(logcounts[which(rownames(logcounts) %in% gene_input[index]),])
}
expr_high_BayesSpace = matrix(NA, dim(expr_high)[1], dim(expr_high)[2])
for(index in 1:length(gene_input)){
  print(index)
  pred_expr =c(sce_pred(index))
  expr_high_BayesSpace[index,] = pred_expr[-id_empty]
}
rownames(expr_high_BayesSpace) = gene_input

save(expr_high_BayesSpace, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/result/expr_high_BayesSpace_result_scenairo_",i,"_rep_",j,".RData"))

coordinate_allsubspot = c()
for(k in 1:dim(location_subspot)[1]){
  x = round(location_subspot[k,1],3)
  y = round(location_subspot[k,2],3)
  coordinate_allsubspot[k] = paste0(x,"_",y)
}
coordinate_BayesSpacesubspot = c()
for(k in 1:dim(BayesSpace_high_ST_result$location)[1]){
  x = round(BayesSpace_high_ST_result$location[k,1],3)
  y = round(BayesSpace_high_ST_result$location[k,2],3)
  coordinate_BayesSpacesubspot[k] = paste0(x,"_",y)
}
ind_match_totalsubspots = match(coordinate_BayesSpacesubspot, coordinate_allsubspot)
real_subspot_count = count_location_subspot$count_spot[,ind_match_totalsubspots]
real_subspot_count = real_subspot_count[,-id_empty]
rownames(real_subspot_count) = paste0("gene",1:5000)

real_subspot_count_sub = real_subspot_count[match(rownames(expr_high),rownames(real_subspot_count)),]
ind_remove = which(is.na(real_subspot_count_sub[1,]))
dat_result = list()
dat_result$real_subspot_count=real_subspot_count
dat_result$real_subspot_count_sub=real_subspot_count_sub[,-ind_remove]
dat_result$expr_high=expr_high[,-ind_remove]
dat_result$expr_high_BayesSpace=expr_high_BayesSpace[,-ind_remove]

save(dat_result, file = paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/NC/Reviewer1/Q3/simulation/prediction/result/Compare_pred_expr_",i,"_rep_",j,".RData"))




```









