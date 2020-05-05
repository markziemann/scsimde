library(SingleCellExperiment)
library(scater)
library(ExperimentHub)
library(muscat)
library(sctransform)
library(mitch)

##########################################################
# Prepare the KANG data to use as a template for simulation 
##########################################################
# Download the KANG data
eh <- ExperimentHub()
query(eh, "Kang")
(sce <- eh[["EH2259"]])

# Filter out low count cells
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
qc <- perCellQCMetrics(sce)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
# Filter out low count genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Perform some normalisations
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assays(sce)$vstresiduals <- vst(counts(sce), show_progress = FALSE)$y

# Tidy up some of the data
sce$id <- paste0(sce$stim, sce$ind)
(sce <- prepSCE(sce, 
                cluster_id = "cell", # subpopulation assignments
                group_id = "stim",   # group IDs (ctrl/stim)
                sample_id = "id",    # sample IDs (ctrl/stim.1234)
                drop = TRUE))        # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# have a look at the structure of the data
# str(sce)
# str(sce@colData@listData$cluster_id)
# sce@colData@listData$cluster_id

##########################################################
# perform the simulation using Muscat simulation functions
##########################################################
# prep. SCE for simulation
# need to see why we are losing some cell types at this step.
sce <- prepSim(sce, verbose = TRUE)

# simulate data with 10% DE
(sim <- simData(sce,
                n_genes = 5000, n_cells = 20000,
                p_dd = c(0.9, 0, 0.1, 0, 0, 0)))

# simulation metadata
head(gi <- metadata(sim)$gene_info)

sim <- sim[rowSums(counts(sim) > 0) > 0, ]

##########################################################
# process the simulated data with Muscat pseudobulk function
##########################################################
pb <- aggregateData(sim,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

res <- pbDS(pb, verbose = TRUE)

tbl <- res$table$B


for (CELL_ID in names(tbl)) {
  
  mytbl <- metadata(sim)$gene_info[which(metadata(sim)$gene_info$cluster_id == CELL_ID),]
  genes_up <- mytbl[which(mytbl$logFC>0),1]
  genes_dn <- mytbl[which(mytbl$logFC<0),1]
  
}

##########################################################
# Extract the ground truth DEGs
##########################################################
y <- vector("list", length(unique(gi$cluster_id)))

n=0
for ( CELL_ID in unique(gi$cluster_id) ) {
  n=n+1
  y1 <- subset(gi,logFC>0 & cluster_id==CELL_ID)$gene
  y2 <- subset(gi,logFC<0 & cluster_id==CELL_ID)$gene

  yy <- list(y1,y2)
  names(yy)=c("up","dn")

  y[n] <- list(yy)
}
names(y) <- unique(gi$cluster_id)

##########################################################
# create random gene sets
##########################################################
randomGeneSets <- function(a) {
  gsets<-sapply( rep(50,1000) , function(x) {
    list(as.character(sample(a,x))) 
  } )
  names(gsets) <- stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
  gsets
}

library(stringi)
mygenes <- unique(gi$gene)
randomSets <- randomGeneSets(mygenes)

##########################################################
# create gene sets with a fraction that are DE
##########################################################
gs <- vector("list", length(unique(gi$cluster_id)))

# gene sets will consist of nde #genes that are de and nee #genes that are not de (ee)
nde=50
nee=30
nsets=10

n=0
for ( CELL_ID in unique(gi$cluster_id) ) {
  n=n+1
  # first make sets of upregulated genes
  i=0
  up <- vector("list", nsets)
  for ( i in seq(nsets) ) {
    deg <- sample(subset(gi,logFC>0 & cluster_id==CELL_ID)$gene,nde)
    ndeg <- sample(unique(gi$gene),nee)
    mygs <- union(deg,ndeg)
    up[i] <- list(mygs)
    i=i+1
  }
  names(up) <- stri_rand_strings(length(up), 15, pattern = "[A-Za-z]")

  # now make sets of downregulated genes
  i=0
  dn <- vector("list", nsets)
  for ( i in seq(nsets) ) {
    deg <- sample(subset(gi,logFC<0 & cluster_id==CELL_ID)$gene,nde)
    ndeg <- sample(unique(gi$gene),nee)
    mygs <- union(deg,ndeg)
    dn[i] <- list(mygs)
    i=i+1  
  }
  names(dn) <- stri_rand_strings(length(up), 15, pattern = "[A-Za-z]")

  cell<-list("up"=up,"dn"=dn)

  gs[n] <- list(cell)

}
names(gs) <- unique(gi$cluster_id)

# flatten the list
desets <- unlist(unlist(gs,recursive=FALSE),recursive=FALSE)

#########################################################
# combine DE genesets with some random genesets
#########################################################

randomGeneSets <- function(a) {
  gsets<-sapply( rep(50,920) , function(x) {
    list(as.character(sample(a,x)))
  } )
  names(gsets) <- stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
  gsets
}

library(stringi)
mygenes <- unique(gi$gene)
randomSets <- randomGeneSets(mygenes)

# merge the random ones and the DE sets
mysets <- c(desets,randomSets)

##########################################################
# analyse with mitch
##########################################################


x <- mitch_import(res$table$B, DEtype = "muscat", geneIDcol = "gene")

mitch_res <- mitch_calc(x,mysets)

head(mitch_res$enrichment_result)

#todo need to perform p-adjust for each cell type and then calculate the sensitivity, specificity and F1 score
sig <- subset(mitch_res$enrichment_result, p.adjustMANOVA <= 0.05)

obs <- sig$set

TP = length(intersect(obs,names(desets)))
FP = length(setdiff(obs,names(desets)))
FN= length(setdiff(names(desets),obs))
TN= length(mysets) - (TP+FP+FN)
sensitivity=TP/(TP+FN)
specificity=TN/(TN+FP)
precision=TP/(TP+FP)
recall=TP/(TP+FN)
F1=2*precision*recall/(precision+recall)
str(F1)
#p.adjust()
save.image("scrna.Rdata")
