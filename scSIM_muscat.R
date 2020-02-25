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
# create random gene sets
##########################################################
randomGeneSets<-function(a){
  gsets<-sapply( rep(50,1000) , function(x) {list(as.character(sample(rownames(a),x))) } )
  names(gsets)<-stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
  gsets
}

mygenes <- unique(gi$gene)

# TODO here need to set the ground truth gene sets

##########################################################
# analyse with mitch
##########################################################

x<-mitch_import(res$table$B,DEtype="muscat",geneIDcol = "gene")



#todo
#- process sim with muscat
t(head(assay(pb)))

metadata(sim)

str(gi)
head(gi,30)
dim(gi)
head(gi$category,30)
table(gi$category)
table(gi$cluster_id)
summary(gi$sim_disp)
str(sim)

table(sim$sample_id)

?metadata
str(metadata(sim))
