# Bradley May 2022
# General pipeline for fine mapping in GTEx - any gene, any site
# PEER analysis of the TMM normalised, INT expression values
# Need to make sure this is done in R version 3.3.0 - Seems to be the only version can get PEER to run on

# IN this instance, doing this for multiple sites at the same time
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/3.3")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
library(optparse)
option_list = list(make_option(c("-s", "--site"), action = "store", default = NA, type ="character", help="cluster number"),
										make_option(c("-g", "--gene"), action = "store", default = NA, type ="character", help="gene name"))
opt = parse_args(OptionParser(option_list=option_list))
s = opt$s;
print(paste("The option for site to study is", s, sep = " "))
gene = opt$g
print(paste("The option for gene to study is", gene, sep = " "))

# Rename the s variable if using any with brackets
s<-gsub("\\(", "", s)
s<-gsub("\\)", "", s)

#~~~~~~ Set up
library("peer", lib="/exports/igmm/eddie/NextGenResources/software/R/x86_64-pc-linux-gnu-library/3.3")
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx")

#~~~~~~ Load the INT expression and design
expr <- read.csv(paste("data", s, gene, "expr_TMM_INT.csv", sep = "/"), row.names = 1)
covs <- read.csv(paste("data", s, gene, "meta_ready.csv",sep = "/"), row.names = 1)
dim(expr)
dim(covs)

#~~~~~~ Start the PEER analysis
#Build the model
set.seed(123)
model=PEER() #Build model
PEER_setPhenoMean(model, as.matrix(expr))
dim(PEER_getPhenoMean(model))

# Set up the parameters
max_fac <- floor(nrow(expr)/4)
if(max_fac > 100){
	max_fac <- 100
}
PEER_setNk(model, max_fac) #Set the maximum number of unobserved factors to model
PEER_setNmax_iterations(model, 10000)

# Set the observed covariates
covs[,1:ncol(covs)]<-sapply(covs[,1:ncol(covs)],as.numeric)
PEER_setCovariates(model, as.matrix(covs)) #Set the observed covariaties

# using the expr and covariates, update the model
PEER_update(model) #Does 1000 iterations. If this is not converged after, and the variance of residuals keeps decreasing, choose a higher value, e.g 10,000: PEER_setNMax_iterations(model, 10000)


#~~~~~~ Save the results of the analysis
if(file.exists(paste("results", gsub("\\ ", "", s), sep = "/")) == F){
	dir.create(paste("results", gsub("\\ ", "", s), sep = "/"))
	if(file.exists(paste("results", gsub("\\ ", "", s), gene, sep = "/")) == F){
		dir.create(paste("results", gsub("\\ ", "", s), gene, sep = "/"))
		dir.create(paste("results", gsub("\\ ", "", s), gene, "PEER", sep = "/"))
	}
}
pathOut <- paste("results", gsub("\\ ", "", s), gene, "PEER", sep = "/")

factors = PEER_getX(model)
dim(factors)
rownames(factors) <- rownames(covs)
write.csv(factors, paste(pathOut, "Peer_factors_age_gender_batch.csv", sep = "/"))

weights = PEER_getW(model)
dim(weights)
write.csv(weights, paste(pathOut, "Peer_weights_age_gender_batch.csv", sep ="/"))

precision = PEER_getAlpha(model)
dim(precision)
write.csv(precision, paste(pathOut, "Peer_precision_age_gender_batch.csv", sep ="/"))


residuals = PEER_getResiduals(model)
rownames(residuals) <- rownames(expr)
colnames(residuals) <- colnames(expr)
dim(residuals)
write.csv(residuals, paste(pathOut, "Peer_residuals_age_gender_batch.csv", sep = "/"))


pdf(file = paste(pathOut, "precision_age_gender_batch.pdf", sep = "/"))
plot(precision)
dev.off()
pdf(paste(pathOut, "PEER_model_age_gender_batch.pdf", sep="/")); PEER_plotModel(model); dev.off()
