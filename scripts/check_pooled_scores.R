library(jsonlite)
root <- normalizePath(commandArgs(trailingOnly=TRUE)[1]);run <- file.path(root,'data/gwas/pooled-measured-20261004')
manifest <- fromJSON(file.path(run,'manifest.private.json'));checks <- list()
for(i in seq_len(nrow(manifest$jobs))){
 folder <- file.path(root,manifest$jobs$folder[i]);info <- fromJSON(file.path(folder,'diagnostics.private.json'))
 if(!isTRUE(info$gmmat_null_converged)) next
 obj <- readRDS(file.path(folder,'gmmat_null.private.rds'));dogs <- read.delim(file.path(folder,'dogs.private.tsv'))
 stopifnot(identical(as.character(obj$id_include),as.character(dogs$id)),all(obj$X[,1]==1),all(obj$X[,2]==dogs$male),all(obj$X[,3]==dogs$gsp))
 g <- read.table(file.path(folder,'genotypes.txt'),header=FALSE,check.names=FALSE);result <- read.table(file.path(folder,'gmmat.score.txt'),header=TRUE)
 selected <- unique(c(1L,ceiling(nrow(g)/2),nrow(g)))
 for(j in selected){
  geno <- as.numeric(g[j,-1]);centered <- geno-mean(geno);score <- as.numeric(crossprod(centered,obj$residuals));variance <- as.numeric(crossprod(centered,obj$P%*%centered));p <- pchisq(score^2/variance,1,lower.tail=FALSE)
  stopifnot(as.character(g[j,1])==as.character(result$SNP[j]),abs(score-result$SCORE[j])<max(1e-6,abs(score)*1e-5),abs(variance-result$VAR[j])<max(1e-6,abs(variance)*1e-5),abs(p-result$PVAL[j])<max(1e-8,abs(p)*1e-5))
 }
 checks[[length(checks)+1L]] <- list(folder=manifest$jobs$folder[i],sample_order_and_covariates_verified=TRUE,score_reconstructions=length(selected))
}
write_json(list(status='pooled_sample_orders_covariates_and_selected_scores_verified',checks=checks),file.path(run,'score-check.private.json'),auto_unbox=TRUE,pretty=TRUE)
cat(length(checks),'pooled chromosome orders, covariates and selected scores verified.\n')
