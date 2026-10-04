library(GMMAT)
library(jsonlite)
root <- normalizePath(commandArgs(trailingOnly=TRUE)[1]);folder <- file.path(root,'results/private/version2/08_cross_breed')
cfg <- fromJSON(file.path(root,'config/cross_breed.json'));checks <- list()
for(scenario in cfg$schemes){
 input <- file.path(folder,paste0(scenario,'_',c('LR','GSP'),'_score.private.tsv'))
 output <- file.path(folder,paste0(scenario,'_native_gmmat_meta.private.tsv'))
 glmm.score.meta(input,output)
 native <- read.delim(output,check.names=FALSE)
 data <- read.delim(gzfile(file.path(folder,paste0(scenario,'.private.tsv.gz'))),check.names=FALSE)
 stopifnot(identical(as.character(native$SNP),as.character(data$marker)),all(native$A1==data$ALT),all(native$A2==data$REF))
 stopifnot(max(abs(native$SCORE-data$score_meta_score))<1e-5,max(abs(native$VAR-data$score_meta_variance))<1e-5,max(abs(native$PVAL-data$score_meta_p))<1e-7)
 # Independent base-R calculation of every effect meta-analysis and Fisher probability.
 w1 <- 1/data$LR_GEMMA_se^2;w2 <- 1/data$GSP_GEMMA_se^2
 beta <- (w1*data$LR_GEMMA_beta+w2*data$GSP_GEMMA_beta)/(w1+w2);se <- sqrt(1/(w1+w2))
 valid <- is.finite(beta)&is.finite(se)&se>0
 stopifnot(max(abs(beta[valid]-data$iv_beta[valid]))<1e-10,max(abs(se[valid]-data$iv_se[valid]))<1e-10,max(abs(2*pnorm(-abs(beta[valid]/se[valid]))-data$iv_p[valid]))<1e-10)
 for(key in cfg$fisher_tests){
  p1 <- data[[paste0('LR_',key)]];p2 <- data[[paste0('GSP_',key)]];valid <- is.finite(p1)&is.finite(p2)&p1>0&p2>0&p1<=1&p2<=1
  expected <- pchisq(-2*(log(p1[valid])+log(p2[valid])),df=4,lower.tail=FALSE)
  stopifnot(max(abs(expected-data[[paste0(key,'_fisher')]][valid]))<1e-10)
 }
 checks[[scenario]] <- list(markers=nrow(data),native_gmmat_score_meta_verified=TRUE,effect_meta_verified=TRUE,fisher_verified=TRUE)
}
write_json(list(status='all_cross_breed_meta_calculations_verified',checks=checks),file.path(folder,'verification.private.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Native GMMAT meta-analysis and independent R effect/Fisher calculations verified.\n')
