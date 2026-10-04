library(GMMAT)
library(jsonlite)
folder <- normalizePath(commandArgs(trailingOnly=TRUE)[1]);setwd(folder)
meta <- fromJSON('input.private.json');dogs <- read.delim('dogs.private.tsv',stringsAsFactors=FALSE)
K <- as.matrix(read.table('kinship.txt'));rownames(K)<-colnames(K)<-dogs$id
warnings <- character()
capture <- function(expr) withCallingHandlers(expr,warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
null <- tryCatch(capture(glmmkin(y~male+gsp,data=dogs,kins=K,id='id',family=binomial(link='logit'),method='REML',method.optim='Brent',maxiter=500)),error=function(e) list(error=conditionMessage(e)))
valid <- isTRUE(null$converged)
if(valid){
 stopifnot(identical(as.character(null$id_include),as.character(dogs$id)))
 saveRDS(null,'gmmat_null.private.rds')
 capture(glmm.score(null,infile='genotypes.txt',outfile='gmmat.score.txt',select=seq_len(nrow(dogs)),infile.sep='\t',infile.ncol.skip=1,infile.header.print='SNP',center=TRUE,MAF.range=c(0,0.5),missing.method='impute2mean',nperbatch=100))
}
status <- system2(meta$gemma,c('-g','gemma_genotypes.txt','-p','phenotype.txt','-k','kinship.txt','-c','covariates.txt','-lmm','4','-maf','0','-miss','1','-o','gemma'),stdout='gemma.console.txt',stderr='gemma.console.txt')
if(status!=0) stop('GEMMA failed: ',status)
write_json(list(gmmat_null_converged=valid,gmmat_null_error=null$error,theta=if(valid) as.numeric(null$theta) else NULL,warnings=unique(warnings),gemma_exit_status=status,versions=list(R=R.version.string,GMMAT=as.character(packageVersion('GMMAT')))),'diagnostics.private.json',auto_unbox=TRUE,pretty=TRUE,na='null',digits=16)
