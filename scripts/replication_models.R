library(GMMAT)
library(jsonlite)
args <- commandArgs(trailingOnly=TRUE)
folder <- normalizePath(args[1]);setwd(folder)
meta <- fromJSON("input.private.json")
dogs <- read.delim("dogs.private.tsv",stringsAsFactors=FALSE)
K <- as.matrix(read.table("kinship.txt"))
rownames(K) <- colnames(K) <- dogs$id
stopifnot(nrow(K)==nrow(dogs), !anyNA(dogs), all(dogs$y %in% 0:1))
captured <- character()
capture <- function(expr) withCallingHandlers(expr,warning=function(w) {
    captured <<- c(captured,conditionMessage(w));invokeRestart("muffleWarning")
})
fit_null <- function(data) {
    obj <- glmmkin(y~male,data=data,kins=K,id="id",
        family=binomial(link="logit"),method="REML",method.optim="Brent",maxiter=500)
    if(!isTRUE(obj$converged)) stop("GMMAT null model did not converge")
    obj
}
run_score <- function(obj,name) {
    glmm.score(obj,infile="genotypes.txt",outfile=name,infile.sep="\t",infile.ncol.skip=1,
        infile.ncol.print=1,infile.header.print="SNP",center=TRUE,select=seq_len(nrow(dogs)))
    read.table(name,header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
}
run_gemma <- function(y,name,full=FALSE) {
    write.table(y,"current_phenotype.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
    argv <- c("-g","gemma_genotypes.txt","-p","current_phenotype.txt","-k","kinship.txt",
        "-c","covariates.txt","-lmm",if(full) "4" else "3","-maf","0","-miss","1","-o",name)
    status <- system2(meta$gemma,argv,stdout=paste0(name,".console.txt"),stderr=paste0(name,".console.txt"))
    if(status!=0) stop("GEMMA failed: ",status)
    result <- read.table(file.path("output",paste0(name,".assoc.txt")),header=TRUE,
        check.names=FALSE,colClasses="character")
    for(column in setdiff(names(result),c("chr","rs","allele1","allele0")))
        result[[column]] <- as.numeric(result[[column]])
    stopifnot(result$allele1==meta$alternate_allele,result$allele0==meta$reference_allele)
    result
}
set.seed(meta$seed)
null <- capture(fit_null(dogs))
saveRDS(null,"gmmat_null.private.rds")
score <- capture(run_score(null,"gmmat_observed.score.txt"))
wald <- tryCatch(capture(glmm.wald(y~male,data=dogs,kins=K,id="id",family=binomial(link="logit"),
    infile="genotypes.txt",snps=meta$marker,method.optim="Brent",infile.sep="\t",infile.ncol.skip=1,
    missing.method="omit",maxiter=500,select=seq_len(nrow(dogs)))),error=function(e) list(error=conditionMessage(e)))
gemma <- run_gemma(dogs$y,"gemma_observed",full=TRUE)
observed_warnings <- unique(captured); captured <- character()
# Conditional parametric null draws preserve K, sex, and the observed case count.
# This is a limited calibration screen, not a tail-accurate replication p-value.
tau <- as.numeric(null$theta[2])
alpha <- as.numeric(null$coefficients)
if(length(alpha)!=2) alpha <- as.numeric(null$alpha)
stopifnot(length(alpha)==2,all(is.finite(alpha)),is.finite(tau),tau>=0)
decomposition <- eigen(K,symmetric=TRUE)
rootK <- decomposition$vectors %*% diag(sqrt(pmax(decomposition$values,0)*tau))
X <- model.matrix(~male, dogs)
bootstrap <- list(); attempted <- 0L; failures <- character()
for (iteration in seq_len(meta$replicates)) {
    repeat {
        attempted <- attempted+1L
        if(attempted>100000L) stop("Conditional null sampler exceeded attempt limit")
        probability <- plogis(as.vector(X%*%alpha+rootK%*%rnorm(nrow(dogs))))
        simulated <- rbinom(nrow(dogs),1,probability)
        if(sum(simulated)==sum(dogs$y)) break
    }
    simulated_dogs <- dogs;simulated_dogs$y <- simulated
    values <- tryCatch({
        model <- capture(fit_null(simulated_dogs))
        gmmat <- capture(run_score(model,"calibration.score.txt"))
        linear <- run_gemma(simulated,"gemma_calibration")
        p <- c(as.numeric(gmmat$PVAL),as.numeric(linear$p_score))
        if(length(p)!=2 || any(!is.finite(p)) || any(p<0 | p>1)) stop("Invalid simulated p-value")
        list(GMMAT=p[1],GEMMA=p[2])
    },error=function(e) {failures <<- c(failures,conditionMessage(e));NULL})
    if(!is.null(values)) bootstrap[[length(bootstrap)+1L]] <- values
    if(iteration%%20==0) cat("Null screen:",iteration,"of",meta$replicates,"\n")
}
result <- list(GMMAT_score=score,GMMAT_wald=wald,GEMMA=gemma,
    null_theta=as.numeric(null$theta),null_coefficients=alpha,
    null_converged=isTRUE(null$converged),null_field_names=names(null),warnings=observed_warnings,
    calibration=list(requested=meta$replicates,completed=length(bootstrap),failed=length(failures),
        failures=failures,attempted_draws=attempted,null_simulation="conditional logistic mixed model, fixed case count",
        warnings=unique(captured),values=bootstrap),versions=list(R=R.version.string,GMMAT=as.character(packageVersion("GMMAT"))))
write_json(result,"result.private.json",auto_unbox=TRUE,pretty=TRUE,na="null",digits=16)
