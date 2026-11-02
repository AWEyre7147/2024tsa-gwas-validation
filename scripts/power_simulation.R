library(GMMAT)
library(jsonlite)
folder <- normalizePath(commandArgs(trailingOnly=TRUE)[1]); setwd(folder)
meta <- fromJSON('input.private.json')
dogs <- read.delim('dogs.private.tsv',stringsAsFactors=FALSE)
K <- as.matrix(read.table('kinship.txt')); rownames(K) <- colnames(K) <- dogs$id
g <- as.numeric(unlist(read.table('genotypes.txt',header=FALSE)[1,-1]))
stopifnot(length(g)==nrow(dogs),all(is.finite(g)),sum(dogs$y)==meta$cases,
          all(g>=0 & g<=2),nrow(K)==nrow(dogs),max(abs(K-t(K)))<1e-10)
original <- readRDS(meta$null_model)
stopifnot(isTRUE(original$converged),identical(as.character(original$id_include),as.character(dogs$id)))
alpha <- as.numeric(original$coefficients)
if(length(alpha)!=2) alpha <- as.numeric(original$alpha)
tau <- as.numeric(original$theta[2]); stopifnot(length(alpha)==2,is.finite(tau),tau>=0)
X <- model.matrix(~male,dogs)
eig <- eigen(K,symmetric=TRUE)
stopifnot(min(eig$values)> -1e-8)
# Standard normal Gauss-Hermite quadrature: expected case count is recentered
# once per scenario, not once per drawn latent effect.
q <- 32L; J <- matrix(0,q,q)
for(i in 1:(q-1)) J[i,i+1] <- J[i+1,i] <- sqrt(i/2)
qe <- eigen(J,symmetric=TRUE); nodes <- sqrt(2)*qe$values; weights <- qe$vectors[1,]^2
warnings <- character()
capture <- function(expr) withCallingHandlers(expr,warning=function(w){
  warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
})
records <- list(); index <- 0L
for(multiplier in meta$variance_multipliers) for(OR in meta$odds_ratios) {
  scenario <- paste0('OR',OR,'_variance',multiplier)
  beta <- log(OR); tau_sim <- tau*multiplier
  root <- eig$vectors %*% diag(sqrt(pmax(eig$values,0)*tau_sim))
  fixed <- X[,2]*alpha[2]+g*beta
  means <- function(intercept) sum(vapply(seq_len(nrow(dogs)),function(i)
      sum(weights*plogis(intercept+fixed[i]+sqrt(max(K[i,i],0)*tau_sim)*nodes)),numeric(1)))
  intercept <- uniroot(function(a) means(a)-meta$cases,c(-100,100),tol=1e-10)$root
  # Separate deterministic RNG stream for every context/scenario.
  set.seed(as.integer(meta$seed+round(multiplier*10000)+round(OR*100)))
  for(iteration in seq_len(meta$replicates)) {
    warnings <- character(); accepted <- FALSE; attempted <- 0L
    while(attempted < meta$maximum_attempts) {
      attempted <- attempted+1L
      y <- rbinom(nrow(dogs),1,plogis(intercept+fixed+as.vector(root%*%rnorm(nrow(dogs)))))
      if(sum(y)==meta$cases) {accepted <- TRUE; break}
    }
    r <- list(context=meta$context,scenario=scenario,OR=OR,variance_multiplier=multiplier,
      replicate=iteration,generating_intercept=intercept,expected_cases=means(intercept),
      generating_tau=tau_sim,attempted_draws=attempted,generated=accepted,outcome_bits='',
      GMMAT_score=NA_real_,GEMMA_wald=NA_real_,GEMMA_lrt=NA_real_,GEMMA_score=NA_real_,
      GMMAT_converged=FALSE,GEMMA_exit=NA_integer_,GEMMA_mle_boundary=NA,
      GEMMA_beta=NA_real_,GEMMA_SE=NA_real_,error='',warnings='')
    if(accepted) {
      r$outcome_bits <- paste0(y,collapse='')
      simulated <- dogs; simulated$y <- y
      null <- tryCatch(capture(glmmkin(y~male,data=simulated,kins=K,id='id',family=binomial(link='logit'),
          method='REML',method.optim='Brent',maxiter=500)),error=function(e)list(error=conditionMessage(e)))
      r$GMMAT_converged <- isTRUE(null$converged)
      if(r$GMMAT_converged) {
        score <- tryCatch({
          stopifnot(identical(as.character(null$id_include),as.character(dogs$id)))
          capture(glmm.score(null,infile='genotypes.txt',outfile='current.score.txt',select=seq_len(nrow(dogs)),
            infile.sep='\t',infile.ncol.skip=1,infile.ncol.print=1,infile.header.print='SNP',center=TRUE,
            MAF.range=c(0,.5),missing.method='impute2mean'))
          read.table('current.score.txt',header=TRUE,stringsAsFactors=FALSE)
        },error=function(e)list(error=conditionMessage(e)))
        if(!is.null(score$PVAL)) r$GMMAT_score <- as.numeric(score$PVAL[1])
        else r$error <- paste(r$error,score$error)
      } else r$error <- paste(r$error,'GMMAT null fit failed',null$error)
      write.table(y,'current_phenotype.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
      r$GEMMA_exit <- system2(meta$gemma,c('-g','gemma_genotypes.txt','-p','current_phenotype.txt','-k','kinship.txt',
          '-c','covariates.txt','-lmm','4','-maf','0','-miss','1','-o','current'),stdout='current.gemma.console.txt',stderr='current.gemma.console.txt')
      if(r$GEMMA_exit==0L) {
        linear <- tryCatch(read.table('output/current.assoc.txt',header=TRUE,colClasses='character'),error=function(e)NULL)
        if(!is.null(linear) && nrow(linear)==1) {
          stopifnot(linear$allele1==meta$ALT,linear$allele0==meta$REF)
          r$GEMMA_wald <- as.numeric(linear$p_wald); r$GEMMA_lrt <- as.numeric(linear$p_lrt)
          r$GEMMA_score <- as.numeric(linear$p_score); r$GEMMA_beta <- as.numeric(linear$beta);r$GEMMA_SE <- as.numeric(linear$se)
          ratio <- as.numeric(linear$l_mle)
          r$GEMMA_mle_boundary <- is.finite(ratio) && (ratio<=1.00001e-5 || ratio>=99999)
        } else r$error <- paste(r$error,'GEMMA association row unavailable')
      } else r$error <- paste(r$error,'GEMMA failed')
    } else r$error <- 'Conditional sampler attempt limit reached'
    r$warnings <- paste(unique(warnings),collapse=' | ')
    index <- index+1L; records[[index]] <- r
  }
  cat(meta$context,scenario,meta$replicates,'replicates completed\n');flush.console()
}
result <- do.call(rbind,lapply(records,as.data.frame,stringsAsFactors=FALSE))
write.table(result,'draws.private.tsv',sep='\t',row.names=FALSE,quote=TRUE,na='NA')
write_json(list(R=R.version.string,GMMAT=as.character(packageVersion('GMMAT')),
  null_sex_coefficient=alpha[2],null_variance=tau,case_count=meta$cases,
  simulation='Scenario-recentered logistic mixed model, conditioned on observed case count',
  source_null=meta$null_model,replicates=nrow(result)), 'simulation_metadata.private.json',pretty=TRUE,auto_unbox=TRUE,digits=16)
