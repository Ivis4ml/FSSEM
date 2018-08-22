source("solver.min.R")
library(optparse)

option_list = list(
  make_option(
    c("-n", "--nobs"),
    type = "numeric",
    default = 100,
    help = "sample number [default = %default]",
    metavar = "numeric"
  ),
  make_option(
    c("-p", "--gene"),
    type = "numeric",
    default = 10,
    help = "gene number [default = %default]",
    metavar = "numeric"
  ),
  make_option(
    c("-k", "--eqtl"),
    type = "numeric",
    default = 3,
    help = "eQTL number [default = %default]",
    metavar = "numeric"
  ),
  make_option(
    c("-c", "--supp"),
    type = "numeric",
    default = 3,
    help = "support coeffcient gene number [default = %default]",
    metavar = "numeric"
  ),
  make_option(
    c("-s", "--sigma"),
    type = "numeric",
    default = 1,
    help = "noise on SEM's observations",
    metavar = "numeric"
  ),
  make_option(
    c("-d", "--seed"),
    type = "numeric",
    default = 1,
    help = "initialized seed for network structure",
    metavar = "numeric"
  ),
  make_option(
    c("-a", "--dag"),
    type = "numeric",
    default = 1,
    help = "dag or dcg for network structure",
    metavar = "numeric"
  )
)

option_parser = OptionParser(option_list = option_list)
args.opt = parse_args(option_parser)

N  = args.opt$nobs
Ng = args.opt$gene
Ne = args.opt$eqtl
Nk = Ng * Ne
sigma = args.opt$sigma
seed  = args.opt$seed
## support genes for each gene
Ns = args.opt$supp
## DAG or DCG
dag = ifelse(args.opt$dag == 1, TRUE, FALSE)

set.seed(seed)
seeds = sample.int(1e8, size = 5, replace = F)


iterator = 1

cv.err   = NULL
cv.1se   = NULL
cv.ll    = NULL
## repeat 3 times
while (iterator <= 1) {
  seed = seeds[iterator]
  cat("Init Seeds = ", seed, "\n")
  set.seed(seed)
  data = getrandsem(
    N = N,
    Ng = Ng,
    Nk = Nk,
    r = Ns / Ng,
    d = 0.1,
    sigma = sigma,
    dag = TRUE
  )
  
  Xs = submatX(data)
  Ys = data$obs[c("Y1", "Y2")]
  
  sigma2.init = getsigma2_L2reg_multi(Xs,
                                      Ys,
                                      nrho = 50,
                                      M = Ng,
                                      N = N)
  params.init = constrained_L2reg_multi(Xs, Ys, sigma2.init$rho.opt, M = Ng, N = N)
  
  cv.params = cv_multiSML(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    nlambda = 10,
    nrho = 10,
    threshold = 1e-4
  )
  
  cvlambda.opt = optimLambda_cv(cv.params, se = FALSE, type = "err")
  cv.fit = multiSML_iPALM(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    lambda = cvlambda.opt$lambda,
    rho = cvlambda.opt$rho,
    maxit = 2000,
    threshold = 1e-6,
    use.strict = FALSE,
    acc = TRUE,
    sparse = T
  )
  
  cv.err = rbind(cv.err,
                 c(
                   (TPR(cv.fit$B[[1]], data$var$B1) + TPR(cv.fit$B[[2]], data$var$B2)) / 2,
                   (FDR(cv.fit$B[[1]], data$var$B1) + FDR(cv.fit$B[[2]], data$var$B2)) / 2,
                   TPR2(cv.fit$B[[1]], cv.fit$B[[2]], data$var$B1, data$var$B2),
                   FDR2(cv.fit$B[[1]], cv.fit$B[[2]], data$var$B1, data$var$B2)
                 ))
  
  ### 1se selection
  cvlambda.1se = optimLambda_cv(cv.params, se = TRUE, type = "err")
  se.fit = multiSML_iPALM(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    lambda = cvlambda.1se$lambda,
    rho = cvlambda.1se$rho,
    maxit = 2000,
    threshold = 1e-6,
    use.strict = FALSE,
    acc = TRUE,
    sparse = T
  )
  
  cv.1se = rbind(cv.1se,
                 c(
                   (TPR(se.fit$B[[1]], data$var$B1) + TPR(se.fit$B[[2]], data$var$B2)) / 2,
                   (FDR(se.fit$B[[1]], data$var$B1) + FDR(se.fit$B[[2]], data$var$B2)) / 2,
                   TPR2(se.fit$B[[1]], se.fit$B[[2]], data$var$B1, data$var$B2),
                   FDR2(se.fit$B[[1]], se.fit$B[[2]], data$var$B1, data$var$B2)
                 ))
  
  cvlambda.opt = optimLambda_cv(cv.params, se = FALSE, type = "loglik")
  cv.fit = multiSML_iPALM(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    lambda = cvlambda.opt$lambda,
    rho = cvlambda.opt$rho,
    maxit = 2000,
    threshold = 1e-6,
    use.strict = FALSE,
    acc = TRUE,
    sparse = T
  )
  
  cv.ll = rbind(cv.ll,
                 c(
                   (TPR(cv.fit$B[[1]], data$var$B1) + TPR(cv.fit$B[[2]], data$var$B2)) / 2,
                   (FDR(cv.fit$B[[1]], data$var$B1) + FDR(cv.fit$B[[2]], data$var$B2)) / 2,
                   TPR2(cv.fit$B[[1]], cv.fit$B[[2]], data$var$B1, data$var$B2),
                   FDR2(cv.fit$B[[1]], cv.fit$B[[2]], data$var$B1, data$var$B2)
                 ))
  
  rm(data, Xs, Ys, sigma2.init, params.init, cv.params, cv.fit, se.fit)
  gc()
  iterator = iterator + 1
} ## replications

colnames(cv.err) = c("TPR", "FDR", "TPR.2", "FDR.2")

saveRDS(
  list(
    error.cv = cv.err,
    error.se = cv.1se,
    LL.cv    = cv.ll,
    params = list(
      N = N,
      Ng = Ng,
      Ne = Ne,
      sigma2 = sigma ^ 2
    ),
    seed = seeds
  ),
  file =  sprintf("n%dg%de%ds%.2f.rds", N, Ng, Ne, sigma)
)
