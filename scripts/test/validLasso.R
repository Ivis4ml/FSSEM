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

set.seed(as.integer(seed))
seeds = sample.int(1e8, size = 5, replace = F)

iterator = 1

cv.err   = NULL
cv.1se   = NULL
## repeat 2 times
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
    dag = dag
  )
  
  Xs = submatX(data)
  Ys = data$obs[c("Y1", "Y2")]
  
  sigma2.init = getsigma2_L2reg_multi(Xs,
                                      Ys,
                                      nrho = 50,
                                      M = Ng,
                                      N = N)
  params.init = constrained_L2reg_multi(Xs, Ys, sigma2.init$rho.opt, M = Ng, N = N)
  cv.params = cv_SMLasso(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    nlambda = 20
  )
  
  opt.lasso = optLasso_cv(cv.params, se = FALSE)
  
  
  cv.fit1 = sparse_maximum_likehood_iPALM(
    B = params.init$B[[1]],
    f = params.init$F[[1]],
    Y = Ys[[1]],
    X = Xs,
    sigma2 = params.init$sigma2[1],
    N = data$var$N,
    Ng = data$var$Ng,
    lambda = opt.lasso,
    maxit = 2000
  )
  
  cv.fit2 = sparse_maximum_likehood_iPALM(
    B = params.init$B[[2]],
    f = params.init$F[[2]],
    Y = Ys[[2]],
    X = Xs,
    sigma2 = params.init$sigma2[1],
    N = data$var$N,
    Ng = data$var$Ng,
    lambda = opt.lasso,
    maxit = 2000
  )
  
  cv.err = rbind(cv.err,
                 c(
                   (TPR(cv.fit1$B, data$var$B1) + TPR(cv.fit2$B, data$var$B2)) / 2,
                   (FDR(cv.fit1$B, data$var$B1) + FDR(cv.fit2$B, data$var$B2)) / 2,
                   TPR2(cv.fit1$B, cv.fit2$B, data$var$B1, data$var$B2),
                   FDR2(cv.fit1$B, cv.fit2$B, data$var$B1, data$var$B2)
                 ))
  
  ###############
  opt.1se = optLasso_cv(cv.params, se = TRUE)
  
  
  cv.1se1 = sparse_maximum_likehood_iPALM(
    B = params.init$B[[1]],
    f = params.init$F[[1]],
    Y = Ys[[1]],
    X = Xs,
    sigma2 = params.init$sigma2[1],
    N = data$var$N,
    Ng = data$var$Ng,
    lambda = opt.1se,
    maxit = 2000
  )
  
  cv.1se2 = sparse_maximum_likehood_iPALM(
    B = params.init$B[[2]],
    f = params.init$F[[2]],
    Y = Ys[[2]],
    X = Xs,
    sigma2 = params.init$sigma2[1],
    N = data$var$N,
    Ng = data$var$Ng,
    lambda = opt.1se,
    maxit = 2000
  )
  
  cv.1se = rbind(cv.1se,
                 c(
                   (TPR(cv.1se1$B, data$var$B1) + TPR(cv.1se2$B, data$var$B2)) / 2,
                   (FDR(cv.1se1$B, data$var$B1) + FDR(cv.1se2$B, data$var$B2)) / 2,
                   TPR2(cv.1se1$B, cv.1se2$B, data$var$B1, data$var$B2),
                   FDR2(cv.1se1$B, cv.1se2$B, data$var$B1, data$var$B2)
                 ))
  
  iterator = iterator + 1
}

saveRDS(
  list(
    error.cv = cv.err,
    error.se = cv.1se,
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
