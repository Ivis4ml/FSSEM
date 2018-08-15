setwd("./src")
source("solver.R")
#### data simulation
seed = 36435431
cat("Init Seeds = ", seed, "\n")
N = 500         ## sample size
Ng = 10         ## gene number
Nk = 10 * 3     ## eQTL number
Ns = 3          ## sparse ratio
sigma = 0.5 
set.seed(seed)
data = getrandsem(
  N = N,
  Ng = Ng,
  Nk = Nk,
  r = Ns / Ng,
  d = 0.1,
  sigma = sigma,
  dag = T
)
Xs = submatX(data)
Ys = data$obs[c("Y1", "Y2")]

### training FSSEM model
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
  nlambda = 20,
  nrho = 20,
  threshold = 1e-4
)

cvlambda.opt = optimLambda_cv(cv.params, se = T, type = "loglik")
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
  use.strict = TRUE,
  acc = TRUE,
  sparse = T
)


### result validation based on learning error from FSSEM
{
  err = sqrt(cv.fit$err)
  metric = c(
    (
      TPR(cv.fit$B[[1]], data$var$B1, precision = err) + TPR(cv.fit$B[[2]], data$var$B2, precision = err)
    ) / 2,
    (
      FDR(cv.fit$B[[1]], data$var$B1, precision = err) + FDR(cv.fit$B[[2]], data$var$B2, precision = err)
    ) / 2,
    TPR2(cv.fit$B[[1]], cv.fit$B[[2]], data$var$B1, data$var$B2, precision = err),
    FDR2(cv.fit$B[[1]], cv.fit$B[[2]], data$var$B1, data$var$B2, precision = err)
  )
  cat("PD\tFDR\tPD_diff\tFDR_diff\n")
  cat(sprintf("%.3f\t%.3f\t%.3f\t%.3f\n", metric[1], metric[2], metric[3], metric[4]))
}

