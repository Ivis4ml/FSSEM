## ridge regression for estimate sigma2 in gene expression
#' @example
#' X = submatX(data)
#' Y = data$obs$Y1
#' B = constrained_L2reg(X, Y, rho = 0.1)
constrained_L2reg = function(X, Y, rho) {
  # gene number(M) & sample number(N)
  std = centralize(X, Y)
  X = std$X
  Y = std$Y
  meanX = std$muX
  meanY = std$muY
  M = ncol(Y)
  N = nrow(Y)
  B = Matrix(0,
             nrow = M,
             ncol = M,
             sparse = T)
  f = list()
  err = 0
  for (i in 1:M) {
    Xi = X[[i]]                   ## n x sk
    Pi = diag(N) - projection(Xi) ## n x n
    yi = Y[, i, drop = F]         ## n x 1
    Yi = Y[,-i, drop = F]         ## n x (p-1)
    ## (Y^TPY + rho)^{-1}Y^TPy
    bi = solve(crossprod(Yi, Pi %*% Yi) + rho * diag(M - 1)) %*% t(Yi) %*% Pi %*% yi
    ## bi = glmnet(Pi %*% Yi, Pi %*% yi, alpha = 0, lambda = rho)[["beta"]][, 1]
    B[i, -i] = bi
    f[[i]] = solve(crossprod(Xi)) %*% t(Xi) %*% (yi - Yi %*% bi)
    err = err + crossprod(yi - Yi %*% bi - Xi %*% f[[i]])
  }
  sigma2 = err / (M * N - 1)
  mu     = (diag(M) - B) %*% meanY - sapply(1:M, function(i) {
    meanX[[i]] %*% f[[i]]
  })
  list(
    B = as.matrix(B),
    F = f,
    sigma2 = sigma2,
    mu = mu
  )
}

## cross-validation on ridge regression to estimate sigma2
#' @param nrho number of L2 penalty's coefficient
#' @param ncv  number of cross-validation
#' @example
#' sigma2 = getsigma2_L2reg(X, Y, nrho = 15, ncv = 5)
getsigma2_L2reg = function(X, Y, nrho = 10, ncv = 5) {
  rho_factors = 10 ** (seq(-6, 2, length.out = nrho))
  N = ncol(Y)
  M = nrow(Y)
  cv.err  = matrix(0, nrow = nrho, ncol = ncv)
  cv.fold = sample(seq(1, ncv), size = N, replace = T)
  irho    = 1
  for (rho in rho_factors) {
    for (cv in 1:ncv) {
      ytrain = Y[, cv.fold != cv]
      xtrain = lapply(X, function(x) {
        x[, cv.fold != cv, drop = F]
      })
      ytest  = Y[, cv.fold == cv]
      xtest  = lapply(X, function(x) {
        x[, cv.fold == cv, drop = F]
      })
      fit    = constrained_L2reg(xtrain, ytrain, rho)
      ftest  = lapply(1:M, function(i) {
        crossprod(fit$F[[i]], xtest[[i]])
      })
      ftest  = do.call(rbind, ftest)
      cv.err[irho, cv] = norm((diag(M) - fit$B) %*% ytest - ftest - fit$mu, type = "f") ^
        2
    }
    irho = irho + 1
  }
  cv.mean = rowMeans(cv.err)
  rho.min = rho_factors[which.min(cv.mean)]
  fit = constrained_L2reg(X, Y, rho.min)
  list(rho.opt = rho.min, sigma2.opt = fit$sigma2)
}



## ridge regression for estimate sigma2 initialization in gene expression
#' @param M number of gene
#' @param N number of sample
#' @example
#' M = data$var$Ng
#' N = data$var$N
#' B = constrained_L2reg_multi(Xs, Ys, sigma2$rho.opt, M, N)
constrained_L2reg_multi = function(X, Ys, rho, M, N) {
  K = length(Ys)
  B = list()
  F = list()
  mu = list()
  err = 0
  df  = 0
  for (i in 1:K) {
    fit = constrained_L2reg(X, Ys[[i]], rho)
    B[[i]]  = as.matrix(fit$B)
    F[[i]]  = fit$F
    mu[[i]] = fit$mu
    err = err + fit$sigma2 * (N * M - 1)
    df  = df + (N * M - 1)
  }
  sigma2 = err / df
  list(
    B = B,
    F = F,
    sigma2 = sigma2,
    mu = mu
  )
}

## cross-validation on ridge regression to estimate sigma2
#' @param nrho number of L2 penalty's coefficient
#' @param ncv  number of cross-validation
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatX(data)
#' M  = data$var$Ng
#' N  = data$var$N
#' sigma2 = getsigma2_L2reg_multi(Xs, Ys, nrho = 20, M = M, N = N)
getsigma2_L2reg_multi = function(X,
                                 Ys,
                                 nrho = 10,
                                 ncv = 5,
                                 M,
                                 N) {
  rho_factors = 10 ** (seq(-6, 2, length.out = nrho))
  cv.err  = matrix(0, nrow = nrho, ncol = ncv)
  cv.fold = sample(seq(1, ncv), size = N, replace = T)
  irho    = 1
  for (rho in rho_factors) {
    for (cv in 1:ncv) {
      ytrain = lapply(Ys, function(y) {
        y[, cv.fold != cv, drop = F]
      })
      xtrain = lapply(X, function(x) {
        x[, cv.fold != cv, drop = F]
      })
      ytest  = lapply(Ys, function(y) {
        y[, cv.fold == cv, drop = F]
      })
      xtest  = lapply(X, function(x) {
        x[, cv.fold == cv, drop = F]
      })
      Ntrain = sum(cv.fold != cv)
      fit    = constrained_L2reg_multi(xtrain, ytrain, rho, M, Ntrain)
      for (k in 1:length(Ys)) {
        ftest  = lapply(1:M, function(i) {
          crossprod(fit$F[[k]][[i]], xtest[[i]])
        })
        ftest  = do.call(rbind, ftest)
        cv.err[irho, cv] = cv.err[irho, cv] + norm((diag(M) - fit$B[[k]]) %*% ytest[[k]] - ftest - fit$mu[[k]], type = "f")
      }
    }
    irho = irho + 1
  }
  cv.mean = rowMeans(cv.err)
  rho.min = rho_factors[which.min(cv.mean)]
  fit = constrained_L2reg_multi(X, Ys, rho.min, M, N)
  list(
    rho.opt = rho.min,
    sigma2.opt = fit$sigma2[1],
    cv.ram = list(rho = rho_factors, cvm = cv.mean)
  )
}


################### 
## ridge regression for estimate sigma2 initialization
## on different gene expressionand different eQTLs
#' @param M number of gene
#' @param N number of sample
#' @example
#' M = data$var$Ng
#' N = data$var$N
#' params.init = constrained_L2reg_gen(Xs, Ys, sigma2$rho.opt, M, N)
constrained_L2reg_gen = function(Xs, Ys, rho, M, N) {
  K = length(Ys)
  B = list()
  F = list()
  mu = list()
  err = 0
  df  = 0
  for (i in 1:K) {
    fit = constrained_L2reg(Xs[[i]], Ys[[i]], rho)
    B[[i]]  = as.matrix(fit$B)
    F[[i]]  = fit$F
    mu[[i]] = fit$mu
    err = err + fit$sigma2 * (N[i] * M - 1)
    df  = df + (N[i] * M - 1)
  }
  sigma2 = err / df
  list(
    B = B,
    F = F,
    sigma2 = sigma2,
    mu = mu
  )
}


## generalized cross-validation on ridge regression to estimate sigma2
## on different gene expressionand different eQTLs
#' @param nrho number of L2 penalty's coefficient
#' @param ncv  number of cross-validation
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatXs(data)
#' M  = data$var$Ng
#' N  = data$var$N
#' sigma2 = getsigma2_L2reg_gen(Xs, Ys, nrho = 20, M = M, N = N)
getsigma2_L2reg_gen = function(Xs,
                               Ys,
                               nrho = 10,
                               ncv = 5,
                               M,
                               N) {
  rho_factors = 10 ** (seq(-6, 2, length.out = nrho))
  cv.err  = matrix(0, nrow = nrho, ncol = ncv)
  cv.fold = list()
  cv.fold[[1]] = sample(seq(1, ncv), size = N[1], replace = T)
  cv.fold[[2]] = sample(seq(1, ncv), size = N[2], replace = T)
  irho    = 1
  for (rho in rho_factors) {
    for (cv in 1:ncv) {
      ytrain = lapply(1:2, function(ix) {
        Ys[[ix]][, cv.fold[[ix]] != cv, drop = F]
      })
      xtrain = lapply(1:2, function(ix) {
        lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] != cv, drop = F]
        })
      })
      ytest  = lapply(1:2, function(ix) {
        Ys[[ix]][, cv.fold[[ix]] == cv, drop = F]
      })
      xtest  = lapply(1:2, function(ix) {
        lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] == cv, drop = F]
        })
      })
      Ntrain = sapply(cv.fold, function(f){ sum(f != cv) })
      fit    = constrained_L2reg_gen(xtrain, ytrain, rho, M, Ntrain)
      for (k in 1:length(Ys)) {
        ftest  = lapply(1:M, function(i) {
          crossprod(fit$F[[k]][[i]], xtest[[k]][[i]])
        })
        ftest  = do.call(rbind, ftest)
        cv.err[irho, cv] = cv.err[irho, cv] + norm((diag(M) - fit$B[[k]]) %*% ytest[[k]] - ftest - fit$mu[[k]], type = "f")
      }
    }
    irho = irho + 1
  }
  cv.mean = rowMeans(cv.err)
  rho.min = rho_factors[which.min(cv.mean)]
  fit = constrained_L2reg_gen(Xs, Ys, rho.min, M, N)
  list(
    rho.opt = rho.min,
    sigma2.opt = fit$sigma2[1],
    cv.ram = list(rho = rho_factors, cvm = cv.mean)
  )
}
