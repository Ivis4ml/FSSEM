#### select subset of hyper-parameters within 1 standard error
############# new version
## ultility functions
## select top k lambda and rho tuple by user define for stability selection
subLambda_ss = function(cvparams,
                        type = c("err", "loglik"),
                        ntop = 10) {
  cvm = if (type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = sapply(cvparams$opt.hyperparams, `[`, 1),
    rho    = sapply(cvparams$opt.hyperparams, `[`, 2),
    cvmean = sapply(cvm, mean),
    cvsd   = sapply(cvm, sd)
  )
  cv.min     = which.min(cvfuns$cvmean)
  cv.1se     = cvfuns[cv.min, 3] + cvfuns[cv.min, 4]
  cvfun.1se  = cvfuns[cvfuns$cvmean < cv.1se, c(1, 2, 3, 4)]
  cvfun.1se[order(cvfun.1se[,3], decreasing = F)[1:10],c(1,2)]
}

####################################################
## stability selection on a class of parameters
##################
##' @title ss_fssem
## stability selection
ss_fssem = function(Bs,
                    fs,
                    Ys,
                    Xs,
                    sigma2,
                    Ng,
                    params = NULL,
                    wBs = inverse(Bs),
                    rB  = flinv(Bs),
                    maxit = 100,
                    acc = TRUE,
                    inertial = inertial_pars("lin"),
                    threshold = 1e-3,
                    sparse = FALSE,
                    use.strict = TRUE,
                    Nbootstrap = 100,
                    Nsample = 0.75,
                    verbose = 2)  {
  N = ncol(Ys[[1]])
  ss.fold = vector("list", Nbootstrap)
  i = 1
  while(i <= Nbootstrap) {
    subs = sample(seq(1, N), ceiling(N * Nsample), replace = F)
    ss.fold[[i]] = sort(subs)
    ss.fold[[i+1]] = setdiff(seq(1, N), subs)
    i = i + 2
  }
  ss.fit = list()
  N.ss = vector("list", 2)
  D.ss = NULL
  for (j in 1:nrow(params)) {
    lambda = params[j, 1]
    rho    = params[j, 2]
    for (i in 1:Nbootstrap) {
      Yss = lapply(Ys, function(Y) {
        Y[, ss.fold[[i]]]
      })
      Xss = list()
      for (k in 1:length(Xs)) {
        Xss[[k]] = lapply(Xs[[k]], function(X) {
          X[, ss.fold[[k]], drop = F]
        })
      }
      ss.fit[[i]] = genSML_iPALM(
        Bs = params.init$B,
        fs = params.init$F,
        Ys = Yss,
        Xs = Xss,
        sigma2 = params.init$sigma2[1],
        Ng = data$var$Ng,
        lambda = lambda,
        rho = rho,
        maxit = maxit,
        threshold = threshold,
        use.strict = use.strict,
        acc = acc,
        sparse = sparse
      )
      err2abs = ss.fit[[i]]$err
      if (is.null(N.ss[[1]])) {
        N.ss[[1]] = ifelse(abs(ss.fit[[i]]$B[[1]]) > err2abs, 1, 0)
      } else {
        N.ss[[1]] = N.ss[[1]] + ifelse(abs(ss.fit[[i]]$B[[1]]) > err2abs, 1, 0)
      }
      if (is.null(N.ss[[2]])) {
        N.ss[[2]] = ifelse(abs(ss.fit[[i]]$B[[2]]) > err2abs, 1, 0)
      } else {
        N.ss[[2]] = N.ss[[2]] + ifelse(abs(ss.fit[[i]]$B[[2]]) > err2abs, 1, 0)
      }
      nonzero = c(as.numeric(ss.fit[[i]]$B[[1]]), as.numeric(ss.fit[[i]]$B[[2]]))
      nonzero = nonzero[nonzero != 0]
      thresh.2 = sort(abs(nonzero))[round(0.2 * length(nonzero)) + 1]
      if (is.null(D.ss)) {
        D.ss = ifelse(
          abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > pmin(abs(ss.fit[[i]]$B[[1]]), abs(ss.fit[[i]]$B[[2]])) &
            abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > thresh.2,
          1,
          0
        )
      } else {
        D.ss = D.ss + ifelse(
          abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > pmin(abs(ss.fit[[i]]$B[[1]]), abs(ss.fit[[i]]$B[[2]])) &
            abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > thresh.2,
          1,
          0
        )
      }
    }
  }
  list(N = N.ss, D = D.ss)
}





