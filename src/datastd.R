## ultility funcitons
center = function(X) {
  apply(X, 1, function(x) {
    x - mean(x)
  })
}

submatX = function(data) {
  ## submatrix for X on eQTL
  sk = data$obs$sk
  X  = data$obs$X
  lapply(sk, function(ix) {
    X[ix, , drop = F]
  })
}

## X(X^TX)^{-1}X^T
projection = function(X) {
  X %*% solve(crossprod(X)) %*% t(X)
}


## use QR more slowly
projection.QR = function(X) {
  qr = qr.default(X, LAPACK = TRUE)
  Q  = qr.qy(qr, diag(1, nrow = nrow(qr$qr), ncol = qr$rank))
  tcrossprod(Q)
}

## centeralized Y (gene expression) and X (eQTL quantitive)
centralize = function(X, Y) {
  meanX = lapply(X, rowMeans)
  meanY = rowMeans(Y)
  X = lapply(X, center)
  Y = center(Y)
  list(X = X,
       Y = Y,
       muX = meanX,
       muY = meanY)
}

###### fused lasso ########
## centeralized Ys (gene expression) and Xs (eQTL quantitive)
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = data$obs$X
centralize_mult = function(Xs, Ys) {
  meanX = lapply(Xs, rowMeans)
  meanY = lapply(Ys, rowMeans)
  Xs = lapply(Xs, center)
  Ys = lapply(Ys, center)
  list(X = Xs,
       Y = Ys,
       muX = meanX,
       muY = meanY)
}

#' @description build sumbatrix for eQTLs observation for subset of corresponding genes
submatXs = function(data) {
  sk = data$obs$sk
  Xs  = list(X1 = data$obs$X1, X2 = data$obs$X2)
  lapply(Xs, function(X) {
    lapply(sk, function(ix) {
      X[ix, , drop = F]
    })
  })
}

## centralized for multiple Ys and Xs
#' @description generalized centralization of Ys and Xs
#'              Xs -> n x sk
#'              Ys -> n x ng
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatXs(data)
centralize_gen = function(Xs, Ys) {
  meanX = lapply(Xs, function(X) {
    lapply(X, rowMeans)
  })
  meanY = lapply(Ys, rowMeans)
  Xs = lapply(Xs, function(X) {
    lapply(X, center)
  })
  Ys = lapply(Ys, center)
  list(X = Xs,
       Y = Ys,
       muX = meanX,
       muY = meanY)
}



