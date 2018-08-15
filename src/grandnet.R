## generate random graph(DAG) of G genes and their correspond differential network
## B_1 & B_2
## @param Ng gene number nodes
## @param e Expected number of edges per node
## @param d Expected ratio of differential edges per node (0.1)
## @param dag DAG or not
require(igraph)
require(Matrix)
require(glmnet)
getrandDAG = function(Ng,
                      e,
                      dag = TRUE,
                      d = 0.1,
                      Bmin = 0.5,
                      Bmax = 1,
                      maxit = Ng * Ng) {
  B1 = matrix(0,
              nrow = Ng,
              ncol = Ng)
  Nc = Ng * Ng
  Ne = rbinom(1, Nc, e / (Ng - 1))
  ## iteration mark
  iter1 = 0
  while (sum(B1) < Ne & iter1 < maxit) {
    edge     = runif(1, min = 1, max = Nc)
    B1[edge] = TRUE
    if (dag) {
      g        = graph_from_adjacency_matrix(B1)
      B1[edge] = is.dag(g)
    }
    iter1 = iter1 + 1
  }
  B2  = B1
  nn  = which(B1 != 0)
  nz  = which(B1 == 0)
  Nd  = ceiling(Ne * d)
  Ndf = rbinom(1, Nd, 0.5)
  while (sum(abs(B1 - B2)) < Ndf) {
    edge     = sample(nn, 1)
    B2[edge] = FALSE
  }
  iter2 = 0
  while (sum(B2) < Ne & iter2 < maxit) {
    edge     = sample(nz, 1)
    B2[edge] = TRUE
    if (dag) {
      g        = graph_from_adjacency_matrix(B2)
      B2[edge] = is.dag(g)
    }
    iter2 = iter2 + 1
  }
  ne  = which(B1 & B2)
  n1  = which(B1 & !(B2))
  n2  = which(!(B1) & B2)
  B1[ne] = B2[ne] = runif(length(ne), min = Bmin, max = Bmax) * sample(c(-1, 1), length(ne), replace = T)
  B1[n1] = runif(length(n1), min = Bmin, max = Bmax) * sample(c(-1, 1), length(n1), replace = T)
  B2[n2] = runif(length(n2), min = Bmin, max = Bmax) * sample(c(-1, 1), length(n2), replace = T)
  
  if(iter1 < maxit & iter2 < maxit & any(B1 != B2)) {
    if(!dag) {
      detIB1 = det(diag(Ng) - B1)
      detIB2 = det(diag(Ng) - B2)
      if(abs(detIB1) > 1e-6 & abs(detIB2) > 1e-6){
        list(B1 = B1, B2 = B2)
      } else {
        NULL
      }
    } else {
      list(B1 = B1, B2 = B2)
    }
  } else {
    NULL
  }
}


#' @param N   number of sample
#' @param Ng  number of gene
#' @param k   number of eQTL
getrandeQTL = function(N, Ng, Nk) {
  step = Nk / Ng
  X = round(2 * matrix(runif(Nk * N), nrow = Nk)) + 1
  G = matrix(0,
             nrow = Ng,
             ncol = Nk)
  ix = lapply(1:Ng,
              function(i) {
                s = seq(0, step - 1) * Ng + i
                G[i, s] <<- 1
                s
              })
  list(G = G, X = X, sk = ix)
}

## randNetinit
## randomly generate regulatory network for fixed seed
randNetinit = function(Ng = 10,
                       Nk = 10,
                       r = 0.3,
                       d = 0.1,
                       ...) {
  B = getrandDAG(Ng,
                 e = Ng * r,
                 dag = dag,
                 d = d,
                 ...)
  while (is.null(B)) {
    B = getrandDAG(Ng,
                   e = Ng * r,
                   dag = dag,
                   d = d,
                   ...)
  }
  B
}

require(mvtnorm)
getrandsem = function(N = 200,
                      Ng = 10,
                      Nk = 10,
                      r = 0.3,
                      d = 0.1,
                      dag = TRUE,
                      sigma = 0.1,
                      B = NULL,
                      ...) {
  if (is.null(B)) {
    B = getrandDAG(Ng,
                   e = Ng * r,
                   dag = dag,
                   d = d,
                   ...)
    while (is.null(B)) {
      B = getrandDAG(Ng,
                     e = Ng * r,
                     dag = dag,
                     d = d,
                     ...)
    }
  }
  Q = getrandeQTL(N, Ng, Nk)
  F = Q[[1]]
  X = Q[[2]]
  sk = Q[[3]]
  E1 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  E2 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  Y1 = solve(diag(Ng) - B[[1]]) %*% (F %*% X + E1)
  Y2 = solve(diag(Ng) - B[[2]]) %*% (F %*% X + E2)
  list(
    obs = list(
      Y1 = Y1,
      Y2 = Y2,
      X = X,
      sk = sk
    ),
    var = list(
      B1 = Matrix(B[[1]], sparse = T),
      B2 = Matrix(B[[2]], sparse = T),
      F = Matrix(F, sparse = T),
      N = N,
      Ng = Ng,
      Nk = Nk
    )
  )
}

## data = getrandsem(N = 200, Ng = 30, Nk = 90, r = 0.1, d = 0.1, sigma = 1, dag = TRUE)
## datn = getrandsem(N = 200, Ng = 30, Nk = 90, r = 0.1, d = 0.1, sigma = 1, dag = FALSE)

## eQTL and gene regulatory network are all different with each other under different conditions
#' @description getdiffeQTL
#' @param N   number of sample
#' @param Ng  number of gene
#' @param k   number of eQTL
#' @param d   differential ratio = 0.1
getdiffeQTL = function(N, Ng, Nk, d = 0.1) {
  step = Nk / Ng
  X = vector("list", 2)
  X[[1]] = round(2 * matrix(runif(Nk * N), nrow = Nk)) + 1
  X[[2]] = apply(X[[1]], 2, function(x) {
    Nd = Nk * d
    dx = sample(1:Nk, Nd, replace = F)
    x[dx] = round(2 * runif(Nd)) + 1
    x
  })
  G = matrix(0,
             nrow = Ng,
             ncol = Nk)
  ix = lapply(1:Ng,
              function(i) {
                s = seq(0, step - 1) * Ng + i
                G[i, s] <<- 1
                s
              })
  list(G = G, X = X, sk = ix)
}


#' @description getdiffsem
#' @details eQTL measurement for different condition are generated with proportional difference.
#' @param f difference proportion of each gene's eQTL measurement, such as SNP
getdiffsem = function(N = 200,
                      Ng = 10,
                      Nk = 10,
                      r = 0.3,
                      d = 0.1,
                      f = 0.1,
                      dag = TRUE,
                      sigma = 0.1) {
  B = getrandDAG(Ng, e = Ng * r, dag = dag, d = d)
  Q = getdiffeQTL(N, Ng, Nk, f)
  F = Q[[1]]
  X = Q[[2]]
  sk = Q[[3]]
  E1 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  E2 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  Y1 = solve(diag(Ng) - B[[1]]) %*% (F %*% X[[1]] + E1)
  Y2 = solve(diag(Ng) - B[[2]]) %*% (F %*% X[[2]] + E2)
  
  list(
    obs = list(
      Y1 = Y1,
      Y2 = Y2,
      X1 = X[[1]],
      X2 = X[[2]],
      sk = sk
    ),
    var = list(
      B1 = Matrix(B[[1]], sparse = T),
      B2 = Matrix(B[[2]], sparse = T),
      F = Matrix(F, sparse = T),
      N = N,
      Ng = Ng,
      Nk = Nk
    )
  )
}

## data = getdiffsem(N = 200, Ng = 30, Nk = 90, r = 0.1, d = 0.1, f = 0.1, sigma = 1, dag = TRUE)
