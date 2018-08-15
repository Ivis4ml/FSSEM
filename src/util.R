## ultility
#' @param x estimated
#' @param b background
#' @example
#' FPR(params.opt8$B[[2]], data$var$B2)
#' FDR(params.opt8$B[[1]], data$var$B1)
#' TPR(params.opt8$B[[1]], data$var$B1)
TPR = function(x, b, precision = 1e-3) {
  x = as.matrix(x)
  b = as.matrix(b)
  sum(abs(x[b != 0]) > precision) / sum(b != 0)
}
FPR = function(x, b, precision = 1e-3) {
  x = as.matrix(x)
  b = as.matrix(b)
  sum(b == 0 & abs(x) > precision) / sum(b == 0)
}
FDR = function(x, b, precision = 1e-3) {
  x = as.matrix(x)
  b = as.matrix(b)
  if(sum(x !=0 ) != 0) {
    sum(abs(x[b == 0]) > precision) / sum(x != 0)
  } else {
    0
  }
}

TPR2 = function(x1, x2, b1, b2, precision = 1e-3) {
  x1 = as.matrix(x1)
  b1 = as.matrix(b1)
  x2 = as.matrix(x2)
  b2 = as.matrix(b2)
  ib = which(b1 != b2)
  ix = which(abs(x1 - x2) > precision)
  sum((ix %in% ib)) / length(ib)
}

FPR2 = function(x1, x2, b1, b2, precision = 1e-3) {
  x1 = as.matrix(x1)
  b1 = as.matrix(b1)
  x2 = as.matrix(x2)
  b2 = as.matrix(b2)
  ib = which(b1 == b2)
  ix = which(abs(x1 - x2) > precision)
  (sum(ix %in% ib) / length(ib))
}

FDR2 = function(x1, x2, b1, b2, precision = 1e-3) {
  x1 = as.matrix(x1)
  b1 = as.matrix(b1)
  x2 = as.matrix(x2)
  b2 = as.matrix(b2)
  ib = which(b1 == b2)
  ix = which(abs(x1 - x2) > precision)
  if(length(ix) > 0) {
    sum((ix %in% ib)) / length(ix)
  } else {
    0 
  }
}
