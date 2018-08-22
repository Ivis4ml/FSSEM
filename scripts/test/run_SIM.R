source("publication.R")
require(latex2exp)
## Large system simulation running on pegasus system
simuGenerator = function(method = c("validLasso", "validSML"), type = c("DAG", "DCG"), f = 1) {
  params.opt = list(
    N = seq(100, 1000, by = 100),
    Ng = c(10, 30),
    Ne = c(1, 3),
    sigma = c(sqrt(0.01), sqrt(0.09), sqrt(0.25))
  )
  method = match.arg(method)
  type   = match.arg(type)
  dag    = ifelse(type == "DAG", 1, 0)
  head = sprintf(
    "#!/bin/bash
#BSUB -o %%J.out
#BSUB -e %%J.err
#BSUB -P cai
#BSUB -W %s
#BSUB -q %s
#BSUB -n %d
#BSUB -B
#BSUB -N
#BSUB -u xxz220@miami.edu
    
cd .
module load R/3.3.1
    ",
    "100:00",
    "general",
    1
  )
  subcmd = c("#!/bin/bash")
  system(sprintf("mkdir test%s/%s_%d", type, method, f))
  system(sprintf("cp solver.min.R test%s/%s_%d", type, method, f))
  system(sprintf("cp %s.R test%s/%s_%d", method, type, method, f))
  cmds = NULL
  ## nonrepeat init seed
  init = sample.int(1e8, size = 2 * 2 * 3, replace = F)
  ix = 1
  for (g in params.opt$Ng) {
    for (k in params.opt$Ne) {
      for (sigma in params.opt$sigma) {
        seed = init[ix]
        ix = ix + 1
        for (n in params.opt$N) {
          job = sprintf("test%s/%s_%d/n%dg%de%ds%.2f.job",
                        type,
                        method,
                        f,
                        n,
                        g,
                        k,
                        sigma)
          cmd = sprintf(
            "Rscript %s.R -n %d -p %d -k %d -s %f -a %s -d %s > n%dg%de%ds%.2f.log \n",
            method,
            n,
            g,
            k,
            sigma,
            dag,
            seed,
            n,
            g,
            k,
            sigma
          )
          cmds = c(cmds, cmd)
          context = paste0(c(head,
                             cmd), collapse = "\n")
          writeLines(context, job, sep = '')
          subcmd = c(subcmd, sprintf(
            "bsub < %s",
            sprintf("n%dg%de%ds%.2f.job", n, g, k, sigma)
          ))
        }
      }
    }
  }
  writeLines(paste0(cmds, collapse = "\n"),
             sprintf("test%s/%s_%d/run.sh", type, method, f))
  writeLines(subcmd, sprintf("test%s/%s_%d/submit.sh", type, method, f), sep = "\n")
  system(sprintf("chmod +x test%s/%s_%d/submit.sh", type, method, f))
}


## Large system simulation processing on pegasus system
simuProcess = function(dir = "SMLg10", Ng = 10) {
  params.opt = list(
    N = seq(100, 500, by = 100),
    Ne = c(1, 3),
    sigma = c(0.01, 0.1, 1)
  )
  
  errCV = NULL
  errBIC = vector("list", 3)
  names(errBIC) = c(0, 0.5, 1)
  for (k in params.opt$Ne) {
    for (sigma in params.opt$sigma) {
      for (n in params.opt$N) {
        cat(sprintf("%s/n%dg%de%ds%.2f.rds", dir, n, Ng, k, sigma), "\n")
        sdat = readRDS(sprintf("%s/n%dg%de%ds%.2f.rds", dir, n, Ng, k, sigma))
        errCV = rbind(errCV, c(n, k, sigma, colMeans(sdat$error.cv[complete.cases(sdat$error.cv), ])))
        errBIC[["0"]]   = rbind(errBIC[["0"]], c(n, k, sigma, colMeans(sdat$error.eBIC$eBIC0[complete.cases(sdat$error.eBIC$eBIC0), ])))
        errBIC[["0.5"]] = rbind(errBIC[["0.5"]], c(n, k, sigma, colMeans(sdat$error.eBIC$eBIC.5[complete.cases(sdat$error.eBIC$eBIC.5), ])))
        errBIC[["1"]]   = rbind(errBIC[["1"]], c(n, k, sigma, colMeans(sdat$error.eBIC$eBIC1[complete.cases(sdat$error.eBIC$eBIC1), ])))
      }
    }
  }
  colnames(errCV) = c("N", "eQTLs", "sigma2", "TPR", "FDR", "TPR2", "FDR2")
  colnames(errBIC[["0"]]) = c("N", "eQTLs", "sigma2", "TPR", "FDR", "TPR2", "FDR2")
  colnames(errBIC[["0.5"]]) = c("N", "eQTLs", "sigma2", "TPR", "FDR", "TPR2", "FDR2")
  colnames(errBIC[["1"]]) = c("N", "eQTLs", "sigma2", "TPR", "FDR", "TPR2", "FDR2")
  
  errCV = as.data.frame(errCV)
  errBIC = lapply(errBIC, as.data.frame)
  
  list(cv = errCV, eBIC = errBIC)
}

res$cv$id = paste("eQtl=", res$cv$eQTLs, ", sigma=", res$cv$sigma2, sep = '')
ggplot(data = res$cv, aes(x = N, y = TPR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res$cv, aes(x = N, y = FDR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res$cv, aes(x = N, y = TPR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res$cv, aes(x = N, y = FDR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()

res$eBIC$`1`$id = paste("eQtl=", res$eBIC$`1`$eQTLs, ", sigma=", res$eBIC$`1`$sigma2, sep = '')
ggplot(data = res$eBIC$`1`, aes(x = N, y = TPR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res$eBIC$`1`, aes(x = N, y = FDR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res$eBIC$`1`, aes(x = N, y = TPR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res$eBIC$`1`, aes(x = N, y = FDR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()




res2$cv$id = paste("eQtl=", res2$cv$eQTLs, ", sigma=", res2$cv$sigma2, sep = '')
ggplot(data = res2$cv, aes(x = N, y = TPR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$cv, aes(x = N, y = FDR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$cv, aes(x = N, y = TPR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$cv, aes(x = N, y = FDR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()

res2$eBIC$`1`$id = paste("eQtl=", res2$eBIC$`1`$eQTLs, ", sigma=", res2$eBIC$`1`$sigma2, sep = '')
ggplot(data = res2$eBIC$`1`, aes(x = N, y = TPR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~http://r.789695.n4.nabble.com/Subscript-td877549.html id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$eBIC$`1`, aes(x = N, y = FDR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$eBIC$`1`, aes(x = N, y = TPR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$eBIC$`1`, aes(x = N, y = FDR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()


res2$eBIC$`0`$id = paste("eQtl=", res2$eBIC$`0`$eQTLs, ", sigma=", res2$eBIC$`0`$sigma2, sep = '')
ggplot(data = res2$eBIC$`0`, aes(x = N, y = TPR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$eBIC$`0`, aes(x = N, y = FDR)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$eBIC$`0`, aes(x = N, y = TPR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()
ggplot(data = res2$eBIC$`0`, aes(x = N, y = FDR2)) + geom_point(color = "black", alpha = 0.5) + geom_line(color = "blue") +
  facet_wrap( ~ id) + theme_Publication() + ylim(c(0, 1)) +
  scale_fill_Publication() + scale_colour_Publication()




### batch create DAG
subcmd = c("#!/bin/bash")
for(i in 1:50) {
  simuGenerator("validLasso", "DAG", i)
  simuGenerator("validSML", "DAG", i)
  subcmd = c(subcmd, sprintf("cd ./validLasso_%d\n./submit.sh\ncd ..", i), sprintf("cd ./validSML_%d\n./submit.sh\ncd ..", i))
}

writeLines(subcmd, "testDAG/submit_all.sh", sep = "\n")

### batch create for DCG
subcmd = c("#!/bin/bash")
for(i in 1:10) {
  simuGenerator("validLasso", "DCG", i)
  simuGenerator("validSML", "DCG", i)
  subcmd = c(subcmd, sprintf("cd ./validLasso_%d\n./submit.sh\ncd ..", i), sprintf("cd ./validSML_%d\n./submit.sh\ncd ..", i))
}

writeLines(subcmd, "testDCG/submit_all.sh", sep = "\n")


#### generation statistic
## Large system simulation processing on pegasus system
preProcessim = function(dir = "testhat", Ng = 10) {
  params.opt = list(
    N = seq(100, 1000, by = 100),
    Ne = c(1, 3),
    sigma = c(sqrt(0.01), sqrt(0.09), sqrt(0.25))
  )
  
  cvLasso = NULL
  cvSML   = NULL
  for (k in params.opt$Ne) {
    for (sigma in params.opt$sigma) {
      for (n in params.opt$N) {
        cat(sprintf("%s/n%dg%de%ds%.2f.rds", dir, n, Ng, k, sigma), "\n")
        params = sprintf("%d_%d_%d_%.2f", n, Ng, k, sigma)
        for(i in 1:10) {
          rname = sprintf("%s/validLasso_%s/n%dg%de%ds%.2f.rds", dir, i, n, Ng, k, sigma)
          if(file.exists(rname)) {
            res = readRDS(rname)
            if(is.null(cvLasso[[params]])) {
              cvLasso[[params]] = res[1:3]
            } else {
              cvLasso[[params]][[1]] = rbind(cvLasso[[params]][[1]], res[[1]])
              cvLasso[[params]][[2]] = rbind(cvLasso[[params]][[2]], res[[2]])
            }
          }
        }
        for(i in 1:10) {
          rname = sprintf("%s/validSML_%s/n%dg%de%ds%.2f.rds", dir, i, n, Ng, k, sigma)
          if(file.exists(rname)) {
            res = readRDS(rname)
            if(is.null(cvSML[[params]])) {
              cvSML[[params]] = res[1:3]
            } else {
              cvSML[[params]][[1]] = rbind(cvSML[[params]][[1]], res[[1]])
              cvSML[[params]][[2]] = rbind(cvSML[[params]][[2]], res[[2]])
            }
          }
        }
      }
    }
  }
  
  list(Lasso = cvLasso, SML = cvSML)
}


ggplotsim = function(dat = NULL, Ng = 10) {
  params.opt = list(
    N = seq(100, 1000, by = 100),
    Ne = c(1, 3),
    sigma = c(sqrt(0.01), sqrt(0.09), sqrt(0.25))
  )
  
  res = NULL
  for (k in params.opt$Ne) {
    for (sigma in params.opt$sigma) {
      for (n in params.opt$N) {
        type = ifelse(sigma < 1, "error.se", "error.cv")
        params = sprintf("%d_%d_%d_%.2f", n, Ng, k, sigma)
        data = dat$Lasso[[params]][[type]]
        data[is.nan(data)] = 0
        res = rbind(res, c(n, Ng, k, sigma, colMeans(data), 0))
        data = dat$SML[[params]][[type]]
        data[is.nan(data)] = 0
        res = rbind(res, c(n, Ng, k, sigma, colMeans(data), 1))
      }
    }
  }
  res = as.data.frame(res)
  colnames(res) = c("N", "Ng", "eQTLs", "sigma2", "TPR", "FDR", "TPR2", "FDR2", "method")
  res$method = ifelse(res$method == 0, "SML", "FSSEM")
  res$id = paste("n==", res$Ng, "~n[e]==", res$eQTLs, "~sigma^2==", (res$sigma2)^2, sep = '')
  res
}

library(latex2exp)
res = ggplotsim(preProcessim("testhat/testDAG", Ng = 10), Ng = 10)
res = res[res$sigma2 != 0.3,]

override.shape = c(16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21)
override.linetype = c()

ggplot(data = res, aes(x = N, y = TPR, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() +
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) +
  ylab(TeX('PD of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDB_10_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDB_10_DAG.pdf /media/xinchou/Code/Plot_Figure/PDB_10_DAG.svg")
ggplot(data = res, aes(x = N, y = FDR, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRB_10_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRB_10_DAG.pdf /media/xinchou/Code/Plot_Figure/FDRB_10_DAG.svg")
ggplot(data = res, aes(x = N, y = TPR2, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDdB_10_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDdB_10_DAG.pdf /media/xinchou/Code/Plot_Figure/PDdB_10_DAG.svg")
ggplot(data = res, aes(x = N, y = FDR2, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRdB_10_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRdB_10_DAG.pdf /media/xinchou/Code/Plot_Figure/FDRdB_10_DAG.svg")

res = ggplotsim(preProcessim("testhat/testDAG", Ng = 30), Ng = 30)
res = res[res$sigma2 != 0.3,]
ggplot(data = res, aes(x = N, y = TPR, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() +
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDB_30_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDB_30_DAG.pdf /media/xinchou/Code/Plot_Figure/PDB_30_DAG.svg")
ggplot(data = res, aes(x = N, y = FDR, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRB_30_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRB_30_DAG.pdf /media/xinchou/Code/Plot_Figure/FDRB_30_DAG.svg")
ggplot(data = res, aes(x = N, y = TPR2, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDdB_30_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDdB_30_DAG.pdf /media/xinchou/Code/Plot_Figure/PDdB_30_DAG.svg")
ggplot(data = res, aes(x = N, y = FDR2, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRdB_30_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRdB_30_DAG.pdf /media/xinchou/Code/Plot_Figure/FDRdB_30_DAG.svg")

####################
require(gridExtra)
res0 = res[res$sigma2 == 0.5 & res$Ng == 30 & res$eQTLs == 3,]
ggplot(data = res0, aes(x = N, y = TPR, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) +
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0.7, 1)) + 
  ylab(TeX('PD of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PD_30_3_0.5_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PD_30_3_0.5_DAG.pdf /media/xinchou/Code/Plot_Figure/PD_30_3_0.5_DAG.svg")
ggplot(data = res0, aes(x = N, y = FDR, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) + 
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0, 0.3)) + 
  ylab(TeX('FDR of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDR_30_3_0.5_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDR_30_3_0.5_DAG.pdf /media/xinchou/Code/Plot_Figure/FDR_30_3_0.5_DAG.svg")
ggplot(data = res0, aes(x = N, y = TPR2, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) + 
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0.7, 1)) + 
  ylab(TeX('PD of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDDB_30_3_0.5_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDDB_30_3_0.5_DAG.pdf /media/xinchou/Code/Plot_Figure/PDDB_30_3_0.5_DAG.svg")
ggplot(data = res0, aes(x = N, y = FDR2, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) + 
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRDB_30_3_0.5_DAG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRDB_30_3_0.5_DAG.pdf /media/xinchou/Code/Plot_Figure/FDRDB_30_3_0.5_DAG.svg")



res = ggplotsim(preProcessim("testhat/testDCG", Ng = 10), Ng = 10)
res = res[res$sigma2 != 0.3,]
ggplot(data = res, aes(x = N, y = TPR, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() +
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDB_10_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDB_10_DCG.pdf /media/xinchou/Code/Plot_Figure/PDB_10_DCG.svg")
ggplot(data = res, aes(x = N, y = FDR, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\mathbf{B}^{(1)}$ and $\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRB_10_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRB_10_DCG.pdf /media/xinchou/Code/Plot_Figure/FDRB_10_DCG.svg")
ggplot(data = res, aes(x = N, y = TPR2, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDdB_10_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDdB_10_DCG.pdf /media/xinchou/Code/Plot_Figure/PDdB_10_DCG.svg")
ggplot(data = res, aes(x = N, y = FDR2, colour = method)) + geom_point(aes(shape = method), size = 2) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(16, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRdB_10_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRdB_10_DCG.pdf /media/xinchou/Code/Plot_Figure/FDRdB_10_DCG.svg")

##########################



res = ggplotsim(preProcessim("testhat/testDCG", Ng = 30), Ng = 30)
res = res[res$sigma2 != 0.3,]
ggplot(data = res, aes(x = N, y = TPR, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line() +
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDB_30_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDB_30_DCG.pdf /media/xinchou/Code/Plot_Figure/PDB_30_DCG.svg")
ggplot(data = res, aes(x = N, y = FDR, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRB_30_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRB_30_DCG.pdf /media/xinchou/Code/Plot_Figure/FDRB_30_DCG.svg")
ggplot(data = res, aes(x = N, y = TPR2, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('PD of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDdB_30_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDdB_30_DCG.pdf /media/xinchou/Code/Plot_Figure/PDdB_30_DCG.svg")
ggplot(data = res, aes(x = N, y = FDR2, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line() + 
  facet_wrap( ~ id, nrow = 2, labeller = label_parsed) + theme_Publication() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRdB_30_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRdB_30_DCG.pdf /media/xinchou/Code/Plot_Figure/FDRdB_30_DCG.svg")

##########
res0 = res[res$sigma2 == 0.5 & res$Ng == 30 & res$eQTLs == 3,]
ggplot(data = res0, aes(x = N, y = TPR, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) +
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0.7, 1)) + 
  ylab(TeX('PD of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PD_30_3_0.5_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PD_30_3_0.5_DCG.pdf /media/xinchou/Code/Plot_Figure/PD_30_3_0.5_DCG.svg")
ggplot(data = res0, aes(x = N, y = FDR, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) + 
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0, 0.3)) + 
  ylab(TeX('FDR of $\\,\\mathbf{B}^{(1)}\\,$ and $\\,\\mathbf{B}^{(2)}$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDR_30_3_0.5_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDR_30_3_0.5_DCG.pdf /media/xinchou/Code/Plot_Figure/FDR_30_3_0.5_DCG.svg")
ggplot(data = res0, aes(x = N, y = TPR2, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) + 
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0.7, 1)) + 
  ylab(TeX('PD of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/PDDB_30_3_0.5_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/PDDB_30_3_0.5_DCG.pdf /media/xinchou/Code/Plot_Figure/PDDB_30_3_0.5_DCG.svg")
ggplot(data = res0, aes(x = N, y = FDR2, colour = method)) + geom_point(aes(shape = method), size = 5) + geom_line(size = 1) + 
  facet_wrap( ~ id) + theme_Bioinformatics() + ylim(c(0, 1)) + 
  ylab(TeX('FDR of $\\,|\\mathbf{B}^{(2)} - \\mathbf{B}^{(1)}|$')) +
  scale_fill_Publication() + scale_colour_Publication() + scale_shape_manual(values = c(2, 4))
ggsave(filename = "/media/xinchou/Code/Plot_Figure/FDRDB_30_3_0.5_DCG.svg", device = "svg", width = 6, height = 6)
system("rsvg-convert -f pdf -o /media/xinchou/Code/Plot_Figure/FDRDB_30_3_0.5_DCG.pdf /media/xinchou/Code/Plot_Figure/FDRDB_30_3_0.5_DCG.svg")

saveRDS(list(pdB, fdrB, pdDB, fdrDB), file = "exp/DCG30Fig.rds")

