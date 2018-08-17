# data preprocessing for the LSCNC data in paper.
### download CCLE data: GSE33356
### GPL570 GPL6801;
library(GEOquery)
gse <- getGEO(filename="exp/GSE33356_family.soft.gz")
gse1 <- GEOquery:::parseGSEMatrix("exp/GSE33356-GPL570_series_matrix.txt.gz",destdir = "exp", AnnotGPL = FALSE, 
    getGPL = F)
gse2 <- GEOquery:::parseGSEMatrix("exp/GSE33356-GPL6801_series_matrix.txt.gz",destdir = "exp", AnnotGPL = FALSE, 
    getGPL = F)

platforms <- lapply(GSMList(gse), function(x) {Meta(x)$platform_id})

## platform GPL570 GPL6801
# exprs = gse1$eset@assayData$exprs
# pbids = rownames(exprs)
# exprs = log2(exprs)
exprlist = BiocGenerics::Filter(function(gsm) {Meta(gsm)$platform_id=='GPL570'},  GSMList(gse))
probes = Table(GPLList(gse)[[1]])$ID
ex.matrix = do.call('cbind', lapply(exprlist, function(x) {
  tab = Table(x)
  mymatch = match(probes, tab$ID_REF)
  tab$VALUE[mymatch]
}))
ex.matrix = apply(ex.matrix, 2, function(x) {as.numeric(as.character(x))})
ex.matrix = log2(ex.matrix)
rownames(ex.matrix) = as.character(probes)

snplist = BiocGenerics::Filter(function(gsm) {Meta(gsm)$platform_id=='GPL6801'},  GSMList(gse))
snpids  = as.character(Table(GPLList(gse)[[2]])$ID)
annot_SNP = read.csv("./exp/GenomeWideSNP_6.na29.annot.csv", comment.char = "#", sep = ",",stringsAsFactors = F)
## annot_SNP = read.csv("./exp/GRCh37_hg19_AffyID2rsnumbers.txt", sep = ",",stringsAsFactors = F)
annot_SNP = annot_SNP[, c(1, 2, 3, 4, 9, 10, 11)]
snp.hash = annot_SNP[,2]
names(snp.hash) = annot_SNP[,1]
## check duplicate
snp.hash = snp.hash[!duplicated(snp.hash)]
snp.matrix = do.call('cbind', lapply(snplist, function(x) {
  tab = Table(x)
  mymatch = match(snpids, tab$ID_REF)
  as.character(tab$VALUE[mymatch])
}))
filterID = intersect(snpids, names(snp.hash))
rownames(snp.matrix) = snpids
snp.matrix = snp.matrix[filterID,]
snp.matrix = snp.matrix[complete.cases(snp.matrix),]
snpids = rownames(snp.matrix)
rssnp  = snp.hash[snpids]
rownames(snp.matrix) = rssnp

snp.matrix[snp.matrix == "AA"] = 0
snp.matrix[snp.matrix == "AB"] = 1
snp.matrix[snp.matrix == "BB"] = 2
snp.matrix[snp.matrix == "NC"] = 3
nr = nrow(snp.matrix)
rname = rownames(snp.matrix)
cname = colnames(snp.matrix)
snp.matrix = matrix(as.numeric(snp.matrix), nrow = nr)
rownames(snp.matrix) = rname
colnames(snp.matrix) = cname
snp.matrix[snp.matrix == 3] = NA
nix = apply(snp.matrix, 1, function(x){
  length(unique(x)) == 1 & is.na(unique(x)[1])
})
snp.matrix = snp.matrix[!nix,]
snp.matrix = apply(snp.matrix, 1, function(x){
  ni = which(is.na(x))
  if(length(ni) > 0) {
    x[ni] = sample(x = x[!is.na(x)], size = length(ni), replace = T)
  }
  x
})
snp.matrix = t(snp.matrix)

library(org.Hs.eg.db)
library(clusterProfiler)

## 54675 affy probes
probes = as.character(probes)
annoAffy  = read.csv("./exp/HG-U133_Plus_2.na36.annot.csv", comment.char = "#", sep = ",",stringsAsFactors = F)
prob2gene = annoAffy$Entrez.Gene
names(prob2gene) = annoAffy$Probe.Set.ID

## interaction sample names
## 42 pairs of data
library(stringr)
pdata1 = phenoData(gse1$eset)
pdata2 = phenoData(gse2$eset)

sample1 = as.character(pdata1@data$title)
sample1 = str_sub(sample1, start = 13, end = -1)
sample2 = as.character(pdata2@data$title)
sample2 = str_extract(sample2, "[\\d]+[T|N]")
cosample = intersect(sample2, sample1)
subid1  = sapply(cosample, function(i){which(sample1 == i)})
subid2  = sapply(cosample, function(i){which(sample2 == i)})
pdata1  = pdata1@data[subid1,]
pdata2  = pdata2@data[subid2,]
ex.matrix = ex.matrix[,subid1]
snp.matrix = snp.matrix[,subid2]


library(hgu133plus2.db)
library(biomaRt)
p2txtable = toTable(hgu133plus2ENSEMBL)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
attributes = listAttributes(ensembl)
mapping = getBM(
  attributes = c("affy_hg_u133_plus_2", "ensembl_transcript_id", "entrezgene", "hgnc_symbol", "ensembl_peptide_id", "refseq_mrna", "ensembl_exon_id"),
  filters = "affy_hg_u133_plus_2",
  values = probes,
  mart = ensembl
)

## mapping probes to transcript remove nonspecific probes
mapping = mapping[complete.cases(mapping),]
ns_probes = names(prob2gene[!(str_detect(prob2gene, "///"))])
mapping = mapping[mapping$affy_hg_u133_plus_2 %in% ns_probes,]
mapix = sapply(1:nrow(mapping), function(ix){
  c = mapping[ix,]
  prob2gene[c[1,1]] == c[1,3]
})
mapping = mapping[mapix,]
mapping = mapping[mapping$ensembl_peptide_id != "",]
ex.smatrix = ex.matrix[ns_probes,]

## download cdf package of hgu133plus2hsentrezgcdf here: 
## http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133plus2hsentrezgprobe_22.0.0.tar.gz
library(hgu133plus2hsentrezgcdf)
library(affy)
filenames = paste0("exp/GSE33356/", colnames(ex.matrix), ".CEL.gz")
rawdata = ReadAffy(filenames = filenames)
rawdata@cdfName = "HGU133Plus2_Hs_ENSG"
rawdata@cdfName = "HGU133Plus2_Hs_ENTREZG"
data = rma(rawdata)
samples = sapply(str_split(colnames(data), "\\."), `[`, 1)
genesid = sapply(str_split(rownames(data), "_"), `[`, 1)
ex.smatrix = data@assayData$exprs
rownames(ex.smatrix) = genesid
colnames(ex.smatrix) = samples
ex.smatrix = ex.smatrix[,colnames(ex.matrix)]
snp.smatrix = snp.matrix

## MatrixEQTL detect eQTL
library(MatrixEQTL)
library("BSgenome.Hsapiens.UCSC.hg19")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

geneloc = getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  mart = ensembl
)
tmp = bitr(unique(rownames(ex.smatrix)), fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")
genenames = tmp$ENTREZID
names(genenames) = tmp$ENSEMBL
geneloc = geneloc[geneloc$ensembl_gene_id %in% names(genenames),]
geneloc$entrezgene = genenames[geneloc$ensembl_gene_id]
geneloc = geneloc[complete.cases(geneloc),]
chr   = sapply(geneloc$entrezgene, function(x){org.Hs.egCHR[[x]]})
geneloc$chromosome_name = sapply(chr, `[`, 1)
geneloc = geneloc[complete.cases(geneloc),]

ex.smatrix = ex.smatrix[as.character(unique(geneloc$entrezgene)),]
exprmat = SlicedData$new()
colnames(ex.smatrix) = NULL
exprmat$initialize(ex.smatrix)
snpmat  = SlicedData$new()
colnames(snp.smatrix) = NULL
snpmat$initialize(snp.smatrix)
geneloc = geneloc[as.character(geneloc$entrezgene) %in% rownames(exprmat),]
geneloc = geneloc[,c(5, 2, 3, 4)]
colnames(geneloc) = c("geneid", "chr", "left", "right")
rownames(geneloc) = NULL
snploc = annot_SNP[,c(2,3,4)]
snploc = unique(snploc)
colnames(snploc) = c("snpid", "chr", "pos")
snploc = snploc[snploc$snpid %in% rownames(snpmat),]
covs = NULL
status = as.character(pdata1$characteristics_ch1)
status[status == "tissue: lung cancer"] = "tumor"
status[status != "tumor"] = "normal"
age = as.numeric(str_sub(as.character(pdata1$characteristics_ch1.2), start = 6, end = 7))
covs = rbind(covs,  status = ifelse(status == "normal", 0, 1), age)
covmat = SlicedData$new()
covmat$initialize(covs)


## cancer type
## considering all tx or specific tx
status = as.character(pdata1$characteristics_ch1)
status[status == "tissue: lung cancer"] = "tumor"
status[status != "tumor"] = "normal"
statvar= model.matrix(~status)
library(limma)
fit = lmFit(ex.smatrix, statvar)
efit = eBayes(fit)
prvals = topTable(efit, coef = 2, number = Inf)
sig.genes =  rownames(prvals)[prvals$adj.P.Val <= 1e-2]

## eQTL analysis and preparation of MatrixEQTL data
snp.smatrix = snp.smatrix[snploc$snpid[snploc$chr %in% c(as.character(seq(1,24)), "X", "Y")], ]
SNPS = data.frame(id = rownames(snp.smatrix), snp.smatrix)
colnames(SNPS) = c("id", paste("Sample_", seq(1, 84), sep=""))
write.table(SNPS, "SNP.txt", quote = F, row.names = F, col.names = TRUE, sep = "\t")
GES  = data.frame(id = rownames(ex.smatrix), ex.smatrix)
colnames(GES) = c("id", paste("Sample_", seq(1, 84), sep=""))
write.table(GES, "GE.txt", quote = F, row.names = F, col.names = TRUE, sep = "\t")
CVs = data.frame(id = rownames(covs), covs)
colnames(CVs) = c("id", paste("Sample_", seq(1, 84), sep=""))
write.table(CVs, "Covariates.txt", quote = F, row.names = F, col.names = TRUE, sep = "\t")

## snp and genepos
saveRDS(snploc, "SNP.rds")
saveRDS(geneloc, "GE.rds")

## run MatrixEQTL
source("run_EQTL.R")

### eQTL data prepration
## read eQTL data
library(entropy)
snp.smatrix = snp.smatrix + 1
normal = which(status == "normal")
tumor  = which(status == "tumor")
filter_by_entropy = apply(snp.smatrix, 1, function(x) {
  entropy(table(x[tumor])) > 0.5 & entropy(table(x[normal])) > 0.5
})
filter_by_MAF = apply(snp.smatrix, 1, function(x) {
  map.t = min(table(x[tumor]) / length(tumor))
  map.n = min(table(x[normal]) / length(normal))
  (map.t >= 0.15 & map.n >= 0.15)
})

eqtl_genes = read.csv("/media/xinchou/Storage/EQTL/SML/cis_eQTL_results_R.txt", sep = "\t", stringsAsFactors = F)
## pick significant
sig_egene = eqtl_genes
filtered_ix = names(filter_by_MAF)[filter_by_MAF]
sig_index = apply(sig_egene, 1, function(x){x[1] %in% filtered_ix})
## sig.egene = eqtl_genes[eqtl_genes$FDR < 1e-2, ]
sig_egene = sig_egene[sig_index,]

eqtl_map = split(sig_egene$SNP, sig_egene$gene)

### remove equal SNP data
eqtl_map = lapply(eqtl_map, function(s) {
  s = rownames(unique(snp.smatrix[s,,drop = F]))
  mat1 = t(center(snp.smatrix[s,tumor, drop = F]))
  mat2 = t(center(snp.smatrix[s,normal, drop = F]))
  s1 = rownames(unique(mat1))
  s2 = rownames(unique(mat2))
  intersect(s1, s2)
})

eqtl_map = eqtl_map[sapply(eqtl_map, function(x){length(x) > 0})]

for(ge in names(eqtl_map)) {
  if(length(eqtl_map[[ge]]) > 1) {
    g = as.numeric(ge)
    fdr = NULL
    for(snp in eqtl_map[[ge]]) {
      fdr = c(fdr, sig_egene[sig_egene$SNP == snp & sig_egene$gene == g, 6])
    }
    min = which.min(fdr)
    eqtl_map[[ge]] = eqtl_map[[ge]][min]
  }
}


## build proper data structrure for FSSEM algorithm
getRealdata = function(log = F, gene_cand = NULL, ...) {
  gene = gene_cand
  normal = which(status == "normal")
  tumor  = which(status == "tumor")
  Y1 = ex.smatrix[gene, normal]
  Y2 = ex.smatrix[gene, tumor]
  eqtl = eqtl_map[gene]
  eqtlmat = snp.smatrix[unique(unlist(eqtl)),]
  N1 = ncol(Y1)
  N2 = ncol(Y2)
  Ng = length(gene)
  Nk = nrow(eqtlmat)
  G  = matrix(0, nrow = Ng, ncol = Nk)
  sk = lapply(1:Ng,
              function(i) {
                s = which(rownames(eqtlmat) %in% eqtl[[i]])
                G[i, s] <<- 1
                s
              })
  X1 = eqtlmat[,normal]
  X2 = eqtlmat[,tumor]
  if(log) {
    Y1 = log2(Y1)
    Y2 = log2(Y2)
  }
  data = list(
    obs = list(
      Y1 = Y1,
      Y2 = Y2,
      X1 = X1,
      X2 = X2,
      sk = sk
    ),
    var = list(
      Names = gene,
      Ns = c(N1, N2),
      Ng = Ng,
      Nk = Nk
    )
  )
}


### data of interest filter; get all high confidence edges
### Please donwload the `lung_top` from https://s3-us-west-2.amazonaws.com/humanbase/networks/lung_top.gz
lung_net = read.csv("lung_top", sep = "\t", header = F)
lung_net = lung_net[lung_net[,3] >= 0.80, ]
lung_net[,1] = as.character(lung_net[,1])
lung_net[,2] = as.character(lung_net[,2])
## candidate nodes
vtx_cand = names(eqtl_map)
cand_fileter = apply(lung_net, 1, function(e){ 
  e[1] %in% vtx_cand & e[2] %in% vtx_cand
})
lung_net = lung_net[cand_fileter, ]
rownames(lung_net) = NULL
lung_net = graph_from_data_frame(lung_net[,1:2], directed = F)
gene_cand = names(V(lung_net))

### Run LSCSN simulation
source("../src/solver.min.R")
require(igraph)
data = getRealdata(log = F, gene_cand = gene_cand)

Xs = submatXs(data)
Ys = data$obs[c("Y1", "Y2")]

sigma2.init = getsigma2_L2reg_gen(Xs,
                                  Ys,
                                  nrho = 50,
                                  M = data$var$Ng,
                                  N = data$var$Ns)
params.init = constrained_L2reg_gen(Xs, Ys, sigma2.init$rho.opt, M = data$var$Ng, N = data$var$Ns)

gen.params = cv_genSML(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    weighted = T,
    nlambda = 20,
    nrho = 20,
    threshold = 1e-4
)

cvlambda.opt = optimLambda_cv1(gen.params, se = F, type = "err")  
allfit = genSML_iPALM(
    Bs = params.init$B,
    fs = params.init$F,
    Ys = Ys,
    Xs = Xs,
    sigma2 = params.init$sigma2[1],
    Ng = data$var$Ng,
    lambda = cvlambda.opt$lambda,
    rho = cvlambda.opt$rho,
    maxit = 2000,
    threshold = 1e-4,
    use.strict = T,
    acc = TRUE,
    sparse = F
)

## run stability selection
ssfit = ssSML_iPALM(
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
  use.strict = T,
  acc = TRUE,
  sparse = F,
  Nbootstrap = 100,
  Nsample = 0.5
)

## modify network for comparison
idmap = bitr(data$var$Names, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
rownames(allfit$B[[1]]) = idmap[,2]
colnames(allfit$B[[1]]) = idmap[,2]
rownames(allfit$B[[2]]) = idmap[,2]
colnames(allfit$B[[2]]) = idmap[,2]
allfit$N = list()
err2abs = sqrt(allfit$err)
allfit$N[[1]] = ifelse((abs(allfit$B[[1]]) > err2abs) | t(abs(allfit$B[[1]]) > err2abs), 1, 0)
allfit$N[[2]] = ifelse((abs(allfit$B[[2]]) > err2abs) | t(abs(allfit$B[[2]]) > err2abs), 1, 0)
nonzero = c(as.numeric(allfit$B[[1]]), as.numeric(allfit$B[[2]]))
nonzero = nonzero[nonzero != 0]
thresh.2 = sort(abs(nonzero))[round(0.2 * length(nonzero))+1]
allfit$D = ifelse(abs(allfit$B[[1]] - allfit$B[[2]]) > pmin(abs(allfit$B[[1]]), abs(allfit$B[[2]])) & 
                  abs(allfit$B[[1]] - allfit$B[[2]]) > thresh.2, 1, 0)


## build differential network 
ssD = ssfit$D / 100
rownames(ssD) = idmap[,2]
colnames(ssD) = idmap[,2]
## cutoff for stability selection is defined at 0.80
diffedges = ifelse(ssD >= 0.8, 1, 0)
differentNet = ssD * diffedges
signDiff = ifelse(allfit$B[[2]] - allfit$B[[1]] > 0, 1, 2) * diffedges
dnet = graph_from_adjacency_matrix(differentNet, weighted = T)
dnet = graph_from_adjacency_matrix(signDiff, weighted = T)
V(dnet)$label.cex = .4
deg = (igraph::degree(dnet, mode = "all") + 1) * 4 + 1
V(dnet)$size = deg
V(dnet)$color = "grey90"
E(dnet)$width = 1
E(dnet)$arrow.size = .2

postscript("network.eps",  fonts=c("serif", "Palatino"))
plot.igraph(dnet,
    layout = norm_coords(layout.fruchterman.reingold(dnet), 
                          ymin = -0.5, ymax = 0.5, 
                          xmin = -0.5, xmax = 0.5), vertex.frame.color = "grey90",
    edge.color = ifelse(E(dnet)$weight == 1, "red", "dark blue"), edge.curved = .1)
dev.off()


