# test RNA-seq datasets
rm(list=ls())
setwd("D:/001big/001Batch/test/GOF_supp")
#1. load package----------------------------------------------------------------
#Load required packages
library(foreach)
reqpkg <- c("phyloseq", "ggplot2", "parallel", "nleqslv", "MASS", "fdrtool", 
            "reshape2", "xtable", "edgeR", "DESeq2", "muStat","EstimationTools")
for(i in reqpkg)
{library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, 
         character.only = TRUE)
};rm(i)
#  necessary functions
for(j in list.files("R")){source(file.path("R",j))};rm(j) #Smooth test
for(j in list.files("auxFun")){source(file.path("auxFun",j))};rm(j)#Auxiliary functions

#2. neuroblastoma dataset-------------------------------------------------------
##2.1 RNA-Seq dataset of neuroblastoma cell line data, treated by nutlin or ethanol.
neuroblastoma0 = readRDS("Datasets/celine_neuroblastoma_data.RData")
neuroblastoma = phyloseq(otu_table(neuroblastoma0$counts[neuroblastoma0$mRNA,], 
                                   taxa_are_rows = TRUE), 
                         sample_data(data.frame(group = neuroblastoma0$group, 
                                                row.names = names(neuroblastoma0$counts))))
neuroblastoma_count <- neuroblastoma0$counts[neuroblastoma0$mRNA,]

### estimate the dispersion-------------------------------------------
neuroblastoma_NB <- estNBparamsPhylo(neuroblastoma,"group")
dispList = neuroblastoma_NB$theta
odDF <- melt(dispList)
Dispersion  <- 1/odDF
mean(Dispersion[,1])

#physeq <- neuroblastoma
#phyVars <- "group"
estAndSave(neuroblastoma,"group","neuroblastoma") # estimate the dispersion and save in Rdata



# 保存parameter到Rdata
a <- load('HPCbootstrap/params/neuroblastomaparams.RData')
# define, reorder, save
rawNames <- "neuroblastoma"
tags = factor(x=rawNames,labels=c("Neuroblastoma (cell line)"),
              levels = rawNames, ordered = TRUE)
vars <- c("group")
Origin <- factor("RNA-seq",levels="RNA-seq",labels = "RNA-seq",ordered = TRUE)
phyList <- list(neuroblastoma)
seqTechs <- factor("Illumina")
names(phyList) =  names(Origin) = names(seqTechs) = names(tags) = tags
nSamples = sapply(phyList, nsamples)
paramList <- lapply(rawNames, function(tag){load(file.path("HPCbootstrap/params",paste0(tag, "params.RData")));
  get(paste0("params", tag))}) # parameter of NB
save(tags, nSamples, phyList, vars, seqTechs, Origin, paramList, rawNames, file ="tags_neuroblastoma.RData")
rm(neuroblastoma)

#2.2 boxplot of dispersion parameters------------------------------------------
load("tags_neuroblastoma.RData")
dispList = lapply(paramList, function(x){x$theta}) # theta parameter
names(dispList) = tags
odDF = melt(dispList)
names(odDF) = c("Dispersion","Dataset")
odDF$Dataset = factor(odDF$Dataset, levels = tags, labels = tags, ordered = TRUE)
odDF$Origin = Origin[odDF$Dataset]
odDFmerge = merge(odDF, data.frame("SequencingTechnology" = seqTechs, Dataset = levels(odDF$Dataset)))
dispPlot = ggplot(aes(y = 1/Dispersion, x =  Dataset, fill = Origin, col = SequencingTechnology), data = odDFmerge) + geom_boxplot(size = 0.4, outlier.size = 0.7) + ylab("Estimated dispersion parameter")  + theme_bw() +  scale_y_continuous(trans ="log10", limits = c(1e-8,1e3), breaks = 10^seq(-8,4,by = 1))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + scale_fill_brewer(palette = "Accent", name = "Origin") + scale_colour_brewer(palette = "Dark2", name = "Sequencing \ntechnology") + geom_hline(yintercept = 10, col = "red", linetype ="dashed") + xlab("")
dispPlot
rm(dispPlot, odDFmerge, odDF, dispList)

#2.3 failed to fit--------------------------------------------------------------
# 2.3.1 Barplots of proportions of features for which the NB distribution failed to fit.
#estTaxa = sapply(paramList, function(x){sum(!is.na(x$coef[1,]))})
#presentTaxa = sapply(phyList, ntaxa)
#nonTrimmedTaxa = sapply(phyList, function(phy){sum(colMeans(otu_table(phy)@.Data==0)<1)})
#dfTrim = data.frame(Dataset = tags, failFrac = 1-estTaxa/nonTrimmedTaxa, Origin = Origin)
#ggplot(aes(x = Dataset, y = failFrac, fill = Origin), data = dfTrim) + geom_bar(stat ="identity")
#+ scale_y_continuous(limits = c(0,0.06), name = "Proportion of failed fits")
#+ theme_bw() + theme(axis.text.x = element_text(angle = 90)) + xlab("")
# 2.3.2  smooth test------------------------------------------------------------
testResultsList = smoothTestPhylo(neuroblastoma,"group")# 进行平滑检验
# 2.3.3 生成基于负二项分布的bootstrap样本----------------------------------------
testbootstrap <- genBootStrap(neuroblastoma_NB,reps = 5, 
                              single = TRUE, xFac = NULL, filePath = NULL, nCores = 1)
# 参数向量params、样本重复次数reps、是否使用单一抽样single、
#是否需要对因子变量进行因子化xFac、结果文件保存路径filePath和计算核心数nCores

#2.3.4 使用平滑检验估计bootsrap 样本参数, bootstrap test statistics-------------
mat <- neuroblastoma0$counts[neuroblastoma0$mRNA,]
testBootRes <- smoothTestMat(t(mat),testbootstrap$X)

#2.3.5 select the proportion of poorly fit features-----------------------------
# 对于每个参数计算其实际结果是否小于等于 bootstrap 重采样结果的比例。
# 最终返回所有标签的 p-value 列表,从而得出比例


B = 1000
filePath = "HPCbootstrap/testBootstrap"
reps = 5
foo = mclapply(tags, mc.cores = 1 , mc.preschedule = FALSE, function(tag){
  filPat = paste0(filePath, "/", tag)
  dir.create(filPat)
  #Generate bootstrapped datasets per feature.
  genBoot = genBootStrap(paramList[[tag]], reps = 5, filePath = filPat)
  assign("X", genBoot$X) # Also save the design matrix
  save(X, file = file.path(filPat, "X.RData"))
})

foo3 = mclapply(tags, mc.cores = 1 , mc.preschedule = FALSE, function(tag){
  #Create a folder to save the results
  dir.create(paste0(filePath, "/", tag, "/results"))})

filPat = paste0(filePath, "/", tags)
dir.create(paste0(filPat, "/testResults"))
sapply(seq_len(reps), function(i){
  if(!file.exists(paste0(filPat, "/testResults/testRes",i,".RData"))){
    load(paste0(filPat, "/rep", i, "LibsPars.RData"))
    load(file.path(filPat, "X.RData"))
    #Obtain bootstrap test statistics
    testBootRes = smoothTestMat(yBase, x = X)
    save(testBootRes, file = paste0(filPat, "/testResults/testRes",i,".RData"))
  }
})

#Now calculate the p-values (preparation)
odTrusted = 0.1
trustID = !tags %in% c("American gut project", "Armpit")
tagsTrusted = tags[trustID]
#rawNamesTrusted = rawNames[trustID]
rawNamesTrusted <- "Neuroblastoma (cell line)"
if(!file.exists(file ="testBootList.RData")){
  testBootList = lapply(rawNamesTrusted, loadTestStats)
  names(testBootList) = rawNamesTrusted
  #Now calculate the p-values
  pValsList = lapply(rawNamesTrusted, function(tag){
    realRes = testResultsList[[which(tag==rawNames)]]["testStat",]
    bootRes = testBootList[[tag]]
    names(bootRes) = names(paramList[[tags[rawNames==tag]]]$theta)[as.integer(names(bootRes))]
    pVals = sapply(names(bootRes), function(name){
      mean(realRes[name] < bootRes[[name]], na.rm = TRUE) #bootstrap p-value
    })
    return(pVals)
  })
  pValsList = lapply(pValsList, function(x){x[!is.na(names(x))]})
  names(pValsList) = tagsTrusted
  pValsListTrusted = lapply(names(pValsList), function(x){pValsList[[x]][paramList[[x]]$theta > odTrusted]})
  names(pValsListTrusted) = tagsTrusted
  save(pValsList, pValsListTrusted, file ="testBootList.RData")
} else {load(file ="testBootList.RData")}


# poor fit features

fdrtool(statistic = "pvalue", plot = FALSE, verbose = FALSE)


#3. GFRN dataset----------------------------------------------------------------\
rm(list=ls())
load("Datasets/GFRN_paperRaw.Rdata")
#3.1 Estimating dispersions
estNBparams = function(Y, X, prevCutOff = 1){
  Y = Y[rowSums(Y)>0, colMeans(Y==0)<prevCutOff] #去除表达值都为0或太低的基因，只保留有效基因
  Libs =  rowSums(Y) #计算每个样本的总表达量，即基因数
  tmpFit = apply(Y, 2, function(y){ # 对每个基因进行负二项分布的拟合，具体实现如下：
    nbFit = glm.nb2(y = y, reg = X, s = Libs) # 使用负二项分布拟合该基因的表达量数据，其中reg和s作为回归变量；
    list(theta = nbFit$theta, coef = nbFit$betas) # 返回拟合的参数，theta表示负二项分布的dispersion参数，coef表示回归系数；
  })
  list(Libs = Libs, theta = sapply(tmpFit, function(x){x$theta}), coef = sapply(tmpFit, function(x){x$coef}), X = X)
  # 返回拟合结果，包括每个样本的总表达量、每个基因的dispersion参数和回归系数，以及协变量矩阵X。
} 
X <- model.matrix(~group_sub)
Y <- cts_sub
# 此function只能计算condition 为两种的情况，目前参数估计结果为NA
raw_NB <- estNBparams(Y, X, prevCutOff = 1)
#----
# 改用edgeR包 glmfit 函数进行估计

# 数据收集
counts <- sim_counts_5
batch <- as.factor(batch)
n_batch <- nlevels(batch)  # number of batches
dge_obj <- DGEList(counts=counts)
batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch 
n_batches <- sapply(batches_ind, length)
n_sample <- sum(n_batches)
batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
mod <- model.matrix(~group)
design <- cbind(batchmod, mod)

## Check for intercept in covariates, and drop if present
check <- apply(design, 2, function(x) all(x == 1))
#if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
design <- as.matrix(design[,!check])
cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')

## Check if the design is confounded
if(qr(design)$rank<ncol(design)){
  #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
  if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
  if(ncol(design)>(n_batch+1)){
    if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
    }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
}

########  Estimate gene-wise dispersions within each batch  ########
cat("Estimating dispersions\n")
## Estimate common dispersion within each batch as an initial value
disp_common <- sapply(1:n_batch, function(i){
  if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
    # not enough residual degree of freedom
    return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
  }else{
    return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
  } 
})

## Estimate gene-wise dispersion within each batch 
genewise_disp_lst <- lapply(1:n_batch, function(j){
  if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
    # not enough residual degrees of freedom - use the common dispersion
    return(rep(disp_common[j], nrow(counts)))
  }else{
    return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                  dispersion=disp_common[j], prior.df=0))
  }
})
names(genewise_disp_lst) <- paste0('batch', levels(batch))

## construct dispersion matrix
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}

phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
for(k in 1:n_batch){
  phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
}

#######  Estimate parameters from NB GLM  ########
cat("Fitting the GLM model\n")
glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
  vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 

gamma_hat <- glm_f2$coefficients[, 1:n_batch]
mu_hat <- glm_f2$fitted.values
phi_hat <- do.call(cbind, genewise_disp_lst)

# ComBat_seq进行校正
sim_counts_5_combatseq <- ComBat_seq(sim_counts_5,batch,group) 

# pearson residuals assess glm model
residuals.glm(glm_f2,tpye="pearson") #  跑不通
mu <- glm_f2$fitted.values

#NBGOF test---------------------------------------------------------------------
# 缺陷：对于common tagwise dispersion model的估计不太准确

# extract quantities:
mu.hat.m = mu_hat   # mu may be close to 0
phi.hat.m = glm_f2$dispersion     # there may be NA's
v = mu.hat.m + phi.hat.m * mu.hat.m^2 # variance of counts
res.m = as.matrix((counts - mu.hat.m) / sqrt(v)) # residuals

# make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
res.m[ is.nan(res.m) ] = 0
res.m[ is.infinite(res.m) ] = 0

# sort res.m with care!
grp.ids = factor(apply(X, 1, function(x){paste(rev(x), collapse = ".")}),  # 为每个样本生成一个ID
                 labels = seq(ncol(X)))
sort.vec = function(x, grp.ids)  ave(x, grp.ids, FUN = sort)
res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
ord.res.v = as.vector(t(res.om))

model_tgc_m_obj = list(mu.hat.mat = mu.hat.m,
                       res.mat = res.m,
                       res.omat = res.om,
                       ord.res.vec = ord.res.v,
                       phi.hat.mat = phi.hat.m
)

mu.hat.mat0 = model_tgc_m_obj$mu.hat.mat
phi.hat.mat0 = model_tgc_m_obj$phi.hat.mat
res.omat0 = model_tgc_m_obj$res.omat
ord.res.vec0 = model_tgc_m_obj$ord.res.vec


#simulate new datasets and re-fit
sim <- 999
m <- dim(counts)[1]
n <- dim(counts)[2]
N <- m*n
seed <- 539
ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
  #setTxtProgressBar(pb, i/sim)
  set.seed(i+seed)
  y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
  dim(y.mat.h) = dim(counts)
  rownames(y.mat.h) = rownames(counts)
  colnames(y.mat.h) = colnames(counts)
  model.edgeR.tagcom(y.mat.h, X, lib.sizes=colSums(y.mat.h), prior.df = 43, method="CoxReid")$ord.res.vec
  #model.edgeR.tagcom(y.mat.h, X, lib.sizes=colSums(y.mat.h), prior.df = 43, method="CoxReid")
}

dimnames(ord.res.sim.mat.tmp) = NULL
ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)

ord.typ.res.sim = median(ord.res.sim.mat[sim,]) # on simulated datasets ONLY!
# subtract the typical residual vector from each row of the ordered residual matrix
dists.mat.res = (sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))^2
# construct new distance matrix D of dimension (R+1)-by-m
grp.vec = ( seq_len( ncol(dists.mat.res) ) - 1 ) %/% n     # grouping vector
dist.mat = t( rowsum(t(dists.mat.res), grp.vec) )    # vertical distance matrix (sim. + obs.)
# THIS dist.mat IS UN-SORTED
pear.mat = t( rowsum(t(ord.res.sim.mat)^2, grp.vec) )  # Pearson stats matrix (sim. + obs.)
## calcualte one p-value for each gene based on dist.mat or pear.mat, so a total of m p-values
## this p-value is the Monte Carlo p-value simply from a single univariate case
v.pvals = numeric(m)     # M.C. p-values based on vertical distances
p.pvals.1s = numeric(m)  # M.C. p-values based on Pearson statistics (1-sided)
for (i in 1:m){
  v.pvals[i] = (sum(dist.mat[1:sim,i] >= dist.mat[(sim+1),i]) + 1) / (sim + 1)
  p.pvals.1s[i] = (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1)
}

# to avoid potential zero p-values
v.pvals[v.pvals == 0] = 1/sim
p.pvals.1s[p.pvals.1s == 0] = 1/sim

# GOF_supp test-----------------------------------------------------------------
# referene: Ghent 
# GFRN dataset
batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
mod <- model.matrix(~group_sub)
design <- cbind(batchmod,mod)
check <- apply(design, 2, function(x) all(x == 1))
design <- as.matrix(design[,!check])

cts_glm <- edgeRdispersion(cts_sub,batch,group_sub)
cts_glm_GOF <- estNBparams(cts_sub,design)

# GOF test: AIC
#Quasi-Negative Binomial model
QNBD<- function(x,theta,beta,log = FALSE){
  
  loglik <- log((((beta-1)/(beta-1+beta*x))*(factorial(beta-1+beta*x))/
                   
                   (factorial(x)*factorial(beta-1+beta*x-x))*
                   
                   ((theta^x)*(1-theta)^(beta*x+beta-1-x))))
  
  if(log == FALSE)
    
    density <- exp(loglik)
  
  else density<-loglik
  
  return(density)
  
}
cts_para <- maxlogL(cts_sub,dist="QNBD")
summary(cts_para)

