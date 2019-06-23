## created: 201906


#--------#--------#--------#--------#--------#--------#--------
# 0: preliminaries --------------------------------------------
#--------#--------#--------#--------#--------#--------#--------

rm(list=ls(all=T)) # clean the environment
set.seed(10)


## root 
root = "/mnt/f/Brinkman group/current/Alice/ai4all2019_bio"
setwd(root)


## directories
input_dir = paste0(root,"/00_input") # raw data directory
feat_dir = paste0(root,"/01_features") # feature directory
model_dir = paste0(root, "/02_models") # model directory
result_dir = paste0(root, "/03_results") # stats/plots directory
sapply(c(feat_dir, model_dir, result_dir), 
       function(x) dir.create(x, showWarnings=F))


## load packages; need to fix according to what model we'll be using
pkgs = c("Rfast", "stringr", "plyr", "dplyr", # var, str_, llply, etc
         "lattice", "ggplot2", # barplot, plots
         "foreach", "doMC", # parallel back-end
         "caret", "e1071", "ranger", "ANN2", "randomForest",
         "glmnet","RSNNS","plsRglm") # ml
# pkgs = c("frbs", "brnn", "monomvn", "Cubist", "elasticnet", "fastICA", "lars", "leaps", "MASS", "RWeka", "neuralnet", "rqPen", "nnls", "penalized", "KRLS", "pls", "quantregForest", "qrnn", "rqPen", "kernlab", "relaxo", "foba", "spikeslab", "superpc", "ipred", "e1071", "logicFS", "earth", "bartMachine", "arm", "mboost", "import", "bst", "party", "partykit", "rpart", "randomGLM", "xgboost", "elmNN", "gam", "mgcv", "h2o", "kknn", "LiblineaR", "LogicReg", "nnet", "monmlp", "RSNNS", "msaenet", "FCNN4R", "keras", "mxnet", "partDSA", "plsRglm", "ranger", "Rborist", "randomForest", "extraTrees", "RRF", "kohonen", "spls", "deepnet", "gbm", "evtree", "nodeHarvest")
pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
sapply(pkgs, require, character.only=T)


## script options
no_cores = detectCores()-1 # number of cores to use in parallel
registerDoMC(no_cores)


## load input files

# sample annotation file columns:
#  SampleID: unique identifier of the sample (matching the name of the .CEL file in HTA20 folder, except for extension .CEL);
#  GA: gestational age as determined by the last menstrual period and or ultrasound;
#  Batch: the batch identifier;
#  Set: name of the source dataset;
#  Train: 1 for samples to be used for training, 0 for samples to be used for test;
meta = read.csv(paste0(input_dir,"/anoSC1_v11_nokey.csv"))

# RNASEQ data: each row is a gene and each column a patient 
#  probeset: gene ID
#  pid: patient id 
data0 = t(get(load(paste0(input_dir,"/HTA20_RMA.RData"))))


## data exploration ----------------------------
range(data0) # range: 0.9365626 14.2836285
gid = colnames(data0); head(gid) # genes IDS
pid = rownames(data0); head(pid) # Patient IDS  

# plot stats: mean count, pearson/spearman corr
datavars = colVars(data0)
meancount = colMeans(data0)
meancounto = order(meancount)

tr_ind0 = which(meta$Train==1)
corpe = apply(data0[tr_ind0,meancounto], 2, function(x) 
  cor(x, meta$GA[tr_ind0], method="pearson"))
corpep = apply(data0[tr_ind0,meancounto], 2, function(x) 
  cor.test(x, meta$GA[tr_ind0], method="pearson")$p.value)

corsp = apply(data0[tr_ind0,meancounto], 2, function(x) 
  cor(x, meta$GA[tr_ind0], method="spearman"))
corspp = apply(data0[tr_ind0,meancounto], 2, function(x) 
  cor.test(x, meta$GA[tr_ind0], method="spearman")$p.value)

# plot stats
png(paste0(result_dir,"/gene_stats.png"), width=800, height=1000)
par(mfcol=c(4,1), mar=c(3,3,3,3))
plot(density(data0), main="count distribution")
plot(log(meancount), datavars, pch=16, cex=.3, 
     main="gene ln(mean count) x variance")
plot(corpe, cex=1-corpep, pch=16, 
     main="gene (asc mean count order) x pearson corr with GA (size=1-pvalue)")
# plot(abs(corpe), cex=1-corpep, pch=16)
plot(corsp, cex=1-corspp, pch=16, 
     main="gene (asc mean count order) x spearman corr with GA (size=1-pvalue)")
# plot(abs(corsp), cex=1-corspp, pch=16)
graphics.off()



#--------#--------#--------#--------#--------#--------#--------
# 1: feature extraction & model testing -----------------------
#--------#--------#--------#--------#--------#--------#--------

# overwrite model?
overwrite = F

# feature
xi = "features_raw" # ra, a, pca

# model
model = "enet"

# list model parameters to test; see https://topepo.github.io/caret/available-models.html
parsi = expand.grid(lambda=10^runif(5, min=-5, 1), fraction=runif(5, min=0, max=1)) # parsi = NULL; if no parameters need to be tested

## 0) load feature ----------------------------------------
m0 = read.csv(paste0(feat_dir,"/", xi))
rownames(m0) = m0[,1]
m0 = as.matrix(m0[,-1])


## 1) prep cvn-fold cross validation & rmse function ----------
cvinds_path = paste0(root,"/cvinds.Rdata")
if (file.exists(cvinds_path)) {
  load(cvinds_path)
} else {
  cvn = 10
  tr_ind0 = which(meta$Train==1)
  te_ind = sample(tr_ind0, ceiling(length(tr_ind0)/11))
  tr_ind = sample(tr_ind0[!tr_ind0%in%te_ind])
  ctr = as.numeric(meta$GA[tr_ind])
  cte = as.numeric(meta$GA[te_ind])
  save(cvn,tr_ind0,te_ind,tr_ind,ctr,cte, file=cvinds_path)
}
fitcv = trainControl(method="cv", number=cvn)


## 2) test regression models ---------------------------------

# best RMSE parameters chosen by default
mtr = m0[tr_ind,]
dir.create(paste0(model_dir,"/",xi), showWarnings=F)
fname = paste0(model_dir,"/",xi,"/",model,".Rdata")
if (!file.exists(fname) | overwrite) { try ({
  t2i = NULL
  if (!is.null(parsi)) {
    t2i = caret::train(y=ctr, x=mtr, model, trControl=fitcv, tuneGrid=parsi)
  } else {
    t2i = caret::train(y=ctr, x=mtr, model, trControl=fitcv)
  }
  if (!is.null(t2i)) save(t2i, file=fname)
}) }
load(fname)

# results to data frame
df = data.frame(
  rmse=t2i$results$RMSE[which.min(t2i$results$RMSE)],
  time=as.numeric(t2i$times$everything[3]),
  model_=t2i$modelInfo$label, 
  feature=xi, model=model, 
  par=paste0( paste0(names(t2i$bestTune), collapse="_"), ": ", paste0(t2i$bestTune, collapse="_") )
  , stringsAsFactors=F)
df



## 4) get test prediction results from models ----------------------
pred = extractPrediction(t2i, testX=m0[te_ind,], testY=cte, unkX=m0[-tr_ind0,])

# plot graph to compare models
wth = 200
dir.create(paste0(result_dir,"/obsVSpred"), showWarnings=F)
png(paste0(result_dir,"/obsVSpred/",xi,"_",model,".png"), width=wth)
pl = plotObsVsPred(pred)
print(pl)
graphics.off()


## 5) get rmse results and plot ------------------------------------
# preds = unlist(preds0,recursive=F)
rmse = function(x,y) sqrt(mean((x-y)^2))
rdf = ldply(unique(pred$model), function(mi) {
  mii = pred$model==mi
  pr = pred$pred[mii]
  data.frame(rmse=c(rmse(pr[1:length(tr_ind)], ctr),
                    rmse(pr[(length(tr_ind)+1):(length(tr_ind)+length(te_ind))], cte),
                    rmse(pr[1:length(tr_ind0)], meta$GA[append(tr_ind,te_ind)])),
             feature=rep(xi,3), model=rep(mi,3), type=c("train","test","all"))
})
rdf



## 7) save final results of one model/feature ------------------

# submission template
class_final = read.csv(paste0(input_dir,"/TeamX_SC1_prediction.csv")) 

class_final$GA = round(pred$pred[is.na(pred$obs) & pred$model==model],1)
write.csv(class_final, file=paste0(result_dir,"/TeamX_SC1_prediction.csv"))
