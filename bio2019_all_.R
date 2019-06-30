## input:
##   - features: (samples x 32830/925032 gene/probeset) RNAseq counts (+extracted features); data has been batch and count normalized
##   - class: (367 train sample) gestational age 8-42 weeks
## output:
##   - class: (368 test sample) gestational age 8-42 weeks rounded to 1 decimal place
## created: 201906
## author: 


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
pkgs = c("Rfast", "stringr", "plyr", "dplyr", "Matrix", # var, str_, llply, etc
         "lattice", "ggplot2", # barplot, plots
         "foreach", "doMC", # parallel back-end
         "caret", "e1071", "ranger", "ANN2", "randomForest",
         "elasticnet", "fastICA", "foba", "glmnet","kernlab", 
         "KRLS", "lars", "leaps", "nnls", "nodeHarvest", 
         "partDSA", "pls", "plsRglm", "rpart", "rqPen",
         "RSNNS", "spikeslab", "xgboost") # ml
pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
sapply(pkgs, require, character.only=T)


## script options
no_cores = detectCores()-1 # number of cores to use in parallel
registerDoMC(no_cores)

overwrite = F # overwrite results?


## load input files

# sample annotation file columns:
#  SampleID: unique identifier of the sample (matching the name of the .CEL file in HTA20 folder, except for extension .CEL);
#  GA: gestational age as determined by the last menstrual period and or ultrasound;
#  Batch: the batch identifier;
#  Set: name of the source dataset;
#  Train: 1 for samples to be used for training, 0 for samples to be used for test;
meta = read.csv(paste0(input_dir,"/anoSC1_v11_nokey.csv"))

# submission template
class_final = read.csv(paste0(input_dir,"/TeamX_SC1_prediction.csv")) 

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
# 1: feature extraction ---------------------------------------
#--------#--------#--------#--------#--------#--------#--------

## temp data prep -------------------------------------------
# save only high variance genes and those with high sig pearson corr with GA
m = data0[,datavars>quantile(datavars,0.7) & abs(corsp)>quantile(abs(corsp),0.7) & corspp<quantile(corspp,0.1) & meancount>quantile(meancount,0.5)]
# # rfe to reduce features random forest
# rfe_res = rfe(m[tr_ind,], meta$GA[tr_ind], sizes=c(1:8), rfeControl=rfeControl(functions=rfFuncs, method="cv", number=10))
# print(rfe_res)
# predictors(rfe_res)
# plot(rfe_res, type=c("g", "o"))
write.csv(m, file=paste0(feat_dir,"/features_raw.csv"))


## 1) Removing genes with low variance ------------------------
eliminate = data.frame(col = c(1:dim(data0)[2]), var_gene = colVars(data0))
# head(eliminate)
eliminate = subset(eliminate, var_gene<quantile(eliminate$var_gene, 0.3))
# dim(data0)
data1 = subset(data0, select = -c(eliminate$col))
# dim(data1)


## 2) Removing elements with low correlation with GA ----------
data2 = data.frame(SampleID=row.names(data1), data1)
row.names(data2) = NULL
# Combining the two datasets 
data2 = merge(sample[,c(1,2,5)], data2, by.x = 'SampleID', by.y = 'SampleID', all = T)
#Using only the train dataset to make the feature importance 
data2 = subset(data2, Train == 1)
data2 = subset(data2, select = -c(Train,SampleID))

eliminate = data.frame(col = c(1:dim(data1)[2]),corr = apply(data2[,-1], 2, cor, x=data2$GA))
eliminate$corr = abs(eliminate$corr)
eliminate = subset(eliminate, corr<quantile(eliminate$corr,0.3))
dim(data2)
data2 = subset(data1, select = -c(eliminate$col))
dim(data2)
# write.csv(data2, paste0(feat_dir,'/features_raw.csv'), row.names=T)


## 3) PCA -----------------------------------------------------
# https://www.datacamp.com/community/tutorials/pca-analysis-r
# Don't need the values to be predicted, only the features
PrePCA = preProcess(data2,method="pca")
feat.pca = predict(PrePCA,data2)
# PrePCA
write.csv(feat.pca, paste0(feat_dir,'/features_pca.csv'), row.names=T)


## 4) Random Forest --------------------------------------------
#https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest
#https://uc-r.github.io/random_forests
#http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
# It depends on the features and the value to be predicted 
# data2 = as.matrix(data2)
metric = "Accuracy"
#Number randomely variable selected is mtry
PreRF = caret::train(y=meta$GA[tr_ind0], x=data2[tr_ind0,], method='ranger',importance='impurity')

# print(PreRF)
PreRF.i = varImp(PreRF)$importance
# PreRF.i = PreRF.i[order(PreRF.i$Overall, decreasing=T)[1:500],]
# PreRF.i = PreRF.i[1:500,]

feat.ra = data2[,order(PreRF.i, decreasing=T)[1:500]]
write.csv(feat.ra, paste0(feat_dir,'/features_ra.csv'), row.names = T)


## 5) Autoencoder --------------------------------------------
#https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/
#https://www.rdocumentation.org/packages/ANN2/versions/1.5/topics/autoencoder
#https://www.rdocumentation.org/packages/ANN2/versions/1.5/topics/autoencoder

#epochs changed to 100
preA = autoencoder(data2,hidden.layers = c(1000, 500, 1000))
feat.A = encode(preA, data2)
write.csv(feat.A, paste0(feat_dir,'/features_a.csv'), row.names = T)

# save.image()



#--------#--------#--------#--------#--------#--------#--------
# 2: regression models ----------------------------------------
#--------#--------#--------#--------#--------#--------#--------

## 0) load features ----------------------------------------
feat_paths = list.files(feat_dir) # feature paths
m0s = llply(feat_paths, function(xi) {
  m0 = read.csv(paste0(feat_dir,"/", xi))
  rownames(m0) = m0[,1]
  m0 = as.matrix(m0[,-1])
})
names(m0s) = gsub(".csv","",feat_paths)


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
# # cross validation function
# rmse = function(x,y) sqrt(mean((x-y)^2))
# cv_inds = split(tr_ind, cut(seq_along(tr_ind), cvn, labels=F))
# cv_class = llply(cv_inds, function(is) meta$GA[is])
# cv = function(data0, cv_inds, cv_class, fun, ...) {
#   llply (1:length(cv_inds), function(i) {
#     mtr = data0[unlist(cv_inds[-i]),]
#     ctr = unlist(cv_class[-i])
#     mte = data0[cv_inds[[i]],]
#     cte_ = unlist(cv_class[-i])
#     pte = fun(mtr, ctr, mte, ...)
#     return(list(pte=pte, rmse=rmse(pte,cte_)))
#   })
# }
# # usage:
# result_10x_rsme_pred = cv(data0, cv_inds, cv_class, function(mtr, ctr, mte) {
#   ...
#   return(pred) # vector of test class prediction
# })
fitcv = trainControl(method="cv", number=cvn)


## 2) test regression models ---------------------------------

# list models to test
models = c(# "ANFIS", # takes too long; RMSE 20
  #"avNNet","bag",
  # "bagEarth", # 8.9
  # "bagEarthGCV", # 8.6; repeat
  #"bam",        "bartMachine",
  # "bayesglm", # 10
  #"blackboost", 
  ## "blasso", # 8.4; takes a bit longer 950
  # "blassoAveraged", # 8.4; repeat
  ## "bridge", # 8.4; takes a bit longer 950
  #"brnn",       "BstLm",      "bstSm",      "bstTree",   
  #"cforest",    "ctree",      "ctree2",     "cubist",    
  #"DENFIS",     "dnn",        "earth",      "elm",       
  "enet", # 8.5      
  # "evtree",     "extraTrees", "FIR.DM", # no pkg  
  "foba", # 8.5      
  # "FS.HGD", # 9
  # "gam",        "gamboost", "gamLoess",   "gamSpline", # error on run
  # "gaussprLinear", # 10
  "gaussprPoly", # 8.3
  "gaussprRadial", # 8.8
  # "gbm", # 9
  # "gbm_h2o",  # error on run
  # "gcvEarth", # 9.5
  # "GFS.FR.MOGUL","GFS.LT.RS",  "GFS.THRIFT", # takes too long
  # "glm", "glm.nb", # 10
  # "glmboost",  # error on run
  "glmnet", # 8.5
  # "glmnet_h2o",  # error on run
  # "glmStepAIC", # 10
  # "HYFIS", # takes too long
  "icr", # 8.4
  "kernelpls", # 8.5
  # "kknn",# 9.6       
  # "knn", # 9       
  ## "krlsPoly", # 8.5; takes a bit longer 900
  "krlsRadial", # 8.5
  # "lars", # 8.5       
  "lars2", # 8.5      
  "lasso", # 8.7     
  "leapBackward", "leapForward", "leapSeq", # 8.5
  # "lm","lmStepAIC", # 10
  # "logicBag",   "logreg", "M5","M5Rules", # no pkg
  # "mlp", # 10
  # "mlpKerasDecay","mlpKerasDropout",  # error on run
  # "mlpSGD", # error on run
  # "mlpWeightDecay", # 10 
  # "mlpWeightDecayML", # 9 
  # "monmlp", # 10
  # "msaenet", # 10   
  # "mxnet",      "mxnetAdam", # no pkg
  # "neuralnet", # error on run
  # "nnet", # 26
  "nnls", # 8.5
  "null", # 8.5
  # "parRF", # 9      
  "partDSA", # 8.5
  # "pcaNNet", # 26   
  # "pcr", # 8.5
  # "penalized", # 9
  "pls",        "plsRglm",    
  # "ppr", "qrf", # 11     
  # "qrnn", "randomGLM", # takes too long
  # "ranger", # 9   
  "rbf", # 8.4       
  # "rbfDDA", # 27
  # "Rborist", # takes too long
  # "relaxo", # 8.7
  # "rf", #9.7
  # "rfRules", # error on run
  # "ridge", # 9.5
  # "rlm", # error on run       
  "rpart", # 8.4     
  # "rpart1SE", # 11  # "rpart2", # 8.6    
  "rqlasso", # 8.4
  # "RRF",        "RRFglobal", # 9 
  # "rvmLinear", # 8.6
  "rvmPoly",    "rvmRadial", # 8.5
  # "SBC", # 12        
  "simpls", "spikeslab",  # "spls", # 8.5; spls takes a bit longer 950
  # "superpc", # lowest score 27
  # "svmBoundrangeString",
  # "svmExpoString", # error on run
  # "svmLinear",  "svmLinear2", # 11 # "svmLinear3", # 9
  "svmPoly",    "svmRadial",  "svmRadialCost","svmRadialSigma", # 8.5
  # "svmSpectrumString", # error on run
  # "treebag", # 9.3   
  "widekernelpls", # 8.5
  # "WM", # 11
  "rqnc", # 8.4
  "nodeHarvest", # 8.6
  "mlpML", # 8.5
  "xgbDART" # extreme gradient boosting is good; "xgbLinear", # takes too long
  # "xgbTree",    "xyf"
)
# models = unique(modelLookup()[modelLookup()$forReg,c(1)])

# list model parameters to test
pars = list(
  # gbm=expand.grid(interaction.depth = c(1:5), #3 # gradient booted machine
  #                 n.trees = (seq(1,30,3))*50,
  #                 shrinkage = 0.1,
  #                 n.minobsinnode = 20),
  # ANFIS=expand.grid(max.iter=c(10,30,60,100)), # adaptive-network-based fuzzy inference system
  # avNNet=expand.grid(size=c(50,100), decay=10^runif(5, min = -5, 1), bag=T),
  #brnn=expand.grid(neurons=c(50,100)),
  enet=expand.grid(lambda=10^runif(5, min=-5, 1), fraction=runif(5, min=0, max=1)),
  # neuralnet=expand.grid(layer1=c(50,100),layer2=c(50,100),layer3=c(50,100)),
  # blackboost=expand.grid(maxdepth=c(1,3,6,10)),
  # xgbDART=expand.grid(nrounds, max_depth, eta, gamma, subsample, colsample_bytree, rate_drop, skip_drop, min_child_weight),
  # xgbLinear=expand.grid(nrounds, lambda, alpha, eta),
  # xgbTree=expand.grid(nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample),
  mlpML=expand.grid(layer1=c(50,100),layer2=c(50,100),layer3=c(50,100))
  # mlpKerasDropout=expand.grid(size=c(50,100), dropout=seq(0, .7, length=3), batch_size=floor(length(tr_ind)/3), lr=c(2e-6, 2e-3,.1,.5), rho=c(.2,.5,.9), decay=c(0,.3), activation=c("relu","softmax","linear")),
  # mlpKerasDecay=expand.grid(size=c(50,100), lambda=seq(0, .7, length=3), batch_size=floor(length(tr_ind)/3), lrlr=c(2e-6, 2e-3,.1,.5), rho=c(.2,.5,.9), decay=c(0,.3), activation=c("relu","softmax","linear"))
)

# use lapply/loop to run everything; best RMSE chosen by default
for (model in models) {
  cat("\n", model, " ------------------------------------------");
  for (xi in names(m0s)) {
    m0 = m0s[[xi]]
    mtr = m0[tr_ind,]
    
    dir.create(paste0(model_dir,"/",xi), showWarnings=F)
    fname = paste0(model_dir,"/",xi,"/",model,".Rdata")
    if (!file.exists(fname) | overwrite) { try ({ cat("\n", xi)
      t2i = NULL
      if (model%in%names(pars)) {
        t2i = caret::train(y=ctr, x=mtr, model, trControl=fitcv, tuneGrid=pars[[model]])
      } else {
        t2i = caret::train(y=ctr, x=mtr, model, trControl=fitcv)#, tuneGrid=pars[[model]])
      }
      if (!is.null(t2i)) save(t2i, file=fname)
    }) }
  }
}



## 3) load models -----------------------------------
feat_dirs = list.dirs(model_dir, full.names=F)
feat_dirs = feat_dirs[!feat_dirs%in%""]
result0 = llply(feat_dirs, function(data_type) { 
  models = gsub(".Rdata","",list.files(paste0(model_dir,"/",data_type)))
  a = llply(models, function(model) 
    get(load(paste0(model_dir,"/", data_type,"/",model,".Rdata"))) )
  names(a) = models
  return(a)
})
names(result0) = feat_dirs

# cat("min rmse's\n")
result = unlist(result0,recursive=F)
# r2 = llply(1:length(result), function(i) {
#   cat(models[i],": ",  round(result[[i]]$results$RMSE[which.min(result[[i]]$results$RMSE)],4),"\t",  result[[i]]$times$everything[3],"\n")
# })

# results to data frame
df1 = ldply (names(result), function(i) {
  fm = str_split(i,"[.]")[[1]]
  data.frame(
    rmse=result[[i]]$results$RMSE[which.min(result[[i]]$results$RMSE)],
    time=as.numeric(result[[i]]$times$everything[3]),
    model_=result[[i]]$modelInfo$label, 
    feature=fm[1], model=fm[2], 
    par=paste0( paste0(names(result[[i]]$bestTune), collapse="_"), ": ", paste0(result[[i]]$bestTune, collapse="_") )
    , stringsAsFactors=F)
})
# df1
write.csv(df1, file=paste0(result_dir,"/rmse_train.csv"))

## 4) print all results as prediction plots -----------------------
# result = unlist(result0, recursive=F)
# scores = ldply(names(result), function(xi) {
#   x = result[[xi]]
#   score = laply(x, function(xi) xi$rmse)
#   return(data.frame(feature.model=rep(xi,length(scores)),
#                     rmse=score, cv=1:length(scores)))
# })
# 
# png(paste0(result_dir, "/10cv.png"))
# pl = barchart(rmse~feature.model, data=scores, groups=cv, 
#               auto.key=list(columns=2),
#               cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
#               main="rmse scores over 10-fold cv's")
# print(pl)
# dev.off()


## 4) get test prediction results from models ----------------------
preds0 = llply(names(result0), function(xi) {
  m0 = m0s[[xi]]
  res = extractPrediction(result0[[xi]], 
    testX=m0[te_ind,], testY=cte, unkX=m0[-tr_ind0,]) # some features don't have test data
  return(res)
}, .parallel=T)
names(preds0) = names(result0)
save(preds0,file=paste0(result_dir,"/preds.Rdata"))

load(paste0(result_dir,"/preds.Rdata"))
wth = 200
dir.create(paste0(result_dir,"/obsVSpred"), showWarnings=F)
for (pred0n in names(preds0)) {
  # plot graph to compare models
  png(paste0(result_dir,"/obsVSpred/",pred0n,".png"), 
      width=wth*length(levels(preds0[[pred0n]]$model)))
  pl = plotObsVsPred(preds0[[pred0n]])
  print(pl)
  graphics.off()
  
  # png(paste0(result_dir,"/",xi,"_rmse.png"))
  # dotplot(caret::resamples(preds0[[xi]]))
  # graphics.off()
}


## 5) get rmse results and plot ------------------------------------
# preds = unlist(preds0,recursive=F)
rmse = function(x,y) sqrt(mean((x-y)^2))
rmsedf = ldply(names(preds0), function(xi) {
  x = preds0[[xi]]
  ldply(unique(x$model), function(mi) {
    mii = x$model==mi
    pr = x$pred[mii]
    data.frame(rmse=c(rmse(pr[1:length(tr_ind)], ctr),
                      rmse(pr[(length(tr_ind)+1):(length(tr_ind)+length(te_ind))], cte),
                      rmse(pr[1:length(tr_ind0)], meta$GA[append(tr_ind,te_ind)])),
               feature=rep(xi,3), model=rep(mi,3), type=c("train","test","all"))
  })
})
write.csv(rmsedf, file=paste0(result_dir,"/rmse.csv"))

wth = 2000
png(paste0(result_dir,"/rmse_all.png"), width=wth)
par(mfrow=c(3,1))
pl = barchart(rmse~model, data=rmsedf[rmsedf$type=="all",,drop=F], groups=feature, 
              auto.key = list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="all rmse for each model grouped by feature type")
print(pl)
graphics.off()
png(paste0(result_dir,"/rmse_test.png"), width=wth)
pl = barchart(rmse~model, data=rmsedf[rmsedf$type=="test",,drop=F], groups=feature, 
              auto.key = list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="test")
print(pl)
graphics.off()
png(paste0(result_dir,"/rmse_train.png"), width=wth)
pl = barchart(rmse~model, data=rmsedf[rmsedf$type=="train",,drop=F], groups=feature, 
              auto.key = list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="train")
print(pl)
graphics.off()


## 6) check for outlier subjects ------------------------------------
trte_diff = ldply(names(preds0), function(xi) {
  x = preds0[[xi]]
  xdf = ldply(unique(x$model), function(mi) {
    mii = x$model==mi
    pr = x$pred[mii]
    ob = x$pred[mii]
    mdf = data.frame(diff=c(pr[1:length(tr_ind0)] - meta$GA[append(tr_ind,te_ind)]), feat.model=paste0(rep(xi,length(tr_ind0)), ".", rep(mi,length(tr_ind0))), sample=meta$SampleID[append(tr_ind,te_ind)])
    return(mdf)
  })
  return(xdf)
})

png(paste0(result_dir,"/outliers.png"), width=5000)
pl = barchart(abs(diff)~sample, data=trte_diff[order(trte_diff$diff),], groups=feat.model, 
              auto.key = list(columns=5),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="all abs(pred-obs) for each sample grouped by feature.model type")
print(pl)
graphics.off()

# xdiffs = ldply(names(result), function(xi) {
#   xdiff = unlist(llply(result[[xi]], function(x) x$pte)) - meta$GA[tr_ind]
#   return(data.frame(feature.model=rep(xi,length(xdiff)),
#                     GAdiff=xdiff, sample=meta$SampleID[tr_ind]))
# })
# 
# png(paste0(result_dir, "/10cv_sample.png"))
# par(mfrow=c(1,2))
# pl = barchart(GAdiff~sample, data=xdiffs, groups=feature.model, 
#               auto.key=list(columns=2),
#               cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
#               main="prediction - ground truth, over sample")
# print(pl)
# pl = barchart(abs(GAdiff)~sample, data=xdiffs, groups=feature.model, 
#               auto.key=list(columns=2),
#               cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
#               main="absolute prediction - ground truth, over sample")
# print(pl)
# dev.off()


## 7) save final results of one model/feature ------------------
xi = "features_raw"
model = "enet"
# finalsol = llply(preds, function(x) x$pred[is.na(x$obs)])
class_final$GA = round(preds0[[xi]]$pred[is.na(preds0[[xi]]$obs) & preds0[[xi]]$model==model],1)
write.csv(class_final, file=paste0(result_dir,"/TeamX_SC1_prediction.csv"))
