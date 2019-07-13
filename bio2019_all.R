## input:
##   - features: (samples x 32830/925032 gene/probeset) RNAseq counts (+extraga_vald features); data has been batch and count normalized
##   - class: (367 train sample) gestational age 8-42 weeks
## output:
##   - class: (368 test sample) gestational age 8-42 weeks rounded to 1 decimal place
## created: 2019-06
## author: alice yue (aya43@sfu.ca); raquel aoki (raquel_aoki@sfu.ca)


#--------#--------#--------#--------#--------#--------#--------
# 0: preliminaries --------------------------------------------
#--------#--------#--------#--------#--------#--------#--------

rm(list=ls(all=T)) # clean the environment
set.seed(10)


## root: define the working directory of the project 
root = "/mnt/f/Brinkman group/current/Alice/ai4all2019_bio"
# root = "C:\\Users\\raque\\Documents\\GitHub\\ai4all2019_bio"
setwd(root)


## directories: create directories to save features, models, results 
input_dir = paste0(root,"/00_input") # raw data directory
feat_dir = paste0(root,"/01_features") # feature directory
model_dir = paste0(root, "/02_models") # model directory
result_dir = paste0(root, "/03_results") # stats/plots directory
sapply(c(input_dir,feat_dir, model_dir, result_dir), 
       function(x) dir.create(x, showWarnings=F))
cvinds_path = paste0(root,"/cvinds.Rdata") # train/test indices


## load packages
pkgs = c("Rfast", "stringr", "plyr", "dplyr", "Matrix", # var, str_, llply, etc
         "lattice", "ggplot2", # barplot, plots
         "foreach", "doMC", # parallel back-end
         "caret", "caretEnsemble", "Metrics",
         "deepnet", "e1071", "ranger", "ANN2", "randomForest",
         "elasticnet", "fastICA", "foba", "glmnet","kernlab", 
         "KRLS", "lars", "leaps", "nnls", "nodeHarvest", 
         "partDSA", "pls", "plsRglm", "rpart", "rqPen",
         "RSNNS", "spikeslab") # ml
pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
sapply(pkgs, require, character.only=T)


## script options
no_cores = detectCores()-1 # no of cores to use in parallel processing
registerDoMC(no_cores)

overwrite = F # overwrite results?


## load workspace for time consuming tasks 
load(paste(input_dir,'script.RData',sep='\\'))


## load input files

# sample annotation file
meta = read.csv(paste0(input_dir,"/anoSC1_v11_nokey.csv"))

# RNASEQ data
data0 = t(get(load(paste0(input_dir,"/HTA20_RMA.RData"))))

# submission template
submission = read.csv(paste0(input_dir,"/TeamX_SC1_prediction.csv")) 



## data exploration (ggplot in .rmd) ----------------------------
cat('Range\n');range(data0);
cat('Genes IDs\n'); gid = colnames(data0); head(gid);
cat('Patients IDs\n'); pid = rownames(data0); head(pid)
cat('Meta Data head and shape\n'); head(meta); dim(meta)
cat('RNASEQ head and shape\n'); head(data0[,c(1:5)]); dim(data0)
cat('Submission head and shape\n'); head(submission); dim(submission)

## plot stats: mean count, var, pearson/spearman corr (use only the train dataset to calculate)
meta_data_match = match(meta$SampleID, rownames(data0))
data1 = data.frame(GA=meta$GA[meta_data_match], data0)
data1 = data1[meta$Train[meta_data_match]==1,]

data_expl = data.frame(col=colnames(data1), 
                       corr_p=apply(data1,2,cor, x=data1$GA, method='pearson'),
                       corr_s=apply(data1,2,cor, x=data1$GA, method='spearman'),
                       variance=c(var(data1$GA),colVars(data0)),
                       mean=c(mean(data1$GA),colMeans(data0))
)
data_expl = data_expl[-1,] # removing GA from the dataset

# plot
p1 <- ggplot(data_expl, aes(x=mean)) + 
  geom_histogram(fill='lightgreen') + 
  xlab('Average Expression') + labs(title='(a)')

p2 <- ggplot(data_expl, aes(x=mean, y=variance)) + 
  geom_point(color='lightgreen')+
  xlab('Average Expression') + 
  ylab('Variance')+ labs(title='(b)')

p3 <- ggplot(data_expl, aes(x=corr_p)) + 
  geom_histogram(fill='lightgreen')+
  xlab('Pearson Correlation between Genes and Gestacional Age')+ labs(title='(c)')

p4 <- ggplot(data_expl, aes(x=corr_p, y=corr_s)) + 
  geom_point(color='lightgreen')+ labs(title='(d)')+
  xlab('Pearson Correlation') + ylab('Spearman Correlation')

grid.arrange(p1, p2, p3, p4, nrow=2)

# # archived version
# datavars = colVars(data0)
# meancount = colMeans(data0)
# meancounto = order(meancount)
# 
# train_index = which(meta$Train==1)
# data1 = data0[train_index,meancounto]
# ga_tr0 = meta$GA[train_index]
# 
# corpe = apply(data1, 2, function(x) cor(x, ga_tr0, method="pearson"))
# corpep = apply(data1, 2, function(x) cor.test(x, ga_tr0, method="pearson")$p.value)
# 
# corsp = apply(data1, 2, function(x) cor(x, ga_tr0, method="spearman"))
# corspp = apply(data1, 2, function(x) cor.test(x, ga_tr0, method="spearman")$p.value)
# 
# # plot stats
# png(paste0(result_dir,"/gene_stats.png"), width=800, height=1000)
# par(mfcol=c(4,1), mar=c(3,3,3,3))
# plot(density(data0), main="count distribution")
# plot(log(meancount), datavars, pch=16, cex=.3, 
#      main="gene ln(mean count) x variance")
# plot(corpe, cex=1-corpep, pch=16, 
#      main="gene (asc mean count order) x pearson corr with GA (size=1-pvalue)")
# # plot(abs(corpe), cex=1-corpep, pch=16)
# plot(corsp, cex=1-corspp, pch=16, 
#      main="gene (asc mean count order) x spearman corr with GA (size=1-pvalue)")
# # plot(abs(corsp), cex=1-corspp, pch=16)
# graphics.off()



#--------#--------#--------#--------#--------#--------#--------
# 1: feature extraction ---------------------------------------
#--------#--------#--------#--------#--------#--------#--------

## temp data prep -------------------------------------------
# save only high variance genes and those with high sig pearson corr with GA
# feat = data0[,datavars>quantile(datavars,0.7) & abs(corsp)>quantile(abs(corsp),0.7) & corspp<quantile(corspp,0.1) & meancount>quantile(meancount,0.5)]
# # rfe to reduce features random forest
# rfe_res = rfe(m[train_index_tr,], meta$GA[train_index_tr], sizes=c(1:8), rfeControl=rfeControl(functions=rfFuncs, method="cv", number=10))
# print(rfe_res)
# predictors(rfe_res)
# plot(rfe_res, type=c("g", "o"))
# write.csv(feat, file=paste0(feat_dir,"/features_raw.csv"))


## 1 and 2) remove genes with low variance and absolute correlation
data_expl$corr_p = abs(data_expl$corr_p)
keep = subset(data_expl, variance>quantile(data_expl$variance, 0.3) & corr_p>quantile(data_expl$corr_p, 0.3) & mean>quantile(mean,0.2))
data2 = subset(data0, select=keep$col)
data2 = scale(data2)

keep1 = subset(data_expl, variance>quantile(data_expl$variance, 0.8) & corr_p>quantile(data_expl$corr_p, 0.8) & mean>quantile(mean, 0.8))
data3 = subset(data0, select=keep1$col)
data3 = scale(data3)
write.csv(data3, paste0(feat_dir,'/features_raw.csv'), row.names=T) # save for future use


## 3) PCA -----------------------------------------------------
# https://www.datacamp.com/community/tutorials/pca-analysis-r
PrePCA = preProcess(data2, method="pca")
feat.pca = predict(PrePCA, data2)
write.csv(feat.pca, paste0(feat_dir,'/features_pca.csv'), row.names=T)


## 4) Random Forest --------------------------------------------
#https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest
#https://uc-r.github.io/random_forests
#http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
data1 = subset(data1, select =keep$col)
data1 = data.frame(ga_tr0, data1)
if (!exists('PreRF')) {
  PreRF = caret::train(y=data1[,1], x=data1[,-1], method='ranger', importance='impurity')
  PreRF.i = varImp(PreRF)$importance
} else {
  PreRF.i = varImp(PreRF)$importance
}

feat.ra = data2[,order(PreRF.i, decreasing=T)[1:500]]
write.csv(feat.ra, paste0(feat_dir,'/features_ra.csv'), row.names = T)


## 5) Autoencoder --------------------------------------------
#https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/
#https://www.rdocumentation.org/packages/ANN2/versions/1.5/topics/autoencoder
#https://www.rdocumentation.org/packages/ANN2/versions/1.5/topics/autoencoder
if(!exists('preA')){
  preA = autoencoder(data2, hidden.layers = c(1000, 500, 1000))
  feat.A = encode(preA, data2)  
}
rownames(feat.A) = rownames(feat.ra)
write.csv(feat.A, paste0(feat_dir,'/features_a.csv'), row.names = T)



#--------#--------#--------#--------#--------#--------#--------
# 2: regression models ----------------------------------------
#--------#--------#--------#--------#--------#--------#--------

## 0) load features ----------------------------------------
feat_paths = list.files(feat_dir) # feature paths
features = llply(feat_paths, function(xi) {
  feature = read.csv(paste0(feat_dir,"/", xi))
  rownames(feature) = feature[,1]
  feature = as.matrix(feature[,-1])
})
names(features) = gsub(".csv","",feat_paths)


## 1) prep cvn-fold cross validation & rmse function ----------
if (!file.exists(cvinds_path)) {
  cvn = 10
  train_index = which(meta$Train==1) #selecting only the training examples
  test_index = which(meta$Train==0)
  train_index_val = sample(train_index, ceiling(length(train_index)/11)) # cross validation set
  train_index_tr = sample(train_index[!train_index%in%train_index_val]) # remove cross validation set from training set 
  ga_val = as.numeric(meta$GA[train_index_val])
  ga_tr = as.numeric(meta$GA[train_index_tr]) 
  save(cvn, train_index, train_index_val, train_index_tr, ga_val, ga_tr, file=cvinds_path)
}
load(cvinds_path)

# 10-fold cv
fitcv = trainControl(method="cv", number=cvn)

# # cross validation function
# rmse = function(x,y) sqrt(mean((x-y)^2))
# cv_inds = split(train_index_tr, cut(seq_along(train_index_tr), cvn, labels=F))
# cv_class = llply(cv_inds, function(is) meta$GA[is])
# cv = function(data0, cv_inds, cv_class, fun, ...) {
#   llply (1:length(cv_inds), function(i) {
#     mtr = data0[unlist(cv_inds[-i]),]
#     ga_tr = unlist(cv_class[-i])
#     mte = data0[cv_inds[[i]],]
#     ga_val_ = unlist(cv_class[-i])
#     pte = fun(mtr, ga_tr, mte, ...)
#     return(list(pte=pte, rmse=rmse(pte,ga_val_)))
#   })
# }
# # # usage:
# # result_10x_rsme_pred = cv(data0, cv_inds, cv_class, function(mtr, ga_tr, mte) {
# #   ...
# #   return(pred) # vector of test class prediction
# # })


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
  #"cforest",    "ga_tree",      "ga_tree2",     "cubist",    
  #"DENFIS",     "earth",      "elm",       
  "dnn", 
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
  "lm", #"lmStepAIC", # 10
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
  "rvmPoly", # 8.5
  "rvmRadial", # 8.5
  # "SBC", # 12        
  "simpls", "spikeslab",  # "spls", # 8.5; spls takes a bit longer 950
  # "superpc", # lowest score 27
  # "svmBoundrangeString",
  # "svmExpoString", # error on run
  # "svmLinear",  "svmLinear2", # 11 # "svmLinear3", # 9
  "svmPoly",    "svmRadial",  "svmRadialCost","svmRadialSigma", # 8.5
  # "svmSpega_trumString", # error on run
  # "treebag", # 9.3   
  "widekernelpls", # 8.5
  # "WM", # 11
  # "rqnc", # 8.4
  "nodeHarvest" # 8.6
  # "mlpML" # 8.5
  # "xgbDART" # extreme gradient boosting is good; "xgbLinear", # takes too long
  # "xgbTree",    "xyf"
)
# models = unique(modelLookup()[modelLookup()$forReg,c(1)])

# list model parameters to test
pars = list(
  dnn=expand.grid(layer1=100, layer2=50, layer3=10, hidden_dropout=c(0,2,5), visible_dropout=c(0,2,5)), # drop out fraction for hidden/input layer. 
  # gbm=expand.grid(interaction.depth = c(1:5), #3 # gradient booted machine
  #                 n.trees = (seq(1,30,3))*50,
  #                 shrinkage = 0.1,
  #                 n.minobsinnode = 20),
  # ANFIS=expand.grid(max.iter=c(10,30,60,100)), # adaptive-network-based fuzzy inference system
  # avNNet=expand.grid(size=c(50,100), decay=10^runif(5, min = -5, 1), bag=T),
  #brnn=expand.grid(neurons=c(50,100)),
  enet=expand.grid(lambda=10^runif(3, min=-5, 1), fraction=runif(3, min=0, max=1)),
  # neuralnet=expand.grid(layer1=c(50,100),layer2=c(50,100),layer3=c(50,100)),
  # blackboost=expand.grid(maxdepth=c(1,3,6,10)),
  # xgbDART=expand.grid(nrounds, max_depth, eta, gamma, subsample, colsample_bytree, rate_drop, skip_drop, min_child_weight),
  # xgbLinear=expand.grid(nrounds, lambda, alpha, eta),
  # xgbTree=expand.grid(nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample),
  mlpML=expand.grid(layer1=c(50,100),layer2=c(50,100),layer3=c(50,100))
  # mlpKerasDropout=expand.grid(size=c(50,100), dropout=seq(0, .7, length=3), batch_size=floor(length(train_index_tr)/3), lr=c(2e-6, 2e-3,.1,.5), rho=c(.2,.5,.9), decay=c(0,.3), activation=c("relu","softmax","linear")),
  # mlpKerasDecay=expand.grid(size=c(50,100), lambda=seq(0, .7, length=3), batch_size=floor(length(train_index_tr)/3), lrlr=c(2e-6, 2e-3,.1,.5), rho=c(.2,.5,.9), decay=c(0,.3), activation=c("relu","softmax","linear"))
)

# use lapply/loop to run everything; best RMSE chosen by default
for (model in models) {
  cat("\n", model, " ------------------------------------------");
  for (xi in names(features)) {
    feature = features[[xi]]
    mtr = feature[train_index_tr,]
    
    dir.create(paste0(model_dir,"/",xi), showWarnings=F)
    fname = paste0(model_dir,"/",xi,"/",model,".Rdata")
    if (!file.exists(fname) | overwrite) { try ({ cat("\n", xi)
      t2i = NULL
      if (model%in%names(pars)) {
        t2i = caret::train(y=ga_tr, x=mtr, model, trControl=fitcv, tuneGrid=pars[[model]])
      } else {
        t2i = caret::train(y=ga_tr, x=mtr, model, trControl=fitcv)#, tuneGrid=pars[[model]])
      }
      if (!is.null(t2i)) save(t2i, file=fname)
    }) }
  }
}



## 3) load models -----------------------------------
fd = list.dirs(model_dir, full.names=F); feat_dirs = fd[!fd%in%""]
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
write.table(df1, file=paste0(result_dir,"/rmse_train.csv"), sep=",")

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
  feature = features[[xi]]
  res = extractPrediction(result0[[xi]], 
                          testX=feature[train_index_val,], testY=ga_val, unkX=feature[-train_index,]) # some features don't have test data
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
    pr1 = pr[1:length(train_index_tr)]
    pr1[pr1<8] = 8; pr1[pr1>42] = 42
    pr2 = pr[(length(train_index_tr)+1):(length(train_index_tr)+length(train_index_val))]
    pr2[pr2<8] = 8; pr2[pr2>42] = 42
    pr3 = pr[1:length(train_index)]
    pr3[pr3<8] = 8; pr3[pr3>42] = 42
    data.frame(rmse=c(rmse(pr1, ga_tr), rmse(pr2, ga_val),
                      rmse(pr3, meta$GA[append(train_index_tr,train_index_val)])),
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
    mdf = data.frame(diff=c(pr[1:length(train_index)] - meta$GA[append(train_index_tr,train_index_val)]), feat.model=paste0(rep(xi,length(train_index)), ".", rep(mi,length(train_index))), sample=meta$SampleID[append(train_index_tr,train_index_val)])
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
#   xdiff = unlist(llply(result[[xi]], function(x) x$pte)) - meta$GA[train_index_tr]
#   return(data.frame(feature.model=rep(xi,length(xdiff)),
#                     GAdiff=xdiff, sample=meta$SampleID[train_index_tr]))
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
submission$GA = round(preds0[[xi]]$pred[is.na(preds0[[xi]]$obs) & preds0[[xi]]$model==model],1)
submission$GA[submission$GA<8] = 8
submission$GA[submission$GA>42] = 42
write.csv(submission, file=paste0(result_dir,"/TeamX_SC1_prediction.csv"), row.names=F)



#--------#--------#--------#--------#--------#--------#--------
# 3: extra - ensembles ----------------------------------------
#--------#--------#--------#--------#--------#--------#--------
load(cvinds_path)
fitcv = trainControl(method="cv", number=10)
rmse = function(x,y) sqrt(mean((x-y)^2))

# load features
xi = "features_pca"
feature = read.csv(paste0(feat_dir,"/", xi,".csv"))
rownames(feature) = feature[,1]; feature = as.matrix(feature[,-1])
mtr = feature[train_index_tr,]
mte = feature[train_index_val,]
mte0 = feature[test_index,]

# EDIT HERE - train models; note some are classification models too
# exercise: pick and choose a few models to ensemble
models_reg = c(
  "enet", # 8.5
  "glmnet", # 8.5
  "kernelpls", # 8.5
  "krlsRadial", # 8.5
  "lars2", # 8.5
  "nnls", # 8.5
  "rbf", # 8.4
  "rqlasso", # 8.4
  "rvmPoly", # 8.5
  "rvmRadial", # 8.5
  "simpls", "spikeslab", # "spls", # 8.5; spls takes a bit longer 950
  "rqnc" # 8.4
)
model_list = caretList(
  y=ga_tr, x=mtr, trControl=fitcv, metric="RMSE",
  methodList=models_reg, continue_on_fail=T)

# weighted linear combination of model predictions
ensemble = caretEnsemble(
  model_list, 
  metric="RMSE",
  trControl=fitcv)
summary(ensemble)

# predict test set GA; get rmse
# exercise: do ensembling methods perform better? why do you think that is?
ens_preds = predict(ensemble, newdata=mte)
rmse(ens_preds, ga_val)

# save as submission
ens_predfinal = predict(ensemble, newdata=mte0)
submission$GA = round(ens_predfinal,1)
submission$GA[submission$GA<8] = 8
submission$GA[submission$GA>42] = 42
write.csv(submission, file=paste0(result_dir,"/submission_ensemble.csv"), row.names=F)
