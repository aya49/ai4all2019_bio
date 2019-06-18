## input:
##   - features: (samples x 32830/925032 gene/probeset) RNAseq counts (+extracted features); data has been batch and count normalized
##   - class: (367 train sample) gestational age 8-42 weeks
## output:
##   - class: (368 test sample) gestational age 8-42 weeks rounded to 1 decimal place
## process:
## created: 20190605
## author: alice yue

set.seed(10)

## root
root = "/mnt/f/Brinkman group/current/Alice/ai4all2019_bio"
setwd(root)

## input
input_dir = paste0(root,"/input")
meta_dir = paste0(input_dir,"/anoSC1_v11_nokey.csv") # meta data
class_dir = paste0(input_dir,"/TeamX_SC1_prediction.csv") # class
data_dir = paste0(root,"/data") # features directory

## ouput
result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)
models_dir = paste0(result_dir, "/models"); dir.create(models_dir, showWarnings=F)
preds_dir = paste0(result_dir, "/preds.Rdata")

overwrite = F

## load packages
libr = function(pkgs) {
  pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
  if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
  sapply(pkgs, require, character.only=T)
}
libr(c("stringr", "plyr", "lattice", "foreach", "doMC", "Rfast", "caret"))
no_cores = detectCores()-6
registerDoMC(no_cores)


## load files
meta = read.csv(meta_dir)
class_final = read.csv(class_dir)
data_paths = list.files(data_dir)
m0 = t(get(load(paste0(input_dir,"/HTA20_RMA.RData"))))
range(m0) # range: 0.9365626 14.2836285
m0vars = colVars(m0) 


## prep 10-fold cross validation & rmse function
cvn = 10
rmse = function(x,y) sqrt(mean((x-y)^2))
tr_ind0 = which(meta$Train==1)
te_ind = sample(tr_ind0, ceiling(length(tr_ind0)/11))
tr_ind = sample(tr_ind0[!tr_ind0%in%te_ind])
# cv_inds = split(tr_ind, cut(seq_along(tr_ind), cvn, labels=F))
# cv_class = llply(cv_inds, function(is) meta$GA[is])
# cv = function(m0, cv_inds, cv_class, fun, ...) {
#   llply (1:length(cv_inds), function(i) {
#     mtr = m0[unlist(cv_inds[-i]),]
#     ctr = unlist(cv_class[-i])
#     mte = m0[cv_inds[[i]],]
#     cte_ = unlist(cv_class[-i])
#     pte = fun(mtr, ctr, mte, ...)
#     return(list(pte=pte, rmse=rmse(pte,cte_)))
#   })
# }
# # usage:
# result_10x_rsme_pred = cv(m0, cv_inds, cv_class, function(mtr, ctr, mte) {
#   ...
#   return(pred) # vector of test class prediction
# })
fitcv = trainControl(method="cv", number=cvn)


## plot stats: mean count, pearson/spearman corr
meancount = colMeans(m0)
meancounto = order(meancount)

corpe = apply(m0[tr_ind0,meancounto], 2, function(x) 
  cor(x, meta$GA[tr_ind0], method="pearson"))
corpep = apply(m0[tr_ind0,meancounto], 2, function(x) 
  cor.test(x, meta$GA[tr_ind0], method="pearson")$p.value)

corsp = apply(m0[tr_ind0,meancounto], 2, function(x) 
  cor(x, meta$GA[tr_ind0], method="spearman"))
corspp = apply(m0[tr_ind0,meancounto], 2, function(x) 
  cor.test(x, meta$GA[tr_ind0], method="spearman")$p.value)

# plot stats
png(paste0(result_dir,"/gene_stats.png"), width=800, height=1000)
par(mfcol=c(4,1), mar=c(3,3,3,3))
plot(density(m0), main="count distribution")
plot(log(meancount), m0vars, pch=16, cex=.3, 
     main="gene ln(mean count) x variance")
plot(corpe, cex=1-corpep, pch=16, 
     main="gene (asc mean count order) x pearson corr with GA (size=1-pvalue)")
# plot(abs(corpe), cex=1-corpep, pch=16)
plot(corsp, cex=1-corspp, pch=16, 
     main="gene (asc mean count order) x spearman corr with GA (size=1-pvalue)")
# plot(abs(corsp), cex=1-corspp, pch=16)
graphics.off()


## temp data prep (rm bottom 10% var genes)
# save only high variance genes
m = m0[,m0vars>quantile(m0vars,0.7) & abs(corsp)>quantile(abs(corsp),0.9) & corsp<quantile(corsp,0.9) & meancount>quantile(meancount,0.3)] # 29547 features, too much?
# # rfe to reduce features random forest (for testing)
# rfe_res = rfe(m[tr_ind,], meta$GA[tr_ind], sizes=c(1:8), rfeControl=rfeControl(functions=rfFuncs, method="cv", number=10))
# print(rfe_res)
# predictors(rfe_res)
# plot(rfe_res, type=c("g", "o"))
save(m, file=paste0(data_dir,"/raw.Rdata"))


## test regression models ---------------------------------
models = c("ANFIS","avNNet","bag","bagEarth",
           "bagEarthGCV","bam",        "bartMachine","bayesglm",  
           "blackboost", "blasso",     "blassoAveraged","bridge",    
           "brnn",       "BstLm",      "bstSm",      "bstTree",   
           "cforest",    "ctree",      "ctree2",     "cubist",    
           "DENFIS",     "dnn",        "earth",      "elm",       
           "enet",       "evtree",     "extraTrees", "FIR.DM",    
           "foba",       "FS.HGD",     "gam",        "gamboost",  
           "gamLoess",   "gamSpline",  "gaussprLinear","gaussprPoly",        
           "gaussprRadial","gbm",        "gbm_h2o",    "gcvEarth",  
           "GFS.FR.MOGUL","GFS.LT.RS",  "GFS.THRIFT", "glm",       
           "glmboost",   "glm.nb",     "glmnet",     "glmnet_h2o",
           "glmStepAIC", "HYFIS",      "icr",        "kernelpls", 
           "kknn",       "knn",        "krlsPoly",   "krlsRadial",
           "lars",       "lars2",      "lasso",      "leapBackward",
           "leapForward","leapSeq",    "lm","lmStepAIC", 
           "logicBag",   "logreg",     "M5","M5Rules",   
           "mlp",        "mlpKerasDecay","mlpKerasDropout","mlpML",     
           "mlpSGD",     "mlpWeightDecay","mlpWeightDecayML","monmlp",    
           "msaenet",    "mxnet",      "mxnetAdam",  "neuralnet", 
           "nnet",       "nnls",       "nodeHarvest","null",      
           "parRF",      "partDSA",    "pcaNNet",    "pcr",       
           "penalized",  "pls",        "plsRglm",    "ppr",       
           "qrf",        "qrnn",       "randomGLM",  "ranger",    
           "rbf",        "rbfDDA",     "Rborist",    "relaxo",    
           "rf","rfRules",    "ridge",      "rlm",       
           "rpart",      "rpart1SE",   "rpart2",     "rqlasso",   
           "rqnc",       "RRF",        "RRFglobal",  "rvmLinear", 
           "rvmPoly",    "rvmRadial",  "SBC",        "simpls",    
           "spikeslab",  "spls",       "superpc",    "svmBoundrangeString",
           "svmExpoString","svmLinear",  "svmLinear2", "svmLinear3",
           "svmPoly",    "svmRadial",  "svmRadialCost","svmRadialSigma",     
           "svmSpectrumString","treebag",    "widekernelpls","WM",        
           "xgbDART",    "xgbLinear",  "xgbTree",    "xyf")
# models <- unique(modelLookup()[modelLookup()$forReg,c(1)])
models = models[!models%in%c("M5Rules", "M5", "logicBag", "bartMachine", "elm", "mlpSGD", "mxnet", "mxnetAdam", "extraTrees")]

pkgs = c("frbs", "brnn", "monomvn", "Cubist", "elasticnet", "fastICA", "lars", "leaps", "MASS", "RWeka", "neuralnet", "rqPen", "nnls", "penalized", "KRLS", "pls", "quantregForest", "qrnn", "rqPen", "kernlab", "relaxo", "foba", "spikeslab", "superpc", "ipred", "e1071", "logicFS", "earth", "bartMachine", "arm", "mboost", "import", "bst", "party", "partykit", "rpart", "randomGLM", "xgboost", "elmNN", "gam", "mgcv", "h2o", "kknn", "LiblineaR", "LogicReg", "nnet", "monmlp", "RSNNS", "msaenet", "FCNN4R", "keras", "mxnet", "partDSA", "plsRglm", "ranger", "Rborist", "randomForest", "extraTrees", "RRF", "kohonen", "spls", "deepnet", "gbm", "evtree", "nodeHarvest")
pkgs_ui = setdiff(pkgs, rownames(installed.packages())) # c("RWeka", "logicFS", "bartMachine", "elmNN", "FCNN4R", "mxnet", "extraTrees")
pkgs = pkgs[!pkgs%in%pkgs_ui]
libr(pkgs)



# can't be access outside ;w;
pars = list(
  gbm=expand.grid(interaction.depth = c(1:5), #3 # gradient booted machine
                  n.trees = (seq(1,30,3))*50,
                  shrinkage = 0.1,
                  n.minobsinnode = 20),
  # ANFIS=expand.grid(max.iter=c(10,30,60,100)), # adaptive-network-based fuzzy inference system
  avNNet=expand.grid(size=c(50,100), decay=10^runif(5, min = -5, 1), bag=T),
  brnn=expand.grid(neurons=c(50,100)),
  enet=expand.grid(lambda=10^runif(5, min=-5, 1), fraction=runif(5, min=0, max=1)),
  neuralnet=expand.grid(layer1=c(50,100),layer2=c(50,100),layer3=c(50,100)),
  blackboost=expand.grid(maxdepth=c(1,3,6,10)),
  # xgbDART=expand.grid(nrounds, max_depth, eta, gamma, subsample, colsample_bytree, rate_drop, skip_drop, min_child_weight),
  # xgbLinear=expand.grid(nrounds, lambda, alpha, eta),
  # xgbTree=expand.grid(nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample),
  mlpML=expand.grid(layer1=c(50,100),layer2=c(50,100),layer3=c(50,100)),
  mlpKerasDropout=expand.grid(size=c(50,100), dropout=seq(0, .7, length=3), batch_size=floor(length(tr_ind)/3), lr=c(2e-6, 2e-3,.1,.5), rho=c(.2,.5,.9), decay=c(0,.3), activation=c("relu","softmax","linear")),
  mlpKerasDecay=expand.grid(size=c(50,100), lambda=seq(0, .7, length=3), batch_size=floor(length(tr_ind)/3), lrlr=c(2e-6, 2e-3,.1,.5), rho=c(.2,.5,.9), decay=c(0,.3), activation=c("relu","softmax","linear"))
)

# models = c("glm","enet") # test

ctr = meta$GA[tr_ind]
# result0 = NULL
for (data_path in gsub(".Rdata","",data_paths)) {
  cat(data_path, " ------------------------------------------\n");
  
  m0 = as.matrix(get(load(paste0(data_dir,"/", data_path,".Rdata"))))
  class(m0)="numeric"
  mtr = m0[tr_ind,]

  # use lapply/loop to run everything; best RMSE chosen by default
  data_path_ = paste0(models_dir,"/",data_path)
  dir.create(data_path_, showWarnings=F)
  for (model in models) {
    fname = paste0(data_path_,"/",model,".Rdata")
    if (!file.exists(fname) | overwrite) {
      try ({
        t2i = NULL
        if (model%in%names(pars)) {
          t2i = caret::train(y=ctr, x=mtr, (model), trControl=fitcv, tuneGrid=pars[[model]])
        } else {
          t2i = caret::train(y=ctr, x=mtr, (model), trControl=fitcv)#, tuneGrid=pars[[model]])
        }
        if (!is.null(t2i)) save(t2i, file=fname)
      })
    }
  }
}





## print training results as table ---------------------------------
data_dirs = list.dirs(models_dir, full.names=F)
data_dirs = data_dirs[!data_dirs%in%""]
result0 = llply(data_dirs, function(data_type) { 
  models = gsub(".Rdata","",list.files(paste0(models_dir,"/",data_type)))
  a = llply(models, function(model) 
    get(load(paste0(models_dir,"/", data_type,"/",model,".Rdata"))) )
  names(a) = models
  return(a)
})
names(result0) = data_dirs

# cat("min rmse's\n")
result = unlist(result0,recursive=F)
# r2 = llply(1:length(result), function(i) {
#   cat(models[i],": ",  round(result[[i]]$results$RMSE[which.min(result[[i]]$results$RMSE)],4),"\t",  result[[i]]$times$everything[3],"\n")
# })

# resuls to data frame
df1 = ldply (names(result), function(i) {
  fm = str_split(i,"[.]")[[1]]
  score = data.frame(
    rmse=result[[i]]$results$RMSE[which.min(result[[i]]$results$RMSE)],
    time=as.numeric(result[[i]]$times$everything[3]),
    model_=result[[i]]$modelInfo$label, 
    feature=fm[1], model=fm[2], 
    par=paste0( paste0(names(result[[i]]$bestTune), collapse="_"), ": ", paste0(result[[i]]$bestTune, collapse="_") )
, stringsAsFactors=F)
})
# df1


## print all results as plots ------------------------------------
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

cte = as.numeric(meta$GA[te_ind])
preds0 = NULL
for (data_path in names(result0)) {
  m0 = as.matrix(get(load(paste0(data_dir,"/", data_path,".Rdata"))))
  class(m0) = "numeric"
  preds0[[data_path]] = extractPrediction(
    result0[[data_path]], 
    testX=m0[te_ind,], testY=cte, unkX=m0[-tr_ind0,])
}
save(preds0,file=preds_dir)

load(preds_dir)
for (pred0n in names(preds0)) {
  # plot graph to compare models
  png(paste0(result_dir,"/",data_path,"_obsvspred.png"))
  plotObsVsPred(preds0[[data_path]])
  graphics.off()
  
  # png(paste0(result_dir,"/",data_path,"_rmse.png"))
  # dotplot(caret::resamples(preds0[[data_path]]))
  # graphics.off()
}

# preds = unlist(preds0,recursive=F)
rmsedf = ldply(names(preds0), function(xi) {
  x = preds0[[xi]]
  ldply(unique(x$model), function(mi) {
    mii = x$model==mi
    pr = x$pred[mii]
    ob = x$pred[mii]
    data.frame(rmse=c(rmse(pr[1:length(tr_ind)], ctr),
                      rmse(pr[(length(tr_ind)+1):(length(tr_ind)+length(te_ind))], cte),
                      rmse(pr[1:length(tr_ind0)], meta$GA[append(tr_ind,te_ind)])),
               feature=rep(xi,3), model=rep(mi,3), type=c("train","test","all"))
  })
})

png(paste0(result_dir,"/rmse_test.png"))
pl = barchart(rmse~model, data=rmsedf[rmsedf$type=="test",,drop=F], groups=feature, 
              auto.key = list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="test rmse for each model grouped by feature type")
print(pl)
graphics.off()
png(paste0(result_dir,"/rmse_train.png"))
pl = barchart(rmse~model, data=rmsedf[rmsedf$type=="train",,drop=F], groups=feature, 
              auto.key = list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="train rmse for each model grouped by feature type")
print(pl)
graphics.off()
png(paste0(result_dir,"/rmse_trte.png"))
pl = barchart(rmse~model, data=rmsedf[rmsedf$type=="all",,drop=F], groups=feature, 
              auto.key = list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="all rmse for each model grouped by feature type")
print(pl)
graphics.off()




## check for outliers ------------------------------------
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

png(paste0(result_dir,"/diff_outliers.png"), width=5000)
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



## save final results of one model/feature ------------------
feature = "raw"
model = "enet"
# finalsol = llply(preds, function(x) x$pred[is.na(x$obs)])
class_final$GA = round(preds0[[feature]]$pred[is.na(preds0[[feature]]$obs) & preds0[[feature]]$model==model],1)
write.csv(class_final, file=paste0(result_dir,"/TeamX_SC1_prediction.csv.csv"))
