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
root = "~/projects/ai4all2019_bio"
setwd(root)

## input
input_dir = paste0(root,"/input")
meta_dir = paste0(input_dir,"/anoSC1_v11_nokey.csv") # meta data
class_dir = paste0(input_dir,"/TeamX_SC1_prediction.csv") # class
data_dir = paste0(root,"/data") # features directory

## ouput
result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)
pred_dir = 


## load packages
pkgs = c("stringr", "plyr", "caret")
pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
sapply(pkgs, require, character.only=T)


# ## temp data prep (rm bottom 10% var genes)
# m0 = t(get(load(paste0(input_dir,"/HTA20_RMA.RData"))))
# range(m0) # range: 0.9365626 14.2836285
# plot(density(m0))
# 
# require(Rfast)
# m0vars = colVars(m0)
# m = m0[,m0vars>quantile(m0vars,0.1)]
# save(m, file=paste0(data_dir,"/raw.Rdata"))


## load files
meta = read.csv(meta_dir)
class_final = read.csv(class_dir)
data_paths = list.files(data_dir)

## prep 10-fold cross validation indices
tr_ind0 = meta$Train==1
tr_ind0r = sample(which(tr_ind0))
cv_inds = split(tr_ind0r, cut(seq_along(tr_ind0r), 10, labels=F))
cv_class = llply(cv_inds, function(is) meta$GA[is])
cv10 = function(m0, cv_inds, cv_class, fun, ...) {
  llply (1:length(cv_inds), function(i) {
    mtr = m0[unlist(cv_inds[-i]),]
    ctr = unlist(cv_class[-i])
    mte = m0[cv_inds[[i]],]
    cte = fun(mtr, ctr, mte, ...)
    return(cte)
  })
}

## test regression models
result0 = llply(data_paths, function(data_path) {
  res = NULL
  
  m0 = get(load(paste0(data_dir,"/", data_path)))
  mte0 = m0[!tr_ind0,]
  
  # random forest; pilot paper used caret recursive feature selection + random forest
  res$rf = cv10(m0, cv_inds, cv_class, function(mtr, ctr, mte) {
    
    return(cte) # vector of test class prediction
  }, ...)
})

## evaluation
rmse = function(x,y) sqrt(mean((x-y)^2))
result = unlist(result0, recursive=F)
scores = ldply(names(result), function(xi) {
  x = result[[xi]]
  score = laply (1:length(x), function(i) rmse(x[[i]], cv_class[[i]]) )
  return(data.frame(feature.model=rep(xi,length(scores)),
                    rmse=score, cv=1:length(scores)))
})

png(paste0(result_dir, "/10cv.png"))
pl = barchart(rmse~feature.model, data=scores, groups=cv, 
              auto.key=list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="rmse scores over 10-fold cv's")
print(pl)
dev.off()


## check for outliers
xdiffs = ldply(names(result), function(xi) {
  xdiff = unlist(result[[xi]]) - meta$GA[tr_ind0r]
  return(data.frame(feature.model=rep(xi,length(xdiff)),
                    GAdiff=xdiff, sample=meta$SampleID[tr_ind0r]))
})

png(paste0(result_dir, "/10cv_sample.png"))
par(mfrow=c(1,2))
pl = barchart(GAdiff~sample, data=xdiffs, groups=feature.model, 
              auto.key=list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="prediction - ground truth, over sample")
print(pl)
pl = barchart(abs(GAdiff)~sample, data=xdiffs, groups=feature.model, 
              auto.key=list(columns=2),
              cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
              main="absolute prediction - ground truth, over sample")
print(pl)
dev.off()

