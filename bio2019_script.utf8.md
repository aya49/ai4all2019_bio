
<!-- rnb-text-begin -->

---
title: "Invent the Future"
output: html_notebook
---

# Input Files:
- feature: (samples x 32830/925032 gene/probeset) RNAseq counts (+extracted feature); data has been batch and count normalized
- Value: (367 train sample) gestational age 8-42 weeks

# Output:
- Value: (368 test sample) gestational age 8-42 weeks rounded to 1 decimal place

# Data Analysis 

##1. Preprocessing  

In the first part we will work on the environment we will use. This means we will clean the environment, declare the directory of the project and install libraries.  



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Next, we will load 3 input files:

1. Sample Annotation file: 

* SampleID: unique identifier of the sample (matching the name of the .CEL file in HTA20 folder, except for extension .CEL);
* GA: gestational age as determined by the last menstrual period and or ultrasound;
* Batch: the batch identifier;
* Set: name of the source dataset;
* Train: 1 for samples to be used for training, 0 for samples to be used for test;

2. RNASEQ Data: each row is a gene and each column a patient 

*  probeset: gene ID
*  pid: patient id 

3. Submission template for the competition


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyMgbG9hZCBpbnB1dCBmaWxlc1xuXG4jIHNhbXBsZSBBbm5vdGF0aW9uIGZpbGVcbm1ldGEgPSByZWFkLmNzdihwYXN0ZTAoaW5wdXRfZGlyLFwiL2Fub1NDMV92MTFfbm9rZXkuY3N2XCIpKVxuXG4jIFJOQVNFUSBkYXRhXG5kYXRhMCA9IHQoZ2V0KGxvYWQocGFzdGUwKGlucHV0X2RpcixcIi9IVEEyMF9STUEuUkRhdGFcIikpKSlcblxuIyBzdWJtaXNzaW9uIHRlbXBsYXRlXG5zdWJtaXNzaW9uID0gcmVhZC5jc3YocGFzdGUwKGlucHV0X2RpcixcIi9UZWFtWF9TQzFfcHJlZGljdGlvbi5jc3ZcIikpIFxuYGBgIn0= -->

```r
## load input files

# sample Annotation file
meta = read.csv(paste0(input_dir,"/anoSC1_v11_nokey.csv"))

# RNASEQ data
data0 = t(get(load(paste0(input_dir,"/HTA20_RMA.RData"))))

# submission template
submission = read.csv(paste0(input_dir,"/TeamX_SC1_prediction.csv")) 
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Data Exploration: 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY2F0KCdSYW5nZVxcbicpO3JhbmdlKGRhdGEwKTtcbmNhdCgnR2VuZXMgSURzXFxuJyk7IGdpZCA9IGNvbG5hbWVzKGRhdGEwKTsgaGVhZChnaWQpO1xuY2F0KCdQYXRpZW50cyBJRHNcXG4nKTsgcGlkID0gcm93bmFtZXMoZGF0YTApOyBoZWFkKHBpZClcbmNhdCgnTWV0YSBEYXRhIGhlYWQgYW5kIHNoYXBlXFxuJyk7IGhlYWQobWV0YSk7IGRpbShtZXRhKVxuY2F0KCdSTkFTRVEgaGVhZCBhbmQgc2hhcGVcXG4nKTsgaGVhZChkYXRhMFssYygxOjUpXSk7IGRpbShkYXRhMClcbmNhdCgnU3VibWlzc2lvbiBoZWFkIGFuZCBzaGFwZVxcbicpOyBoZWFkKHN1Ym1pc3Npb24pOyBkaW0oc3VibWlzc2lvbilcbmBgYCJ9 -->

```r
cat('Range\n');range(data0);
cat('Genes IDs\n'); gid = colnames(data0); head(gid);
cat('Patients IDs\n'); pid = rownames(data0); head(pid)
cat('Meta Data head and shape\n'); head(meta); dim(meta)
cat('RNASEQ head and shape\n'); head(data0[,c(1:5)]); dim(data0)
cat('Submission head and shape\n'); head(submission); dim(submission)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Now we will explore the data with some plots. First we need to prepare the data and summarize the informatoin we will use in the plots. 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBwbG90IHN0YXRzOiBtZWFuIGNvdW50LCBwZWFyc29uL3NwZWFybWFuIGNvcnJcblxuIyB1c2Ugb25seSB0aGUgdHJhaW4gZGF0YXNldCB0byBjYWxjdWxhdGUgdGhlIGNvcnJlbGF0aW9uLCBtZWFuIGFuZCB2YXJpYW5jZSBcbmRhdGExID0gZGF0YS5mcmFtZSgnU2FtcGxlSUQnID0gcm93bmFtZXMoZGF0YTApLCBkYXRhMClcbmRhdGExID0gbWVyZ2UobWV0YVssYygxLDIsNSldLGRhdGExLGJ5LnggPSAnU2FtcGxlSUQnLGJ5LnkgPSAnU2FtcGxlSUQnLCBhbGwgPSBUKVxuXG4jIHNwbGl0IHRyYWluaW5nIHNldCBhbmQgdGVzdGluZyBzZXQgXG5kYXRhMSA9IHN1YnNldChkYXRhMSwgVHJhaW4gPT0gMSkgXG5kYXRhMSA9IHN1YnNldChkYXRhMSwgc2VsZWN0ID0gLWMoVHJhaW4sU2FtcGxlSUQpKVxuXG4jdHJhaW4gc2V0IFxuZGF0YV9leHBsID0gZGF0YS5mcmFtZShjb2wgPSBuYW1lcyhkYXRhMSksY29ycl9wID0gYXBwbHkoZGF0YTEsIDIsIGNvciwgeCA9IGRhdGExJEdBLCBtZXRob2QgPSAncGVhcnNvbicpKVxuZGF0YV9leHBsID0gZGF0YS5mcmFtZShkYXRhX2V4cGwsIGNvcnJfcyA9IGFwcGx5KGRhdGExLCAyLCBjb3IsIHggPSBkYXRhMSRHQSwgbWV0aG9kID0gJ3NwZWFybWFuJykpXG5yb3duYW1lcyhkYXRhX2V4cGwpID0gTlVMTCBcblxuXG4jRGF0YXNldCB3aXRoIGFsbCBpbmZvcm1hdGlvbiBcbmRhdGFfZXhwbCA9IGRhdGEuZnJhbWUoZGF0YV9leHBsLCB2YXJpYW5jZSA9IGModmFyKGRhdGExJEdBKSxjb2xWYXJzKGRhdGEwKSksIG1lYW4gPSBjKG1lYW4oZGF0YTEkR0EpLGNvbE1lYW5zKGRhdGEwKSkpIFxuZGF0YV9leHBsID0gZGF0YV9leHBsWy0xLF0gIyByZW1vdmluZyBHQSBmcm9tIHRoZSBkYXRhc2V0XG5gYGAifQ== -->

```r
# plot stats: mean count, pearson/spearman corr

# use only the train dataset to calculate the correlation, mean and variance 
data1 = data.frame('SampleID' = rownames(data0), data0)
data1 = merge(meta[,c(1,2,5)],data1,by.x = 'SampleID',by.y = 'SampleID', all = T)

# split training set and testing set 
data1 = subset(data1, Train == 1) 
data1 = subset(data1, select = -c(Train,SampleID))

#train set 
data_expl = data.frame(col = names(data1),corr_p = apply(data1, 2, cor, x = data1$GA, method = 'pearson'))
data_expl = data.frame(data_expl, corr_s = apply(data1, 2, cor, x = data1$GA, method = 'spearman'))
rownames(data_expl) = NULL 


#Dataset with all information 
data_expl = data.frame(data_expl, variance = c(var(data1$GA),colVars(data0)), mean = c(mean(data1$GA),colMeans(data0))) 
data_expl = data_expl[-1,] # removing GA from the dataset
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


To make the plots, we will use a library called ggplot2. Here are shown some plots we can make using this library. 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBwbG90IHN0YXRzXG5wMSA8LSBnZ3Bsb3QoZGF0YV9leHBsLCBhZXMoeD1tZWFuKSkgKyBcbiAgZ2VvbV9oaXN0b2dyYW0oZmlsbCA9ICdsaWdodGdyZWVuJykgKyBcbiAgeGxhYignQXZlcmFnZSBFeHByZXNzaW9uJykgKyBsYWJzKHRpdGxlPScoYSknKVxuXG5wMiA8LSBnZ3Bsb3QoZGF0YV9leHBsLCBhZXMoeD1tZWFuLCB5ID0gdmFyaWFuY2UpKSArIFxuICBnZW9tX3BvaW50KGNvbG9yID0gJ2xpZ2h0Z3JlZW4nKStcbiAgeGxhYignQXZlcmFnZSBFeHByZXNzaW9uJykgKyBcbiAgeWxhYignVmFyaWFuY2UnKSsgbGFicyh0aXRsZT0nKGIpJylcblxucDMgPC0gZ2dwbG90KGRhdGFfZXhwbCwgYWVzKHg9Y29ycl9wKSkgKyBcbiAgZ2VvbV9oaXN0b2dyYW0oZmlsbCA9ICdsaWdodGdyZWVuJykrXG4gIHhsYWIoJ1BlYXJzb24gQ29ycmVsYXRpb24gYmV0d2VlbiBHZW5lcyBhbmQgR2VzdGFjaW9uYWwgQWdlJykrIGxhYnModGl0bGU9JyhjKScpXG5cbnA0IDwtIGdncGxvdChkYXRhX2V4cGwsIGFlcyh4PWNvcnJfcCwgeSA9IGNvcnJfcykpICsgXG4gIGdlb21fcG9pbnQoY29sb3IgPSAnbGlnaHRncmVlbicpKyBsYWJzKHRpdGxlPScoZCknKStcbiAgeGxhYignUGVhcnNvbiBDb3JyZWxhdGlvbicpICsgeWxhYignU3BlYXJtYW4gQ29ycmVsYXRpb24nKVxuXG5ncmlkLmFycmFuZ2UocDEsIHAyLCBwMywgcDQsIG5yb3cgPSAyKVxuYGBgIn0= -->

```r
# plot stats
p1 <- ggplot(data_expl, aes(x=mean)) + 
  geom_histogram(fill = 'lightgreen') + 
  xlab('Average Expression') + labs(title='(a)')

p2 <- ggplot(data_expl, aes(x=mean, y = variance)) + 
  geom_point(color = 'lightgreen')+
  xlab('Average Expression') + 
  ylab('Variance')+ labs(title='(b)')

p3 <- ggplot(data_expl, aes(x=corr_p)) + 
  geom_histogram(fill = 'lightgreen')+
  xlab('Pearson Correlation between Genes and Gestacional Age')+ labs(title='(c)')

p4 <- ggplot(data_expl, aes(x=corr_p, y = corr_s)) + 
  geom_point(color = 'lightgreen')+ labs(title='(d)')+
  xlab('Pearson Correlation') + ylab('Spearman Correlation')

grid.arrange(p1, p2, p3, p4, nrow = 2)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


### Exercise 1: 
1. What can we conclude from these plots above? 
2. Why we don't have the correlation between the gene expression and GA for the testing set? 
3. Can you think in any other intesresting plots or analysis? 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI2luc2VydCBoZXJlIGEgbmV3IHBsb3Qgb3IgYW5hbHlzaXMgeW91IHRoaW5rIGlzIHdvcnRoIGluc3RlcmluZ1xuXG5gYGAifQ== -->

```r
#insert here a new plot or analysis you think is worth instering

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




Answer (delete): 
1. a) The average expression is around 5. There are a few genes that are more expressed. b) There is no strong correlation between variance and average gene expression. c) The person correlation between the genes and the Gestacional Age looks like a normal distributio with mean equals 0. This means that most of genes will have a association equals 0. d) The correlation between the Spearman Correlation and Pearson correlation is very strong, meaning that we can use either Spearman or Pearson that the results will be very similar. 
2. We don't make these correlation because we don't know the GA for the testing set. In fact, this is the exacly information that we are looking for. 

##2. Features Extraction 

The original dataset has 32830 genes expressions. However, from the plot of the correlation we already know that some of these genes aren't associate with GA. Besides, genes with a very low variance might don't contribute for the prediction of gestacional age. Also, large datasets might add noise to machine leanring models. 

We will use 5 methods do extract feature. 
1. Elimination of 30% of genes with lowest variance; 
2. Elimination of 30% of genes with lowest correlation with GA; 
3. PCA
4. Random Forest 
5. Autoenconder 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyMgMykgUENBXG5QcmVQQ0EgPSBwcmVQcm9jZXNzKGRhdGEyLG1ldGhvZD1cInBjYVwiKVxuZmVhdC5wY2EgPSBwcmVkaWN0KFByZVBDQSxkYXRhMilcbndyaXRlLmNzdihmZWF0LnBjYSwgcGFzdGUwKGZlYXRfZGlyLCcvZmVhdHVyZV90cnBjYS5jc3YnKSwgcm93Lm5hbWVzPVQpXG5cbmBgYCJ9 -->

```r
## 3) PCA
PrePCA = preProcess(data2,method="pca")
feat.pca = predict(PrePCA,data2)
write.csv(feat.pca, paste0(feat_dir,'/feature_trpca.csv'), row.names=T)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyMgNCkgUmFuZG9tIEZvcmVzdCBcbmRhdGExID0gc3Vic2V0KGRhdGExLCBzZWxlY3QgPWtlZXAkY29sKVxuZGF0YTEgPSBkYXRhLmZyYW1lKEdBMSwgZGF0YTEpXG5pZighZXhpc3RzKCdQcmVSRicpKXtcbiAgUHJlUkYgPSBjYXJldDo6dHJhaW4oeT1HQTEsIHg9ZGF0YTFbLC0xXSwgIG1ldGhvZD0ncmFuZ2VyJyxpbXBvcnRhbmNlPSdpbXB1cml0eScpXG4gIFByZVJGLmkgPSB2YXJJbXAoUHJlUkYpJGltcG9ydGFuY2Vcbn1lbHNle1xuICBQcmVSRi5pID0gdmFySW1wKFByZVJGKSRpbXBvcnRhbmNlXG59XG5cbmZlYXQucmEgPSBkYXRhMlssb3JkZXIoUHJlUkYuaSwgZGVjcmVhc2luZz1UKVsxOjUwMF1dXG53cml0ZS5jc3YoZmVhdC5yYSwgcGFzdGUwKGZlYXRfZGlyLCcvZmVhdHVyZV90cnJhLmNzdicpLCByb3cubmFtZXMgPSBUKVxuXG5gYGAifQ== -->

```r
## 4) Random Forest 
data1 = subset(data1, select =keep$col)
data1 = data.frame(GA1, data1)
if(!exists('PreRF')){
  PreRF = caret::train(y=GA1, x=data1[,-1],  method='ranger',importance='impurity')
  PreRF.i = varImp(PreRF)$importance
}else{
  PreRF.i = varImp(PreRF)$importance
}

feat.ra = data2[,order(PreRF.i, decreasing=T)[1:500]]
write.csv(feat.ra, paste0(feat_dir,'/feature_trra.csv'), row.names = T)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyMgNSkgQXV0b2VuY29kZXIgXG5pZighZXhpc3RzKCdwcmVBJykpe1xuICBwcmVBID0gYXV0b2VuY29kZXIoZGF0YTIsIGhpZGRlbi5sYXllcnMgPSBjKDEwMDAsIDUwMCwgMTAwMCkpXG4gIGZlYXQuQSA9IGVuY29kZShwcmVBLCBkYXRhMikgIFxufVxucm93bmFtZXMoZmVhdC5BKSA9IHJvd25hbWVzKGZlYXQucmEpXG53cml0ZS5jc3YoZmVhdC5BLCBwYXN0ZTAoZmVhdF9kaXIsJy9mZWF0dXJlX3RyYS5jc3YnKSwgcm93Lm5hbWVzID0gVClcbmBgYCJ9 -->

```r
## 5) Autoencoder 
if(!exists('preA')){
  preA = autoencoder(data2, hidden.layers = c(1000, 500, 1000))
  feat.A = encode(preA, data2)  
}
rownames(feat.A) = rownames(feat.ra)
write.csv(feat.A, paste0(feat_dir,'/feature_tra.csv'), row.names = T)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


### Exercise 2: 
1. What parameters do you think you could change? 
2. Choose between PCA and Random Forest. Try to save a new feature set with a different number of feature. 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI2ZlYXQudW5pcXVlID0gcHJlZGljdChQcmVQQ0EsZGF0YTIpICNmb3IgcGNhIE9SXG4jZmVhdC51bmlxdWUgPSBkYXRhMltdICNmb3IgcmFuZG9tIGZvcmVzdFxuI3dyaXRlLmNzdihmZWF0LnVuaXF1ZSwgKVxuXG5gYGAifQ== -->

```r
#feat.unique = predict(PrePCA,data2) #for pca OR
#feat.unique = data2[] #for random forest
#write.csv(feat.unique, )

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



##3. Machine Learning models

In this section we will work which the Machine Learning Models. We have 5 different sets of feature (RAW, PCA, RA, A, UNIQUE) and we will first work with 3 models: Linear Regression, RVM and Random Forest. The best model will be selected based on the difference between the values we observed on the training set for GA (Gestacional Age) and the predicted values for the combination of MODEL+FEATURE. The lowest is this difference on the validation set, the better. 

Before we begin to test the models, we need to set up an experimental framework to so that we can evaluate our results once they come out. One such technique is a resampling procedure called ```cvn```-fold cross validation (here we set ```cvn``` to 10).

Typically, data sets are split up into a train and test set; we extract these indices into the ```tr_ind``` and ```te_ind``` variables respectively. The former set of examples is used to train our model, which is later tested on the test set, a set of samples our model has not seen before. However, since we assume we do not have the GA for our test set, we would not know how our model performs on the test set.

Therefore, we further split the train set into ```cvn```=10 equal sized chunks. Since we know the Ga to all our train samples, we can evaluate the model a total of 10 times, each time training the model 9 of the chunks and testing the model on 1 of the chunks. Thereafter, these metrics can be combined (e.g. via a mean) to produce an appropriate evaluation of the model.

### Exercise 3:
1. Why do we want to perform cross validation? Why would we cross validating the model ```cvn``` times over just 1 time?

Answer (delete): 
1. Evaluation becomes more robust with resampling



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyMgMSkgcHJlcCBjdm4tZm9sZCBjcm9zcyB2YWxpZGF0aW9uICYgcm1zZSBmdW5jdGlvblxuY3ZpbmRzX3BhdGggPSBwYXN0ZTAocm9vdCxcIi9jdmluZHMuUmRhdGFcIilcbmlmKCFmaWxlLmV4aXN0cyhjdmluZHNfcGF0aCkpe1xuICBjdm4gPSAxMFxuICB0cmFpbl9pbmRleCA9IHdoaWNoKG1ldGEkVHJhaW49PTEpICNzZWxlY3Rpbmcgb25seSB0aGUgdHJhaW5pbmcgZXhhbXBsZXNcbiAgdGVzdF9pbmRleCA9IHdoaWNoKG1ldGEkVHJhaW49PTApXG4gIHRyYWluX2luZGV4X3ZhbCA9IHNhbXBsZSh0cmFpbl9pbmRleCwgY2VpbGluZyhsZW5ndGgodHJhaW5faW5kZXgpLzExKSkgIyBjcm9zcyB2YWxpZGF0aW9uIHNldFxuICB0cmFpbl9pbmRleF90ciA9IHNhbXBsZSh0cmFpbl9pbmRleFshdHJhaW5faW5kZXglaW4ldHJhaW5faW5kZXhfdmFsXSkgIyByZW1vdmUgY3Jvc3MgdmFsaWRhdGlvbiBzZXQgZnJvbSB0cmFpbmluZyBzZXQgXG4gIGdhX3ZhbCA9IGFzLm51bWVyaWMobWV0YSRHQVt0cmFpbl9pbmRleF92YWxdKVxuICBnYV90ciA9IGFzLm51bWVyaWMobWV0YSRHQVt0cmFpbl9pbmRleF90cl0pIFxuICBzYXZlKGN2biwgdHJhaW5faW5kZXgsIHRyYWluX2luZGV4X3ZhbCwgdHJhaW5faW5kZXhfdHIsIGdhX3ZhbCwgZ2FfdHIsIGZpbGU9Y3ZpbmRzX3BhdGgpXG59XG5sb2FkKGN2aW5kc19wYXRoKVxuXG5gYGAifQ== -->

```r
## 1) prep cvn-fold cross validation & rmse function
cvinds_path = paste0(root,"/cvinds.Rdata")
if(!file.exists(cvinds_path)){
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

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


The RMSE (root mean squared error) is a metric that indicates how different our models' predicted GA are compared to the true GA. In other words, the smaller the RMSE, the better.This will be the metric used to compare the models we will be working on. 

###3.1 Linear Regression 

Remember $y=ax+b$? That's linear regression in a nutshell. Linear regression assumes that there is a linear relationship between the given multidimensional RNAseq train data $x$ and the corresponding variable we want to predict GA $y$. $a$ represents the slope and $b$ represents the y-intercept. In the real world though, there often isn't a perfect linear relationship, so below, we try to estimate the best line (i.e. $a$ and $b$ parameters) given $x$ and $y$.

Linear regression is fast and intuitive, but there are many data sets that do not conform with this linear assumption.


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZmVhdHVyZV90eXBlID0gJ3BjYScgI29wdGlvbnMgYXJlICdwY2EnLCAncmF3JywgJ3JhJyBmb3IgcmFuZG9tIGZvcmVzdCBhbmQgJ2EnIGZvciBhdXRvZW5jb2RlclxuIyMgMCkgbG9hZCBmZWF0dXJlIFxuZmVhdHVyZSA9IHJlYWQuY3N2KHBhc3RlKGZlYXRfZGlyLFwiL2ZlYXR1cmVfdHJcIixmZWF0dXJlX3R5cGUsJy5jc3YnLCBzZXAgPSAnJykpXG5yb3duYW1lcyhmZWF0dXJlKSA9IGZlYXR1cmVbLDFdXG5mZWF0dXJlID0gYXMubWF0cml4KGZlYXR1cmVbLC0xXSlcblxuI1NQbGl0aW5nIGZlYXR1cmUgb2YgdHJhaW5pbmcgYW5kIHZhbGlkYXRpb24gc2V0XG5mZWF0dXJlX3RyID0gZmVhdHVyZVt0cmFpbl9pbmRleF90cixdXG5mZWF0dXJlX3ZhbCA9IGZlYXR1cmVbdHJhaW5faW5kZXhfdmFsLF1cbnRyYWluXyA9IGRhdGEuZnJhbWUoZ2FfdHIsZmVhdHVyZV90cilcblxuI0xpbmVhciBSRWdyZXNzaW9uIG1vZGVsXG5tb2RlbDEgPSBsbShnYV90ciB+IC4sZGF0YSA9IHRyYWluXylcbiNQcmVkaWN0aW9uc1xucHJlZGljdGlvbnNfbTEgPSBwcmVkaWN0KG1vZGVsMSwgbmV3ZGF0YSA9IHRyYWluX1ssLTFdKVxuI0Vycm9yIGZyb20gdGhlIG9ic2VydmVkIHZhbHVlcyBhbmQgZml0dGVkIHZhbHVlc1xucm1zZShnYV90ciwgcHJlZGljdGlvbnNfbTEpXG5cbmBgYCJ9 -->

```r
feature_type = 'pca' #options are 'pca', 'raw', 'ra' for random forest and 'a' for autoencoder
## 0) load feature 
feature = read.csv(paste(feat_dir,"/feature_tr",feature_type,'.csv', sep = ''))
rownames(feature) = feature[,1]
feature = as.matrix(feature[,-1])

#SPliting feature of training and validation set
feature_tr = feature[train_index_tr,]
feature_val = feature[train_index_val,]
train_ = data.frame(ga_tr,feature_tr)

#Linear REgression model
model1 = lm(ga_tr ~ .,data = train_)
#Predictions
predictions_m1 = predict(model1, newdata = train_[,-1])
#Error from the observed values and fitted values
rmse(ga_tr, predictions_m1)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGxvdGRhdGEgPSBkYXRhLmZyYW1lKGdhX3RyLHByZWRpY3Rpb25zX20xKVxuZ2dwbG90KHBsb3RkYXRhLCBhZXMoeCA9IGdhX3RyLHkgPSBwcmVkaWN0aW9uc19tMSkpK1xuICBnZW9tX3BvaW50KGNvbG9yID0gJ2xpZ2h0Z3JlZW4nKStcbiAgeGxhYignTGluZWFyIFJlZ3Jlc3Npb24gTW9kZWwgLSBUcmFpbmluZyBTZXQnKStcbiAgeWxhYignUHJlZGljdGVkIFZhbHVlcycpXG5cbmBgYCJ9 -->

```r
plotdata = data.frame(ga_tr,predictions_m1)
ggplot(plotdata, aes(x = ga_tr,y = predictions_m1))+
  geom_point(color = 'lightgreen')+
  xlab('Linear Regression Model - Training Set')+
  ylab('Predicted Values')

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


### Exercise 4: 
Calculate the predictions for feature_val and ga_val using linear regression. With the predicted values, make a plot predictions_m1_val and ga_val

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3ByZWRpY3Rpb25zX20xX3ZhbCA9IHByZWRpY3QoKVxuI3Jtc2UoKVxuYGBgIn0= -->

```r
#predictions_m1_val = predict()
#rmse()
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3Bsb3RkYXRhID0gZGF0YS5mcmFtZSgpXG4jZ2dwbG90KHBsb3RkYXRhLCBhZXMoeCA9ICx5ID0gKSkrZ2VvbV9wb2ludChjb2xvciA9ICcnKSt4bGFiKCcnKSt5bGFiKCcnKVxuXG5gYGAifQ== -->

```r
#plotdata = data.frame()
#ggplot(plotdata, aes(x = ,y = ))+geom_point(color = '')+xlab('')+ylab('')

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




###3.2 RVM


RVM (relevance vector machines) is an application of the Bayesian treatment of general linear models to SVMs (support vector machines). 

SVMs is a model that draws support vectors (i.e. a line) to separate data in an effort to predict their corresponding discrete classes (e.g. early GA, late GA) as opposed to predicting continuous values (e.g. 10 week GA, 32 week GA). SVMs extract this line from the train data by first calculating a distance between all the data points. These distances are used to define a best line, one that not only separates the data of different classes, but are also the largest distance from data points of both classes (i.e. the line sites in the middle of the space between the two classes, rather than being closer to one or the other).

RVM is functionally identicle to SVM. The Bayesian treatment is a fancy way of saying that instead of categorizing all data points on one side of the line as a single class, let's give them continuous class values that indicates how far the data points are from the line.

In both cases, the user defined formua that calculates distances between data points are called "kernels".

RVMs are fast to train and easy to use, but it can become less effective if the data set used has too many feature because it depends on possibly dimension sensative kernel functions.


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5tb2RlbDIgPSBydm0oeCA9IGZlYXR1cmVfdHIsIGdhX3RyLCB0eXBlPVwicmVncmVzc2lvblwiKVxucHJlZGljdGlvbnNfbTIgPSBwcmVkaWN0KG1vZGVsMixkYXRhID0gYXMuZGF0YS5mcmFtZShmZWF0dXJlX3RyKSlcbnJtc2UoZ2FfdHIscHJlZGljdGlvbnNfbTIpXG5cblxucGxvdGRhdGEgPSBkYXRhLmZyYW1lKGdhX3RyLHByZWRpY3Rpb25zX20yKVxuZ2dwbG90KHBsb3RkYXRhLCBhZXMoeCA9IGdhX3RyLHkgPSBwcmVkaWN0aW9uc19tMikpK1xuICBnZW9tX3BvaW50KGNvbG9yID0gJ2xpZ2h0Z3JlZW4nKStcbiAgeWxhYignUHJlZGljdGVkIFZhbHVlcycpK1xuICB4bGFiKCdSVk0gLSBUcmFpbmluZyBTZXQnKVxuXG5gYGAifQ== -->

```r

model2 = rvm(x = feature_tr, ga_tr, type="regression")
predictions_m2 = predict(model2,data = as.data.frame(feature_tr))
rmse(ga_tr,predictions_m2)


plotdata = data.frame(ga_tr,predictions_m2)
ggplot(plotdata, aes(x = ga_tr,y = predictions_m2))+
  geom_point(color = 'lightgreen')+
  ylab('Predicted Values')+
  xlab('RVM - Training Set')

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

### Exercise 5: 
Calculate the predictions for feature_val and ga_val using rvm. With the predicted values, make a plot predictions_m1_val and ga_val

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3ByZWRpY3Rpb25zX20yX3ZhbCA9IHByZWRpY3QoKVxuI3Jtc2UoKVxuYGBgIn0= -->

```r
#predictions_m2_val = predict()
#rmse()
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3Bsb3RkYXRhID0gZGF0YS5mcmFtZSgpXG4jZ2dwbG90KHBsb3RkYXRhLCBhZXMoeCA9ICx5ID0gKSkrZ2VvbV9wb2ludChjb2xvciA9ICcnKSt4bGFiKCcnKSt5bGFiKCcnKVxuXG5gYGAifQ== -->

```r
#plotdata = data.frame()
#ggplot(plotdata, aes(x = ,y = ))+geom_point(color = '')+xlab('')+ylab('')

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



###3.3 Random Forest

A random forest is an aggregate of decision trees formed via model training. A decision tree is a decision making template in the form of "if this then that"; in our case, this decision would be "if this gene is expressed often, then the woman is in her 32nd week of gestation". 

Decision trees and the rules they contain are interpretable and easy to understand. As well combining many decision trees together yield more robust results. However, it can be slow if there are too many feature and it is more suited for predicting categorical results.


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubW9kZWwzID0gcmFuZ2VyKGdhX3Ryfi4sZGF0YSA9IHRyYWluXylcbnByZWRpY3Rpb25zX20zID0gcHJlZGljdChtb2RlbDMsIGRhdGEgPSBhcy5kYXRhLmZyYW1lKGZlYXR1cmVfdHIpKSBcbnJtc2UoZ2FfdHIscHJlZGljdGlvbnNfbTMkcHJlZGljdGlvbnMpXG5cbnBsb3RkYXRhID0gZGF0YS5mcmFtZShnYV90cixwcmVkaWN0aW9uc19tMyRwcmVkaWN0aW9ucylcbmdncGxvdChwbG90ZGF0YSwgYWVzKHggPSBnYV90cix5ID0gcHJlZGljdGlvbnNfbTMucHJlZGljdGlvbnMpKStcbiAgZ2VvbV9wb2ludChjb2xvciA9ICdsaWdodGdyZWVuJykrXG4gIHlsYWIoJ1ByZWRpY3RlZCBWYWx1ZXMnKStcbiAgeGxhYignUmFuZG9tIEZvcmVzdCAtIFRyYWluaW5nIFNldCcpXG5cbmBgYCJ9 -->

```r
model3 = ranger(ga_tr~.,data = train_)
predictions_m3 = predict(model3, data = as.data.frame(feature_tr)) 
rmse(ga_tr,predictions_m3$predictions)

plotdata = data.frame(ga_tr,predictions_m3$predictions)
ggplot(plotdata, aes(x = ga_tr,y = predictions_m3.predictions))+
  geom_point(color = 'lightgreen')+
  ylab('Predicted Values')+
  xlab('Random Forest - Training Set')

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

### Exercise 6: 
Calculate the predictions for feature_val and ga_val using Random Forest. With the predicted values, make a plot predictions_m1_val and ga_val

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3ByZWRpY3Rpb25zX20zX3ZhbCA9IHByZWRpY3QoKVxuI3Jtc2UoKVxuYGBgIn0= -->

```r
#predictions_m3_val = predict()
#rmse()
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3Bsb3RkYXRhID0gXG4jZ2dwbG90KHBsb3RkYXRhLCBhZXMoeCA9ICx5ID0gKSkrZ2VvbV9wb2ludChjb2xvciA9ICcnKSt4bGFiKCcnKSt5bGFiKCcnKVxuXG5gYGAifQ== -->

```r
#plotdata = 
#ggplot(plotdata, aes(x = ,y = ))+geom_point(color = '')+xlab('')+ylab('')

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



### Exercise 7:
Repite the Linear Regression, SCM and Random Forest with 2 other sets of feature: 'raw', 'ra' for random forest or 'a' for autoencoder and the 'unique' set of feature you made on  Exercise 2. From the set of feature you choose, which one have the best result? Why?

SET OF FEATURES 1

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIy0tLS0tLS0tLS0tLS0tLS0tLS0tLS0gICBDSEFOR0UgSEVSRVxuI2ZlYXR1cmVfdHlwZSA9ICAjb3B0aW9ucyBhcmUgJ3BjYScsICdyYXcnLCAncmEnIGZvciByYW5kb20gZm9yZXN0IGFuZCAnYScgZm9yIGF1dG9lbmNvZGVyXG5cbiMjIDApIGxvYWQgZmVhdHVyZSBcbmZlYXR1cmUgPSByZWFkLmNzdihwYXN0ZShmZWF0X2RpcixcIi9mZWF0dXJlX3RyXCIsZmVhdHVyZV90eXBlLCcuY3N2Jywgc2VwID0gJycpKVxucm93bmFtZXMoZmVhdHVyZSkgPSBmZWF0dXJlWywxXVxuZmVhdHVyZSA9IGFzLm1hdHJpeChmZWF0dXJlWywtMV0pXG5mZWF0dXJlX3RyID0gZmVhdHVyZVt0cmFpbl9pbmRleF90cixdXG5mZWF0dXJlX3ZhbCA9IGZlYXR1cmVbdHJhaW5faW5kZXhfdmFsLF1cbnRyYWluXyA9IGRhdGEuZnJhbWUoZ2FfdHIsZmVhdHVyZV90cilcblxuIy0tLS0tLS0tLS0tLS0tLS0tLS0tLS0gICBMaW5lYXIgUmVncmVzc2lvblxubW9kZWwxID0gbG0oZ2FfdHIgfiAuLGRhdGEgPSB0cmFpbl8pXG5wcmVkaWN0aW9uc19tMSA9IHByZWRpY3QobW9kZWwxLCBuZXdkYXRhID0gdHJhaW5fWywtMV0pXG5ybXNlKGdhX3RyLCBwcmVkaWN0aW9uc19tMSlcblxuI1BMT1QgXG5wbG90ZGF0YSA9IGRhdGEuZnJhbWUoZ2FfdHIscHJlZGljdGlvbnNfbTEpXG5nZ3Bsb3QocGxvdGRhdGEsIGFlcyh4ID0gZ2FfdHIseSA9IHByZWRpY3Rpb25zX20xKSkrXG4gIGdlb21fcG9pbnQoY29sb3IgPSAnbGlnaHRncmVlbicpK1xuICB4bGFiKCdMaW5lYXIgUmVncmVzc2lvbiBNb2RlbCAtIFRyYWluaW5nIFNldCcpK1xuICB5bGFiKCdQcmVkaWN0ZWQgVmFsdWVzJylcblxuIy0tLS0tLS0tLS0tLS0tLS0tLS0tLS0gICBDSEFOR0UgSEVSRTogQUREIExJTkVBUiBSRUdSRVNTSU9OIEZPUiBUSEUgVkFMSURBVElPTiBTRVQgXG5cbiMtLS0tLS0tLS0tLS0tLS0tLS0tLS0tICAgQUREIFNWTSBGT1IgVFJBSU5JTkcgQU5EIFZBTElEQVRJT04gU0VUIFxuXG4jLS0tLS0tLS0tLS0tLS0tLS0tLS0tLSAgIEFERCBSQU0gRk9SIFRSQUlOSU5HIEFORCBWQUxJREFUSU9OIFNFVCBcblxuXG5cblxuXG5cbmBgYCJ9 -->

```r
#----------------------   CHANGE HERE
#feature_type =  #options are 'pca', 'raw', 'ra' for random forest and 'a' for autoencoder

## 0) load feature 
feature = read.csv(paste(feat_dir,"/feature_tr",feature_type,'.csv', sep = ''))
rownames(feature) = feature[,1]
feature = as.matrix(feature[,-1])
feature_tr = feature[train_index_tr,]
feature_val = feature[train_index_val,]
train_ = data.frame(ga_tr,feature_tr)

#----------------------   Linear Regression
model1 = lm(ga_tr ~ .,data = train_)
predictions_m1 = predict(model1, newdata = train_[,-1])
rmse(ga_tr, predictions_m1)

#PLOT 
plotdata = data.frame(ga_tr,predictions_m1)
ggplot(plotdata, aes(x = ga_tr,y = predictions_m1))+
  geom_point(color = 'lightgreen')+
  xlab('Linear Regression Model - Training Set')+
  ylab('Predicted Values')

#----------------------   CHANGE HERE: ADD LINEAR REGRESSION FOR THE VALIDATION SET 

#----------------------   ADD SVM FOR TRAINING AND VALIDATION SET 

#----------------------   ADD RAM FOR TRAINING AND VALIDATION SET 

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



SET OF FEATURES unique

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIy0tLS0tLS0tLS0tLS0tLS0tLS0tLS0gICBDSEFOR0UgSEVSRVxuI2ZlYXR1cmVfdHlwZSA9ICAjb3B0aW9ucyBhcmUgJ3BjYScsICdyYXcnLCAncmEnIGZvciByYW5kb20gZm9yZXN0IGFuZCAnYScgZm9yIGF1dG9lbmNvZGVyXG5cbiMjIDApIGxvYWQgZmVhdHVyZSBcbmZlYXR1cmUgPSByZWFkLmNzdihwYXN0ZShmZWF0X2RpcixcIi9mZWF0dXJlX3RyXCIsZmVhdHVyZV90eXBlLCcuY3N2Jywgc2VwID0gJycpKVxucm93bmFtZXMoZmVhdHVyZSkgPSBmZWF0dXJlWywxXVxuZmVhdHVyZSA9IGFzLm1hdHJpeChmZWF0dXJlWywtMV0pXG5mZWF0dXJlX3RyID0gZmVhdHVyZVt0cmFpbl9pbmRleF90cixdXG5mZWF0dXJlX3ZhbCA9IGZlYXR1cmVbdHJhaW5faW5kZXhfdmFsLF1cbnRyYWluXyA9IGRhdGEuZnJhbWUoZ2FfdHIsZmVhdHVyZV90cilcblxuIy0tLS0tLS0tLS0tLS0tLS0tLS0tLS0gICBMaW5lYXIgUmVncmVzc2lvblxubW9kZWwxID0gbG0oZ2FfdHIgfiAuLGRhdGEgPSB0cmFpbl8pXG5wcmVkaWN0aW9uc19tMSA9IHByZWRpY3QobW9kZWwxLCBuZXdkYXRhID0gdHJhaW5fWywtMV0pXG5ybXNlKGdhX3RyLCBwcmVkaWN0aW9uc19tMSlcblxuI1BMT1QgXG5wbG90ZGF0YSA9IGRhdGEuZnJhbWUoZ2FfdHIscHJlZGljdGlvbnNfbTEpXG5nZ3Bsb3QocGxvdGRhdGEsIGFlcyh4ID0gZ2FfdHIseSA9IHByZWRpY3Rpb25zX20xKSkrXG4gIGdlb21fcG9pbnQoY29sb3IgPSAnbGlnaHRncmVlbicpK1xuICB4bGFiKCdMaW5lYXIgUmVncmVzc2lvbiBNb2RlbCAtIFRyYWluaW5nIFNldCcpK1xuICB5bGFiKCdQcmVkaWN0ZWQgVmFsdWVzJylcblxuIy0tLS0tLS0tLS0tLS0tLS0tLS0tLS0gICBDSEFOR0UgSEVSRTogQUREIExJTkVBUiBSRUdSRVNTSU9OIEZPUiBUSEUgVkFMSURBVElPTiBTRVQgXG5cbiMtLS0tLS0tLS0tLS0tLS0tLS0tLS0tICAgQUREIFNWTSBGT1IgVFJBSU5JTkcgQU5EIFZBTElEQVRJT04gU0VUIFxuXG4jLS0tLS0tLS0tLS0tLS0tLS0tLS0tLSAgIEFERCBSQU0gRk9SIFRSQUlOSU5HIEFORCBWQUxJREFUSU9OIFNFVCBcblxuXG5cblxuXG5cbmBgYCJ9 -->

```r
#----------------------   CHANGE HERE
#feature_type =  #options are 'pca', 'raw', 'ra' for random forest and 'a' for autoencoder

## 0) load feature 
feature = read.csv(paste(feat_dir,"/feature_tr",feature_type,'.csv', sep = ''))
rownames(feature) = feature[,1]
feature = as.matrix(feature[,-1])
feature_tr = feature[train_index_tr,]
feature_val = feature[train_index_val,]
train_ = data.frame(ga_tr,feature_tr)

#----------------------   Linear Regression
model1 = lm(ga_tr ~ .,data = train_)
predictions_m1 = predict(model1, newdata = train_[,-1])
rmse(ga_tr, predictions_m1)

#PLOT 
plotdata = data.frame(ga_tr,predictions_m1)
ggplot(plotdata, aes(x = ga_tr,y = predictions_m1))+
  geom_point(color = 'lightgreen')+
  xlab('Linear Regression Model - Training Set')+
  ylab('Predicted Values')

#----------------------   CHANGE HERE: ADD LINEAR REGRESSION FOR THE VALIDATION SET 

#----------------------   ADD SVM FOR TRAINING AND VALIDATION SET 

#----------------------   ADD RAM FOR TRAINING AND VALIDATION SET 

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


###Exercise 8: 
Create a plot to identify which model or set of feature is the best between so far. 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI1NhdmluZyB0aGUgcmVzdWx0cyBvZiB0aGUgcHJlZGljdGlvbnMgd2hlbiB0aGUgLnVuaXF1ZSBmZWF0dXJlIHdlcmUgdXNlZCBcbm91dHB1dF91bmlxdWUgPSBkYXRhLmZyYW1lKFxuICBtb2RlbCA9IGMoJ0xSJywnUlZNJywnUkYnKSxcbiAgZmVhdHVyZSA9IGMocmVwKCd1bmlxdWUnLDMpICksXG4gIHJtc2VfdCA9IGMocm1zZShnYV90ciwgcHJlZGljdGlvbnNfbTEpLHJtc2UoZ2FfdHIscHJlZGljdGlvbnNfbTIpLHJtc2UoZ2FfdHIscHJlZGljdGlvbnNfbTMkcHJlZGljdGlvbnMpICksXG4gIHJtc2VfdmFsID0gYyhybXNlKGdhX3ZhbCwgcHJlZGljdGlvbnNfbTF2KSxybXNlKGdhX3ZhbCxwcmVkaWN0aW9uc19tMnYpLHJtc2UoZ2FfdmFsLHByZWRpY3Rpb25zX20zdiRwcmVkaWN0aW9ucykgKVxuKVxuXG5vdXRwdXQgPSByZWFkLnRhYmxlKHBhc3RlMChpbnB1dF9kaXIsJy8vb3V0cHV0X2Jhc2ljX21vZGVscy5jc3YnKSxzZXA9JywnLCBoZWFkZXI9VClcbm91dHB1dCA9IHJiaW5kKG91dHB1dCxvdXRwdXRfdW5pcXVlKVxuYGBgIn0= -->

```r
#Saving the results of the predictions when the .unique feature were used 
output_unique = data.frame(
  model = c('LR','RVM','RF'),
  feature = c(rep('unique',3) ),
  rmse_t = c(rmse(ga_tr, predictions_m1),rmse(ga_tr,predictions_m2),rmse(ga_tr,predictions_m3$predictions) ),
  rmse_val = c(rmse(ga_val, predictions_m1v),rmse(ga_val,predictions_m2v),rmse(ga_val,predictions_m3v$predictions) )
)

output = read.table(paste0(input_dir,'//output_basic_models.csv'),sep=',', header=T)
output = rbind(output,output_unique)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI3A3PC0gZ2dwbG90KGRhdGEgPSBvdXRwdXQsIGFlcyh4ID0gLCB5ID0gKSkgK1xuI3A4PC0gZ2dwbG90KGRhdGEgPSBvdXRwdXQsIGFlcyh4ID0gLCB5ID0gKSkgKyBcbmBgYCJ9 -->

```r
#p7<- ggplot(data = output, aes(x = , y = )) +
#p8<- ggplot(data = output, aes(x = , y = )) + 
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




## Extra Models
Basides the models shown in the prevous section, there are many other models or variations of that models we could use to explore this problem. However, some of them are more time consuming, thus we trained some extra models and we saved their results, so now all we need is make the predictions and check their RMSE. 


First we define the feature and the model we want to explore:

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBvdmVyd3JpdGUgbW9kZWw/XG5vdmVyd3JpdGUgPSBGXG5cbiMgZmVhdHVyZVxuZmVhdHVyZV9uYW1lID0gXCJmZWF0dXJlX3RyYVwiICMgcmEsIGEsIHBjYSwgcmF3XG5cbiMjbG9hZCB0aGUgZmVhdHVyZVxuZmVhdHVyZSA9IHJlYWQuY3N2KHBhc3RlMChmZWF0X2RpcixcIi9cIiwgZmVhdHVyZV9uYW1lLCcuY3N2JykpXG5yb3duYW1lcyhmZWF0dXJlKSA9IGZlYXR1cmVbLDFdXG5mZWF0dXJlID0gYXMubWF0cml4KGZlYXR1cmVbLC0xXSlcblxuI3NwbGl0aW5nIHRoZSB2YWxpZGF0aW9uIGFuZCB0ZXN0aW5nIHNldFxuZmVhdHVyZV90ciA9IGZlYXR1cmVbdHJhaW5faW5kZXhfdHIsXVxuZmVhdHVyZV92YWwgPSBmZWF0dXJlW3RyYWluX2luZGV4X3ZhbCxdXG5gYGAifQ== -->

```r
# overwrite model?
overwrite = F

# feature
feature_name = "feature_tra" # ra, a, pca, raw

##load the feature
feature = read.csv(paste0(feat_dir,"/", feature_name,'.csv'))
rownames(feature) = feature[,1]
feature = as.matrix(feature[,-1])

#spliting the validation and testing set
feature_tr = feature[train_index_tr,]
feature_val = feature[train_index_val,]
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Then we run this chunk that will load the model case the model is available or run the model with the feature and model specified.

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBtb2RlbCBvcHRpb25zOlxuI2VuZXQsIGZvYmEsIGdhdXNzcHJQb2x5ZywgZ2F1c3NwclJhZGlhbCxnbG1uZXQsIGljciwga2VybmVscGxzLCBrcmxzUmFkaWFsLCBsYXJzMiwgbGFzc28sIGxlYXBCYWNrd2FyZCwgbGVhcEZvcndhcmQsIGxlYXBTZXEsIG5ubHMsIHBhcnREU0EsIHBscyAsIHBsc1JnbG0sIHJiZiwgcnBhcnQsIHJxbGFzc28sIHJ2bVBvbHksIHJ2bVJhZGlhbCwgc2ltcGxzLCBzcGlrZXNsYWIsIHNwbHMsIHN2bVBvbHksIHN2bVJhZGlhbCwgc3ZtUmFkaWFsQ29zdCwgc3ZtUmFkaWFsU2lnbWEsIHdpZGVrZXJuZWxwbHMsIHJxbmMsIG5vZGVIYXJ2ZXN0LCBtbHBNTCwgeGdiREFSVFxuI3NlZSBodHRwczovL3RvcGVwby5naXRodWIuaW8vY2FyZXQvYXZhaWxhYmxlLW1vZGVscy5odG1sIGZvciBtb3JlIGluZm9ybWF0aW9uIGFib3V0IGVhY2ggbW9kZWxcblxubW9kZWwgPSBcImVuZXRcIlxucGFyc2kgPSBleHBhbmQuZ3JpZChsYW1iZGE9MTBecnVuaWYoNSwgbWluPS01LCAxKSwgZnJhY3Rpb249cnVuaWYoNSwgbWluPTAsIG1heD0xKSkgIyBwYXJzaSA9IE5VTEw7IGlmIG5vIHBhcmFtZXRlcnMgbmVlZCB0byBiZSB0ZXN0ZWRcbmZpdGN2ID0gdHJhaW5Db250cm9sKG1ldGhvZD1cImN2XCIsIG51bWJlcj1jdm4pXG5cblxuZm5hbWUgPSBwYXN0ZTAobW9kZWxfZGlyLFwiL1wiLGZlYXR1cmVfbmFtZSxcIi9cIixtb2RlbCxcIi5SZGF0YVwiKVxuaWYgKCFmaWxlLmV4aXN0cyhmbmFtZSkgfCBvdmVyd3JpdGUpIHsgdHJ5ICh7XG4gIHQyaSA9IE5VTExcbiAgaWYgKCFpcy5udWxsKHBhcnNpKSkge1xuICAgIHQyaSA9IGNhcmV0Ojp0cmFpbih5PWdhX3RyLCB4PWZlYXR1cmVfdHIsIG1vZGVsLCB0ckNvbnRyb2w9Zml0Y3YsIHR1bmVHcmlkPXBhcnNpKVxuICB9IGVsc2Uge1xuICAgIHQyaSA9IGNhcmV0Ojp0cmFpbih5PWdhX3RyLCB4PWZlYXR1cmVfdHIsIG1vZGVsLCB0ckNvbnRyb2w9Zml0Y3YpXG4gIH1cbiAgaWYgKCFpcy5udWxsKHQyaSkpIHNhdmUodDJpLCBmaWxlPWZuYW1lKVxufSkgfVxubG9hZChmbmFtZSlcblxuIyByZXN1bHRzIHRvIGRhdGEgZnJhbWVcbmRmID0gZGF0YS5mcmFtZShcbiAgcm1zZT10MmkkcmVzdWx0cyRSTVNFW3doaWNoLm1pbih0MmkkcmVzdWx0cyRSTVNFKV0sXG4gIHRpbWU9YXMubnVtZXJpYyh0MmkkdGltZXMkZXZlcnl0aGluZ1szXSksXG4gIG1vZGVsXz10MmkkbW9kZWxJbmZvJGxhYmVsLCBcbiAgZmVhdHVyZT1mZWF0dXJlX25hbWUsIG1vZGVsPW1vZGVsLCBcbiAgcGFyPXBhc3RlMCggcGFzdGUwKG5hbWVzKHQyaSRiZXN0VHVuZSksIGNvbGxhcHNlPVwiX1wiKSwgXCI6IFwiLCBwYXN0ZTAodDJpJGJlc3RUdW5lLCBjb2xsYXBzZT1cIl9cIikgKVxuICAsIHN0cmluZ3NBc0ZhY3RvcnM9RilcbmRmXG5cbmBgYCJ9 -->

```r
# model options:
#enet, foba, gaussprPolyg, gaussprRadial,glmnet, icr, kernelpls, krlsRadial, lars2, lasso, leapBackward, leapForward, leapSeq, nnls, partDSA, pls , plsRglm, rbf, rpart, rqlasso, rvmPoly, rvmRadial, simpls, spikeslab, spls, svmPoly, svmRadial, svmRadialCost, svmRadialSigma, widekernelpls, rqnc, nodeHarvest, mlpML, xgbDART
#see https://topepo.github.io/caret/available-models.html for more information about each model

model = "enet"
parsi = expand.grid(lambda=10^runif(5, min=-5, 1), fraction=runif(5, min=0, max=1)) # parsi = NULL; if no parameters need to be tested
fitcv = trainControl(method="cv", number=cvn)


fname = paste0(model_dir,"/",feature_name,"/",model,".Rdata")
if (!file.exists(fname) | overwrite) { try ({
  t2i = NULL
  if (!is.null(parsi)) {
    t2i = caret::train(y=ga_tr, x=feature_tr, model, trControl=fitcv, tuneGrid=parsi)
  } else {
    t2i = caret::train(y=ga_tr, x=feature_tr, model, trControl=fitcv)
  }
  if (!is.null(t2i)) save(t2i, file=fname)
}) }
load(fname)

# results to data frame
df = data.frame(
  rmse=t2i$results$RMSE[which.min(t2i$results$RMSE)],
  time=as.numeric(t2i$times$everything[3]),
  model_=t2i$modelInfo$label, 
  feature=feature_name, model=model, 
  par=paste0( paste0(names(t2i$bestTune), collapse="_"), ": ", paste0(t2i$bestTune, collapse="_") )
  , stringsAsFactors=F)
df

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

Now we will explore the predictions of these models+feature:


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



### Exercise 9: 
Explore more combinations of models and feature. Keep the results of the top 3 models you found. Discuss with your colleagues their results. Answer: 
* Are they similar? How you defined which models were the best? 
* Can you think in a more efficient way to compare all the models available? 


####Model 1 

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBmZWF0dXJlXG5mZWF0dXJlX25hbWUgPSBcImZlYXR1cmVfdHJhXCIgIyByYSwgYSwgcGNhLCByYXdcblxuIyNsb2FkIHRoZSBmZWF0dXJlXG5mZWF0dXJlID0gcmVhZC5jc3YocGFzdGUwKGZlYXRfZGlyLFwiL1wiLCBmZWF0dXJlX25hbWUsJy5jc3YnKSlcbnJvd25hbWVzKGZlYXR1cmUpID0gZmVhdHVyZVssMV1cbmZlYXR1cmUgPSBhcy5tYXRyaXgoZmVhdHVyZVssLTFdKVxuXG4jc3BsaXRpbmcgdGhlIHZhbGlkYXRpb24gYW5kIHRlc3Rpbmcgc2V0XG5mZWF0dXJlX3RyID0gZmVhdHVyZVt0cmFpbl9pbmRleF90cixdXG5mZWF0dXJlX3ZhbCA9IGZlYXR1cmVbdHJhaW5faW5kZXhfdmFsLF1cblxuIyBtb2RlbCBvcHRpb25zOlxuI2VuZXQsIGZvYmEsIGdhdXNzcHJQb2x5ZywgZ2F1c3NwclJhZGlhbCxnbG1uZXQsIGljciwga2VybmVscGxzLCBrcmxzUmFkaWFsLCBsYXJzMiwgbGFzc28sIGxlYXBCYWNrd2FyZCwgbGVhcEZvcndhcmQsIGxlYXBTZXEsIG5ubHMsIHBhcnREU0EsIHBscyAsIHBsc1JnbG0sIHJiZiwgcnBhcnQsIHJxbGFzc28sIHJ2bVBvbHksIHJ2bVJhZGlhbCwgc2ltcGxzLCBzcGlrZXNsYWIsIHNwbHMsIHN2bVBvbHksIHN2bVJhZGlhbCwgc3ZtUmFkaWFsQ29zdCwgc3ZtUmFkaWFsU2lnbWEsIHdpZGVrZXJuZWxwbHMsIHJxbmMsIG5vZGVIYXJ2ZXN0LCBtbHBNTCwgeGdiREFSVFxuI3NlZSBodHRwczovL3RvcGVwby5naXRodWIuaW8vY2FyZXQvYXZhaWxhYmxlLW1vZGVscy5odG1sIGZvciBtb3JlIGluZm9ybWF0aW9uIGFib3V0IGVhY2ggbW9kZWxcblxubW9kZWwgPSBcImVuZXRcIlxucGFyc2kgPSBleHBhbmQuZ3JpZChsYW1iZGE9MTBecnVuaWYoNSwgbWluPS01LCAxKSwgZnJhY3Rpb249cnVuaWYoNSwgbWluPTAsIG1heD0xKSkgIyBwYXJzaSA9IE5VTEw7IGlmIG5vIHBhcmFtZXRlcnMgbmVlZCB0byBiZSB0ZXN0ZWRcbmZpdGN2ID0gdHJhaW5Db250cm9sKG1ldGhvZD1cImN2XCIsIG51bWJlcj1jdm4pXG5cbmZuYW1lID0gcGFzdGUwKG1vZGVsX2RpcixcIi9cIixmZWF0dXJlX25hbWUsXCIvXCIsbW9kZWwsXCIuUmRhdGFcIilcbmlmICghZmlsZS5leGlzdHMoZm5hbWUpIHwgb3ZlcndyaXRlKSB7IHRyeSAoe1xuICB0MmkgPSBOVUxMXG4gIGlmICghaXMubnVsbChwYXJzaSkpIHtcbiAgICB0MmkgPSBjYXJldDo6dHJhaW4oeT1nYV90ciwgeD1mZWF0dXJlX3RyLCBtb2RlbCwgdHJDb250cm9sPWZpdGN2LCB0dW5lR3JpZD1wYXJzaSlcbiAgfSBlbHNlIHtcbiAgICB0MmkgPSBjYXJldDo6dHJhaW4oeT1nYV90ciwgeD1mZWF0dXJlX3RyLCBtb2RlbCwgdHJDb250cm9sPWZpdGN2KVxuICB9XG4gIGlmICghaXMubnVsbCh0MmkpKSBzYXZlKHQyaSwgZmlsZT1mbmFtZSlcbn0pIH1cbmxvYWQoZm5hbWUpXG5cbiMgcmVzdWx0cyB0byBkYXRhIGZyYW1lXG5kZiA9IGRhdGEuZnJhbWUoXG4gIHJtc2U9dDJpJHJlc3VsdHMkUk1TRVt3aGljaC5taW4odDJpJHJlc3VsdHMkUk1TRSldLFxuICB0aW1lPWFzLm51bWVyaWModDJpJHRpbWVzJGV2ZXJ5dGhpbmdbM10pLFxuICBtb2RlbF89dDJpJG1vZGVsSW5mbyRsYWJlbCwgXG4gIGZlYXR1cmU9ZmVhdHVyZV9uYW1lLCBtb2RlbD1tb2RlbCwgXG4gIHBhcj1wYXN0ZTAoIHBhc3RlMChuYW1lcyh0MmkkYmVzdFR1bmUpLCBjb2xsYXBzZT1cIl9cIiksIFwiOiBcIiwgcGFzdGUwKHQyaSRiZXN0VHVuZSwgY29sbGFwc2U9XCJfXCIpIClcbiAgLCBzdHJpbmdzQXNGYWN0b3JzPUYpXG5kZlxuXG4jI2dldCB0ZXN0IHByZWRpY3Rpb24gcmVzdWx0cyBmcm9tIG1vZGVscyAtLS0tLS0tLS0tLS0tLS0tLS0tLS0tXG5wcmVkID0gcHJlZGljdCh0MmksbmV3ZGF0YSA9IGFzLmRhdGEuZnJhbWUoZmVhdHVyZV90cikpXG5wcmVkX3ZhbCA9IHByZWRpY3QodDJpLG5ld2RhdGEgPSBhcy5kYXRhLmZyYW1lKGZlYXR1cmVfdmFsKSlcbnRpdGxlID0gcGFzdGUoZmVhdHVyZV9uYW1lLG1vZGVsLCBzZXAgPSAnLScpXG5cbiMgcGxvdCBncmFwaCB0byBjb21wYXJlIG1vZGVsc1xucGxvdDEgPSBkYXRhLmZyYW1lKGdhX3RyLHByZWQpXG5wbG90MiA9IGRhdGEuZnJhbWUoZ2FfdmFsLCBwcmVkX3ZhbClcbnByZWQxIDwtIGdncGxvdChwbG90MSwgYWVzKHggPSBnYV90cix5PXByZWQpKSArIGdlb21fcG9pbnQoY29sb3I9J29yYW5nZScpK1xuICB5bGFiKCdQcmVkaWN0ZWQnKSt4bGFiKCdUcmFpbmluZyBzZXQnKStsYWJzKHRpdGxlPXRpdGxlKVxuXG5wcmVkMiA8LSBnZ3Bsb3QocGxvdDIsIGFlcyh4ID0gZ2FfdmFsLHk9cHJlZF92YWwpKSArIGdlb21fcG9pbnQoY29sb3I9J29yYW5nZScpK1xuICB5bGFiKCdQcmVkaWN0ZWQnKSt4bGFiKCdWYWxpZGF0aW9uIHNldCcpXG5cbmdyaWQxIDwtIGdyaWQuYXJyYW5nZShwcmVkMSwgcHJlZDIsIG5yb3cgPSAxKVxuZ3JpZDFcbmdnc2F2ZShmaWxlbmFtZSA9IHBhc3RlKHJlc3VsdF9kaXIsJy8nLCdwcmVkaWN0aW9ucy0nLHRpdGxlLCcucG5nJyxzZXA9JycpLGdyaWQxKVxuXG5gYGAifQ== -->

```r
# feature
feature_name = "feature_tra" # ra, a, pca, raw

##load the feature
feature = read.csv(paste0(feat_dir,"/", feature_name,'.csv'))
rownames(feature) = feature[,1]
feature = as.matrix(feature[,-1])

#spliting the validation and testing set
feature_tr = feature[train_index_tr,]
feature_val = feature[train_index_val,]

# model options:
#enet, foba, gaussprPolyg, gaussprRadial,glmnet, icr, kernelpls, krlsRadial, lars2, lasso, leapBackward, leapForward, leapSeq, nnls, partDSA, pls , plsRglm, rbf, rpart, rqlasso, rvmPoly, rvmRadial, simpls, spikeslab, spls, svmPoly, svmRadial, svmRadialCost, svmRadialSigma, widekernelpls, rqnc, nodeHarvest, mlpML, xgbDART
#see https://topepo.github.io/caret/available-models.html for more information about each model

model = "enet"
parsi = expand.grid(lambda=10^runif(5, min=-5, 1), fraction=runif(5, min=0, max=1)) # parsi = NULL; if no parameters need to be tested
fitcv = trainControl(method="cv", number=cvn)

fname = paste0(model_dir,"/",feature_name,"/",model,".Rdata")
if (!file.exists(fname) | overwrite) { try ({
  t2i = NULL
  if (!is.null(parsi)) {
    t2i = caret::train(y=ga_tr, x=feature_tr, model, trControl=fitcv, tuneGrid=parsi)
  } else {
    t2i = caret::train(y=ga_tr, x=feature_tr, model, trControl=fitcv)
  }
  if (!is.null(t2i)) save(t2i, file=fname)
}) }
load(fname)

# results to data frame
df = data.frame(
  rmse=t2i$results$RMSE[which.min(t2i$results$RMSE)],
  time=as.numeric(t2i$times$everything[3]),
  model_=t2i$modelInfo$label, 
  feature=feature_name, model=model, 
  par=paste0( paste0(names(t2i$bestTune), collapse="_"), ": ", paste0(t2i$bestTune, collapse="_") )
  , stringsAsFactors=F)
df

##get test prediction results from models ----------------------
pred = predict(t2i,newdata = as.data.frame(feature_tr))
pred_val = predict(t2i,newdata = as.data.frame(feature_val))
title = paste(feature_name,model, sep = '-')

# plot graph to compare models
plot1 = data.frame(ga_tr,pred)
plot2 = data.frame(ga_val, pred_val)
pred1 <- ggplot(plot1, aes(x = ga_tr,y=pred)) + geom_point(color='orange')+
  ylab('Predicted')+xlab('Training set')+labs(title=title)

pred2 <- ggplot(plot2, aes(x = ga_val,y=pred_val)) + geom_point(color='orange')+
  ylab('Predicted')+xlab('Validation set')

grid1 <- grid.arrange(pred1, pred2, nrow = 1)
grid1
ggsave(filename = paste(result_dir,'/','predictions-',title,'.png',sep=''),grid1)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


####model 2 

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


####model 3

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



#Prediction For the Competition 

As all Data SCience/Machine Learning project, we have spend most of the time so far procesing the data, exploring the dataset, creating/selecting feature, exploring different models to define which one is the best option for our main task: predict the Gestacional Age of our testing set. 

Now, we will finale create the predictions for our test set. 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI1JVTiBIRVJFIE9OQ0UgTU9SRSBUSEUgTU9ERUwgWU9VIFRISU5LIElTIFRIRSBCRVNUXG4jUkVNRU1CRVIgVE8gTE9BRCBUSEUgRkVBVFVSRVMgQU5EIFRSQUlOIFRIRSBNT0RFTCBVU0lORyBUSEUgVFJBSU5JTkcgU0VUXG5gYGAifQ== -->

```r
#RUN HERE ONCE MORE THE MODEL YOU THINK IS THE BEST
#REMEMBER TO LOAD THE FEATURES AND TRAIN THE MODEL USING THE TRAINING SET
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZmVhdHVyZV90cnByZWQgPSBmZWF0dXJlW3Rlc3RfaW5kZXgsXVxuXG4jSU5TRVJUIEhFUkUgTU9ERUxcbiNwcmVkaWN0aW9ucyA9IHByZWRpY3QoLU1PREVMLSwgZGF0YSA9IGFzLmRhdGEuZnJhbWUoZmVhdHVyZV90cnByZWQpKSBcblxuIyBzdWJtaXNzaW9uIHRlbXBsYXRlXG5jbGFzc19maW5hbCA9IHJlYWQuY3N2KHBhc3RlMChpbnB1dF9kaXIsXCIvVGVhbVhfU0MxX3ByZWRpY3Rpb24uY3N2XCIpKSBcbmNsYXNzX2ZpbmFsJEdBID0gcm91bmQocHJlZGljdGlvbnMsMSlcbndyaXRlLmNzdihjbGFzc19maW5hbCwgZmlsZT1wYXN0ZTAocmVzdWx0X2RpcixcIi9zdWJtaXNzaW9uX2ZpbGUuY3N2XCIpKVxuYGBgIn0= -->

```r
feature_trpred = feature[test_index,]

#INSERT HERE MODEL
#predictions = predict(-MODEL-, data = as.data.frame(feature_trpred)) 

# submission template
class_final = read.csv(paste0(input_dir,"/TeamX_SC1_prediction.csv")) 
class_final$GA = round(predictions,1)
write.csv(class_final, file=paste0(result_dir,"/submission_file.csv"))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->








<!-- rnb-text-end -->

