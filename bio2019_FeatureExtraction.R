#--------#--------#--------#--------#--------#--------#--------#--------#--------#
#--------#--------#--------#--------#--------#--------#--------#--------#--------#
# Author: Raquel Aoki
# Data: 2019/06 
#--------#--------#--------#--------#--------#--------#--------#--------#--------#
#--------#--------#--------#--------#--------#--------#--------#--------#--------#


#--------#--------#--------#--------#--------#--------#--------#--------#--------#
# I: Preprossing the data 
#--------#--------#--------#--------#--------#--------#--------#--------#--------#

#-------- Cleaning the environment 
rm(list=ls(all=T))


#-------- Packages 
if (!require("Rfast")) install.packages("Rfast")
if (!require("ggplot2")) install.packages("ggplot2")
#if (!require("randomForest")) install.packages("randomForest")
if (!require("caret")) install.packages("caret")
if (!require("e1071")) install.packages("e1071")
if (!require("ranger")) install.packages("ranger")
if (!require("dplyr")) install.packages("dplyr")
if (!require("deepnet")) install.packages("deepnet")


require(Rfast) #Variance by column s
require(ggplot2) #Plots
require(caret) #Random Forest 
#require(randomForest)
require(e1071)
require(ranger)
require(dplyr)
require(deepnet)


#-------- Loading the dataset
#RNASEQ data: each row is a gene and each column a patient 
#probeset: gene ID
#pid: patient id 
load("Data/HTA20_RMA.RData")

#Reading sample annotation file 
#SampleID: unique identifier of the sample (matching the name of the .CEL file in HTA20 folder, 
#except for extension .CEL);
#GA: gestational age as determined by the last menstrual period and or ultrasound;
#Batch: the batch identifier;
#Set: name of the source dataset;
#Train: 1 for samples to be used for training, 0 for samples to be used for test;
#Platform: gene expression platform used to generate the cell files.
sample = read.csv('Data/anoSC1_v11_nokey.csv', sep = ',', header=T)


#--------  Checking datasets dimentions 
dim(eset_HTA20)
dim(sample)
  

#-------- Saving the genes IDS and Patient IDS  
gid = rownames(eset_HTA20)
pid = colnames(eset_HTA20)
#Checking the ids 
head(gid)
head(pid)


#-------- Transforming gene dataset: now rows are patients and columns are genes
data = t(eset_HTA20)
rm(eset_HTA20)


#--------#--------#--------#--------#--------#--------#--------#--------#--------#
# II: Feature Extraction  
#--------#--------#--------#--------#--------#--------#--------#--------#--------#

#-------- 1) Removing genes with lowest variance 
eliminate = data.frame(col = c(1:dim(data)[2]),var_gene = colVars(data))
head(eliminate)
eliminate = subset(eliminate, var_gene<quantile(eliminate$var_gene,0.3))
dim(data)
data = subset(data, select = -c(eliminate$col))
dim(data)

#--------2) PCA
#https://www.datacamp.com/community/tutorials/pca-analysis-r
#Don't need the values to be predicted, only the features

#Option 1: slower 
#genes.pca <- prcomp(data, center = TRUE,scale. = TRUE)
#vars <- colVars(genes.pca$x)  
#props <- data.frame(pc = 1:length(vars),cumprob = cumsum(vars / sum(vars)))
#pca <- train(Class~., data=dataset, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
#set.seed(12)
#ggplot(data = props, mapping = aes(x = props$pc,y=props$cumprob))+
#  geom_point()+xlab('PC')+ylab('Cumulative Proportion')+ 
#  theme(panel.background = element_rect(fill = "slategray2"))

#Option 2: faster
PrePCA <- preProcess(data,method="pca")
feat.pca <- predict(PrePCA,data)
PrePCA

write.csv(feat.pca, 'Data/features_pca.csv', row.names = T)
#--------3) Removing elements with low correlation with GA
data1 = data.frame(SampleID=row.names(data), data)
row.names(data1) = NULL
#Combining the two datasets 
data1 = merge(sample[,c(1,2,5)],data1,by.x = 'SampleID',by.y = 'SampleID', all = T)
#Using only the train dataset to make the feature importance 
data1 = subset(data1, Train == 1)
data1 = subset(data1, select = -c(Train,SampleID))

eliminate = data.frame(col = c(1:dim(data)[2]),corr = apply(data1[,-1], 2, cor, x =data1$GA))
eliminate$corr = abs(eliminate$corr)
eliminate = subset(eliminate, corr<quantile(eliminate$corr,0.3))
dim(data1)
data1 = subset(data1, select = -c(eliminate$col))
dim(data1)



#--------4) Random Forest 
#https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest
#https://uc-r.github.io/random_forests
#http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
#It depends on the features and the value to be predicted 
#data1 = as.matrix(data1)
metric <- "Accuracy"
set.seed(123)
#Number randomely variable selected is mtry
PreRF <- train(GA~., data=data1, method='ranger')

print(PreRF)



#5) Autoencoder
#https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/
#  https://www.rdocumentation.org/packages/ANN2/versions/1.5/topics/autoencoder
PreA <- train(GA~.,  data=data1, method = 'dnn',hidden = c(1),
              activation = 'linear', output = 'linear',
              numepochs = 3, batchsize = 100)
print(PreA)


