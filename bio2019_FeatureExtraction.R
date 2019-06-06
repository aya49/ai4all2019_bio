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
if (!require("randomForest")) install.packages("randomForest")

require(Rfast) #Variance by column s
require(ggplot2) #Plots
require(randomForest) #Random Forest 

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

#-------- 1) Removing genes with many missing values
#There is no missing values
#na_count <-sapply(data, function(y) sum(length(which(is.na(y)))))
#sum(na_count)

#-------- 2) Removing genes with lowest variance 
eliminate = data.frame(col = c(1:dim(data)[2]),var_gene = colVars(data))
head(eliminate)
eliminate = subset(eliminate, var_gene<quantile(eliminate$var_gene,0.1))
dim(data)
data = subset(data, select = -c(eliminate$col))
dim(data)

#--------3) PCA
#https://www.datacamp.com/community/tutorials/pca-analysis-r
#Don't need the values to be predicted, only the features

genes.pca <- prcomp(data, center = TRUE,scale. = TRUE)
vars <- colVars(genes.pca$x)  
props <- data.frame(pc = 1:length(vars),cumprob = cumsum(vars / sum(vars)))


ggplot(data = props, mapping = aes(x = props$pc,y=props$cumprob))+
  geom_point()+xlab('PC')+ylab('Cumulative Proportion')+ 
  theme(panel.background = element_rect(fill = "slategray2"))

write.csv(genes.pca$x[,c(1:500)], 'Data/features_pca.csv', row.names = T)

#--------4) Random Forest 
#https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest
#https://uc-r.github.io/random_forests
#It depends on the features and the value to be predicted 
data1 = data.frame(SampleID=row.names(data), data)
row.names(data1) = NULL
#Combining the two datasets 
data1 = merge(data1, sample[,c(1,2,5)],by.x = 'SampleID',by.y = 'SampleID', all = T)
#Using only the train dataset to make the feature importance 
data1 = subset(data1, Train == 1)
data1 = subset(data1, select = -c(Train,SampleID))
data1 = as.matrix(data1)
genes.rf = randomForest(formula = GA~., data = data1, ntrees = 50, importance = T )



5) Autoencoder
#https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/
#  https://www.rdocumentation.org/packages/ANN2/versions/1.5/topics/autoencoder


