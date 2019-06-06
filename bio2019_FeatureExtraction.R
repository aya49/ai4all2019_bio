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
require(Rfast) #Variance by column s
require(ggplot2) #Plots


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

#Feature Extraction 
#1) Removing genes with many missing values
#There is no missing values
#na_count <-sapply(data, function(y) sum(length(which(is.na(y)))))
#sum(na_count)

#2) Removing genes with lowest variance 


#var_gene = colVars(data)
eliminate = data.frame(col = c(1:dim(data)[2]),var_gene = colVars(data))
head(eliminate)
eliminate = subset(eliminate, var_gene<quantile(eliminate$var_gene,0.1))
dim(data)
data = subset(data, select = -c(eliminate$col))
dim(data)

#3) PCA
genes.pca <- prcomp(data, center = TRUE,scale. = TRUE)
#https://www.datacamp.com/community/tutorials/pca-analysis-r
vars <- colVars(genes.pca$x)  
props <- cumsum(vars / sum(vars))


plot(props, xlab = 'Components (1-735')


#4) Random Forest 
#5) Autoencoder 


