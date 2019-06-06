#Deleting previous files
rm(list=ls(all=T))

#Loading the dataset
load("Data/HTA20_RMA.RData")
dim(eset_HTA20)
  
  
#probeset is compoised of entrez-gene ID
#pid has the patient id 
probeset = rownames(eset_HTA20)
pid = colnames(eset_HTA20)
head(probeset,20)
head(pid)


#Transforming dataset: now rows are patients and columns are genes
data = t(eset_HTA20)
head(eset_HTA20[,c(1:30)])
rm(eset_HTA20)

#Reading sample annotation file 
#SampleID: unique identifier of the sample (matching the name of the .CEL file in HTA20 folder, 
#except for extension .CEL);
#GA: gestational age as determined by the last menstrual period and or ultrasound;
#Batch: the batch identifier;
#Set: name of the source dataset;
#Train: 1 for samples to be used for training, 0 for samples to be used for test;
#Platform: gene expression platform used to generate the cell files.
sample = read.csv('Data/anoSC1_v11_nokey.csv', sep = ',', header=T)


#Feature Extraction 
#1) Removing genes with many missing values
#There is no missing values
#na_count <-sapply(data, function(y) sum(length(which(is.na(y)))))
#sum(na_count)

#2) Removing genes with lowest variance 
require(Rfast)
#var_gene = colVars(data)
eliminate = data.frame(col = c(1:dim(data)[2]),var_gene = colVars(data))
head(eliminate)
eliminate = subset(eliminate, var_gene<quantile(eliminate$var_gene,0.1))

da
#3) PCA



#4) Random Forest 
#5) Autoencoder 


