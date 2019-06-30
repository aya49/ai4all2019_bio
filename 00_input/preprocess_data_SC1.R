#preprocess HTA2.0 at gene level

library(pd.hta20.hs.entrezg)
library(oligo)
anoSC1=read.csv("anoSC1_v11_nokey.csv")

pile=list()
blist=list(s1=c(1,2,3,4),s2=c(5,6,7,8),s3=c(9,10,32))
for(myset in names(blist)){
ano=anoSC1[anoSC1$Batch%in%blist[[myset]],]

list = list.celfiles("./CELL/",full.names=TRUE)[dir("./CELL/")%in%paste(ano$SampleID,".CEL",sep="")]
data = read.celfiles(list,pkgname="pd.hta20.hs.entrezg")
# = read.celfiles(list)

ppData <- rma(data) 
eset=exprs(ppData) 
eset=data.frame(eset)
colnames(eset)=gsub(".CEL","",colnames(eset))
pile[[myset]]<-eset

}
eset_HTA20=do.call(cbind,pile)
colnames(eset_HTA20)<-unlist(lapply(strsplit(colnames(eset_HTA20),split="\\."),"[[",2))
batch=anoSC1$Batch[match(colnames(eset_HTA20),anoSC1$SampleID)]

library(preprocessCore)
library(limma)
x=normalize.quantiles(as.matrix(eset_HTA20))
rownames(x)<-rownames(eset_HTA20)
colnames(x)<-colnames(eset_HTA20)
eset_HTA20 =  removeBatchEffect(eset_HTA20,batch=factor(batch))

save(eset_HTA20,file="HTA20_RMA.RData")


#preprocess HTA2.0 at probeset level

rm(list=ls())

library("hta20stprobeset.db")
library(oligo)
anoSC1=read.csv("anoSC1_v11_nokey.csv")

pile=list()
blist=list(s1=c(1,2),s2=c(3,4),s3=c(5,6),s4=c(7,8),s5=c(9,10,32))
for(myset in names(blist)){
  ano=anoSC1[anoSC1$Batch%in%blist[[myset]],]
  
  list = list.celfiles("./CELL/",full.names=TRUE)[dir("./CELL/")%in%paste(ano$SampleID,".CEL",sep="")]
  data = read.celfiles(list)
  # = read.celfiles(list)
  
  ppData <- rma(data,target="probeset") 
  eset=exprs(ppData) 
  eset=data.frame(eset)
  colnames(eset)=gsub(".CEL","",colnames(eset))
  pile[[myset]]<-eset
  
}
eset_HTA20=do.call(cbind,pile)

colnames(eset_HTA20)<-unlist(lapply(strsplit(colnames(eset_HTA20),split="\\."),"[[",2))
batch=anoSC1$Batch[match(colnames(eset_HTA20),anoSC1$SampleID)]


library(preprocessCore)
library(limma)
x=normalize.quantiles(as.matrix(eset_HTA20))
rownames(x)<-rownames(eset_HTA20)
colnames(x)<-colnames(eset_HTA20)
eset_HTA20_probeset =  removeBatchEffect(eset_HTA20,batch=factor(batch))

save(eset_HTA20_probeset,file="HTA20_RMA_probeset.RData")
#

