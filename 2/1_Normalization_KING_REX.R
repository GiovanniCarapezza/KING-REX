library(DESeq2)
library(stringi)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2/"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/2",sep=""))
#path to raw counts
path_input_data=paste(sep="",path_to_results,"/Results/raw_counts_KING-REX/")
#path reformat rawcount for deseq
path=paste("raw_deseq_KING-REX",sep="")
if(!dir.exists(path)){dir.create(path,recursive = T)}
path_output=paste("normalized_counts_KING-REX",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}


#read kinase list
kin_list=read.table(paste(path_to_results,"/Results/1/kinase_list_and_info.txt",sep=""),sep="\t",stringsAsFactors=F,row.names = 1,header = T)
kin_list=kin_list[sort(rownames(kin_list)),]

#read annotation file
input=read.table(paste(path_to_results,"/Results/annotations_KING-REX.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)

#reformat for deseq
for(i in 1:dim(input)[1]){
  a=read.table(file=paste(path_input_data,"/",input[i,1],".genes.results",sep=""),sep="\t",row.names=1,header = T,stringsAsFactors = F)
  change_names=which(rownames(a)=="SGK110")
  rownames(a)[change_names]="SBK3"
  change_names=which(rownames(a)=="SGK196")
  rownames(a)[change_names]="POMK"
  change_names=which(rownames(a)=="KIAA1804")
  rownames(a)[change_names]="MAP3K21"
  a2=a[rownames(kin_list),]
  write.table(file=paste(path,"/",input[i,1],".txt",sep=""),
              x=cbind(rownames(a2),round(a2$expected_count,0)),sep="\t",row.names=F,col.names=F,quote=F,append = F)
}


#create sample table for deseq
sampleName=c()
sampleFiles=c()
sampleCondition=c()
sampleCondition2=c()
names=c()
for(i in 1:dim(input)[1]){
  name=input$Sample[i]
  sampleName=c(sampleName,name)
  sampleFiles=c(sampleFiles,paste(name,".txt",sep=""))
  sampleCondition=c(sampleCondition,input[i,"Sample_Name"])
  names=c(names,input$Sample_Name[i])
}
sampleTable<-data.frame(sampleName=sampleName, fileName=sampleFiles, condition=sampleCondition)

#create deseq object
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=path, design=~condition)
#Estimate the size factors
esf=estimateSizeFactors(ddsHTSeq)
sf=sizeFactors(esf)
#extract normalized data and log2
count_norm=counts(esf,normalize=T)
count_norm_log2=log2(count_norm+1)

#write normalized counts
write.table(count_norm,paste(path_output,"/","normalized_counts.txt",sep=""),quote=FALSE,sep="\t")
write.table(count_norm_log2,paste(path_output,"/","normalized_counts_log2.txt",sep=""),quote=FALSE,sep="\t")
#

#plots
##################################
distscn <- dist(t(count_norm_log2))
mat <- as.matrix(distscn)
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
annotation<-as.data.frame(cbind(sampleCondition),row.names=sampleName)
pheatmap(mat,fontsize = 24,
         clustering_distance_rows=distscn,
         annotation_col=annotation,
         clustering_distance_cols=distscn,
         col = rev(hmcol))


melt_norm=melt(count_norm_log2)
colnames(melt_norm)[3]="Log2 Norm Count"
colnames(melt_norm)[2]="Samples"
p=ggplot(melt_norm, 
         aes(x=Samples, y=`Log2 Norm Count`)) + coord_cartesian(ylim=c(0,25)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),text = element_text(size=24)) + geom_boxplot()
print(p)

unlink(path,recursive = T)






