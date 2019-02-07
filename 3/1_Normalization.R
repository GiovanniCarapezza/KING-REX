library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(stringi)

rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/3",sep=""))
#path to raw counts
path_input_data=paste(sep="",path_to_results,"/Results/raw_counts_KING-REX/")
#path reformat rawcount for deseq
path=paste("raw_counts_IN_OUT",sep="")
if(dir.exists(path)){
  unlink(path,recursive = T)
  dir.create(path,recursive = T)
}else{
  dir.create(path,recursive = T)
}
#path output
path_output=paste("normalized_counts_IN_OUT",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}
path_output2=paste("edgeR_input",sep="")
if(!dir.exists(path_output2)){dir.create(path_output2,recursive = T)}


#read kinase list
kin_list=read.table(paste(path_to_results,"/Results/1/kinase_list_and_info.txt",sep=""),sep="\t",stringsAsFactors=F,row.names = 1,header = T)
kin_list=kin_list[sort(rownames(kin_list)),]

#read annotation file
input=read.table(paste(path_to_results,"/Results/annotations_KING-REX.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)

#reformat for deseq and select only the 319 kinase for differential expression analysis
for(i in 1:dim(input)[1]) {
  a=read.table(file=paste(path_input_data,"/",input[i,1],"_trex_intervals_new.cov",sep=""),sep="\t",row.names=1)
  for (j in 1:dim(kin_list)[1]) {
    kin=rownames(kin_list)[j]
    kin_in_v=a[paste(kin,"_IN",sep=""),"V2"]
    kin_out_v=a[paste(kin,"_OUT",sep=""),"V2"]
    
    if(kin_list$DIFF[j]=="yes"){
        write.table(file=paste(path,"/",input[i,1],"_IN.txt",sep=""),x=cbind(kin,kin_in_v),sep="\t",row.names=F,col.names=F,quote=F,append = T)
        write.table(file=paste(path,"/",input[i,1],"_OUT.txt",sep=""),x=cbind(kin,kin_out_v),sep="\t",row.names=F,col.names=F,quote=F,append = T)
    }
  }
}

#create sample table for deseq
sampleName=c()
sampleFiles=c()
sampleCondition=c()
sampleCondition2=c()
names=c()
for(i in 1:dim(input)[1]){
  name=input$Sample[i]
  sampleName=c(sampleName,paste(name,"_IN",sep=""),paste(name,"_OUT",sep=""))
  sampleFiles=c(sampleFiles,paste(name,"_IN.txt",sep=""),paste(name,"_OUT.txt",sep=""))
  sampleCondition=c(sampleCondition,"IN","OUT")
  names=c(names,input$Sample_Name[i],input$Sample_Name[i])
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

############Further step of normailzation######################
#sort
count_norm_log2=count_norm_log2[rownames(kin_list)[which(kin_list$DIFF=="yes")],]
count_norm_log2_old=count_norm_log2
#calculate scaling factor
eff_med=c()
for(kin_a in rownames(count_norm_log2_old)){
  eff_medd_app=c()
  #discard dilution samples and low quality
  for(sample in input$Sample[c(1:24,35:42,45:48)]){
    val1=count_norm_log2_old[kin_a,paste(sample,"_IN",sep="")]
    val2=count_norm_log2_old[kin_a,paste(sample,"_OUT",sep="")]
    eff_medd_app=c(eff_medd_app,(val1-val2))
  }
  
  eff_med=c(eff_med,median(eff_medd_app))
}
#apply scaling factor
count_norm_log2[,stri_detect_fixed(str = colnames(count_norm_log2),pattern = "_OUT")]=
  count_norm_log2[,stri_detect_fixed(str = colnames(count_norm_log2),pattern = "_OUT")]+
  eff_med
count_norm_log2[count_norm_log2<0]=0


#write normalized data only for folder 3 data
write.table(count_norm_log2_old,paste(path_output,"/","prenormalized_IN_OUT_counts_log2.txt",sep=""),quote=FALSE,sep="\t")
write.table(count_norm_log2,paste(path_output,"/","normalized_IN_OUT_counts_log2.txt",sep=""),quote=FALSE,sep="\t")

#prepare edgeR file
for(i in 1:dim(count_norm_log2)[2]) {
  counts_col=round(as.numeric((2^count_norm_log2[,i])-1),0)
  counts_col[counts_col<0]=0
  write.table(file=paste(path_output2,"/",colnames(count_norm_log2)[i],".txt",sep=""),x=cbind("Gene","counts"),sep="\t",row.names=F,col.names=F,quote=F,append = F)  
  write.table(file=paste(path_output2,"/",colnames(count_norm_log2)[i],".txt",sep=""),x=cbind(rownames(count_norm_log2),counts_col),sep="\t",row.names=F,col.names=F,quote=F,append = T)  
}
write.table(cbind(names,sampleName,sampleFiles,sampleCondition),file =paste(path_output2,"/", "condition",sep=""),quote=FALSE,sep="\t",row.names = F)


#distance matrix pre and post normalization
annotation=as.data.frame(cbind(names,sampleName,sampleFiles,sampleCondition),stringsAsFactors=F)
annotation=annotation[annotation$names %in% c("KARPAS299","KM12","LC2AD","U118MG","NCIH716"),]
annotation[,"Tissue"]=""
for(row in 1:nrow(annotation)){
  annotation$Tissue[row]=input$Tissue[which(input$Sample_Name==annotation$names[row])[1]]
}

trex_IO=count_norm_log2_old[,annotation$sampleName]
distscn <- dist(t(trex_IO))
mat <- as.matrix(distscn)
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
annotation2<-as.data.frame(cbind(annotation$sampleCondition),row.names=annotation$sampleName)
ann_colors=list(
  V1=c(IN="green",OUT="orange")
)
wd = 35*nrow(annotation2)
jpeg(paste(path_output,"Fig2A_distance_matrix_prenormalized.jpeg",sep="/"),width = wd+400,height = wd+100)
pheatmap(mat,fontsize = 24,cellwidth = 20,cellheight = 20,
         clustering_distance_rows=distscn,
         annotation_col=annotation2,
         annotation_colors = ann_colors,
         clustering_distance_cols=distscn,
         col = rev(hmcol))
dev.off()


trex_IO=count_norm_log2[,annotation$sampleName]
distscn <- dist(t(trex_IO))
mat <- as.matrix(distscn)
jpeg(paste(path_output,"Fig2B_distance_matrix_normalized.jpeg",sep="/"),width = wd+400,height = wd+100)
pheatmap(mat,fontsize = 24,cellwidth = 20,cellheight = 20,
         clustering_distance_rows=distscn,
         annotation_col=annotation2,
         annotation_colors = ann_colors,
         clustering_distance_cols=distscn,
         col = rev(hmcol))
dev.off()

unlink(path,recursive = T)
