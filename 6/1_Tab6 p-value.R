#to excecute this script, you must run first script in folder 3
library(edgeR)
library(stringi)
library(pheatmap)
rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/6",sep=""))

path_input=paste(path_to_results,"/Results/3/edgeR_input",sep="")
path_output=paste("Imbalance analysis",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}

#sample to test ordered from T=0 to T=60 minutes
samples_to_test=c("KM12T0", "KM12T10", "KM12T20", "KM12T60")

#output file
filter_pval=1E-70
filter_FC=2
out_fusion=paste("fusion","_pv_",filter_pval,"_fc_",filter_FC,".txt",sep="")
write(x = c("Sample","Kinase","FC","p-value","IN_A","IN_B","OUT_A","OUT_B"),file = paste(path_output,out_fusion,sep="/"),ncolumns = 8,append = F,sep = "\t")

#detect imbalanced kinases
for(sample in samples_to_test){
  
  #create sample table for edgeR 
  sampleTable=read.table(paste(sep="",path_input,"/condition"),header=T,stringsAsFactors = F)
  samples=unique(sampleTable[,1])
  dir_sample=paste(path_output,"/",sample,sep="")
  if(!dir.exists(dir_sample)){dir.create(dir_sample,recursive = T)}
  sel=sampleTable[,1]==sample
  sampleTable2=sampleTable[sel,c(2,3,4)]
  sampleTable2[,3]=as.factor(sampleTable2[,3])
  rownames(sampleTable2)=c(1:dim(sampleTable2)[1])
  colnames(sampleTable2)=c("sampleName","fileName","condition")
  targets=sampleTable2
  targets$fileName <- paste(path_input, targets$fileName, sep='/')
  colnames(targets)=c("labels","files","description")
  targets$group=as.character(targets$description)
  targets$group[targets$group=="IN"]=1
  targets$group[targets$group=="OUT"]=2
  targets=targets[,c(2,4,3,1)]
  
  ###DE analysis with EdgeR
  d <- readDGE(targets, skip=0, comment.char = "!")
  d <- estimateCommonDisp(d, verbose = F)
  d <- estimateTagwiseDisp(d)
  et <- exactTest(d)
  # select top 
  top <- topTags(et,p.value = filter_pval,adjust.method ="BH",n=Inf)
  
  
  #extract count 
  counts_norm_log2=as.data.frame(log2(d$counts+1))
  colnames(counts_norm_log2)=d$samples$labels
  
  
  if(nrow(top)>0){
    #filter FC
    top$table$logFC=-top$table$logFC
    sel1=top$table$logFC>filter_FC
    top$table=top$table[sel1,]
    
    
    ##################
    ####Heatmap#######
    df <- as.data.frame(d$samples$description,row.names =d$samples$labels )
    colnames(df)[1]="condition"
    list_col=list(condition=c(IN="black",OUT="blue"))
    hmcol=colorRampPalette(c(
      rep("green",4),
      rep("darkgreen",2),
      rep("yellow",2),
      rep("orange",2),
      rep("orangered2",2),
      rep("red",8)
    ),space="rgb")(20)
    brs=0:20
    kins=rownames(top$table)
    kins2=intersect(kins,rownames(counts_norm_log2))
    data_log1=counts_norm_log2[kins2,]
    byFC=rownames(top$table)[order((top$table$logFC),decreasing = T)]
    data_log1=data_log1[byFC,]
    data_log1=data_log1[,c(which(stri_detect_fixed(colnames(data_log1),"IN")),which(stri_detect_fixed(colnames(data_log1),"OUT")))]
    wd = 15*length(colnames(data_log1))+300
    ht = 12*length(kins)+250
    if(length(kins)>30){cw=15;ch=8;fn=8}#;wd=600;ht=1200}
    if(length(kins)<=30){cw=30;ch=25;fn=10}#;wd=600;ht=1024}
    if(length(kins)<=10){cw=30;ch=30;fn=10}#;wd=500;ht=500}
    if(length(kins)<=5){cw=30;ch=30;fn=10}#;wd=360;ht=360}
    png(paste(dir_sample,"/Heatmap_counts2.png",sep=""),width=wd,height = ht)
    c=pheatmap(data_log1,annotation_col=df,color =hmcol ,scale="none",breaks = brs,
               fontsize_number=fn,
               border_color="white",
               number_color="#2F2F2F",
               cellwidth =cw,
               cellheight =ch,
               annotation_colors=list_col,display_numbers=round(data_log1,1),
               cluster_rows=F,cluster_cols = T)
    dev.off()
    
    
    #write results on file
    for(row_top in 1:nrow(top$table)){
      write(x = c(sample,
                  rownames(top$table)[row_top],
                  round(top$table$logFC[row_top],2),
                  top$table$PValue[row_top],
                  as.numeric(data_log1[rownames(top$table)[row_top],])),
            file = paste(path_output,out_fusion,sep="/"),ncolumns = 8,append = T,sep = "\t")
      x=as.data.frame(cbind(sample,
                            rownames(top$table)[row_top],
                            round(top$table$logFC[row_top],2),
                            top$table$PValue[row_top],
                            round(as.numeric(data_log1[rownames(top$table)[row_top],1]),1),
                            round(as.numeric(data_log1[rownames(top$table)[row_top],2]),1),
                            round(as.numeric(data_log1[rownames(top$table)[row_top],3]),1),
                            round(as.numeric(data_log1[rownames(top$table)[row_top],4]),1)))
      colnames(x)=c("Sample","Kinase","FC","p-value","IN_A","IN_B","OUT_A","OUT_B")
      print(x)
      print("")
    }
  }else{
    write(x = c(sample,"","","","","","",""),
          file = paste(path_output,out_fusion,sep="/"),ncolumns = 8,append = T,sep = "\t")
    
  }
}

