library(stringi)
library(xlsx)
# library(ggplot2)
# library(reshape2)

rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2/"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/2",sep=""))
#path to raw counts
path_input_data=paste(sep="",path_to_results,"/Results/raw_counts_KING-REX/")
#path to raw counts
path_input_data2=paste(sep="",path_to_results,"/Results/raw_counts_transcriptome/")
#path to raw counts
path_input_data3=paste(sep="",path_to_results,"/Results/raw_counts_CCLE/")


#load kingrex raw count 
#read annotation file
kingrex_annot=read.table(paste(path_to_results,"/Results/annotations_KING-REX.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)
kingrex_data=matrix(data = 0,nrow = 28517,ncol=length(unique(kingrex_annot$Sample_Name)))
colnames(kingrex_data)=unique(kingrex_annot$Sample_Name)
for(i in 1:length(unique(kingrex_annot$Sample_Name))){
  n=unique(kingrex_annot$Sample_Name)[i]
  #print(n)
  a=read.table(file = paste(sep="",path_input_data,n,"_A.genes.results"),header = T,sep = "\t",stringsAsFactors = F)
  kingrex_data[,n]=a$expected_count
  if(i==1){
    rownames(kingrex_data)=a$gene_id
  }
}
kingrex_data=as.data.frame(kingrex_data)

#load transcriptome raw count 
transcriptome_annot=read.table(paste(path_to_results,"/Results/annotations_Transcriptome.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)
transcriptome_data=matrix(data = 0,nrow = 28517,ncol=length(unique(transcriptome_annot$Sample_Name)))
colnames(transcriptome_data)=unique(transcriptome_annot$Sample_Name)
for(i in 1:length(unique(transcriptome_annot$Sample_Name))){
  n=unique(transcriptome_annot$Sample_Name)[i]
  #print(n)
  a=read.table(file = paste(sep="",path_input_data2,n,".genes.results"),header = T,sep = "\t",stringsAsFactors = F)
  transcriptome_data[,n]=a$expected_count
  if(i==1){
    rownames(transcriptome_data)=a$gene_id
  }
}
transcriptome_data=as.data.frame(transcriptome_data)

#load ccle data
ccle_data=read.csv(file = paste(sep="",path_input_data3,"/CCLE_RNAseq_15_Aug_2017.reads.gct2"),header = T,check.names = F,sep = "\t",stringsAsFactors = F,skip = 2)


#################################################################################
#read kinase list
list_kin=read.table(paste(path_to_results,"/Results/1/kinase_list_and_info.txt",sep=""),sep="\t",stringsAsFactors=F,row.names = 1,header = T)
list_kin=list_kin[sort(rownames(list_kin)),]

#change name
list_kin=list_kin[sort(rownames(list_kin)),]
change_name=which(rownames(kingrex_data)=="SGK110")
rownames(kingrex_data)[change_name]="SBK3"
change_name=which(rownames(kingrex_data)=="SGK196")
rownames(kingrex_data)[change_name]="POMK"
change_name=which(rownames(kingrex_data)=="KIAA1804")
rownames(kingrex_data)[change_name]="MAP3K21"
rownames(transcriptome_data)=rownames(kingrex_data)

#extract only kinase and sort
transcriptome_data_kin=transcriptome_data[rownames(list_kin),]
kingrex_data_kin=kingrex_data[rownames(list_kin),]
ccle_data_kin=ccle_data[ccle_data$Description %in% rownames(list_kin),]
rownames(ccle_data_kin)=ccle_data_kin$Description
#setdiff(rownames(kingrex_data_kin),rownames(ccle_data_kin))

ccle_data_kin=ccle_data[ccle_data$Description %in% c(rownames(list_kin),"MLTK","MLK4"),]
rownames(ccle_data_kin)=ccle_data_kin$Description
change_name=which(rownames(ccle_data_kin)=="MLTK")
rownames(ccle_data_kin)[change_name]="ZAK"
change_name=which(rownames(ccle_data_kin)=="MLK4")
rownames(ccle_data_kin)[change_name]="MAP3K21"

ccle_data_kin=ccle_data_kin[sort(rownames(ccle_data_kin)),3:ncol(ccle_data_kin)]
kingrex_data_kin=kingrex_data_kin[rownames(ccle_data_kin),sort(colnames(kingrex_data_kin))]
transcriptome_data_kin=transcriptome_data_kin[rownames(ccle_data_kin),sort(colnames(transcriptome_data_kin))]

#normaliztion function
UQ1000=function(a){
  b=a[a!=0]
  c=quantile(b,probs = 0.75)
  d=1000/c
  e=a*d
  return(e)
}

##Normalization
transcriptome_data_norm=transcriptome_data_kin
for(i in 1:ncol(transcriptome_data_kin)){
  to_norm=as.numeric(transcriptome_data_kin[,i])
  norm=UQ1000(to_norm)
  transcriptome_data_norm[,i]=log(x = norm+1,base = 2)
}

ccle_data_norm=ccle_data_kin
ccle_data_raw=ccle_data_kin
for(i in 1:ncol(ccle_data_kin)){
  to_norm=as.numeric(ccle_data_kin[,i])
  ccle_data_raw[,i]=log(x = to_norm+1,base = 2)
  norm=UQ1000(to_norm)
  ccle_data_norm[,i]=log(x = norm+1,base = 2)
}

kingrex_data_norm=kingrex_data_kin
for(i in 1:ncol(kingrex_data_kin)){
  to_norm=as.numeric(kingrex_data_kin[,i])
  norm=UQ1000(to_norm)
  kingrex_data_norm[,i]=log(x = norm+1,base = 2)
}


#select large intestine cell line from transcriptome and kingrex
name2=c("COLO205","COLO678","HCT116","HCT15","KM12","LS180","RKO","SW1417","SW480","SW948")
KINGREX_LI=kingrex_data_norm[,name2]
transcriptome_LI=transcriptome_data_norm[,name2]

#select cell line from ccle and king-rex
indx=c()
name=c()
name2=c()
for(i in 1:ncol(kingrex_data_norm)){
  n=colnames(kingrex_data_norm)[i]
  a=which(stri_detect_fixed(colnames(ccle_data_norm),n))
  if(length(a)>0){
    indx=c(indx,a)
    name=c(name,n)
  }else{
    name2=c(name2,n)
  }
}
ccle_data_norm=ccle_data_norm[,indx]
colnames(ccle_data_norm)=name

#large intestine cell lines to select
name2=c("COLO678","HCT116","HCT15","KM12","LS180","RKO","SW1417","SW480","SW948")

KINGREX_LI2=kingrex_data_norm[,name2]
CCLE_LI=ccle_data_norm[,name2]

#mix cell lines to select
name2=c("BT474","HPAC","HUH7","K562","KARPAS299","NCIH716","SNU1079","U118MG")
KINGREX_Mix=kingrex_data_norm[,name2]
CCLE_Mix=ccle_data_norm[,name2]

#remove all, except for dataset 
a=ls()
b=c("KINGREX_LI","KINGREX_LI2","KINGREX_Mix","CCLE_LI","CCLE_Mix","transcriptome_LI")
a=a[!(a %in% b)]
rm(list=a)
rm(b)


#function to calulate sentivity specificity and accuracy for each kinase
test_metrics=function(trsh,KINGREX,Reference){
  
  test_P_N=matrix("",nrow = length(rownames(KINGREX)),ncol=length(colnames(KINGREX)))
  rownames(test_P_N)=rownames(KINGREX)
  colnames(test_P_N)=colnames(KINGREX)
  for(i in 1:length(rownames(KINGREX))){
    kin=rownames(test_P_N)[i]
    for(j in 1:length(colnames(KINGREX))){
      cl=colnames(test_P_N)[j]
      if(Reference[kin,cl]>trsh){test_P_N[i,j]="P"}else{test_P_N[i,j]="N"}
    }
  }
  
  test_TP_TN_FP_FN=matrix("",nrow = length(rownames(KINGREX)),ncol=length(colnames(KINGREX)))
  rownames(test_TP_TN_FP_FN)=rownames(KINGREX)
  colnames(test_TP_TN_FP_FN)=colnames(KINGREX)
  for(i in 1:length(rownames(KINGREX))){
    kin=rownames(test_TP_TN_FP_FN)[i]
    for(j in 1:length(colnames(KINGREX))){
      cl=colnames(test_TP_TN_FP_FN)[j]
      if(Reference[kin,cl]>trsh & KINGREX[kin,cl]>trsh){
        test_TP_TN_FP_FN[i,j]="TP"
      }
      if(Reference[kin,cl]<=trsh & KINGREX[kin,cl]<=trsh){
        test_TP_TN_FP_FN[i,j]="TN"
      }
      if(Reference[kin,cl]<=trsh & KINGREX[kin,cl]>trsh){
        test_TP_TN_FP_FN[i,j]="FP"
      }
      if(Reference[kin,cl]>trsh & KINGREX[kin,cl]<=trsh){
        test_TP_TN_FP_FN[i,j]="FN"
      }
    }
  }
  
  
  
  test=matrix(0,nrow = length(rownames(KINGREX)),ncol=6)
  rownames(test)=rownames(KINGREX)
  colnames(test)=c("P","N","TP","FP","TN","FN")
  for(i in 1:length(rownames(KINGREX))){
    test[i,"P"]=length(which(test_P_N[i,]=="P"))
    test[i,"N"]=length(which(test_P_N[i,]=="N"))
    test[i,"TP"]=length(which(test_TP_TN_FP_FN[i,]=="TP"))
    test[i,"FP"]=length(which(test_TP_TN_FP_FN[i,]=="FP"))
    test[i,"TN"]=length(which(test_TP_TN_FP_FN[i,]=="TN"))
    test[i,"FN"]=length(which(test_TP_TN_FP_FN[i,]=="FN"))
  }
  test=as.data.frame(test)
  test[,"Threshold"]=trsh
  test[,"Recall"]=round((test$TP/test$P)*100,1)
  test[,"Precision"]=round((test$TP/(test$TP+test$FP))*100,1)
  test[,"Fmeasure"]=round((2*test$TP/(2*test$TP+test$FP+test$FN))*100,1)
  #######Note sensitivity == recall
  #######
  test=cbind(rownames(test),test)
  colnames(test)[1]="Kinases"
  ###############################################################
  ############add sheet with details for each Threshold##########à
  #################################################################
  ###########################################################à##
  sheet  <- createSheet(wb, sheetName=paste("Threshold_",trsh,sep=""))
  addDataFrame(test, sheet, startRow=1, startColumn=1,row.names = F,colnamesStyle = TABLE_COLNAMES_STYLE,characterNA ="NaN")
  rownames(test)=NULL
  return(test) 
}



#function to write average table 
to_wr=function(test_to_wr_app){
  a=test_to_wr_app$Threshold[1]
  b=round(mean(test_to_wr_app$Recall[!is.nan(test_to_wr_app$Recall)]),1)
  c=round(mean(test_to_wr_app$Precision[!is.nan(test_to_wr_app$Precision)]),1)
  d=round(mean(test_to_wr_app$Fmeasure[!is.nan(test_to_wr_app$Fmeasure)]),1)
  return(c(a,b,c,d))
}

#save Tab1 e tab2
wb2 <- createWorkbook("xlsx")
TABLE_COLNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 

#save details of table 1.A to excel file
wb <- createWorkbook("xlsx")
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 
#tab1.a 
print("Tab.1: KING-REX vs. In-house transcriptome data on a panel of CRC cancer cell lines")
test_05=test_metrics(0.5,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_1=test_metrics(1,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_2=test_metrics(2,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_3=test_metrics(3,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_4=test_metrics(4,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_5=test_metrics(5,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_6=test_metrics(6,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
test_7=test_metrics(7,KINGREX = KINGREX_LI,Reference = transcriptome_LI)
#save excel
saveWorkbook(wb, "Tab_S4_details_of_tab1a.xlsx")
rm(wb)

test_to_wr=matrix(data=0,nrow = 8,ncol = 4)
colnames(test_to_wr)=c("Threshold","Recall","Precision","Fmeasure")
test_to_wr[1,]=to_wr(test_05)
test_to_wr[2,]=to_wr(test_1)
test_to_wr[3,]=to_wr(test_2)
test_to_wr[4,]=to_wr(test_3)
test_to_wr[5,]=to_wr(test_4)
test_to_wr[6,]=to_wr(test_5)
test_to_wr[7,]=to_wr(test_6)
test_to_wr[8,]=to_wr(test_7)
print(test_to_wr)
sheet  <- createSheet(wb2, sheetName="Tab_1A")
addDataFrame(test_to_wr, sheet, startRow=1, startColumn=1,row.names = F,colnamesStyle = TABLE_COLNAMES_STYLE2,characterNA ="NaN")


wb <- createWorkbook("xlsx")
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 
print("Tab.1: KING-REX vs. CCLE transcriptome data on a panel of CRC cancer cell lines")
test_05=test_metrics(trsh = 0.5,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_1=test_metrics(1,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_2=test_metrics(2,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_3=test_metrics(3,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_4=test_metrics(4,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_5=test_metrics(5,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_6=test_metrics(6,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
test_7=test_metrics(7,KINGREX = KINGREX_LI2,Reference = CCLE_LI)
saveWorkbook(wb, "Tab_S5_details_of_tab1B.xlsx")
rm(wb)

test_to_wr=matrix(data=0,nrow = 8,ncol = 4)
colnames(test_to_wr)=c("Threshold","Recall","Precision","Fmeasure")
test_to_wr[1,]=to_wr(test_05)
test_to_wr[2,]=to_wr(test_1)
test_to_wr[3,]=to_wr(test_2)
test_to_wr[4,]=to_wr(test_3)
test_to_wr[5,]=to_wr(test_4)
test_to_wr[6,]=to_wr(test_5)
test_to_wr[7,]=to_wr(test_6)
test_to_wr[8,]=to_wr(test_7)
print(test_to_wr)
sheet  <- createSheet(wb2, sheetName="Tab_1B")
addDataFrame(test_to_wr, sheet, startRow=1, startColumn=1,row.names = F,colnamesStyle = TABLE_COLNAMES_STYLE2,characterNA ="NaN")

wb <- createWorkbook("xlsx")
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 
print("Tab.2: KING-REX vs. CCLE transcriptome data on a panel of heterogeneous cancer cell lines")
test_05=test_metrics(trsh = 0.5,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_1=test_metrics(1,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_2=test_metrics(2,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_3=test_metrics(3,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_4=test_metrics(4,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_5=test_metrics(5,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_6=test_metrics(6,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
test_7=test_metrics(7,KINGREX = KINGREX_Mix,Reference = CCLE_Mix)
saveWorkbook(wb, "Tab_S6_details_of_tab2.xlsx")
rm(wb)

test_to_wr=matrix(data=0,nrow = 8,ncol = 4)
colnames(test_to_wr)=c("Threshold","Recall","Precision","Fmeasure")
test_to_wr[1,]=to_wr(test_05)
test_to_wr[2,]=to_wr(test_1)
test_to_wr[3,]=to_wr(test_2)
test_to_wr[4,]=to_wr(test_3)
test_to_wr[5,]=to_wr(test_4)
test_to_wr[6,]=to_wr(test_5)
test_to_wr[7,]=to_wr(test_6)
test_to_wr[8,]=to_wr(test_7)
print(test_to_wr)
sheet  <- createSheet(wb2, sheetName="Tab_2")
addDataFrame(test_to_wr, sheet, startRow=1, startColumn=1,row.names = F,colnamesStyle = TABLE_COLNAMES_STYLE2,characterNA ="NaN")
saveWorkbook(wb2,file = "Tab1A_1B_2.xlsx")