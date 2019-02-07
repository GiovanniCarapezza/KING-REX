library(stringi)
#to excecute this script, you must run first script in folder 2
rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC review/"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/6",sep=""))
#path to norm counts
path_input_data=paste(sep="",path_to_results,"/Results/2/normalized_counts_KING-REX")
#path output
path_output=paste("R2",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}

#read annotation file
input=read.table(paste(path_to_results,"/Results/annotations_KING-REX.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)

#load normalized counts
kingrex_genes=read.table(file=paste(path_input_data,"normalized_counts_log2.txt",sep="/"),sep="\t",header=T,row.names=1)

#function to calculate R2 and scatterplot
R2=function(a,b,x,y,mn){
  plot(a,b,pch=20,col="blue",xlab = x,ylab = y,main = mn,cex.main=2,cex.axis=2,cex.lab=1.6,xlim =c(0,17.5),ylim =c(0,17.5))
  abline(lm(b~a),col="red")
  r2_1=summary(lm(b~a))
  text(x = max(a)-2,y=min(b)+0.2,expression(paste("",R^2,": ",sep="")),cex = 2)
  text(x = max(a)-0.8,y=min(b)+0.1,paste(round(r2_1$r.squared,2)),cex=2)
  print(mn)
  print(paste("R^2:",round(r2_1$r.squared,2)))
  print("")
  
}


#function to pass duplicate mean to R2 function
myR2_mean=function(s1,s2,mn){
  sel_kingrex_genes=stri_detect_fixed(colnames(kingrex_genes),s1)
  kingrex_genes2=kingrex_genes[,sel_kingrex_genes]
  sel_kingrex_genes=stri_detect_fixed(colnames(kingrex_genes),s2)
  kingrex_genes3=kingrex_genes[,sel_kingrex_genes]
  jpeg(filename = paste(path_output,"/",s1[1],"_",s2[1],".jpeg",sep=""),width = 700,height = 600)
  R2(rowMeans(kingrex_genes2),rowMeans(kingrex_genes3),paste("Mean","_",s1[1],"_and_B",sep=""),paste("Mean","_",s2[1],"_and_B",sep=""),mn)
  dev.off()
}


#T0 vs T10
s1=c("KM12T0_A","KM12T0_B")
s2=c("KM12T10_A","KM12T10_B")
mn="T0 vs. T10"
myR2_mean(s1,s2,mn)
#T0 vs T20
s1=c("KM12T0_A","KM12T0_B")
s2=c("KM12T20_A","KM12T20_B")
mn="T0 vs. T20"
myR2_mean(s1,s2,mn)
#T0 vs T60
s1=c("KM12T0_A","KM12T0_B")
s2=c("KM12T60_A","KM12T60_B")
mn="T0 vs. T60"
myR2_mean(s1,s2,mn)

