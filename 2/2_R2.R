library(stringi)

rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/2",sep=""))

#scatterplot and rsquared calculation 
R2=function(a,b,x,y,mn){
  plot(a,b,pch=20,col="blue",xlab = x,ylab = y,main = mn,cex.main=6,xlim =c(0,17.5),ylim =c(0,17.5))
  abline(lm(b~a),col="red")
  r2_1=summary(lm(b~a))
  text(x = max(a)-2,y=min(b)+0.2,expression(paste("",R^2,": ",sep="")))
  text(x = max(a)-0.8,y=min(b)+0.1,paste(round(r2_1$r.squared,2)))
  return(round(r2_1$r.squared,2))
}

#
large_intestine_cell_lines=c("COLO205","COLO678","HCT116","HCT15",
                             "KM12","LS180","RKO","SW1417","SW480","SW948")


#load king-rex data
kingrex_genes=read.table(file="normalized_counts_KING-REX/normalized_counts_log2.txt",sep="\t",header=T,row.names=1)
#select large intestine cell line
sel=colnames(kingrex_genes) %in% paste(large_intestine_cell_lines,"A",sep="_")
kingrex_genes2=kingrex_genes[,sel]
#sort columns
kingrex_genes2=kingrex_genes2[,sort(colnames(kingrex_genes2))]
#remove _A from name 
colnames(kingrex_genes2)=stri_replace_all_fixed(colnames(kingrex_genes2),"_A","")

#load transcriptome data 
transcriptome_genes=read.table(file="normalized_counts_transcriptome/normalized_counts_log2.txt",sep="\t",header=T,row.names=1)
#select kinase names
transcriptome_genes2=transcriptome_genes[rownames(kingrex_genes),]

if(!dir.exists("R2")){
  dir.create("R2")
}

R2_matrix=as.data.frame(matrix(0,nrow = length(large_intestine_cell_lines),ncol = 1))
colnames(R2_matrix)="R^2"
rownames(R2_matrix)=large_intestine_cell_lines
for(cl in large_intestine_cell_lines){
  a=kingrex_genes2[,cl]
  b=transcriptome_genes2[,cl]
  x="KING-REX"
  y="Transcriptome"
  mn = cl
  jpeg(filename = paste("R2/",cl,".jpeg",sep=""))
  R2_matrix[cl,]=R2(a,b,x,y,mn)
  dev.off()
}

print(R2_matrix)
print("Average R squared value:")
print(round(mean(R2_matrix$`R^2`),1))

