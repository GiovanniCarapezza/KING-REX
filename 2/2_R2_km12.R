library(stringi)

rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2/"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/2",sep=""))

#scatterplot and rsquared calculation 
R2=function(a,b,x,y,mn){
  plot(a,b,pch=20,col="blue",xlab = x,ylab = y,main = mn,cex.main=6,xlim =c(0,17.5),ylim =c(0,17.5))
  abline(lm(b~a),col="red")
  r2_1=summary(lm(b~a))
  text(x = max(a)-2,y=min(b)+0.2,expression(paste("",R^2,": ",sep="")))
  text(x = max(a)-0.8,y=min(b)+0.1,paste(round(r2_1$r.squared,3)))
  return(round(r2_1$r.squared,4))
}

#
KM12s=c("KM12_A","KM12_B","KM12T0_A","KM12T0_B")


#load king-rex data
kingrex_genes=read.table(file="normalized_counts_KING-REX/normalized_counts_log2.txt",sep="\t",header=T,row.names=1)
#select KM12s
sel=colnames(kingrex_genes) %in% KM12s
kingrex_genes2=kingrex_genes[,sel]

if(!dir.exists("R2_km12")){
  dir.create("R2_km12")
}

R2_matrix=as.data.frame(matrix(0,nrow = length(KM12s),ncol = length(KM12s)))
colnames(R2_matrix)=KM12s
rownames(R2_matrix)=KM12s
for(cl in KM12s){
  for(cl2 in KM12s){
    a=kingrex_genes2[,cl]
    b=kingrex_genes2[,cl2]
    x=cl
    y=cl2
    mn = ""
    jpeg(filename = paste("R2_km12/",cl,"_",cl2,".jpeg",sep=""))
    R2_matrix[cl,cl2]=R2(a,b,x,y,mn)
    dev.off() 
  }
}

print(R2_matrix)
print("Average R squared value for KM12 cell lines:")
print(round(mean(R2_matrix[R2_matrix<1]),2))

