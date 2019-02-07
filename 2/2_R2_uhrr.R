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
  text(x = max(a)-0.8,y=min(b)+0.1,paste(round(r2_1$r.squared,2)))
  return(round(r2_1$r.squared,2))
}

#
UHRRs=c("UHRR_A","UHRR_B","UHRR_C","UHRR_D","UHRR_E")


#load king-rex data
kingrex_genes=read.table(file="normalized_counts_KING-REX/normalized_counts_log2.txt",sep="\t",header=T,row.names=1)
#select UHRRS
sel=colnames(kingrex_genes) %in% UHRRs
kingrex_genes2=kingrex_genes[,sel]

if(!dir.exists("R2_uhrr")){
  dir.create("R2_uhrr")
}

R2_matrix=as.data.frame(matrix(0,nrow = length(UHRRs),ncol = length(UHRRs)))
colnames(R2_matrix)=UHRRs
rownames(R2_matrix)=UHRRs
for(cl in UHRRs){
  for(cl2 in UHRRs){
    a=kingrex_genes2[,cl]
    b=kingrex_genes2[,cl2]
    x=cl
    y=cl2
    mn = ""
    jpeg(filename = paste("R2_uhrr/",cl,"_",cl2,".jpeg",sep=""))
    R2_matrix[cl,cl2]=R2(a,b,x,y,mn)
    dev.off() 
  }
}

print(R2_matrix)
print("Average R squared value for UHRR cell lines:")
print(round(mean(R2_matrix[R2_matrix<1]),2))

