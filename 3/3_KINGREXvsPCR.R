library(stringi)
library(reshape2)
library(ggplot2)
library(xlsx)

rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/3",sep=""))
#path out
path_output=paste("Fig_S1_KR_vs_PCR/",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}



##
norm_in_out=read.csv(file = "normalized_counts_IN_OUT/normalized_IN_OUT_counts_log2.txt",header = T,sep = "\t",row.names = 1)

#sel genes
genes=c("ALK","NTRK1","RET","ROS1")
norm_in_out2=norm_in_out[genes,]

#sel cell lines 
cls=c("KARPAS299","U118MG","KM12","LC2AD","UHRR")
col_app=stri_split_fixed(colnames(norm_in_out2),"_",simplify = T)[,1]
norm_in_out3=norm_in_out2[,col_app %in% cls]

#fusions
fusions=c("KARPAS299","KM12","LC2AD","U118MG")
names(fusions)=genes
#sort cell lines
col_app2=stri_split_fixed(colnames(norm_in_out3),"_",simplify = T)[,1]
index_sort=c()
for(i in 1:length(cls)){
  index_sort=c(index_sort,which(col_app2==cls[i]))
}
norm_in_out4=norm_in_out3[,index_sort[1:(ncol(norm_in_out3)-6)]]
#calcolo medie duplicat
norm_in_out5=matrix(0,nrow = length(genes),ncol = length(cls)*2)
colnames(norm_in_out5)=paste(rep(cls,each=2),c("_IN","_OUT"),sep="")
rownames(norm_in_out5)=genes
for(cl in cls){
  for(g in genes){
    app1=(2^(norm_in_out4[g,paste(cl,"_A_IN",sep="")])-1+(2^(norm_in_out4[g,paste(cl,"_B_IN",sep="")])-1))/2
    norm_in_out5[g,paste(cl,"_IN",sep="")]=log(x = (app1+1),base = 2)
    app2=(2^(norm_in_out4[g,paste(cl,"_A_OUT",sep="")])-1+(2^(norm_in_out4[g,paste(cl,"_B_OUT",sep="")])-1))/2
    norm_in_out5[g,paste(cl,"_OUT",sep="")]=log(x = (app2+1),base = 2)
  }
}
#relativ quant vs. UHRR
norm_in_out6=matrix(0,nrow = length(genes),ncol = length(cls)*2)
colnames(norm_in_out6)=paste(rep(cls,each=2),c("_IN","_OUT"),sep="")
rownames(norm_in_out6)=genes
for(cl in cls){
  for(g in genes){
    app1=2^((norm_in_out5[g,paste(cl,"_IN",sep="")]-norm_in_out5[g,"UHRR_IN"]))
    norm_in_out6[g,paste(cl,"_IN",sep="")]=app1
    app2=2^((norm_in_out5[g,paste(cl,"_OUT",sep="")]-norm_in_out5[g,"UHRR_OUT"]))
    norm_in_out6[g,paste(cl,"_OUT",sep="")]=app2
  }
}

#read pcr 
pcr_data=read.xlsx(file = "Summary qPCR expression data KING-REX validation.xlsx",sheetIndex = 1,as.data.frame = T)
pcr_data2=t(pcr_data)
colnames(pcr_data2)=pcr_data2[1,]
pcr_data2=pcr_data2[2:11,]
pcr_data3=pcr_data2[c(2,6,8,10),c(1,2,4,5,7)]
colnames(pcr_data3)=paste(cls,"_IN",sep="")
rownames(pcr_data3)=genes
pcr_data4=pcr_data2[c(1,5,7,9),c(1,2,4,5,7)]
colnames(pcr_data4)=paste(cls,"_OUT",sep="")
rownames(pcr_data4)=genes
pcr_data5=cbind(pcr_data3,pcr_data4)
pcr_data6=pcr_data5[,paste(rep(cls,each=2),c("_IN","_OUT"),sep="")]

melt_pcr_data=melt(pcr_data6)
melt_pcr_data[,c("Cell_line","Assay")]=stri_split_fixed(melt_pcr_data$Var2,"_",simplify = T)
colnames(melt_pcr_data)[c(1,3)]=c("Gene","Relative Quantification (UHRR)")
melt_pcr_data=melt_pcr_data[,c(1,4,5,3)]

melt_kingrex_data=melt(norm_in_out6)
melt_kingrex_data[,c("Cell_line","Assay")]=stri_split_fixed(melt_kingrex_data$Var2,"_",simplify = T)
colnames(melt_kingrex_data)[c(1,3)]=c("Gene","Relative Quantification (UHRR)")
melt_kingrex_data=melt_kingrex_data[,c(1,4,5,3)]

for(gene in genes){
  #format plot 
  ggbase1=  geom_bar(position="dodge",stat="identity",width=0.7) 
  ggbase2=scale_fill_manual(name="Assay",labels=c("IN","OUT"),values=c("royalblue4","red"))
  ggbase3= theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
                 legend.position="right", 
                 text=element_text(size=24),
                 plot.title = element_text(hjust = 0.5)#,
                 #axis.title.x=element_text(size=0)
  )
  # ggbase4=coord_cartesian(ylim=c(0,20))
  
  #Fig. 4B
  melt_dataset=melt_kingrex_data[which(melt_kingrex_data$Gene==gene),]
  melt_dataset=melt_dataset[which(melt_dataset$Cell_line %in% c(fusions[gene],"UHRR")),]
  melt_dataset$Cell_line=factor(x = melt_dataset$Cell_line,levels = cls,ordered = T)
  melt_dataset$Assay=factor(x = melt_dataset$Assay,levels = c("IN","OUT"),ordered = T)
  gg3<-ggplot(melt_dataset, aes(x=Cell_line, y=`Relative Quantification (UHRR)`, fill=Assay)) + 
    ggbase1 + ggbase2 + ggbase3  + ggtitle(paste("KINGREX:",gene)) 
  
  # p=p+geom_hline(aes(yintercept = 7.7,linetype="200 reads"),col="green")
  # p=p+geom_hline(aes(yintercept = 6.7,linetype="100 reads"),col="orange")
  # p=p+geom_hline(aes(yintercept = 5.7,linetype="50 reads"),col="red")
  # p=p+scale_linetype_manual(name="",labels=c("200 reads","100 reads","50 reads") ,values = c(1,1,1), guide = guide_legend(override.aes = list(color = c("green","orange","red"))))
  # gg3=gg3+geom_hline(aes(yintercept = 5,linetype="Threshold"),col="green")
  # 
  #save plot
  jpeg(paste(path_output,gene,"_kingrex.jpeg",sep=""),width = 700,height = 500)
  print(gg3)
  dev.off()
  
}


for(gene in genes){
  #format plot 
  ggbase1=  geom_bar(position="dodge",stat="identity",width=0.7) 
  ggbase2=scale_fill_manual(name="Assay",labels=c("IN","OUT"),values=c("royalblue4","red"))
  ggbase3= theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
                 legend.position="right", 
                 text=element_text(size=24),
                 plot.title = element_text(hjust = 0.5)#,
                 #axis.title.x=element_text(size=0)
  )
  # ggbase4=coord_cartesian(ylim=c(0,20))
  
  #Fig. 4B
  melt_dataset=melt_pcr_data[which(melt_pcr_data$Gene==gene),]
  melt_dataset=melt_dataset[which(melt_dataset$Cell_line %in% c(fusions[gene],"UHRR")),]
  melt_dataset$Cell_line=factor(x = melt_dataset$Cell_line,levels = cls,ordered = T)
  melt_dataset$Assay=factor(x = melt_dataset$Assay,levels = c("IN","OUT"),ordered = T)
  melt_dataset$`Relative Quantification (UHRR)`=as.numeric(as.character(melt_dataset$`Relative Quantification (UHRR)`))
  gg3<-ggplot(melt_dataset, aes(x=Cell_line, y=`Relative Quantification (UHRR)`, fill=Assay))+ 
    ggbase1 + ggbase2 + ggbase3 + ggtitle(paste("qPCR:",gene)) 
  
  
  
  #save plot
  jpeg(paste(path_output,gene,"_pcr.jpeg",sep=""),width = 700,height = 500)
  print(gg3)
  dev.off()
  
}

b=as.numeric(as.character(melt_pcr_data$`Relative Quantification (UHRR)`))
a=melt_kingrex_data$`Relative Quantification (UHRR)`
r2_1=summary(lm(b~a))
print("R squared value KINGREX vs PCR data:")
print(round(r2_1$r.squared,4))

