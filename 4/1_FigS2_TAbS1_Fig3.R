library(stringi)
library(pheatmap)
library(RColorBrewer)
#to excecute this script, you must run first script in folder 2
rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/4",sep=""))
#path to norm counts
path_input_data=paste(sep="",path_to_results,"/Results/2/normalized_counts_KING-REX")
#path output
path_output=paste("output",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}

#read kinase list
kin_list=read.table(paste(path_to_results,"/Results/1/kinase_list_and_info.txt",sep=""),sep="\t",stringsAsFactors=F,row.names = 1,header = T)
kin_list=kin_list[sort(rownames(kin_list)),]

#read annotation file
input=read.table(paste(path_to_results,"/Results/annotations_KING-REX.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)

#load normalized counts
kingrex_genes=read.table(file=paste(path_input_data,"normalized_counts_log2.txt",sep="/"),sep="\t",header=T,row.names=1)
kingrex_genes_not_log2=read.table(file=paste(path_input_data,"normalized_counts.txt",sep="/"),sep="\t",header=T,row.names=1)

#annotation
annotation<-as.data.frame(cbind(input$Tissue),row.names=input$Sample,stringsAsFactors = F)
#selct dilution
indexes=c(which(stri_startswith_fixed(rownames(annotation),"KAR")),
          which(stri_startswith_fixed(rownames(annotation),"U118")),
          which(stri_startswith_fixed(rownames(annotation),"KU")),
          which(stri_startswith_fixed(rownames(annotation),"UK")))
annotation<-as.data.frame(cbind(input$Tissue[indexes]),row.names=input$Sample[indexes],stringsAsFactors = F)
ann_colors=list(
  V1=c(Mix="blue",
       NERVOUS_SYSTEM="purple",LYMPHOMA="lightblue")
)
kingrex_genes=kingrex_genes[,colnames(kingrex_genes) %in% rownames(annotation)]
kingrex_genes_not_log2=kingrex_genes_not_log2[,colnames(kingrex_genes_not_log2) %in% rownames(annotation)]

#distance matrix of figure S2 on dilutions
distscn <- dist(t(kingrex_genes))
mat <- as.matrix(distscn)
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
wd = 35*nrow(annotation)
jpeg(paste(path_output,"FigS2_distance_matrix_dilutions.jpeg",sep="/"),width = wd+400,height = wd+150)
pheatmap(mat,fontsize = 24,cellwidth = 20,cellheight = 20,
         clustering_distance_rows=distscn,
         annotation_col=annotation,
         annotation_colors = ann_colors,
         clustering_distance_cols=distscn,
         col = rev(hmcol))
dev.off()


###select genes not expressed in karpas
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
sel_0=apply(kingrex_genes[,c("KARPAS299_A","KARPAS299_B")],MARGIN = 1,mean)<1 
U118MG=kingrex_genes[sel_0,]
#select genes with expression >5 in u118mg
sel_1=apply(U118MG[,c("U118MG_A","U118MG_B")],MARGIN = 1,mean)>5
U118MG=U118MG[sel_1,]
U118MG=U118MG[order(U118MG$U118MG_A,decreasing = T),]
U118MG=U118MG[,c("U118MG_A","U118MG_B",
                 "UK8712_A","UK8712_B",
                 "UK7525_A","UK7525_B",
                 "KU50_A","KU50_B",
                 "KU7525_A","KU7525_B",
                 "KU8712_A","KU8712_B",
                 "KARPAS299_A","KARPAS299_B")]

#heatmap
hg = 30*nrow(U118MG)+100
wd = 30*ncol(U118MG)+100
jpeg(paste(path_output,"TabS1_A.jpeg",sep="/"),width = wd,height = hg)
pheatmap(U118MG,show_rownames = T,
         color = hmcol,
         clustering_distance_rows = distanza,
         clustering_distance_cols = distanza,
         scale="none",
         display_numbers = round(U118MG,1),
         fontsize_number=10,
         border_color="white",
         number_color="black",
         cellwidth =22,
         cellheight = 22,
         cluster_rows=F,
         cluster_cols=F)
dev.off()

###select genes not expressed in u118mg
sel_0=apply(kingrex_genes[,c("U118MG_A","U118MG_B")],MARGIN = 1,mean)<1
KARPAS=kingrex_genes[sel_0,]
##select genes with expression >5 in Karpas
sel_1=apply(KARPAS[,c("KARPAS299_A","KARPAS299_B")],MARGIN = 1,mean)>5
KARPAS=KARPAS[sel_1,]
KARPAS=KARPAS[order(KARPAS$KARPAS299_A,decreasing = T),]
KARPAS=KARPAS[,c("KARPAS299_A","KARPAS299_B",
                 "KU8712_A","KU8712_B",
                 "KU7525_A","KU7525_B",
                 "KU50_A","KU50_B",
                 "UK7525_A","UK7525_B",
                 "UK8712_A","UK8712_B",
                 "U118MG_A","U118MG_B")]

hg = 30*nrow(KARPAS)+100
wd = 30*ncol(KARPAS)+100
jpeg(paste(path_output,"TabS1_B.jpeg",sep="/"),width = wd,height = hg)
pheatmap(KARPAS,show_rownames = T,
         color = hmcol,
         clustering_distance_rows = distanza,
         clustering_distance_cols = distanza,
         scale="none",
         display_numbers = round(KARPAS,1),
         fontsize_number=10,
         border_color="white",
         number_color="black",
         cellwidth =22,
         cellheight = 22,
         cluster_rows=F,
         cluster_cols=F)
dev.off()
#write to file
write.table(x = U118MG,file = paste(path_output,"TabS1_A.txt",sep="/"),append = F,quote = F,sep = "\t",row.names = T,col.names = T)
write.table(x = KARPAS,file = paste(path_output,"TabS1_B.txt",sep="/"),append = F,quote = F,sep = "\t",row.names = T,col.names = T)


###R2 Figure 3 Theoretical versus measured gene expression value
#function to scatterplot and R2 calculation
R2=function(a,b,x,y,mn){
  plot(a,b,pch=20,col="blue",xlab = x,ylab = y,main = mn,cex.main=2,cex.axis=2,cex.lab=1.6,xlim =c(0,17.5),ylim =c(0,17.5))
  abline(lm(b~a),col="red")
  r2_1=summary(lm(b~a))
  text(x = max(a)-2,y=min(b)+0.2,expression(paste("",R^2,": ",sep="")),cex = 2)
  text(x = max(a)-0.8,y=min(b)+0.1,paste(round(r2_1$r.squared,2)),cex=2)
}

#function to calculate theoretical values
calcola_expr=function(a,b,perc1){
  perc2=1-perc1
  prev=a*perc1+b*perc2
  return(prev)
}

#mean of measured values 
KARPAS299=rowMeans(kingrex_genes_not_log2[,c("KARPAS299_A","KARPAS299_B")])
U118MG=rowMeans(kingrex_genes_not_log2[,c("U118MG_A","U118MG_B")])
KU8712=rowMeans(kingrex_genes_not_log2[,c("KU8712_A","KU8712_B")])
KU7525=rowMeans(kingrex_genes_not_log2[,c("KU7525_A","KU7525_B")])
KU50=rowMeans(kingrex_genes_not_log2[,c("KU50_A","KU50_B")])
UK7525=rowMeans(kingrex_genes_not_log2[,c("UK7525_A","UK7525_B")])
UK8712=rowMeans(kingrex_genes_not_log2[,c("UK8712_A","UK8712_B")])
#log2 of measured values
x=c(log(x=(KU8712+1),base=2),
    log(x=(KU7525+1),base=2),
    log(x=(KU50+1),base=2),
    log(x=(UK7525+1),base=2),
    log(x=(UK8712+1),base=2))

#theoretical value calculation and log2 transformation
y=c(log(x=(calcola_expr(KARPAS299,U118MG,0.875)+1),base=2),
    log(x=(calcola_expr(KARPAS299,U118MG,0.75)+1),base=2),
    log(x=(calcola_expr(KARPAS299,U118MG,0.5)+1),base=2),
    log(x=(calcola_expr(KARPAS299,U118MG,0.25)+1),base=2),
    log(x=(calcola_expr(KARPAS299,U118MG,0.125)+1),base=2))

#scatterplot and R2
s1="Measured"
s2="Theoretical"
mn = "Measured vs. Theoretical"
jpeg(filename = paste(path_output,"/Fig3_",s1,"_",s2,"_scatterplot.jpeg",sep=""),width = 700,height = 600)
R2(x,y,s1,s2,mn)
dev.off()

