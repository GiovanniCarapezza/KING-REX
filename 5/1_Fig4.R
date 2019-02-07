#to excecute this script, you must run first script in folder 3

library(stringi)
library(reshape2)
library(ggplot2)
rm(list = ls())
#substitute your path to results
path_to_results="~/Scrivania/KING-REX data and scripts for BMC2"
#set working directory to this source file location
setwd(paste(path_to_results,"/Results/5",sep=""))
#path to norm counts
path_input_data=paste(sep="",path_to_results,"/Results/3/normalized_counts_IN_OUT")
#path output
path_output=paste("output",sep="")
if(!dir.exists(path_output)){dir.create(path_output,recursive = T)}

#read kinase list
kin_list=read.table(paste(path_to_results,"/Results/1/kinase_list_and_info.txt",sep=""),sep="\t",stringsAsFactors=F,row.names = 1,header = T)
kin_list=kin_list[sort(rownames(kin_list)),]

#read annotation file
input=read.table(paste(path_to_results,"/Results/annotations_KING-REX.txt",sep=""),sep="\t",stringsAsFactors=F,header = T)
input2=read.table(paste(path_to_results,"/Results/3/edgeR_input/condition",sep=""),sep="\t",stringsAsFactors=F,header = T)

#load normalized counts for assay in e out
kingrex_InOut=read.table(file=paste(path_input_data,"normalized_IN_OUT_counts_log2.txt",sep="/"),sep="\t",header=T,row.names=1)

#Add tissue to input2 from input
input2[,"Tissue"]=""
for(row in 1:nrow(input2)){
  input2$Tissue[row]=input$Tissue[which(input$Sample_Name==input2$names[row])[1]]
}

#extract ros1 expression
U118MG=kingrex_InOut["ROS1",]
#select dilutions and order from 100%u118mg to 0% and IN-OUT
U118MG=U118MG[,c("U118MG_A_IN","U118MG_B_IN","U118MG_A_OUT","U118MG_B_OUT",
                 "UK8712_A_IN","UK8712_B_IN","UK8712_A_OUT","UK8712_B_OUT",
                 "UK7525_A_IN","UK7525_B_IN","UK7525_A_OUT","UK7525_B_OUT",
                 "KU50_A_IN","KU50_B_IN","KU50_A_OUT","KU50_B_OUT",
                 "KU7525_A_IN","KU7525_B_IN","KU7525_A_OUT","KU7525_B_OUT",
                 "KU8712_A_IN","KU8712_B_IN","KU8712_A_OUT","KU8712_B_OUT",
                 "KARPAS299_A_IN","KARPAS299_B_IN","KARPAS299_A_OUT","KARPAS299_B_OUT")]

#extract alk expression
KARPAS=kingrex_InOut["ALK",]
#select dilutions and order from 100%karpas299 to 0% and IN-OUT
KARPAS=KARPAS[,c("KARPAS299_A_IN","KARPAS299_B_IN","KARPAS299_A_OUT","KARPAS299_B_OUT",
                 "KU8712_A_IN","KU8712_B_IN","KU8712_A_OUT","KU8712_B_OUT",
                 "KU7525_A_IN","KU7525_B_IN","KU7525_A_OUT","KU7525_B_OUT",
                 "KU50_A_IN","KU50_B_IN","KU50_A_OUT","KU50_B_OUT",
                 "UK7525_A_IN","UK7525_B_IN","UK7525_A_OUT","UK7525_B_OUT",
                 "UK8712_A_IN","UK8712_B_IN","UK8712_A_OUT","UK8712_B_OUT",
                 "U118MG_A_IN","U118MG_B_IN","U118MG_A_OUT","U118MG_B_OUT"
)]


####fig 4 
#divide assay in from assay out
IN=KARPAS[1,stri_detect_fixed(colnames(KARPAS),"_IN")]
colnames(IN)=stri_replace_all_fixed(colnames(IN),"_IN","")
OUT=KARPAS[1,stri_detect_fixed(colnames(KARPAS),"_OUT")]
colnames(OUT)=stri_replace_all_fixed(colnames(IN),"_OUT","")
KARPAS2=rbind(IN,OUT)
rownames(KARPAS2)=c("IN","OUT")
#calculate mean fc for each percentage
FC=IN-OUT
FC2=c(mean(as.numeric(FC[1:2])),
      mean(as.numeric(FC[3:4])),
      mean(as.numeric(FC[5:6])),
      mean(as.numeric(FC[7:8])),
      mean(as.numeric(FC[9:10])),
      mean(as.numeric(FC[11:12])),
      mean(as.numeric(FC[13:14])))

#rename columns according to kaprpas299 percentage
colnames(KARPAS2)=c("100%_A","100%_B","87.5%_A","87.5%_B","75%_A","75%_B","50%_A","50%_B","25%_A","25%_B","12.5%_A","12.5%_B","0%_A","0%_B")
KARPAS2=as.data.frame(t(KARPAS2))
KARPAS2[,"Perc"]=c(rep("100%",2),rep("87.5%",2),rep("75%",2),rep("50%",2),rep("25%",2),rep("12.5%",2),rep("0%",2))
#add duplicate info A or B
KARPAS2[,"Dup"]=rep(c("A","B"),7)

#format data for ggplot2
melt_dataset=melt(KARPAS2)
colnames(melt_dataset)=c("Percentage","Dup","Assay","Log2(NC)")
melt_dataset$Percentage=factor(x =melt_dataset$Percentage,levels = c("100%","87.5%","75%","50%","25%","12.5%","0%"),ordered =T )
#format plot 
ggbase1=  geom_bar(position="dodge",stat="identity",width=0.7) 
ggbase2=scale_fill_manual(name="Assay",labels=c("IN","OUT"),values=c("royalblue","royalblue4"))
ggbase5=scale_color_manual(name="Duplicate",labels=c("A","B"),values=c("gold","red"))
ggbase3= theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
               legend.position="right", 
               text=element_text(size=38),
               plot.title = element_text(hjust = 0.5)#,
               #axis.title.x=element_text(size=0)
)
ggbase4=coord_cartesian(ylim=c(0,20))

#Fig. 4B
gg3<-ggplot(melt_dataset, aes(x=Percentage, y=`Log2(NC)`, fill=Assay,col=Dup)) + 
  ggbase1 + ggbase2 + ggbase3 + ggbase4 + ggbase5+
  ggtitle("ALK") 

#add fold change text
gg3= gg3 + annotate("text",size=8, x=(1:7), y=rep(15,7),
                    label=as.character(round(FC2,1)),
                    colour="blue")

gg3= gg3 + annotate("text",size=8, x=(1:7), y=rep(17,7),
                    label=as.character(rep("FC",7)),
                    colour="blue")

#save plot
jpeg(paste(path_output,"Fig.4B_ALK_gene_expression.jpeg",sep="/"),width = 700,height = 500)
print(gg3)
dev.off()

#divide assay in from assay out
IN=U118MG[1,stri_detect_fixed(colnames(U118MG),"_IN")]
colnames(IN)=stri_replace_all_fixed(colnames(IN),"_IN","")
OUT=U118MG[1,stri_detect_fixed(colnames(U118MG),"_OUT")]
colnames(OUT)=stri_replace_all_fixed(colnames(IN),"_OUT","")
U118MG2=rbind(IN,OUT)
rownames(U118MG2)=c("IN","OUT")
#calculate mean fc for each percentage
FC=IN-OUT
FC2=c(mean(as.numeric(FC[1:2])),
      mean(as.numeric(FC[3:4])),
      mean(as.numeric(FC[5:6])),
      mean(as.numeric(FC[7:8])),
      mean(as.numeric(FC[9:10])),
      mean(as.numeric(FC[11:12])),
      mean(as.numeric(FC[13:14])))
#rename columns according to u118mg percentage
colnames(U118MG2)=c("100%_A","100%_B","87.5%_A","87.5%_B","75%_A","75%_B","50%_A","50%_B","25%_A","25%_B","12.5%_A","12.5%_B","0%_A","0%_B")
U118MG2=as.data.frame(t(U118MG2))
U118MG2[,"Perc"]=c(rep("100%",2),rep("87.5%",2),rep("75%",2),rep("50%",2),rep("25%",2),rep("12.5%",2),rep("0%",2))
#add duplicate info A or B
U118MG2[,"Dup"]=rep(c("A","B"),7)

#format data for ggplot2
melt_dataset=melt(U118MG2)
colnames(melt_dataset)=c("Percentage","Dup","Assay","Log2(NC)")
melt_dataset$Percentage=factor(x =melt_dataset$Percentage,levels = c("100%","87.5%","75%","50%","25%","12.5%","0%"),ordered =T )

#Fig4.A
gg3<-ggplot(melt_dataset, aes(x=Percentage, y=`Log2(NC)`, fill=Assay,col=Dup)) + 
  ggbase1 + ggbase2 + ggbase3 + ggbase4 + ggbase5+
  ggtitle("ROS1") 

#add fold change text
gg3= gg3 + annotate("text",size=8, x=(1:7), y=rep(15,7),
                    label=as.character(round(FC2,1)),
                    colour="blue")

gg3= gg3 + annotate("text",size=8, x=(1:7), y=rep(17,7),
                    label=as.character(rep("FC",7)),
                    colour="blue")

#save plot
jpeg(paste(path_output,"Fig.4A_ROS1_gene_expression.jpeg",sep="/"),width = 700,height = 500)
print(gg3)
dev.off()

