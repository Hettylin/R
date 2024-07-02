#library(GEOquery)
#library(dplyr)
#library(tidyverse)
#library (data.table)
#library(limma)
#library(DESeq2)
setwd("C:/R/workplace/AD")
Sys.setlocale('LC_ALL','C')
#获取GPL文件
GPL_table = read.table('C:\\R\\workplace\\AD\\GPL13534-11288.txt',sep = "\t",
                       comment.char = "#", stringsAsFactors = F,
                       header = T, fill = TRUE, quote = "")  #下载好的数据
#获取GSE文件
#gset <- getGEO("GSE48350",destdir = ".",AnnotGPL = F,getGPL = F)
#gset[["GSE48350_series_matrix.txt.gz"]]@annotation
#gset[[1]]
#exp<-exprs(gset[[1]]) #样本数据
#cli<-pData(gset[[1]])	## 临床数据
gse <- read.table('C:\\R\\workplace\\AD\\GSE111006\\GSE111006_series_matrix.txt',sep = "\t",
                      comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)  #下载好的数据
ID_Sybmol = GPL_table[,c(1,22)] #获取DI和Sybmol
colnames(ID_Sybmol)[2]="Symbol" #改列名
#c(1,11)取1和11列数据，根据需求修改
#合并ID
Exp = merge(ID_Sybmol,gse,by.x = "ID",by.y = "ID_REF",all=T) #将ID_Sybmol的ID合并到gse 
Exp = Exp[,-1]
#Symbol中存在///的数据，保留第一个
Exp$`Symbol`<-data.frame(sapply(Exp$`Symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
#数据过滤，对多余的基因进行删除
#Exp = Exp[Exp$Symbol != "",] #对空白赋值
Exp = na.omit(Exp)   #删除NA数据
#table(duplicated(Exp$`Symbol`)) #查看有多少重复值
#EXP1<-avereps(Exp[,-c(1,ncol(Exp))],ID=Exp$`Symbol`)#去重同时转标准化矩阵
Exp1 <-Exp[!duplicated(Exp$Symbol),]  #去重
#write.table(Exp1,"GSE109887.txt",row.names = F,quote = F,sep="\t")
rownames(Exp1)<-Exp1[,1] #将数据框的第一列作为行名
Exp1<-Exp1[,-1] #将数据框的第一列删除，只留下剩余的列作为数据
#Exp1 <- log2(Exp1+1)#对于没有标准化处理过的数据需要做log2转换,有的gene表达量为0，会影响对数转换，所以对所有数据+1，保证结果不会报错
write.csv(Exp1, "C:\\R\\workplace\\AD\\GSE99624.csv ")



#####另一种处理ID合并的方法
setwd("C:\\R\\workplace\\AD")	
options(stringAsFactors = F)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*600)	# 设置内存量
library(GEOquery)
library(limma)
library(affy)
#解除60s限时
#library(SeuratData)  
#getOption('timeout')
#options(timeout=10000)
#提高下载速度
#options( 'download.file.method.GEOquery' = 'libcurl' )
gset <- getGEO('GSE111006', destdir=".",
                AnnotGPL = T,     ## 注释数据
                getGPL = T)       ## 临床数据
gset <- getGEO("GSE48350",destdir = ".",AnnotGPL = F,getGPL = T)
gset[["GSE48350_series_matrix.txt.gz"]]@annotation
gset[[1]]
exp<-exprs(gset[[1]]) #表达矩阵
cli<-pData(gset[[1]])	#### 获取临床信息
group<-c(rep("control",3),rep("hht",3))	## 查看分组信息，根据需求修改
#用下面的方式可以避免超时问题，但需要下载好数据压缩包
#gset <- getGEO("GSE48350",destdir = ".",AnnotGPL = F,getGPL = F)
#gset[["GSE48350_series_matrix.txt.gz"]]@annotation
#gset[[1]]
#exp<-exprs(gset[[1]]) 
#cli<-pData(gset[[1]])	
gpl <- data.table::fread("GPL570-55999.txt", sep = "\t", header = TRUE)
#GPL <- gpl[gpl$`Gene Symbol`!= '', ]##???????????????
GPL<-fData(gset[[1]])	## 获取平台数据
gpl<-GPL[,c(1,3)]  #取1和3列，根据需要修改
##有两个Gene symbol与之对应并用“///”分隔开，只取第一个Gene symbol来达到去重效果
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
exp<-as.data.frame(exp) #转化为数据矩阵
exp$ID<-rownames(exp)	# 提取列名
exp_symbol<-merge(exp,gpl,by="ID")#合并
exp_symbol<-na.omit(exp_symbol)#删除空值
table(duplicated(exp_symbol$`Gene symbol`)) #查看重复量
exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`Gene symbol`)#去重同时转化为标准化矩阵
write.table(exp_unique,"C:\\R\\workplace\\AD\\GSE48350.txt")

###自动ID转换
library(org.Hs.eg.db)
#基因ID转换#
keytypes(org.Hs.eg.db) #查看所有可转化类型
darkgrey <-read.csv(file = "GSE111006_allSamplesCounts_htseqcov1_hss_forGEO.csv",row.names = 1)
darkgrey$X<-rownames(darkgrey)
symb = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
              keys = darkgrey$X, #将输入的gene_name列进行数据转换
              keytype = "ENSEMBL", #输入数据的类型
              column = "SYMBOL")#输出数据的类型
symb  = na.omit(symb)  #na省略entrezid_all中不是一一对应的数据情况
symb = data.frame(symb) #将entrezid_all变成数据框格式
symb <-symb[!duplicated(symb$symb),]#去重
symb$X<-rownames(symb)
data1 <- merge(darkgrey,symb,by= "X")
row.names(data1)<-data1$symb
data1<-data1[,-41]
write.csv(data1,"GSE111006_idchange.csv")
getwd()
##中位数筛选
library(dplyr)
exp2 <- mutate(exp,probe_id=rownames(exp))
table(exp2$probe_id %in% ids2$probe_id) 
exp2<- merge(x = exp2,y = ids2, by="probe_id")
exp2=exp2[order(exp2$symbol,exp2$median,decreasing = T),]
exp2 <- exp2[!duplicated(exp2$symbol),] 
rownames(exp2) <- exp2$symbol
write.csv(exp2,"symbol_GSEq12288.csv")
