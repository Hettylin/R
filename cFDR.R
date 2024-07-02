#install.packages("randomForest")
library(vcfR) 
library(randomForest)
###IEU数据库数据处理：VCF文件的转换，转成CSV
setwd("C:\\R\\workplace\\cFDR")
vcf<-read.vcfR("ebi-a-GCST006368.vcf") #读数据
##获取基因信息
gt<-vcf@gt
gt<-as.data.frame(gt)#转格式
colnames(gt)<-c("aa","bb")#重命名
gt = subset(gt,gt$aa == "ES:SE:LP:AF:ID")#将有缺失的值删除，如果不是五个数据都有下面的步骤会报错
BBMI<-gsub(c(":",":"),",", gt$bb)#分割数据，转成按列呈现
BBBMI<-read.table(text = BBMI,sep=",",header = F)
colnames(BBBMI)<-c("beta","se","pval","eaf","rsids")#列命名
BBBMI$pval<-10^-(BBBMI$pval)#P值的转换，将不是P值的数据转为P值（对应log的情况）（LP不是P值）
gt1 = subset(BBBMI, select = c(rsids,pval))#提取列
write.csv(gt1,"ebi-a-GCST006368.csv")
#获取前八列数据中的位置信息
fix<-vcf@fix
fix1 <-subset(fix, select = c(CHROM, POS,ID))#提取列
fix1<-as.data.frame(fix1)#转格式
colnames(fix1)<-c("chr","pos","rsids")#改列名
data1<-data.frame(gt1,fix1)#合并数据
data1<-data5[,-1]
write.csv(data1,"ebi-a-GCST90000025_1.csv")
#芬兰数据库数据处理
library(data.table)
data<-fread('finngen_R5_F5_ALZHDEMENT',head=T)
#colnames(data)[1]="chr"#改列名
data2 = subset(data, select = c(rsids,pval,chr,pos))
data2 <- as.data.frame(data2)#转格式
#write.csv(data2,"finn-b-F5_ALZHDEMENT_1.csv")

###找共同的SNP数据
library(data.table)
library(dplyr)
data1<-read.csv("ebi-a-GCST005348.csv")#读取处理好的数据
data1<-data1[,-1]#删掉读取csv后产生的第一列序列号，直接从上面获取的不需要
data2<-read.csv("prot-a-114.csv")
data2<-data2[,-1]
#data3<-intersect(data1$rsids,data2$rsids)#两个数据的rsids列取交集（需要values形式而不是数据矩阵）
data3 <- merge(data1,data2,by= "rsids")
#data3 <- data3[,c("rsids","pval.x","pval.y",)]
write.csv(data3,"ebi-a-GCST006368&prot-a-114.csv")

#解除60s限时
#getOption('timeout')
#options(timeout=10000)
###LD分析
#colnames(data1_)<-c("x","chr.exposure","pos.exposure","SNP","pval.exposure")#列命名
#colnames(data2_)<-c("x","SNP","chr.exposure","pos.exposure","pval.exposure")#列命名
#data1_3 <-clump_data(data1_,clump_r2=0.001,clump_kb=10000)
#data1_3 <-clump_data(data2_,clump_r2=0.2,clump_kb=200)
library(TwoSampleMR)
#if(!require("pacman"))install.packages("pacman",update = F,ask = F)
library("pacman")
#install.packages("yulab-smu/yulab.uti1s")
library(yulab.utils)
p_load(data.table,dplyr,tidyr)##对包进行加载安装
library("ieugwasr")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
data3 <- read.csv(file = "ukb-a-248&prot-a-114.csv",row.names = 1)
#data3<-data3[,-1]#删掉读取csv后产生的第一列序列号，直接从上面获取的不需要
data3_=data3[,c("rsids","pval.x")]
colnames(data3_)<-c("rsid","pval")#列命名
#dat2_1 <- subset(dat2_, pval <= 1)#删除p值大于1的（不知道为什么但不这样会报错）；也有可能是有空值或者SNP检索不到；后面发现是UKB数据库的LP没有转换成P值，转换之后就不会有问题
expo_clump <- ieugwasr::ld_clump_local(dat = data3_,
                                       clump_r2 = 0.2,
                                       clump_p = 1,
                                       clump_kb = 200,
                                       plink_bin = "C:\\plink\\plink",
                                       bfile = "C:\\R\\workplace\\cFDR\\data_maf0.01_rs_ref\\data_maf0.01_rs_ref"
)#bfile上写的是之前下载的data_maf0.01_rs_ref的工作路径，设置完后再加一个data_maf0.01_rs_ref（代表文件夹里的三个文件）。同理，plink_bin上写的之前下载的plink工作路径，再在后面加上plink
#pink link软件：https://www.cog-genomics.org/plink/1.9/
#data_maf0.01_rs_ref下载地址：http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz
df <- merge(expo_clump,data3,by.x= "rsid",by.y= "rsids")
df <- df[-2]#删除无用列
write.csv(df,"ebi-a-GCST005348&prot-a-114_1.csv")

#P1和P2分别对应两个表型的P值
#install.packages("devtools")
#devtools::install_github("KehaoWu/GWAScFDR")
library(GWAScFDR) 
#cf1=cfdr::cfdr(df[,p1],df[,p2])
#cf2=cfdr::cfdr(df[,p2],df[,p1])
#cf1=GWAScFDR::cFDR(p1,p2)
#cf2=GWAScFDR::cFDR(p2,p1)
cf1=cFDR(df[,"pval.x"],df[,"pval.y"])
cf2=cFDR(df[,"pval.y"],df[,"pval.x"])
cf=data.frame(cFDR1=cf1,cFDR2=cf2)
cf$ccFDR=apply(cf, 1, max)
#ccf=condFDR::ccFDR(df,p1,p2,p_threshold  = 0.1,mc.cores = 8)
res=data.frame(df,cf)
#res=plyr::rename(res,c("cFDR1"=paste0(toupper(t[1]),"|",toupper(t[2]),".cFDR"),"cFDR2"=paste0(toupper(t[2]),"|",toupper(t[1]),".cFDR"),"ccFDR"=paste0(toupper(t[1]),"|",toupper(t[2]),".ccFDR")))
write.csv(res,"ebi-a-GCST005348&prot-a-114_2.csv")
#画图
#install.packages("qqman")
library(qqman)           #加载qqman包
res<-read.csv("ukb-a-379&finn-b-F5_ALZHDEMENT_2.csv")
P = manhattan(res,p = "ccFDR",bp = "pos",chr = "chr",snp =  "rsid",#读取数据
              main="manhattan Plot",#标题命名
              xlab = "Chromosomal Location",#横坐标命名
              ylab = "-log10(ccFDR)",#纵坐标命名
              suggestiveline = -log10(5e-02),#－log10(5e－5)处添加"suggestive"横线
              #genomewideline = -log10(5e-08),#－log10(5e－8)处添加"genome-wide sigificant"横线
              #annotatePval = 0.05,#标记p值小于0.05的点
              #annotateTop = F,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
              #highlight = snpsOfInterest,#内置高亮的snp数据， 也可以对snpOfInterest进行设置
              #col = c("blue4", "orange3","red"), #设置散点的颜色，默认是黑灰配色
              ylim = c(0, 5), #设置y轴范围
              #cex = 0.6, #设置点的大小
              #cex.axis = 0.9, #设置坐标轴字体大小
              #suggestiveline = F, genomewideline = F, #移除标记线
              #chrlabs = c(paste0("chr",c(1:20)),"P","Q") #设置x轴染色体标签名
              )
#注意chr列中不能有X或Y染色体，可尝试转化为23或将其删除，并将该列转换为int形式 
#res1<-res1[-440393,] 
#res1$chr<-as.integer(res1$chr)

qq(res$pval.x, #选择数据
   main = "Q-Q plot", #标题命名
   xlab = "Empirical -log (q[BMD|Aβ])",ylab = "Nominal -log p[BMD]", #轴命名
   #xlim = c(0, 7), ylim = c(0, 12), #设置x,y轴范围
   pch = 18, col = "skyblue", cex = 1, las = 1)    #qq图只需要一列p值的数据
