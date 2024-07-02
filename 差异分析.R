##########step1:下载数据########
library(GEOquery)
setwd("C:/R/workplace/AD/GSE111006")
gse="GSE111006"
eSet <- getGEO(gse,
               destdir = '.',
               getGPL = F)#获取gse 的相关信息，表达矩阵，注释信息，平台信息
exp <- exprs(eSet[[1]]) #获取表达矩阵
exp[1:4,1:4]
#exp <- read.csv(file = "C:\\R\\workplace\\AD\\GSE85426\\GSE85426_series_matrix.csv",row.names = 1)#直接导入
exp = log2(exp+1) #根据实际情况判断是否须有加上log2，如果太大了就需要

##########step2：获取样本临床信息#########
pd <- pData(eSet[[1]]) #获取样本信息，年龄、性别、患病状态
pd[1:4,1:4]
gpl <- eSet[[1]]@annotation #获取平台号，要从这里得到注释文件，得到基因symbol
p = identical(rownames(pd),colnames(exp)) #identical（）判断样本名称的顺序是否一致
#保存上面提取的数据到文件中
save(gse,exp,pd,gpl,file = "step2_output.Rdata")

############step3：生成分组信息##########
#清空当前环境下的所有变量
rm(list = ls())
#加载Step2中保存的数据
load("step2_output.Rdata")
#加载stringr包
library(stringr)
#查看表达矩阵的列名
table(colnames(exp)) #第一行或者说第一行的上一行
#建立分组
group_list = ifelse(str_detect(pd$sarcopenia,"1"),"test","normal")#选中pd对应的分类列，识别疾病关键词标记为test
group_list=factor(group_list,
                  levels = c("test","normal"))
#查看分组信息
group_list
table(group_list)
#获取注释信息
if(T){
  #a = getGEO(gpl,destdir = ".")
  #b = a@dataTable@table
  b = read.table('C:\\R\\workplace\\AD\\GPL81-57556.txt',sep = "\t",
                 comment.char = "#", stringsAsFactors = F,
                 header = T, fill = TRUE, quote = "") #直接导入下载好的
  colnames(b)
  ids2 = b[,c("ID","Gene.Symbol")]#对应的修改为对应列名,提取出ID和symbol两列
  colnames(ids2) = c("probe_id","symbol")
  ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),] #去除symbol是空的以及含有两个的
}
#保存上述提取的数据
save(group_list,ids2,exp,pd,file = "step3_output.Rdata")

##########Step4：主成分分析##########
#接下来，使用 FactoMineR、factoextra 包执行主成分分析 (PCA)，以确定治疗组和对照组之间是否存在显着分组。
#清空当前环境下的变量
rm(list = ls())
#加载数据
load("step2_output.Rdata")
load("step3_output.Rdata")
#对表达矩阵进行转置操作(使原先的行变为列，列变为行)并转换成数据框格式
dat=as.data.frame(t(exp))#每一行为一个样本，每一列为一个基因
#加载包
#install.packages("factoextra")
library(FactoMineR)
library(ggplot2)
library(factoextra)
#执行PCA
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point",
                         col.ind = group_list,
                         palette = c("#E7B800", "#00AFBB"),#颜色
                                      addEllipses = TRUE,#添加椭圆
                                      legend.title = "Groups" #图例名
                         )
#查看PCA图
print(pca_plot)
#保存图片和数据
ggsave(plot = pca_plot, filename = paste0(gse,"PCA.png"))
save(pca_plot, file = "pca_plot. Rdata")

#######Step5:差异基因表达分析##########
#清除环境变量并加载Step2和3中保存的数据
rm(list = ls())
load("step2_output.Rdata")
load("step3_output.Rdata")
#加载limma包
library(limma)
#根据 group_list 创建一个设计矩阵，用于差异比较的设计
design=model.matrix(~group_list)
#将表达矩阵 exp 和设计矩阵 design 进行拟合
fit=lmFit(exp, design)
#对拟合对象'fit'进行贝叶斯估计
fit=eBayes(fit)
#从拟合对象 fit 中获取差异基因的结果
deg=topTable(fit, coef=2, number = Inf)
#查看差异基因结果的前几行
head(deg)
#向deg数据框添加几列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg)) 
#将deg数据框中的行名作为新的列 probe_id 添加到 deg
#tibble::rownames_to_column(deg) #和上行代码等价
head(deg)
#colnames(ids2)[1]="ID" #因为前面已经将第一列命名为probe_id所以要先改掉
#ids2 <- ids2 %>%
#    mutate(probe_id = symbol) #复制symbol列并命名为probe_id
#合并表
table(deg$probe_id %in% ids2$probe_id) 
#比较deg数据框中的probe_id列和ids数据框中的probe_id列的匹配情况，并生成一个计数表，显示匹配和不匹配的数量。
#deg <- inner_join(deg,ids,by="probe_id") #和上行代码等价
deg <- merge(x = deg,y = ids2, by="probe_id")
#将deg数据框和ids2数据框按照probe_id列进行合并
deg <-deg[!duplicated(deg$symbol),] 
#找到deg数据框中的重复行，并使用逻辑索引!duplicated(deg$symbol) 来删除重复行
dim(deg) #显示deg 数据框的维度
head(deg)
#为后续火山图绘制增加一列(上调或下调)
logFC_t=0.5  #如果结果不好可以调整一下1或者0.5之类的，参考其他文献
PValue=0.05
#logFC_t=mean(deg$logFC)+2*sd(deg$logFC) #和上行代码作用一样，都是设置logFC 的阈值
change=ifelse(deg$adj.P.Val > PValue,'stable',  #如果数据结果不好可以考虑用P值直接做
              ifelse( deg$logFC >logFC_t,'up',
                      ifelse( deg$logFC < -logFC_t,'down','stable') )
) #根据条件判断，如果P.Value大于0.05,则赋值为'stable';如果logFC大于logFC_t，则赋值为'up'；如果logFC小于-logFC_t，则赋值为'down'；否则赋值为'stable'
deg <- mutate(deg, change) #将change列添加到deg数据框中。
head(deg) #显示deg的前几行
table(deg$change) #根据 deg 数据框的change列中的值创建频数表
#添加 ENTREZID 列，后续将在富集分析中使用
library(ggplot2) #载入ggplot2软件包，用于数据可视化
#if (!require("BiocManager", quietly = TRUE))
#install.packages("clusterProfiler")
#BiocManager::install("clusterProfiler")
library(clusterProfiler) #载入clusterProfiler软件包，用于生物学注释和富集分析
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db) #载入org.Hs.eg.db软件包，该软件包包含了人类基因组的注释信息
s2e <- bitr(unique(deg$symbol), fromType = "SYMBOL",
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db) 
#将deg数据框中的symbol列的基因符号（SYMBOL）转换为对应的Entrez ID
head(s2e) #打印 s2e 对象的前几行，以查看基因符号到Entrez ID的转换结果
head(deg) #打印 deg 数据框的前几行，以查看处理前的结果
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL")) 
#将deg数据框和s2e对象根据symbol列和SYMBOL列进行内连接（内部匹配），即将具有相同基因符号的行合并
head(deg)#查看基因注释后的前几个结果
dim(deg)
save(logFC_t,deg,file = "step5_output.Rdata") #保存数据，以便后续的分析使用

##########Step6：可视化火山图和热图#########
#通常可以通过火山图和热图来可视化差异差异基因表达数据
rm(list = ls())
load("step2_output.Rdata")
load("step3_output.Rdata")
load("step5_output.Rdata")
library(dplyr)
#绘制火山图
dat <-mutate(deg,v=-log10(P.Value))
head(dat)
if(T){
  for_label <- dat %>%
    filter(symbol %in% c("RUNX2","FN1"))
}
if(F){
  for_label <- dat %>% head(5)
}
if(T) {
  x1 = dat %>%
    filter(change == "up") %>%
    head(3)
  x2 = dat %>%
    filter(change == "down") %>%
    head(3)
  for_label = rbind(x1,x2)
}
p <- ggplot(data = dat,
            aes(x = logFC,
                y = v)) +
  geom_point(alpha=0.4, size=3.5,
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) + #geom_vline添加垂线
  geom_hline(yintercept = -log10(0.05), lty=4, col="black", lwd=0.8) + #geom_hline添加水平线
  theme_bw()
p
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )#添加for_label包含的基因名字
volcano_plot
ggsave(plot = volcano_plot, filename = paste0(gse,"volcano.png"))
#绘制热图
cg=names(tail(sort(apply(exp,1,sd)),1000))
#这行代码是为了筛选前1000的基因
#apply(X数据对象，1表示行/2表示列，FUN函数(sd标准差))
#sort()排序，排序结果不可逆转,默认是升序.
#tail() 函数用于获取向量,矩阵,表,DataFrame 或函数的最后部分,tail(x, n)x:指定的数据类型,n:需要打印的行数
#names()R语言中的函数用于获取或设置对象的名称
n=exp[cg,] #n为exp中的前1000的基因数据
annotation_col=data.frame(group=group_list)
rownames(annotation_col) = colnames(n)
library(pheatmap)
heatmap_plot <- pheatmap(n,#表达矩阵
                         show_colnames=F,#是否展示列名
                         show_rownames = F,#是否展示行名
                         annotation_col = annotation_col, #列的分组信息
                         cutree_rows=2,cutree_cols=2,## cutree_rows, cutree_cols可以根据行列的聚类数将热图分隔开
                          scale = "row")#表示值均一化的方向，或者按照行或列，或者没有，值可以是"row", “column” 或者"none"
#保存结果
library(ggplot2)
png(file = paste0(gse,"heatmap.png"))
ggsave(plot = heatmap_plot, filename = paste0(gse,"heatmap.png"))
dev.off()
write.csv(deg,"deg.csv")
################step7:GO和KEGG富集分析#############
#用到的文档就是画火山图的数据
rm(list = ls())
load("step5_output.Rdata")
dim(deg)
deg <-deg[!duplicated(deg$symbol),] 
dim(deg)
rownames(deg) = deg$symbol
#.libPaths()
#.libPaths("C:/Users/lxt/AppData/Local/R/win-library/4.3")#永久定义R包安装路径
install.packages("Rgraphviz")
BiocManager::install("topGO")
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(ggplot2)
library(stringr)
#数据处理-差异基因筛选#
vector = abs(deg$logFC) > 1.5 & deg$P.Value < 0.05 & deg$symbol !="" 
#abs绝对值;通常logFC> 1和PValue< 0.05条件进行筛选；deg$gene_name != ""表示gene_name不为空白
#deg$gene_name<-str_to_title(deg$gene_name)#用stringr将基因名称的第一个字母大写（小鼠首字母为大写）
data_sgni= deg[vector,]#筛选差异基因
head(data_sgni)
dim(data_sgni)
#All_gene <- rownames(deg) # 提取所有基因基因名
#已有OrgDb的常见物种#
library(org.Hs.eg.db)
#基因ID转换#
keytypes(org.Hs.eg.db) #查看所有可转化类型
entrezid_all = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
                      keys = data_sgni$symbol, #将输入的gene_name列进行数据转换
                      keytype = "SYMBOL", #输入数据的类型
                      column = "ENTREZID")#输出数据的类型
entrezid_all  = na.omit(entrezid_all)  #na省略entrezid_all中不是一一对应的数据情况
entrezid_all = data.frame(entrezid_all) #将entrezid_all变成数据框格式
head(entrezid_all)
###GO富集分析###
GO_enrich = enrichGO(gene = entrezid_all[,1], #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", #输入数据的类型
                     ont = "ALL", #可以输入CC/MF/BP/ALL
                     #universe = 背景数据集 # 表示背景基因，无参的物种选择组装出来的全部unigenes作为背景基因；有参背景基因则不需要。
                     pvalueCutoff = 1,qvalueCutoff = 1, #表示筛选的阈值，阈值设置太严格可导致筛选不到基因。可指定 1 以输出全部
                     readable = T) #是否将基因ID映射到基因名称。
GO_enrich  = data.frame(GO_enrich) #将GO_enrich导成数据框格式

table(GO_enrich$ONTOLOGY)
#数据保存#
###KEGG富集分析###
install.packages("R.utils")
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], #即待富集的基因列表
                         keyType = "kegg",
                         pAdjustMethod = 'fdr',  #指定p值校正方法
                         organism= "human",  #hsa，可根据你自己要研究的物种更改，可在https://www.kegg.jp/brite/br08611中寻找
                         qvalueCutoff = 1, #指定 p 值阈值（可指定 1 以输出全部）
                         pvalueCutoff=1) #指定 q 值阈值（可指定 1 以输出全部）
KEGG_enrich  = data.frame(KEGG_enrich)
#输出的GO/KEGG富集结果各列内容：
#ONTOLOGY：GO的BP（生物学过程）、CC（细胞组分）或MF（分子功能）三个方面内容；
#ID：富集到的GO term/KEGG term；
#Description：对GO term/KEGG term的生物学功能和意义进行描述；
#GeneRatio：富集到该GO term/KEGG term中的基因数目/给定基因的总数目；
#BgRatio：该GO term/KEGG term中背景基因总数目/该物种所有已知GO功能基因的数目；
#pvalue、p.adjust和qvalue：p值、校正后p值和q值；
#geneID和Count：富集到该GO term/KEGG term中的基因名称和数目。
#write.csv(GO_enrich,'GO_enrich.csv') #数据导出
#保存富集分析的结果
save(GO_enrich,KEGG_enrich,entrezid_all,file = "step7_output.Rdata")
################step8:可视化GO和KEGG富集结果#############
#数据载入与处理#
rm(list = ls())
load("step2_output.Rdata")
load("step7_output.Rdata")
library(ggplot2)
#library(openxlsx)
#go_enrich = read.csv("GO_enrich.csv",sheet= "ONTOLOGY",sep=',')  
go_enrich <- GO_enrich
go_enrich$term <- paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term,ordered = T)
#go_enrich <- go_enrich[order(go_enrich$Count, decreasing = T),]
table(go_enrich$ONTOLOGY)
#分开三个goterm:BP,CC,MF
go_enrich_BP <- as.data.frame(go_enrich[go_enrich$ONTOLOGY== "BP",])
go_enrich_CC <- as.data.frame(go_enrich[go_enrich$ONTOLOGY== "CC",])
go_enrich_MF <- as.data.frame(go_enrich[go_enrich$ONTOLOGY== "MF",])
#针对count数排序
go_enrich_BP <- go_enrich_BP[order(go_enrich_BP$Count, decreasing = T),]
go_enrich_CC <- go_enrich_CC[order(go_enrich_CC$Count, decreasing = T),]
go_enrich_MF <- go_enrich_MF[order(go_enrich_MF$Count, decreasing = T),]
#选取展示的个数
display_number = c(15, 10, 5)  ##这三个数字分别代表选取的BP、CC、MF的数量
go_enrich_BP = as.data.frame(go_enrich_BP)[1:display_number[1], ]
go_enrich_CC = as.data.frame(go_enrich_CC)[1:display_number[2], ]
go_enrich_MF = as.data.frame(go_enrich_MF)[1:display_number[3], ]
#将提取的各组数据进行整合
go_enrich2 <- rbind(go_enrich_BP, go_enrich_CC)
go_enrich2 <- rbind(go_enrich2, go_enrich_MF)
go_enrich2$term <- factor(go_enrich2$Description,levels=go_enrich2$Description,ordered = T)
head(go_enrich2)
#横向柱状图#
goplotR <- ggplot(go_enrich2, 
       aes(x=term,y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  xlab("GO term") + #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
goplotR
ggsave(plot = goplotR, 
       filename = paste0(gse,"goR.png"),
       limitsize = FALSE,
       width = 10,             # 宽
       height = 10,            # 高
       units = "in",          # 单位
       dpi = 300)
#纵向柱状图#
goplotC <- ggplot(go_enrich2, 
       aes(x=term,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()
#help(theme) #查阅这个函数其他具体格式
goplotC
ggsave(plot = goplotC, 
       filename = paste0(gse,"goC.png"),
       limitsize = FALSE,
       width = 10,             # 宽
       height = 10,            # 高
       units = "in",          # 单位
       dpi = 300)
#气泡图#
goplotqipao <- ggplot(go_enrich2,
       aes(y=term,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()
goplotqipao
ggsave(plot = goplotqipao, 
       filename = paste0(gse,"goqipao.png"),
       limitsize = FALSE,
       width = 10,             # 宽
       height = 10,            # 高
       units = "in",          # 单位
       dpi = 300)
###KEGG可视化###
#数据导入#
kk_result <- KEGG_enrich
#排序
kk_result <- kk_result[order(kk_result$Count, decreasing = T),]
#数据处理#
display_number = 30#显示数量设置
kk_result = as.data.frame(kk_result)[1:display_number[1], ]
kk = as.data.frame(kk_result)
rownames(kk) = 1:nrow(kk)
kk$order=factor(rev(as.integer(rownames(kk))),labels = rev(kk$Description))
#柱状图#
keggzhu <-  ggplot(kk,aes(y=order,x=Count,fill=pvalue))+
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "KEGG Pathways Enrichment",  #设置标题、x轴和Y轴名称
       x = "Gene number",
       y = "Pathway")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
keggzhu
ggsave(plot = keggzhu, 
       filename = paste0(gse,"keggzhu.png"),
       limitsize = FALSE,
       width = 10,             # 宽
       height = 10,            # 高
       units = "in",          # 单位
       dpi = 300)
#气泡图#
keggqipao <- ggplot(kk,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"),
       x="Gene Ratio",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
keggqipao
ggsave(plot = keggqipao, 
       filename = paste0(gse,"keggqipao.png"),
       limitsize = FALSE,
       width = 10,             # 宽
       height = 10,            # 高
       units = "in",          # 单位
       dpi = 300)