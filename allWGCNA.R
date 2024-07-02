rm(list = ls())
library(WGCNA)  #加载WGCNA包
options(stringsAsFactors = FALSE)  #开启多线程
########################数据准备######################
setwd("C:/R/workplace/AD/GSE111006")
expr <- read.csv(file = "GSE111006.csv",row.names = 1)
dim(expr)
#femData <- expr[order(apply(expr,1,mad), decreasing = T)[1:5000],]# 取表达量前5000
datExpr0 = as.data.frame(t(expr))
#datExpr0就是一个以每行为样本，一列为一个基因的数据框
#b <- rownames(datExpr0[53,])
#c <- rownames(datExpr0[108,])#判断不需要的数据的位置
##去除离群值
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 90, col = "red") #划定需要剪切的枝长h
clust = cutreeStatic(sampleTree, cutHeight = 90, minSize = 1)#保留红线下的部分也就是clust==1#cutHeight跟h需要同步修改,如果样本量很小注意修改最小样本量
table(clust) 
keepSamples = (clust!=0)  #保留非离群(clust!=0)的样本
keepSamples = (clust==1) 
datExpr = datExpr0[keepSamples, ]  #去除离群值后的数据
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#构建样本分类矩阵，分为疾病和控制组
datTraits = read.csv(file = "GSE111006_datTraits.csv",row.names = 1)#导入构建好的分类矩阵
#datTraits <- data.frame(case=c(rep(1,25),rep(0,9)), control= c(rep(0,25),rep(1,9)))#根据原始数据建分类矩阵，1表示是对应的分组
rownames(datTraits) <- rownames(datExpr0)#将样本名作为行名
datTraits <- datTraits[rownames(datTraits) %in% rownames(datExpr), ]#比较行名，如果原来的样本在删减后的样本中则保留
#write.csv(datTraits,"GSE111006_datTraits.csv")#新构建的分类矩阵保存一下方便下次使用

####################网络构建##############
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
pdf("1Threshold.pdf",width = 9,height = 5)
par(mfrow = c(1,2))
cex1 = 0.9;
#无标度拓扑拟合指数
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))+
     text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")+
     abline(h=0.9,col="red")  #根据软阈值设置h,，可以改变高度值
#平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))+
     text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#运行下面的代码，如果有合适的软阈值，系统会自动推荐给你
sft$powerEstimate

#2.2 一步法构建网络和模块检测
net = blockwiseModules(datExpr, power = 5,
                       TOMType = "unsigned", minModuleSize = 150,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
# power = 6是刚才选择的软阈值
#minModuleSize：模块中最少的基因数
#mergeCutHeight ：模块合并阈值，阈值越大，模块越少（重要）
#saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM"保存TOM矩阵，名字为"femaleMouseTOM"
#net$colors 包含模块分配，net$MEs 包含模块的模块特征基因
#查看划分的模块数和每个模块里面包含的基因个数
table(net$colors);
#模块标识的层次聚类树状图，可以使用以下代码将树状图与颜色分配一起显示
#sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf("new_module.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#保存分配模块和模块包含的基因信息。
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
text <- unique(moduleColors)
for(i in 1:length(text)){
  y = t(assign(paste(text[i],"expr",sep = "."),
               datExpr[moduleColors == text[i]]))
  write.csv(y,paste(text[i],"csv",sep="."),quote = F)
}
####################3表型模块相关性分析###############
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#相关性
moduleTraitCor = cor(MEs, datTraits, use = "p");
#显著性 P值
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
textMatrix = paste(signif(moduleTraitCor, 2), 
                   "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("4module-trait.pdf",width = 5.2,height = 7)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits), cex.lab.x=0.4,
               yLabels = names(MEs), cex.lab.y=0.35,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(200),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("Module-trait relationships"))
dev.off()

########3.2 基因与表型数据的关系、重要模块：基因显著性和模块成员########
#基因都参与什么样的功能，可能进行功能富集
colnames(datTraits)[1] <-"case"
case = as.data.frame(datTraits$case);
names(case) = "case";
modNames = substring(names(MEs), 3)#从第三个字母开始取
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, case, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(case), sep="");
names(GSPvalue) = paste("p.GS.", names(case), sep="");

#########3.3 模块内分析：鉴定具有高GS和高MM的基因###########
module = "purple"  #选择显著的模块
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
##根据MM>0.8且GS>0.2筛选出核心基因
MM = abs(geneModuleMembership[moduleGenes, column])
GS = abs(geneTraitSignificance[moduleGenes, 1])
MMGS<-as.data.frame(cbind(MM,GS))
module<-read.csv(file="purple.csv",row.names = 1)#输入最显著模块,这里是项目自动保存的位置
rownames(MMGS)=rownames(module)
hub_b<-abs(MMGS$MM)>0.8&abs(MMGS$GS)>0.2
table(hub_b)
hub_a<-subset(MMGS, abs(MMGS$MM)>0.8&abs(MMGS$GS)>0.2)
write.csv(hub_a, "hubgene_MMGS.csv")#将筛选出来的hub基因保存文档
