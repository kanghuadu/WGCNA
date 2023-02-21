加权基因共表达网络分析 (WGCNA, Weighted correlation network analysis)是用来描述不同样品之间基因关联模式的系统生物学方法，可以用来鉴定高度协同变化的基因集, 并根据基因集的内连性和基因集与表型之间的关联鉴定候补生物标记基因或治疗靶点。

该分析方法旨在寻找协同表达的基因模块(module)，并探索基因网络与关注的表型之间的关联关系，以及网络中的核心基因。

适用于复杂的数据模式，推荐`5组`(或者`15个样品`)以上的数据。一般可应用的研究方向有：不同器官或组织类型发育调控、同一组织不同发育调控、非生物胁迫不同时间点应答、病原菌侵染后不同时间点应答。

----
**输出结果**
![6](https://files.mdnice.com/user/40177/99d2bd7e-5597-43e1-8300-ffc1c17db0e0.png)!
[7](https://files.mdnice.com/user/40177/1c47960f-7f27-4b45-b1b1-24f997bd1bc3.png)!
[10](https://files.mdnice.com/user/40177/181e82a1-2bef-4472-b7b4-63827b8f8753.png)

### 分析代码
**关于WGCNA分析，如果你的数据量较大，建议使用服务期直接分析，本地分析可能导致R崩掉。**
## 设置文件位置
```
setwd("~/00_WGCNA/20230217_WGCNA/WGCNA_01")
```
加载分析所需的安装包
```
install.packages("WGCNA")
#BiocManager::install('WGCNA')
library(WGCNA)
options(stringsAsFactors = FALSE)
```
注意，如果你想打开多线程分析，可以使用一下代码
```
enableWGCNAThreads() 
```
![取决于你的电脑线程数量](https://files.mdnice.com/user/40177/3758bf0d-d904-4aa4-a3e3-1f09ad050a95.png)
# 一、导入基因表达量数据
```
## 读取txt文件格式数据
WGCNA.fpkm = read.table("ExpData_WGCNA.txt",header=T,
                        comment.char = "",
                        check.names=F)
###############
# 读取csv文件格式
WGCNA.fpkm = read.csv("ExpData_WGCNA.csv", header = T, check.names = F)
```
![输入数据格式](https://files.mdnice.com/user/40177/42042c36-530b-4f75-af9f-4d233477a887.png)

## 数据处理
```
dim(WGCNA.fpkm)
names(WGCNA.fpkm)
datExpr0 = as.data.frame(t(WGCNA.fpkm[,-1]))
names(datExpr0) = WGCNA.fpkm$sample;##########如果第一行不是ID命名，就写成fpkm[,1]
rownames(datExpr0) = names(WGCNA.fpkm[,-1])
```
### 过滤数据
```
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
```
**过滤低于设定的值的基因**
```
##filter
meanFPKM=0.5  ###--过滤标准，可以修改
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]
# for meanFpkm in row n+1 and it must be above what you set--select meanFpkm>opt$meanFpkm(by rp)
filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
write.table(filtered_fpkm, file="mRNA.filter.txt",
            row.names=F, col.names=T,quote=FALSE,sep="\t")
```
### Sample cluster
```
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1.sampleClustering.pdf", width = 15, height = 8)
par(cex = 0.6)
par(mar = c(0,6,6,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2)
### Plot a line to show the cut
#abline(h = 180, col = "red")##剪切高度不确定，故无红线
dev.off()
```
-----
### 不过滤数据
**如果你的数据不进行过滤直接进行一下操作,此步与前面的操作相同，任选异种即可。**
```
## 不过滤
## Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr0 = datExpr0[keepSamples, ]
write.table(datExpr0, file="mRNA.symbol.uniq.filter.sample.txt",
            row.names=T, col.names=T,quote=FALSE,sep="\t")

###
#############Sample cluster###########
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1.sampleClustering.filter.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
### Plot a line to show the cut
#abline(h = 50000, col = "red")##剪切高度不确定，故无红线
dev.off()
```
---
# 二、导入性状数据
```
traitData = read.table("TraitData.txt",row.names=1,header=T,comment.char = "",check.names=F)
allTraits = traitData
dim(allTraits)
names(allTraits)
```
![](https://files.mdnice.com/user/40177/e2227a55-d511-40b5-9b05-522aac2493d2.png)
```
## 形成一个类似于表达数据的数据框架
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits)
collectGarbage()
```
![](https://files.mdnice.com/user/40177/1f4d924c-d38f-4ff8-a23f-3f90d3c0ceba.png)
### 再次样本聚类
```
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
```
输出样本聚类图
```
pdf(file="2.Sample_dendrogram_and_trait_heatmap.pdf",width=20,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2)
dev.off()
```
![](https://files.mdnice.com/user/40177/b21a415f-4876-4d77-bb26-b98989f51879.png)

# 三、WGCNA分析（后面都是重点）
## 筛选软阈值
```
enableWGCNAThreads()
# 设置soft-thresholding powers的数量
powers = c(1:30)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
```
此步骤是比较耗费时间的，静静等待即可。
![](https://files.mdnice.com/user/40177/b5fe7a93-f9f4-4f4f-bc8d-4d0041ade183.png)
### 绘制soft Threshold plot
```
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
```
![](https://files.mdnice.com/user/40177/b9dacde2-e24e-4eea-913e-a0c1103d98fe.png)
## 选择softpower
选择softpower是一个`玄学`的过程，可以直接使用软件自己认为是最好的softpower值，但是不一定你要获得最好结果；其次，我们自己选择自己认为比较好的softpower值，但是，需要自己不断的筛选。因此，从这里开始WGCNA的分析结果就开始受到不同的影响。
```
## 选择软件认为是最好的softpower值
#softPower =sft$powerEstimate
---
# 自己设定softpower值
softPower = 9
```
继续分析
```
adjacency = adjacency(datExpr0, power = softPower)
```
## 将邻接转化为拓扑重叠
这一步建议去服务器上跑，后面的步骤就在服务器上跑吧，数据量太大；如果你的数据量较小，本地也就可以
```
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
```
```
geneTree = hclust(as.dist(dissTOM), method = "average");
```
### 绘制聚类树(树状图)
```
pdf(file="4_Gene clustering on TOM-based dissimilarity.pdf",width=24,height=18)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
```
![](https://files.mdnice.com/user/40177/1aa89c91-2ced-46f6-866f-82d307d894f2.png)
### 加入模块
```
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="5_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
```
![](https://files.mdnice.com/user/40177/7fa4b9f3-f9ae-4bc3-b1cf-c8dc17e3ba6c.png)
## 合并模块
做出的WGCNA分析中，具有较多的模块，但是在我们后续的分析中，是使用不到这么多的模块，以及模块越多对我们的分析越困难，那么就必须合并模块信息。具体操作如下。
```
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

######剪切高度可修改
MEDissThres = 0.4  
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
```
![](https://files.mdnice.com/user/40177/bddbfcb4-0eaa-4742-b693-c001945c8def.png)
**合并及绘图**
```
 = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
table(mergedColors)

#sizeGrWindow(12, 9)
pdf(file="7_merged dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
```
![](https://files.mdnice.com/user/40177/0539c905-d51a-427a-b27d-8b7d5088ffc6.png)
## Rename to moduleColors
```
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
```
**详细代码教程可看：**[WGCNA分析 | 全流程分析代码 | 代码一](https://mp.weixin.qq.com/s/M0LAlE-61f2ZfpMiWN-iQg)
![image](https://user-images.githubusercontent.com/62915930/220238568-cc598afe-5083-402a-9c29-b3f03aaa440a.png)

> 小杜的生信筆記 ，主要发表或收录生物信息学的教程，以及基于R的分析和可视化（包括数据分析，图形绘制等）；分享感兴趣的文献和学习资料!!
