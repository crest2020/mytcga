if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install(c("Seurat"))
also installing the dependencies â€˜GenomeInfoDbDataâ€?, â€˜lambda.râ€?, â€˜futile.optionsâ€?, â€˜digestâ€?, â€˜gtableâ€?, â€˜lazyevalâ€?, â€˜reshape2â€?, â€˜rlangâ€?, â€˜scalesâ€?, â€˜tibbleâ€?, â€˜viridisLiteâ€?, â€˜withrâ€?, â€˜GenomicRangesâ€?, â€˜Biobaseâ€?, â€˜IRangesâ€?, â€˜GenomeInfoDbâ€?, â€˜beeswarmâ€?, â€˜viporâ€?, â€˜matrixStatsâ€?, â€˜HDF5Arrayâ€?, â€˜gridExtraâ€?, â€˜RcppAnnoyâ€?, â€˜RcppHNSWâ€?, â€˜irlbaâ€?, â€˜rsvdâ€?, â€˜futile.loggerâ€?, â€˜snowâ€?, â€˜BHâ€?, â€˜SingleCellExperimentâ€?, â€˜ggplot2â€?, â€˜BiocGenericsâ€?, â€˜SummarizedExperimentâ€?, â€˜ggbeeswarmâ€?, â€˜DelayedArrayâ€?, â€˜DelayedMatrixStatsâ€?, â€˜S4Vectorsâ€?, â€˜viridisâ€?, â€˜Rcppâ€?, â€˜BiocNeighborsâ€?, â€˜BiocSingularâ€?, â€˜BiocParallelâ€?, â€˜beachmatâ€?

BiocManager::install(c("multtest"))
install.packages()
library(scater)
library(ggplot2)
library(Seurat)
install.packages("Seurat")
Package which is only available in source form, and may need compilation of C/C++/Fortran: â€˜Seuratâ€?

BiocManager::install("sva")
?BiocManager::install

library(Seurat)
?seurat
example(Read10X())
?CreateAssayObject()
getwd()
install.packages("OSCAUtils")
BiocManager::install("ggpubr")
library(OSCAUtils)

library(scRNAseq)
sce.416b=LunSpikeInData(which="416b")
?LunSpikeInData
sce <- LunSpikeInData()
library(AnnotationHub)
an=AnnotationHub()
cd=c(1,2,3,45)
grepl(3,cd)
BiocManager::install("pheatmap")
library(gridExtra)
?plotExpression
library(scater)
?grid.arrange

options(download.file.method="libcurl")
options(url.method="libcurl")
library(GEOquery)
BiocManager::install("GEOquery_2.54.1.zip",site_repository = NULL)
getwd()
.libPaths()
?install
install.packages("GEOquery_2.54.1.zip",repos = NULL)
library(GEOquery)
eset1=getGEO("GSE83521",destdir = '.',getGPL = F)
eset2=getGEO("GSE89143",destdir = '.',getGPL=F)
options(stringsAsFactors = F)
eset1
eset2
experimentData(eset1[[1]])
exp1=exprs(eset1[[1]])
exp1[1:4,1:4]
boxplot(exp1,las=2)
dim(exp1)
pd1=pData(eset1[[1]])

gpl1=eset1[[1]]@annotation
head(gpl1)

exp2=exprs(eset2[[1]])
exp2[1:4,1:4]
dim(exp2)
boxplot(exp2,las=2)
library(limma)
exp2=log2(exp2+1)
boxplot(exp2,las=2)
exp2=normalizeBetweenArrays(exp2)
boxplot(exp2,las=2)
pd2=pData(eset2[[1]])
gpl2=eset2[[1]]@annotation
gpl2
gpl1
index=sort.int(pd2$characteristics_ch1,index.return = T)
index
index1=sort(pd2$characteristics_ch1,index.return = T)
index1
?sort.int
pd2=pd2[index$ix,]
dim(exp2)
rownames(pd2)
match(rownames(pd2),colnames(exp2))
exp2=exp2[,match(rownames(pd2),colnames(exp2))]
head(exp2)
##probe annotation
library(GEOquery)
gpl=getGEO('GPL19978',destdir = '.')
dim(gpl)
datatable=gpl@dataTable
datatable@columns
table=(datatable@table)
colnames(Table(gpl))
gpl
?Table
ids=Table(gpl)
ids=ids[-c(1:15),c(1,2)]
x1=exp1[rownames(exp1)%in%ids$ID,]
x2=exp2[rownames(exp2)%in%ids$ID,]
boxplot(x1,las=2)
boxplot(x2,las=2)
cg=intersect(rownames(x1),rownames(x2))
head(cg)
x_merge=cbind(x1[cg,],x2[cg,])
head(x_merge)
boxplot(x_merge,las=2)

##remove batch effect
pd1$title
pd2$characteristics_ch1
group_list1=c(rep('tumour',6),rep('normal',6),rep(c('tumour','normal'),each=3))
group_list
gse=c(rep('GSE83527',12),rep('GSE89143',6))
gse
table(group_list,gse)
dat=x_merge

library(sva)
library(limma)
dat[1:4,1:4]
batch=gse
batch
model.matrix(~group_list1)
design
group_list1
ex_b_limma=removeBatchEffect(dat,batch = batch,design = design)
dim(ex_b_limma)
boxplot(ex_b_limma,las=2)
##differential analysis
library(limma)
fit=lmFit(ex_b_limma,design)
fit
fit=eBayes(fit)
fit
options(digits = 4)
topTable(fit,coef = 2,adjust.method = "BH")
design
deg=topTable(fit,coef=2,adjust="BH",number = Inf)
head(deg)
##volcano plot
nrdeg=deg
head(nrdeg)
attach(nrdeg)
plot(logFC,-log10(P.Value))

library(ggpubr)
nrdeg$v=-log10(P.Value)
ggscatter(nrdeg,x="logFC",y="v",size=0.5)
nrdeg$g=ifelse(nrdeg$P.Value>0.05,'stable',ifelse(nrdeg$logFC>1,"up",ifelse(nrdeg$logFC<(-1),"down","stable")))
table(nrdeg$g)
nrdeg$name=rownames(nrdeg)
head(nrdeg)
ggscatter(nrdeg,x="logFC",y="v",size=0.5,color = 'g')
ggscatter(nrdeg,x="logFC",y="v",size=0.5,color = 'g',label = 'name',repel=T,label.select = head(rownames(nrdeg)),palette = c("#00AFBB", "#E7B800", "#FC4E07"))

##pheatmap
up=nrdeg[nrdeg$g=='up',]
down=nrdeg[nrdeg$g=='down',]
x=rbind(up,down)
x$ID=rownames(x)
x
dim(x)
y=merge(x,ids,by="ID")
y
ex_b_limma2=ex_b_limma[match(y$ID,rownames(ex_b_limma)),]
head(y$ID)
head(rownames(ex_b_limma))
rownames(ex_b_limma2)=y$circRNA
library(pheatmap)
pheatmap(ex_b_limma2,show_colnames = T,show_rownames = T)
n=t(scale(t(ex_b_limma2)))
head(n)
n[n>2]=2
n[n<(-2)]=-2
pheatmap(n,show_colnames = F,show_rownames = F)
ac=data.frame(sampleType=group_list,GeoDatabase=gse)
ac
rownames(ac)=colnames(n)
pheatmap(n,show_colnames = T,show_rownames = F,cluster_cols = T,annotation_col = ac,fontsize = 8)
