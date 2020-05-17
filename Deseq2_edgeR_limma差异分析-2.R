options(stringsAsFactors = F)
load(file = "TCGA-KIRC-miRNA-example.Rdata")
sessionInfo()
dim(expr)
dim(meta)
##DESeq2
group_list=ifelse(as.numeric(substr(colnames(expr),14,15))<10,"tumour","normal")
head(group_list)
table(group_list)
exprset=na.omit(expr)
dim(exprset)
View(exprset[1:10,1:10])
library(DESeq2)
coldata=data.frame(row.names = colnames(exprset),group_list=group_list)
View(coldata)
dds=DESeqDataSetFromMatrix(countData = exprset,colData = coldata,design = ~group_list)
dds=DESeq(dds)
dim(exprset)
browseVignettes("DESeq2")
boxplot(exprset,las=2)
res=results(dds)
res
resordered=res[order(res$padj),]
deg=as.data.frame(resordered)
head(deg)
deseq_deg=na.omit(deg)
dim(resordered)
dim(deg)
dim(deseq_deg)
nrdeg=deseq_deg[,c(2,6)]
colnames(nrdeg)=c("log2foldchange","pvalue")
#boxplot(exprset,las=2)
#boxplot(assays(dds),las=2)
#boxplot(log10(assays(dds)),las=2)
#edgeR
library(edgeR)

d=DGEList(counts = exprset,group = factor(group_list))
keep=rowSums(cpm(d)>1)>=2
d=d[keep, ,keep.lib.sizes=F]
d$samples$lib.size=colSums(d$counts)
d=calcNormFactors(d)
d
d$samples
dge=d
View(colnames(dge))
design=model.matrix(~0+factor(group_list))
rownames(design)=colnames(dge)

colnames(design)=levels(factor(group_list))
head(design)
head(colnames(exprset))

dge=estimateGLMCommonDisp(dge,design)
dge=estimateGLMTrendedDisp(dge, design)
dge=estimateGLMTagwiseDisp(dge, design)

fit=glmFit(dge, design)
lrt=glmLRT(fit,contrast = c(-1,1))
nrdeg=topTags(lrt,n=nrow(dge))
nrdeg=as.data.frame(nrdeg)
head(nrdeg)
edge_deg=nrdeg
nrdeg=edge_deg[,c(1,5)]
colnames(nrdeg)=c("log2foldchange","padj")

#limma

library(limma)
design=model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprset)
design
dge=DGEList(counts=exprset)
dge=calcNormFactors(dge)
logcpm=cpm(dge,log=T,prior.count = 3)
v=voom(dge,design,plot=T,normalize.method = "quantile")
fit=lmFit(v,design)
cont.matrix=makeContrasts(contrasts = c("tumour-normal"),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
tempoutput=topTable(fit2,coef="tumour-normal",n=Inf)
deg_limma_voom=na.omit(tempoutput)
head(deg_limma_voom)
nrdeg=deg_limma_voom[,c(1,4)]
colnames(nrdeg)=c("log2foldchange","pvalue")
save.image(file="deseq2-edgeR-limma.Rdata")
