library(survival)
library(survminer)

head(phe$ID)
head(colnames(exprset))
all(phe$ID==substr(colnames(exprset),1,12))
BiocManager::install("randomForest")
library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
head(colnames(exprset))

x=t(log2(exprset+1))
y=phe$event
View(x[1:10,1:6])
table(y)
tmp=as.vector(table(y))
View(tmp)

num_class=length(tmp)
num_class

min_size=tmp[order(tmp,decreasing = F)[1]]
min_size
order(tmp,decreasing = F)
samplesize=rep(min_size,num_class)
samplesize

rf_output=randomForest(x=x,y=y,importance = T,ntree = 10001,proximity = T)
head(rf_output)
View(rf_output)
str(rf_output)

rf_importance=importance(rf_output,scale=F)
str(rf_importance)
head(rf_importance)
varImpPlot(rf_output,type=2,n.var=30,scale=F,main="top 30 predictors",cex=0.7)

target_labels=as.vector(y)
table(target_labels)
MDSplot(rf_output,factor(y),k=2,xlab="",ylab="",palette = c("red","blue"),main="MDS plot")
View(target_labels)
?MDSplot
dim(y)
head(y)
class(y)

head(rf_importance)
?order

choose_gene=rownames(tail(rf_importance[order(rf_importance[,2]),],50))
head(choose_gene)
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]

n=t(scale(t(log2(choose_matrix+1))))
n[n>2]=2
n[n<(-2)]=-2
n[1:4,1:4]

group_list=ifelse(as.numeric(substr(colnames(expr),14,15))<10,"tumour","normal")

annotation_col=data.frame(group_list=group_list)
rownames(annotation_col)=colnames(expr)

pheatmap(n,show_colnames = F,annotation_col = annotation_col)

library(ggfortify)
df=as.data.frame(t(choose_matrix))
View(df[1:4,1:4])
df$group=group_list
autoplot(prcomp(df[,1:(ncol(df)-1)]),data=df,colour="group")+theme_bw()


library(FactoMineR)
library(factoextra)
dat.pca=PCA(t(choose_matrix),graph = F)
fviz_pca_ind(dat.pca,rep=T,geom.ind = "point",col.ind = group_list,addEllipses = T,legend.title="Groups")
save.image(file="KIRC-random-forest.Rdata")
