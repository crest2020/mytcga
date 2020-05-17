options(stringsAsFactors = F)
load('TCGA-KIRC-miRNA-example.Rdata')
dim(expr)
dim(meta)
group_list=ifelse(as.numeric(substr(colnames(expr),14,15))<10,"tumour","normal")
table(group_list)
exprset=na.omit(expr)
library(survival)
library(survminer)

exprset=exprset[,group_list=="tumour"]
head(meta)
colnames(meta)
meta[,3][is.na(meta[,3])]=0
meta[,4][is.na(meta[,4])]=0
meta$days=as.numeric(meta[,3])+as.numeric(meta[,4])
View(meta[,c(2,3,4)])
colnames(meta)
meta=meta[,c(1:2,5:9)]
colnames(meta)=c("ID","event","race","age","gender","stage","days")
meta
colnames(meta)

meta$event=ifelse(meta$event=="alive",0,1)
meta$age=as.numeric(meta$age)
library(stringr)
View(meta[1:10,])
meta$stage=str_split(meta$stage,' ',simplify = T)[,2]
table(meta$stage)
meta$age_group=ifelse(meta$age>median(meta$age),"older","younger")
colnames(meta)
meta$time=meta$days/30

# 用my.surv <- surv(OS_MONTHS,OS_STATUS=='DECEASED')构建生存曲线。
# 用kmfit2 <- survfit(my.surv~TUMOR_STAGE_2009)来做某一个因子的KM生存曲线。
# 用 survdiff(my.surv~type, data=dat)来看看这个因子的不同水平是否有显著差异，其中默认用是的logrank test 方法。
# 用coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung) 来检测自己感兴趣的因子是否受其它因子(age,gender等等)的影响。

meta$ID=toupper(meta$ID)
phe=meta

dim(phe)
phe=phe[match(substr(colnames(exprset),1,12),phe$ID),]
head(phe)
dim(phe)
head(meta$ID)
head(colnames(exprset))
exprset[1:4,1:4]
sfit=survfit(Surv(time,event)~gender,data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval=T)
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
sfit$time

sfit1=survfit(Surv(time, event)~gender, data=phe)
sfit2=survfit(Surv(time, event)~age_group, data=phe)

arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1,pval =TRUE, data = phe, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2,pval =TRUE, data = phe, risk.table = TRUE)
dev.off()

## 挑选感兴趣的基因做生存分析
tmp=as.data.frame(row.names(exprset))
g1='hsa-mir-21'
g2='hsa-mir-143'
g3='hsa-mir-192'
g4='hsa-mir-183'
g5='hsa-mir-10b'
gs=c(g1,g2,g3,g4,g5)
gs
splots=lapply(gs,function(g){
  phe$gene=ifelse(exprset[g,]>median(exprset[g,]),"high","low")
  table(phe$gene)
  sfit1=survfit(Surv(time,event)~gene,data=phe)
  ggsurvplot(sfit1,pval=T,risk.table = T,data=phe)
})
View(phe[1:4,])
arrange_ggsurvplots(splots,print=T,ncol=2,nrow=3,risk.table.height = 0.4)
dev.off()

##批量生存分析
View(exprset[1:6,1:6])
mysurv=with(phe,Surv(time,event))
log_rank_p=apply(exprset,1,function(gene){
  phe$group=ifelse(gene>median(gene),"high","low")
  data.survdiff=survdiff(mysurv~group,data=phe)
  p.val=1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  return(p.val)
})
head(log_rank_p)
View(log_rank_p)
log_rank_p=sort(log_rank_p)
boxplot(log_rank_p)
table(log_rank_p<0.01)
log_rank_p[log_rank_p<0.01]
gs %in% names(log_rank_p[log_rank_p<0.01])

library(pheatmap)
choose_gene=names(log_rank_p[log_rank_p<0.01])
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]

choose=log2(choose_matrix+1)
?scale
#scale 函数默认对列进行标准化，所以两次转置
View(choose[1:6,])
n=t(scale(t(log2(choose_matrix+1))))
dim(n)
View(n[1:5,])
n[n>2]=2
n[n<(-2)]=-2
n[1:4,1:4]
head(group_list)
annotation_col=data.frame(group_list=group_list)
rownames(annotation_col)=colnames(expr)
pheatmap(n,show_colnames = F,annotation_col = annotation_col,filename = "log_rank_genes_selected.png")
BiocManager::install("ggfortify")
library(ggfortify)
df=as.data.frame(t(choose_matrix))
View(df[1:10,1:10])
df$group=group_list
autoplot(prcomp(df[,1:(ncol(df)-1)]),data=df,colour='group')+theme_bw()
dev.off()


library(FactoMineR)
library(factoextra)
dat.pca=PCA(t(choose_matrix),graph=F)
fviz_pca_ind(dat.pca,repel = T,geom.ind = "point",col.ind = group_list,addEllipses = T,legend.title="Groups")

dim(exprset)
dim(phe)
dim(meta)
all(colnames(exprset)==rownames(phe))
?sort

order

rownames(phe)=phe$ID
phe_sort=phe[order(rownames(phe)),]
exprset_sort=exprset[,order(colnames(exprset))]
all(colnames(exprset_sort)==rownames(phe_sort))
dim(exprset_sort)
dim(phe_sort)
head(colnames(exprset_sort))
head(rownames(phe_sort))
View(phe[1:5])
table(duplicated(phe$ID))
duplicated(phe$ID)

View(exprset[1:9,1:9])
View(phe)
View(meta[1:10,])
###cox analysis

cox_result=apply(exprset,1,function(gene){
  #gene=exprset[1,]
  group=ifelse(gene>median(gene),"high","low")
  survival_dat=data.frame(group=group,stage=phe$stage,age=phe$age,gender=phe$gender,stringsAsFactors = F)
  m=coxph(mysurv~gender+age+stage+group,data = survival_dat)
  beta=coef(m)
  se=sqrt(diag(vcov(m)))
  HR=exp(beta)
  HRse=HR*se
  
  
  tmp=round(cbind(coef=beta,se=se,z=beta/se,p=1-pchisq((beta/se)^2,1),HR=HR,HRse=HRse,HRz=(HR-1)/HRse,HRp=1-pchisq(((HR-1)/HRse)^2,1),HRCILL=exp(beta-qnorm(.975,0,1)*se),HRCIUL=exp(beta+qnorm(.975,0,1)*se)),3)
  return(tmp['grouplow',])
})
###dim(tmp)
head(tmp)
dim(cox_result)
View(cox_result)
m
coef(m)
se
HR
class(HRse)
dim(tmp)
View(tmp)
View(HRse)###


cox_result=t(cox_result)
table(cox_result[,4]<0.05)
View(cox_result)

##logrank test
mysurv=with(phe,Surv(time,event))
log_rank_p=apply(exprset,1,function(gene){
  #gene=exprset[1,]
  phe$group=ifelse(gene>median(gene),"high","low")
  data.survdiff=survdiff(mysurv~group,data=phe)
  p.val=1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  return(p.val)
})
install.packages("VennDiagram")
library("VennDiagram")
venn.list=list(cox=rownames(cox_result[cox_result[,4]<0.05,]),log=names(log_rank_p[log_rank_p<0.05]))
venn.plot=venn.diagram(venn.list,NULL,fill=c("darkmagenta", "darkblue"),alpha=c(0.5,0.5),cex=2,cat.fontface=4,main="overlap of coxph and log-rank test")
grid.draw(venn.plot)


head(colnames(exprset))
head(phe$ID)
all(substring(colnames(exprset),1,12)==phe$ID)

a=data.frame(name=c("a","b","c"),value=c(4,7,5))
a
b=data.frame(n=c("c","b","a"),v=c(100,200,300))
b
b[match(a$name,b$n),]
match(a$name,b$n)

##coxph-ggforst
library(survival)
library(survminer)
sfit=survfit(Surv(time,event)~gender,data=phe)
View(summary(sfit))
ggsurvplot(sfit,conf.int = F,pval = T)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),risk.table = F,pval=T,conf.int = T,xlab="Time in months",ggtheme = theme_light(),ncensor.plot=T)

sfit1=survfit(Surv(time,event)~gender,data=phe)
sfit2=survfit(Surv(time,event)~age_group,data = phe)
splots=list()
splots[[1]]=ggsurvplot(sfit1,pval=T,data = phe,risk.table = T)
splots[[2]]=ggsurvplot(sfit2,pval=T,data=phe,risk.table = T)
arrange_ggsurvplots(splots,print=T,ncol=2,nrow=1,risk.table.height = 0.4)
dev.off()
