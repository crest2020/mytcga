options(stringsAsFactors = F)
load("TCGA-KIRC-miRNA-example.Rdata")
dim(expr)
dim(meta)

library(survival)
library(survminer)
group_list=ifelse(as.numeric(substr(colnames(expr),14,15))<10,"tumour","normal")
table(group_list)
exprset=na.omit(expr)
dim(exprset)
exprset=expr[,group_list=='tumour']

meta[,3][is.na(meta[,3])]=0
meta[,4][is.na(meta[,4])]=0
meta$days=as.numeric(meta[,3])+as.numeric(meta[,4])
meta=meta[,c(1:2,5:9)]
colnames(meta)=c('ID','event','race','age','gender','stage',"days")
meta$event=ifelse(meta$event=='alive',0,1)
meta$age=as.numeric(meta$age)
View(meta[1:10,])
library(stringr)
meta$stage=str_split(meta$stage,' ',simplify = T)[,2]
table(meta$stage)
meta$age_group=ifelse(meta$age>median(meta$age),'older','younger')
table(meta$race)
meta$time=meta$days/30
phe=meta
phe$ID=toupper(phe$ID)
?toupper
View(phe[1:10,])
phe=phe[match(substr(colnames(exprset),1,12),phe$ID),]
View(phe[1:10,])
View(exprset[1:10,1:10])
head(phe$ID)
head(colnames(exprset))
dim(exprset)
dim(phe)

phe=phe[match(substr(colnames(exprset),1,12),phe$ID),]
all(phe$ID==substr(colnames(exprset),1,12))

## 挑选感兴趣的基因构建coxph模型 
e=t(exprset[c('hsa-mir-21','hsa-mir-143','hsa-mir-10b','hsa-mir-192','hsa-mir-183'),])
dim(e)
e=log2(e)
colnames(e)=c('miR21','miR143','miR10b','miR192','miR183')
dat=cbind(phe,e)
View(dat[1:10,])
dat$gender=factor(dat$gender)
dat$stage=factor(dat$stage)
colnames(dat)
s=Surv(time, event) ~ miR21+miR143+miR10b+miR192+miR183
s
model=coxph(s,data=dat)
summary(model,data=dat)
options(scipen=1)
ggforest(model,data=dat,main="Hazard ratio",cpositions = c(0.1,0.22,0.4),fontsize = 1,refLabel = "1",noDigits = 4)


fp=predict(model)
summary(model,data=dat)
library(Hmisc)
options(scipen = 200)
with(dat,rcorr.cens(fp,Surv(time,event)))
?rcorr.cens
Surv(time,event)
?Surv
with(dat,Surv(time,event))
save(exprset,phe,file="TCGA-KIRC-miRNA-survival_input.Rdata")
save.image(file="coxph-ggforest-6.Rdata")
