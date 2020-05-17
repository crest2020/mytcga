options(stringsAsFactors =F)

library(RTCGA.miRNASeq)
??RTCGA.miRNASeq
s=rownames(LUAD.miRNASeq)[seq(1,nrow(LUAD.miRNASeq),by=3)]
head(s)

expr=expressionsTCGA(LUAD.miRNASeq)
dim(expr)
expr[1:10,1:4]
expr=as.data.frame(expr[seq(1,nrow(expr),by=3),3:ncol(expr)])
dim(expr)
mi=colnames(expr)
expr=apply(expr,1,as.numeric)
colnames(expr)=s
rownames(expr)=mi
expr[1:4,1:4]
expr=na.omit(expr)
dim(expr)
expr=expr[apply(expr,1,function(x){sum(x>1)>10}),]
dim(expr)

library(RTCGA.clinical)
??RTCGA.clinical
meta=LUAD.clinical
#meta2=LUAD.clinical
View(meta[1:10,1:10])
tmp=as.data.frame(colnames(meta))
View(tmp)
dim(meta)

grep('patient.bcr_patient_barcode',colnames(meta))
grep('patient.days_to_last_followup',colnames(meta))
grep("patient.days_to_death",colnames(meta))
grep("patient.vital_status",colnames(meta))
grep("patient.race",colnames(meta))
grep("patient.age_at_initial_pathologic_diagnosis",colnames(meta))
grep("patient.gender",colnames(meta))
grep("patient.stage_event.pathologic_stage",colnames(meta))
View(meta[,"patient.stage_event.pathologic_stage"])
meta=as.data.frame(meta[c('patient.bcr_patient_barcode',"patient.vital_status","patient.days_to_death",'patient.days_to_last_followup',"patient.race","patient.age_at_initial_pathologic_diagnosis","patient.gender","patient.stage_event.pathologic_stage")])

View(meta[1:10,])
group_list=ifelse(as.numeric(substr(colnames(expr),14,15))<10,"tumour","normal")
table(group_list)

exprset=na.omit(expr)
dim(exprset)
dim(expr)

exprset=exprset[,group_list=="tumour"]
exprset[1:4,1:4]


meta[,3][is.na(meta[,3])]=0
meta[,4][is.na(meta[,4])]=0
meta$days=as.numeric(meta[,3])+as.numeric(meta[,4])
meta=meta[,c(1:2,5:9)]
colnames(meta)
colnames(meta)=c('ID','event','race','age','gender','stage',"days") 
View(meta[1:10,])

meta$event=ifelse(meta$event=="alive",0,1)
meta$age=as.numeric(meta$age)
library(stringr)
#meta$stage=
a=str_split(meta$stage[1]," ",simplify = T)
meta$stage[1]
meta$stage[2]
a
View(meta[1:10,])
dim(meta2)
dim(meta)
all(rownames(meta2)==rownames(meta))
#meta$stage=meta2[,"patient.stage_event.pathologic_stage"]
View(meta[1:10,])

meta$stage=str_split(meta$stage," ",simplify = T)[,2]
table(meta$stage)
meta$age_group=ifelse(meta$age>median(meta$age,na.rm = T),"older","younger")
table(meta$race)
meta$time=meta$days/30
phe=meta
save(expr,exprset,meta,phe,group_list,file="TCGA-LUAD-miRNA-example.Rdata")

phe$ID=toupper(phe$ID)
phe=phe[match(substr(colnames(exprset),1,12),phe$ID),]
View(phe[1:10,])

save(expr,exprset,meta,phe,group_list,file="TCGA-LUAD-miRNA-example.Rdata")

## 挑选感兴趣的基因构建coxph模型 
# miR-31, miR-196b, miR-766, miR-519a-1, miR-375, miR-187, miR-331 and miR-101-1
# hsa-mir-31,hsa-mir-196b,hsa-mir-766,hsa-mir-519a-1,hsa-mir-375,hsa-mir-187,hsa-mir-331,hsa-mir-101-1
View(exprset[1:6,1:6])
e=t(exprset[c('hsa-mir-31','hsa-mir-196b','hsa-mir-766','hsa-mir-519a-1','hsa-mir-375','hsa-mir-187','hsa-mir-331','hsa-mir-101-1'),])
View(e[1:6,1:6])
e=log2(e+1)
colnames(e)=c('miR31','miR196b','miR766','miR519a1','miR375','miR187','miR331','miR101')
dat=cbind(phe,e)

dat$gender=factor(dat$gender)
dat$stage=factor(dat$stage)
colnames(dat)

##survival

library(survival)
library(survminer)
s=Surv(time,event)~miR31+miR196b+miR766+miR519a1+miR375+miR187+miR331+miR101
s
model=coxph(s,data=dat)
summary(model,data=dat)

options(scipen = 1)
ggforest(model,data=dat,main="Hazard ratio",cpositions = c(0.1,0.22,0.4),fontsize = 1,refLabel = "1",noDigits = 4)

new_dat=dat
fp=predict(model,new_dat,type="risk")
boxplot(fp)
fp=predict(model,new_dat,type="expected")
boxplot(fp)
fp=predict(model,new_dat)
boxplot(fp)

basehaz(model)
library(Hmisc)
options(scipen = 200)
with(new_dat,rcorr.cens(fp,Surv(time,event)))

head(fp)
library(cowplot)
library(pheatmap)
fp_dat=data.frame(s=1:length(fp),v=as.numeric(sort(fp)))
dim(fp_dat)
head(fp_dat)
sur_dat=data.frame(s=1:length(fp),t=phe[names(sort(fp)),"time"],e=phe[names(sort(fp)),"event"])
head(sur_dat)
sur_dat$e=ifelse(sur_dat$e==0,"alive","death")
exp_dat=new_dat[names(sort(fp)),10:17]
head(sur_dat)
plot.point=ggplot(fp_dat,aes(x=s,y=v))+geom_point()
print(plot.point)

plot.surv=ggplot(sur_dat,aes(x=s,y=t))+geom_point()
print(plot.surv)
mycolors=colorRampPalette(c("black","green","red"),bias=1.2)(100)
head(mycolors)
tmp=t(scale(exp_dat))
tmp[tmp>1]=1
tmp[tmp<(-1)]=-1
plot.h=pheatmap(tmp,col=mycolors,show_colnames = F,cluster_cols = T)
plot.h=pheatmap(tmp,col=mycolors,show_colnames = F,cluster_cols = F)
plot_grid(plot.point,plot.surv,plot.h$gtable,labels = c("A","B","C"),align="v",ncol=1)

View(dat[1:6,])
save.image(file="TCGA-LUAD-survival_input.Rdata-step07.Rdata")
