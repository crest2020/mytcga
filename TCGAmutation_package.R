library(TCGAmutations)
tcga_available()
tcga_load(study="LGG")
tcaga_mc3=tcga_lgg_mc3
laml=tcaga_mc3
phe=as.data.frame(laml@clinical.data)

str(laml)
View(laml@data[1:10,])

table(phe$vital_status)
phe=phe[phe$vital_status%in%c("Alive","Dead"),]
colnames(phe)
phe$time=0
pos=which(phe$vital_status=="Dead")
phe$time[pos]=phe$days_to_death[pos]
pos=which(phe$vital_status=="Alive")
phe$time[pos]=phe$days_to_last_followup[pos]
phe=phe[,c("Tumor_Sample_Barcode","gender","vital_status",'time')]
View(phe)
phe$vital_status=ifelse(phe$vital_status=="Alive",0,1)
table(nchar(as.character(phe$Tumor_Sample_Barcode))==12)
length(phe$Tumor_Sample_Barcode)
length("ahgie")
table(duplicated(phe$Tumor_Sample_Barcode))


library(survival)
library(survminer)
mysurv=with(phe,Surv(time,vital_status))
sfit=survfit(mysurv~gender,data=phe)
p=ggsurvplot(sfit,data=phe,pval=T,xlab="Time in days",ggtheme = theme_light())

print(p)
save.image(file="TCGAmutations-survival.Rdata")
