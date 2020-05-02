options(stringsAsFactors = F)
library(AnnoProbe)
library(Biobase)
gset=geoChina("GSE20685")
gset
eset=exprs(gset[[1]])
?exprs
head(eset)
pheno=pData(gset[[1]])
head(pheno)
dim(pheno)
dim(eset)
all(colnames(eset)==rownames(pheno))
checkGPL(gset[[1]]@annotation)
gset[[1]]@annotation
e2g=idmap(gset[[1]]@annotation)
?filterEM
eset=filterEM(eset,e2g)
head(eset,n=2)
gset
View(pheno)
dat <- pheno[,c("characteristics_ch1.3","characteristics_ch1.4")]
library(stringr)
colnames(dat)=c("status","OS")
dat$status=as.numeric(str_split(dat$status,":")[[1]][2])
?sub
View(dat)
dat$OS=as.numeric(lapply(str_split(dat$OS,":")[[1]][2]))
table(dat$status)
table(dat$OS)
str_split(dat$OS,":")[[1]][2]

unlist(str_split(dat$status[1:10],":"))[2]


View(dat)
length(dat$status[1])
nchar(dat$status[1])
dat$status=as.numeric(substr(dat$status,14,14))
nchar(dat$OS)
dat$OS=as.numeric(substr(dat$OS,29,31))
View(dat)
table(dat$status)
boxplot(dat$OS~dat$status)

library(survival)
cg=array(dim=nrow(eset))
head(cg)
cg
survdata=dat
my.surv=Surv(survdata$OS,survdata$status==0)
eset[1:4,1:4]
dim(eset)
head(dat)
#for(i in 1:nrow(eset)){
  gene_exprs=eset[1,]
  dim(gene_exprs)
  rownames(gene_exprs)
  gene_exprs
  gene_exprs=as.data.frame(t(gene_exprs))
  gene_exprs[,1]=as.numeric(gene_exprs[,1])
  gene_exprs$g=ifelse(gene_exprs[,1]>=median(gene_exprs[,1]),"high","low")
  gene_exprs$sample_name=rownames(gene_exprs)
  gene_surv=survdata
  View(gene_surv)
  gene_surv$sample_name=rownames(gene_surv)
  all(gene_surv$sample_name==gene_exprs$sample_name)
  gene_surv=cbind(gene_surv,gene_exprs)
  View(gene_surv)
  ?merge
  
  
  x=data.frame(name=c(1,2,6,3),value=c(6,7,8,23))
  y=data.frame(name=c(3,2,1,6),value2=c(8,6,3,0))
  x
  y

  
  merge(x,y,"name",sort=F)
  
  my.surv
  kmdif=survdiff(my.surv~gene_surv$g,data=gene_surv)
  kmdif$chisq
  length(kmdif$n)
  p.value=1-pchisq(kmdif$chisq,length(kmdif$n)-1)
p.value  
kmdif$call
  dim(cg)
  cg[1]=(p.value<0.05)
cg  
View(gene_surv)


browseVignettes("survival")
