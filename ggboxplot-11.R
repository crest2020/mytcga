library(maftools)
table(group_list)
library(survival)
library(survminer)
(exprset['hsa-mir-10b',])
dim(expr)
expr[1:4,1:4]
head(phe$ID)
head(colnames(exprset))
head(mut)
dat=data.frame(gene=log2(exprset['hsa-mir-10b',]+1),stage=phe$stage)

head(dat)
boxplot(dat$gene~dat$stage)

library(ggpubr)
p=ggboxplot(dat,x="stage",y="gene",color="stage",palette = "jco",add="jitter")
p+stat_compare_means()
BiocManager::install("ggstatsplot")
library(ggstatsplot)
ggbetweenstats(data=dat,x=stage,y=gene)
?ggbetweenstats
dev.off()

packageVersion("ggplot2")

library(ggplot2)
install.packages("ggplot2", dependencies=TRUE)

res.aov=aov(gene~stage,data=dat)
summary(res.aov)
TukeyHSD(res.aov)

VHL_mut=substr(as.character(as.data.frame(mut[mut$Hugo_Symbol=="VHL","Tumor_Sample_Barcode"])[,1]),1,12)
as.data.frame(mut[mut$Hugo_Symbol=="VHL","Tumor_Sample_Barcode"])


dat1=data.frame(gene=log2(exprset["hsa-mir-10b",]),mut=substr(colnames(exprset),1,12)%in%VHL_mut)
View(dat1)

library(ggpubr)
p=ggboxplot(dat1,x="mut",y="gene",color="mut",palette = "jco",add="jitter")
p+stat_compare_means(method="t.test")

ggplot(dat1,aes(x=mut,y=gene))+geom_boxplot()+geom_jitter()+geom_violin()+theme_bw()


head(phe)
head(mut)

dat2=data.frame(gene1=log2(exprset['hsa-mir-10b',]+1),gene2=log2(exprset['hsa-mir-143',]+1),stage=phe$stage)

dat2$stage=as.factor(dat2$stage)
sp=ggscatter(dat2,x="gene1",y="gene2",add="reg.line",add.params = list(color="blue",fill="lightgray"),conf.int = T)
sp
sp+stat_cor(method = "pearson",label.x = 15,label.y=20)

sp=ggscatter(dat2,x="gene1",y="gene2",color = "stage",palette = "jco",add="reg.line",conf.int = T)
sp+stat_cor(aes(color=stage),label.x = 15)

##split cohort
dim(expr)
set.seed(1234567890)
k=sample(1:593,300)
t_exp=expr[,k]
v_exp=expr[,-k]
table(ifelse(substr(colnames(t_exp),14,15)=='01','tumor','normal'))
table(ifelse(substr(colnames(v_exp),14,15)=='01','tumor','normal'))
head(colnames(exprset))
table(substr(colnames(t_exp),14,15)=='02')
lapply(1:10,function(x){table(as.numeric(substr(colnames(expr),14,15))==x)})


t_tumour=t_exp[,substr(colnames(t_exp),14,15)=='01']
v_tumour=v_exp[,substr(colnames(v_exp),14,15)=='01']
t_phe=phe[match(substr(colnames(t_tumour),1,12),phe$ID),]
v_phe=phe[match(substr(colnames(v_tumour),1,12),phe$ID),]
table(t_phe$stage)
table(v_phe$stage)
install.packages("caret",dependencies = T)
library(caret)
set.seed(123456789)
sam=createDataPartition(phe$event,p=0.5,list=F)
sam
train=phe[sam,]
test=phe[-sam,]

prop.table(table(train$stage))
prop.table(table(test$stage))
t_tumour[1:4,1:4]

install.packages(c("lars","glmnet"),dependencies = T)
library(lars)
library(glmnet)

x=t(log2(t_tumour+1))
y=t_phe$event
cv_fit=cv.glmnet(x=x,y=y,alpha=1,nlambda=1000)
plot.cv.glmnet(cv_fit)

cv_fit$lambda.min
cv_fit$lambda.1se

model_lasso=glmnet(x=x,y=y,alpha=1,lambda = cv_fit$lambda.1se)
lasso.prob=predict(cv_fit,newx = t(log2(v_tumour+1)),s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(v_phe$event,lasso.prob)
head(re)
View(lasso.prob)
library(ROCR)
?prediction
pred=prediction(re[,2],re[,1])
pred
View(re[,2])
View(re[,1])
View(re)
perf=performance(pred,"tpr","fpr")
performance(pred,"auc")
plot(perf,colorize=F,col="black")
lines(c(0,1),c(0,1),col="gray",lty=4)
library(survival)
new_dat=phe
library(timeROC)
new_dat$fp=as.numeric(lasso.prob[,1])
head(new_dat$fp)
?timeROC
with(new_dat,ROC<<-timeROC(T=time,delta = event,marker = fp,cause = 1,weighting = "marginal",times =c(60,100),ROC=T,iid=T))
plot(ROC,time=60,col="blue",add=F)

plot(ROC,time=100,col="red",add=F)
ROC$AUC

library(ROCR)
library(glmnet)
library(caret)

save.image(file="timeROC.Rdata")
