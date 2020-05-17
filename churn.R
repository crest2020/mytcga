library(ggplot2)
library(ALL)
data("ALL")
library(limma)
eset=ALL[,ALL$mol.biol%in%c("BCR/ABL","ALL1/AF4")]
head(eset)
f=factor(as.character(eset$mol.biol))
f
design=model.matrix(~f)
eset$mol.biol
design
fit=eBayes(lmFit(eset,design))
head(fit)
summary(fit)
head(fit$p.value)
selected=p.adjust(fit$p.value[,2])<0.00
head(selected)
esetSel=eset[selected,]
dim(esetSel)
library(hgu95av2.db)
library(BiocInstaller)
biocLite("hgu95av2.db")
data=exprs(esetSel)
probes=rownames(data)
symbol=mget(probes,hgu95av2SYMBOL,ifnotfound = NA)
head(symbol)
symbol=do.call(rbind,symbol)
head(symbol)
symbol[is.na(symbol[,1]),1]=rownames(symbol)[is.na(symbol[,1])]
rownames(data)=symbol[probes,1]
head(data)
library(plotrix)
head(symbol)
heatmap(data,cexRow = 0.5)
design
attach(mtcars)
layout(matrix(c(1,1,2,3),2,2,byrow=T))
hist(wt)
hist(mog)
hist(mpg)
hist(disp)
layout(matrix(c(1,1,1,1,2,2,2,3,2,2,2,3),3,4,byrow=T))
matrix(c(1,1,1,1,2,2,2,3,2,2,2,3),3,4,byrow=T)
hist(wt)
hist(mpg)
hist(disp)
layout(matrix(c(1,1,2,1,1,1),nrow = 2,byrow=T))
hist(wt)
hist(mpg)
par(fig=c(0,0.85,0,0.85), new=F)
plot(mtcars$wt, mtcars$mpg, xlab="Miles Per Gallon",
        ylab="Car Weight")
par(fig=c(0,0.8,0.55,0.8), new=TRUE)
boxplot(mtcars$wt, horizontal=TRUE, axes=FALSE)
> par(fig=c(0.65,1,0,0.8),new=TRUE)
> boxplot(mtcars$mpg, axes=FALSE)
> mtext("Enhanced Scatterplot", side=3, outer=TRUE, line=-3)
> par(fig=c(0.4,0.75,0.4,0.7),new=TRUE,mar=c(2,2,0,0),mgp=c(1,.4,0),cex=1,cex.lab=0.7,cex.axis=0.7)
> hist(mtcars$mpg, main="")
par(fig=c(0,0.85,0,0.85),new=F)
plot(mtcars$wt,mtcars$mpg,xlab="miles per gallon",ylab="car weight")
par(fig=c(0,0.8,0.55,0.90),new=T)
boxplot(mtcars$wt,horizontal=T,axes=F)
set.seed(1)
x=c(rnorm(100,0,1))
x
x[23]=100
boxplot(x)
cat("The average is",mean(x),"and the sd is",sd(x))
median(x)
mad(x)
y=c(rnorm(100,0,1))
y[23]=84
library(rafalib)
mypar()
plot(x,y,main=paste0("corellation=",round(cor(x,y),3)),pch=21,bg=1,xlim=c(-3,100),ylim=c(-3,100))
abline(0,1)
mypar(1,2)
plot(x,y,main=paste0("correlation=",round(cor(x,y),3)),pch=21,bg=1,xlim=c(-3,100),ylim=c(-3,100))
plot(rank(x),rank(y),main=paste0("correlation=",round(cor(x,y,method="spearman"),3)),pch=21,bg=1,xlim=c(-3,100),ylim=c(-3,100))
abline(0,1)
rank(x)
x=2^(rnorm(100))
y=2^(rnorm(100))
ratios=x/y
ratios
mypar(1,2)
hist(ratios)
logratios=log2(ratios)
hist(logratios)
hist(rnorm(100)-rnorm(100))
wilcox.test(x,y)$p.value
t.test(x,y)$p.value
data("ChickWeight")
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
chick=reshape(ChickWeight,idvar=c("Chick","Diet"),timevar = "Time",direction = "wide")
head(chick)
?reshape
library(reshape)
md=melt(mydata,id=c("id","time"))
ID=c(1,1,2,2)
Time=c(1,2,1,2)
X1=c(5,3,6,2)
X2=c(6,5,1,4)
mydata=data.frame(ID,Time,X1,X2)
mydata
md=melt(mydata,id=c("ID","Time"))
md
cast(md,ID+variable~Time)
cast(md,ID~variable+Time)
md
cast(md,ID~variable,mean)
cast(md,Time~variable,mean)
cast(md,ID~Time,mean)
data(mtcars)
head(mtcars)
mtcars$car=rownames(mtcars)
head(mtcars)
mtcarsMelt=melt(mtcars)
head(mtcarsMelt)
?melt
cast(mtcarsMelt,car~variable,mean)
melt(mtcars,id.vars=c('cyl','gear'),variable.name='carVariable',variable.value='carValue')
data("iris")
table.iris=table(iris$Species)
table.iris
pie(table.iris)
hist(iris$Sepal.Length)
boxplot(Petal.Width~Species,data=iris)
plot(x=iris$Petal.Length,y=iris$Sepal.Width,col=iris$Species)
boxplot(iris$Sepal.Length~iris$Species)
plot(x=iris$Petal.Width,y=iris$Petal.Length,col=iris$Species)
pairs(iris[1:4],main="Edgar Anderson's Iris Data",pch=21,bg=c("red","green3","blue")[unclass(iris$Species)])
unclass(iris$Species)
head(iris[1:4])
sample(1:10)
sample(1:10,size=5)
sample.int(20,12)
dnorm(0)
dnorm(0,mean=3,sd=5)
curve(dnorm,-3,3)
pnorm(0.5)
pnorm(1.5)
pnorm(1.5,lower.tail = F)
curve(pnorm(x),-3,3)
qnorm(0.5)
qnorm(pnorm(0))
set.seed(50)
x=rnorm(100,mean=3,sd=5)
x
hist(x)
set.seed(50)
y=runif(100,0,5)
hist(y)
shapiro.test(x)
shapiro.test(y)
help("TDist")
help(Binomial)
help("distributions")
data(mtcars)
range(mtcars$mpg)
length(mtcars$mpg)
mean(mtcars$mpg)
median(mtcars$mpg)
sd(mtcars$mpg)
var(mtcars$mpg)
IQR(mtcars$mpg)
quantile(mtcars$mpg)
max(mtcars$mpg)
min(mtcars$mpg)
cummax(mtcars$mpg)
mtcars$mpg
cummin(mtcars$mpg)
summary(mtcars$mpg)
summary(mtcars)
table(mtcars$cyl)
stem(mtcars$mpg)
?stem
mtcars$mpg
length(mtcars$mpg)
library(ggplot2)
qplot(mtcars$mpg,binwidth=2)
cov(mtcars[1:3])
cor(mtcars[1:3])
library(reshape2)
qplot(x=Var1,y=Var2,data=melt(cor(mtcars[1:3])),fill=value,geom="tile")
melt(cor(mtcars[1:3]))
lmfit=lm(mtcars$mpg~mtcars$cyl)
lmfit
summary(lmfit)
plot(mtcars$cyl,mtcars$mpg)
abline(lmfit)
plot(mtcars$mpg,mtcars$cyl)
abline(lmfit)
?lm
binom.test(x=92,n=315,p=1/6)
boxplot(mtcars$mpg,mtcars$mpg[mtcars$am==0],ylab="mpg",names=c("overall","automobile"))
abline(h=mean(mtcars$mpg),lwd=2,col="red")
abline(h=mean(mtcars$mpg[mtcars$am==0]),lwd=2,col="blue")
mpg.mu=mean(mtcars$mpg)
mpg.am=mtcars$mpg[mtcars$am==0]
t.test(mpg.am,mu=mpg.mu)
boxplot(mtcars$mpg~mtcars$am,ylab="mpg",names=c("automatic","manual")
abline(h=mean(mtcars$mpg[mtcars$am==0]),lwd=2,col="blue")
abline(h=mean(mtcars$mpg[mtcars$am==1]),lwd=2,col="red")
t.test(mtcars$mpg~mtcars$am)
library(rms)
library(car)
install.packages("car")
olsfit=ols(log(wages)~age+sex+education,data=SLID,x=T,y=T)

summary(olsfit)
robcov(olsfit)
lmfit=lm(wages~age+sex+education,data=SLID)
summary(lmfit)
lmfit1=glm(wages~age+sex+education,family = gaussian,data=SLID)
summary(lmfit1)
lmfit2=lm(wages~age+sex+education,data=SLID)
summary(lmfit2)
anova(lmfit1,lmfit2)
data("warpbreaks")
head(warpbreaks)
rs1=glm(breaks~tension,data = warpbreaks,family = "poisson")
summary(rs1)
table(warpbreaks$tension)
head(mtcars$vs)
lm1=glm(vs~hp+mpg+gear,data=mtcars,family = "binomial")
summary(lm1)
library(mgcv)
library(MASS)
attach(Boston)
str(Boston)
fit=gam(dis~s(nox))
summary(fit)
?gam
plot(nox,dis)
x=seq(0,1,length=500)
y=predict(fit,data.frame(nox=x))
lines(x,y,col="red",lwd=2)
plot(fit)
fit2=gam(medv~crim+zn+crim:zn,data=Boston)
vis.gam(fit2)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
data=read.table("clipboard",header = F,sep='\t')
head(data)
sym=select(org.Hs.eg.db,keys=data2,columns = c("ENSEMBL","SYMBOL","GENENAME"),keytype = "ENSEMBL")
sym
sum(table(sym$SYMBOL)==1)
head(sym)
head(data)
data1=as.vector(data)
data2=as.character(as.matrix(data1))
head(data2)
head(sym)
write.csv(sym,file="sym.csv")
getwd()


attach(mtcars)
plot(wt,mpg)
abline(lm(mpg~wt))
lm(mpg~wt)
title("Regression of MPG on Weight")
detach(mtcars)
dose=c(20,30,40,45,60)
drugA=c(16,20,27,40,60)
drugB=c(15,18,25,31,40)
plot(dose,drugA,type="b")
opar=par(no.readonly = T)
par(lty=2,pch=17)
plot(dose,drugA,type="b")
par(opar)
plot(dose,drugA,type="b",lty=3,lwd=3,pch=15,cex=2)
opar=par(no.readonly = T)
par(pin=c(2,3))
par(lwd=2,cex=1.5)
par(cex.axis=.75,font.axis=3)
plot(dose,drugA,type="b",pch=19,lty=2,col="red")
plot(dose,drugB,type="b",pch=23,lty=6,col="blue",bg="green")
par(opar)
plot(dose,drugA,type="b",col="red",lty=2,pch=2,lwd=2,main="Clinical Treials for Drug A",sub="This is hypothetical data",xlab="Dosage",ylab="Drug Response",xlim=c(0,60),ylim=c(0,70))
hist(mtcars$wt,main="Histogram of wt")
boxplot(mtcars$wt,mian="boxplot of wt")
attach(mtcars)
opar=par(no.readonly = T)
par(mfrow=c(3,1))
hist(wt)
hist(mpg)
hist(disp)
par(opar)
detach(mtcars)
attach(mtcars)
layout(matrix(c(1,1,2,3),2,2,byrow=T))
hist(wt)
hist(mpg)
boxplot(disp)
par(opar)
layout(matrix(c(1,1,2,3),2,2,byrow=T))
hist(mpg)
hist(wt)
hist(disp)
detach(mtcars)
set.seed(1)
population = unlist( read.csv("femaleControlsPopulation.csv") )
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
library(devtools)
install_github("genomicsclass/dagdata")
library(dagdata)
dir=system.file(package = "dagdata")
dir
list.files(file.path(dir,"extdata"))
filename=file.path(dir,"extdata/femaleControlsPopulation.csv")
filename
population=unlist(read.csv(filename))
set.seed(1)
N=12
B=10000
pvals=replicate(B,{
  control=sample(population,N)
  treatment=sample(population,N)
  t.test(control,treatment)$p.val
  
  
})

hist(pvals)
library(GSE5859Subset)
data(GSE5859Subset)
g=sampleInfo$group
g
e=geneExpression[25,]
library(rafalib)
mypar(1,2)
qqnorm(e[g==1])
qqline(e[g==1])
e
qqnorm(e[g==0])
qqline(e[g==0])
t.test(e[g==1],e[g==0])$p.val
myttest=function(x) t.test(x[g==1],x[g==0],var.equal = T)$p.value
pvals=apply(geneExpression,1,myttest)
hist(pvals)
sum(pvals<0.05)
library(genefilter)
prowtest=rowttests(geneExpression,factor(g))
?rowttests
sum(prowtest$p.value<0.05)

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
library(BiocInstaller)
biocLite(c("GO.db", "preprocessCore", "impute"))
install.packages("WGCNA")
library(WGCNA)
install.packages("scales")
library(dynamicTreeCut)
library(lazyeval)
update.packages("lazyeval")
library(lazyeval)
install.packages("lazyeval")
remove.packages("lazyeval")
install.packages("lazyeval")
library(lazyeval)
library(WGCNA)
library(GSE5859Subset)
data(GSE5859Subset)
dim(geneExpression)
dim(geneAnnotation)
head(geneExpression)
head(geneAnnotation)
dim(sampleInfo)
head(sampleInfo)
all.equal(sampleInfo$filename,colnames(geneExpression))
all.equal(geneAnnotation$PROBEID,rownames(geneExpression))
options(digits = 2)
cbind(sampleInfo[1:3,],colnames(geneExpression)[1:3],t(geneExpression)[1:3,1:4])
rownames(sampleInfo)=sampleInfo$filename
head(sampleInfo$filename)
rownames(geneAnnotation)=geneAnnotation$PROBEID
head(geneAnnotation[,1])
library(Biobase)
es5859=ExpressionSet(assayData = geneExpression,)
es5859
pData(es5859)=sampleInfo
fData(es5859)=geneAnnotation
es5859[which(fData(es5859)$CHR=="chrY"
             methods(class="ExpressionSet")
library(annotate)           
mi=pmid2MIAME("17206142")
experimentData(es5859)=mi
es5859
annotation(es5859)="hgfocus.db"
es5859
abstract(es5859)
substr(abstract(es5859),1,50)
library(GEOquery)
glioMA=getGEO("GSE78703")[[1]]
glioMA
head(pData(glioMA))
fData(glioMA)
library(BiocInstaller)
biocLite("")

library(GEOquery)
gset=getGEO("GSE13535",GSEMatrix = T,AnnotGPL = T)
get(gset)
show(gset)
fdata=fData(gset[[1]])
head(fdata)
dim(fdata)
colnames(fdata)

exprSet=(exprs(gset[[1]]))
pdata=pData(gset[[1]])
head(pdata,n=1)
view(pdata)
library(DT)
view(pdata)
??view
library(grid)
View(pdata)
sample=pdata$geo_accession
treat_time=rep(c('2h','18h'),each=11)
nrow(sample)
dim(sample)
table(sample)
length(sample)

treat_type=rep(rep(c("vehicle_control","PE1.3_embolized","PE2.0_embolized"),c(3,4,4)),times=2)
treat_type
design_df=data.frame(sample,treat_time,treat_type)
design_df
TS=paste(design_df$treat_time,design_df$treat_type,sep=".")
TS
TS=factor(TS,levels = unique(TS))
TS
unique(TS)
design=model.matrix(~0+TS)
design
View(design)
fit=lmFit(exprSet,design)
library(limma)
fit
limmaUsersGuide(view=TRUE)
design
cont.matrix=makeContrasts(vs1=TS18h.vehicle_control-TS2h.vehicle_control,vs2=TS18h.PE2.0_embolized-TS2h.PE2.0_embolized,vs3=TS18h.PE1.3_embolized-TS2h.PE1.3_embolized,diff=(TS18h.PE2.0_embolized-TS18h.vehicle_control)-(TS18h.PE1.3_embolized-TS18h.vehicle_control),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2
results=decideTests(fit2)
vennDiagram(results)
fit3=ebayes(fit2)
top=topTable(fit3)
cont.matrix1=makeContrasts(TS18h.vehicle_control-TS2h.vehicle_control,TS18h.PE2.0_embolized-TS2h.PE2.0_embolized,diff=(TS18h.PE2.0_embolized-TS18h.vehicle_control)-(TS18h.PE1.3_embolized-TS18h.vehicle_control),levels = design)
fite=contrasts.fit(fit,cont.matrix1)
fitee=eBayes(fite)
topgene=topTable(fitee,coef = 1,adjust.method = "BH")
colnames(design)
design
topgene
colnames(design)
View(topgene)

summary(fit)
fit$coefficients
show(fit)
?eBayes
?topTable
names=rownames(topgene)
gene_symbol=fdata[names,]$`Gene symbol`
gene_symbol
?makeContrasts
targets=c("MU","MU","WT","WT","WT")
targets
factor(targets,levels=c("WT","MU"))
model.matrix(~0+targets)


library(tissuesGeneExpression)
data(tissuesGeneExpression)
d=dist(t(e))
library(rafalib)
mypar()
hc=hclust(d)
hc
plot(hc,labels=tissue,cex=0.5)
myplclust(hc,labels = tissue,lab.col = as.fumeric(tissue),cex=0.5)
myplclust(hc,labels=tissue,lab.col = as.fumeric(tissue),cex=0.5)
abline(h=120)
hclusters=cutree(hc,h=120)
table(true=tissue,cluster=hclusters)
hclu=cutree(hc,k=8)
table(true=tissue,cluster=hclu)

set.seed(1)
km=kmeans(t(e[1:2,]),centers = 7)
km
names(km)
mypar(1,2)
plot(e[1,],e[2,],col=as.fumeric(tissue),pch=16)
plot(e[1,],e[2,],col=km$cluster,pch=16)
table(true=tissue,cluster=km$cluster)
km=kmeans(t(e),centers=7)
mds=cmdscale(d)
?cmdscale

plot(mds[,1],mds[,2])
plot(mds[,1],mds[,2],col=km$cluster,pch=16)
table(true=tissue,cluster=km$cluster)
library(RColorBrewer)
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol
library(genefilter)
rv=rowVars(e)
idx=order(-rv)[1:40]
library(gplots)
cols=palette(brewer.pal(8,"Dark2"))[as.fumeric(tissue)]
head(cbind(colnames(e),cols))
heatmap.2(e[idx,],labCol = tissue,trace="none",ColSideColors = cols,col=hmcol)
library(Biobase)
library(SpikeIn)
library(BiocInstaller)
biocLite("hgu95acdf")

n
library(hgu95acdf)
data("SpikeIn95")
library(rafalib)
mypar()
plot(X,Y)
detach("package:AnnotationDbi")

library(BiocInstaller)
biocLite("cometExactTest")
source("https://bioconductor.org/biocLite.R")
remove.packages("BiocInstaller")
source("https://bioconductor.org/biocLite.R")
library(maftools)
install.packages("cometExactTest")
biocLite()
install.packages("ggplot2")
data(mtcars)
df=mtcars[,c("mpg","cyl","wt")]
df$cyl=as.factor(df$cyl)
head(df)
qplot(x,y=null,data,geom="auto")
library(ggplot2)
qplot(x=mpg,y=wt,data=df,geom=c("point","smooth"))
qplot(x=mpg,y=wt,data=df,color=cyl,shape=cyl)
set.seed(1234)
wdata=data.frame(sex=factor(rep(c("F","M"),each=200)),weight=c(rnorm(200,55),rnorm(200,58)))
head(wdata)
boxplot(weight~sex,wdata)
qplot(sex,weight,data=wdata,geom="boxplot",fill=sex)
qplot(sex,weight,data=wdata,geom="violin")
qplot(sex,weight,data=wdata,geom="dotplot",stackdir="center",binaxis="y",dotsize=0.5,color=sex)
qplot(weight,data=wdata,geom="histogram",fill=sex)
qplot(weight,data=wdata,geom="density",color=sex,linetype=sex)
qplot(x=wt,y=mpg,data=df,geom="point")
ggplot(data=df,aes(x=mpg,y=wt))+geom_point()
ggplot(data=df,aes(x=mpg,y=wt))+geom_point(color="blue",size=2,shape=23)
ggplot(wdata,aes(x=weight))+geom_density()
ggplot(wdata,aes(x=weight))+stat_density()
library(plyr)
mu=ddply(wdata,"sex",summarise,grp.mean=mean(weight))
mu
?apply
head(wdata)
colnames(wdata)
head(wdata$sex)
subset(wdata,sex="F",select=weight)
biocLite("affy")
3


library(ggplot2)
data("mtcars")
qplot(x=mpg,y=wt,data=df,colour=cyl,shape=cyl)
qplot(x=mpg, y=wt, data = df, colour=cyl, shape=cyl)
df <- mtcars[, c("mpg","cyl","wt")]
df$cyl <- as.factor(df$cyl)
qplot(x=mpg,y=wt,data=df,geom=c("point","smooth"))
set.seed(1234)
head(wdata)
wdata <- data.frame(
  sex=factor(rep(c("F", "M"), each=200)),
  weight=c(rnorm(200, 55), rnorm(200, 58))
)
head(wdata)
qplot(sex,weight,data=wdata,geom="boxplot",fill=sex)
qplot(sex,weight,data=wdata,geom="violin")
qplot(sex,weight,data=wdata,geom="dotplot",stackdir="center",binaxis="y",dotsize=0.5,color=sex)
qplot(weight,data=wdata,geom="histogram",fill=sex)
qplot(weight,data=wdata,geom="density",color=sex,linetype=sex)
ggplot(data=df,aes(x=mpg,y=wt))+geom_point()
ggplot(data=df,aes(x=mpg,y=wt))+geom_point(color="blue",size=2,shape=23)
ggplot(wdata,aes(x=weight))+geom_density()
ggplot(wdata,aes(x=weight))+stat_density()
library(plyr)
mu=ddply(wdata,"sex",summarise,grp.mean=mean(weight))
mu
a=ggplot(wdata,aes(x=weight))
a
a+geom_area()
a+geom_density()
a+geom_histogram()
library(ggplot2)
ggplot(data=df,aes(x=mpg,y=wt))+geom_point()
ggplot(data=df,aes(x=mpg,y=wt))+geom_point(color="blue",size=2,shape=23)
ggplot(wdata,aes(x=weight))+geom_density()
ggplot(wdata,aes(x=weight))+stat_density()
a+geom_area(aes(y=..density..),stat="bin")
a+geom_density()
a+geom_density(aes(color=sex))
a+geom_density(aes(fill=sex),alpha=0.4)
a+geom_density(aes(color=sex))+geom_vline(data=mu,aes(xintercept=grp.mean,color=sex),linetype="dashed")
a+geom_dotplot(aes(fill=sex))+scale_fill_manual(values=c("#999999","#E69F00"))
a+geom_freqpoly(aes(color=sex))
a+geom_freqpoly(aes(color=sex,linetype=sex))+theme_minimal()
a+geom_histogram()
a+geom_histogram(aes(color=sex),fill="white",position="dodge")+theme_classic()
a+stat_ecdf()
ggplot(data=mtcars,aes(sample=mpg))+stat_qq()
data(mpg)
b=ggplot(mpg,aes(x=fl))
b+geom_bar(fill="steelblue",color="black")+theme_classic()
head(mpg)
rownames(mpg)
colnames(mpg)
b=ggplot(data=mtcars,aes(x=wt,y=mpg))
b+geom_point()
b+geom_smooth()
b+geom_point()
b+geom_point(aes(color=factor(cyl),shape=factor(cyl)))+scale_color_manual(values=c("#999999","#E69F00","#56B4E9"))+theme_classic()
b+geom_point()+geom_smooth(method="loess",se=F)
b+geom_point(aes(color=factor(cyl),shape=factor(cyl)))+geom_smooth(aes(color=factor(cyl),shape=factor(cyl)),method="lm",se=F,fullrange=T)
ggplot(data=mpg,aes(cty,hwy))+geom_point()+geom_quantile()+theme_minimal()
install.packages("quantreg")
ggplot(data=faithful,aes(x=eruptions,y=waiting))+geom_point()+geom_rug()
p=ggplot(data=mpg,aes(displ,hwy))
p+geom_point()
p+geom_jitter(width=0.5,height=0.5)
b+geom_text(aes(label=rownames(mtcars)))
library(Biostrings)
library(BiocInstaller)
biocLite("Biostrings")
n
letters
LETTERS
sample(LETTERS[c(1,3,7,20)],size=10,replace=T)
DNA_ALPHABET
seq=sample(DNA_ALPHABET[1:4],size=10,replace = T)
seq
str(seq)
paste(seq,collapse = "")
bstring=BString("I am a BString object")
bsstring
bstring
length(bstring)
dnastring=DNAString("TTGCAA-CTC-N")
dnastring
dnastring
length(dnastring)
str(dnastring)
slotNames(dnastring)
alphabetFrequency(dnastring)
alphabetFrequency(dnastring,baseOnly=T,as.prob=T)
letterFrequency(dnastring,"A",as.prob=T)

reverseComplement(dnastring)
dnastring[2]
dnastring[7:12]
view=Views(dnastring,start=3:0,end=5:8)
view
length(view)
view[4:2]
view[[2]]
set=NULL
for(i in 1:4)
  set=c(set,paste(sample(DNA_ALPHABET[1:4],10,replace = T),collapse = ""))
set
length(set)
reverseComplement(set)
alphabetFrequency(set,baseOnly=T,as.prob=T)

letterFrequency(set,"C",as.prob=T)
width(set)
subset(set,1,3)
subseq(set,1,3)

library(rmarkdown)
install.packages("rmarkdown")
library(car)
install.packages
install.packages("viridis")
 library(viridis)
plot(pressure)
library(BiocInstaller)
biocLite("GenomicRanges")
library(GenomicRanges)
exon=GRanges(seqnames=c("chr1"),strand="+",IRanges(start=c(1:5),end=c(6:11)))
exon
names(exon)=c("a","a","a","b","c","c")
exon
b=c("a","a","a","b","c","c")
mcols(exon)=b
unique(b)
exon
?unique
width(reduce(exon[exon$X=="a",]))
length=data.frame(name=unique(b))
length=data.frame()
for(i in c("a","b","c")){
length[i,"length"]=width(reduce(exon[exon$X==i,]))}

length
dim(length)



data(iris)
table.iris=table(iris$Species)
table.iris
pie(table.iris)
hist(iris$Sepal.Length)
boxplot(Petal.Width~Species,data=iris)
plot(x=iris$Sepal.Length,y=iris$Sepal.Width,col=iris$Species)
install.packages("C50")
data(churn)
library(C50)
str(churnTrain)
churnTrain=churnTest[,!names(churnTrain)%in%c("state","area_code","account_length")]
set.seed(2)
ind=sample(2,nrow(churnTrain),replace = T,prob=c(0.7,0.3))
table(ind)
trainset=churnTrain[ind==1,]
testset=churnTrain[ind==2,]
install.packages("e1071")
library(e1071)
model=svm(churn~.,data=trainset,kernel="radial",cost=1,gamma=1/ncol(trainset))
model
summary(model)
colnames(trainset)
dim(trainset)
svm.pred=predict(model,testset)
summary(svm.pred)
svm.table=table(svm.pred,testset$churn)
svm.table
install.packages("neuralnet")
library(neuralnet)
churn_neural=neuralnet(yes+no~number_vmail_messages+total_day_minutes+total_day_calls+total_day_charge+total_eve_calls,data=trainset,hidden=2)

?neuralnet()
colnames(trainset)
head(trainset)
colnames(iris)
table(iris$Species)
table(trainset$churn)
trainset$yes=trainset$churn=="yes"
trainset$yes
trainset$no=trainset$churn=="no"

churn_neural$result.matrix
net.predict=compute(churn_neural,testset[,c(3,4,5,6,7)])
head(net.predict$net.result)
index=apply(net.predict$net.result,1,which.max)
head(index)
index
table(index==2)
net_pred=c("yes","no")[index]
table(net_pred)
table(testset$churn,net_pred)
library(scater)
library()