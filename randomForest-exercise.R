url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data"

data <- read.csv(url, header=FALSE)
View(data)

# 取回的数据没有列名
head(data)
# 手动修改列名
colnames(data) <- c(
  "age",  "sex",  "cp",   "trestbps", 
  "chol",   "fbs",    "restecg", 
  "thalach",   "exang",    "oldpeak", 
  "slope",  "ca",   "thal",  "hd" 
) # 共14个

head(data)
str(data)
data[data=="?"]=NA
data$sex=factor(data$sex,levels = c(0,1),labels = c("F","M"))
head(data$sex)

data$cp=as.factor(data$cp)
data$fbs=as.factor(data$fbs)
data$restecg=as.factor(data$restecg)
data$exang=as.factor(data$exang)
data$slope=as.factor(data$slope)

str(data$ca)
table(data$thal)
View(data$thal)
data$ca=as.integer(data$ca)
data$ca=as.factor(data$ca)
data$thal=as.integer(data$thal)
data$thal=as.factor(data$thal)

data$hd=ifelse(data$hd==0,"Healthy","Unhealthy")
data$hd=as.factor(data$hd)

library(randomForest)
data.imputed=rfImpute(hd~.,data=data,iter=6)
model=randomForest(hd~.,data=data.imputed,proximity=T)

oob.error.data=data.frame(trees=rep(1:nrow(model$err.rate),times=3),Type=rep(c("OOB","Healthy","Unhealthy"),each=nrow(model$err.rate)),Error=c(model$err.rate[,"OOB"],model$err.rate[,"Healthy"],model$err.rate[,"Unhealthy"]))
View(oob.error.data)                        

library(ggplot2)
ggplot(oob.error.data,aes(x=trees,y=Error))+geom_line(aes(color=Type))

oob.values=vector(length=10)
for(i in 1:10){
  temp.model=randomForest(hd~.,data=data.imputed,mtry=i)
  oob.values[i]=temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values
