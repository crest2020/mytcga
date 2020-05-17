head(phe)
exprset[1:4,1:4]
set.seed(1234567890)
dim(expr)
k=sample(1:593,300)
t_exp=expr[,k]
v_exp=expr[,-k]
head(t_phe)
x=t(log2(t_tumour+1))
y=t_phe$event
View(x)
all(t_phe$ID==colnames(t_exp))
dim(expr)
dim(phe)

data=cbind(x,y)
colnames(data)[ncol(data)]="event"
data=as.data.frame(data)
data$event=as.factor(data$event)
View(data)
library(e1071)
model=svm(formula=event~.,data=data,kernel="linear")
summary(model)
pred=predict(model,data)
table(pred,data$event)


x=t(log2(v_tumour+1))
y=v_phe$event
data=cbind(x,y)
colnames(data)[ncol(data)]="event"
data=as.data.frame(data)
data$event=as.factor(data$event)
pred=predict(model,data)
table(pred,data$event)
save.image(file="svm-predict.Rdata")
