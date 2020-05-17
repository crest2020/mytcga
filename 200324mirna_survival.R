options(stringsAsFactors = F)
library(RTCGA.miRNASeq)
BRCA.miRNASeq[1:10,1:6]
expr=expressionsTCGA(BRCA.miRNASeq)
dim(expr)
expr[1:10,1:6]
expr1<-expr[seq(1,nrow(expr),by=3),3:ncol(expr)]
dim(expr1)
expr1[1:6,1:6]
str(expr1)
expr1<-apply(expr1,2,as.numeric)
expr1[1:6,1:6]
rownames(expr1)=rownames(BRCA.miRNASeq)[seq(1,nrow(BRCA.miRNASeq),by=3)]
expr1[1:6,1:6]
expr1<-t(expr1)
dim(expr1)
table(apply(expr1,1,is.na)) 
Expr<-na.omit(expr1)
table(apply(Expr,1,is.na))
dim(Expr)
Expr<-Expr[apply(Expr, 1, function(x){
  sum(x>1)>10}),]
dim(Expr)
Expr[1:6,1:6]

library(RTCGA.clinical)
BRCA_clinical<-BRCA.clinical
dim(BRCA_clinical)
View(BRCA_clinical[1,])
write.csv(BRCA_clinical,quote = F,"BRCA_clinical_Data.csv")
BRCA_clinicaldata<-BRCA_clinical[c('patient.bcr_patient_barcode',
                                   'patient.vital_status',
                                   'patient.days_to_death',
                                   'patient.days_to_last_followup',
                                   'patient.race',
                                   'patient.age_at_initial_pathologic_diagnosis',
                                   'patient.gender',
                                   'patient.stage_event.pathologic_stage'
)]
str(BRCA_clinicaldata)
head(BRCA_clinicaldata)
BRCA_clinicaldata[1:4,1:4]
rownames(BRCA_clinicaldata) <- NULL
head(BRCA_clinicaldata)
BRCA_clinicaldata<-tibble::column_to_rownames(BRCA_clinicaldata,var="patient.bcr_patient_barcode")
head(BRCA_clinicaldata) 
head(Expr)
rownames(BRCA_clinicaldata)<-stringr::str_to_upper(rownames(BRCA_clinicaldata))
head(Expr)
head(BRCA_clinicaldata)
dim(BRCA_clinicaldata)
dim(Expr)
Expr[1:4,1:4]
save(Expr,BRCA_clinicaldata,file="step1_TCGA_BRCA_miRNAexpr_clinical.Rdata")

library(edgeR)
library(stringr)
group_list<-ifelse(as.numeric(substr(colnames(Expr),14,15))<10, 'Tumor','Normal')
group_list <- factor(group_list,levels = c("Normal","Tumor"))
group_list
table(group_list)
dge<-DGEList(counts=Expr,group=group_list)
keep <- rowSums(cpm(dge)>1) >= 3
table(keep)
dge <- dge[keep,keep.lib.sizes=FALSE]
dge
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge, method = 'TMM')
dge$samples
logCPM <- cpm(dge, log=TRUE, prior.count=2)
design<-model.matrix(~0+group_list)
colnames(design)<-levels(group_list)
rownames(design)<-colnames(dge)
head(design)
dge <-estimateDisp(dge,design)
fit <- glmFit(dge, design) 
fit2 <- glmLRT(fit, contrast=c(-1,1))
allDEM<-topTags(fit2, n=nrow(dge))
allDEM<-as.data.frame(allDEM)
head(allDEM)
write.csv(allDEM,file="allDEM.csv") 
FDR_cutoff<-0.05
logFC_cutoff<-1
allDEM$change<-as.factor(ifelse(allDEM$FDR>FDR_cutoff,'No', #如果FDR>0.05，则定义为'No'，否则进行下一步的ifelese
                                ifelse(allDEM$logFC>logFC_cutoff,'Up', #在FDR<=0.05前提下，logFC>0.05则定义'Up',否则进行下一步的ifelse
                                       ifelse(allDEM$logFC< -logFC_cutoff,'Down','No'))))
table(allDEM$change)
SigDEMname<-rownames(allDEM)[allDEM$change !="No"]
SigDEM<-allDEM[SigDEMname,]
write.csv(SigDEM,file="SigDEM.csv") #输出数据
Up_DEMname<-rownames(allDEM)[allDEM$change=="Up"]


Up_DEM<-allDEM[Up_DEMname,]
write.csv(Up_DEM,file="Up_DEM.csv")
Down_DEMname<-rownames(allDEM)[allDEM$change=="Down"]
Down_DEM<-allDEM[Down_DEMname,]
write.csv(Down_DEM,file="Down_DEM.csv")
head(allDEM)
Expr['hsa-mir-133a-1',]
library(ggplot2)
volcano_plot <- ggplot(data = allDEM,
                       aes(x = logFC,
                           y = -log10(FDR))) +
  geom_point(alpha=0.4, size=3.5,
             aes(color=change))
volcano_plot

library(survival)
library(survminer)
Up_DEM_expr<-Expr[rownames(Up_DEM),]
Up_DEM_expr[1:4,1:4]

dim(Up_DEM_expr)
dim(BRCA_clinicaldata)
table(group_list)
library("FactoMineR")
library("factoextra")

Exprset_PCA=t(Up_DEM_expr)
Exprset_PCA<-as.data.frame(Exprset_PCA)
Exprset_PCA.pca <- PCA(Exprset_PCA, graph = F)
Exprset_PCA=cbind(Exprset_PCA,group_list)
fviz_pca_ind(Exprset_PCA.pca,
             title = "Principal Component Analysis",
             geom.ind = "point", ###  show points only (nbut not "text"),这里c("point", "text)2选1
             col.ind = Exprset_PCA$group_list, ###  color by groups
             ###  palette = c("#00AFBB", "#E7B800"),自定义颜色
             addEllipses = TRUE, ###  Concentration ellipses加圆圈
             legend.title = "Groups")
head(Exprset_PCA.pca$ind)
dim(Exprset_PCA.pca$eig)

class(Exprset_PCA.pca)
(Exprset_PCA.pca$svd[1])
length(Exprset_PCA.pca$var)
View(Up_DEM_expr[1:4,1:4])
head(group_list)
Up_DEM_expr<-Up_DEM_expr[,group_list=='Tumor']
dim(Up_DEM_expr)
dim(BRCA_clinicaldata)
View(BRCA_clinicaldata[1,]) 
colnames(BRCA_clinicaldata)<-c("event","days_to_death","days_to_last_followup","race","age","gender","stage")
head(BRCA_clinicaldata)
head(Up_DEM_expr)
library(stringr)
BRCA_clinicaldata=BRCA_clinicaldata[match(substr(colnames(Up_DEM_expr),1,12),rownames(BRCA_clinicaldata)),]
dim(Up_DEM_expr)
dim(BRCA_clinicaldata)
head(colnames(Up_DEM_expr))
head(BRCA_clinicaldata)
table((rownames(BRCA_clinicaldata))==(substr(colnames(Up_DEM_expr),1,12)))
rownames(BRCA_clinicaldata)=colnames(Up_DEM_expr)
str(BRCA_clinicaldata)
View(BRCA_clinicaldata[1:4,])
BRCA_clinicaldata[,2][is.na(BRCA_clinicaldata[,2])]<-0
BRCA_clinicaldata[,3][is.na(BRCA_clinicaldata[,3])]<-0
str(BRCA_clinicaldata)
BRCA_clinicaldata$time_days<-as.numeric(BRCA_clinicaldata[,2])+as.numeric(BRCA_clinicaldata[,3])
View(BRCA_clinicaldata[1:4,])
table(BRCA_clinicaldata$time_days<= 0)
BRCA_clinicaldata<-BRCA_clinicaldata[BRCA_clinicaldata$time_days>0,]  
table(BRCA_clinicaldata$time_days<= 0)
BRCA_clinicaldata$time_year<-BRCA_clinicaldata$time_days/365
View(BRCA_clinicaldata[1:5,])
BRCA_clinicaldata$event<-ifelse(BRCA_clinicaldata$event=='alive',0,1)
str(BRCA_clinicaldata)
table(BRCA_clinicaldata$event)
dim(BRCA_clinicaldata)
dim(Up_DEM_expr)
Up_DEM_expr<-Up_DEM_expr[,rownames(BRCA_clinicaldata)]
dim(Up_DEM_expr)
View(Up_DEM_expr[1:4,1:4])
all(colnames(Up_DEM_expr)==rownames(BRCA_clinicaldata))

library(survival)
Surv1=with(BRCA_clinicaldata,Surv(time_year,event))
?Surv
Up_DEM_expr<-t(as.data.frame(Up_DEM_expr))
dim(Up_DEM_expr)
class(Up_DEM_expr)
Up_DEM_expr<-as.data.frame(Up_DEM_expr)
class(Up_DEM_expr)
res.cox=coxph(Surv1~Up_DEM_expr[,1],data = Up_DEM_expr)
res.cox
res.cox_summary=summary(res.cox)
res.cox_summary
View(BRCA_clinicaldata)
mySurv<-with(BRCA_clinicaldata,Surv(time_year,event))
View(Up_DEM_expr[1:4,1:4])
mySurv
View(mySurv[1:4,1:3])
univ_models=apply(Up_DEM_expr,2,function(gene){res.cox=coxph(mySurv~gene,data=Up_DEM_expr)})
options(scipen=200)
class(univ_models)
head(univ_models)
univcox_results=lapply(univ_models,function(x){
  x=summary(x)
  P.value=round(x$wald["pvalue"],digits = 4)
  wald.test=round(x$coef["test"],digits = 4)
  beta<-round(x$coef[1],digits=4)
  HR<-round(x$coef[2],digits=4)
  HR_CI_lower<-round(x$conf.int[,"lower .95"],4)
  HR_CI_upper<-round(x$conf.int[,"upper .95"],4)
  HR<-paste0(HR,"(",HR_CI_lower,"-",HR_CI_upper,")")
  z<-round(x$coef[4],digits=4)
  res<-c(beta, HR, z,wald.test, P.value)
  names(res)<-c("beta","HR (95% CI)","z","wald.test","P-value")
  return(res)
})
class(univcox_results)
results <- t(as.data.frame(univcox_results, check.names = FALSE))
class(results)
results<-as.data.frame(results)
View(results)
results<-results[order(results$`P-value`,decreasing = F),]
View(results)
write.csv(results,file="univariate COX analysis_log2(x+1).csv")

options(stringsAsFactors = F)

length(Up_univariate.res) 
Up_univariate.res<-results[unclass(results$`P-value`)< 0.05,]
Up_univariate.res
head(as.numeric(results$`P-value`))
results=as.data.frame(results,stringsAsFactors=F)
View(results)
class(results)
?as.data.frame
d=data.frame(a=c(1,3,4),b=c(3,4,5))
View(d)
as.character(results$`P-value`)
Up_univariate.res<-results[as.numeric(as.character(results$`P-value`))< 0.05,]
dim(Up_univariate.res)
View(Up_univariate.res)
multi_input.expr<-Up_DEM_expr[,rownames(Up_univariate.res)]
multi_input.expr[1:4,1:4]
mySurv<-with(BRCA_clinicaldata,Surv(time_year,event))
multi_COX<-coxph(mySurv ~ ., data=multi_input.expr)  
summary(multi_COX)
step.multi_COX=step(multi_COX,direction = "both")
step.multi_COX
RiskScore=predict(step.multi_COX,type="risk",newdata = multi_input.expr)
RiskScore
RiskGroup=ifelse(RiskScore>=median(RiskScore),"high","low")
View(RiskScore)
View(RiskGroup)
library(dplyr)
library(pheatmap)
ten_miRNA_names=c('hsa-mir-1293','hsa-mir-466', 'hsa-mir-3923','hsa-mir-2276','hsa-mir-3927', 'hsa-mir-9-3','hsa-mir-3929','hsa-mir-556','hsa-mir-3609','hsa-mir-20b')
miRNA_heatmap_input<-Up_DEM_expr[,ten_miRNA_names]
miRNA_heatmap_input[1:4,1:4]
risk_group<-as.data.frame(RiskGroup)
all(rownames(risk_group)==rownames(miRNA_heatmap_input))
heatmap.dat<-cbind(id=rownames(miRNA_heatmap_input),risk_group,miRNA_heatmap_input)
head(factor(heatmap.dat$RiskGroup,levels=c("low","high"),ordered = TRUE))
dim(heatmap.dat)
View(heatmap.dat[1:10,])
head(heatmap.dat$RiskGroup)
table(heatmap.dat$RiskGroup)
heatmap_order=(heatmap.dat[order(heatmap.dat$RiskGroup,decreasing = F),])
library(pheatmap)
miRNA_heatmap_input<-heatmap_order[,3:ncol(heatmap_order)]
miRNA_heatmap_input<-t(miRNA_heatmap_input)
ten_miRNA_heatmap<-pheatmap(log2(miRNA_heatmap_input+1), #准备的数据
                            color = colorRampPalette(c("green", "black", "red"))(50), #这里用这个颜色和文章颜色差不多
                            annotation_col = risk_group, #分组注释
                            show_colnames =FALSE, #不显示列名
                            show_rownames = TRUE, #显示行名
                            cluster_cols=FALSE,  #这里不对列聚类，因为我们已经自己分好组了
                            cluster_rows=TRUE)   #对行聚类
all(substr(names(RiskScore),1,12)==substr(rownames(BRCA_clinicaldata),1,12))
all(substr(names(RiskGroup),1,12)==substr(rownames(BRCA_clinicaldata),1,12))
dim(BRCA_clinicaldata)
length(RiskGroup)
KM.input<-cbind(BRCA_clinicaldata[,c("event","time_year")],RiskScore,RiskGroup)
str(KM.input)
fit<-survfit(Surv(time_year,event) ~RiskGroup , data=KM.input)
print(fit) 
summary(fit)
library(survminer)
KMsurvival_plot<-ggsurvplot(fit,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            conf.int = TRUE, #添加置信区间
                            conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = FALSE,###  show bars instead of names in text annotations in legend of risk table.不显示注释名字
                            xlab = "Time in years",   ###  customize X axis label.自定义x的标签为time in years
                            surv.median.line = "hv", #添加中位生存时间的线
                            ncensor.plot = FALSE, #我这里不显示删失的图，TRUE就显示
                            legend.labs =
                              c("high risk", "low risk"),    ###  对legend的标签重新命名
                            palette = c("#E7B800", "#2E9FDF"), ###  自定义颜色
                            ggtheme = theme_light()#绘图主题
)
KMsurvival_plot
KMres.sum  <- surv_summary(fit)
head(KMres.sum)

surv_diff <- survdiff(Surv(time_year, event) ~ RiskGroup, data = KM.input)
surv_diff
install.packages("survivalROC")
library(survivalROC)
Survival_ROC_input<-KM.input
survival_ROC=survivalROC(Stime = Survival_ROC_input$time_year,status = Survival_ROC_input$event,marker = Survival_ROC_input$RiskScore,predict.time = 5,method = "KM")
survival_ROC_AUC<-round(survival_ROC$AUC,3)
plot(survival_ROC$FP,survival_ROC$TP)
plot(survival_ROC$FP,survival_ROC$TP,type="l",xlim=c(0,1),ylim=c(0,1),
     xlab="False positive rate",  
     ylab="True positive rate",
     main=paste0("ROC Curve", " (", "AUC = ",survival_ROC_AUC," )"),  #标题
     cex.main=1.5,#放大标题
     cex.lab=1.3,#坐标轴标签（名称）的缩放倍数
     cex.axis=1.3, font=1.2, #放大轴上的数字
     lwd=1.5, #指定线条宽度
     col="red"  #红色
)
abline(a=0,b=1)
install.packages("timeROC")
library(timeROC)
time_ROC_input<-KM.input
time_ROC<-timeROC(T=time_ROC_input$time_year, #生存时间(dead和alive的生存时间).
                  delta=time_ROC_input$event, #生存结局，Censored的样本必须用0表示
                  marker=time_ROC_input$RiskScore, #预测的变量，这里是风险评分，在没有特殊情况下，值越大，越容易发生危险事件
                  cause=1, #阳性结局的赋值（必须是1或2），也就是dead的赋值，这里dead是1表示的
                  weighting = "marginal", #权重计算方法，这是默认方法，选择KM计算删失分布，weighting="aalen" [选用COX]，weighting="aalen" [选用Aalen]
                  times = c(3,5,10), #计算3、5、10年的ROC曲线
                  ROC=TRUE,
                  iid=TRUE #计算AUC
)
time_ROC
summary(time_ROC)
time_ROC$AUC
plot(time_ROC,time=3,col="red",title=FALSE,lwd=2)
plot(time_ROC,time=5,col="blue",add=TRUE,title=FALSE,lwd=2)
plot(time_ROC,time=10,col="black",add=TRUE,title=FALSE,lwd=2)
legend("bottomright", #图例设置在右下角
       c(paste0("AUC at 3 years = ", round(time_ROC$AUC[1],3)),
         paste0("AUC at 5 years = ", round(time_ROC$AUC[2],3)),
         paste0("AUC at 10 years = ", round(time_ROC$AUC[3],3))),
       col=c("red","blue","black"),lwd=2,bty="n" #或者bty+“o"
)
library(ggplot2)
time_ROC$TP
summary(time_ROC)
time_ROC.res<-data.frame(TP_3year=time_ROC$TP[,1], #获取3年的ROC的TP
                         FP_3year=time_ROC$FP[,1],  #获取3年的ROC的FP
                         TP_5year=time_ROC$TP[,2],  #获取5年的ROC的TP
                         FP_5year=time_ROC$FP[,2], #获取5年的ROC的FP
                         TP_10year=time_ROC$TP[,3], #获取10年的ROC的TP
                         FP_10year=time_ROC$FP[,3]) #获取10年的ROC的FP
time_ROC$AUC
TimeROC_plot<-ggplot()+
  geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=1,color="red")+
  geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=1,color="blue")+
  geom_line(data=time_ROC.res,aes(x=FP_10year,y=TP_10year),size=1,color="black")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2 #添加虚线
  )+
  theme_bw()+
  annotate("text",x=0.75,y=0.25,size=4.5,
           label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[1]],3)),color="red")+
  annotate("text",x=0.75,y=0.15,size=4.5,
           label=paste0("AUC at 5 years = ",round(time_ROC$AUC[[2]],3)),color="blue")+
  annotate("text",x=0.75,y=0.05,size=4.5,
           label=paste0("AUC at 10 years = ",round(time_ROC$AUC[[3]],3)),color="black")+
  labs(x="False positive rate",y="True positive rate")+
  theme(axis.text=element_text(face="bold", size=11,  color="black"),#加粗刻度标签
        axis.title=element_text(face="bold", size=14, color="black")) #加粗xy轴标签名字
TimeROC_plot
