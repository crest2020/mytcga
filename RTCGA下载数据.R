library(RTCGA)
library(RTCGA.miRNASeq)
browseVignettes("RTCGA.miRNASeq")
checkTCGA('Dates')
data(package="RTCGA.miRNASeq")
rownames(RTCGA.miRNASeq)
s=rownames(KIRC.miRNASeq)[seq(1,nrow(KIRC.miRNASeq),by=3)]
head(s)
expr=expressionsTCGA(KIRC.miRNASeq)
View(expr[,1:4])
mi=colnames(expr)
head(mi)
expr=apply(expr,)
dim(expr)

View(KIRC.miRNASeq[,1:3])
expr=as.data.frame(expr[seq(1,nrow(expr),by=3),3:ncol(expr)])
View(expr[,1:3])
mi=colnames(expr)
head(mi)
expr=apply(expr,1,as.numeric)
warnings()
#as.numeric矩阵发生了转置
colnames(expr)=s
rownames(expr)=mi
expr=na.omit(expr)
expr=expr[apply(expr,1,function(x){sum(x>1)>10}),]
dim(expr)
#临床信息
library(RTCGA.clinical)
meta=KIRC.clinical
tmp=as.data.frame(colnames(meta))
View(meta[1:10,1:10])
grepl(1,c(3,2,1))
head(meta[grepl("patient.bcr_patient_barcode",colnames(meta))])
head(meta[(grepl('patient.days_to_last_followup',colnames(meta)))])
head(meta[(grepl('patient.days_to_death',colnames(meta)))])
head( meta[(grepl('patient.vital_status',colnames(meta)))])
meta=as.data.frame(meta[c('patient.bcr_patient_barcode','patient.vital_status',
                          'patient.days_to_death','patient.days_to_last_followup',
                          'patient.race',
                          'patient.age_at_initial_pathologic_diagnosis',
                          'patient.gender' ,
                          'patient.stage_event.pathologic_stage')])
View(meta[1:10,])
save(expr,meta,file = "TCGA-KIRC-miRNA-example.Rdata")
