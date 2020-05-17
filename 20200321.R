a=read.table('GSE113143_Normal_Tumor_Expression.tab.gz',sep='\t',quote="",fill=T,comment.char = "!",header = T)
dim(a)
rownames(a)=a[,1]
a=a[,-1]
exp=a
fpkmtoTpm=function(fpkm){
  exp(log(fpkm)-log(sum(fpkm))+log(1e6))
}
tpms=apply(exp,2,fpkmtoTpm)
tpms[1:3,]
colSums(tpms)
group_list=c(rep(c("normal","tumour"),each=3))
?rep
group_list
group
##differential analysis
group_list=factor(group_list,levels=c("normal","tumour"),ordered = F)
group_list

exprset=tpms
boxplot(exprset,outline=F,notch=T,col=group_list,las=2)
exprset=log2(exprset+1)
dat=exprset
design=model.matrix(~factor(group_list))
design
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
warnings()
topTable(fit,coef=2,adjust="BH")
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p=ggboxplot(df,x="stage",y="gene")
}

deg=topTable(fit,coef=2,adjust="BH",number = Inf)
head(deg)
logFC_t=1.5

BiocManager::install("org.Mm.eg.db")
install.packages("org.Mm.eg.db_3.10.0.tar.gz",repos = NULL)
library(clusterProfiler)
library(org.Mm.eg.db)
deg$g=ifelse(deg$P.Value>0.05,'stable',ifelse(deg$logFC>logFC_t,'up',ifelse(deg$logFC<(-logFC_t),'down','stable')))
table(deg$g)
deg$symbol=rownames(deg)

library(ggplot2)
df=bitr(unique(deg$symbol),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
head(df)

DEG=merge(deg,df,by.x="symbol",by.y="SYMBOL")
head(DEG)
save(file="")
?save.image
gene_up=DEG[DEG$g=='up','ENTREZID']
gene_down=DEG[DEG$g=='down','ENTREZID']
enrich=enrichKEGG(gene=c(gene_down,gene_up),organism = "mmu")
head(enrichkk)[,1:6]
head(gene_up)
head(gene_down)
gene_down1=gene_down[1:1000]
length(gene_up)
length(gene_down)
browseKEGG(enrichkk, 'hsa04512')
dotplot(enrich
?enrichKEGG
enrichKK=DOSE::setReadable(enrich, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
enrichKK 
gene_down[1:1000]
gene_up
gene_diff=DEG[DEG$g!="stable",]
genelist=gene_diff$logFC
names(genelist)=gene_diff$ENTREZID
genelist=sort(genelist,decreasing = T)
head(genelist)
?cnetplot
cnetplot(enrich,foldChange = genelist,colorEdge=T,circular=T)
emapplot(enrich)
heatplot(enrich)
ego_bp=enrichGO(gene=DEG$ENTREZID,OrgDb = org.Mm.eg.db,keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.05)
goplot(ego_bp)

##RNAseq数据，下载GEO中的FPKM文件后该怎么下游分析
##https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247492196&idx=1&sn=1bcdef5e19cdde2b0d8cc947c759fac6&chksm=9b4ba2dfac3c2bc94329084cc201f937cf0c906ddc7e45a083aaad3c33e2f84b961577756b3d&mpshare=1&scene=24&srcid=&sharer_sharetime=1582970081206&sharer_shareid=c784b3dab32e04ddf5a7f3787407e657&key=b75e9a2bedf85391a26fac3d9c9ca9b6b1d645e5ecb03cd43e703aba91246c5ef7a60461fe599ddf8ec13936dbf0d0875b13e743a4aa3a6921520766e341af915bcacdef886d38f8933fb098f74511da&ascene=14&uin=MjE3NzczNTgyMQ%3D%3D&devicetype=Windows+7&version=6208006f&lang=zh_CN&exportkey=Ae0tPb4agovv9rV4njZ3QxI%3D&pass_ticket=hclFePS1EtLvIIr4jqXeCeKl1AJVGqcT8MFD9%2FrpY3dEQ1NF6JwtXahkdrpJUPmp


library(stringr)
libray(ggplot2)
barplot(ego_bp,showCategory = 16,title = "The go_bp enrichment analysis of all deg")+scale_size(range=c(2,12))+scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))
save.image(file="clusterProfiler.RData")

BiocManager::install("biomaRt")
library(biomaRt)
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
head(listDatasets(ensembl))
searchDatasets(mart=ensembl,pattern = "hsa.*")
mart=useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")



exp=data.frame(a=c(1,2,3),b=c(4,5,6),d=c(7,8,9))
exp
fun=function(i){return(data.frame(row=1:))}
i=2
data.frame(row=1:nrow(exp),col=i,e=exp[,i])

BiocManager::install("RTCGA")
BiocManager::install("RTCGA.miRNASeq")
BiocManager::install("RTCGA.clinical")
library(RTCGA.miRNASeq)
?miRNASeq
View(BRCA.miRNASeq[1:10,1:6])
expr=expressionsTCGA(BRCA.miRNASeq)

dim(expr)

expr1=expr[seq(1,nrow(expr),by=3),3:ncol(expr)]
dim(expr1)
str(expr1)
expr1=as.data.frame(expr1)
rownames(expr1)=rownames(BRCA.miRNASeq)[seq(1,nrow(expr),by=3)]

View(expr1[1:10,1:6])
expr1=apply(expr1,2,as.numeric)
expr1[1:6,1:6]
expr1=t(expr1)
dim(expr1)
table(apply(expr1,1,is.na))
expr2=na.omit(expr1)
dim(expr2)
table(apply(expr2,1,is.na))
expr2=expr2[apply(expr2,1,function(x){sum(x>1)>10}),]
dim(expr2)
View(expr2[1:10,1:6])
1046-677
##clinical

library(RTCGA.clinical)
brca_clinical=BRCA.clinical
dim(brca_clinical)
View(brca_clinical[1:4,])
brca_clinicaldata<-brca_clinical[c('patient.bcr_patient_barcode',
                                   'patient.vital_status',
                                   'patient.days_to_death',
                                   'patient.days_to_last_followup',
                                   'patient.race',
                                   'patient.age_at_initial_pathologic_diagnosis',
                                   'patient.gender',
                                   'patient.stage_event.pathologic_stage'
)]
str(brca_clinicaldata)
brca_clinicaldata[1:4,1:4]
rownames(brca_clinicaldata)=NULL
rownames(brca_clinicaldata)=brca_clinicaldata$patient.bcr_patient_barcode
View(brca_clinicaldata[1:6,1:6])
##differential analysis
BiocManager::install("edgeR")
library(edgeR)
library(stringr)
group_list=ifelse(as.numeric(substr(colnames(expr2),14,15))<10,"Tumour","Normal")
table(group_list)
head(rownames(expr2))

View(BRCA.miRNASeq)
?row.names



