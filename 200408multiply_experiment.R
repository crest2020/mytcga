BiocManager::install("gam")
library(BiocParallel)
library(SingleCellExperiment)
library(clusterExperiment)
library(scone)
library(zinbwave)
library(slingshot)

library(gam)
library(RColorBrewer)
library(devtools)
install_github("drisso/fletcher2017data")
install_github("drisso/fletcher2017data")
library(fletcher2017data)
data("fletcher")
fletcher2017data::fletcher
fletcher
detable=matrix(c(28,142,501,12000),nrow=2,dimnames = list(c("DE","Not.De"),c("In.gene.set","Not.In.gene.set")))
detable
fisher.test(detable,alternative = "greater")
BiocManager::install("ALL")
library(ALL)
data("ALL")
ind.bs=grep("^B",ALL$BT)
ind.bs
head(ind.bs)
head(ALL$BT)
ALL
ind.mut=which(ALL$mol.biol%in%factor(c("BCR/ABL","NEG")))
head(ind.mut)
factor(c("BCR/ABL","NEG"))
View(ALL$mol.biol)
ind.mut
sset=intersect(ind.bs,ind.mut)
all.eset=ALL[,sset]
dim(all.eset)
exprs(all.eset)[1:4,1:4]
all
ALL
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl")
allse=probe2gene(all.eset)
?probe2gene
library(hgu95av2.db)
BiocManager::install("regioneR")
library(EnrichmentBrowser)
library(regioneR)
allse=probe2gene(all.eset)
head(names(allse))
library(airway)
data("airway")
airse=airway[grep("^ENSG",names(airway))]
dim(airse)
assay(airse)[1:4,1:4]
allse$GROUP=ifelse(allse$mol.biol=="BCR/ABL",1,0)
table(allse$group)
airse$GROUP=ifelse(colData(airway)$dex=="trt",1,0)
head(airse$group)
table(airse$group)
airse$block=airway$cell
table(airse$block)
allse=deAna(allse)
?deAna
rowData(allse,use.names = T)
head(airse)
airse$GROUP

airse$GROUP
airse$GROUP
airse
head(assay(airse))
allse
assay$GROUP=c(0,1,0,1,0,1,0,1)
head(assay)
assay=assay(airse)
allse
airse
airse=deAna(airse)

yes
assay(airse)

airse$GROUP
airse1=deAna(airse)
?deAna
kegg.gs=getGenesets(org="hsa",db="kegg")
go.gs=getGenesets(org="hsa",db="go",go.onto = "BP",go.mode = "GO.db")
data.dir=system.file("extdata",package = "EnrichmentBrowser")
gmt.file=file.path(data.dir,"hsa_kegg_gs.gmt")
hsa.gs=getGenesets(gmt.file)
hsa.gs[1:2]
ora.all=sbea(method="ora",se=allse,gs=hsa.gs,perm=0,alpha=0.2)
gsRanking(ora.all)
eaBrowse(ora.all)
gsea.all=sbea(method = "gsea",se=allse,gs=hsa.gs,perm=1)
airse
colData(airse)
colData(allse)
?makeSummarizedExperimentFromExpressionSet

hsa.grn=compileGRN(org="hsa",db="kegg")
library(regioneR)

re=read.table("model-based-cpg-islands-hg19.txt",header = T)

head(re)
cpgHMM <- toGRanges(re)
cpgHMM=filterChromosomes(cpgHMM,chr.type = "canonical")
cpgHMM=sort(cpgHMM)
cpgHMM
bed=read.table("UCSC.promoters.hg19.bed",header = T)
head(bed)
promoters=toGRanges(bed)
promoters=filterChromosomes(promoters,chr.type = "canonical")
promoters=sort(promoters)
promoters
cpg=cpgHMM[seqnames(cpgHMM)%in%c("chr21","chr22")]
prom=promoters[seqnames(promoters)%in%c("chr21","chr22")]
pt=overlapPermTest(cpg,prom,genome="hg19",ntimes=100,per.chromosome=T,count.once=T)
pt
pt[[1]]
summary(pt[[1]]$permuted)
BiocManager::install("EnsDb.Hsapiens.v86")

library(MultiAssayExperiment)
library(GenomicRanges)
library(RaggedExperiment)
library(curatedTCGAData)
library(GenomicDataCommons)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(TCGAutils)
library(UpSetR)
library(mirbase.db)
library(AnnotationFilter)
library(EnsDb.Hsapiens.v86)
library(survival)
library(survminer)
library(pheatmap)
extends("SingleCellExperiment")
showClass("RaggedExperiment")
sample1=GRanges(c(A="chr1:1-10:-",B="chr1:8-14:+",C="chr1:15-18:+"),score=3:5,type=c("germline","somatic","germline"))
sample2=GRanges(c(D="chr1:1-10:-",E="chr1:11-18:+"),score=11:12,type=c("germline","somatic"))
sample1
sample2
coldata=DataFrame(id=1:2,status=factor(c("control","case")))
coldata
ragexp=RaggedExperiment(sample1=sample1,sample2=sample2,colData = coldata)
ragexp
grl=GRangesList(sample1=sample1,sample2=sample2)
ragexp2=RaggedExperiment(grl,colData = coldata)
identical(ragexp,ragexp2)
rowRanges(ragexp)
sparseAssay(ragexp)
unlist(grl)
assay(ragexp,"score")
assay(ragexp,"type")
compactAssay(ragexp)
unlist(grl)
compactAssay(ragexp,"type")
disjoinAssay(ragexp,simplifyDisjoin = mean)
grl
query=GRanges(c("chr1:1-14:-","chr1:15-18:+"))
query

weightedmean=function(scores,ranges,qranges){
  isects=pintersect(ranges,qranges)
  sum(scores*width(isects))/sum(width(isects))
}
qreduceAssay(ragexp,query,simplifyReduce = weightedmean)
data("miniACC")
miniACC
colData(miniACC)[1:4,1:4]
table(miniACC$race)
experiments(miniACC)
sampleMap(miniACC)
metadata(miniACC)
miniACC[c("MAPK14","IGFBP2"),]
miniACC[,miniACC$pathologic_stage=="stage iv",]
miniACC[,,"RNASeq2GeneNorm"]
miniACC[[1]]
summary(complete.cases(miniACC))
accmatched=intersectColumns(miniACC)
accmatched
colnames(accmatched)
accmatched=intersectRows(miniACC[,,c("RNASeq2GeneNorm","gistict","Mutations")])
rownames(accmatched)
class(assay(miniACC))
assays(miniACC)
longFormat(miniACC[c("TP53","CTNNB1"),,],colDataCols = c("vital_status","days_to_death"))
wideFormat(miniACC[c("TP53","CTNNB1"),,],colDataCols = c("vital_status","days_to_death"))

MultiAssayExperiment(experiments = experiments(miniACC),colData = colData(miniACC),sampleMap = sampleMap(miniACC),metadata = metadata(miniACC))
experiments(miniACC)
colData(miniACC)
sampleMap(miniACC)
metadata(miniACC)
library(curatedTCGAData)
curatedTCGAData("ACC")
acc=curatedTCGAData("ACC",assays=c("miRNASeqGene", "RPPAArray", "Mutation", "RNASeq2GeneNorm", "CNVSNP"),dry.run = F)
acc
?curatedTCGAData

simpa=TCGAutils::simplifyTCGA(miniACC)
upsetSamples(miniACC)

Surv(miniACC$days_to_death,miniACC$vital_status)
survivial=miniACC[,complete.cases(miniACC$days_to_death,miniACC$vital_status),]
survivial

fit=survfit(Surv(days_to_death,vital_status)~pathology_N_stage,data=colData(survivial))
ggsurvplot(fit,data=colData(survivial),risk.table = T)
wideacc=wideFormat(miniACC["EZH2",,],colDataCols = c("vital_status","days_to_death","pathology_N_stage"))
browseVignettes("curatedTCGAData")
wideacc
boxplot(wideacc$RNASeq2GeneNorm_EZH2)
wideacc$y=Surv(wideacc$days_to_death,wideacc$vital_status)
head(wideacc)
coxph(Surv(days_to_death,vital_status)~gistict_EZH2+log2(RNASeq2GeneNorm_EZH2)+pathology_N_stage,data=wideacc)

fit2=survfit(Surv(days_to_death,vital_status)~gistict_EZH2+log2(RNASeq2GeneNorm_EZH2)+pathology_N_stage,data=wideacc)
ggsurvplot(fit2, data = wideacc, risk.table = TRUE)


multi_cox=coxph(Surv(days_to_death,vital_status)~gistict_EZH2+log2(RNASeq2GeneNorm_EZH2)+pathology_N_stage,data=wideacc)
step_multi=step(multi_cox,direction = "both")
wideacc
multi_cox
?step
subacc=miniACC[,,c("RNASeq2GeneNorm", "gistict")]
subacc
subacc=intersectColumns(subacc)
subacc=intersectRows(subacc)
subacc
subacc.list=assays(subacc)
subacc.list
subacc.list[[1]]=log2(subacc.list[[1]]+1)
subacc.list=lapply(subacc.list, t)
corres=cor(subacc.list[[1]],subacc.list[[2]])
hist(diag(corres))
hist(corres[upper.tri(corres)])
?diag
View(corres)
which.max(diag(corres))
df=wideFormat(subacc["EIF4E",,])
head(df)
boxplot(RNASeq2GeneNorm_EIF4E~gistict_EIF4E,data=df,varwidth=T)
?prcomp
pc=prcomp(USArrests, scale = TRUE)
pc$rotation
rownames(miniACC)[[5]]
experiments(miniACC)
browseVignettes("ensembldb")
ebb=EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
genes(ebb)
