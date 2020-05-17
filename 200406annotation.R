library(GEOquery)
gset=getGEO("GSE2513",destdir = '.',AnnotGPL = F,getGPL = F)
?getGEO
a=gset[[1]]
dat=exprs(a)
boxplot(dat,las=2)
BiocManager::install("oligo")
library(oligo)
BiocManager::install("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
dir="GSE2513_RAW/"
od=getwd()
setwd(dir)
celfiles=list.celfiles(listGzipped = T)
celfiles
affyraw=read.celfiles(celfiles)
setwd(od)
eset=rma(affyraw)
a=eset
dat=exprs(a)
dim(dat)
dat[1:5,1:5]
boxplot(dat,las=2)
?boxplot
cols=rainbow(ncol(dat)*1.2)

boxplot(dat,col=cols,las=2)
library(limma)
dat=normalizeBetweenArrays(dat)
boxplot(dat,col=cols,las=2)
dd=dat
dd=t(dd)
dd=as.data.frame(dd)
dim(dd)
library(stringr)
dim(dd)
?str_match

data=data.table::fread()
rownames(dd)
group=c(rep("conjunctiva",6),rep("pterygium",6))
group_list=data.frame(rownames(dd),group)
group_list
group=as.factor(group)
group
dd=cbind(dd,group)
library(FactoMineR)
library(factoextra)
dim(dat)
dim(dd)
dd.pca=PCA(dd[,-ncol(dd)],graph = F)
fviz_pca_ind(dd.pca,geom.ind = "point",col.ind = group,palette=c('#00AFBB','#E7B800'),addEllipses = T)
colnames(group)="group"

library(limma)
design=model.matrix(~0+group)
design
contrast.matrix=makeContrasts('ptergium-conjunctiva',levels=design)
class(group)
group
design
colnames(design)=levels(group)
design
rownames(design)=colnames(dat)
design
contrast.matrix=makeContrasts('pterygium-conjunctiva',levels=design)
fit=lmFit(dat,design)
fit2=contrasts.fit(fit,contrast.matrix)
fit2=eBayes(fit2)
tem=topTable(fit2,coef=1,n=Inf)
nrDEG=na.omit(tem)

dim(nrDEG)
head(nrDEG)
?topTable
fit2
attach(nrDEG)
plot(logFC,-log10(P.Value))
library(ggpubr)
df=nrDEG
df$v=-log10(P.Value)
ggscatter(df,x="logFC",y='v',size=0.5,color = col)
col=ifelse(df$P.Value > 0.05,'black',
                 ifelse( df$logFC > 0.5849625,'red',
                         ifelse( df$logFC < -0.5849625,'blue','black') )
) 
length(df$P.Value)
length(col)

BiocManager::install("cgdsr")
nn
install.packages("DT")

library(cgdsr)
library(DT)
mycgds=CGDS("http://www.cbioportal.org/")
all_TCGA=getCancerStudies(mycgds)
DT::datatable(all_TCGA)
browseVignettes("cgdsr")
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
browseVignettes("edgeR")
limmaUsersGuide()
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,3,3,3)))
design
View(design)
library(edgeR)
edgeRUsersGuide()
group=factor(c(1,1,2,2))
group
model.matrix(~0+group)
BiocManager::install("DESeq2")
browseVignettes("DESeq2")
?read.csv
install.packages("Bioc2018Anno")
BiocManager::install("Bioc2018Anno")
BiocManager::install("Bioconductor/BiocWorkshops")
install.packages("devtools")
library(devtools)
install_github("jmacdon/Bioc2018Anno")
BiocManager::install("GenomicFeatures")
n

install.packages("D:/Rdocument/TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2.tar.gz", repos = NULL, type = "source")
BiocManager::install("EnsDb.Hsapiens.v79")
BiocManager::install("Homo.sapiens")
library(org.Hs.eg.db)
BiocManager::install("hgu95av2.db")

BSgenome.Hsapiens.UCSC.hg19
hugene20sttranscriptcluster.db
EnsDb.Mmusculus.v79
Organism.dplyr

library(hgu95av2.db)
ls("package:hgu95av2.db")
hgu95av2.db
columns(hgu95av2.db)
help("SYMBOL")
head(keys(hgu95av2.db))

k=head(keys(hgu95av2.db,keytype = "PROBEID"))   
select(hgu95av2.db,keys=k,columns = c("SYMBOL","GENENAME"),keytype = "PROBEID")
mapIds(hgu95av2.db,keys=k,column = c("GENENAME"),keytype = "PROBEID")
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
unikeys=head(keys(org.Hs.eg.db,keytype = "UNIPROT"))
unikeys
cols=c("SYMBOL","PATH")
cols
select(org.Hs.eg.db,keys = unikeys,keytype = "UNIPROT",columns = c("SYMBOL","PATH","UNIPROT"))
help("PATH")
load(system.file("extdata","resultTable.Rda",package = "AnnotationDbi"))
head(resultTable)
annots=select(org.Hs.eg.db,keys=rownames(resultTable),keytype = "ENTREZID",columns = c("SYMBOL","GENENAME","ENSEMBL"))
head(annots)

?grepl
name=keys(hgu95av2.db,keytype = "GENENAME")
index=(grepl(pattern = "lnc",name))
name[index]
?merge
resultTable1=merge(resultTable,annots,by.x=0,by.y="ENTREZID")
head(resultTable1)
library(GO.db)
GOIDS=c("GO:0042254","GO:0044183")
select(GO.db,keys=GOIDS,columns = "DEFINITION",keytype = "GOID")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
columns(txdb)
keytypes(txdb)
keys=head(keys(txdb,keytype = "GENEID"))
keys
cols=c("TXID","TXSTART","TXNAME")
cols
select(txdb,keys = keys,keytype = "GENEID",columns = cols)
library(EnsDb.Hsapiens.v79)
edb=EnsDb.Hsapiens.v79
edb
columns(edb)
keytypes(edb)
keys=head(keys(edb,keytype="GENEID"))

select(edb,keys=keys,columns = c("TXID","TXSEQSTART","TXBIOTYPE"),keytype = "GENEID")
filter(edb)
?filter
AnnotationFilterList()
linky=keys(edb,filter=list(GeneBiotypeFilter("lincRNA"),SeqNameFilter("Y")))
library(AnnotationFilter)
browseVignettes("AnnotationFilter")
supportedFilters()
length(linky)
txs=select(edb,keys=linky,keytype = "GENEID",columns = c("TXID"))
txs
select(edb,keys=list(GeneBiotypeFilter("lincRNA"),SeqNameFilter("Y")),columns = c("TXID"))
load(system.file("extdata/eset.Rdata",package = "Bioc2018Anno"))
eset
head(exprs(eset))
head(pData(phenoData(eset)))
head(pData(eset))
head(fData(eset))
library(hugene20sttranscriptcluster.db)
set.seed(12345)
ids=featureNames(eset)[sample(1:25000,5)]
ids
select(hugene20sttranscriptcluster.db,ids,"SYMBOL")
keytypes(hugene20sttranscriptcluster.db)
columns(hugene20sttranscriptcluster.db)
select(hugene20sttranscriptcluster.db,ids,c("SYMBOL","MAP"))
mapIds(hugene20sttranscriptcluster.db,ids,"SYMBOL","PROBEID",multiVals = "filter")
filter()
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
select(TxDb.Hsapiens.UCSC.hg19.knownGene,c("1","10"),c("TXNAME","TXCHROM"),"GENEID")
library(EnsDb.Hsapiens.v79)
select(EnsDb.Hsapiens.v79,c("1","10"),c("GENEID","GENENAME"),"ENTREZID")
gns=genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gns
txs=transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs
txs[txs%over%gns[1:2,]]
library(Homo.sapiens)
Homo.sapiens
head(genes(Homo.sapiens,columns=c("ENTREZID","ALIAS","UNIPROT")),4)
library(Organism.dplyr)
src=src_organism(dbpath = hg38light())
src
options(width=120)
promoters(src)
tbl(src,"id")
library(BSgenome)
head(available.genomes())
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
getSeq(Hsapiens,"chr1")
getSeq(Hsapiens,gns["5467",])
library(AnnotationHub)
hub=AnnotationHub()
library(biomaRt)
listMarts()
mart=useMart("ENSEMBL_MART_ENSEMBL")
mart
grep(pattern = "hsapien",listDatasets(mart)$dataset)
listDatasets(mart)[78,]
grep("1",c(1,2,3,1))
mart=useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl")
df <- data.frame(
  seqnames = rep(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  start = c(101, 105, 125, 132, 134, 152, 153, 160, 166, 170),
  end = c(104, 120, 133, 132, 155, 154, 159, 166, 171, 190),
  strand = rep(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10),
  row.names = head(letters, 10))
gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
gr
seq(6,0,length=10)
seqnames(gr)
ranges(gr)
granges(gr)
strand(gr)
seqnames(gr)
start(gr)
end(gr)
width(gr)
mcols(gr)
mcols(gr)$score
score(gr)
seqinfo(gr)=Seqinfo(genome="hg38")
seqinfo(gr)
names(gr)
length(gr)
gr[2:3,"GC"]
subset(gr,strand=="+"&score>5,select=score)
grmod=gr
grmod[2]=gr[1]
grmod
rep(gr[2],times=3)
gr
rev(gt)
rev(gr)
window(gr,start=2,end=4)
gr[IRanges(start=c(2,7),end=c(3,9))]
IRanges(start=c(2,7), end=c(3,9))
gr
seqinfo(gr)
sp=split(gr,rep(1:2,each=5))
sp
rep(1:2,each=5)
split(gr,~strand)
c(sp[[1]],sp[[2]])
stack(sp,index.var="group")
stack(sp)
aggregate(gr,score~strand,mean)
aggregate(gr,~strand,n_score=lengths(score),mean_score=mean(score))
g
flank(gr,10,start = F)
promoters(gr)
flank(unstrand(gr),10)
shift(gr,5)
gr
resize(gr,30)
reduce(gr)
reduce(gr,ignore.strand=T)
gaps(gr)
disjoin(gr)
cov=coverage(gr)
cov
gr
cov_gr=GRanges(cov)
cov_gr
cov_chr1=cov[[1]]
cov_chr1
Lengths(cov_chr1)
cov_chr1$Values
sum(runLength(cov_chr1)*runValue(cov_chr1))
runValue(cov_chr1)
cov=coverage(cov_gr,weight = "score")
cov
GPos(cov[1:3])
rg=reduce(gr,with.revmap=T)
rg
rg$score=mean(extractList(gr$score,rg$revmap))
rg
g2=head(gr,n=2)
g2
union(gr,g2)
intersect(gr,g2)
setdiff(gr,g2)
g3=gr[1:2]
g3
ranges(g3[1])=IRanges(start=105,end=112)
g3
punion(g2,g3)

set.seed(66+105+111+99+49+56)

pos <- sample(1:200, size = 30L)
size <- 10L
end <- size + pos - 1L
chrom <- sample(paste0("chr", 1:3), size = 30L, replace = TRUE)
query_df <- data.frame(chrom = chrom, 
                       start = pos,
                       end = end)
query_dfs <- split(query_df, 1:3)
q1 <- rename(query_dfs[[1L]], start = "pos")
q2 <- rename(query_dfs[[2L]], chrom = "ch", start = "st")
q3 <- rename(query_dfs[[3L]], end = "last")

q1
q2
q3
q1 <- makeGRangesFromDataFrame(q1, start.field = "pos")
q1
q2 <- makeGRangesFromDataFrame(q2, seqnames.field = "ch",
                               start.field = "st")
q3 <- makeGRangesFromDataFrame(q3, end.field = "last")
q3
query=mstack(q1,q2,q3,.index.var="replicate")
sort(query,by=~start)
subject=gr
subsetByOverlaps(query,subject,ignore.strand=T)
hits=findOverlaps(query,subject,ignore.strand=T)
hits
joined=query[queryHits(hits)]
joined
joined$score=subject$score[subjectHits(hits)]
joined
hitsbyquery=as(hits,"List")
hitsbyquery
counts=lengths(hitsbyquery)
counts
sum(countQueryHits(hits))
countOverlaps(query,subject,ignore.strand=T)
hitsbyquery
library(tools)
BiocManager::install("airway")
library(airway)
bams=list_files_with_exts(system.file("extdata",package="airway"),"bam")
bams
basename(bams)
names(bams)=sub("_[^_]+$","",basename(bams))
bams
View(bams)
library(Rsamtools)
bams=BamFileList(bams)
bams
first_bam=bams[[1]]
first_bam
first_bam_cvg=coverage(first_bam)
head(table(first_bam_cvg)[1,])
library(GenomicAlignments)
reads=grglist(readGAlignments(first_bam))
reads[lengths(reads)>2]
library(GenomicFeatures)
library(EnsDb.Hsapiens.v79)
tx=exonsBy(EnsDb.Hsapiens.v79,"gene")
tx
reads=keepStandardChromosomes(reads)
counts=countOverlaps(tx,reads,ignore.strand=T)
airway
library(airway)
data(package="airway")
airway

library(GEOquery)
gse=getGEO("GSE103512")[[1]]
BiocManager::install("GenomicDataCommons")
library("GenomicDataCommons")
GenomicDataCommons::status()
ge_manifest=files()%>%filter(~cases.project.project_id=="TCGA-OV"&type=="gene_expression"&analysis.workflow_type=="HTSeq - Counts")%>%manifest()
ge_manifest
library(tcltk)
BiocManager::install("tcltk")
library(PharmacoGx)
'magicaxis', 'downloader', 'lsa'