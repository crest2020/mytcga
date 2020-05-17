library(GenomicRanges)

gr=GRanges(seqnames = c("chr1","chr2","chr2"),ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),strand = c("+","-","-"))
gr
gr[1:2,]
gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),
                          end=c(100,200,300)),
           names=c("id1","id3","id2"),
           scores=c(100,90,50)
)
gr
mcols(gr)=DataFrame(name2=c("pax6","meis1","zic4"),score2=c(1,2,3))
gr
mcols(gr)=cbind(mcols(gr),
                DataFrame(name2=c("pax6","meis1","zic4")) )
gr
elementMetadata(gr)
values(gr)
gr$name3 = c("A","C", "B")
gr
?system.file
filepath=system.file("extdata","cpi.hg19.chr21.bed",package="compGenomRData")
filepath
library(compGenomRData)
devtools::install_github("compgenomr/compGenomRData")
?file.path
filepa)th=file.path("extdata","cpgi.hg19.chr21.bed")
filepath
cpgi.df = read.table(filepath, header = FALSE,
                     stringsAsFactors=FALSE) 
head(cpgi.df)
table(grep("_",cpgi.df)
cpgi.df=cpgi.df[grep("_",cpgi.df[,1],invert = T),]
dim(cpgi.df)
grep("_",cpgi.df
)
table(grep("_",cpgi.df[,1]))
grep("_",cpgi.df[,1])     
seq.gr=GRanges(seqnames = cpgi.df[,1],ranges=IRanges(start=cpgi.df$V2,end=cpgi.df$V3))
seq.gr
ref.filepath=file.path("extdata","refseq.hg19.chr21.bed")
ref.filepath
ref.df = read.table(ref.filepath, header = FALSE,
                    stringsAsFactors=FALSE) 
head(ref.df)
ref.gr=GRanges(seqnames=ref.df[,1],
               ranges=IRanges(start=ref.df[,2],
                              end=ref.df[,3]),
               strand=ref.df[,6],name=ref.df[,4])
ref.gr
tss.gr=ref.gr
end(tss.gr[strand(tss.gr)=="+",])  =start(tss.gr[strand(tss.gr)=="+",])
start(tss.gr[strand(tss.gr)=="-",])=end(tss.gr[strand(tss.gr)=="-",])
tss.gr=tss.gr[!duplicated(tss.gr),]
tss.gr
BiocManager::install("Gviz")
library(genomation)

library(ggpubr)
library(ggbio)
library(Gviz)
library(rtracklayer)
import.bed(ref.filepath)
session <- browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
session
peak.filepath=file.path("extdata","wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz")
pk1.gr=readBroadPeak(peak.filepath)
pk1.gr
subsetByOverlaps(pk1.gr,seq.gr)
countOverlaps(pk1.gr,seq.gr)
findOverlaps(pk1.gr,seq.gr)
n.ind=nearest(seq.gr,tss.gr)
dists=distanceToNearest(seq.gr,tss.gr,select="arbitrary")
dists
boxplot(log10(mcols(dists)[,1]+1),las=2)

promoter.gr=tss.gr
start(promoter.gr)=start(promoter.gr)-1000
end(promoter.gr)=end(promoter.gr)+1000
promoter.gr=promoter.gr[seqnames(promoter.gr)=="chr21"]
promoter.gr
library(Rsamtools)
bamfilepath=file.path("extdata","wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam")
bamfilepath
param=ScanBamParam(which=promoter.gr)
counts=countBam(bamfilepath,param = param)
counts
head(counts)
library(GenomicAlignments)
alns=readGAlignments(bamfilepath,param = param)
alns
covs=coverage(alns)
covs
covs=coverage(bamfilepath,param=param)
covs
bwfile=file.path("extdata","wgEncodeHaibTfbsA549.chr21.bw")
bwfile
bw.gr=import(bwfile,which=promoter.gr,as="RleList")
bw.gr

library(rtracklayer)
bwfile

nrows=200
ncols=6
counts=matrix(runif(nrows*ncols,1,1e4),nrows)
head(counts)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=paste0("gene", 1:200))
rowRanges
coldata=DataFrame(timepoint=1:6,row.names = LETTERS[1:6])
coldata
library(GenomicRanges)
library(SummarizedExperiment)
se=SummarizedExperiment(assays=list(counts=counts),
                        rowRanges=rowRanges, colData=coldata)
head(counts)
?SummarizedExperiment
se
rowData(se)
colData(se)
fdata(se)
assays(se)
se[1:5,1:3]
se[,se$timepoint==1]
ranges(se)
library(Gviz)
cpgi.track=AnnotationTrack(seq.gr,name="CpG")
gene.track=BiomartGeneRegionTrack(genome="hg19",chromosome = "chr21",start=27698681,end=28083310,name="ENSEMBL")
gene.track
chipseqfile=file.path("extdata","wgEncodeHaibTfbsA549.chr21.bw")
?DataTrack
cov.track=DataTrack(chipseqfile,type="l",name="coverge")
tracklist=list(cpgi.track,cov.track)
plotTracks(tracklist,from=27698681,28083310,chromosome = "chr21")
BiocManager::install("UCSC")
library(genomation)
transcriptFile=file.path("extdata","refseq.hg19.chr20.bed")
feat=readTranscriptFeatures(transcriptFile,remove.unusual = T,up.flank = 500,down.flank = 500)
feat
unlist(feat)
names(feat)
prom=feat$promoters
h3k4file=file.path("extdata","H1.ESC.H3K4me3.chr20.bw")
h3k4file
sm=ScoreMatrix(h3k4file,prom,type="bigWig",strand.aware = T)
library(ggbio)
data(ideoCyto, package = "biovizBase")
p <- autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")
p
cpgfile=file.path("extdata","CpGi.hg19.table.txt")
cpgfile
cpgi.gr=readGeneric(cpgfile,chr = 1, start = 2, end = 3,header=TRUE, 
                    keep.all.metadata =TRUE,remove.unusual=TRUE)
p+layout_karyogram(cpgi.gr)
p+layout_karyogram(cpgi.gr,aes(x=start,y=obsExp),geom="point",ylim = c(2,50), color = "red",
                   size=0.1,rect.height=1)
p=ggplot()+layout_circle(ideoCyto$hg19,geom="ideo",fill="white",colour="white",cytoband=T,radius=39,trackWidth=2)
p
p=p+layout_circle(cpgi.gr,geom="point",grid=T,size=0.01,aes(y=obsExp),color="red",radius=42,trackWidth=10)
p
p=p+layout_circle(as(seqinfo(ideoCyto$hg19),"GRanges"),geom="text",aes(label=seqnames),vjust=0,radius=55,trackWidth=7,size=3)
p
BiocManager::install("MotifDb")
library(MotifDb)
motifs=query(query(MotifDb,"Hsapiens"),"CTCF")
motifs
ctcf_motif=motifs[[9]]
ctcf_motif
BiocManager::install("TFBSTools")
library(seqLogo)
seqLogo(ctcf_motif)
library(TFBSTools)


library(AnnoProbe)
BiocManager::install("TCGAmutations")
devtools::install_github("jmzeng1314/AnnoProbe")
?cel
library(affy)
library(TCGAbiolinks)
browseVignettes("TCGAbiolinks")
devtools::install_github(repo = "PoisonAlien/TCGAmutations")
