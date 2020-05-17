library(maftools)
options(stringsAsFactors = F)
laml=read.maf("TCGA.KIRC.mutect.somatic.maf.gz")
project="TCGA_KIRC"
library(BSgenome)
library("deconstructSigs")

library(BSgenome.Hsapiens.UCSC.hg38)
mut=laml@data
mut=mut[mut$Variant_Type=="SNP",]
a=mut[,c(16,5,6,12,13)]
colnames(a)=c("Sample","chr", "pos","ref",  "alt")
View(a)
a$Sample=as.character(a$Sample)
plot(table(a$Sample),las=2)
sigs.input=mut.to.sigs.input(mut.ref = a,sample.id = "Sample",chr="chr",pos="pos",ref="ref",alt="alt",bsg=BSgenome.Hsapiens.UCSC.hg38)
View(sigs.input[,])
class(sigs.input)
barplot(as.numeric(sigs.input[1,]))

w=lapply(unique(a$Sample)[1:10],function(i){
  sample=whichSignatures(tumor.ref = sigs.input,signatures.ref = signatures.cosmic,sample.id = i,contexts.needed = T,tri.counts.method = "default")
  print(i)
  return(sample$weights)
})
w
View(w[[1]])
w=do.call(rbind,w)
library(pheatmap)
pheatmap(w,cluster_rows = T,cluster_cols =F)


View(w)

browseVignettes("deconstructSigs")
??`deconstructSigs-package`
