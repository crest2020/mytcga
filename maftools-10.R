library(RTCGA)
data(package="RTCGA")
browseVignettes("RTCGA")

BiocManager::install("TCGAmutations")
BiocManager::install("maftools")
devtools::install_github("PoisonAlien/TCGAmutations")
demo()
library(maftools)
laml=read.maf(maf="TCGA.KIRC.mutect.somatic.maf.gz")

laml@data=laml@data
View(laml@data)

(grep("^MUC",laml@data$Hugo_Symbol))
laml@data$Hugo_Symbol[1066]
laml@data$t_vaf=laml@data$t_alt_count/laml@data$t_depth

length(unique(laml@data$Tumor_Sample_Barcode))
View(getSampleSummary(laml))
View(getGeneSummary(laml))
View(getFields(laml))
plotmafSummary(maf=laml,rmOutlier = T,showBarcodes = T,addStat = "median",dashboard = T,titvRaw = F)

oncoplot(maf=laml,top=30,fontSize = 12,showTumorSampleBarcodes = F)

table(laml@data$Tumor_Sample_Barcode)
project="TCGA_KIRC"
oncoplot(maf=laml,top=15,fontSize = 12,clinicalFeatures = c("subtype"),sortByAnnotation = T)

dev.off()
laml.mutload=tcgaCompare(maf=laml,cohortName = project)


het=inferHeterogeneity(maf=laml,tsb="TCGA-CZ-5466-01A-01D-1501-10",vafCol = "t_vaf")
print(het$clusterMeans)
plotClusters(clusters = het)
mut=laml@data[,c("Hugo_Symbol","Chromosome","Start_Position","Tumor_Sample_Barcode","t_vaf")]
mut$pos=paste(mut$Chromosome,mut$Start_Position,sep=":")
tail(sort(table(mut$Hugo_Symbol)))
table(mut$Hugo_Symbol)
View(mut[order(mut$Hugo_Symbol),])
View(mut)
save.image(file="TCGA_KIRC_mut.Rdata")
