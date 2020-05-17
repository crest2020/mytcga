options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/"))
options(download.file.method = 'libcurl')
options(url.method='libcurl')
BiocManager::install("multiMiR")

library(topGO)
library(miRNAtap)
library(miRNAtap.db)
library(org.Hs.eg.db)
mir="miR-10b"
mir
predictions=getPredictedTargets(mir,species = "hsa",method="geom",min_src = 2)
head(predictions)


mir="hsa-miR-18a-3p"
tmp1=getPredictedTargets(mir,species = "hsa",method = "geom",min_src = 2)
head(rownames(tmp1))

library(multiMiR)
example1=get_multimir(mirna = mir,summary = T)
tmp2=example1@data
View(tmp2)
setdiff(rownames(tmp1),tmp2$target_entrez)

mir_mu=getPredictedTargets(mirna = "mmu-miR-9-5p",species = "mmu",method = "geom",min_src = 2)

example1=get_multimir(mirna="mmu-miR-9-5p",org="mmu",summary = T)
tmp2=example1@data
table(tmp2$database)
intersect(rownames(mir_mu),tmp2$target_entrez)
