options(stringsAsFactors = F)
##下载数据
# ./gdc-client download -m ~/Documents/Nutstore/github/TCGA-KIRC-miRNA-example/GDC/gdc_manifest.2018-08-05-clinical.txt -d clinical

# ./gdc-client download -m ~/Documents/Nutstore/github/TCGA-KIRC-miRNA-example/GDC/gdc_manifest.2018-08-05-LUAD-miRNA-seq.txt -d miRNAseq

library(XML)
library(methods)
dir1=file.path(getwd(),"tcga/clinical")
##读取临床信息
all_files=list.files(path=dir1,pattern = "*.xml$",recursive = T)
all_files

x=all_files[1]
result=xmlParse(file=file.path(dir1,x))
result
rootnode=xmlRoot(result)
rootnode
xmldataframe=xmlToDataFrame(rootnode[2])
xmldataframe
head(rootnode)
rootnode
browseVignettes("XML")
??XML
head(rootnode[2])
View(t(xmldataframe))

c1=lapply(all_files,function(x){
  result=xmlParse(file=file.path(dir1,x))
  rootnode=xmlRoot(result)
  xmldataframe=xmlToDataFrame(rootnode[2])
  return(t(xmldataframe))
})
head(c1)
View(c1)
c1_df=t(do.call(cbind,c1))
View(c1_df)

#读取miRNA表达量
dir2=file.path(getwd(),'tcga/miRNA')
dir2
file2=list.files(path =dir2,pattern="*.mirnas.quantification.txt$",recursive = T)
head(file2)
length(file2)
result=read.table(file=file.path(dir2,file2[1]),sep="\t",header = T)

mi2=lapply(file2,function(x){
  result=read.table(file=file.path(dir2,x),sep="\t",header = T)[,1:2]
  return(result)
})
View(mi)
identical(mi,mi2)
mi_df=t(do.call(cbind,mi))

dim(mi_df)
mi_df[1:4,1:3]

colnames(mi_df)=mi_df[1,]
mi_df=mi_df[seq(2,nrow(mi_df),by=2),]
mi_df[1:4,1:4]
View(as.numeric(mi_df))
class(mi_df)
mi_df_dataframe=as.data.frame(mi_df)
head(mi_df_dataframe[1,])

##
library(survival)
?survexp
ggforst