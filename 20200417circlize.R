library(PharmacoGx)
browseVignettes("survival")
BiocManager::install("TENxBrainData")
library("TENxBrainData")
tenx=TENxBrainData()

library(DelayedArray)
mat=matrix(rep(1:20,1:20),ncol=2)
mat
de_mat=DelayedArray(seed=mat)
de_mat
library(Matrix)
Mat=Matrix(mat)
de_mat=DelayedArray(seed = mat)
de_mat
df=as.data.frame(mat)
df
da_df=DelayedArray(seed=df)
da_df
library(tibble)
tbl=as_tibble(mat)
tbl
da_tbl=DelayedArray(seed=tbl)
da_tbl
da_rle=RleArray(data=Rle(mat),dim=dim(mat))
?RleArray
da_rle
is(da_rle,"DelayedArray")
showClass(getClass("RleMatrix",where="DelayedArray"))
library(SummarizedExperiment)
se=SummarizedExperiment(da_rle,colData = DataFrame(row.names = c("A","B")))
se
assay(se)
library(HDF5Array)
hdf5_file=file.path("500_Effectively_Using_the_DelayedArray_Framework", "hdf5_mat.h5")
hdf5_file
rhdf5::h5ls(hdf5_file)
showtree(da_rle)

library(survival)
fit1=survfit(Surv(futime,fustat)~resid.ds,data=ovarian)
data(package="survival")
fit1
data("ovarian")
ovarian
fit1
print(fit1,rmean=730)
?print
summary(fit1,times=(0:4)*182.5,scale=365)
plot(fit1,col=1:2,xscale=365.25,lwd=2,mark.time=T,xlab="Years since study entry",ylab="Survival")
legend(450,1.5,c("No residual disease","Redisual disease"),col=1:2,lwd=2,bty="n")
?legend
fit2=survfit(Surv(time,status)~sex+ph.ecog,data=lung)
fit2
fit2[1:3]
plot(fit2[1:3],lty=1:3,lwd=2,xscale=365)
crdata <- data.frame(time= c(1:8, 6:8),
                     endpoint=factor(c(1,1,2,0,1,1,3,0,2,3,0),
                                     labels=c("censor", "a", "b", "c")),
                     istate=rep("entry", 11),
                     id= LETTERS[1:11])
crdata
tfit=survfit(Surv(time,endpoint)~1,data=crdata,id=id,istate = istate)
dim(tfit)
summary(tfit)
plot(tfit,col=1:4,lty=1:4,lwd=2,ylab="Probability in state")
myeloid[1:5,]
sfit0=survfit(Surv(futime,death)~trt,myeloid)
plot(sfit0,xscale=365,xaxs='r',col=1:2,lwd=2)
library(survival)
strata(c(0,1,0,1,1,1,0))
?t.test
library(mvtnorm)
set.seed(42)
n=50
mu=c(2,1)
sigma=diag(c(1,1.5))
sigma
alpha=0.95
mu_diff=1.2
x=rmvnorm(n,mu,sigma)
x
t.test(x=x[,1],y=x[,2],mu=mu_diff,paired = T,alternative = "two.sided")


sq=function(x,y){
  result=x^2+y^2
  return(result)
}

sq
class(sq)
sq(2,3)
sqprint=function(x,y){
  result=x^2+y^2
  cat("there is the result:",result,"\n")
}
sqprint(2,4)
sqprint

cpg=function(bedrow){
  cpglen=bedrow[3]-bedrow[2]+1
  if(cpglen>1500){
    cat("there is large\n")
  }
  else if(cpglen<=1500&cpglen>700){
    cat("there is normal\n")
  }
  else{
    cat("there is short\n")
  }
}
class(cpg)
cpg
for(i in 1:10){print("hello world")}
mat=rnorm(30,c(5,6))
mat
class(mat)
mat=array(rnorm(30),c(5,6))
mat
mat=as.matrix(mat)
mat
sapply(mat, sum)
class(mat)
BiocManager::install("compGenomRData")
library(devtool)
install_github("compgenomr/compGenomRData")
paste0("patient",1:4)
dd=data.frame("TRXX4"=c(11,13,2,1),"OCT4"=c(10,13,4,3),"PAX6"=c(1,3,10,9))
dd
rownames(dd)=paste0("patient",1:4)
dd
dist(dd,method = "manhattan")
dist(dd,method="euclidean")
as.dist(1-cor(t(dd)))
cor(dd)
cor(t(dd))
t(dd)
d=dist(dd)
hc=hclust(d,method="complete")
plot(hc)
BiocManager::install("circlize")
library(circlize)
set.seed(999)
n=1000
df=data.frame(factors=sample(letters[1:8],n,replace = T),x=rnorm(n),y=runif(n))
head(df)
circos.par("track.height"=0.1)
circos.initialize(factors = df$factors,x=df$x)
circos.track(factors=df$factors,y=df$y,panel.fun=function(x,y){
  circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2]+uy(5,"mm"),CELL_META$sector.index)
  circos.axis(labels.cex = 0.6)
})
CELL_META$name
col=rep(c("#FF0000", "#00FF00"),4)
circos.trackPoints(df$factors,df$x,df$y,col=col,pch=16,cex=0.5)
circos.text(-1,0.5,"text",sector.index = "a",track.index = 1)
bgcol=rep(c("#EFEFEF","#CCCCCC"),4)
circos.trackHist(df$factors,df$x,bin.size = 0.2,bg.col = bgcol,col=NA)
CELL_META
names(CELL_META)
?uy
circos.track(factors=df$factors,x=df$x,y=df$y,panel.fun=function(x,y){
  ind=sample(length(x),10)
  x2=x[ind]
  y2=y[ind]
  od=order(x2)
  circos.lines(x2[od],y2[od])
})
circos.update(sector.index="d",track.index=2,bg.col="#FF8080",bg.border="black")
circos.points(x=-2:2,y=rep(0.5,5),col="white")
circos.text(CELL_META$xcenter,CELL_META$ycenter,"updated",col="white")
circos.track(ylim=c(0,1),panel.fun=function(x,y){
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  breaks=seq(xlim[1],xlim[2],by=0.1)
  n_breaks=length(breaks)
  circos.rect(breaks[-n_breaks],rep(ylim[1],n_breaks-1),breaks[-1],rep(ylim[2],n_breaks-1),col=rand_color(n_breaks),border=NA)
})
?circos.rect
breaks
xlim=CELL_META$xlim
ylim=CELL_META$ylim
breaks=seq(xlim[1],xlim[2],by=0.1)
head(breaks)
xlim
length(breaks)
c(breaks[-n_breaks],rep(ylim[1],n_breaks-1),breaks[-1],rep(ylim[2],n_breaks-1))
xlim(CELL_META)
class(CELL_META)
CELL_META

circos.link("a",0,"b",0,h=0.4)
circos.link("c",c(-0.5,0.5),"d",c(-0.5,0.5),col="red")
circos.link("e",0,"g",c(-1,1),col="green",border = "black",lwd=2,lty=2)
circos.clear()
fa = c("d", "f", "e", "c", "g", "b", "a")
fa
f1=factor(fa)
f1
circos.initialize(factors = f1,xlim=c(0,1))
f2=factor(fa,levels=fa)

f2
circos.initialize(factors = f2,xlim=c(0,1))
factors=c("a","a","a","b","b")
x=1:5
y=5:1
circos.track(factors=factors,x=x,y=y,panel.fun=function(x,y){circos.points(x,y)})
get.cell.meta.data()
factors=c("a","b")
circos.initialize(factors,xlim=c(0,1))
circos.track(ylim=c(0,1))
circlize(0.5,0.5,sector.index = "a",track.index = 1)
reverse.circlize(90,0.9,sector.index = "a",track.index = 1)

circos.clear()
BiocManager::install(c("yaml","EBImage"))
library(yaml)
data = yaml.load_file("https://raw.githubusercontent.com/Templarian/slack-emoji-pokemon/master/pokemon.yaml")
data
set.seed(123)
pokemon_list=data$emojis[sample(length(data$emojis),40)]
pokemon_name=sapply(pokemon_list, function(x) x$name)
pokemon_name
pokemon_src=sapply(pokemon_list, function(x) x$src)
pokemon_src
library(EBImage)
circos.clear()
circos.par("points.overflow.warning"=F)
circos.initialize(pokemon_name,xlim=c(0,1))
circos.track(ylim=c(0,1),panel.fun=function(x,y){
  pos=(circlize(CELL_META$xcenter,CELL_META$ycenter))
  image=EBImage::readImage(pokemon_src[CELL_META$sector.numeric.index])
  circos.text(CELL_META$xcenter,CELL_META$cell.ylim[1]-uy(2,"mm"),CELL_META$sector.index,facing="clockwise",niceFacing = T,adj=c(1,0.5),cex=0.6)
  rasterImage(image,xleft = pos[1, 1] - 0.05, ybottom = pos[1, 2] - 0.05,
              xright = pos[1, 1] + 0.05, ytop = pos[1, 2]+ 0.05)
},bg.border=1,track.height=0.5)


circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  image = EBImage::readImage(pokemon_src[CELL_META$sector.numeric.index])
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - uy(2, "mm"),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  rasterImage(image, 
              xleft = pos[1, 1] - 0.05, ybottom = pos[1, 2] - 0.05,
              xright = pos[1, 1] + 0.05, ytop = pos[1, 2]+ 0.05)
}, bg.border = 1, track.height = 0.15)

fa = letters[1:10]
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
circos.initialize(fa, xlim = cbind(rep(0, 10), runif(10, 0.5, 1.5)))
circos.track(ylim = c(0, 1), track.height = uh(5, "mm"),
             panel.fun = function(x, y) {
               circos.lines(c(0, 0 + ux(5, "mm")), c(0.5, 0.5), col = "blue")
             })
circos.track(ylim = c(0, 1), track.height = uh(1, "cm"),
             track.margin = c(0, uh(2, "mm")),
             panel.fun = function(x, y) {
               xcenter = get.cell.meta.data("xcenter")
               circos.lines(c(xcenter, xcenter), c(0, uy(1, "cm")), col = "red")
             })
circos.track(ylim = c(0, 1), track.height = uh(5, "mm"),
             track.margin = c(0, uh(5, "mm")),
             panel.fun = function(x, y) {
               line_length_on_x = ux(1*sqrt(2)/2, "mm")
               line_length_on_y = uy(1*sqrt(2)/2, "mm")
               circos.lines(c(0, line_length_on_x), c(0, line_length_on_y), col = "orange")})
 

factors = letters[1:3]              
circos.initialize(factors = factors, xlim = c(1, 2))
circos.info()
circos.track(ylim = c(0, 1))
circos.info(sector.index = "a", track.index = 1)
circos.clear()


BiocManager::install("ComplexHeatmap")
library(circlize)
col_fun=colorRamp2(c(-2,0,2),c("green","yellow","red"))
col_fun
circlize_plot=function(){
  set.seed(12345)
  fa=letters[1:10]
  circos.initialize(fa,xlim=c(0,1))
  circos.track(ylim=c(0,1),panel.fun=function(x,y){
    circos.points(runif(20),runif(20),cex=0.5,pch=16,col=2)
    circos.points(runif(20),runif(20),cex=0.5,pch=16,col=3)
  })
  circos.track(ylim=c(0,1),panel.fun=function(x,y){
    circos.lines(sort(runif(20)),runif(20),col=4)
    circos.lines(sort(runif(20)),runif(20),col=5)
  })
  for(i in 1:10){
    circos.link(sample(fa,1),sort(runif(10))[1:2],
    sample(fa,1),sort(runif(10))[1:2],col="orange")
  }
  
}
circlize_plot()
library(ComplexHeatmap)
lgd_points=Legend(at=c("label1","label2"),type="points",legend_gp = gpar(col=2:3),title_position = "topleft",title="Track1")
lgd_lines = Legend(at = c("label3", "label4"), type = "lines", 
                   legend_gp = gpar(col = 4:5, lwd = 2), title_position = "topleft", 
                   title = "Track2")
circos.clear()
?circos.track
lgd_lines = Legend(at = c("label3", "label4"), type = "lines", 
                   legend_gp = gpar(col = 4:5, lwd = 2), title_position = "topleft", 
                   title = "Track2")
lgd_links = Legend(at = c(-2, -1, 0, 1, 2), col_fun = col_fun, 
                   title_position = "topleft", title = "Links")
lgd_list_vertical = packLegend(lgd_points, lgd_lines, lgd_links)
lgd_list_vertical
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
set.seed(999)
bed = generateRandomBed()
dim(bed)
head(bed)
bed=generateRandomBed(nr=200,nc=4)
head(bed)
bed=generateRandomBed(nc=2,fun=function(k)sample(letters,k,replace = T))
head(bed)
circos.initializeWithIdeogram()
text(0,0,"default",cex=1)
circos.info()
circos.initializeWithIdeogram(species = "mm10")
cytoband.file = system.file(package = "circlize", "extdata", "cytoBand.txt")
head(cytoband.file)
df=read.table(cytoband.file, colClasses = c("character", "numeric",
                                            "numeric", "character", "character"), sep = "\t")
head(df)
circos.initializeWithIdeogram(df)
circos.initializeWithIdeogram(chromosome.index = paste0("chr",c(3,5,2,8)))
circos.initializeWithIdeogram(plotType = c("axis", "labels"))
circos.initializeWithIdeogram(plotType = NULL)
circos.clear()
circos.par("start.degree" = 90)
circos.initializeWithIdeogram()
circos.clear()
circos.par("gap.degree" = rep(c(2, 4), 12))
rep(c(2, 4), 12)
circos.initializeWithIdeogram()
circos.clear()
set.seed(123)
circos.initializeWithIdeogram(plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)
df = data.frame(
  name  = c("TP53",  "TP63",    "TP73"),
  start = c(7565097, 189349205, 3569084),
  end   = c(7590856, 189615068, 3652765))
df
circos.genomicInitialize(df)
tp_family = readRDS(system.file(package = "circlize", "extdata", "tp_family_df.rds"))
head(tp_family)
circos.genomicInitialize(tp_family)
circos.track(ylim = c(0, 1), 
             bg.col = c("#FF000040", "#00FF0040", "#0000FF40"), 
             bg.border = NA, track.height = 0.05)
n = max(tapply(tp_family$transcript, tp_family$gene, function(x) length(unique(x))))
head(n)
circos.genomicTrack(tp_family, ylim = c(0.5, n + 0.5), 
                    panel.fun = function(region, value, ...) {
                      all_tx = unique(value$transcript)
                      for(i in seq_along(all_tx)) {
                        l = value$transcript == all_tx[i]
                        # for each transcript
                        current_tx_start = min(region[l, 1])
                        current_tx_end = max(region[l, 2])
                        circos.lines(c(current_tx_start, current_tx_end), 
                                     c(n - i + 1, n - i + 1), col = "#CCCCCC")
                        circos.genomicRect(region[l, , drop = FALSE], ytop = n - i + 1 + 0.4, 
                                           ybottom = n - i + 1 - 0.4, col = "orange", border = NA)
                      }
                    }, bg.border = NA, track.height = 0.4)
circos.clear()
extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}
set.seed(123)
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)
circos.initializeWithIdeogram(plotType = c("labels", "axis"))
circos.track(ylim = c(0, 1))
circos.genomicIdeogram()
circos.genomicIdeogram(track.height = 0.2)
circos.clear()
circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 100, nc = 4)
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.genomicHeatmap(bed, col = col_fun, side = "inside", border = "white")
col_fun
circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicHeatmap(bed, col = col_fun, side = "outside",
                      line_col = as.numeric(factor(bed[[1]])))
circos.genomicIdeogram()
circos.clear()
circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 50, fun = function(k) sample(letters, k, replace = TRUE))
head(bed)
bed[1, 4] = "aaaaa"
circos.genomicLabels(bed, labels.column = 4, side = "inside")
circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicLabels(bed, labels.column = 4, side = "outside",
                     col = as.numeric(factor(bed[[1]])), line_col = as.numeric(factor(bed[[1]])))
circos.genomicIdeogram()
circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicIdeogram()
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top")
})
circos.track(ylim = c(0, 1), track.height = 0.1)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "bottom", direction = "inside")
})
circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
circos.genoimcRainfall(bed)
circos.genomicRainfall(bed)
circos.genomicRainfall(bed_list, col = c("red", "green"))

load(system.file(package = "circlize", "extdata", "DMR.RData"))
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
bed_list = list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(DMR_hypo, col = c("#0000FF80"), track.height = 0.1)
circos.clear()
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
load(system.file(package = "circlize", "extdata", "tagments_WGBS_DMR.RData"))
head(tagments, n = 4)
head(DMR1, n = 4)
head(correspondance, n = 4)


mat = matrix(1:9, 3)
mat
rownames(mat) = letters[1:3]
mat
df = data.frame(from = letters[1:3], to = LETTERS[1:3], value = 1:3);df
set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
mat
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
mat
df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df
library(circlize)
chordDiagram(mat)
circos.clear()
chordDiagram(df)
circos.clear()
chordDiagram(mat, order = c("S2", "S1", "S3", "E4", "E1", "E5", "E2", "E6", "E3"))
circos.clear()
circos.par(gap.after = c(rep(5, nrow(mat)-1), 15, rep(5, ncol(mat)-1), 15))
chordDiagram(mat)
circos.clear()
circos.par(gap.after = c(rep(5, length(unique(df[[1]]))-1), 15, 
                         rep(5, length(unique(df[[2]]))-1), 15))
chordDiagram(df)
circos.clear()
chordDiagram(mat, big.gap = 10)
circos.clear()
circos.par(start.degree = 90, clock.wise = FALSE)
chordDiagram(mat)
circos.clear()
grid.col = c(S1 = "red", S2 = "green", S3 = "blue",
             E1 = "grey", E2 = "grey", E3 = "grey", E4 = "grey", E5 = "grey", E6 = "grey")
chordDiagram(mat, grid.col = grid.col)
chordDiagram(t(mat), grid.col = grid.col)
chordDiagram(mat, grid.col = grid.col, transparency = 0)
col_mat = rand_color(length(mat), transparency = 0.5)
dim(col_mat) = dim(mat)
chordDiagram(mat, grid.col = grid.col, col = col_mat)
col = rand_color(nrow(df))
chordDiagram(df, grid.col = grid.col, col = col)
col_fun = colorRamp2(range(mat), c("#FFEEEE", "#FF0000"), transparency = 0.5)
chordDiagram(mat, grid.col = grid.col)
grid.col
chordDiagram(df, grid.col = grid.col, col = col_fun)
chordDiagram(mat, grid.col = grid.col, row.col = 1:3)
chordDiagram(mat, grid.col = grid.col, column.col = 1:6)
chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.lty = 2, link.border = "red")
lwd_mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
lwd_mat[mat > 12] = 2
border_mat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
border_mat[mat > 12] = "red"
chordDiagram(mat, grid.col = grid.col, link.lwd = lwd_mat, link.border = border_mat)
border_mat2 = matrix("black", nrow = 1, ncol = ncol(mat))
rownames(border_mat2) = rownames(mat)[2]
colnames(border_mat2) = colnames(mat)
border_mat2
chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.border = border_mat2)
lty_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(1, 2, 3))
lwd_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(2, 2, 2))
border_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(1, 1, 1))
chordDiagram(mat, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df,
             link.border = border_df)
chordDiagram(mat, grid.col = grid.col, row.col = c("#FF000080", "#00FF0010", "#0000FF10"))
mat
col_mat[mat < 12] = "#00000000"
chordDiagram(mat, grid.col = grid.col, col = col_mat) 
col_fun = function(x) ifelse(x < 12, "#00000000", "#FF000080");
chordDiagram(mat, grid.col = grid.col, col = col_fun)
col_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), 
                    c("#FF000080", "#00FF0080", "#0000FF80"))
chordDiagram(mat, grid.col = grid.col, col = col_df)

col = rand_color(nrow(df))
col[df[[3]] < 10] = "#00000000"
chordDiagram(df, grid.col = grid.col, col = col)
col = rand_color(nrow(df))
chordDiagram(df, grid.col = grid.col, link.visible = df[[3]] >= 10)
chordDiagram(mat, grid.col = grid.col, link.sort = TRUE, link.decreasing = TRUE)
title("link.sort = TRUE, link.decreasing = TRUE", cex = 0.8)
chordDiagram(mat, grid.col = grid.col, link.sort = TRUE, link.decreasing = FALSE)
chordDiagram(mat, grid.col = grid.col, transparency = 0)
chordDiagram(mat, grid.col = grid.col, transparency = 0, link.rank = rank(mat))
chordDiagram(df, grid.col = grid.col, transparency = 0, link.rank = rank(df[[3]]))
df2 = data.frame(start = c("a", "b", "c", "a"), end = c("a", "a", "b", "c"))
chordDiagram(df2, grid.col = 1:3, self.link = 1)
chordDiagram(df2, grid.col = 1:3, self.link = 2)
mat3 = matrix(rnorm(25), 5)
colnames(mat3) = letters[1:5]
cor_mat = cor(mat3)
col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
chordDiagram(cor_mat, grid.col = 1:5, symmetric = TRUE, col = col_fun)
chordDiagram(cor_mat, grid.col = 1:5, col = col_fun)
chordDiagram(mat, grid.col = grid.col, directional = 1)
chordDiagram(mat, grid.col = grid.col, directional = 1, diffHeight = uh(5, "mm"))
chordDiagram(mat, grid.col = grid.col, directional = -1)
mat2 = matrix(sample(100, 35), nrow = 5)
rownames(mat2) = letters[1:5]
colnames(mat2) = letters[1:7]
mat2
chordDiagram(mat2, grid.col = 1:7, directional = 1, row.col = 1:5)
mat3 = mat2
for(cn in intersect(rownames(mat3), colnames(mat3))) {
  mat3[cn, cn] = 0
}
mat3
chordDiagram(mat3, grid.col = 1:7, directional = 1, row.col = 1:5)
arr.col = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), 
                     c("black", "black", "black"))
chordDiagram(mat, grid.col = grid.col, directional = 1, direction.type = "arrows",
             link.arr.col = arr.col, link.arr.length = 0.2)
arr.col = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), 
                     c("black", "black", "black"))
chordDiagram(mat, grid.col = grid.col, directional = 1, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.col = arr.col, link.arr.length = 0.2)
chordDiagram(matx, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
chordDiagram(matx, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", diffHeight = -uh(2, "mm"))
chordDiagram(df, directional = 1)
set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
grid.col = c(S1 = "red", S2 = "green", S3 = "blue",
             E1 = "grey", E2 = "grey", E3 = "grey", E4 = "grey", E5 = "grey", E6 = "grey")
par(mfrow = c(1, 2))
chordDiagram(mat, grid.col = grid.col)
chordDiagram(mat, grid.col = grid.col, scale = TRUE)

mat = matrix(rnorm(36), 6, 6)
rownames(mat) = paste0("R", 1:6)
colnames(mat) = paste0("C", 1:6)
mat[2, ] = 1e-10
mat[, 3] = 1e-10

chordDiagram(mat)
circos.info()

df = expand.grid(letters[1:3], LETTERS[1:4])
df1 = df
df1$value = sample(10, nrow(df), replace = TRUE)
df2 = df
df2$value = -sample(10, nrow(df), replace = TRUE)
df = rbind(df1, df2)
df
grid.col = structure(1:7, names = c(letters[1:3], LETTERS[1:4]))
grid.col
chordDiagram(df, col = ifelse(df$value > 0, "red", "green"), grid.col = grid.col)
chordDiagram(df, col = "red", link.visible = df$value > 0, grid.col = grid.col)
chordDiagram(df, col = "green", link.visible = df$value < 0, grid.col = grid.col)
df = expand.grid(letters[1:3], LETTERS[1:4])
df$value = 1
df$value2 = 3 

chordDiagram(df[, 1:3], grid.col = grid.col)
chordDiagram(df[, 1:4], grid.col = grid.col) 

chordDiagram(mat)
circos.info()
chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid")
chordDiagram(mat, grid.col = grid.col, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.01))
chordDiagram(mat, grid.col = grid.col, annotationTrack = NULL)
chordDiagram(mat, preAllocateTracks = 2)
circos.info()
chordDiagram(mat, annotationTrack = NULL,
             preAllocateTracks = list(track.height = 0.3))
chordDiagram(mat, annotationTrack = NULL,
             preAllocateTracks = list(list(track.height = 0.1),
                                      list(bg.border = "black")))
chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) 
chordDiagram(mat, grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(5, "mm"))
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "white")
}

set.seed(123)
mat2 = matrix(rnorm(100), 10)
chordDiagram(mat2, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    xplot = get.cell.meta.data("xplot")
    
    circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3) # dotted line
    by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.2, 0.5)
    for(p in seq(by, 1, by = by)) {
      circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.1, 
                  paste0(p*100, "%"), cex = 0.3, adj = c(0.5, 0), niceFacing = TRUE)
    }
    
    circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, 0))
  }, bg.border = NA)
  
  
  chordDiagram(mat, grid.col = grid.col, big.gap = 20)
  abline(h = 0, lty = 2, col = "#00000080")  
  abline(v = 0, lty = 2, col = "#00000080")  

  
  mat1 = matrix(sample(20, 25, replace = TRUE), 5)  
  chordDiagram(mat1, directional = 1, grid.col = rep(1:5, 2), transparency = 0.5,
               big.gap = 10, small.gap = 1)  
  mat2 = mat1 / 2  
  gap = calc_gap(mat1, mat2, big.gap = 10, small.gap = 1)
  chordDiagram(mat2, directional = 1, grid.col = rep(1:5, 2), transparency = 0.5,
               big.gap = gap, small.gap = 1)  

  options(digits = 2)
  mat1 = matrix(rnorm(25), nrow = 5)
  rownames(mat1) = paste0("A", 1:5)
  colnames(mat1) = paste0("B", 1:5)
  
  mat2 = matrix(rnorm(25), nrow = 5)
  rownames(mat2) = paste0("A", 1:5)
  colnames(mat2) = paste0("C", 1:5)
  
  mat3 = matrix(rnorm(25), nrow = 5)
  rownames(mat3) = paste0("B", 1:5)
  colnames(mat3) = paste0("C", 1:5)
  
  
  mat = matrix(0, nrow = 10, ncol = 10)
  rownames(mat) = c(rownames(mat2), rownames(mat3))
  colnames(mat) = c(colnames(mat1), colnames(mat2))
  mat[rownames(mat1), colnames(mat1)] = mat1
  mat[rownames(mat2), colnames(mat2)] = mat2
  mat[rownames(mat3), colnames(mat3)] = mat3
  mat
  circos.clear()
  circos.par(gap.after = rep(c(rep(1, 4), 8), 3))
  chordDiagram(mat, annotationTrack = c("grid", "axis"),
               preAllocateTracks = list(
                 track.height = uh(4, "mm"),
                 track.margin = c(uh(4, "mm"), 0)
               ))  

  
  circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, niceFacing = TRUE)
  }, bg.border = NA)  
  
  
  highlight.sector(rownames(mat1), track.index = 1, col = "red", 
                   text = "A", cex = 0.8, text.col = "white", niceFacing = TRUE)
  highlight.sector(colnames(mat1), track.index = 1, col = "green", 
                   text = "B", cex = 0.8, text.col = "white", niceFacing = TRUE)
  highlight.sector(colnames(mat2), track.index = 1, col = "blue", 
                   text = "C", cex = 0.8, text.col = "white", niceFacing = TRUE)  
  circos.clear()  
  chordDiagram(mat, order = c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5)))  
  library(reshape2)  
mat1  
melt(mat1)
df2 = do.call("rbind", list(melt(mat1), melt(mat2), melt(mat3)))
chordDiagram(df2, order = c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5)))
circos.clear()
