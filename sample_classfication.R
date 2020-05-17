mdsvals = cmdscale(dist(t(dat)))
mdsvals = as.data.frame(mdsvals)
mdsvals$Type=factor(colData(se)[,'cancer.type.ch1'])
mdsvals$Normal = factor(colData(se)[,'normal.ch1'])
head(mdsvals)
#>                   V1        V2 Type Normal
#> GSM2772660  8.531331 -18.57115   BC     no
#> GSM2772661  8.991591 -13.63764   BC     no
#> GSM2772662 10.788973 -13.48403   BC     no
#> GSM2772663  3.127105 -19.13529   BC     no
#> GSM2772664 13.056599 -13.88711   BC     no
#> GSM2772665  7.903717 -13.24731   BC     no
And do the plot.

library(ggplot2)
ggplot(mdsvals, aes(x=V1,y=V2,shape=Normal,color=Type)) + 
  geom_point( alpha=0.6) + theme(text=element_text(size = 18))