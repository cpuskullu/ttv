setwd("D:/Akademik/TEZ/DR/ExoTTV/Qatar1b")
veri <- read.table("qatar1bTTV_vonEssen2013.dat", header=FALSE, sep="\t")

xekseni <- 'frequency'  # frequency / period
 
# Frekans/Period Aralýðý
bas <- 0.00  # Varsayýlan: frekans=0, period=1000
son <- 0.025

library(lomb)

veri <- veri[order(veri[1]),]

rlspcozum <- randlsp(repeats=50, veri, type = xekseni, from = bas, to = son, ofac = 100, plot=TRUE, alpha=0.05)


summary(rlspcozum)
qchisq(1-rlspcozum$p.value, df=rlspcozum$n-1)/(rlspcozum$n-1)