#########################################################################################
#                      -- MINoku.R - 20.09.2020 --                                      #
# GIRDI:: .MIN verileri                                                                 #
# ISLEM:: AVE minimum dosyalarýnýn verilerini okur, tek bir dosyaya yazdirir.           #
# CIZIM:: YOK                                                                           #
# CIKTI:: *.min:                                                                        #
#########################################################################################

### --BUYUK HARFLERLE-- YAZILAN BASLIKLAR ICINDEKI GIRDILER DUZENLENMELIDIR.      ###
### Kucuk Harflerle yazilan basliklar altinda islemler yapilir.                   ###

####################
### Kutuphaneler ###
library(lomb)   # lsp()
library(scales) # renk saydamligi -alpha() icin kullaniliyor
library(quantmod) # findPeaks()
library(grid)
library(lattice)

###################
### --AYARLAR-- ###
setwd("J:/Akademik/Araþtýrma Projeleri/2019 - BAP (Ortak)/Veri/KOI22 (K422b)/Extracted Data/AVE")
hedefAdi <- "Kepler422b"
t0 <- (2455010.25005 - 2454833)
P <- 7.8914483

dosyalar <- list.files(path = ".", pattern = "*.MIN", all.files = FALSE)
dosyasayisi <- length(dosyalar)
onJD <- 2454833

# DOSYA ISIMLERI
minimumlar <- paste0(hedefAdi,".min")
oc <- paste0(hedefAdi,".oc")
if(file.exists(minimumlar)){file.remove(minimumlar)}
if(file.exists(oc)){file.remove(oc)}

#veri_tum <- c()

for(i in 1:dosyasayisi){
veri <- read.table(dosyalar[i], header=FALSE, sep=" ", skip=0)

write.table(veri, file=minimumlar, col.names=FALSE, row.names=FALSE, quote=FALSE, append = TRUE, sep = "\t")

veri$V3 <- veri$V2
veri$V2 <- (t0 + round((veri$V1 - t0)/P, digits = 0)*P)-veri$V1
veri$V1 <- round((veri$V1 - t0)/P, digits = 0)

#veri_tum$V1 <- append(veri_tum, veri$V1)
#veri_tum$V2 <- append(veri_tum, veri$V2)

write.table(veri, file=oc, col.names=FALSE, row.names=FALSE, quote=FALSE, append = TRUE, sep = "\t")
}

# Grafik Cizimi
#par(mfrow = c(dosyasayisi/6, 6), oma=c(0,0,0,0))
#plot(veri_tum$V1, veri_tum$V2, xlab="", ylab="", axes=F, col="black", main="")



