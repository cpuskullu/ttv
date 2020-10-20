#####################################################################################
#                              O-C.R - 05.12.2015                                   #
#####################################################################################
### -- AYARLAR -- ###                                                               
setwd("D:/Akademik/TEZ/DR/ExoTTV/Qatar1b")
veri <- read.table("qatar1bTTV.dat", header=FALSE, sep="\t")
model <- read.table("qatar1bTTVmodel.dat", header=FALSE, sep="\t")

###########################
### -- BICIMLENDIRME -- ###
# Metin
bicem <- "serif"  # "serif", "sans", "mono", "symbol" 
gorun <- 1 # 1:plain, 2:bold, 3: italic

# Grafik
noktaboyut <- 1
noktatur <- 16
noktarenk <- "#00000090"

simgeler <- c(21,22,23,24,25) # matches levels of category column (in alphabetical order)
renkler <- c("red","green","blue","darkgrey","orange") # yada rainbow(<sayi>)
kategorisutunu <- veri$V6  # Kategorileri iceren sutun
calismaetiketi <- "This work"

# Eksen
xi <- -1000
xs <- 500
yi <- -0.003
ys <- 0.003
uznciz <- 0.6
kisciz <- 0.3

######################
### -- ÝÞLEMLER -- ###

# O-C Grafigi
par(mgp=c(4,1,0), mar=c(6,6,3,6), font=gorun, family=bicem, bg=NA, las=1, lwd=2, ps=18)
plot(veri$V1, veri$V2, type="n", #axes=F,
     xlab="Epoch", ylab="O–C (day)", 
     pch=simgeler[kategorisutunu], col="black", bg=renkler[kategorisutunu], lwd=1,
     xlim=c(xi, xs), ylim=c(yi, ys), tcl=uznciz, xaxt='n',
     #at = seq(-1000, 1000, by = 100),
     cex=noktaboyut, cex.lab=1.0, cex.axis=1.0, cex.main=0.5, cex.sub=0.1)

# Hata Çubuklarý: noktalardan once
arrows(veri$V1, veri$V2-veri$V3, 
       veri$V1, veri$V2+veri$V3, code=3, angle=90, length=0.05, lwd=1, col="#00000099")

# Noktalar: hata cubuklarindan once
points(veri$V1, veri$V2, pch=simgeler[kategorisutunu], col="black", bg=renkler[kategorisutunu], lwd=1,
     cex=noktaboyut, cex.lab=1.0, cex.axis=1.0, cex.main=0.5, cex.sub=0.1)

etiket_listesi <- levels(kategorisutunu)
etiket_sayisi <- length(etiket_listesi)
# Calisma etiketini sirasinin en alta at
for(i in 1:etiket_sayisi){
     if(etiket_listesi[i] == calismaetiketi) {
        etiket_listesi[i] <- etiket_listesi[etiket_sayisi]
        etiket_listesi[etiket_sayisi] <- calismaetiketi
        
        temp <- simgeler[i]
        simgeler[i] <- simgeler[etiket_sayisi]
        simgeler[etiket_sayisi] <- temp
        
        temp <- renkler[i]
        renkler[i] <- renkler[etiket_sayisi]
        renkler[etiket_sayisi] <- temp
        }
}

# Gosterge
op <- par(family = "sans")
legend("bottomleft", inset=.04, etiket_listesi, bty = "n", box.lty=0,     
     cex=0.6, text.font=1, text.col=renkler,  #bg='lightyellow',
     pch=simgeler, col="black", pt.bg=renkler, pt.lwd=1
)
par(op)

if(T){
axis(side=1, at = seq(xi, xs, by = 100), labels=TRUE, font=gorun, family=bicem, tcl=uznciz)
axis(side=1, at = seq(xi, xs, by = 25), labels=FALSE, tcl=kisciz)
}

# MODEL
lines(model$V1, model$V2, col="#00000050", lwd="2.0")
sigma = 0.000782163 
lines(c(-1700,1200),c(sigma,sigma), col="#00000080", lty="dashed", lwd="1.4")
lines(c(-1700,1200),c(-sigma,-sigma), col="#00000080", lty="dashed", lwd="1.4")
       
par(new=T, las=1, font=gorun, family=bicem)
plot(veri$V4, veri$V5, axes=F, type="n", xlab="", ylab="",
     #font=gorun, family=bicem, pch=noktatur, col=noktarenk,
     xlim=c(-96.6111746, 2033.426016), ylim=c(-4.32, 4.32), tcl=uznciz,     
     cex=noktaboyut, cex.lab=1.0, cex.axis=1.0, cex.main=0.5, cex.sub=0.1)

axis(side=3, at = seq(-90, 2000, by = 100), labels=TRUE,
     font=gorun, family=bicem, tcl=uznciz)
mtext("BJD – 2455500", side=3, line=3, las=3)      
axis(side=4, at = seq(-4, 4, by = 1), labels=TRUE, 
     font=gorun, family=bicem, tcl=uznciz)
mtext("O–C (min)", side=4, line=3, las=3)      
 
### DOSYA SONU