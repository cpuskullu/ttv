#########################################################################################
#                      -- LombScargle_wOC.R - 06.08.2017 --                             #
# GIRDI:: O-C verileri                                                                  #
# ISLEM:: Tek veya birden cok zaman serisi verisinin Lomb-Scargle frekans analizini     #
#         yapar, secilen frekanslari (nf) kullanarak O-C dagilimina en iyi sinus        #
#         egrisini uydurur.                                                             #
# CIZIM:: evrelendirilmis O-C, cikti *.png: LSP                                         #
# CIKTI:: *.out: LSP verisi, SINUS FIT ve LS analiz sonucu                              #
# KAYNAK:: http://menugget.blogspot.com.tr/2013/01/                                     #
#                 lomb-scargle-periodogram-for-unevenly.html                            #
#########################################################################################

### --BUYUK HARFLERLE-- YAZILAN BASLIKLAR ICINDEKI GIRDILER DUZENLENMELIDIR.      ###
### Kucuk Harflerle yazilan basliklar altinda islemler yapilir.                   ###

####################
### Kutuphaneler ###
library(lomb)   # lsp()
library(scales) # renk saydamligi -alpha() icin kullaniliyor
library(quantmod) # findPeaks()

###################
### --AYARLAR-- ###
setwd("D:/Akademik/TEZ/DR/ExoLCA-TEZ/Kepler491 (KOI201)/O-C/Solutions LC")
dosyaAdi <- "Kepler491bLC_O-C-[1].dat"
hedefAdi <- "Kepler491b"

cozumLSP <- T # LSP cozumu icin T; RANDLSP cozumu icin F
xekseni <- 'frequency'  # frequency / period 
nyq <- T # Tayf penceresi Nyquist frekansýna kadar açýlsýn mý?  F ise aralik belirle!
nyq_turu <- 'mean' # Nyquist frekansýnda dt'nin hesabý: min / mean
# ÖNEMLÝ: Nyquist frekansýnda dt'nin en küçük ve ortalama deðeri, sýfýrlar dýþarýda býrakýlarak  hesaplanýr.
nf <- 2 # Evrelendirmede kullanýlacak frekansýn sýrasý, baskýn frekans icin 1: lspcozum$peak.at ile belirlenir, sonraki frekanslar icin 'tepeler' degisenine bak: ilk 5 freq icin head(tepeler,5) 
alfa <- 0.05 # Güvenirlik sýnýrý (FAP) varsayýlan: 0.05 (%95) - Tepenin (Pnmax) degeri P-value.
dil <- 'en'  # en / tr

# Baslik: \1 Nesne adi t0 P \2 Sutun Basliklari
# Girdi Sutunlari: Epoch O-C(d) Error Reference
veri <- read.table(dosyaAdi, header=FALSE, sep="\t", skip=2)
# DOSYA ISIMLERI
ozetCiktiDosyasi <- paste0(hedefAdi,"-nf",nf,"-Ozet.out")
if(file.exists(ozetCiktiDosyasi)){file.remove(ozetCiktiDosyasi)}
LSPCiktiDosyasi <- paste0(hedefAdi,"-LSP.out")

# Frekans/Period Araligi
bas <- 0.0  # Varsayýlan: f=0, P=1
son <- 0.90 # Varsayýlan: f=0.1, P=1000
trials <- 1000 # RANDLSP acilirsa gecerlidir. Varsayýlan: R=1000 (Deneme: repeats=10^5 ve ofac=100 iken 10saate sadece 5800e kadar ulasiyor.)

ocmodel <- F
# Eksen hesabi ve O-C sinus modeli icin Efemeris
bilgiSatiri <- strsplit(readLines(dosyaAdi, n=1), split=" ")
options(digits=15)
# Bilgi satiri okunanacak
T0 <- as.double(strsplit(bilgiSatiri[[1]][2], split="=")[[1]][2])
P <-  as.double(strsplit(bilgiSatiri[[1]][4], split="=")[[1]][2])
# Bilgi satiri yoksa:
# T0 <- 2456082.496115
# P <- 1.306186483

### --AYARLAR-- ###
###################

# Veriyi zamana göre sýrala
veri <- veri[order(veri[1]),]
        
###################
### --GORUNUM-- ###
# Metin
bicem <- "serif"  # "serif", "sans", "mono", "symbol" 
gorun <- 1 # 1:plain, 2:bold, 3: italic
metinboyut <- 32

# Grafik
icerikboyut <- 1
kalinlik <- 3
noktaboyut <- 1.8
noktatur <- 16
noktarenk <- "#00000090"
modelkalinlik <- 2
modelrenk <- alpha("gray20",1.0)
hatacubukciz <- T
hatacubukrenk <- alpha("gray20",1.0)
hatacubukkalinlik <- 1

kategorisutunu <- veri$V4  # Kategorileri iceren sutun
simgeler <- c(9,0,1,21,6,23,12,22,13) # matches levels of category column (in alphabetical order)
# TrES3b icin: c(23,22,21,24,25,22,21,23,24)
renkler <- c("red","green","blue","black","red","orange","gray","blue","green") # yada rainbow(<sayi>)
 # TrES3b icin: c("red","green","blue","gray","orange","purple","black","cyan","darkgreen")
calismaetiketi <- "This study"
onJD <- 2450000

# Eksen
eksenJD <- F  # Ikincil yatay eksen JD ise T; Yil ise F
eksenaciklik <- 100
anaciz <- 0.8
araciz <- 0.4
yuvarla <- function(x) {round(x+5,-2)}
xi <- yuvarla(min(veri$V1)-eksenaciklik)
xs <- yuvarla(max(veri$V1)+eksenaciklik)
yi <- round(min(veri$V2),3)-0.001
ys <- round(max(veri$V2),3)+0.001
kesisme <- "i" # Eksenlerin kesisme durumu:  "r" (regular), "i" (internal)

# Dil
xbaslik1 <- c("Frequency (cycl/P)", "Frekans (çevrim/P)") 
xbaslik2 <- c("Period (cycl)", "Dönem (çevrim)")
xbaslik3 <- c("Phase", "Evre")
ybaslik1 <- c("Normalized Power", "Normalize Þiddet")
ybaslik2 <- c("O-C", "O-C")

kenaralt <- 6
kenarsol <- 9
kenarust <- 6
kenarsag <- 6

### --GORUNUM-- ###
###################

par(mar=c(kenaralt,6,3,4), mgp=c(4, 1.3, 0),
     cex=icerikboyut, ps=metinboyut, font=gorun, family=bicem, 
     pch=noktatur, bg=NA, las=1, lwd=kalinlik, xaxs=kesisme, yaxs=kesisme)

########################
### Nyquist frekansý ###
zamanfarki <- (veri[nrow(veri),1] - veri[nrow(veri)-1,1])
for(i in (nrow(veri)-1):2){
      zamanfarki <- c(zamanfarki, (veri[i,1] - veri[i-1,1])) }

if(nyq_turu == 'min'){
   if(xekseni == 'frequency'){zf <- 1/(2*min(zamanfarki[zamanfarki > 0]))} else{zf <- (2*min(zamanfarki[zamanfarki > 0]))}
} else if(nyq_turu == 'mean') {
   if(xekseni == 'frequency'){zf <- 1/(2*mean(zamanfarki[zamanfarki > 0]))} else{zf <- (2*mean(zamanfarki[zamanfarki > 0]))}
}
if(nyq){ son <- zf }
### Nyquist frekansý ###
########################

##############################
### LombScarglePeriodagram ###
if(cozumLSP){
# LSP:      
lspcozum <- lsp(veri$V2,times=veri$V1,type=xekseni,from=bas,to=son,ofac=50,plot=F,alpha=alfa)
} else {
# RANDLSP: $sig.level:NULL olacaktir. Yardim metni ornegide ayni sonucu veriyor
set.seed(444)
rand.times <- sample(1:length(veri$V1),length(veri$V1))
lspcozum <- randlsp(repeats=trials, rand.times, times=veri$V1, type=xekseni, from=bas, to=son, ofac=10, plot=TRUE, alpha=alfa)
}

lspverisi <- data.frame(lspcozum$scanned,lspcozum$power)
write.table(lspverisi, file=LSPCiktiDosyasi, sep="\t", quote=FALSE)

if(dil == 'en'){i=1}else{i=2}

# Grafiði Biçimlendirme 
if(xekseni == 'frequency'){xbaslik <- xbaslik1[i]} else {xbaslik <- xbaslik2[i]}

plot.lsp(lspcozum, main= "", type="l", level=TRUE, 
	xlab=xbaslik, ylab=ybaslik1[i], #axes=F, 
	tcl=0.8, lwd.ticks = par("lwd"), xlim=c(bas+0.0005, son)  #, ylim=c(0,8)
)

if(xekseni == 'period'){
ticks <- seq(0, 3, by=1)
label <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1, at=c(1, 10, 100, 1000), labels=label, tcl=araciz)
box()
}         
# Tepe noktalarý bulma
tepeler <- lspverisi[findPeaks(lspverisi[,2], 0),]
tepeler <- tepeler[with(tepeler, order(-tepeler[,2])),] 

if(FALSE){
# En buyuk tepe noktasi cizgisi
abline(h=lspcozum$peak, col="gray", lty = "dotted")
# Tepe noktalarini isaretle
points(tepeler, col="red") 
}

### LombScarglePeriodagram ###
##############################

################################
### Periodagram Evre Grafigi ###
dev.new()
if(xekseni == 'frequency'){ donem <- 1/ifelse(nf==1,lspcozum$peak.at[2],tepeler[nf-1,1]) } else { donem <- ifelse(nf==1,lspcozum$peak.at[2],tepeler[nf-1,1]) }
plot((veri[,1]-veri[1,1])/donem-as.integer((veri[,1]-veri[1,1])/donem),veri[,2], 
  xlab=xbaslik3[i], ylab=ybaslik2[i])

fit <- lm(veri[,2] ~ poly(((veri[,1]-veri[1,1])/donem-as.integer((veri[,1]-veri[1,1])/donem)), 3))
lines(sort((veri[,1]-veri[1,1])/donem-as.integer((veri[,1]-veri[1,1])/donem)), fitted(fit)[order((veri[,1]-veri[1,1])/donem-as.integer((veri[,1]-veri[1,1])/donem))], col="red", lwd=2)

# Olasýlýk Daðýlýmý (P-deðeri Histogramý)
if(FALSE){
library(RobPer)
dev.new()
betavalues <- betaCvMfit(lspcozum$power)
crit.val <- qbeta((0.05)^(bas:son),shape1=betavalues[1], shape2=betavalues[2])
hist(lspcozum$power, breaks=100, freq=FALSE, col=8, main ="")
betafun <- function(x) dbeta(x, shape1=betavalues[1], shape2=betavalues[2])
curve(betafun, add=TRUE, lwd=2)
abline(v=crit.val, lwd=1, lty=2, col="red")
}
### Periodagram Evre Grafigi ###
################################

########################
### --Sinus Egrisi-- ###
P0 <- c(); sinegrisiOzet <- c()
for(k in 1:(nf)){
# Baslangic degerleri; P0, A0 ve E0 icin
P0[k] <- ifelse(k==1,lspcozum$peak.at[2],1/tepeler[k-1,1])  # Eger nf, 1'den farkliysa sin egrisi sonraki harmonikler icin cozulur. 
A0 <- 0.00180
E0 <- 20

xdata <- veri$V1
ydata <- veri$V2

OC <- function(x,A0,E0,P0) {A0*sin(2*pi*(x-E0)/P0)}

library(minpack.lm)
# EnKK yontemiyle egriyi uydur
sinfit <- nlsLM(formula = ydata ~ OC(xdata,A0,E0,P0[1]), start=list(A0 = A0, E0 = E0))#, weights=veri$V3)  # weights must be positive, otherwise R will produce an error)

# ozet
sinegrisiOzet[[k]] <- summary(sinfit)
} # DONGU SONU: k

# Draw the fit on the plot by getting the prediction from the fit at 200 x-coordinates across the range of xdata
sinegrisi <- data.frame(xdata = seq(xi,xs,len=2000))
sinegrisi$ydata <- predict(sinfit,newdata=sinegrisi)

# the sum of squared residuals
sum(resid(sinfit)^2)
# chi-square of y
qchisq(1-summary(sinfit)$coefficients[1,4], df=summary(sinfit)$df[2])/(summary(sinfit)$df[2])
# the parameter confidence intervals
#confint(fit, level=0.95)

### --Sinus Egrisi-- ###
########################

#######################
### --O-C Grafigi-- ###

dev.new(width=1280, height=720)
par(mar=c(kenaralt,kenarsol,kenarust,kenarsag), mgp=c(0, 1.3, 0),
 cex=icerikboyut, ps=metinboyut, font=gorun, family=bicem, pch=noktatur, bg=NA, las=1, lwd=kalinlik, xaxs=kesisme, yaxs=kesisme)

plot(veri$V1, veri$V2, xlab="", ylab="", type="n", axes=F,     
     pch=simgeler[kategorisutunu], col="black", bg=renkler[kategorisutunu], 
     xlim=c(xi, xs), ylim=c(yi, ys))

axis(side=1, at = seq(xi, xs, by = 200), tcl=anaciz, lwd=par("lwd"), labels=TRUE)
axis(side=1, at = seq(xi, xs, by = 25), tcl=araciz, lwd=par("lwd"), labels=FALSE)
axis(side=2, at = seq(yi, ys, by = 0.002), tcl=anaciz, lwd=par("lwd"), labels=TRUE)
axis(side=2, at = seq(yi, ys, by = 0.0005), tcl=araciz, lwd=par("lwd"), labels=FALSE)

mtext("Epoch", side=1, line=4)      
mtext("O–C (day)", side=2, line=7, las=3)

box()

# Sifir cizgisi
#lines(c(xi,xs),c(0,0), col="black", lty="solid", lwd=1)

if(hatacubukciz){
# Hata Çubuklarý: noktalardan once
arrows(veri$V1, veri$V2-veri$V3,veri$V1, veri$V2+veri$V3, 
     code=3, angle=90, length=0.05, col=hatacubukrenk, lwd=hatacubukkalinlik)
     }

# Noktalar: hata cubuklarindan sonra
points(veri$V1, veri$V2, pch=simgeler[kategorisutunu], cex=noktaboyut, lwd=1,
     col="black", bg=alpha(renkler[kategorisutunu],0.9))

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
legend("bottomleft", inset=.025, etiket_listesi, bty = "n", box.lty=0,     
     cex=0.3, text.font=1, text.col="black",  bg='gray50',   #text.col=renkler
     y.intersp=1.8, x.intersp = 2.5,  
     pch=simgeler, pt.cex=1, pt.lwd=1, col="black", pt.bg=renkler 
)
par(op)

if(ocmodel){
# O-C modeli ve hata sinir cizgileri cizimi 
#model <- read.table("qatar1bTTVmodel.dat", header=FALSE, sep="\t")
#lines(model$V1, model$V2, col="#00000050", lwd="2.0")
lines(sinegrisi$xdata, sinegrisi$ydata, col=modelrenk, lwd=modelkalinlik, lty="solid")

sigma = summary(sinfit)$parameters[1,2] 
#lines(c(xi,xs),c(sigma,sigma), col="black", lty="dashed", lwd=1)
#lines(c(xi,xs),c(-sigma,-sigma), col="black", lty="dashed", lwd=1)
}

#if(F){       
par(new=TRUE)
# Ikincil yatay eksen JD yada Yil olarak cizilir
if(eksenJD){
plot(veri$V4, veri$V5, xlab="", ylab="", axes=F, type="n", 
     xlim=c((T0+xi*P)-onJD,(T0+xs*P)-onJD), ylim=c(yi*24*60, ys*24*60))

axis(side=3, at=seq(yuvarla((T0+xi*P)-onJD),yuvarla((T0+xs*P)-onJD),by=200), 
     tcl=anaciz, lwd=par("lwd"), labels=TRUE)
axis(side=3, at=seq(yuvarla((T0+xi*P)-onJD),yuvarla((T0+xs*P)-onJD),by=100), 
     tcl=araciz, lwd=par("lwd"), labels=FALSE)
mtext(paste0("BJD – ",onJD), side=3, line=3)

} else {                                            
plot(veri$V4, veri$V5, xlab="", ylab="", axes=F, type="n", 
     xlim=c(as.numeric(format(as.Date((T0+xi*P)-2415018.5, origin="1899-12-30"), "%Y")),as.numeric(format(as.Date((T0+xs*P)-2415018.5, origin="1899-12-30"), "%Y"))), ylim=c(yi*24*60, ys*24*60))

axis(side=3, at=seq(as.numeric(format(as.Date((T0+xi*P)-2415018.5, origin="1899-12-30"), "%Y")),as.numeric(format(as.Date((T0+xs*P)-2415018.5, origin="1899-12-30"), "%Y")),by=1), 
     tcl=anaciz, lwd=par("lwd"), labels=TRUE)
axis(side=3, at=seq(as.numeric(format(as.Date((T0+xi*P)-2415018.5, origin="1899-12-30"), "%Y")),as.numeric(format(as.Date((T0+xs*P)-2415018.5, origin="1899-12-30"), "%Y")),by=0.5), 
     tcl=araciz, lwd=par("lwd"), labels=FALSE)
mtext(paste0("Year"), side=3, line=4)
}
# Ikincil dikey eksen O-C (dk) olarak cizilir
axis(side=4, at=seq(round(yi*24*60,0),round(ys*24*60,0), by = 0.5), 
     tcl=araciz, lwd=par("lwd"), labels=FALSE)
axis(side=4, at=seq(round(yi*24*60,0),round(ys*24*60,0), by = 2), 
     tcl=anaciz, lwd=par("lwd"), labels=TRUE)      
mtext("O–C (min)", side=4, line=4, las=3) 
#}

### --O-C Grafigi-- ###
#######################

####################
### Bilgi Kutusu ###
print(tepeler)
lspOzet <- summary(lspcozum)
lspOzet
print(lspcozum$alpha)
print(lspcozum$sig.level)
qchisq(1-lspcozum$p.value, df=lspcozum$n-2)/(lspcozum$n-2)
### Bilgi Kutusu ###
####################

# DOSYAYA YAZDIR
write(dosyaAdi, file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("-- SINUS CURVE FIT TO O-C -----------------", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
for(k in 1:nf){
write(paste("Harmonic = ",k), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("A0 = ",A0,"\t","P = ",P0[k],"\t","E0 = ",E0), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
capture.output(sinegrisiOzet[k], file=ozetCiktiDosyasi, append=TRUE)
write("-- ---------------------- -----------------", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
}
write("-- LOMB-SCARGLE FREQ ANALYSIS TO O-C ------", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("Nyquist frequency:",zf), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("LSP peak at frequency:",lspcozum$peak.at[1]), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("LSP peak at period:",lspcozum$peak.at[2]), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("LSP peak power:",lspcozum$peak), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
capture.output(lspOzet, file=ozetCiktiDosyasi, append=TRUE)
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("LSP alpha:",lspcozum$alpha), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("LSP sig.level:",lspcozum$sig.level), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write.table(tepeler, file=ozetCiktiDosyasi, append=TRUE)
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("-- --------------------------------- ------", file=ozetCiktiDosyasi, append=TRUE)

# False Alarm Probability (FAP) EKSÝK: MC iterasyon
# KAYNAK: Searching for Periodic Signals in Time Series Data presentation, Uni.Oulu
# N <- -6.362 + 1.193*lspcozum$n + 0.00098*(lspcozum$n^2)
# fapN <- 1-(1-exp(-lspcozum$peak))^N  # Unknown period
# fapP <- exp(-lspcozum$peak)  # Known period
# fapN
# fapP

####################
### --KOD SONU-- ###
winDialog("ok", paste0("Kod Sonlandi!"))
### --KOD SONU-- ###
####################
### DOSYA SONU