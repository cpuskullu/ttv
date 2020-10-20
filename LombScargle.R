#########################################################################################
#                          -- LombScargle.R - 20.07.2015 --                             #
# GIRDI:: O-C verileri                                                                  #
# ISLEM:: Tek veya birden cok zaman serisi verisinin Lomb-Scargle frekans analizini     #
#         yapar, secilen frekanslari (nf) kullanarak O-C dagilimina en iyi sinus        #
#         egrisini uydurur.                                                             #
# CIZIM:: evrelendirilmis O-C, cikti *.png: LSP                                         #
# CIKTI:: *.out: LSP verisi, SINUS FIT ve LS analiz sonucu                              #
# KAYNAK:: http://menugget.blogspot.com.tr/2013/01/                                     #
#                 lomb-scargle-periodogram-for-unevenly.html                            #
#########################################################################################

### --BUYUK HARFLERLE-- YAZILAN BASLIKLAR ICINDEKI GIRDILER DUZENLENMELIDIR.      
### Kucuk Harflerle yazilan basliklar altinda islemler yapilir.                   

####################
### Kutuphaneler ###
library(lomb)
library(quantmod) # findPeaks()

###################
### --AYARLAR-- ###                                                             
setwd("D:/Akademik/TEZ/DR/ExoLCA-TEZ/Kepler491 (KOI201)/LCV/Time Data")
#dosyalar <- list.files(path = ".", pattern = "*.DIF", all.files = FALSE, include.dirs = FALSE)
dosyalar <- c()
dosyalar[1] <- "test.ctbl"
hedefAdi <- "Kepler491bSC"
LSP <- T # LSP cozumu icin T; RANDLSP cozumu icin F
xekseni <- 'frequency'  # frequency / period 
nyq <- F # Tayf penceresi Nyquist frekansýna kadar açýlsýn mý?
nyq_turu <- 'mean' # Nyquist frekansýnda dt'nin hesabý: min / mean
# ÖNEMLÝ: Nyquist frekansýnda dt'nin en küçük ve ortalama deðeri, sýfýrlar dýþarýda býrakýlarak  hesaplanýr.
nf <- 50 # Evrelendirmede kullanýlacak frekansýn sýrasý, baskýn frekans: lsp$peak.at, sonraki frekanslar icin 'tepeler' degisenine bak. 
alfa <- 0.05 # Güvenirlik sýnýrý (FAP) varsayýlan: 0.05 (%95) - Tepenin (Pnmax) degeri P-value.
dil <- 'en'  # en / tr

# Frekans/Period Aralýðý
bas <- 0.0  # Varsayýlan: f=0, P=1
son <- 5.0 # Varsayýlan: f=0.1, P=1000
trials <- 1000 # RANDLSP acilirsa gecerlidir. Varsayýlan: R=1000 (Deneme: repeats=10^5 ve ofac=100 iken 10saate sadece 5800e kadar ulasiyor.)

### --AYARLAR-- ###
###################

###################
### --GORUNUM-- ###
# Metin
bicem <- "serif"  # "serif", "sans", "mono", "symbol" 
gorun <- 1 # 1:plain, 2:bold, 3: italic

# Grafik
noktaboyut <- 1
noktatur <- 16
noktarenk <- "#00000090"

# Eksen
xi <- -1000
xs <- 500
yi <- -0.003
ys <- 0.003
uznciz <- 0.6
kisciz <- 0.3

kenarust <- 1 
kenarsol <- 5
kenarsag <- 2
kenaralt <- 5
### --GORUNUM-- ###
###################


# DONGU: i - dosyalar

for(j in 1:length(dosyalar)){
veri <- read.table(dosyalar[j], header=FALSE, sep="\t")

# DOSYA ISIMLERI
ciktiAdi <- strsplit(basename(dosyalar[j]), "\\.")[[1]][1]
ozetCiktiDosyasi <- paste0(ciktiAdi,"-nf",nf,"-Ozet.out")
LSPCiktiDosyasi <- paste0(ciktiAdi,"-LSP.out")
if(file.exists(ozetCiktiDosyasi)){file.remove(ozetCiktiDosyasi)}
                        
# Veriyi zamana göre sýrala
veri <- veri[order(veri[1]),]

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

dev.new()
if(LSP){
# LSP
lspcozum <- lsp(veri$V2, times=veri$V1, type = xekseni, from = bas, to = son, ofac = 14, plot=F, alpha=alfa)
} else{
# RANDLSP: $sig.level:NULL olacaktir. Yardim metni ornegide ayni sonucu veriyor
set.seed(444)
rand.times <- sample(1:length(veri$V1),length(veri$V1))
lspcozum <- randlsp(repeats=trials, rand.times, times=veri$V1, type = xekseni, from = bas, to = son, ofac = 10, plot=TRUE, alpha=alfa)
}

grafikverisi <- data.frame(lspcozum$scanned,lspcozum$power)
write.table(grafikverisi, file=LSPCiktiDosyasi, sep="\t", quote=FALSE)

# Tepe noktalarý bulma
tepeler <- grafikverisi[findPeaks(grafikverisi[,2], 0),]
tepeler <- tepeler[with(tepeler, order(-tepeler[,2])),]
#points(tepeler, col="red") 

# DÝL AYARLARI
xbaslik1 <- c("Frequency (cycl/P)", "Frekans (çevrim/P)") 
xbaslik2 <- c("Period (cycl)", "Dönem (çevrim)")
xbaslik3 <- c("Phase", "Evre")
ybaslik1 <- c("Normalized Power", "Normalize Þiddet")
ybaslik2 <- c("O-C", "O-C")

if(dil == 'en'){i=1}else{i=2}
dev.off()
##########################
### --GRAFIGI YAZDIR-- ###
# pdf(), png(), postscript(): saydamligi desteklemez: pdf'i eps'ye cevir.
yazdir <- 1 
if(yazdir){
png(paste0(ciktiAdi,".png"))#, width = 40, height = 30)
# En altta dev.off komutunu kapatmali; bu nedenle 'yazdir' degiseni ile kontrol ediliyor.
}
### --GRAFIGI YAZDIR-- ###
##########################

# Grafiði Biçimlendirme 
if(xekseni == 'frequency'){xbaslik <- xbaslik1[i]} else {xbaslik <- xbaslik2[i]}
#dev.new()
par(bg=NA, las=1, ps=20, lwd=2, mar=c(kenaralt,kenarsol,kenarust,kenarsag))
plot.lsp(lspcozum, main= "", type="l", level=TRUE, 
	xlab=xbaslik, ylab=ybaslik1[i], font=gorun, family=bicem, 
  xaxs = "i", yaxs = "i",  #axes=F,
	tcl=0.2, xlim=c(bas, son) #, ylim=c(0, 8)
)
if(xekseni == 'period'){
ticks <- seq(0, 3, by=1)
label <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1, at=c(1, 10, 100, 1000), labels=label, tcl=kisciz)

box()
}
# Tepe noktasý cizgisi
# abline(h=lspcozum$peak, col="grey", lty = "dotted")

if(yazdir){
dev.off() }
### LombScarglePeriodagram ###
##############################


################################
### Periodagram Evre Grafigi ###
dev.new()
if(xekseni == 'frequency'){ donem <- 1/ifelse(nf==1,lspcozum$peak.at[2],tepeler[nf-1,1]) } else { donem <- ifelse(nf==1,lspcozum$peak.at[2],tepeler[nf-1,1]) }
plot((veri[,1]-veri[1,1])/donem-as.integer((veri[,1]-veri[1,1])/donem),veri[,2], 
  xlab=xbaslik3[i], ylab=ybaslik2[i])#, ylim=c(38940,39040))

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
A0 <- 0.00080
E0 <- 60

xdata <- veri$V1
ydata <- veri$V2

OC <- function(x,A0,E0,P0) {A0*sin(2*pi*(x-E0)/P0)}

library(minpack.lm)
# EnKK yontemiyle egriyi uydur
sinfit <- nlsLM(formula = ydata ~ OC(xdata,A0,E0,P0[k]), start=list(A0 = A0, E0 = E0))#, weights=veri$V3)  # weights must be positive, otherwise R will produce an error)

# ozet
sinegrisiOzet[[k]] <- summary(sinfit)
}

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

# BÝLGÝ ÇIKTISI
# print(tepeler)
lspOzet <- summary(lspcozum)
#print(lspcozum$alpha)
#print(lspcozum$sig.level)
qchisq(1-lspcozum$p.value, df=lspcozum$n-2)/(lspcozum$n-2)


# DOSYAYA YAZDIR

write(ciktiAdi, file=ozetCiktiDosyasi, append=TRUE, sep="\n")
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
write(paste("LSP sig.level:",lspcozum$sig.level), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write(paste("LSP alpha:",lspcozum$alpha), file=ozetCiktiDosyasi, append=TRUE, sep="\n")
write("", file=ozetCiktiDosyasi, append=TRUE, sep="\n")
capture.output(lspOzet, file=ozetCiktiDosyasi, append=TRUE)
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

} # DONGU SONU: j

### DOSYA SONU