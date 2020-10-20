#########################################################################################
#                         -- KWMinima.R - 10.08.2017 --                                 #
# GIRDI:: .min dosyalari                                                                #
# ISLEM:: Kwee - van Woerden minima detection, 1956                                     #
# CIZIM:: 1) sT fonksiyonunu, polinom fiti ile birlikte cizer                           #
#         2) Minimum eðrisini çizer                                                     #
# CIKTI:: #YOK                                                                          #
# KAYNAK:: https://www.r-bloggers.com/fitting-polynomial-regression-in-r/               #
#          https://rpubs.com/kikihatzistavrou/80124                                     #
#          http://www.talkstats.com/threads/error-in-parabolic-regression-of-data-points.14476/
#          https://stackoverflow.com/questions/29999900/poly-in-lm-difference-between-raw-vs-orthogonal
#          https://courses.lumenlearning.com/boundless-algebra/chapter/graphs-of-quadratic-functions/
#                                                                                       #
#########################################################################################

### --BUYUK HARFLERLE-- YAZILAN BASLIKLAR ICINDEKI GIRDILER DUZENLENMELIDIR. Kucuk Harflerle yazilan basliklar altinda islemler yapilir. ###

####################
### Kutuphaneler ###
library(stats)
library(quantmod)

###################
### --AYARLAR-- ###
setwd("D:/Akademik/TEZ/DR/ExoLCA-TEZ/Kepler491 (KOI201)/O-C/KWMinima/deneme2")
# Girdi Sutunlari: ObsTime ObsFlux/Mag 
dosyalar <- list.files(path = ".", pattern = "*.min", all.files = FALSE)

Intp <- F # branch'larý ayri ayri interpole et
yuvarla <- 5 # JD kesirli yuvarlama hanesi
### --AYARLAR-- ###
###################

####################
### Fonksiyonlar ###
# Constructing Quadratic Formula
QuadRoot <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
        x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
        x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
        QuadRoot = c(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
        x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
      b^2-4*a*c
} 
### Fonksiyonlar ###
####################

####################
### --ISLEMLER-- ###  
for(i in 1:length(dosyalar)){
veri <- read.table(dosyalar[i], header=FALSE, skip=0, sep = "", stringsAsFactors=FALSE)

Tbas_tam <- as.integer(veri[1,1])
veri[,1] <- veri[,1] - as.integer(Tbas_tam)

N <- nrow(veri)
dT <- c()
for(n in 1:N-1){
  dT[n] <- veri[n+1,1] - veri[n,1]
}
dT <- mean(dT)
                         
Ti <- veri[1,1] 
Ts <- veri[N,1]
T1 <- Ti + (Ts - Ti)/2

##############
# min Kontrolu
#if(veri[which.min(veri[,2]),1] < T1){ T1 <- veri[which.min(veri[,2]),1] }
# min Kontrolu
##############

#sN <- (T1-Ti)/(dT/2)

Nn <- 10000

mag_b1 <- c(); mag_b2 <- c(); dmag <- c()
veri <- as.data.frame(approx(veri, method="linear", n=Nn))

############
# sT Dongusu
repeat{
branch1 <- veri[round(veri[,1],yuvarla) <= round(T1,yuvarla), ]
branch2 <- veri[round(veri[,1],yuvarla) >= round(T1,yuvarla), ]

if(Intp){ n = 10000 # AVE'nin tahmini degeri: 3640 
branch1 <- as.data.frame(approx(branch1, method="linear", n=n))
branch2 <- as.data.frame(approx(branch2, method="linear", n=n))
}

if(nrow(branch1) != nrow(branch2)){
  
  if(nrow(branch1) > nrow(branch2)){
    n <- nrow(branch2)
    b1 <- branch1[order(-branch1[,1]), ]
    b1 <- as.data.frame(lapply(b1, "[", c(1:n)))
    b2 <- branch2
  }
  
  if(nrow(branch2) > nrow(branch1)){
    n <- nrow(branch1)
    b2 <- as.data.frame(lapply(branch2, "[", c(1:n)))
    b2 <- b2[order(-b2[,1]), ]
    b1 <- branch1  
  }                                    
    
} else {

  n <- nrow(branch1)
  b1 <- branch1
  b2 <- branch2[order(-branch2[,1]), ]

} 

  for(j in 1:n){
    mag_b1[j] <- b1[j,2]
    mag_b2[j] <- b2[j,2]
    dmag[j] = mag_b1[j] - mag_b2[j]
  }
  
  sT1 <- sum(dmag^2)

  #Shifting dT for test
  # T1 + 0.5dT  
  T1art <- T1 + 0.5*dT
  branch1 <- veri[round(veri[,1],yuvarla) <= round(T1art,yuvarla), ]
  branch2 <- veri[round(veri[,1],yuvarla) >= round(T1art,yuvarla), ]

  if(Intp){ 
  branch1 <- as.data.frame(approx(branch1, method="linear", n=n))
  branch2 <- as.data.frame(approx(branch2, method="linear", n=n))
  }

  if(nrow(branch1) != nrow(branch2)){
    
    if(nrow(branch1) > nrow(branch2)){
      n <- nrow(branch2)
      b1 <- branch1[order(-branch1[,1]), ]
      b1 <- as.data.frame(lapply(b1, "[", c(1:n)))
      b2 <- branch2
    }
    
    if(nrow(branch2) > nrow(branch1)){
      n <- nrow(branch1)
      b2 <- as.data.frame(lapply(branch2, "[", c(1:n)))
      b2 <- b2[order(-b2[,1]), ]
      b1 <- branch1  
    }                                    
      
  } else {
  
    n <- nrow(branch1)
    b1 <- branch1
    b2 <- branch2[order(-branch2[,1]), ]
  
  }
  
  for(j in 1:n){
    mag_b1[j] <- b1[j,2]
    mag_b2[j] <- b2[j,2]
    dmag[j] = mag_b1[j] - mag_b2[j]
  }    
  sT1art <- sum(dmag^2)
  
  if(sT1art < sT1){
    T1 <- T1art
  } else {
  
    # T1 - 0.5dT           
    T1eks <- T1 - 0.5*dT
    branch1 <- veri[round(veri[,1],yuvarla) <= round(T1eks,yuvarla), ]
    branch2 <- veri[round(veri[,1],yuvarla) >= round(T1eks,yuvarla), ] 

  if(Intp){ 
  branch1 <- as.data.frame(approx(branch1, method="linear", n=n))
  branch2 <- as.data.frame(approx(branch2, method="linear", n=n))
  }

  if(nrow(branch1) != nrow(branch2)){
    
    if(nrow(branch1) > nrow(branch2)){
      n <- nrow(branch2)
      b1 <- branch1[order(-branch1[,1]), ]
      b1 <- as.data.frame(lapply(branch1, "[", c(1:n)))
      b2 <- branch2
    }
    
    if(nrow(branch2) > nrow(branch1)){
      n <- nrow(branch1)
      b2 <- as.data.frame(lapply(branch2, "[", c(1:n)))
      b2 <- b2[order(-b2[,1]), ]
      b1 <- branch1  
    }                                    
      
  } else {
  
    n <- nrow(branch1)
    b1 <- branch1
    b2 <- branch2[order(-branch2[,1]), ]
  
  }

    for(j in 1:n){
      mag_b1[j] <- b1[j,2]
      mag_b2[j] <- b2[j,2]
      dmag[j] = mag_b1[j] - mag_b2[j]
    }                                     
    sT1eks <- sum(dmag^2)
  
    if(sT1eks < sT1){
      T1 <- T1eks
    } else {
    break
    }
  }
}
# sT Dongusu  
############

Time <- c(T1eks, T1, T1art)
sT <- c(sT1eks, sT1, sT1art)

##############
# Polinom Fiti

# NOT
# R polinom fit uygulamasinda asagidaki iki yaklasim ayni olmakla birlikte 
# degisenlerin 'correlated' olmasini dikkate alir. 
# sTpoly <- lm(sT ~ poly(Time, 2, raw=F))
# sTpoly <- lm(sT ~ Time + I(Time^2))
# Note that when the explanatory variables are strongly correlated,
# the individual confidence intervals will usually underestimate
# the uncertainty in the parameter estimates, as indicated by the confidence region
# The rules of error propagation are different for correlated and uncorrelated errors
# -----------------------------------------------------------------------------------

sTfun <- function(Time) { Time + I(Time^2) }
optimize(sTfun, interval=c(Ti,Ts), maximum=F)

# Polinom Fiti: uncorrelated variables (orthogonal polynomials) 
sTpoly <- lm(sT ~ poly(Time, 2))
sTa <- summary(sTpoly)$coefficients[3,1]
sTb <- summary(sTpoly)$coefficients[2,1]
sTc <- summary(sTpoly)$coefficients[1,1]
delta(sTa,sTb,sTc)
sTpolyroot <- QuadRoot(sTa,sTb,sTc)

yeniTime <- seq(Ti, Ts, 10^-yuvarla)
yenidf <- data.frame(Time=yeniTime)
yenisT <- predict(sTpoly, newdata=yenidf)

grafikverisi <- data.frame(yeniTime,yenisT)
minima <- grafikverisi[findValleys(grafikverisi[,2], 0),]

# sT fonk. grafigi
plot(yeniTime, yenisT, type='l', col="grey", cex=0.5)#, ylim=c(-0.02,0.02)) #, xlim=c(1.04,1.06))
points(Time, sT, col="blue", pch=16)
lines(Time, predict(sTpoly), col="red")

# Polinom Fiti: correlated variables
sTpoly2 <- lm(yenisT ~ yeniTime + I(yeniTime^2))
sTa2 <- summary(sTpoly2)$coefficients[3,1]
sTb2 <- summary(sTpoly2)$coefficients[2,1]
sTc2 <- summary(sTpoly2)$coefficients[1,1]
#approx(x=Time, y=sT, xout=T0)$y
cf = coef(sTpoly2)

#y <- sTa2^2*T0 + sTb2*T0 + sTc2

delta(sTa2,sTb2,sTc2)
sTpoly2root <- QuadRoot(sTa2,sTb2,sTc2)
sTpoly2root

# Minimum time
T0 <- (-sTb2/(2*sTa2))

# Minima mean error over Z
Z <- 1/2*N

# Minima mean error: poly(Time, 2)
sigT0 = sqrt((4*sTa*sTc-(sTb^2))/(4*(sTa^2)*(Z-1)))
#sigT0 = sqrt(abs((4*sTa2*sTc2-(sTb2^2)))/(4*(sTa2^2)*(Z-1)))

#Get our covariance matrix
v <- vcov(sTpoly2)
b <- coef(sTpoly2)
#use delta method to calculate variance
xminvar <- (1/(2*b[3]))^2*v[2,2] + (b[2]/(2*b[3]^2))^2*v[3,3] - (b[2]/(2*b[3]^3))*v[2,3]

# Minima mean error: Time + I(Time^2)
#sTpolyErr <- lm(sT ~ Time + I(Time^2))
#sTa <- summary(sTpolyErr)$coefficients[3,1]
#sTb <- summary(sTpolyErr)$coefficients[2,1]
#sTc <- summary(sTpolyErr)$coefficients[1,1]
#sigT02 = sqrt((4*sTa*sTc-(sTb^2))/(4*(sTa^2)*(Z-1)))

# Polinom Fiti
##############

T0 #+Tbas_tam
minima
sigT0
n
Time
sT
sTpoly

yaz <- paste0(dosyalar[i],"\t", T0+Tbas_tam,"\t", sigT0)
write(yaz, file="#cikti.dat", append=TRUE)

} # DONGU SONU: i
                                    
### --ISLEMLER-- ###  
####################

##########################
### --GRAFIGI YAZDIR-- ###
# pdf(), png(), postscript(): saydamligi desteklemez: pdf'i eps'ye cevir.
yazdir <- 0 
if(yazdir){
pdf(paste0(dosyalar[i],".pdf"))#, width = 40, height = 30)
# En altta dev.off komutunu kapatmali; bu nedenle 'yazdir' degiseni ile kontrol ediliyor.
}
### --GRAFIGI YAZDIR-- ###
##########################

##########################
##### --GRAFIK CIZ-- #####
dev.new()
plot(veri) #, ylim=c(0.990, 1.002))
#lines(approx(veri[,1], veri[,2], method="linear", n=n))
#lines(aradeger)
lines(b1)
lines(b2)
abline(v=T0, col="red", lwd="2")

#lines(spline(branch1[,1], branch1[,2]), df=70, method = "natural")

##### --GRAFIK CIZ-- #####
##########################

if(yazdir){
dev.off() }

####################
### --KOD SONU-- ###
winDialog("ok", paste0("Kod Sonlandý!","\n\nBÝLGÝ: ",
  "\nÝþlenen Dosya Top.Sayýsý = ", length(dosyalar))) 
  #"\nAtýlan Satýr Top.Sayýsý = ", GecersizSayaci))
### --KOD SONU-- ###
####################


#########################################################################################
#                          KWMinima.R - 10.08.2017                                      #
#                         ** Geliþtirme Notlari **                                      #
#                                                                                       #
#  s#: ...                                                                              #
#########################################################################################