### Alpine plant species distribution models with abiotic, plant and microbe predictor variables
# By Cliff Bueno de Mesquita, Fall 2014, Spring 2015
# Paper published in Ecography, 2015

### Setup
library(AICcmodavg)
library(modEvA)
library(car)
library(MASS)
library(minpack.lm)
library(rgl)
library(robustbase)
library(Matrix)
library(qpcR)
library(boot)
source("~/Desktop/Functions/logisticPseudoR2s.R") # Code for pseudo R2 values
setwd("~/Desktop/CU/2Research/SDM")
data <- read.csv("GLVDATA.csv", header = TRUE)
# Remove 10 plots. Now all below 100 stems
plant.pa<-data[-c(1,19,31,43,49,51,52,56,60,61),]

# Make terms to test linearity
plant.pa$logALT<-log(plant.pa$ALT)*plant.pa$ALT
plant.pa$logSOIL_H2O<-log(plant.pa$SOIL_H2O)*plant.pa$SOIL_H2O
plant.pa$logMeanSnow<-log(plant.pa$MeanSnow)*plant.pa$MeanSnow
plant.pa$logSAND<-log(plant.pa$SAND)*plant.pa$SAND
plant.pa$logPH<-log(plant.pa$PH)*plant.pa$PH
plant.pa$logTDN<-log(plant.pa$TDN)*plant.pa$TDN
plant.pa$logDOC<-log(plant.pa$DOC)*plant.pa$DOC
plant.pa$logDptotal<-log(plant.pa$Dptotal)*plant.pa$Dptotal
plant.pa$logDpinorg<-log(plant.pa$Dpinorg)*plant.pa$Dpinorg
plant.pa$logcirsco<-log(plant.pa$cirscoA)*plant.pa$cirscoA
plant.pa$logsilaca<-log(plant.pa$silacaA)*plant.pa$silacaA
plant.pa$logtrillu<-log(plant.pa$trilluA)*plant.pa$trilluA
plant.pa$logphlsib<-log(plant.pa$phlsibA)*plant.pa$phlsibA
plant.pa$loganggra<-log(plant.pa$anggraA)*plant.pa$anggraA
plant.pa$logsenfre<-log(plant.pa$senfreA)*plant.pa$senfreA
plant.pa$loghymgra<-log(plant.pa$hymgraA)*plant.pa$hymgraA
plant.pa$logeriper<-log(plant.pa$eriperA)*plant.pa$eriperA
plant.pa$logoxydig<-log(plant.pa$oxydigA)*plant.pa$oxydigA
plant.pa$logbisbis<-log(plant.pa$bisbisA)*plant.pa$bisbisA
plant.pa$logantalp<-log(plant.pa$antalpA)*plant.pa$antalpA
plant.pa$loggenalg<-log(plant.pa$genalgA)*plant.pa$genalgA
plant.pa$logmoss<-log(plant.pa$mossA)*plant.pa$mossA
plant.pa$loggeuros<-log(plant.pa$geurosA)*plant.pa$geurosA
plant.pa$logelyscr<-log(plant.pa$elyscrA)*plant.pa$elyscrA
plant.pa$logtrispi<-log(plant.pa$trispiA)*plant.pa$trispiA
plant.pa$logfesrub<-log(plant.pa$fesrubA)*plant.pa$fesrubA
plant.pa$logdescae<-log(plant.pa$descaeA)*plant.pa$descaeA
plant.pa$logkobmyo<-log(plant.pa$kobmyoA)*plant.pa$kobmyoA
plant.pa$logcarnar<-log(plant.pa$carnarA)*plant.pa$carnarA
plant.pa$logcarper<-log(plant.pa$carperA)*plant.pa$carperA
plant.pa$logcarpha<-log(plant.pa$carphaA)*plant.pa$carphaA
plant.pa$logcarnig<-log(plant.pa$carnigA)*plant.pa$carnigA
plant.pa$logligfil<-log(plant.pa$ligfilA)*plant.pa$ligfilA

# Code to compare nested models 
#ModelChi <- model$deviance - model2$deviance
#chidf <- model$df.residual - model2$df.residual
#chisq.prob <- 1 - pchisq(ModelChi, chidf)
#chisq.prob

################################### Models ################################################
# Add different sets of variables manually. Forward and backward selection. Select based on AICc

# Carex nardina 
carnar<-glm(carnar ~ ALT + SOIL_H2O + MeanSnow + Dpinorg + DOC, family = binomial, data = plant.pa)
carnar
carnar2<-glm(carnar ~ ALT + SOIL_H2O + MeanSnow + Dpinorg + DOC + descaeA + hymgraA + kobmyoA + mossA, family = binomial, data = plant.pa)
carnar2
carnar3<-glm(carnar ~ ALT + SOIL_H2O + MeanSnow + Dpinorg + DOC + acidGP3 + rhodo, family = binomial, data = plant.pa)
carnar3
carnar4<-glm(carnar ~ ALT + SOIL_H2O + MeanSnow + Dpinorg + DOC + descaeA + hymgraA + kobmyoA + mossA + rhodo, family = binomial, data = plant.pa)
carnar4

AICc(carnar)
AICc(carnar2)
AICc(carnar3)
AICc(carnar4)

logisticPseudoR2s(carnar)
logisticPseudoR2s(carnar2)
logisticPseudoR2s(carnar3)
logisticPseudoR2s(carnar4)

Dsquared(model = carnar, adjust = TRUE)
Dsquared(model = carnar2, adjust = TRUE)
Dsquared(model = carnar3, adjust = TRUE)
Dsquared(model = carnar4, adjust = TRUE)

dwt(carnar4)
vif(carnar4)
1/vif(carnar4)
mean(vif(carnar4))

carnaraics<-c(80.4296, 69.1591, 78.4370, 69.5050)
akaike.weights(carnaraics)

carnarCV<-cv.glm(data=plant.pa,glmfit=carnar,K=10)
carnarAUC<-AUC(model=carnar)

carnar2CV<-cv.glm(data=plant.pa,glmfit=carnar2,K=10)
carnar2AUC<-AUC(model=carnar2)

carnar3CV<-cv.glm(data=plant.pa,glmfit=carnar3,K=10)
carnar3AUC<-AUC(model=carnar3)

carnar4CV<-cv.glm(data=plant.pa,glmfit=carnar4,K=10)
carnar4AUC<-AUC(model=carnar4)

carnarCV$delta
carnar2CV$delta
carnar3CV$delta
carnar4CV$delta

carnarAUC$AUC
carnar2AUC$AUC
carnar3AUC$AUC
carnar4AUC$AUC

# Carex phaeocephala
carpha<-glm(carpha ~ Dptotal, family = binomial, data = plant.pa)
carpha
carpha2<-glm(carpha ~ Dptotal + carperA + senfreA, family = binomial, data = plant.pa)
carpha2
carpha3<-glm(carpha ~ Dptotal + acidGP1, family = binomial, data = plant.pa)
carpha3
carpha4<-glm(carpha ~ Dptotal + carperA + senfreA + acidGP1, family = binomial, data = plant.pa)
carpha4

AICc(carpha)
AICc(carpha2)
AICc(carpha3)
AICc(carpha4)

logisticPseudoR2s(carpha)
logisticPseudoR2s(carpha2)
logisticPseudoR2s(carpha3)
logisticPseudoR2s(carpha4)

Dsquared(model = carpha, adjust = TRUE)
Dsquared(model = carpha2, adjust = TRUE)
Dsquared(model = carpha3, adjust = TRUE)
Dsquared(model = carpha4, adjust = TRUE)

dwt(carpha4)
vif(carpha4)
1/vif(carpha4)
mean(vif(carpha4))

carphaaics<-c(51.5697,49.2826,49.9701,48.5399)
akaike.weights(carphaaics)

carphaCV<-cv.glm(data=plant.pa,glmfit=carpha,K=10)
carphaAUC<-AUC(model=carpha)

carpha2CV<-cv.glm(data=plant.pa,glmfit=carpha2,K=10)
carpha2AUC<-AUC(model=carpha2)

carpha3CV<-cv.glm(data=plant.pa,glmfit=carpha3,K=10)
carpha3AUC<-AUC(model=carpha3)

carpha4CV<-cv.glm(data=plant.pa,glmfit=carpha4,K=10)
carpha4AUC<-AUC(model=carpha4)

carphaCV$delta
carpha2CV$delta
carpha3CV$delta
carpha4CV$delta

carphaAUC$AUC
carpha2AUC$AUC
carpha3AUC$AUC
carpha4AUC$AUC

# Deschampsia caespitosa
descae<-glm(descae ~ PH + SOIL_H2O + Dptotal, family = binomial, data = plant.pa)
descae
descae2<-glm(descae ~ PH + SOIL_H2O + carnarA + senfreA, family = binomial, data = plant.pa)
descae2
descae3<-glm(descae ~ PH + SOIL_H2O + Dptotal + acidGP1 + acidGP7, family = binomial, data = plant.pa)
descae3
descae4<-glm(descae ~ PH + SOIL_H2O + carnarA + senfreA + acidGP7 + acidGP1, family = binomial, data = plant.pa)
descae4

AICc(descae)
AICc(descae2)
AICc(descae3)
AICc(descae4)

logisticPseudoR2s(descae)
logisticPseudoR2s(descae2)
logisticPseudoR2s(descae3)
logisticPseudoR2s(descae4)

Dsquared(model = descae, adjust = TRUE)
Dsquared(model = descae2, adjust = TRUE)
Dsquared(model = descae3, adjust = TRUE)
Dsquared(model = descae4, adjust = TRUE)

dwt(descae4)
vif(descae4)
1/vif(descae4)
mean(vif(descae4))

descaeaics<-c(80.8759,74.6818,77.3984,69.0363)
akaike.weights(descaeaics)

descaeCV<-cv.glm(data=plant.pa,glmfit=descae,K=10)
descaeAUC<-AUC(model=descae)

descae2CV<-cv.glm(data=plant.pa,glmfit=descae2,K=10)
descae2AUC<-AUC(model=descae2)

descae3CV<-cv.glm(data=plant.pa,glmfit=descae3,K=10)
descae3AUC<-AUC(model=descae3)

descae4CV<-cv.glm(data=plant.pa,glmfit=descae4,K=10)
descae4AUC<-AUC(model=descae4)

descaeCV$delta
descae2CV$delta
descae3CV$delta
descae4CV$delta

descaeAUC$AUC
descae2AUC$AUC
descae3AUC$AUC
descae4AUC$AUC

# Elymus scriberneri
elyscr<-glm(elyscr ~ MeanSnow, family = binomial, data = plant.pa)
elyscr

elyscr2<-glm(elyscr ~ MeanSnow + anggraA, family = binomial, data = plant.pa)
elyscr2

elyscr3<-glm(elyscr ~ MeanSnow + delta, family = binomial, data = plant.pa)
elyscr3

elyscr4<-glm(elyscr ~ MeanSnow + anggraA + delta + acidGP1, family = binomial, data = plant.pa)
elyscr4

AICc(elyscr)
AICc(elyscr2)
AICc(elyscr3)
AICc(elyscr4)

logisticPseudoR2s(elyscr)
logisticPseudoR2s(elyscr2)
logisticPseudoR2s(elyscr3)
logisticPseudoR2s(elyscr4)

Dsquared(model = elyscr, adjust = TRUE)
Dsquared(model = elyscr2, adjust = TRUE)
Dsquared(model = elyscr3, adjust = TRUE)
Dsquared(model = elyscr4, adjust = TRUE)

dwt(elyscr4)
vif(elyscr4)
1/vif(elyscr4)
mean(vif(elyscr4))

elyscraics<-c(46.7872, 46.3103, 45.2885, 44.7362)
akaike.weights(elyscraics)

elyscrCV<-cv.glm(data=plant.pa,glmfit=elyscr,K=10)
elyscrAUC<-AUC(model=elyscr)

elyscr2CV<-cv.glm(data=plant.pa,glmfit=elyscr2,K=10)
elyscr2AUC<-AUC(model=elyscr2)

elyscr3CV<-cv.glm(data=plant.pa,glmfit=elyscr3,K=10)
elyscr3AUC<-AUC(model=elyscr3)

elyscr4CV<-cv.glm(data=plant.pa,glmfit=elyscr4,K=10)
elyscr4AUC<-AUC(model=elyscr4)

elyscrCV$delta
elyscr2CV$delta
elyscr3CV$delta
elyscr4CV$delta

elyscrAUC$AUC
elyscr2AUC$AUC
elyscr3AUC$AUC
elyscr4AUC$AUC

# Festuca Rubra
fesrub<-glm(fesrub ~ ALT + TDN, data = plant.pa, family = binomial)
fesrub

fesrub2<-glm(fesrub ~ ALT + TDN + trispiA + carphaA + geurosA + kobmyoA + elyscrA, family = binomial, data = plant.pa)
fesrub2

fesrub3<-glm(fesrub ~ ALT + TDN + oxalo, family = binomial, data = plant.pa)
fesrub3

fesrub4<-glm(fesrub ~ ALT + TDN + trispiA + carphaA + geurosA + kobmyoA + elyscrA + acidGP3, family = binomial, data = plant.pa)
fesrub4

AICc(fesrub)
AICc(fesrub2)
AICc(fesrub3)
AICc(fesrub4)

logisticPseudoR2s(fesrub)
logisticPseudoR2s(fesrub2)
logisticPseudoR2s(fesrub3)
logisticPseudoR2s(fesrub4)

Dsquared(model = fesrub, adjust = TRUE)
Dsquared(model = fesrub2, adjust = TRUE)
Dsquared(model = fesrub3, adjust = TRUE)
Dsquared(model = fesrub4, adjust = TRUE)

dwt(fesrub4)
vif(fesrub4)
1/vif(fesrub4)
mean(vif(fesrub4))

fesrubaics<-c(80.2030, 74.9179, 80.9339, 76.4520)
akaike.weights(fesrubaics)

fesrubCV<-cv.glm(data=plant.pa,glmfit=fesrub,K=10)
fesrubAUC<-AUC(model=fesrub)

fesrub2CV<-cv.glm(data=plant.pa,glmfit=fesrub2,K=10)
fesrub2AUC<-AUC(model=fesrub2)

fesrub3CV<-cv.glm(data=plant.pa,glmfit=fesrub3,K=10)
fesrub3AUC<-AUC(model=fesrub3)

fesrub4CV<-cv.glm(data=plant.pa,glmfit=fesrub4,K=10)
fesrub4AUC<-AUC(model=fesrub4)

fesrubCV$delta
fesrub2CV$delta
fesrub3CV$delta
fesrub4CV$delta

fesrubAUC$AUC
fesrub2AUC$AUC
fesrub3AUC$AUC
fesrub4AUC$AUC

# Kobresia myosuroides
kobmyo<-glm(kobmyo ~ ALT, family = binomial, data = plant.pa)
kobmyo

kobmyo2<-glm(kobmyo ~ ALT + fesrubA, family = binomial, data = plant.pa)
kobmyo2

kobmyo3<-glm(kobmyo ~ ALT + acidGP3, family = binomial, data = plant.pa)
kobmyo3

kobmyo4<-glm(kobmyo ~ ALT + fesrubA + acidGP3 + sphingo, family = binomial, data = plant.pa)
kobmyo4

AICc(kobmyo)
AICc(kobmyo2)
AICc(kobmyo3)
AICc(kobmyo4)

logisticPseudoR2s(kobmyo)
logisticPseudoR2s(kobmyo2)
logisticPseudoR2s(kobmyo3)
logisticPseudoR2s(kobmyo4)

Dsquared(model = kobmyo, adjust = TRUE)
Dsquared(model = kobmyo2, adjust = TRUE)
Dsquared(model = kobmyo3, adjust = TRUE)
Dsquared(model = kobmyo4, adjust = TRUE)

dwt(kobmyo4)
vif(kobmyo4)
1/vif(kobmyo4)
mean(vif(kobmyo4))

kobmyoaics<-c(71.6094, 71.7590, 69.3965, 69.5209)
akaike.weights(kobmyoaics)

kobmyoCV<-cv.glm(data=plant.pa,glmfit=kobmyo,K=10)
kobmyoAUC<-AUC(model=kobmyo)

kobmyo2CV<-cv.glm(data=plant.pa,glmfit=kobmyo2,K=10)
kobmyo2AUC<-AUC(model=kobmyo2)

kobmyo3CV<-cv.glm(data=plant.pa,glmfit=kobmyo3,K=10)
kobmyo3AUC<-AUC(model=kobmyo3)

kobmyo4CV<-cv.glm(data=plant.pa,glmfit=kobmyo4,K=10)
kobmyo4AUC<-AUC(model=kobmyo4)

kobmyoCV$delta
kobmyo2CV$delta
kobmyo3CV$delta
kobmyo4CV$delta

kobmyoAUC$AUC
kobmyo2AUC$AUC
kobmyo3AUC$AUC
kobmyo4AUC$AUC

# Trisetum spicatum
trispi<-glm(trispi ~ MeanSnow + SAND, family = binomial, data = plant.pa)
AICc(trispi)

trispi2<-glm(trispi ~ MeanSnow + SAND + antalpA + anggraA, family = binomial, data = plant.pa)
trispi2

trispi3<-glm(trispi ~ MeanSnow + SAND + pseudo + tm7, family = binomial, data = plant.pa)
trispi3

trispi4<-glm(trispi ~ MeanSnow + SAND + antalpA + anggraA + tm7 + pseudo, family = binomial, data = plant.pa)
trispi4

AICc(trispi)
AICc(trispi2)
AICc(trispi3)
AICc(trispi4)

logisticPseudoR2s(trispi)
logisticPseudoR2s(trispi2)
logisticPseudoR2s(trispi3)
logisticPseudoR2s(trispi4)

Dsquared(model = trispi, adjust = TRUE)
Dsquared(model = trispi2, adjust = TRUE)
Dsquared(model = trispi3, adjust = TRUE)
Dsquared(model = trispi4, adjust = TRUE)

dwt(trispi4)
vif(trispi4)
1/vif(trispi4)
mean(vif(trispi4))

trispiaics<-c(82.2695,79.7835,79.8116,78.3555)
akaike.weights(trispiaics)

trispiCV<-cv.glm(data=plant.pa,glmfit=trispi,K=10)
trispiAUC<-AUC(model=trispi)

trispi2CV<-cv.glm(data=plant.pa,glmfit=trispi2,K=10)
trispi2AUC<-AUC(model=trispi2)

trispi3CV<-cv.glm(data=plant.pa,glmfit=trispi3,K=10)
trispi3AUC<-AUC(model=trispi3)

trispi4CV<-cv.glm(data=plant.pa,glmfit=trispi4,K=10)
trispi4AUC<-AUC(model=trispi4)

trispiCV$delta
trispi2CV$delta
trispi3CV$delta
trispi4CV$delta

trispiAUC$AUC
trispi2AUC$AUC
trispi3AUC$AUC
trispi4AUC$AUC

# Cirsium scopulorum
cirsco<-glm(cirsco ~ MeanSnow + ALT + SOIL_H2O, family = binomial, data = plant.pa)
cirsco
cirsco2<-glm(cirsco ~ MeanSnow + ALT + SOIL_H2O + carphaA + phlsibA + silacaA, family = binomial, data = plant.pa)
cirsco2
cirsco3<-glm(cirsco ~ MeanSnow + ALT + SOIL_H2O + delta + tm7, family = binomial, data = plant.pa)
cirsco3
cirsco4<-glm(cirsco ~ MeanSnow + ALT + SOIL_H2O + carphaA + phlsibA + silacaA + delta, family = binomial, data = plant.pa)
cirsco4

AICc(cirsco)
AICc(cirsco2)
AICc(cirsco3)
AICc(cirsco4)

logisticPseudoR2s(cirsco)
logisticPseudoR2s(cirsco2)
logisticPseudoR2s(cirsco3)
logisticPseudoR2s(cirsco4)

Dsquared(model = cirsco, adjust = TRUE)
Dsquared(model = cirsco2, adjust = TRUE)
Dsquared(model = cirsco3, adjust = TRUE)
Dsquared(model = cirsco4, adjust = TRUE)

dwt(cirsco4)
vif(cirsco4)
1/vif(cirsco4)
mean(vif(cirsco4))

cirscoaics<-c(52.4776, 48.6046, 49.9694,43.3514)
akaike.weights(cirscoaics)

cirscoCV<-cv.glm(data=plant.pa,glmfit=cirsco,K=10)
cirscoAUC<-AUC(model=cirsco)

cirsco2CV<-cv.glm(data=plant.pa,glmfit=cirsco2,K=10)
cirsco2AUC<-AUC(model=cirsco2)

cirsco3CV<-cv.glm(data=plant.pa,glmfit=cirsco3,K=10)
cirsco3AUC<-AUC(model=cirsco3)

cirsco4CV<-cv.glm(data=plant.pa,glmfit=cirsco4,K=10)
cirsco4AUC<-AUC(model=cirsco4)

cirscoCV$delta
cirsco2CV$delta
cirsco3CV$delta
cirsco4CV$delta

cirscoAUC$AUC
cirsco2AUC$AUC
cirsco3AUC$AUC
cirsco4AUC$AUC

# Geum rossii
geuros<-glm(geuros ~ DOC + Dpinorg + Dptotal, family = binomial, data = plant.pa)
geuros

geuros2<-glm(geuros ~ DOC + Dpinorg + Dptotal + carphaA + hymgraA, family = binomial, data = plant.pa)
geuros2

geuros3<-glm(geuros ~ DOC + Dpinorg + Dptotal + cyano, family = binomial, data = plant.pa)
geuros3

geuros4<-glm(geuros ~ DOC + Dpinorg + Dptotal + carphaA + hymgraA + cyano + delta, family = binomial, data = plant.pa)
geuros4

AICc(geuros)
AICc(geuros2)
AICc(geuros3)
AICc(geuros4)

logisticPseudoR2s(geuros)
logisticPseudoR2s(geuros2)
logisticPseudoR2s(geuros3)
logisticPseudoR2s(geuros4)

Dsquared(model = geuros, adjust = TRUE)
Dsquared(model = geuros2, adjust = TRUE)
Dsquared(model = geuros3, adjust = TRUE)
Dsquared(model = geuros4, adjust = TRUE)

dwt(geuros4)
vif(geuros4)
1/vif(geuros4)
mean(vif(geuros4))

geurosaics<-c(58.3605, 51.1933, 58.2956, 50.1202)
akaike.weights(geurosaics)

geurosCV<-cv.glm(data=plant.pa,glmfit=geuros,K=10)
geurosAUC<-AUC(model=geuros)

geuros2CV<-cv.glm(data=plant.pa,glmfit=geuros2,K=10)
geuros2AUC<-AUC(model=geuros2)

geuros3CV<-cv.glm(data=plant.pa,glmfit=geuros3,K=10)
geuros3AUC<-AUC(model=geuros3)

geuros4CV<-cv.glm(data=plant.pa,glmfit=geuros4,K=10)
geuros4AUC<-AUC(model=geuros4)

geurosCV$delta
geuros2CV$delta
geuros3CV$delta
geuros4CV$delta

geurosAUC$AUC
geuros2AUC$AUC
geuros3AUC$AUC
geuros4AUC$AUC

# Oxyria digyna
oxydig<-glm(oxydig ~ MeanSnow + SOIL_H2O + PH + Dpinorg, family = binomial, data = plant.pa)
oxydig

oxydig2<-glm(oxydig ~ MeanSnow + SOIL_H2O + PH + Dpinorg + carphaA + fesrubA, family = binomial, data = plant.pa)
oxydig2

oxydig3<-glm(oxydig ~ MeanSnow + SOIL_H2O + PH + Dpinorg + acidGP7 + delta + pseudo, family = binomial, data = plant.pa)
oxydig3

oxydig4<-glm(oxydig ~ MeanSnow + SOIL_H2O + PH + Dpinorg + carphaA + fesrubA + pseudo + delta, family = binomial, data = plant.pa)
oxydig4

AICc(oxydig)
AICc(oxydig2)
AICc(oxydig3)
AICc(oxydig4)	

logisticPseudoR2s(oxydig)
logisticPseudoR2s(oxydig2)
logisticPseudoR2s(oxydig3)
logisticPseudoR2s(oxydig4)

Dsquared(model = oxydig, adjust = TRUE)
Dsquared(model = oxydig2, adjust = TRUE)
Dsquared(model = oxydig3, adjust = TRUE)
Dsquared(model = oxydig4, adjust = TRUE)

dwt(oxydig3)
vif(oxydig3)
1/vif(oxydig3)
mean(vif(oxydig3))

oxydigaics<-c(47.0146, 44.2358, 43.3568, 41.9309)
akaike.weights(oxydigaics)

oxydigCV<-cv.glm(data=plant.pa,glmfit=oxydig,K=10)
oxydigAUC<-AUC(model=oxydig)

oxydig2CV<-cv.glm(data=plant.pa,glmfit=oxydig2,K=10)
oxydig2AUC<-AUC(model=oxydig2)

oxydig3CV<-cv.glm(data=plant.pa,glmfit=oxydig3,K=10)
oxydig3AUC<-AUC(model=oxydig3)

oxydig4CV<-cv.glm(data=plant.pa,glmfit=oxydig4,K=10)
oxydig4AUC<-AUC(model=oxydig4)

oxydigCV$delta
oxydig2CV$delta
oxydig3CV$delta
oxydig4CV$delta

oxydigAUC$AUC
oxydig2AUC$AUC
oxydig3AUC$AUC
oxydig4AUC$AUC

# Senecio fremontii
senfre<-glm(senfre ~ PH, family = binomial, data = plant.pa)
senfre
senfre2<-glm(senfre ~ PH + mossA, family = binomial, data = plant.pa)
senfre2
senfre3<-glm(senfre ~ PH + oxalo + burk + ktedo, family = binomial, data = plant.pa)
senfre3
senfre4<-glm(senfre ~ PH + mossA + oxalo + burk + ktedo, family = binomial, data = plant.pa)
senfre4

AICc(senfre)
AICc(senfre2)
AICc(senfre3)
AICc(senfre4)

logisticPseudoR2s(senfre)
logisticPseudoR2s(senfre2)
logisticPseudoR2s(senfre3)
logisticPseudoR2s(senfre4)

Dsquared(model = senfre, adjust = TRUE)
Dsquared(model = senfre2, adjust = TRUE)
Dsquared(model = senfre3, adjust = TRUE)
Dsquared(model = senfre4, adjust = TRUE)

dwt(senfre4)
vif(senfre4)
1/vif(senfre4)
mean(vif(senfre4))

senfreaics<-c(75.9411, 74.4281, 66.1767, 68.0320)
akaike.weights(senfreaics)

senfreCV<-cv.glm(data=plant.pa,glmfit=senfre,K=10)
senfreAUC<-AUC(model=senfre)

senfre2CV<-cv.glm(data=plant.pa,glmfit=senfre2,K=10)
senfre2AUC<-AUC(model=senfre2)

senfre3CV<-cv.glm(data=plant.pa,glmfit=senfre3,K=10)
senfre3AUC<-AUC(model=senfre3)

senfre4CV<-cv.glm(data=plant.pa,glmfit=senfre4,K=10)
senfre4AUC<-AUC(model=senfre4)

senfreCV$delta
senfre2CV$delta
senfre3CV$delta
senfre4CV$delta

senfreAUC$AUC
senfre2AUC$AUC
senfre3AUC$AUC
senfre4AUC$AUC

# Silene acaulis
silaca<-glm(silaca ~ DOC + MeanSnow + ALT, family = binomial, data = plant.pa)
silaca

silaca2<-glm(silaca ~ DOC + MeanSnow + ALT + geurosA + anggraA + trispiA + fesrubA, family = binomial, data = plant.pa)
silaca2

silaca3<-glm(silaca ~ DOC + MeanSnow + ALT + acidGP1 + actinomyc, family = binomial, data = plant.pa)
silaca3

silaca4<-glm(silaca ~ DOC + MeanSnow + ALT + geurosA + anggraA + trispiA + fesrubA + burk, family = binomial, data = plant.pa)
silaca4

AICc(silaca)
AICc(silaca2)
AICc(silaca3)
AICc(silaca4)

logisticPseudoR2s(silaca)
logisticPseudoR2s(silaca2)
logisticPseudoR2s(silaca3)
logisticPseudoR2s(silaca4)

Dsquared(model = silaca, adjust = TRUE)
Dsquared(model = silaca2, adjust = TRUE)
Dsquared(model = silaca3, adjust = TRUE)
Dsquared(model = silaca4, adjust = TRUE)

dwt(silaca4)
vif(silaca4)
1/vif(silaca4)
mean(vif(silaca4))

silacaaics<-c(52.8054, 41.5426, 45.6006, 36.5904)
akaike.weights(silacaaics)

silacaCV<-cv.glm(data=plant.pa,glmfit=silaca,K=10)
silacaAUC<-AUC(model=silaca)

silaca2CV<-cv.glm(data=plant.pa,glmfit=silaca2,K=10)
silaca2AUC<-AUC(model=silaca2)

silaca3CV<-cv.glm(data=plant.pa,glmfit=silaca3,K=10)
silaca3AUC<-AUC(model=silaca3)

silaca4CV<-cv.glm(data=plant.pa,glmfit=silaca4,K=10)
silaca4AUC<-AUC(model=silaca4)

silacaCV$delta
silaca2CV$delta
silaca3CV$delta
silaca4CV$delta

silacaAUC$AUC
silaca2AUC$AUC
silaca3AUC$AUC
silaca4AUC$AUC

# Moss
moss<-glm(moss ~ MeanSnow + ALT, family = binomial, data = plant.pa)
moss

moss2<-glm(moss ~ MeanSnow + ALT + carperA + kobmyoA, family = binomial, data = plant.pa)
moss2

moss3<-glm(moss ~ MeanSnow + ALT + rhodo + oxalo, family = binomial, data = plant.pa)
moss3

moss4<-glm(moss ~ MeanSnow + ALT + carperA + kobmyoA + rhodo, family = binomial, data = plant.pa)
moss4

AICc(moss)
AICc(moss2)
AICc(moss3)
AICc(moss4)

logisticPseudoR2s(moss)
logisticPseudoR2s(moss2)
logisticPseudoR2s(moss3)
logisticPseudoR2s(moss4)

Dsquared(model = moss, adjust = TRUE)
Dsquared(model = moss2, adjust = TRUE)
Dsquared(model = moss3, adjust = TRUE)
Dsquared(model = moss4, adjust = TRUE)

dwt(moss4)
vif(moss4)
1/vif(moss4)
mean(vif(moss4))

mossaics<-c(85.0860, 82.4269, 82.6510, 81.5441)
akaike.weights(mossaics)

mossCV<-cv.glm(data=plant.pa,glmfit=moss,K=10)
mossAUC<-AUC(model=moss)

moss2CV<-cv.glm(data=plant.pa,glmfit=moss2,K=10)
moss2AUC<-AUC(model=moss2)

moss3CV<-cv.glm(data=plant.pa,glmfit=moss3,K=10)
moss3AUC<-AUC(model=moss3)

moss4CV<-cv.glm(data=plant.pa,glmfit=moss4,K=10)
moss4AUC<-AUC(model=moss4)

mossCV$delta
moss2CV$delta
moss3CV$delta
moss4CV$delta

mossAUC$AUC
moss2AUC$AUC
moss3AUC$AUC
moss4AUC$AUC

###################################### Predicted Probabilities ############################
carnardata<-subset(plant.pa,select=c(31,2,3,1,11,9,72))
carnardata$PP<-fitted(carnar3)
shapiro.test(carnardata$PP)

carphadata<-subset(plant.pa,select=c(33,10,81))
carphadata$PP<-fitted(carpha3)
shapiro.test(carphadata$PP)

descaedata<-subset(plant.pa,select=c(29,1,7,10,81))
descaedata$PP<-fitted(descae3)
shapiro.test(descaedata$PP)

elyscrdata<-subset(plant.pa,select=c(26,3,9,10,76))
elyscrdata$PP<-fitted(elyscr3)
shapiro.test(elyscrdata$PP)

fesrubdata<-subset(plant.pa,select=c(28,1,8,84))
fesrubdata$PP<-fitted(fesrub3)
shapiro.test(fesrubdata$PP)

kobmyodata<-subset(plant.pa,select=c(30,1,7,72))
kobmyodata$PP<-fitted(kobmyo3)
shapiro.test(kobmyodata$PP)

trispidata<-subset(plant.pa,select=c(27,3,4,75))
trispidata$PP<-fitted(trispi3)
shapiro.test(trispidata$PP)

cirscodata<-subset(plant.pa,select=c(12,3,1,2,76))
cirscodata$PP<-fitted(cirsco3)
shapiro.test(cirscodata$PP)

geurosdata<-subset(plant.pa,select=c(25,9,3,2,80))
geurosdata$PP<-fitted(geuros3)
shapiro.test(geurosdata$PP)

oxydigdata<-subset(plant.pa,select=c(20,7,2,3,11,9,74))
oxydigdata$PP<-fitted(oxydig3)
shapiro.test(oxydigdata$PP)

senfredata<-subset(plant.pa,select=c(17,7,84))
senfredata$PP<-fitted(senfre3)
shapiro.test(senfredata$PP)

silacadata<-subset(plant.pa,select=c(13,3,9,81))
silacadata$PP<-fitted(silaca3)
shapiro.test(silacadata$PP)

# Correlations between predicted prob. and bacterial abundances
cor.test(carnardata$acidGP3, carnardata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(carphadata$acidGP1, carphadata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(descaedata$acidGP1, descaedata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(elyscrdata$delta, elyscrdata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(fesrubdata$oxalo, trispidata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(kobmyodata$acidGP3, kobmyodata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(trispidata$pseudo, trispidata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(cirscodata$delta, cirscodata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(geurosdata$cyano, geurosdata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(oxydigdata$acidGP7, oxydigdata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(senfredata$oxalo, senfredata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)
cor.test(silacadata$acidGP1, silacadata$PP, alternative = "two.sided", method = "kendall", conf.level = 0.95)

#8/12 significant, 1 marginally significant

line1<-lm(PP ~ acidGP3, data = carnardata)
line2<-lm(PP ~ acidGP1, data = carphadata)
line3<-lm(PP ~ acidGP1, data = descaedata)
line4<-lm(PP ~ delta, data = elyscrdata)
line5<-lm(PP ~ oxalo, data = fesrubdata)
line6<-lm(PP ~ acidGP3, data = kobmyodata)
line12<-lm(PP ~ pseudo, data = trispidata)
line7<-lm(PP ~ delta, data = cirscodata)
line8<-lm(PP ~ cyano, data = geurosdata)
line9<-lm(PP ~ acidGP7, data = oxydigdata)
line10<-lm(PP ~ oxalo, data = senfredata)
line11<-lm(PP ~ acidGP1, data = silacadata)

# Figure 3
par(mfrow=c(4,3), oma=c(3,3,1,1) +0.1,mar=c(2,0.75,2,0.75) +0.1,mgp=c(2,1,0),xpd=NA)
par(xaxs="i")

plot(carnardata$acidGP3, carnardata$PP,ylab="",xlab="Acidobacteria Gp3", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="C. nardina", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=TRUE); lines(carnardata$acidGP3, fitted(line1), col="purple"); text(x=0.11, y=0.35, "τ = 0.38, p<0.01", cex=0.75)

plot(carphadata$acidGP1, carphadata$PP,ylab="",xlab="Acidobacteria Gp1", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="C. phaeocephala", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); text(x=0.11, y=0.35, "τ = 0.09, p>0.05", cex=0.75)

plot(descaedata$acidGP1, descaedata$PP,ylab="", xlab="Acidobacteria Gp1", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="D. cespitosa", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); text (x=0.11, y=0.35, "τ = 0.13, p>0.05", cex=0.75)

plot(elyscrdata$delta, elyscrdata$PP,ylab="",xlab="Deltaproteobacteria", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="E. scriberneri", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=TRUE); lines(elyscrdata$delta, fitted(line4), col="purple"); text(x=0.11, y=0.35, "τ = -0.50, p<0.01", cex=0.75)

plot(fesrubdata$oxalo, fesrubdata$PP,ylab="", xlab="Oxalobacteraceae", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="F. rubra", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); text (x=0.11, y=0.35, "τ = -0.12, p>0.05", cex=0.75)

plot(kobmyodata$acidGP3, kobmyodata$PP,ylab="", xlab="Acidobacteria Gp3", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="K. myosuroides", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); lines(kobmyodata$acidGP3, fitted(line6), col="purple"); text (x=0.11, y=0.35, "τ = 0.53, p<0.01", cex=0.75)

plot(trispidata$pseudo, trispidata$PP, ylab="", xlab="Pseudonocardiaceae", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="T. spicatum", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=TRUE); lines(trispidata$pseudo, fitted(line12), col="purple"); text(x=0.11, y=0.35, "τ = 0.58, p<0.01", cex=0.75)

plot(cirscodata$delta, cirscodata$PP,ylab="",xlab="Deltaproteobacteria", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="C. scopulorum", font.main=3, adj=0.99, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); text(x=0.11, y=0.75, "τ = 0.07, p>0.05", cex=0.75)

plot(geurosdata$cyano, geurosdata$PP,ylab="", xlab="Cyanobacteria", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="G. rossii", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); lines(geurosdata$cyano, fitted(line8), col="purple"); text (x=0.11, y=0.35, "τ = -0.38, p<0.01", cex=0.75)

plot(oxydigdata$acidGP7, oxydigdata$PP,ylab="", xlab="Acidobacteria Gp7", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="O. digyna", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=TRUE); lines(oxydigdata$acidGP7, fitted(line9), col="purple"); text (x=0.11, y=0.35, "τ = 0.28, p<0.01", cex=0.75)

plot(senfredata$oxalo, senfredata$PP,ylab="", xlab="Oxalobacteraceae", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="S. fremontii", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); lines(senfredata$oxalo, fitted(line10), col="purple"); text (x=0.11, y=0.35, "τ = 0.55, p<0.01", cex=0.75)

plot(silacadata$acidGP1, silacadata$PP,ylab="", xlab="Acidobacteria Gp1", cex.main=1, ylim=c(0,1),xlim=c(0,0.14),yaxt='n');title(main="S. acaulis", font.main=3, adj=0.85, cex.main=0.90, line=-0.69); axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), labels=FALSE); lines(silacadata$acidGP1, fitted(line11), col="purple"); text (x=0.11, y=0.35, "τ = -0.24, p<0.01", cex=0.75)

mtext("Predicted Probability of Occurrence",side=2,outer=TRUE,cex=0.8,line=1.3)
mtext("Bacteria Relative Abundance",side=1,outer=TRUE,cex=0.8,line=1.25)

################################# Deviance Explained ######################################
# D2 Values for Abiotic, Abiotic + Plant, and Full models
carnar<-c(0.171,0.401,0.425)
carpha<-c(0.143,0.241,0.287)
descae<-c(0.099,0.192,0.304)
elyscr<-c(0.177,0.217,0.316)
fesrub<-c(0.165,0.307,0.310)
kobmyo<-c(0.063,0.077,0.147)
trispi<-c(0.135,0.190,0.240)
cirsco<-c(0.433,0.555,0.656)
geuros<-c(0.124,0.312,0.405)
oxydig<-c(0.437,0.543,0.654)
senfre<-c(0.158,0.189,0.316)
silaca<-c(0.293,0.606,0.732)

# Figure 2
par(mfrow=c(4,3), oma=c(3,3,1,1) +0.1,mar=c(0.75,0.75,0.75,0.75) +0.1,mgp=c(2,1,0),xpd=NA)
par(xaxs="i")

barplot(carnar, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3));title(main="C. nardina", font.main=3,adj=0.05,cex.main=0.90,line=-1); text(x=0.155,y=carnar[2], "*", pos=3, offset=0.1, cex=1.5)

barplot(carpha, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="C. phaeocephala", font.main=3, adj=0.05, cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=carpha[3], "*", pos=3, offset=0.1, cex=1.5)

barplot(descae, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="D. cespitosa", font.main=3,adj=0.05, cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=descae[3], "**", pos=3, offset=0.1, cex=1.5)

barplot(elyscr, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3));title(main="E. scriberneri", font.main=3,adj=0.05, cex.main=0.90,line=-1); text(x=0.255,y=elyscr[3], "*", pos=3, offset=0.1, cex=1.5)

barplot(fesrub, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="F. rubra", font.main=3,adj=0.05, cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.155,y=fesrub[2], "**", pos=3, offset=0.1, cex=1.5)

barplot(kobmyo, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="K. myosuroides", font.main=3,adj=0.05,cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=kobmyo[3], "**", pos=3, offset=0.1, cex=1.5)

barplot(trispi, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, cex.names=0.75, ylim=c(0,0.8),xlim=c(0,0.3));title(main="T. spicatum", font.main=3,adj=0.05,cex.main=0.90,line=-1); text(x=0.255,y=trispi[3], "**", pos=3, offset=0.1, cex=1.5)

barplot(cirsco, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="C. scopulorum", font.main=3,adj=0.05, cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=cirsco[3], "**", pos=3, offset=0.1, cex=1.5)

barplot(geuros, space=0.01, width=0.1, col=c("white","gray70","gray50"), cex.main=1, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="G. rossii", font.main=3,adj=0.05, cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=geuros[3], "*", pos=3, offset=0.1, cex=1.5)

barplot(oxydig, space=0.01, width=0.1,names.arg=c("A","A+P","FULL"), col=c("white","gray70","gray50"), cex.main=1, cex.names=0.75, ylim=c(0,0.8),xlim=c(0,0.3));title(main="O. digyna", font.main=3,adj=0.05, cex.main=0.90,line=-1); text(x=0.255,y=oxydig[3], "**", pos=3, offset=0.1, cex=1.5)

barplot(senfre, space=0.01, width=0.1, names.arg=c("A","A+P","FULL"), col=c("white","gray70","gray50"), cex.main=1,cex.names=0.75, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="S. fremontii", font.main=3,adj=0.05, cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=senfre[3], "**", pos=3, offset=0.1, cex=1.5)

barplot(silaca, space=0.01, width=0.1, names.arg=c("A","A+P","FULL"), col=c("white","gray70","gray50"), cex.main=1, cex.names=0.75, ylim=c(0,0.8),xlim=c(0,0.3),yaxt='n');title(main="S. acaulis", font.main=3,adj=0.05,cex.main=0.90,line=-1); axis(side=2, at=c(0,0.2,0.4,0.6,0.8), labels=FALSE); text(x=0.255,y=silaca[3], "**", pos=3, offset=0.1, cex=1.5)

mtext("Model",side=1,outer=TRUE,cex=0.8,line=1.6)
mtext("Adj. D Squared Value",side=2,outer=TRUE,cex=0.8,line=1.2)
