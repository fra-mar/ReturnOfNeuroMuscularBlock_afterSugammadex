# Goal: make model suitable for Uppmax
library(dplyr)

# Conversion constants
# Rocu and suggis mg/L to microM
roc2microM = 1.887505
sug2microM = 0.499451
# Fraction of parent RocBromide which is rocuronium
# rocBromide2roc = 529.8/609.7 # I considered this in initial versions  

buildParams = function(AGE, BW, CR, RAC, SEV, 
                       bolusROC, i_ROC, durInf,
                       bolusSUG, nSub) {
  
  d0ROC = bolusROC * BW            #induction dose ROC in mg (mgKg-1 *Kg)
  d0ROC = d0ROC  * roc2microM #conversion to micromoles
  infROC = i_ROC * BW / 60#* durInf     # amount ROC per minute
  infROC = infROC  * roc2microM #conversion to micromoles
  d0SUG = bolusSUG * BW            # Reversal dose SUG in mg (mgKg-1 *Kg)
  d0SUG = d0SUG * sug2microM  # Conversion mg to micromoles
  tSUG  = durInf                # timestamp when SUG administered (hours)
  #simDur = (durInf+3) * 60    #infusion duration in minutes + n hour follow up
  
  forEvents = tibble(d0ROC=d0ROC, infROC=infROC, durInf=durInf, 
                d0SUG=d0SUG, tSUG=tSUG)#, simDur=simDur)
  
  CLage   <- 1 + -0.00678 * (AGE - 43)
  tvCLro    <- CLage* 0.269 *(BW/70)^0.75
  V1cr    <- exp(-0.00143 * (CR-119))
  tvV1ro    <- V1cr * 4.73 * (BW/70)^1
  Q2roRac <- ifelse(RAC=='nonAsian',1, 1-0.212)
  Q2ro    <- Q2roRac * 0.279 * (BW/70)^0.75
  V2roAge <- exp(0.00613 * (AGE-43))
  tvV2ro    <- V2roAge * 6.76 * (BW/70)^1
  
  # for SUG
  CLsuBW <-   1 + 0.00378 * (BW-74.5)
  REN <-      (2 * CR/(CR+119))^1.29
  tvCLsu <-   CLsuBW * (REN * 0.093)
  V1suBW <-   1 + -0.00354 *(BW- 74.5)
  V1suRAC <- ifelse(RAC=='nonAsian',1, 1-0.16)
  tvV1su <-   V1suBW * V1suRAC * 4.42 * (BW/70)^1
  tvQ2su <-   0.206 * (BW/70)^0.75
  V2suCR <-   exp( -0.00305 * (CR-119))
  tvV2su <-   V2suCR * 6.35 * (BW/70)^1
  
  # PD model ROC
  ke0Sev <- ifelse(SEV==TRUE, 1-0.567, 1) # 1 no sevo,  1 +  -0.567 if sevo
  tvke0    <- ke0Sev * 0.134 * (BW/70)^-0.25
  EC50sev <- ifelse(SEV==TRUE, 1-0.395, 1) # 1 no sevo,  1 +  -0.395 if sevo
  tvEC50   <- EC50sev * 1.62 # microM
  #PD model SUG
  tvks <- exp(-3.43) * (BW/70)^-0.25
  tvHill <- 7.52
  tvE0 <- 100 #1.04 * 100 #Arbitrary set to 100
  
  # eta
  eta.E0 <- 0.01224571 # CV 11.1%
  eta.V1ro <- 0.05509785 # CV% 23.8%
  eta.CLro <- 0.09982362 # CV% 32.4%
  eta.V2ro <- 0.09865367 # CV% 32.2%
  eta.ke0 <- 0.16032217 #CV % 41.7%
  eta.EC50 <- 0.06015486 # CV 24.9%
  eta.Hill <- 0.15608110 # CV 41.1%
  eta.CLsu <- 0.04895777 # CV 22.4%
  eta.ks   <- 0.83273519 # CV 114%
  
  # Population estimates
  parsPop = tibble(
    CLro = tvCLro,
    Q2ro = Q2ro,
    V1ro = tvV1ro,
    V2ro = tvV2ro,
    V1su = tvV1su,
    V2su = tvV2su,
    CLsu = tvCLsu,
    Q2su = tvQ2su,
    ke0 = tvke0,
    EC50 = tvEC50,
    Hill = tvHill,
    E0 = tvE0,
    ks = tvks
  )
  # Log-normal random distributed parameters
  # To understand sqrt(eta_XXX): 4.1.1 https://nmhelp.tingjieguo.com/V/4
  nSub = nSub - 1
  parsRnd = tibble(
    CLro = tvCLro * exp(rnorm(nSub,mean=0, sd=sqrt(eta.CLro))),
    Q2ro = rep(Q2ro, nSub),
    V1ro = tvV1ro * exp(rnorm(nSub,mean=0, sd= sqrt(eta.V1ro))),
    V2ro = tvV2ro * exp(rnorm(nSub,mean=0, sd= sqrt(eta.V2ro))),
    V1su = rep(tvV1su, nSub),
    V2su = rep(tvV2su, nSub),
    CLsu = tvCLsu * exp(rnorm(nSub,mean=0, sd=sqrt(eta.CLsu))),
    Q2su = rep(tvQ2su, nSub),
    ke0 = tvke0 * exp(rnorm(nSub,mean=0, sd=sqrt(eta.ke0))),
    EC50 = tvEC50 * exp(rnorm(nSub,mean=0, sd=sqrt(eta.EC50))),
    Hill = tvHill * exp(rnorm(nSub,mean=0, sd=sqrt(eta.Hill))),
    E0 = tvE0, # * exp(rnorm(nSub,mean=0, sd=sqrt(eta.E0))),
    ks = tvks * exp(rnorm(nSub,mean=0, sd=sqrt(eta.ks)))
  )
  
  params = rbind(parsPop, parsRnd)
  toReturn = list(forEvents=forEvents, 
                  params=params)
  
  return(toReturn)
}
