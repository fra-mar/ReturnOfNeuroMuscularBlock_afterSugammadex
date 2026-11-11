# This script contains strictly the model
library(rxode2)

mod <- rxode2({
  
  C1ro <- rocCentral / V1ro
  C2ro <- rocPerif   / V2ro
  #SUG, PK      
  C1su <- sugCentral / V1su
  C2su <- sugPerif / V2su
  #Complex, pK
  C1cx <- cpxCentral / V1su #complex same Kel and V1,2 as sug
  C2cx <- cpxPerif   / V2su
  k2cx <- 0.03404745 # in paper log(e) K2 = -3.38, ie exp(-3.38)
  k1cx <- 0.6090779 # in paper k1cx=k2cx/kd
  # ODEs
  #k1ro <- k10ro + k12ro + ke0
  #k1su <- k10su + k12su
  
  d/dt(rocCentral) <- C2ro * Q2ro - 
                      C1ro * CLro - C1ro * Q2ro +
                      C1cx * k2cx - C1ro * C1su * k1cx  # Interaction model
  d/dt(rocPerif)   <- C1ro * Q2ro - C2ro * Q2ro
  
  d/dt(sugCentral) <- C2su * Q2su -
                      C1su * CLsu - C1su * Q2su +
                      C1cx * k2cx - C1ro * C1su * k1cx   #Interaction model
  d/dt(sugPerif)   <- C1su * Q2su - C2su * Q2su
  
  d/dt(cpxCentral) <- C2cx * Q2su - C1cx * CLsu - C1cx *Q2su +
                      C1ro * C1su * k1cx - C1cx * k2cx   #Interaction model
  d/dt(cpxPerif)   <- C1cx * Q2su - C2cx * Q2su
  
  d/dt(Ceff)       <- ke0 * (C1ro - Ceff) - ks * C1su
  TOFr             <- E0 - (E0 * Ceff^Hill / (EC50^Hill + Ceff^Hill))
})
