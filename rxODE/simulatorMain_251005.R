# Goal: make model suitable for Uppmax
library(rxode2)
library(dplyr)

set.seed(1234)
setRxThreads(percent=75)
print(paste0('RxThreads: ', getRxThreads()))

rm(list=ls()) #Uncomment to erase all objects from environment
pwd = getwd()
print(paste('Working on PWD:',pwd))
# Create a folder for outputs
output_dir = "simulations"
if (!file.exists(output_dir)) {dir.create(output_dir)}

# Load the model, the parameters builder and scenarios
source('kleijnModel_251011.R')           
source('kleijnParamBuilder_251011.R')
scenarios = tibble(read.csv('scenarios.csv'))

# ---------------LOAD scenario and NUMBER of SUBJECTS --------------------
args = as.numeric(commandArgs(trailingOnly = TRUE))
i = args[1]
nSub = args[2]

s = scenarios[i,]
print(paste0('Working on scenario: ', s$sceName))

dirName = paste0(pwd,"/",output_dir,"/","scenario_",s$sceName)
if (!file.exists(dirName)) {dir.create(dirName)}

# Calculate parameters, both for event table and model params
allParams = buildParams(AGE = s$AGE, BW = s$BW, CR = s$CR, RAC = s$RAC, 
                        SEV = s$SEV, bolusROC = s$d0ROC, i_ROC = s$infROC, 
                        durInf = s$infDur, bolusSUG = 0, nSub=nSub)

# fe, for events. params, model parameters
fe      = allParams$forEvents
params  = allParams$params
follUp  = 3   # follow up after SUG admin, in HOURS
shrinkSamp = .1 # reduction sampling to save memory: 60*.1 = .1 Hz
lengthoutROC = as.integer(fe$durInf*60 * 60 * shrinkSamp)
lengthoutSUG = as.integer(follUp*60    * 60 * shrinkSamp)

#Create event table
ev <- et( amountUnits = 'micromol',timeUnits = "minutes") %>%
      et( amt = fe$d0ROC, cmt = "rocCentral") %>%    # Start dose ROC
      et( amt = fe$infROC * fe$durInf * 60,  # InfusionROC mg/min * total min
          cmt = "rocCentral", 
          rate = fe$infROC) %>% # Rate as mg/min
  
  et(0, fe$durInf*60, length.out = lengthoutROC) #simDur*60, in 10 seconds

# Simulate from ROC admin just before SUG admin
simROC = rxSolve(mod, params, ev)


#--------------------------------------------------2nd part, SUG
# Get the last state
lastState = simROC |> select(c(sim.id, rocCentral,rocPerif,sugCentral,sugPerif,
                              cpxCentral,cpxPerif,Ceff,TOFr  )) |>
  group_by(sim.id) |> slice(n())

# Intialize final tibble
sim = tibble(simROC)
# Spare calculated SUG doses
sugDoses = c()

# Iterate through subjects 
for (id in lastState$sim.id) {
  idState = lastState |> filter(sim.id == id)
  idPar = params[id,]
  lastTOFr = idState$TOFr
  #Decide about SUG dose  (mg/Kg)
  if (lastTOFr > 90) {
    d0SUG = 0}
  else {
    if (length(grep("fixed", s$d0SUG)) == 1) {
      d0SUG = as.numeric(substr(s$d0SUG,6,nchar(s$d0SUG))) / s$BW
    }
    else if (s$d0SUG == "var"){
      d0SUG = ifelse(lastTOFr < 10,4.0,2.0)
    }
    else{
      print('found nothing about SUG dose')
    }
  }
   
  
 
  sugDoses = c(sugDoses, d0SUG)
  # Conversion total mg to micromol
  d0SUG = d0SUG * s$BW * 0.499451 
  
  # define initial values for mod$state
  inits = c(rocCentral = idState$rocCentral,
            rocPerif   = idState$rocPerif,
            sugCentral = idState$sugCentral,
            sugPerif   = idState$sugPerif,
            cpxCentral = idState$cpxCentral,
            cpxPerif   = idState$cpxPerif,
            Ceff       = idState$Ceff)
  
  evSUG <-  et( amountUnits = 'micromol',timeUnits = "minutes") %>%
            et(amt = d0SUG, cmt = "sugCentral") %>%
            et(0, follUp*60, length.out = lengthoutSUG) # follow up 3 hours
  
  # Run simulation for one subject with ID = id
  simSUG = rxSolve(mod,idPar,evSUG, inits=inits)
  
  simulatedSUG = tibble(simSUG)
  # Reset times sin simSUG to last time before SUG (last time in simROC)
  simulatedSUG$time = simulatedSUG$time + max(simROC$time)
  # Add sim.id as first column
  simulatedSUG = tibble(cbind(data.frame(sim.id = rep(id,follUp*6)), 
                              simulatedSUG))
  # Concatenate to general sim tibble
  sim = rbind(sim, simulatedSUG)
}

# Correct negatives Ceff values and NaN TOFr, clean cols and...
#... put sim.id and times in right order
sim = sim |> mutate(Ceff = replace(Ceff, Ceff <=0, 1E-5),
                    TOFr = replace(TOFr, is.na(TOFr), 100)) |>
  select(c("sim.id", "time", "C1ro", "C2ro", 
           "C1su", "C2su", "C1cx", "C2cx", "Ceff", "TOFr")) |>
  arrange(sim.id, time)

# save relevant data________________________________________________
patientRegime = list(scen = s, forEvents = fe)
saveRDS(patientRegime, paste0(dirName,'/','patientRegime.RDS'))

write.csv(simROC$params,  paste0(dirName,'/','params.csv'), 
          row.names = FALSE )
write.csv(sim, paste0(dirName,'/','simData.csv'), 
          row.names = FALSE )



  

