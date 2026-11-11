library(dplyr)
rm(list = ls())


#Chose a test scenario
scenDirs = list.dirs('simulations', recursive = FALSE)
for (scenDir in scenDirs)  {
  print(paste0("Working on ",scenDir))
  
  patData = readRDS(paste0(scenDir,'/','patientRegime.RDS'))
  simData = tibble(read.csv(paste0(scenDir,'/','simData.csv')))
  
  timeSUG = patData$forEvents$durInf * 60 # conversion hours -> minutes
  
  simData = simData |> filter((time>timeSUG))
  
  # Initialize a tibble for results
  subjectsIDs = unique(simData$sim.id)
  
  #Iterate through subjects
  for (i in c(1:length(subjectsIDs))){
    #print(paste0('Subject: ', i))
    subData = simData |> filter(sim.id==i) |>
      mutate(time = time - timeSUG)

    spntRecover = ifelse(subData$TOFr[1] >=90, TRUE, FALSE)
    
    if (all(subData$TOFr<90)){
      didRecover = FALSE
      t90         = -1
      isRNB       = FALSE
      tIniRNB     = -1
      durRNB      = -1
      lowestTOF   = -1
    }
    
    else {
      didRecover = TRUE
      # When TOFr 90% was first reached
      t90 = subData |> filter(TOFr>=90)
      t90 = as.numeric(t90[1,'time'])
      
      if (any( subData[subData$time > t90,]$TOFr<90) ) {
        isRNB     = TRUE
        under90   = subData |> rowwise() |>
                    filter((time>t90) & (TOFr<90))
        tIniRNB   = as.numeric(under90[1,'time'])
        durRNB    = max(under90$time) - min(under90$time)
        lowestTOF = min(under90$TOFr)
      }
      else {
        isRNB     = FALSE
        tIniRNB   = -1
        durRNB    = -1
        lowestTOF = -1}
    }
    
    #Assess for dips
    dips = subData |> select(c(time, TOFr)) |>
      mutate(dtof = c(0, diff(TOFr)))
    tempDip = dips |> filter( (TOFr<90) & (dtof<0) )
    tdip = ifelse(dim(tempDip)[1]==0,-1,tempDip$time[1])
    if (tdip!=-1){
      if (isRNB==FALSE){
        lowestDIP = as.numeric(tempDip$TOFr[dim(tempDip)[1]])
        isDIP = TRUE
      }
    }
    else{
      lowestDIP = -1
      isDIP = FALSE
    }
    
    if (!exists("calculated")){
      calculated = tibble(sim.id = i,
                          spntRecover = spntRecover,
                          didRecover = didRecover,
                          t90 = t90,
                          isRNB = isRNB,
                          tIniRNB = tIniRNB,
                          durRNB = durRNB,
                          lowestTOF = lowestTOF,
                          isDip = isDIP,
                          tDip = tdip,
                          lowestDIP = lowestDIP)}
    else{
      temp = tibble(sim.id = i,
                    spntRecover = spntRecover,
                    didRecover = didRecover,
                    t90 = t90,
                    isRNB = isRNB,
                    tIniRNB = tIniRNB,
                    durRNB = durRNB,
                    lowestTOF = lowestTOF,
                    isDip = isDIP,
                    tDip = tdip,
                    lowestDIP = lowestDIP)
      calculated = rbind(calculated, temp)
    }
    
    
    
  } #End of for loop
  write.csv(calculated, paste0(scenDir,'/','analyzedTOF.csv'), 
            row.names = FALSE)
  rm(calculated, simData)

}

