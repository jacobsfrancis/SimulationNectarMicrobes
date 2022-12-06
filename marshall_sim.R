#this script is an individual-based simulation model for nectar microbe dispersal and growth
#the aim of the script is to test a hypothesis:
    #expanding landscapes 1) allow greater gamma diversity
                        # 2) favor dispersive species
                        # 3) increase beta diversity

    #AND
    #contracting landscapes 1) allow less gamma diversity
                          # 2) favor species tolerant of others OR that exclude others
                          # 3) decrease beta diversity

#the model simulates flowers, insects, and microbes - iterates dispersal over days
#and each day a 24 hour long hourly generalized Lotka-Volterra model is simulated for within
  #flower dynamics

#the script was written by Marshall McMunn - mmcmunn@gmail.com - 2021




#community size parameters for simulation
    days = 50
    hours = days*24
    flowers = 2000
    insectMax = 40
    insectMin = 10
#MICROBE PARAMETERS
    #number of yeast strains
    nYeast <- 2
    #carrying capacity (density per uL)
    yeast_K_peruL <- 10000
    #rate of growth
    yeast_r <- 0.1
    
    #percapita impact of each yeast species on floral signal
    yeast_s <- c(1/10000,1/20000)  
    if(nYeast != length(yeast_s)) {
      stop("Every yeast species must have param for yeast_s")
    }
    
    
    #number of bacteria strains
    nBact <- 3
    #carrying capacity (density per uL)
    bact_K_peruL <- 100000
    #rate of growth
    bact_r <- 0.15
    
    #percapita impact of each bact species on floral signal
    bact_s <- c(1/10000,1/20000)  
    if(nBact != length(bact_s)) {
      stop("Every bacteria species must have param for bact_s")
    }
    
    
    #draw r - growth rates for microbes
    rVector<-abs(c(rnorm(nYeast, yeast_r , 0.02  ) , rnorm(nBact , bact_r, 0.02)))
    
    #create interaction matrix
    nSpecies <- nYeast + nBact
    
    #draw lotka-volterra alpha parameters
    interactionMatrix <- matrix(data = rnorm(nSpecies^2 , 1, 0.25) , nrow = nSpecies, ncol = nSpecies)
    
    #adjust per capita interaction matrix for interkingdom yeast effect on bacteria x10 effect and 
    #bacteria effect on yeast /10 effects
    interactionMatrix[1:nYeast ,(nYeast+1):nSpecies ]<-interactionMatrix[1:nYeast ,(nYeast+1):nSpecies ]*10
    interactionMatrix[(nYeast+1):nSpecies  ,1:nYeast]<-interactionMatrix[(nYeast+1):nSpecies  ,1:nYeast]/10
    #set intraspecifc interactions to 1
    diag(interactionMatrix) <- 1
    
    
    
    #physical traits of microbes for stick, vector survival, release from vector
      #all stored in one dataframe
        microbeParams<-data.frame(microbeID = LETTERS[1:nSpecies],
          yeastBact = c(rep("yeast", nYeast) , rep("bacteria", nBact)) ,
        microbeProbStick = c(rep(0.3 , nYeast) , rep(0.8, nBact)),
        microbeVectorSurvive = c(rep(0.5 , nYeast) , rep(0.2, nBact)),
        microbeProbRelease = 1- c(rep(0.3 , nYeast) , rep(0.8, nBact)),
        r = rVector,
        K_peruL = c(rep(yeast_K_peruL , nYeast) , rep(bact_K_peruL, nBact))
        )
    

    
#INSECT PARAMETERS
      #number of insect species
      nInsect <- 2
      
      insectParams<-t(data.frame(
            #how many flowers do they each visit?
            insectVisitRate = c(3 , 6) ,
            #how relatively common are they (sum to 1)
            relativeOccurrences = c(0.3 , 0.7),
            #proportion left after visit
            visitKnockDown= c(0.8 , 0.2)))
            colnames(insectParams) <- c("thrips" , "bee")
    
    
#INSECT X MICROBE PREFERENCE
        #a dataframe with dimensions rows = # microbes, cols = # insects
        #empty flower visit prob 1, if deterrence  microbe value less than 1
        #if attraction, greater than 1
        prefMatrix <- data.frame(   c(1, 1, 0.8, 0.8, 0.8) , c(2, 2, 0.2, 0.2, 0.2) )
        colnames(prefMatrix) <- colnames(insectParams)
        rownames(prefMatrix) <- LETTERS[1:nSpecies] 
  
    
    
#FLOWER PARAMETERS
      flowerParams <- t(data.frame(
              flowerLifespan = c(6 , 3),
              #relative abundance
              flowerRelOccur = c(0.8 , 0.2),
              #nectar Volume
              flowerStartVol = c(50 , 20)))
              colnames(flowerParams) <- c("epilobium" , "stachys") 



    
#INSECT X PLANT PREFERENCE
    #a dataframe with dimensions rows = # plants, cols = # insects
    #if attraction, greater than 1
    prefMatrixFlow <- data.frame( c(0.4 , 0.6) ,  c(0.6, 0.2))
    colnames(prefMatrixFlow) <- colnames(insectParams)
    rownames(prefMatrixFlow) <- colnames(flowerParams)
    
    
#DERIVED PARAMETERS AND SETUP
    #flower opening functions
    #linear functions, derived from
        #total flowers
        #minimum flowers per day
        #duration of simulation
    #magnitude of slope
      #based on integral and minimum
    
    #################
    #function for flower change - name decrease OR increase OR constant
    
    #minimum for flowers per day
    flowerMinPerDay <- 20
    ##################
    #set increasing or decreasing
    x <- "increase"
    ####################
          areaOfTriangle<-flowers-flowerMinPerDay*days
          maxFlowerPerDay <- areaOfTriangle / (0.5*days)
    if(x == "increase"){
     flowerFunc<-function(t){ floor(((t/days)*(maxFlowerPerDay)) + flowerMinPerDay) }
    }
          
    if(x == "decrease"){
      flowerFunc<-function(t){ floor((maxFlowerPerDay) -
          ((t/days)*(maxFlowerPerDay)) + flowerMinPerDay)}
    }
          
          if(x == "constant"){
            flowerFunc<-function(t){floor(((flowers)/days))+t-t}
          }
#ctest and look at number of flowers per day
flowerFunc(1:100)

    #the number of insects in time changes proportionally with flowers
    #if this is not the case, then dilution effect is dominant, and there are less visitors per flower
    insectVector <- ceiling(flowerFunc(1:100)/max(flowerFunc(1:100))*(insectMax-insectMin))+insectMin


    #derived insect parameters (no need to edit)
    insectNames <- sample(colnames(insectParams) , replace = TRUE , prob = insectParams["relativeOccurrences" , ], size = insectMax)
    
    #a number that well exceeds a possible poisson draw within reason
    maxFlwrsPerDay<-max(rpois(10000, max(insectParams["insectVisitRate" , ])))*2
    
    #each flower is an slice in an array and contains a matrix where row = time and column = microbe
    #begins filled with 0's
    flow_MicrobeComm_HOUR<-array(data = 0, dim = c(hours, nSpecies, flowers))
    colnames(flow_MicrobeComm_HOUR) <- LETTERS[1:nSpecies]
    
    flow_MicrobeComm_DAY<-array(data = 0, dim = c(days, nSpecies, flowers))
    colnames(flow_MicrobeComm_DAY) <- LETTERS[1:nSpecies]
    
    #each insect is an element in the list and contains a dataframe where row = time and column = microbe
    insect_MicrobeComm<-array(data = 0, dim = c(days, nSpecies, insectMax))
    colnames(insect_MicrobeComm) <- LETTERS[1:nSpecies]
        
    #flower species vector
        flow_PlantID<-sample(colnames(flowerParams) , size = flowers , replace = TRUE, prob = flowerParams[ "flowerRelOccur", ])
      

    #flower age - a dataframe where row = time and column = flowers 
        #set to 0 until flower opens (then 1). Set to NA upon death
        flow_Age <- matrix(data = 0, nrow = days , ncol = flowers  )
        
    #create nectar volume data.frame - fill with 0's
        flow_nectarVol <- data.frame(matrix(0 , nrow = days, ncol = flowers))
        colnames(flow_nectarVol)<-1:flowers
        
        
    #which flowers start open
        startOpen <- sample(1:ncol(flow_Age), size = ceiling(flowerFunc(1)) )
        
      #add microbes randomly to X % of flowers (set here)
        startWithMicrobes <- sample(startOpen,floor(length(startOpen)*0.6))
    
      #open the flowers
        flow_Age[1 , startOpen] <- 1

    #add the microbes - everything starts at 100
              addThese<-sample(1:5, replace = TRUE, length(startWithMicrobes))
           for(f in 1:length(startWithMicrobes)){ 
             
             flow_MicrobeComm_DAY[1,addThese[f],startWithMicrobes[f] ] <-100
           }
          for(f in 1:length(startWithMicrobes)){ 
                
                flow_MicrobeComm_HOUR[1,addThese[f],startWithMicrobes[f] ] <-100
          }
              
      #nectar volumes
              flow_nectarVol <- flowerParams[ "flowerStartVol",flow_PlantID ]

              
              
#begin a loop, every iteration is one timestep (a day)
# 1, microbes grow
# 1, flowers age, microbes set to NA if lifespan is reached
# 3, insects draw a poisson number of ordered visits
# 4, microbes are reduced when visited and then deposited at the next flowers that day

for(i in 2:days){
  #if a flower is alive (above 0 and non-NA), add 1 day to age, if not, remain NA
#before dispersal
  #index columns of flowers that were occupied in i-1

  occupiedFlowersPrevious<-which(colSums(flow_MicrobeComm_DAY[i-1,,])>0)
  
  
  #microbes grew since the last day ended - generalized lotka-volterra logistic growth
for(f in occupiedFlowersPrevious){
  for(h in 1:24){
    for(s in 1:nSpecies){
      flow_MicrobeComm_HOUR[(i-2)*24+h+1 , s, f] <- flow_MicrobeComm_HOUR[(i-2)*24+h  ,s, f ] + flow_MicrobeComm_HOUR[(i-2)*24+h  ,s , f] * rVector[s] * (1- (sum((interactionMatrix[, s] * flow_MicrobeComm_HOUR[(i-2)*24+h  , , f ])) / (microbeParams[s,"K_peruL"]*flow_nectarVol[f])) )
      if(flow_MicrobeComm_HOUR[(i-2)*24+h+1 , s, f]<1){
        flow_MicrobeComm_HOUR[(i-2)*24+h+1 , s, f]<-0
      }
      }
  }
  flow_MicrobeComm_DAY[i , ,f ]<-flow_MicrobeComm_HOUR[(i-1)*24 + 1 , , f]
  }
  
  #some convenient checks to troubleshoot the LV model
  #sum(colSums(flow_MicrobeComm_HOUR[1,,]))
  #sum(colSums(flow_MicrobeComm_HOUR[25,,]))
  #sum(colSums(flow_MicrobeComm_DAY[2,,])>0)
  #flow_MicrobeComm_HOUR[,,occupiedFlowersPrevious[1]]
  
  
  #flowers age
  flow_Age[i,] <- ifelse(flow_Age[i-1,]==0 , 0, flow_Age[i-1,]+1 )

  #flowers die
  
  areDead <- which(flow_Age[i , ] > flowerParams[ "flowerLifespan", flow_PlantID])
  flow_Age[i,areDead]  <-  NA
  flow_MicrobeComm_DAY[i,,areDead]  <- 0
  
  #new flowers open
  newFlowers <- sample( which(flow_Age[i,]==0) , size = flowerFunc(i) , replace = FALSE)
  flow_Age[i , newFlowers] <- 1
  flow_nectarVol[newFlowers] <- flowerParams[ "flowerStartVol", flow_PlantID[newFlowers] ]

  #insects visit - 
  #this loop is once for each insect active that day (i)
  
  #get weighted average of preference by microbe abundance
  availFlowers <- which(flow_Age[i , ]>0)
  
  #thrips microbe preference weights for day i
  microbeVisitWeightsThrips<-list()
  for(f in availFlowers){
    temp <- data.frame(abund = flow_MicrobeComm_DAY[i, , f], species = LETTERS[1:nSpecies] )
    if(sum(temp$abund)>0){
      microbeVisitWeightsThrips[[f]]<-weighted.mean( x = prefMatrix[ ,"thrips"], w = temp$abund   )}
    if(sum(temp$abund)==0){
      microbeVisitWeightsThrips[[f]]<-1
    }
  }
  microbeVisitWeightsThrips<-unlist(microbeVisitWeightsThrips)
  
  #bee microbe preference weights for day i
  microbeVisitWeightsBee<-list()
  for(f in availFlowers){
    temp <- data.frame(abund = flow_MicrobeComm_DAY[i, , f], species = LETTERS[1:nSpecies] )
    if(sum(temp$abund)>0){
      microbeVisitWeightsBee[[f]]<-weighted.mean( x = prefMatrix[ ,"bee"], w = temp$abund   )}
    if(sum(temp$abund)==0){
      microbeVisitWeightsBee[[f]]<-1
    }
  }
  microbeVisitWeightsBee<-unlist(microbeVisitWeightsBee)
  
  
  #insect j on day i
  #this loop is once for each insect active
  for(j in 1:insectVector[i]){
    
    
    #draw the number of flowers visited by each insect on this day
    NumVisits <- rpois(1, insectParams["insectVisitRate" , insectNames[j]])
    
    if(NumVisits>0){
      #loop through each visit and do things - knockdown existing microbes if they are there
                                              #label future flowers visited according to current flower and persistence
                                              #add intitial dispersers if the flower was empty
      #this is for each visit the insect does
      for(k in 1:NumVisits){
        
        #from above, when daily microbe-based flower preferences were calculated
        #this was computationally expensive inside the visit loop...
        microbeVisitWeights<-ifelse(insectNames[j]=="bee" , microbeVisitWeightsBee, microbeVisitWeightsThrips)
        
        #generate weights for plant species preference
        plantIndex<-flow_PlantID[availFlowers]
      
        plantVisitWeights <- prefMatrixFlow[plantIndex ,insectNames[j] ]
        finalWeights <- plantVisitWeights*microbeVisitWeights
        
        
        #this selects the flower for this particular visit
        thisVisit <- sample(x =  which(flow_Age[i , ]>0), size = 1,prob = finalWeights )
        
        #step 1 - pickup microbes
        insect_MicrobeComm[ i , , j ] <- insect_MicrobeComm[ i , , j ] + (flow_MicrobeComm_DAY[i , ,thisVisit ]/flow_nectarVol[thisVisit]) * microbeParams[ ,"microbeProbStick" ]
        
        #step 2 - knockdown microbes in flower
        flow_MicrobeComm_DAY[i , ,thisVisit ] <- insectParams[ "visitKnockDown" , insectNames[j] ]*flow_MicrobeComm_DAY[i , ,thisVisit ]
        
        #step 3 - dropoff microbes to flower
        flow_MicrobeComm_DAY[i , ,thisVisit ] <- flow_MicrobeComm_DAY[i , ,thisVisit ] + (insect_MicrobeComm[ i , , j ] * microbeParams[ ,"microbeProbRelease" ])
        
        #step 4 - microbes left over after visit (picked up, remaining from previous) experience mortality
        insect_MicrobeComm[ i , , j ] <- insect_MicrobeComm[ i , , j ] * microbeParams[ , "microbeVectorSurvive"]
        
        #step 5 - set anything that falls below a cell each visit to 0
        insect_MicrobeComm[ i , , j ][insect_MicrobeComm[ i , , j ]<1] <- 0
      }}}
  #copy over the new microbe totals after visitation
  for(f in 1:flowers){
  flow_MicrobeComm_HOUR[(i-1)*24+1 , ,f ] <- flow_MicrobeComm_DAY[i , ,f ]
  
  }
  
}
            

    
    library(ggplot2)
    library(reshape2)
    library(viridis)
              
              #flower microbe composition
              compLongDF<-melt(flow_MicrobeComm_DAY)
              colnames(compLongDF) <- c("day" , "microbe" , "plant" , "abundance")
              compLongDF$microbePresent<-compLongDF$abundance>0
              microbeTable<-with(compLongDF , table(day, microbe, microbePresent))[,,2]
              
              propComp<-microbeTable/rowSums(flow_Age>0, na.rm = TRUE)
              propCompLong<-melt(propComp)
              colnames(propCompLong)[3]<- "proportionOfOpen"
              ggplot(data = propCompLong, aes(x = day, y = proportionOfOpen  , color = microbe )) + geom_line() + 
                theme_bw()+ggtitle("microbe propotion occupancy") + ylab("proportion flowers occuppied")
              ggsave("figures/propOccupancy_byMicrobe.pdf")
              
              longData<-melt(flow_Age)
              longData<-longData[longData$value!=0,]
              
              #flower ages
              ggplot(longData, aes(x = Var2, y = Var1)) + 
                geom_tile(aes(fill=value)) + 
                scale_fill_gradient(low="grey90", high="red") +
                labs(x="flowerID", y="day", title="Matrix") +
                theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                                   axis.text.y=element_text(size=9),
                                   plot.title=element_text(size=11))+scale_fill_viridis()+ylim(days, 0)+
                theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
                ggtitle("flower age through time")
              ggsave("flower_age.png", path = "figures/", width = 48, height = 34, units = "in")
              
              
              #a look inside a few individual flowers during their open period
              
              flowersWithMicrobes<-which(apply(flow_MicrobeComm_HOUR[ , , ] , 3, sum)>0)
              plottingFlowers<-sample(flowersWithMicrobes, 20)
              
              abunds<-list()
              for(p in plottingFlowers){
                tempFlower<- flow_MicrobeComm_HOUR[ , , p]
                times<-which(rowSums(tempFlower)>0)
                
                if(length(times)>1){
                  abunds<-data.frame(tempFlower[rowSums(tempFlower)>0 , ])
                  abunds<-cbind(times, melt(abunds))
                }
                
                if(length(times)==1){
                  abunds<-t(data.frame(tempFlower[rowSums(tempFlower)>0 , ] , c(0,0,0,0,0)))
                  rownames(abunds)<-1:2
                  abunds<-cbind(times, melt(abunds))
                  abunds<-abunds[,-2]
                }          
                
                
                colnames(abunds) <- c("time_hour" , "species" , "abundance")
                tempName<-paste("figures/flower_",p, ".png",sep="" )
                
                ggplot(abunds, aes(x = time_hour, y = abundance, colour = species))+geom_line()
                ggsave(tempName)
                
              }
              
              
              #compare plant species
              #flower microbe composition
              compLongDF<-melt(flow_MicrobeComm_DAY)
              colnames(compLongDF) <- c("day" , "microbe" , "plant" , "abundance")
              compLongDF$plantSpecies<-flow_PlantID[match(compLongDF$plant , 1:length(flow_PlantID))]
              compLongDF$microbePresent<-compLongDF$abundance>0
              microbeTable<-with(compLongDF , table(day, microbe,plantSpecies,microbePresent))[,,,2]
              propComp<-microbeTable
              spCounts<-matrix(nrow = days, ncol = length(unique(flow_PlantID)))
              for(s in 1:length(unique(flow_PlantID))){
                species<-unique(flow_PlantID)[s]
                spCounts[,s]<-rowSums(flow_Age[ , flow_PlantID==species]>0, na.rm = TRUE)
                propComp[,,s]<-microbeTable[,,s]/spCounts[,s]
              }
              
              propComp<-microbeTable/rowSums(flow_Age>0, na.rm = TRUE)
              
              propCompLong<-melt(propComp)
              colnames(propCompLong)[4]<- "proportionOfOpen"
              ggplot(data = propCompLong, aes(x = day, y = proportionOfOpen  , color = microbe)) + geom_line() + 
                theme_bw()+ggtitle("microbe propotion occupancy - by plant species") +facet_grid(~plantSpecies) + ylab("proportion flowers occuppied")
              ggsave("figures/propOccupancy_byMicrobe_byPlant.pdf")
              
              
              #gamma, alpha, beta diversity 
              #gamma
              compLongDF<-melt(flow_MicrobeComm_DAY)
              colnames(compLongDF) <- c("day" , "microbe" , "plant" , "abundance")
              compLongDF$microbePresent<-compLongDF$abundance>0
              gammaDiv<-rowSums(with(compLongDF , table(day,microbe, microbePresent ))[,,2]>1)
              
              #avgAlpha for each day
              temp<-apply(flow_MicrobeComm_DAY , 2, rbind)
              temp <- cbind(1:days,sort(rep(1:flowers, days), decreasing=FALSE), temp)
              colnames(temp) <- c("day" , "plant" , LETTERS[1:nSpecies])
              temp<-data.frame(temp)
              flow_Age_Long<-melt(flow_Age)
              
              temp$flowAge<-flow_Age_Long[match(paste(flow_Age_Long[,1], flow_Age_Long[,2]) , paste(temp$day , temp$plant)) ,3 ]
              temp[is.na(temp$flowAge) , "flowAge" ] <-0
              justOpen<-temp[ temp$flowAge>0 , ]
              justOpen$alphaDiv<-rowSums(justOpen[,3:(3+(nSpecies-1))]>0 )
              alphaDiv<-with(justOpen , tapply( alphaDiv, day ,mean ))
              
              #occupancy
              propOccupy<-with(justOpen , tapply(alphaDiv>0 ,day ,mean ))
              
              #beta div for each day
              library(vegan)
              brayDist_byDay<-list()
              for(d in 1:days){
                
                oneDayOpen<-justOpen[justOpen$day==d,3:(3+(nSpecies-1))]
                distBray<-vegdist(oneDayOpen)
                brayDist_byDay[[d]] <-mean(distBray, na.rm = TRUE)
              }
              betaDiv<-unlist(brayDist_byDay)
              
              

              
              #plot all 4
              divDF<-data.frame(cbind(gammaDiv , alphaDiv , betaDiv, propOccupy))
              
              divPlot<-melt(divDF)
              divPlot <- cbind(1:days , divPlot)
              colnames(divPlot)<- c("day" , "statistic" , "value")
              ggplot(data = divPlot , aes(x = day, y = value, colour = statistic )) +
                geom_line() +facet_grid(~statistic)
              ggsave("figures/diversity_throughTime.png")
              
              #correlation tests 
              #this currently only works when 5 species##
              #abundance
              cor(temp[,3:7])
              #presence
              cor(temp[,3:7]>0)
              #actual tests
              cor.test(temp[,7], temp[,3])
              cor.test(temp[,7], temp[,4])
              cor.test(temp[,7], temp[,5])
              cor.test(temp[,7], temp[,6])
              
              
              #plot the microbe abundance on the vectors each day
              #choose 20 randomly to plot
              insect_Mic_Long<-melt(insect_MicrobeComm)
              colnames(insect_Mic_Long) <- c("day" , "microbe" , "insect" , "abundance")
              for(g in sample(1:insectMax, 20)){
                thisInsect<-insect_Mic_Long[insect_Mic_Long$insect==g , ]
                
                
                ggplot(data = thisInsect , aes(x = day, y = abundance, fill = microbe))+
                  geom_bar(stat = "identity", position = "dodge")+ggtitle(paste(insectNames[g] , g))
                
                tempName<-paste("figures/insect_",g, ".png",sep="" )
                
                ggsave(tempName)
              }
              
      