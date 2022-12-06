# Agent based model of pollination and microbial dispersal

# This model is designed to test how microbes which change the signal of a flower
# and change its quality might shape bees decisions to visit flowers, and shape their 
# own dispersal.  It also will eventually include parameters for floral defenses that 
# change microbial growth and bee preference.

# There are two types of agents:
  # Bees - Params:
        # 1) their decision threshold (i.e. signal above which they will visit a flower) 
        # 2) the quantity of microbes they are carying
        # 3) the quantity of pollen they are carying (for plant perspective questions)
        # 4) their innate preferences for nectar sugars (and alkaloids for later)
  # Flowers - Params:
        # 1) Their species (for plant perspective questions later)
        # 2) Their signal this is a linear measure of a signal (TODO think of multiple signal axes?)       
        # 3) their microbial carrying capacity
        # 4) Their nectar chemistry (sugars and alkaloids) 
        # 5) Their pollen load (for plant perspective questions later)
        # 6) Their pollen deposition rate (for plant perspective questions later)

  # Microbes - Params:
        
        # 1) their growth rate - r in traditional logistic growth equations
        # 2) their change to the signal, this is multiplied by the population number so it should be small
        # 3) their change to the sugar, this is multiplied by the population number so it should be small
        # 4) their change to the alkaloid, this is multiplied by the population number so it should be small

# At each time step microbes grow using a logarithmic growth function, then
# bees approach a flower (randomly) and decide whether to visit it
# If they visit it they remove nectar and microbes, deposit microbes, 
#   assess flower quality, and update their threshold.

# Some requirements - there must always be more flowers than pollinators!


#clear the work space (important for debugging code in production)

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(plyr)
library(parallel)
library(emmeans)
library(ggpubr)

#Critical functions ####

#Make agents####
BeeInit <- function(thresholdPar = c(0,100),        # we are going to specify the min and max bees have a threshold for each flower species for now just one 
                    microbesPar  = c(0,200),        # the microbial load for now just one species
                    pollenLoadPar = 0,       # the pollen load for now just one species
                    preferencesPar = c(30,.2)) # the initial floral preferences for sugar level and alkaloid
{
 
    time <- 1
    thresholdF1   <- runif(1,thresholdPar[1],thresholdPar[2])
    microbe1      <- runif(1,microbesPar[1],microbesPar[2])
    pollenLoadF1  <- pollenLoadPar[1]
    sugarPref     <- preferencesPar[1]
    alkaloidPref  <- preferencesPar[2]
    volPref       <- 20
    visit         <- 0
    ID <- paste(sample(size=10, x=letters),sep="",collapse="")
    LearnSig <- thresholdF1   
    LearnRew <- sugarPref
    LearnVol <- volPref
    LearnValue <- assess(flower=data.frame(necVol=volPref,sugar=sugarPref),
                         bee= data.frame(volPref=volPref,sugarPref=sugarPref))
      # this is the value function of what they expect for their standard flower, should be 1...
  
  
  return( data.frame(
    time = time,
    thresholdF1   = thresholdF1,
    microbe1      = microbe1,
    pollenLoadF1  = pollenLoadF1,
    sugarPref     = sugarPref,
    alkaloidPref  = alkaloidPref,
    volPref       = volPref,
    visit         = visit,
    ID = ID,
    LearnSig = thresholdF1,    
    LearnRew = sugarPref,
    LearnVol = LearnVol,
    LearnValue = LearnValue))                                             #these are the bees previous experiences.
}
#BeeInit()   # wont work without assess... it is defined later...        

# Initialize a flower
FlowerInit <- function(species ="flower1",             # flower species is coded as a character string
                       signalPar = c(0,100),              # this is a linear measure of a signal (TODO think of multiple signal axes?)       
                       microbesP = c(0,0),            #mean and SD of microbes starting
                       sugarPar = c(30,10),
                       alkaloidPar = c(.2,.01),  # mean and standard deviation for sugar then alkaloid
                       pollenLoadPar = 0,       # starting pollen load
                       pollenDepRate = 10,      # the amount of pollen deposited by a plant
                       maxNectar = 15,
                       refillRate = .025)             
{
  necVol <- runif(1,0,maxNectar)
  return( data.frame(
    time = 1,
    species = species,
    signal = runif(1,min = signalPar[1], max = signalPar[2]),
    microbesP1 = rnorm(1,microbesP[1],microbesP[2]),
    microbeK = 900 * necVol, #900 cells per ul # CHECK!
    necVol = necVol,
    maxNectar = maxNectar,
    refillRate=refillRate,
    sugar = rnorm(1, mean = sugarPar[1], sd = sugarPar[2]),
    alkaloid = rnorm(1, mean= alkaloidPar[1],alkaloidPar[2]),
    pollenF1 = pollenLoadPar[1],
    pollenDepRate = pollenDepRate,
    visit=0,
    ID = paste(sample(size=10, x=letters),sep="",collapse="")
  )
  )
}
FlowerInit()
#Initialize a microbe
MicrobeInit <- function(species = "microbe1",
                        growthRate = 1.2,
                        signalDelt = .002,
                        sugarDelt = .002,
                        alkaloidDelt = .001){
  return(data.frame (
    species = species,
    growthRate = growthRate,
    signalDelt = signalDelt,
    sugarDelt = sugarDelt,
    alkaloidDelt = alkaloidDelt
  ))
}

MicrobeInit()

#Time step functions####
# microbial growth function

# #for debugging functions...
# microbe <- MicrobeInit()
# flower <- FlowerInit()
# bee <- BeeInit()

grow <- function(microbe,flower){
              # define some times
              prev <- max(flower$time)
              current <- prev + 1
              flower[current,]<-flower[prev,]
              flower$time[current]<- current
              
              # flowers secrete nectar if they are below max and we calculate capacity
              
              
              if(flower[prev,"necVol"] < flower[prev,"maxNectar"]){
                flower[current,"necVol"] <- flower[prev,"necVol"]+flower[prev,"refillRate"]}
              if(flower[current,"necVol"] > flower[current,"maxNectar"]){flower[current,"necVol"] <- flower[current,"maxNectar"]}
              flower[current,"microbeK"] <- flower[current,"necVol"] * 500
              
              #first microbes grow using standard logistic growth 
              r <- microbe$growthRate
              P <- flower[prev,"microbesP1"] # REMEMBER TO CHANGE IF WE ADD MORE MICROBES!
              K <- flower[current,"microbeK"]
              flower[current,"microbeP1"] <- r*P*(1-(P/K))
              
              
              # then microbes change the signal of the plant
              flower[current,"signal"]<- flower$signal[1] + microbe$signalDelt * flower[current,"microbesP1"]
              
              #then microbes change the quality of the plant
              flower[current,"sugar"] <- flower$sugar[1] + microbe$sugarDelt* flower[current,"microbesP1"]
              flower[current,"alkaloid"] <- flower$alkaloid[1] + microbe$alkaloidDelt * flower[current,"microbesP1"]
              flower$visit[current] <- 0
              #then return the changed flower
              return(flower)
           }

# Value function
# calculate a reward function, here bees discount a flower if it has less than their 
# total capacity ie volPref very steeply. They have a linear increase in preference to
# as sugar increases from 0 - 45, then a steep decrease to 0 at 55
# they put twice as much weight on nectar sugar as they do volume

VolumeValFun <- function(bee,flower){
  # this gives us a value of ~0-1, it saturates towards 1 volume approaches the bees preference
  # the percentage before the expression in the exponent controls how fast that preference changes
  1/(1+exp(0.4*(-1*flower$necVol+bee$volPref-6)))
     } 
VolumeValFun(bee=data.frame(volPref=20), flower=data.frame(necVol=8))

SugarValFun  <- function(bee,flower){
  #this is linear increase in preference up to 45% sucrose, then steep drop to 0 preference at 55
  pref <-  bee$sugarPref[nrow(flower)]
  avail <- flower$sugar[nrow(flower)]
  
  if(avail <= pref & avail > 0)  {return(1/pref*avail)}
  if(avail > pref & avail <= pref+10) {return((-1/10)*(avail-(pref+10)))}
  if(avail<pref+10)             {return(0)}
  if(avail < 0 )               {return(0)}
}

SugarValFun(bee=data.frame(sugarPref=12), flower=data.frame(sugar=13))

assess <- function (bee,flower){
  # Syntax Value <- assess(bee=bee[currentTime,],flower=flower[currentTime,])
  
  volumeValue <- VolumeValFun(bee=bee[nrow(bee),], flower=flower[nrow(flower),])
  sugarValue  <- SugarValFun(bee=bee[nrow(bee),],flower=flower[nrow(flower),])
  
  # they put twice as much weight on nectar sugar as they do volume (i.e. rep 2)
  return(mean(c(rep(sugarValue,2),volumeValue)))}
  #return(mean(c(rep(sugarValue,2))))} #took out volume




# decision function
  decide <- function(bee, flower){
              # bee senses a flower, decides whether to visit based on its threshold
                # flowsp <- flower$species # may need this when we add more species
                currentTime <- max(flower$time)
              
                if(runif(1,0,100)<=8){
                  #make a small chance (8%) that bees visit flowers that are below their threshold
                  bee[bee$time==currentTime-1,"thresholdF1"]<-flower[flower$time==currentTime, "signal"]-1}
                if(bee[bee$time==currentTime-1,"thresholdF1"] < flower[flower$time==currentTime, "signal"]){
                  
                  return (visit(bee=bee,flower=flower))
 
                }
                else{
                  noChngBee <- bee[currentTime-1,]
                  noChngBee$time <- noChngBee$time+1
                  noChngBee$LearnRew<- NA
                  noChngBee$LearnSig<-NA
                  bee[currentTime,] <- noChngBee
                  
                  
                  flower$visit[currentTime] <- 0
                  bee$visit[currentTime] <- 0
                  
                 # no need to change flower because we already updated it...
                  return(list(bee=bee,flower=flower))
                           
                            }
  }
  
  # visit function 
bee <- BeeInit()
  visit <- function(bee,flower){
    
    #bee gets microbes, assess flower, updates threshold
      
      ## flowsp <- flower[["species"]]   # may need if we add more species of flowers
      currentTime <-  max(flower$time)
      prevBee <- bee[currentTime-1,]
      prevFlow <- flower[currentTime,]
      
      bee[currentTime,] <- bee[currentTime-1,]
      bee[currentTime,"time"]  <- currentTime
      bee[currentTime,"visit"] <- 1
      flower[currentTime,"visit"] <- 1
      
      
      bee[currentTime,"LearnRew"]    <- prevFlow$sugar
      bee[currentTime,"LearnVol"]    <- prevFlow$necVol
      bee[currentTime,"LearnValue"]  <- assess(bee=prevBee,flower=prevFlow)
      bee[currentTime,"LearnSig"]    <- prevFlow$signal
      #transfer microbes 
      
      bee[currentTime,"microbe1"] <-  .98 * prevBee$microbe1 + .2 * prevFlow$microbesP1 # amount bee will get, notice they eat about 92% of microbes so not all are avail for dispersal
      flower[currentTime,"microbesP1"] <- .8 * prevFlow$microbesP1 + .02 * prevBee$microbe1  # amount bee will deposit notice that these decrease by the amount transfered in other steps (e.g. 1-val here)
      flower[currentTime,"necVol"] <- .2 * prevFlow$necVol # pollinator takes 80% of the nectar
      
      #Bees assess flower (todo add volume) and update threshold
      
      #check out https://www.sciencedirect.com/science/article/pii/S0896627319308402 for bayesian decision model
      #TODO make a function for how threshold moves i.e. bayesian updating
       
      #step 1 - look back at the bees experience. 
      #right now we will take up to the 20 most recent visits, but this could change
      #either based on literature or w/e

      experience <- bee[bee$visit==1|bee$time==1,] |> tail(20)
      
   
      # see https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1439-0310.2006.01174.x?casa_token=o_kiaZSxTLQAAAAA:ojoLikVyLYlcP81oB2JOf69ZH1pBiqheJSLgCEcKD2cZdL-qRTuXfnALXk5CE3EOo-sb4s3C741qCu4
      #step 2, bees build a linear model about the world based on what they know
      # and the information they just got. They weight previous experience higher if it was more 
      # recent (3X more), and their most recent experience is even higher weighted (2X highest previous)
        
        mod<-lm(data=experience, LearnSig~LearnValue, weights=c(seq(1,3,length.out=nrow(experience)-1),6))
        
      # then we predict what signal they should use as their threshold based on the model
        
       bee[currentTime,"thresholdF1"]<- predict(object = mod, 
                                              newdata = data.frame(LearnValue=bee$LearnValue[1]))
        
        
        
  #transfer pollen 
       bee[currentTime,"pollenLoadF1"] <- prevBee$pollenLoadF1 + prevFlow$pollenDepRate - .3 * prevBee$pollenLoadF1
       flower[currentTime,"pollenF1"] <- prevFlow$pollenF1     + .3 * prevBee$pollenLoadF1
      
      
      return(list(bee=bee,flower=flower))      
       
  }
  

  
#
  # time step 

  daytimestep <- function(bee = NULL, flower=NULL, microbe=NULL){
    # microbes grow 
    
   
    flower1 <- grow(flower = flower, microbe = microbe)
    
    # bees forage
    
    if(is.null(bee) == F){
    return(decide(bee = bee, flower = flower1))
    }
    
    if(is.null(bee)==T){
      return(flower1)
    }
    ###NOT FOR NOW 
    # flowers get closer to death
    # new flowers open
  }
  
  

  
  nighttimestep <- function(bee, flower){
    #microbes grow 
    #flowers get closer to death
    #new flowers open
  }
  
  
  ####end functions

  
#Test the functions ####

  testFlower <- FlowerInit( signalPar = c(30.1,100))
  testBee <- BeeInit(microbesPar = c(0,200) , thresholdPar = c(0,30))
  testMicrobe <- MicrobeInit(growthRate = 1.002)
  
  for (i in 1:100){
  daytimestep(bee = testBee,flower = testFlower, microbe = testMicrobe) -> tmp
  testBee <- tmp[["bee"]]
  testFlower <- tmp[["flower"]]
  }
  
  plot <- ggplot(testBee[1:100,], aes(x=time,y=thresholdF1))
  plot+geom_point()
  
  plot <- ggplot(testBee[1:100,], aes(x=time,y=visit))
  plot+geom_point()
  
  plot <- ggplot(testBee[2:100,], aes(x=time,y=LearnRew))
  plot+geom_point()
  
  plot <- ggplot(testBee[2:100,], aes(x=time,y=LearnSig))
  plot+geom_point()
  
  
  plot <- ggplot(testBee, aes(x=time,y=pollenLoadF1))
  plot+geom_point()
  
  plot <- ggplot(testFlower, aes(x=time,y=pollenF1))
  plot+geom_point()
  
  plot <- ggplot(testFlower, aes(x=time,y=microbesP1))
  plot+geom_line() +geom_line(aes(x=time,y=necVol))
  
  plot <- ggplot(testFlower, aes(x=time,y=microbeK))
  plot+geom_line()
  
  
  # so for this to update we use the decide function and then assign the lists to each party
  
#lets see if a bee learns####
 
 
 
  
  
  #heres how to do it with FOREACH, but it doesnt really handle generating random numbers in the callwell...####
  # didshelearn <- foreach(z = 1:40,
  #                        beta = betas,
  #                        .packages=c("ggplot2","ggpubr"),
  #                        .combine='rbind') %dopar% {
  #     Learner <- list(BeeInit(microbesPar = c(100,100))) #bee with 100 microbes
  # 
  #     Field <- list()
  # 
  #     for(i in 1:100){
  #         Field[[i]] <- FlowerInit(signalPar = c(0,100),sugarPar = c(30,10))
  #         Field[[i]]$sugar <-abs( Field[[i]]$signal*beta +rnorm(1,0,8))
  #      }
  # 
  #     ti<- do.call(what = rbind,Field)
  #     ggplot(ti,aes(x=signal,y=sugar))+geom_point()
  # 
  #     PassiveMicrobe <- MicrobeInit(sugarDelt = 0,signalDelt = 0)
  # 
  # 
  # 
  #     for(i in 1:100){
  #           #first figure out which flower bees will go to by randomly sampling as many flower as there are bees
  #           # print(i) not needed for parallel
  #           encounters <- sample(1:length(Field), size = length(Learner), replace=F)
  #           encounterList <- vector("list",length(Field))
  #           encounterList[encounters] <- Learner
  # 
  #           test <- mapply(FUN=function(x,y){
  #             tmp <- daytimestep(flower=x,bee=y,microbe=PassiveMicrobe)
  #           },Field,encounterList,SIMPLIFY=F)
  # 
  #           Learner <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
  #           Field <- lapply(X=test, FUN = function(x){
  #             if(!is.null(x$flower)){return(x$flower)}
  #             else(return(x))})
  #     }
  # 
  # 
  #     Field <- lapply(X=Field, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
  #       return(x)})
  # 
  # 
  #     # ggplot(Learner[[1]],aes(x=time,y=visit))+geom_point()
  #     # ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()
  # 
  #     field<-do.call(rbind,Field)
  # 
  #     # ggplot(field,aes(x=time,y=necVol,group=ID)) + geom_line()
  #     # ggplot(field,aes(x=time,y=microbesP1,group=ID)) + geom_line()
  #     # ggplot(field,aes(x=time,y=microbeK,group=ID)) + geom_line()
  #     # ggplot(field,aes(x=time,y=signal,group=ID)) + geom_line()
  #     # ggplot(field,aes(x=time,y=cumulativeVisits,group=ID,color=signal)) + geom_line()
  #     # ggplot(Learner[[1]],aes(x=time,y=LearnValue,size=LearnSig))+geom_point()
  #     #
  # 
  # 
  #     # didshelearn[z,]<-c(NA,NA)
  #     # didshelearn$learnedT[z] <- mean(tail(na.omit(Learner[[1]]$thresholdF1),30))
  #     # didshelearn$realT[z] <- predict(lm(data=ti, signal~sugar),newdata = data.frame(sugar=30))
  #     pdf(paste("Plots/",z,"plot.pdf",sep=""),w=6,h=3)
  #     print(
  #     ggarrange(
  #         ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
  #         ggplot(ti,aes(y=signal,x=sugar))+geom_point()+geom_smooth(method="lm")+scale_y_continuous(limits=c(0,max(ti$signal)))+geom_vline(xintercept=30)
  #       )
  #     )
  #     dev.off()
  #     return(
  #       data.frame(Threshold = mean(tail(na.omit(Learner[[1]]$thresholdF1),30)),
  #                  TrueVal=predict(lm(data=ti, signal~sugar),newdata = data.frame(sugar=30))))
  # 
  # }
  
  
  #instead lets do it with mcapply####
  
  #step 1 designate a function 
  
  simulation <- function(nbees) { # we could name some parameters we wanted to study in this line and then specify them in the initation calls!
                           
    #make a field of random flowers
                           Field <- list()
                           beta<-abs(rnorm(1,0,.4))
                           signals<-NA
                           for(i in 1:100){
                             Field[[i]] <- FlowerInit(signalPar = c(0,100),sugarPar = c(30,10))
                             Field[[i]]$signal <-abs( Field[[i]]$sugar*beta +rnorm(1,0,8)+30)
                             signals[i] <- Field[[i]]$signal
                           }
                           Learner<-list()
                           for(i in 1:nbees){
                             Learner[i] <- list(BeeInit(microbesPar = c(100,100),thresholdPar = c(median(signals)-rnorm(1,10,2),median(signals)+rnorm(1,10,2)))) #bee with 100 microbes
                           }
                           ti<- do.call(what = rbind,Field)
                           ggplot(ti,aes(x=signal,y=sugar))+geom_point()+geom_smooth(method="lm")
                           
                           PassiveMicrobe <- MicrobeInit(sugarDelt = 0,signalDelt = 0,growthRate = 40)
                           
                           
                           
                           for(i in 1:100){
                             #first figure out which flower bees will go to by randomly sampling as many flower as there are bees
                             # print(i) not needed for parallel
                             encounters <- sample(1:length(Field), size = length(Learner), replace=F)
                             encounterList <- vector("list",length(Field))
                             encounterList[encounters] <- Learner
                             
                             test <- mapply(FUN=function(x,y){
                               tmp <- daytimestep(flower=x,bee=y,microbe=PassiveMicrobe)
                             },Field,encounterList,SIMPLIFY=F)
                             
                             Learner <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
                             Field <- lapply(X=test, FUN = function(x){
                               if(!is.null(x$flower)){return(x$flower)}
                               else(return(x))})
                           }
                           
                           
                           Field <- lapply(X=Field, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
                           return(x)})
                           
                           
                           ggplot(Learner[[1]],aes(x=time,y=visit))+geom_point()
                           ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()

                           field<-do.call(rbind,Field)
                           
                            ggplot(field,aes(x=time,y=necVol,group=ID)) + geom_line()
                            ggplot(field,aes(x=time,y=microbesP1,group=ID)) + geom_line()
                            ggplot(field,aes(x=time,y=microbeK,group=ID)) + geom_line()
                            ggplot(field,aes(x=time,y=signal,group=ID)) + geom_line()
                            ggplot(field,aes(x=time,y=cumulativeVisits,group=ID,color=signal)) + geom_line()
                            ggplot(Learner[[1]],aes(x=time,y=LearnValue,size=LearnSig))+geom_point()
                           
                           
                           pdf(paste("Plots/",nbees,"plot.pdf",sep=""),w=6,h=3)
                           print(
                             ggarrange(ncol=2,nrow = 2,
                               ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
                               ggplot(ti,aes(y=signal,x=sugar))+geom_point()+geom_smooth(method="lm")+scale_y_continuous(limits=c(0,max(ti$signal)))+geom_vline(xintercept=30),
                               ggplot(Learner[[1]],aes(x=LearnValue,y=LearnSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
                               ggplot(Learner[[1]],aes(x=LearnRew,y=LearnSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal)))                             )
                           )
                           dev.off()
                           return(list(
                             bees=Learner,
                             flowers=Field,
                             learndata= data.frame(Threshold = mean(tail(na.omit(Learner[[1]]$thresholdF1),30)),
                                        TrueVal=predict(lm(data=ti, signal~sugar),newdata = data.frame(sugar=30))))
                           
                           )
                           
                         }
  
  beePopParms <- rep(seq(from=1,to=5),10)
  simulations <- mclapply(X=beePopParms,FUN=function(x){simulation(nbees =x)}, mc.cores = 11)
  
  
  
  bees<-data.frame(simulations[[1]]$bees)
  bees$beePop <- 1
  bees$simnum <- 1
  for(i in 2:length(beePopParms)){
    df <- do.call(rbind,simulations[[i]]$bees)
    df$beePop <- beePopParms[i]
    df$simnum <-i
    bees<-rbind(bees,df)
  }
  ggplot(bees,aes(x=time,y=thresholdF1,group=ID,color=beePop))+
    geom_line(alpha=.2)+
    facet_grid(rows="beePop")
  
  
  
  flowers<-do.call(rbind,simulations[[1]]$flowers)
  flowers$beePop <- 1
  flowers$simnum <-1
  for(i in 2:length(beePopParms)){
    df <- do.call(rbind,simulations[[i]]$flowers)
    df$beePop <- beePopParms[i]
    df$simnum<-i
    flowers<-rbind(flowers,df)
  }
  
  ggplot(flowers,aes(x=time,y=cumulativeVisits,group=ID,color=beePop))+
    geom_line()+
    facet_grid(rows="beePop")
  
  ggplot(flowers,aes(x=time,y=microbeP1,group=ID,color=beePop))+
    geom_line()+
    facet_grid(rows="beePop")
  
  ggplot(flowers,aes(x=time,y=necVol,group=ID,color=beePop))+
    geom_line()+
    facet_grid(rows="beePop")
  
  finalstate <- flowers[flowers$time==max(flowers$time),]
  finalstate$beePop <- as.factor(finalstate$beePop)
  pct<-function(x){length(x[x$cumulativeVisits==0,])/length(x)}
  ggplot(finalstate,aes(x=as.factor(beePop),y=cumulativeVisits)) +geom_point(position=position_jitter(.3))
  pct<-function(x){nrow(x[x$cumulativeVisits==0,])/nrow(x)}
  
  finalSum <- summarize
        
  
  ggplot(finalSum,aes(x=beePop,y=percentVisit, color=simnum)) +geom_point(size=3,position=position_jitter(width=.2))
  
  # lets look at some microbe parameters...####
  #we are going to use a starting point of 1 pollinator per hundred plants
  #because that gave us something like a 30% visitation rate
  # this matches yeast colonization in the field : https://academic.oup.com/aob/article/103/9/1415/145505
  
  simulation2 <- function(nbees,signalDelt,sugarDelt) { # we could name some parameters we wanted to study in this line and then specify them in the initation calls!
    
    #make a field of random flowers
    Field <- list()
    beta<-abs(rnorm(1,0,.4))
    signals<-NA
    for(i in 1:100){
      Field[[i]] <- FlowerInit(signalPar = c(0,100),sugarPar = c(30,10))
      signals[i] <- Field[[i]]$signal
    }
    
    Learner<-list()
    for(i in 1:nbees){
      Learner[i] <- list(BeeInit(microbesPar = c(100,100),thresholdPar = c(median(signals)-rnorm(1,10,2),median(signals)+rnorm(1,10,2)))) #bee with 100 microbes
    }
    ti<- do.call(what = rbind,Field)
    ggplot(ti,aes(x=signal,y=sugar))+geom_point()+geom_smooth(method="lm")
    
    PassiveMicrobe <- MicrobeInit(sugarDelt = sugarDelt,signalDelt = signalDelt ,growthRate = 20)
    
    
    
    for(i in 1:100){
      #first figure out which flower bees will go to by randomly sampling as many flower as there are bees
      # print(i) not needed for parallel
      encounters <- sample(1:length(Field), size = length(Learner), replace=F)
      encounterList <- vector("list",length(Field))
      encounterList[encounters] <- Learner
      
      test <- mapply(FUN=function(x,y){
        tmp <- daytimestep(flower=x,bee=y,microbe=PassiveMicrobe)
      },Field,encounterList,SIMPLIFY=F)
      
      Learner <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
      Field <- lapply(X=test, FUN = function(x){
        if(!is.null(x$flower)){return(x$flower)}
        else(return(x))})
    }
    
    
    Field <- lapply(X=Field, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
    return(x)})
    
    
    ggplot(Learner[[1]],aes(x=time,y=visit))+geom_point()
    ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()
    
    field<-do.call(rbind,Field)
    
    ggplot(field,aes(x=time,y=necVol,group=ID)) + geom_line()
    ggplot(field,aes(x=time,y=microbesP1,group=ID)) + geom_line()
    ggplot(field,aes(x=time,y=microbeK,group=ID)) + geom_line()
    ggplot(field,aes(x=time,y=signal,group=ID)) + geom_line()
    ggplot(field,aes(x=time,y=cumulativeVisits,group=ID,color=signal)) + geom_line()
    ggplot(Learner[[1]],aes(x=time,y=LearnValue,size=LearnSig))+geom_point()
    
    
    pdf(paste("Plots/",nbees,"plot.pdf",sep=""),w=6,h=3)
    print(
      ggarrange(ncol=2,nrow = 2,
                ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
                ggplot(ti,aes(y=signal,x=sugar))+geom_point()+geom_smooth(method="lm")+scale_y_continuous(limits=c(0,max(ti$signal)))+geom_vline(xintercept=30),
                ggplot(Learner[[1]],aes(x=LearnValue,y=LearnSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
                ggplot(Learner[[1]],aes(x=LearnRew,y=LearnSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal)))                             )
    )
    dev.off()
    return(list(
      bees=Learner,
      flowers=Field,
      learndata= data.frame(Threshold = mean(tail(na.omit(Learner[[1]]$thresholdF1),30)),
                            TrueVal=predict(lm(data=ti, signal~sugar),newdata = data.frame(sugar=30))))
      
    )
    
  }
  
  
  simulations <- mclapply(X=c(1,10,20,30,40,50),FUN=function(x){simulation(nbees =x)}, mc.cores = 11)
  
  
  
  
  simout <-do.call(rbind,simulations)
  ggplot(simout, aes(x=TrueVal,y=Threshold))+geom_point()+geom_abline(intercept=0,slope=1)
  
  
  

  
    #Sensitivity Analyses ####
  
  SugDparm <- seq(from=-.01, to =.01, length=3)
  SigDparm <- seq(from=-.01, to=.01,  length =3)             
  
  Parms <- data.frame(expand.grid(SugDparm,SigDparm))
  ParmsList <- split(Parms,seq(nrow(Parms)))

  
  
  
  onerun <-   function(sigDel, sugDel){
                
              # make a community of pollinators none of which have microbes! 
              #This is a list of bees. Each bee is a data frame
                nPol    <- 4   #number of indiv per species
                Pollinators <- list()
                for(i in 1:nPol){
                  Pollinators[[i]] <- BeeInit(thresholdPar = c(0,1),
                                              microbesPar = c(0,200),
                                              preferencesPar = c(10,100)
                  )
                }  
              
                # make a community of sterile flowers. this is a list of flowers and each is a data frame
                nFlow     <- 20
                Flowers <- list()
                for(i in 1:nFlow){
                  Flowers[[i]] <- FlowerInit(signalPar = c(0,100))
                }
                
                # make a microbe with the parameters provided. we only have a single microbe right now
                ourMicrobe <- MicrobeInit(growthRate = 1.05,signalDelt = sigDel, sugarDelt = sugDel)
                
                for(i in 1:100){
                  #first figure out which flowers bees will go to by randomly sampling as many flowers as there are bees
                  encounters <- sample(1:length(Flowers), size = length(Pollinators), replace=F)
                  encounterList <- vector("list",length(Flowers))
                  encounterList[encounters] <- Pollinators
                  
                  test <- mapply(FUN=function(x,y){
                    tmp <- daytimestep(flower=x,bee=y,microbe=ourMicrobe)
                    },Flowers,encounterList,SIMPLIFY=F)
                  
                  Pollinators <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
                  Flowers <- lapply(X=test, FUN = function(x){
                    if(!is.null(x$flower)){return(x$flower)} 
                      else(return(x))})
                  
                }  
                
                #calculate cumulative visits!
                Flowers <- lapply(X=Flowers, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
                                                              x$signaldelt <- sigDel
                                                              x$sugardelt <- sugDel
                                                                    return(x)})
                
               Pollinators <- lapply(X=Pollinators, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
                x$signaldelt <- sigDel
                x$sugardelt <- sugDel
                return(x)})
                
                updatedFlowers<- do.call(rbind,Flowers)
                updatedBees <- do.call(rbind,Pollinators)
                return(list(flowers=updatedFlowers,bees=updatedBees))
         }               
   onerun(.1,.1)
  
  #try parallel processing
  detectCores()
  Sys.time()
  test<- mclapply(X = ParmsList, 
                  FUN = function(x){onerun(sigDel = x[1,1],sugDel = x[1,2])},
                  mc.cores = 11)
  Sys.time()
  
  flr<-test[[1]]$flowers
  for(i in 2:length(test)){flr<-rbind(flr,test[[i]]$flowers)}
  
  bbz<-test[[1]]$bees
  for(i in 2:length(test)){bbz<-rbind(bbz,test[[i]]$bees)}
  dev.new()      
  
  flr$ID <- as.factor(flr$ID)
  
  timetofirstvisit <-flr[flr$cumulativeVisits==1,]
  timetofirstvisit <- ddply(.data = timetofirstvisit, .variables = .(ID),  function(x){x[x$time==min(x$time),]})
  signalbytimetovis <- ggplot(timetofirstvisit, aes(x=signaldelt,y=time)) 
  signalbytimetovis + geom_point(alpha=.2,position=position_jitter(width=.00002))
  sugarbytimetovis <- ggplot(timetofirstvisit, aes(x=sugardelt,y=time)) 
  sugarbytimetovis + geom_point(alpha=.2,position=position_jitter(width=.00002))
  
  summary(timetofirstvisit)
  
  firstVisSig <- glm(data=timetofirstvisit, time~signaldelt*sugardelt, family="poisson")
  summary(firstVisSig)
  
  summaryTtoV <- ddply(.data = timetofirstvisit, .variables = .(sugardelt,signaldelt), summarize, time=mean(time) )
  ggplot(summaryTtoV, aes(x=signaldelt,y=time,color=(sugardelt),group=as.factor(sugardelt))) + geom_line()
  
  
  
  
  names(Parms) <- c("signaldelt","sugardelt")
  Parms$TtoFirst <- predict(firstVisSig,newdata = Parms,se.fit = T, type="response")$fit
  Parms$TtoFirstSE <- predict(firstVisSig,newdata = Parms,se.fit = T,type="response")$se.fit
  ggplot(Parms, aes(x=signaldelt,y=sugardelt,color=(TtoFirst),fill=TtoFirst)) +geom_tile()
  
  
  
  ggplot(Parms[Parms$sugardelt==-0.010,], aes(x=signaldelt,y=TtoFirst)) + geom_point() +
    geom_errorbar(aes(ymin=TtoFirst-TtoFirstSE, ymax=TtoFirst+TtoFirstSE))
  
  ggplot(Parms, aes(color=as.factor(sugardelt),x=signaldelt,y=TtoFirst)) + geom_point() 
  
  maxvisits <- flr[flr$time==max(flr$time),]
  
  maxvisitsplot <- ggplot(maxvisits, aes(x=signaldelt,y=cumulativeVisits, color=signaldelt))
  maxvisitsplot + geom_point(alpha=.2,position=position_jitter(width=.0002))
  
  summary(maxvisits)
  table(maxvisits$cumulativeVisits)
  signal <- glm(data=maxvisits, cumulativeVisits~signaldelt*sugardelt, family = "poisson")
  summary(signal)    
  
  ggplot(maxvisits, aes(x=signaldelt,y=cumulativeVisits, color=as.factor(sugardelt)))+
    geom_point()+
    facet_grid(vars(sugardelt))
  
  
  Parms$Max <- predict(signal,newdata = Parms[,c("signaldelt","sugardelt")], type="response")
  ggplot(Parms, aes(x=signaldelt,y=sugardelt,color=(Max),fill=Max)) +geom_tile()
  
  
  pollenplot <- ggplot(maxvisits, aes(x=signaldelt,y=signaldelt, color=pollenF1))
  pollenplot + geom_point(alpha=.2,position=position_jitter(width=.0002))
  
  
  
  pollenmod <- lm(data=maxvisits, pollenF1~as.factor(signaldelt)*as.factor(sugardelt))
  summary(pollenmod)
  
  plotdata <- data.frame(emmeans(pollenmod, specs=c("signaldelt","sugardelt")))
  ggplot(plotdata,aes(x=as.factor(signaldelt), color=as.factor(sugardelt),y=emmean)) +
    geom_point(position=position_dodge(.4)) +
    geom_line(aes(group=as.factor(sugardelt)),position=position_dodge(.4)) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position=position_dodge(.4))+
    scale_y_continuous("Predicted pollen receipt")
  
  
  
  
  Parms$Pollen <- predict(pollenmod,newdata = Parms[,c("signaldelt","sugardelt")], type="response")
  ggplot(Parms, aes(x=signaldelt,y=Pollen,color=as.factor(sugardelt),fill=Pollen)) + geom_point()+geom_line(aes(group=as.factor(sugardelt)))
  
  
  
   
  