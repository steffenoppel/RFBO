##########################################################################
#
# RED FOOTED BOOBY POPULATION MODEL FOR ALDABRA
#
##########################################################################

# population model adapted from Cubaynes et al. 2010: https://royalsocietypublishing.org/doi/full/10.1098/rsbl.2010.0778
# implemented in JAGS based on Kery and Schaub 2012
# written by steffen.oppel@vogelwarte.ch in November 2023

## NEED TO DO: NARROW DOWN PRIORS FOR PRODUCTIVITY AND BREEDING PROPENSITY
## model struggles to converge because too many things can change at once
## reverted to version of 29 Oct that seemed to produce sensible output

## dreaded 'slicer stuck at value with infinite density' may occur because there is no variation: https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/1b1ea77d/

### REVISED APPROACH 9 NOV 2023: FIXED SURVIVAL and prop.breeding to reduce ambiguity

library(tidyverse)
library(jagsUI)
library(data.table)
library(lubridate)
library(popbio)
library(janitor)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# 1. ENTER THE DATA  DATA FOR ALDABRA
#########################################################################
## sent by Michelle Risi in October 2023 (WhatsApp)
countdata<-data.frame(Year=seq(1969,2023,1),RFBO=NA)
countdata[countdata$Year==1969,2]<-2277
countdata[countdata$Year==2000,2]<-4095
countdata[countdata$Year==2023,2]<-22607
countdata$phase<-1
countdata$phase[countdata$Year>2000]<-2

## Breeding success (%) at five monitoring colonies in the 2022 NW season
productivity<-c(55.3,61.8,43.6,34.3,47.6)



#########################################################################
# 2. Specify BASIC POPULATION MODEL WITH TWO SCENARIOS
#########################################################################

### DEMOGRAPHIC PARAMETERS 

#Juvenile survival: 	0.85+0.0527 	from Cubaynes et al. 2010
#Adult survival: 	0.92+0.0028 	from Cubaynes et al. 2010
#Age at maturity: 	3 	from	various web sources
#breeding propensity: 0.56+0.0361 adapted from Cubaynes et al. 2010

hist(rbeta(1000,92,8))
hist(rbeta(1000,85,17))
hist(rbeta(1000,52,45))
hist(rbeta(1000,50,85))
hist(rbeta(1000,95,10))


### Calculation of stable age distribution 
### CREATING THE POPULATION MATRIX ###
# 
# seabird.matrix<-matrix(c(
#   0,0,0.519*0.5*0.56,
#   0.85,0,0,
#   0,0.92,0.92),ncol=3, byrow=T)
# stable.stage(seabird.matrix)



setwd("C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/PeripheralProjects/RFBO")
sink("RFBO_popmod_2phase_free_imm.jags")
cat("
  
  
  model {
    #-------------------------------------------------
    # - population model for the Red-Footed Booby population on Aldabra
    # - age structured model with 4 age classes 
    # - all survival based on literature
    # - productivity based on data from 2023 only
    # - goal is to examine whether the increase could have occurred from local productivity or not
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    #mean.prop.good ~ dunif(0.1,0.9)        ## proportion of years that is good or bad (to allow past variation when good years were more common)
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    mean.fec[1] ~ dbeta(50,85) T(0.01,0.99)        ## uninformative prior for breeding success
    mean.fec[2] ~ dbeta(52,45) T(0.01,0.99)        ## uninformative prior for breeding success
    #mean.fec[1] ~ dunif(0.1,0.9)         ## uninformative prior for breeding success
    #mean.fec[2] ~ dunif(0.3,0.9)         ## uninformative prior for breeding success

    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    mean.ad.surv[1] ~ dbeta(92, 8) T(0.01,0.99)            # Prior for mean survival
    mean.ad.surv[2] ~ dbeta(93, 8) T(0.01,0.99)            # Prior for mean survival
    #mean.ad.surv[1] <- 0.92             # Prior for mean survival
    #mean.ad.surv[2] <- 0.92             # Prior for mean survival
    mean.juv.surv[1] ~ dbeta(85,17) T(0.01,0.99)    ## 
    mean.juv.surv[2] ~ dbeta(87,17) T(0.01,0.99)   ##
    #mean.juv.surv[2] <- 0.85  ## avoid different rates for juveniles
    #mean.juv.surv[1] <- 0.85  ## avoid different rates for juveniles
    breed.prop[1] ~ dbeta(90,10) T(0.01,0.99)
    breed.prop[2] ~ dbeta(95,10) T(0.01,0.99)
    #breed.prop[1] ~ dunif(0.5,1)
    #breed.prop[2] ~ dunif(0.5,1)
    #breed.prop[1] <- 0.9
    #breed.prop[2] <- 0.9
    mean.imm[1]<-0
    mean.imm[2]<-50

    
    # -------------------------------------------------        
    # 1.3. Priors FOR POPULATION COUNT ERROR
    # -------------------------------------------------
    # sigma.obs ~ dunif(0,1)  #Prior for SD of observation process (variation in detectability)
    # tau.obs<-pow(sigma.obs,-2)

    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------

      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1
      JUV[1]<-round(Nad.breed[1]*0.5*(mean.fec[1])*breed.prop[1])
      N1[1]<-round(Nad.breed[1]*0.5*(mean.fec[1])*breed.prop[1]*mean.juv.surv[phase[1]])
      N2[1]<-round(Nad.breed[1]*0.5*(mean.fec[1])*breed.prop[1]*mean.juv.surv[phase[1]]*mean.ad.surv[1])
      Nad.breed[1] ~ dunif(2000,2500)         # initial value of population size
      Nad.nonbreed[1] <- Nad.breed[1] * (1-breed.prop[1])         # initial value of non-breeder size
      #prop.good[1] ~ dbern(mean.prop.good)
      ann.imm[1]<-0
      
      for (tt in 2:n.years){
      
        ## RANDOMLY DRAW ANNUAL IMMIGRANTS
        ann.imm[tt] ~ dpois(mean.imm[phase[tt]])
        

        ## THE PRE-BREEDERS ##
        JUV[tt] ~ dbin(mean.fec[phase[tt]],round(0.5 * Nad.breed[tt]*breed.prop[phase[tt]]))                                   ### number of locally produced FEMALE chicks
        N1[tt]  ~ dbin(mean.juv.surv[phase[tt]], max(2,round(JUV[tt-1])))                                               ### number of 1-year old survivors 
        N2[tt] ~ dbin(mean.ad.surv[phase[tt]], max(2,round(N1[tt-1])))                                                      ### number of 2-year old survivors

        ## THE BREEDERS ##
        Nad.surv[tt] ~ dbin(mean.ad.surv[phase[tt]], max(2,round(N2[tt-1]+Nad.breed[tt-1]+Nad.nonbreed[tt-1])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        Nad.breed[tt] ~ dbin(breed.prop[phase[tt]], max(2,Nad.surv[tt]+ann.imm[tt]))                            ### the annual number of breeding birds is the proportion of adult survivors
        Nad.nonbreed[tt] <- max(2,Nad.surv[tt]-Nad.breed[tt])                            ### the annual number of nonbreeding birds is the remainder of adult survivors
      } # tt
      

    # -------------------------------------------------        
    # 2.2. Likelihood for fecundity: binomial regression from the number of surveyed colonies
    # -------------------------------------------------
    for (n in 1:(n.col)){      ### N of colonies where breeding success was measured
      J[n] ~ dbin(mean.fec[2],BP[n])
    } #	close loop over every site
    
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for population.counts
    # -------------------------------------------------
      ## Observation process
    
      # for (t in 1:n.years){
      #   Nad.count[t] ~ dnorm(Nad.breed[t], tau.obs)# Distribution for random error in observed numbers (counts)
      # }# run this loop over t nyears
    
    
    # -------------------------------------------------        
    # 4. DERIVED PARAMETERS
    # -------------------------------------------------

    # ## DERIVED POPULATION GROWTH RATE 
    #   for (tt in 2:n.years){
    #     lambda[tt]<-Nad.breed[tt]/max(1,Nad.breed[tt-1])
    #     loglam[tt]<-log(lambda[tt])
    #   } ## end of tt
    #   
    #   #growth.rate <- exp((1/(n.years-1))*sum(loglam[1:(n.years-1)]))  ### geometric mean growth rate
      
    
  }  ## END MODEL LOOP
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# 4. SET UP AND RUN INTEGRATED POPULATION MODEL
#########################################################################

### INTRODUCE RANDOM NUMBERS IN COUNT DATA TO AVOID SLICER STUCK ERROR DUE TO 0 VARIATION
#countdata$RFBO[17]<-as.integer(runif(1,countdata$RFBO[1],countdata$RFBO[32]))
#countdata$RFBO[46]<-as.integer(runif(1,countdata$RFBO[32],countdata$RFBO[55]))

# Bundle data
jags.data <- list(Nad.count=countdata$RFBO,
                  n.years=length(countdata$RFBO),
                  phase=countdata$phase,
                  J=as.integer(productivity),
                  BP=rep(100,length(productivity)),
                  n.col=length(productivity))

# Initial values 
inits <- function(){list(#Nad.breed=c(runif(1,2275,2280),rep(NA,length(countdata$RFBO)-1)),
                         mean.juv.surv = rbeta(2, 85, 17),
                         mean.ad.surv = rbeta(2, 92, 8),
                         mean.fec=rbeta(2,45,55))}  ### adjusted for REV1 as frequency of good years


# Parameters monitored
parameters <- c("mean.imm","breed.prop","mean.fec","mean.juv.surv","mean.ad.surv","Nad.breed")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 25000
nc <- 3

# Call JAGS from R (model created below)
RFBO_IPM <- jags(jags.data, inits, model.file="C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/PeripheralProjects/RFBO/RFBO_popmod_2phase_free_imm.jags",  ## changed from v4 to v6 on 10 Aug
                     parameters.to.save=parameters,
                 n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T,n.cores=nc,
                 #Rhat.limit=1.5,iter.increment=1000,max.iter=1000000) ### for autojags call
                 n.iter = ni)   ## for normal jags call







#########################################################################
# 5. SUMMARISE OUTPUT AND PLOT POPULATION TRAJECTORY
#########################################################################
## compile output
out<-as.data.frame(RFBO_IPM$summary)
out$parameter<-row.names(RFBO_IPM$summary)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(out)

TABLE1<-out %>% filter(parameter %in% c('mean.imm[1]','mean.imm[2]','mean.fec[1]','mean.fec[2]','breed.prop[1]','breed.prop[2]','mean.ad.surv[1]','mean.ad.surv[2]','mean.juv.surv[1]','mean.juv.surv[2]')) %>%
  select(parameter,c(5,3,7))

names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("immigrants","immigrants","proportion of breeders","proportion of breeders","productivity","productivity","first year survival probability","first year survival probability","annual adult survival probability","annual adult survival probability")
TABLE1$Period<-rep(c("1969-2000","2000-2022"), 5)

#fwrite(TABLE1,"RFBO_demographic_parameter_estimates_REV1.csv")

## FORMAT TABLE FOR MANUSCRIPT

TABLE1<-TABLE1 %>% mutate(MED=paste(round(Median,3)," (",round(lowerCL,3)," - ", round(upperCL,3),")", sep="")) %>%
  select(Parameter,MED, Period) %>%
  spread(key=Period, value=MED)
TABLE1
#fwrite(TABLE1,"TABLE1.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 1: POPULATION TRAJECTORY UNDER BOTH SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
RFBOpop<-out[(grep("Nad.breed\\[",out$parameter)),c(12,1,3,7)] %>%
  mutate(Year=seq(1969,2023)) %>%
  rename(parm=parameter,lcl=`2.5%`,ucl=`97.5%`) %>%
  dplyr::select(parm,Year,mean,lcl,ucl)


### CREATE PLOT FOR BASELINE TRAJECTORY

ggplot()+
  geom_line(data=RFBOpop, aes(x=Year, y=mean), linewidth=1)+
  geom_ribbon(data=RFBOpop,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_point(data=countdata, aes(x=Year, y=RFBO), size=3, colour="firebrick")+
  
  ## format axis ticks
  scale_y_continuous(name="Red-footed Booby pairs", limits=c(0,100000),breaks=seq(0,100000,10000))+
  scale_x_continuous(name="Year", limits=c(1969,2023), breaks=seq(1969,2023,5))+

  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.position = c(0.73,0.89),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
ggsave("RFBO_population_projection_with_immigration.jpg", width=9, height=6)






save.image("RFBO_popmod.RData")

### save model workspace
#save.image("RFBO_IPM_converged.RData")
load("RFBO_IPM_REV1.RData")



