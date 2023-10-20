##########################################################################
#
# RED FOOTED BOOBY POPULATION MODEL FOR ALDABRA
#
##########################################################################

# population model adapted from Cubaynes et al. 2010: https://royalsocietypublishing.org/doi/full/10.1098/rsbl.2010.0778
# implemented in JAGS based on Kery and Schaub 2012
# written by steffen.oppel@vogelwarte.ch in November 2023

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



### Calculation of stable age distribution 
### CREATING THE POPULATION MATRIX ###
# 
# seabird.matrix<-matrix(c(
#   0,0,0.519*0.5*0.56,
#   0.85,0,0,
#   0,0.92,0.92),ncol=3, byrow=T)
# stable.stage(seabird.matrix)



setwd("C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/PeripheralProjects/RFBO")
sink("RFBO_popmod.jags")
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
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    mean.fec ~ dunif(0.3,0.7)         ## uninformative prior for breeding success

    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    mean.ad.surv ~ dunif(0.9, 0.94)             # Prior for mean survival
    mean.juv.surv ~ dunif(0.80,0.85)    ## based on juvenile survival for Balearic shearwaters in the Med.
    breed.prop ~ dunif(0.5,1)
    
    # -------------------------------------------------        
    # 1.3. Priors FOR POPULATION COUNT ERROR
    # -------------------------------------------------
    sigma.obs ~ dunif(0,10)  #Prior for SD of observation process (variation in detectability)
    tau.obs<-pow(sigma.obs,-2)

    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------

      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1
      JUV[1]<-round(Nad.breed[1]*0.5*(mean.fec)*breed.prop)
      N1[1]<-round(Nad.breed[1]*0.5*(mean.fec)*breed.prop*mean.juv.surv)
      N2[1]<-round(Nad.breed[1]*0.5*(mean.fec)*breed.prop*mean.juv.surv*mean.ad.surv)
      Nad.breed[1] ~ dunif(2000,2500)         # initial value of population size
      Nad.nonbreed[1] <- Nad.breed[1] * (1-breed.prop)         # initial value of non-breeder size
      
      for (tt in 2:n.years){

        ## THE PRE-BREEDERS ##
        JUV[tt] ~ dbin(mean.fec,round(0.5 * Nad.breed[tt]*breed.prop))                                   ### number of locally produced FEMALE chicks
        N1[tt]  ~ dbin(mean.juv.surv, max(2,round(JUV[tt-1])))                                               ### number of 1-year old survivors 
        N2[tt] ~ dbin(mean.ad.surv, max(2,round(N1[tt-1])))                                                      ### number of 2-year old survivors

        ## THE BREEDERS ##
        Nad.surv[tt] ~ dbin(mean.ad.surv, max(2,round(N2[tt-1]+Nad.breed[tt-1]+Nad.nonbreed[tt-1])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        Nad.breed[tt] ~ dbin(breed.prop, max(2,Nad.surv[tt]))                            ### the annual number of breeding birds is the proportion of adult survivors
        Nad.nonbreed[tt] <- max(2,Nad.surv[tt]-Nad.breed[tt])                            ### the annual number of nonbreeding birds is the remainder of adult survivors
      } # tt
      

    # -------------------------------------------------        
    # 2.2. Likelihood for fecundity: binomial regression from the number of surveyed colonies
    # -------------------------------------------------
    for (n in 1:(n.col)){      ### N of colonies where breeding success was measured
      J[n] ~ dbin(mean.fec,BP[n])
    } #	close loop over every site
    
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for population.counts
    # -------------------------------------------------
      ## Observation process
    
      for (t in 1:n.years){
        Nad.count[t] ~ dnorm(Nad.breed[t], tau.obs)# Distribution for random error in observed numbers (counts)
      }# run this loop over t= nyears
    
    
    # -------------------------------------------------        
    # 4. DERIVED PARAMETERS
    # -------------------------------------------------

    ## DERIVED POPULATION GROWTH RATE 
      for (tt in 2:n.years){
        lambda[tt]<-Nad.breed[tt]/max(1,Nad.breed[tt-1])
        loglam[tt]<-log(lambda[tt])
      } ## end of tt
      
      #growth.rate <- exp((1/(n.years-1))*sum(loglam[1:(n.years-1)]))  ### geometric mean growth rate
      
    
  }  ## END MODEL LOOP
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# 4. SET UP AND RUN INTEGRATED POPULATION MODEL
#########################################################################

# Bundle data
jags.data <- list(Nad.count=countdata$RFBO,
                  n.years=length(countdata$RFBO),
                  J=as.integer(productivity),
                  BP=rep(100,length(productivity)),
                  n.col=length(productivity))

# Initial values 
inits <- function(){list(mean.ad.surv = runif(1, 0.9, 0.94),
                         mean.fec=runif(1,0.3,0.7))}  ### adjusted for REV1 as frequency of good years


# Parameters monitored
parameters <- c("sigma.obs","breed.prop","mean.fec","mean.juv.surv","mean.ad.surv","Nad.breed")

# MCMC settings
ni <- 1500
nt <- 10
nb <- 500
nc <- 3

# Call JAGS from R (model created below)
RFBO_IPM <- jags(jags.data, inits, model.file="C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/PeripheralProjects/RFBO/RFBO_popmod.jags",  ## changed from v4 to v6 on 10 Aug
                     parameters.to.save=parameters,
                 n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T,n.iter = ni)


### save model workspace
#save.image("RFBO_IPM_REV1.RData")
load("RFBO_IPM_REV1.RData")







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

TABLE1<-out %>% filter(parameter %in% c('mean.fec[1]','mean.fec[2]','orig.fec','prop.good','fec.decrease','mean.ad.surv','mean.juv.surv','growth.rate[1]','growth.rate[2]')) %>%
  select(parameter,c(5,3,7))

names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("original proportion of good breeding year (1956)","productivity (poor year)","productivity (good year)","annual decline in frequency of good breeding year","current proportion of good breeding year (2014-2019)","first year survival probability","annual adult survival probability","annual population growth rate (no eradication)","annual population growth rate (with eradication)")
TABLE1
#fwrite(TABLE1,"RFBO_demographic_parameter_estimates_REV1.csv")

## FORMAT TABLE FOR MANUSCRIPT

TABLE1<-TABLE1 %>% mutate(MED=paste(round(Median,3)," (",round(lowerCL,3)," - ", round(upperCL,3),")", sep="")) %>%
  select(Parameter,MED) %>%
  rename(`Median (95% credible interval)`=MED)
setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RFBO_pop_model")
#fwrite(TABLE1,"TABLE1.csv")

## REPORT QUANTITIES FOR RESULTS SECTION
sum(succ$R)
dim(CH)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 1: POPULATION TRAJECTORY UNDER BOTH SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
RFBOpop<-out[(grep("Nad.breed\\[",out$parameter)),c(12,5,4,6)] %>%
  mutate(Year=seq(1969,2023)) %>%
  rename(parm=parameter,median=`50%`,lcl=`25%`,ucl=`75%`) %>%
  dplyr::select(parm,Year,median,lcl,ucl)


### CREATE PLOT FOR BASELINE TRAJECTORY

ggplot()+
  geom_line(data=RFBOpop, aes(x=Year, y=median), linewidth=1)+
  geom_ribbon(data=RFBOpop,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  
  ## format axis ticks
  scale_y_continuous(name="Red-footed Booby pairs", limits=c(0,25000),breaks=seq(0,25000,5000))+
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
ggsave("Figure1.tif", width = 9, height = 6, device='tiff', dpi=300)
ggsave("RFBO_population_projection_REV1_CI50.jpg", width=9, height=6)


### CREATE INSET PLOT FOR FUTURE PROJECTION
RFBOpop$ucl[RFBOpop$ucl>1000000]<-999999

ggplot()+
  geom_line(data=RFBOpop[RFBOpop$Year>2019,], aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=RFBOpop[RFBOpop$Year>2019,],aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  
  ## format axis ticks
  scale_y_continuous(name="MacGillivray's Prion pairs (millions)", limits=c(0,1000000),breaks=seq(0,1000000,100000),labels=seq(0,1,0.1))+
  scale_x_continuous(name="Year", limits=c(2020,2057), breaks=seq(2021,2056,5), labels=as.character(seq(2021,2056,5)))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.position = c(0.2,0.85),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
ggsave("RFBO_population_projection_REV1_CI50_INSET.jpg", width=9, height=6)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 2: EXTINCTION PROBABILITY OVER TIME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### EXTRACT AND COMBINE DATA FROM ALL CHAINS 
selcol<-grep("Nad.breed",dimnames(RFBO_IPM$samples[[1]])[[2]])   ## FIND COLUMS WE NEED
allchainsamples <- data.frame()
for(chain in 1:3) {
  samplesout<-as.data.frame(RFBO_IPM$samples[[chain]][,selcol]) %>% gather(key="parm", value="value")
  allchainsamples <- rbind(allchainsamples,as.data.frame(samplesout))
}

  
### CALCULATE EXTINCTION PROBABILITY
extprop <- allchainsamples %>%
    mutate(scen.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
    mutate(Scenario=ifelse(scen.index==1,"no eradication","with eradication")) %>%
    mutate(Year=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,1])+1955) %>%
    
    mutate(n=1, inc=ifelse(value<500,1,0)) %>%   ### DEFINE EXTINCTION PROBABILITY HERE
    group_by(Scenario,Year) %>%
    summarise(ext.prob=sum(inc)/sum(n)) %>%
    filter(Year>2019)
  
head(allchainsamples)
head(extprop)
dim(extprop)


## CREATE A COLOUR PALETTE FOR THE NUMBER OF CHICKS RELEASED
colfunc <- colorRampPalette(c("firebrick","cornflowerblue"))


ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=Scenario), size=1)+
  
  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,0.35),breaks=seq(0,0.35,0.05), labels=as.character(seq(0,35,5)))+
  scale_x_continuous(name="Year", breaks=seq(2020,2055,5), labels=as.character(seq(2020,2055,5)))+
  guides(color=guide_legend(title="Scenario"),fill=guide_legend(title="Scenario"))+
  scale_colour_manual(palette=colfunc)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black",angle=45, vjust = 1, hjust=1),
        axis.title=element_text(size=18),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        legend.position = c(0.15,0.90),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
ggsave("Figure2.tif", width = 9, height = 6, device='tiff', dpi=300)
ggsave("RFBO_extinction_probability_REV1_250.jpg", width=9, height=6)



### CREATE TABLE 2 FOR MANUSCRIPT ###

head(extprop)
TABLE2<- extprop %>%
  filter(Year==2057) %>%
  select(ext.prob,Scenario) %>%
  mutate(ext.prob=ext.prob*100) %>%
  mutate(ext.prob=ifelse(ext.prob<1,"< 1%",paste0(ext.prob,"%")))
fwrite(TABLE2,"TABLE2.csv")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROPORTION OF SIMULATIONS WITH NEGATIVE GROWTH RATE AFTER MOUSE ERADICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### EXTRACT AND COMBINE DATA FROM ALL CHAINS 
futlamsamples <- data.frame()
for(chain in 1:3) {
  samplesout<-as.data.frame(RFBO_IPM$samples[[chain]][,9]) %>% gather(key="parm", value="value")
  futlamsamples <- rbind(futlamsamples,as.data.frame(samplesout))
}
head(futlamsamples)

futlamsamples %>% mutate(decline=ifelse(value<0.95,1,0)) %>%
  tabyl(decline)




