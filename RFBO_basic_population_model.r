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
countdata[Year==1969,2]<-2277
countdata[Year==2000,2]<-4095
countdata[Year==2023,2]<-22607

## Breeding success (%) at five monitoring colonies in the 2022 NW season
procuctivity<-c(55.3,61.8,43.6,34.3,47.6)



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



setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR")
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
    
    mean.fec ~ dunif(0,1)         ## uninformative prior for BAD YEARS

    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    mean.phi ~ dunif(0, 1)             # Prior for mean survival
    mean.juv.surv ~ dunif(0.63,0.78)    ## based on juvenile survival for Balearic shearwaters in the Med.

    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    for (scen in 1:2){

      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
      fec.proj[scen,1]<-mean.fec[year.prop.good[scen,1]+1]    ## takes good or bad year fecundity
      year.prop.good[scen,1] ~ dbern(orig.fec)
      
      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
      JUV[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1]))
      N1[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1])*mean.juv.surv)
      N2[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1])*mean.juv.surv*mean.phi)
      N3[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1])*mean.juv.surv*mean.phi*mean.phi)
      Ntot.breed[1,scen] ~ dunif(2000000,5000000)         # initial value of population size
      
      for (tt in 2:65){

        ## LINEARLY DECREASING PROBABILITY OF A GOOD YEAR FROM 1956 to 2014
        year.fec.prop[scen,tt]<- max(0,min(1,(orig.fec + fec.decrease*tt))) ## calculate yearly proportion of good breeding year, but constrain to 0-1 to avoid invalid parent value
        year.prop.good[scen,tt] ~ dbern(year.fec.prop[scen,tt])
        fec.proj[scen,tt]<-mean.fec[year.prop.good[scen,tt]+1]    ## takes good or bad year fecundity
        breed.prop[scen,tt] ~ dunif(0.85,0.95)    ## breeding propensity
        
        ## THE PRE-BREEDERS ##
        JUV[tt,scen] ~ dbin(fec.proj[scen,tt],round(0.5 * Ntot.breed[tt,scen]*breed.prop[scen,tt]))                                   ### number of locally produced FEMALE chicks
        N1[tt,scen]  ~ dbin(mean.juv.surv, max(2,round(JUV[tt-1,scen])))                                               ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(mean.phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(mean.phi, max(2,round(N2[tt-1,scen])))                                                      ### number of 3-year old survivors
        
        ## THE BREEDERS ##
        Ntot.breed[tt,scen] ~ dbin(mean.phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        
      } # tt
      
      for (tt in 66:PROJ){
        
        ## SELECT GOOD OR BAD OR RODENT FREE FECUNDITY FOR FUTURE
        year.fec.prop[scen,tt]<- min(1,max(0,(orig.fec + fec.decrease*tt))) ## calculate yearly proportion of good breeding year, but constrain to 0-1 to avoid invalid parent value
        year.prop.good[scen,tt] ~ dbern(year.fec.prop[scen,tt])
        fec.proj[scen,tt]<-max(mean.fec[year.prop.good[scen,tt]+1],(scen-1)*full.fec)    ## takes current fecundity for scenario 1 and full fecundity for scenario 2
        breed.prop[scen,tt] ~ dunif(0.85,0.95)    ## breeding propensity
        
        ## THE PRE-BREEDERS ##
        JUV[tt,scen] ~ dbin(fec.proj[scen,tt],round(0.5 * Ntot.breed[tt,scen]*breed.prop[scen,tt]))                                   ### need a discrete number otherwise dbin will fail, dpois must be >0
        N1[tt,scen]  ~ dbin(mean.juv.surv, max(2,round(JUV[tt-1,scen])))                                               ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(mean.phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(mean.phi, max(2,round(N2[tt-1,scen])))                                                      ### number of 3-year old survivors
        
        ## THE BREEDERS ##
        Ntot.breed[tt,scen] ~ dbin(mean.phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        
      } # tt
      
    } # scen    
    


    
    # -------------------------------------------------        
    # 2.2. Likelihood for fecundity: Poisson regression from the number of surveyed broods
    # -------------------------------------------------
    for (t in 1:(T.fec)){      ### T-1 or not
      J[t] ~ dpois(rho.fec[t])
      rho.fec[t] <- R[t]*mean.fec[goodyear[t]+1]
      goodyear[t] ~ dbern(prop.good)
    } #	close loop over every year in which we have fecundity data
    
    
    
    
    # -------------------------------------------------        
    # 4. DERIVED PARAMETERS
    # -------------------------------------------------

    ## DERIVED POPULATION GROWTH RATE 
    for (scen in 1:2){
      for (tt in 1:33){
        lambda[tt,scen]<-Ntot.breed[tt+67,scen]/max(1,Ntot.breed[tt+66,scen])
        loglam[tt,scen]<-log(lambda[tt,scen])
      } ## end of tt
      
      growth.rate[scen] <- exp((1/(33))*sum(loglam[1:(33),scen]))  ### geometric mean growth rate
      
    } ## end of scen
    
  }  ## END MODEL LOOP
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# 4. SET UP AND RUN INTEGRATED POPULATION MODEL
#########################################################################

# Bundle data
jags.data <- list(## survival
  y = CH,
  f = f,
  n.occasions = dim(CH)[2],
  nind = dim(CH)[1],
  #mean.juv.surv.prop= 0.728/0.894,  ## juvenile survival based on proportion of adult survival from Jiguet 2007
  
  ## fecundity
  R =succ$R,
  J=succ$J,
  T.fec=length(succ$J),
  goodyear=c(0,0,1,0,0,0),
  
  ## population process
  Ntot.obs=matrix(c(3500000,800000,3500000,800000), ncol=2), ### adjusted for v6 to 2-5 million
  PROJ=66+36)  ### adjusted for v6 to 103 years

# Initial values 
inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1),
                         orig.fec= runif(1, 0.75, 0.85))}  ### adjusted for REV1 as frequency of good years


# Parameters monitored
#parameters <- c("orig.fec","mean.fec","fec.decrease","fec.drop","mean.juv.surv","mean.phi","growth.rate","lambda","Ntot.breed")
parameters <- c("orig.fec","mean.fec","fec.decrease","prop.good","mean.juv.surv","mean.phi","growth.rate","lambda","Ntot.breed")

# MCMC settings
ni <- 150000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (model created below)
MAPR_IPM <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_IPM_REV1_varfec.jags",  ## changed from v4 to v6 on 10 Aug
                     n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T,n.iter = ni)


### save model workspace
setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_press\\MAPR_pop_model")
#save.image("MAPR_IPM_REV1.RData")
load("MAPR_IPM_REV1.RData")







#########################################################################
# 5. SUMMARISE OUTPUT AND PLOT POPULATION TRAJECTORY
#########################################################################
## compile output
out<-as.data.frame(MAPR_IPM$summary)
out$parameter<-row.names(MAPR_IPM$summary)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(out)

TABLE1<-out %>% filter(parameter %in% c('mean.fec[1]','mean.fec[2]','orig.fec','prop.good','fec.decrease','mean.phi','mean.juv.surv','growth.rate[1]','growth.rate[2]')) %>%
  select(parameter,c(5,3,7))

names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("original proportion of good breeding year (1956)","productivity (poor year)","productivity (good year)","annual decline in frequency of good breeding year","current proportion of good breeding year (2014-2019)","first year survival probability","annual adult survival probability","annual population growth rate (no eradication)","annual population growth rate (with eradication)")
TABLE1
#fwrite(TABLE1,"MAPR_demographic_parameter_estimates_REV1.csv")

## FORMAT TABLE FOR MANUSCRIPT

TABLE1<-TABLE1 %>% mutate(MED=paste(round(Median,3)," (",round(lowerCL,3)," - ", round(upperCL,3),")", sep="")) %>%
  select(Parameter,MED) %>%
  rename(`Median (95% credible interval)`=MED)
setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\MAPR_pop_model")
#fwrite(TABLE1,"TABLE1.csv")

## REPORT QUANTITIES FOR RESULTS SECTION
sum(succ$R)
dim(CH)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 1: POPULATION TRAJECTORY UNDER BOTH SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
MAPRpop<-out[(grep("Ntot.breed\\[",out$parameter)),c(12,5,4,6)] %>%
  mutate(Year=rep(seq(1956,2057),2)) %>%
  mutate(scenario=as.numeric(str_extract_all(parameter,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
  mutate(Scenario=ifelse(scenario==1,"no eradication","with eradication")) %>%
  filter(!(Scenario=="with eradication" & Year<2024)) %>%
  #filter((Scenario=="with eradication")) %>%
  #rename(parm=parameter,median=`50%`,lcl=`2.5%`,ucl=`97.5%`) %>%
  rename(parm=parameter,median=`50%`,lcl=`25%`,ucl=`75%`) %>%
  dplyr::select(parm,Scenario,Year,median,lcl,ucl)


### summary for manuscript
MAPRpop %>% filter(Year==2020)
MAPRpop %>% filter(Year==1956)
174684/3502083

### CREATE PLOT FOR BASELINE TRAJECTORY
MAPRpop$ucl[MAPRpop$ucl>5000000]<-4999999

ggplot()+
  geom_line(data=MAPRpop, aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=MAPRpop,aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  
  ## format axis ticks
  scale_y_continuous(name="MacGillivray's Prion pairs (millions)", limits=c(0,5000000),breaks=seq(0,5000000,500000),labels=seq(0,5,0.5))+
  scale_x_continuous(name="Year", limits=c(1956,2057), breaks=seq(1956,2056,20), labels=as.character(seq(1956,2056,20)))+
  
  ## add count data
  geom_segment(aes(x=1956, xend=1956,y=0.4*5000000,yend=0.5*10000000),lineend = "round", size=2, colour="darkblue") +
  geom_segment(aes(x=2000, xend=2000,y=0.4*1500000,yend=0.5*2000000),lineend = "round", size=2, colour="darkblue") +
  
  
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
ggsave("MAPR_population_projection_REV1_CI50.jpg", width=9, height=6)


### CREATE INSET PLOT FOR FUTURE PROJECTION
MAPRpop$ucl[MAPRpop$ucl>1000000]<-999999

ggplot()+
  geom_line(data=MAPRpop[MAPRpop$Year>2019,], aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=MAPRpop[MAPRpop$Year>2019,],aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  
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
ggsave("MAPR_population_projection_REV1_CI50_INSET.jpg", width=9, height=6)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 2: EXTINCTION PROBABILITY OVER TIME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### EXTRACT AND COMBINE DATA FROM ALL CHAINS 
selcol<-grep("Ntot.breed",dimnames(MAPR_IPM$samples[[1]])[[2]])   ## FIND COLUMS WE NEED
allchainsamples <- data.frame()
for(chain in 1:3) {
  samplesout<-as.data.frame(MAPR_IPM$samples[[chain]][,selcol]) %>% gather(key="parm", value="value")
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
ggsave("MAPR_extinction_probability_REV1_250.jpg", width=9, height=6)



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
  samplesout<-as.data.frame(MAPR_IPM$samples[[chain]][,9]) %>% gather(key="parm", value="value")
  futlamsamples <- rbind(futlamsamples,as.data.frame(samplesout))
}
head(futlamsamples)

futlamsamples %>% mutate(decline=ifelse(value<0.95,1,0)) %>%
  tabyl(decline)






##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## SUMMARISE WIND AND RAIN FOR MAPR BREED SUCCESS YEARS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## count number of precipitation-free days

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Weather"), silent=T)
weather<-fread("Gough_weather_data_2000_2020.csv")
head(weather)
weathersummary<-weather %>% group_by(Year) %>%
  summarise(mean_wind=mean(Windspeed, na.rm=T),max_wind=max(Windspeed, na.rm=T),rain=sum(Rain)) %>%
  filter(Year %in% succ$Year)
weathersummary

### McClelland et al. 2018 used 'rain free days' but that metric is suspicious in 2014!
rainsummary<-weather %>% group_by(Year, Month, Day) %>%
  summarise(rain=sum(Rain)) %>%
  mutate(rainfree=ifelse(rain==0,1,0)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(rainfreedays=sum(rainfree)) %>%
  filter(Year %in% succ$Year)
rainsummary



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## PRODUCE WEATHER CHART
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

ggplot(weathersummary, mapping = aes(x = Year, y = mean_wind, group = 1)) + 
  geom_bar(mapping = aes(y = rain/100), stat = "identity", color="cornflowerblue", fill="cornflowerblue", width = 0.5) + 
  geom_line(color="indianred", size=1.5) +
  

  ## format axis ticks
  scale_y_continuous("Mean wind speed (m/s)", sec.axis = sec_axis(~.*100, name = "Precipitation (mm)")) +
  scale_x_continuous(name="Year", limits=c(2013.6,2019.4), breaks=seq(2014,2019,1)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=18))
