
####  Multiple stressor YRC Bayes hierarchical SR model ####

library(tidyverse)
library(postpack)
library(R2jags)
library(mcmcplots)
library(bayesboot)


#Note that a lot of this code is copied/adapted from https://github.com/curryc2/Cook_Inlet_Chinook


#Define Workflow Paths ====================================================
wd <- file.path(getwd())
dir.figs <- file.path(wd,"Plots")
dir.output <- file.path(wd,"Output")

# Read in and process spawner recruit data ------------------------------------------------


## Question 1: 
# Still unsure if this is the correct data. Using the RS_df2 dataframe from Make Figures.R
# but it samples a subset of the posterior (e.g. Line 140 from "Make Figures.R") and I'm wondering why this step was done
# and if I should instead calculate exact medians by removing the subsampling line?
# BC: Yes it is the right data, but you are correct inference was only based on 1000 samples, and the indexing was wonky. A new version 
#     of brood table with proper indexing and based on full set of posteriors is now read in as "brood_table"  

#Question 2:
#Indexing doesn't match online Supplement B. Have adjusted here to match but this needs checking
# BC: new file has correct indexing so can comment all of this out
#RS_df2 <- read.csv("RS_df2.csv")

brood_table <- read.csv("brood_table.csv")

#remove dummy data for 1981

#RS_df2_int1 <- filter(RS_df2,BroodYear>1981)

#add three years to match indexing. This creates dataframe from 1985 to 2019

#RS_df2_int2 <- mutate(RS_df2_int1,BroodYear=(BroodYear+3))

#remove years without recruitment estimates i.e. 2016 and higher

#data <- filter(RS_df2_int2,R_med!="NA")

#Question 3
#is it okay to use years 2013 - 2015 even though they were not estimated based on full recruitment?
# BC: The model that generated the estimates derives complete recruitment estimates through 2015 based on 
#      incomplete observations of brood year recruitment for brood years 2013-2015. So you could use these, but for 2013 onward they 
#      are based on increasingly less data. So I would only used brood years 1985-2012, and have done so for now.  

brood_table <- filter(brood_table,BroodYear<2013)

#remove dummy data for 1981

#RS_df2_int1 <- filter(RS_df2,BroodYear>1981)

#sort into matrices

pops = sort(unique(brood_table$population))
n.pops = length(pops)
n.years = vector(length=n.pops)
years <- matrix(nrow=n.pops,ncol=28)

p <- 1
for(p in 1:n.pops) {
  n.years[p] <- length(unique(brood_table$BroodYear[brood_table$population==pops[p]]))
  years[p,1:n.years[p]] <- sort(unique(brood_table$BroodYear[brood_table$population==pops[p]]))
}#next i

# Spawners and Recruits
Spawners <- matrix(nrow=n.pops,ncol=max(n.years))
Recruits <- matrix(nrow=n.pops,ncol=max(n.years))
ln.Recruits <- matrix(nrow=n.pops,ncol=max(n.years))

p <- 1
for(p in 1:n.pops) {
  y <- 1
  for(y in 1:n.years[p]) {
    year <- years[p,y]
    
    #Spawners
    Spawners[p,y] <- brood_table$S_med[brood_table$population==pops[p] & brood_table$BroodYear==year]
    #Recruits
    Recruits[p,y] <- brood_table$R_med[brood_table$population==pops[p] & brood_table$BroodYear==year]
    ln.Recruits[p,y] <- log(brood_table$R_med[brood_table$population==pops[p] & brood_table$BroodYear==year])
  }#next y
}#next p


#popn order is Carmacks, Lower Main, Middle Main, Pelly, Stewart, Teslin, Upper Lakes, White-Donjek




# Add covariate data ------------------------------------------------------

#A few general covariates for starter model below


##### Ice break up date at Dawson in outmigration year (t = +2) ----------
#data from https://yukonriverbreakup.com/statistics

Ice_int1 <- read.csv("Dawson Ice break up.csv")

#index to outmigration year

Ice_int1 <- Ice_int1 %>% 
  mutate(Year2 = Year-2) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2015
Ice_int2 <- filter(Ice_int1,Year>1984&Year<2013)

#standardize
Ice_int3 <- t(scale(Ice_int2$Julian_day))

#copy rows down
Ice_out <- Ice_int3[rep(seq_len(nrow(Ice_int3)),each=8),]



####### Migration temperature using water temp data from Emmonak (t = 0)##########################################################################

#data from https://www.adfg.alaska.gov/CF_R3/external/sites/aykdbms_website/DataSelection.aspx
#search for Lower Yukon Test Fishing data

Emmonak_water_int1 <- read.csv("Emmonak Water Temp Data.csv")

#sampling locations across year change, look at summary

Emmonak_water_summary1  <- Emmonak_water_int1 %>% 
 group_by(Project.Year,Location) %>% count(Location)

#two years missing middle mouth and one missing big eddy

#ggplot(Emmonak_water_int1_summary,aes(fill=Location,y=n,x=Project.Year))+
  #geom_bar(stat="identity",position="dodge")+scale_y_continuous(trans="log2")

ggplot(Emmonak_water_summary1,aes(fill=Location,y=n,x=Project.Year))+
  geom_bar(stat="identity",position="dodge")+scale_y_continuous(trans="log2")

#sample timing differs too. before 2004 temperatures were only taken at 8 am and 8 pm
#how many have 8 am and 8 pm?

# BC: I get an error running the code below, think it has somehting to do with lubridate::as_date(), 
#     not going to bother trying to figure out, and will instead just generate a dummy covariate @ end of section
Emmonak_water_int2 <- Emmonak_water_int1 %>% 
  select(-Time) %>% 
  rename(Date_Time=Date) %>% 
  separate(Date_Time, into = c("Date", "Time"), sep = " ", remove = FALSE) %>%
  mutate(Date = lubridate::as_date(Date, format = "%Y-%m-%d"),
         Time = hms::as_hms(str_c(Time, ":00")))

Emmonak_water_summary2  <- Emmonak_water_int2 %>% 
  group_by(Project.Year,Location) %>% count(Time)


Emmonak_water_summary3 <- filter(Emmonak_water_summary2,Project.Year<2004)

ggplot(Emmonak_water_summary3,aes(fill=Location,y=n,x=Time))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#missing a few PM samples and only one AM sample in 1995. 
#Use AM for now and sub in PM for 1995 until this can be adjusted

Emmonak_water_int2 <- as.character(Emmonak_water_int2$Time)

Emmonak_water_int2$Time <- format(as.POSIXct(Emmonak_water_int2$Time), "%H:%M:%S")
Emmonak_water_int3 <- filter(Emmonak_water_int2,Time=="08:00:00")

ggplot(Emmonak_water_int3,aes(x=Project.Year,y=Temperature..avg.))+
  geom_point()

Emmonak_water_summary4  <- Emmonak_water_int3 %>% 
  group_by(Project.Year,Location) %>% count(Time)

ggplot(Emmonak_water_summary4,aes(fill=Location,y=n,x=Location))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#best bet for now is to use Big Eddy data at 8 AM. Will be missing for three years:
#1990 was only middle mouth, 1995 not sampled at 8 AM, 2003 also only sampled at MM
#For now sub in MM for BE in 1990 and 2003 and use 8 PM for 1995 from Big Eddy

#extract Big Eddy
Emmonak_water_int4 <- filter(Emmonak_water_int3,Location=="Big Eddy")

#extract MM for 1990 and 2003
Emmonak_water_int5 <- filter(Emmonak_water_int3,Location=="Middle Mouth"&Project.Year=="1990"|
                               Location=="Middle Mouth"&Project.Year=="2003")

#extract 8 PM for Big Eddy in 1995
Emmonak_water_int6 <- filter(Emmonak_water_int2,Time=="20:00:00"&Location=="Big Eddy"&Project.Year=="1995")

#bind dataframes together
Emmonak_water_int7 <- bind_rows(Emmonak_water_int4,Emmonak_water_int5,Emmonak_water_int6)

#add julian date
Emmonak_water_int8 <- mutate(Emmonak_water_int7,julian_date=yday(Date))

ggplot(Emmonak_water_int8,aes(y=Temperature..avg.,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

#figure out general timing window for passing Emmonak for Can-origin 
#travel time is approx one month between river mouth and the border
#estimate late May to early Aug but pick window based on comparability of data to start

#looks like data missing before July 15 in 2000
#1999 started a bit late - June 15
#can maybe pull data from Pilot for this year if needed?
#check Jones paper - does the model just infer absences? or is better to estimate inputs

Emmonak_water_int9 <- filter(Emmonak_water_int8,Temperature..avg.!="NA")

Emmonak_water_summary5  <- Emmonak_water_int9 %>% 
  group_by(Project.Year) %>% summarise(
    count = n(),
    min = min(julian_date),
    max=max(julian_date))

#most years have data from June 9 - Aug 1. use this window for now

Emmonak_water_int10 <- filter(Emmonak_water_int9,julian_date>159&julian_date<214)

ggplot(Emmonak_water_int10,aes(y=Temperature..avg.,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

Emmonak_water_summary6  <- Emmonak_water_int10 %>% 
  group_by(Project.Year) %>% summarise(
    count = n(),
    min = min(julian_date),
    max=max(julian_date))


#calculate mean daily water temp during migration as temporary covariate

Emmonak_water_int11 <- Emmonak_water_int10 %>% 
  group_by(Project.Year) %>% summarise(water_temp=mean(Temperature..avg.))

#remove extra years i.e. <1985 and >2015
Emmonak_water_int12 <- filter(Emmonak_water_int11,Project.Year>1984&Project.Year<2016)

ggplot(Emmonak_water_int12,aes(fill=water_temp,y=water_temp,x=Project.Year))+
  geom_bar(stat="identity")

Migration_temp_int1 <- select(Emmonak_water_int12,water_temp)

#standardize
Migration_temp_int2 <- t(scale(Migration_temp_int1$water_temp))

#copy rows down
Migration_temp <- Migration_temp_int2[rep(seq_len(nrow(Migration_temp_int2)),each=8),]

#To do; check for missing data within each year and determine if needs estimation
#To do: check for obvious outliers/ measurement errors
#TO DO: derive weekly maximums?
#TO DO: estimate run timing past Emmonak and adjust migration temp by population and maybe even by year
#To do: confirm that Emmonak temps are outside of delta influence
#To do: figure out if there is a way to use migration temps for returning years also - this will matter if returns are too stressed to hit the border count

#other possible variables to include:
#Mean summer air temperature (June - Aug) in basin for rearing juveniles (t = +1)
#Mean summer precipitation (June - Aug) in basin for rearing juveniles (t = +1)
#Mean precipitation during spawning and early incubation (Aug - Nov) in basin (t = 0)
#Temperature for spawning/incubation in basin (t = 0)
#Sea surface temp (t = +2, +3)
#climate drivers (t = +2, +3)
#hatchery fish (t = +2, +3)



# Merge covariate data ----------------------------------------------------

# Define covariate names
names.covars <- c("Ice_out","Migration_temp")
n.covars <- length(names.covars)

covars <- array(data=NA,dim=c(n.pops, max(n.years), n.covars))

# BC: creating dummy data for migration temp covariate
Migration_temp <- matrix(data=rnorm(length(Ice_out),0,0.25),
                         nrow = 8,
                         ncol = 28)

##put covariates into covars array for model

library(abind)

covars <- abind(Ice_out,Migration_temp,along=3)
print(covars)

#TO DO: remember to test for correlations between covariates



# Make a list for Jags ----------------------------------------------------

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
)


# Write model code --------------------------------------------------------


jags_model = function(){
  
  #Priors

  #hyperparameters

#question 4: 
#should I be adding a hyperparameter for alpha?
# BC: To start I would say this is not needed as this would tend to reduce among populaiton estimates of intrinsic productivity (b/c it imposes
#     assumption of alphas being drawn from common distribution). If you had alot of variability in amount of data among populaitons then might 
  #   make more sense.  
  
    for(c in 1:n.covars) {
      mu.coef[c] ~ dnorm(0, pow(25,-2))
      sigma.coef[c] ~ dnorm(0, pow(5,-2));T(1e-3,100)
      dist.coef[c] ~ dnorm(mu.coef[c], pow(sigma.coef[c],-2))
    }#next c
  
  #population-specific
  
  for(p in 1:n.pops) {
    exp.alpha[p] ~ dunif(0,25)
    alpha[p] <- log(exp.alpha[p])
    beta[p] ~ dnorm(0,pow(0.1,-2));T(0,100)
    sigma.oe[p] ~ dnorm(0, pow(1,-2));T(1e-3,100)
    
    #Covariate Effects
    for(c in 1:n.covars) {
      coef[p,c] ~ dnorm(mu.coef[c],pow(sigma.coef[c],-2))
    }
  }#next p
  
# question 5:

#should I test with and without autocorrelation term?
# BC: Yes that would be a good idea.   

  #possible priors for autocorrelation if needed:
  #Estimated to be shared among all populations.
  #phi ~ dunif(-0.99, 0.99)
  #log.resid.0 <- 0  # multi-stock model fixes this at zero
  
  #PREDICTIONS
  for(p in 1:n.pops) {
    for(y in 1:n.years[p]) {
      
      for(c in 1:n.covars) {
        cov.eff[p,y,c] <- coef[p,c]*covars[p,y,c]
      }
      
      pred.rec[p,y] <- Sobs[p,y]*exp(alpha[p] - Sobs[p,y]*beta[p] + sum(cov.eff[p,y,1:n.covars]))
      
        }#next y
    
    #Fill in Ragged Components- pretty sure I don't need this step because years for each popn are the same?
    # BC: correct. 
    # for(y in (n.years[p]+1):max(n.years)) {
    #  pred.rec[p,y] <- 0
    #}
  }#next p
  
  #LIKELIHOODS
  for(p in 1:n.pops) {
    for(y in 1:n.years[p]) {
      lnRobs[p,y] ~ dnorm(log(pred.rec[p,y]), pow(sigma.oe[p],-2))
      
    }#next y
  }#next p
  
}
  
#### Specify initial values ####

jags_inits = function() {
  exp.alpha <- runif(n.pops,1,2)
  beta <- runif(n.pops,0,0.001)
  obs.sigma <- runif(n.pops,0.1,1)
  
#question 6: 
#code copied from Curry might be wrong: should this be simga.oe instead of obs.sigma?
# BC: Yeah I think so, but probably does not matter if you are not initializing the 
#    sampler with values drawn form obs.sigma (uniform distribution between 0.1-1) b/c
#    I think JAGS auto initializes if you do not specify something. 

  mu.coef <- rnorm(n.covars, 0, 1)
  sigma.coef <- rgamma(n.covars, 1, 1)
  
  Return <- list(exp.alpha=exp.alpha, beta=beta, obs.sigma=obs.sigma,
                 mu.coef=mu.coef, sigma.coef=sigma.coef)
  return(Return)
}


#### Set nodes to monitor ####
# BC: There were some errors here with naming of parameters b/c some had an underscore in them but in model all use a period instead

#jags_params = c("alpha","exp_alpha", "beta", "sigma.oe","mu_coef","sigma_coef",
#                "coef","cov_eff","dist_coef","pred.rec")

jags_params = c("alpha","exp.alpha", "beta", "sigma.oe","mu.coef","sigma.coef",
                "coef","cov.eff","dist.coef","pred.rec")

##### Specify MCMC Dimensions #####
# BC: I reduced # of iterations to make it quicker to assess model but have not looked at convergence diagnostics etc.
jags_dims = c(
  ni = 1e4,  # number of post-burn-in samples per chain
  nb = 1e4,  # number of burn-in samples
  nt = 200,     # thinning rate
  nc = 3      # number of chains
)


##### Run the model with Jags #####

#try this one
#started 12:06
out <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=200, 
  n.iter=2e4, 
  n.burnin=1e4)  

#Save
saveRDS(out, file=file.path(dir.output,"out.rds"))

#Write Output File for Diagnostics
write.csv(out$BUGSoutput$summary, file=file.path(dir.figs,"out_summary.csv"))

#Visualize results

#Hyper means for covariates
par(mfcol=c(1,2), mar=c(5,0,1,0), oma=c(1,1,3,1))
c <- 1
for(c in 1:n.covars) {
  plotPost(out$BUGSoutput$sims.list$mu.coef[,c], showCurve=TRUE, main='', xlab=names.covars[c],
           xlim=c(-0.5,0.5), rope=0)
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0,alpha=0.5))
}


#Population specific covariate "effects"
par(mfcol=c(1,2), mar=c(2,5,3,1), oma=c(2,2,1,1))
c <- 1
for(c in 1:n.covars) {
  caterplot(out$BUGSoutput$sims.list$coef[,,c],
            labels=pops, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), style='plain', col='blue')
  mtext(names.covars[c], side=3, outer=FALSE, line=1)
  caterpoints(apply(out$BUGSoutput$sims.list$coef[,,c],2,median), reorder=FALSE, pch=21, col='red', bg='orange')
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
}
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5)
mtext('Population', side=2, outer=TRUE, font=2, line=0.5)
