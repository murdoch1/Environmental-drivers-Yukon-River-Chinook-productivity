
####  Multiple stressor YRC Bayes hierarchical SR model ####

library(tidyverse)
library(postpack)
library(R2jags)
library(mcmcplots)
library(bayesboot)
library(lubridate)


#Note that a lot of this code is copied/adapted from https://github.com/curryc2/Cook_Inlet_Chinook


#Define Workflow Paths ====================================================
wd <- file.path(getwd())
dir.figs <- file.path(wd,"Plots")
dir.output <- file.path(wd,"Output")
dir.data <- file.path(wd,"Data")

# Read in and process spawner recruit data ------------------------------------------------

brood_table <- read.csv(file.path(dir.data,"brood_table.csv"))

#removing years without full recruitment estimates
brood_table <- filter(brood_table,BroodYear<2013)

pops = sort(unique(brood_table$population))
n.pops = length(pops)
n.years = vector(length=n.pops)
years <- matrix(nrow=n.pops,ncol=28)

p <- 1
for(p in 1:n.pops) {
  n.years[p] <- length(unique(brood_table$BroodYear[brood_table$population==pops[p]]))
  years[p,1:n.years[p]] <- sort(unique(brood_table$BroodYear[brood_table$population==pops[p]]))
}



# Merge covariate data ----------------------------------------------------

Ice_out <- read.csv(file.path(dir.data,"/Environmental data/Processed/Ice_out.csv"))
Migration_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/Weekly_max_migration_temp_t0.csv"))
rearing_GDD5 <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_GDD.csv"))
#rearing_GDD5 <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_GDD_all.csv"))
rearing_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_prcp_max.csv"))
#rearing_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_prcp_max_all.csv"))
annual_snowpack <- read.csv(file.path(dir.data,"/Environmental data/Processed/snowpack_spawning_mean_Apr1.csv"))
#spawning_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_temp_air.csv"))
spawning_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_est_water_temp1.csv"))
spawning_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_prcp_max.csv"))
SST_summer <- read.csv(file.path(dir.data,"/Environmental data/Processed/SST_summer.csv"))
SST_winter <- read.csv(file.path(dir.data,"/Environmental data/Processed/SST_winter_new.csv"))
Pink_comp <- read.csv(file.path(dir.data,"/Environmental data/Processed/Pink.csv"))
Chum_comp <- read.csv(file.path(dir.data,"/Environmental data/Processed/Chum.csv"))
First_winter_FW <- read.csv(file.path(dir.data,"/Environmental data/Processed/Winter_temp.csv"))

# Define covariate names
names.covars <- c("Ice_out","Migration_temp","rearing_GDD5","rearing_prcp",
                  "annual_snowpack","spawning_temp","spawning_prcp","SST_summer",
                  "SST_winter","Pink_comp","Chum_comp","First_winter_FW")



n.covars <- length(names.covars)

covars <- array(data=NA,dim=c(n.pops, max(n.years), n.covars))



##put covariates into covars array for model

library(abind)

covars <- abind(Ice_out,Migration_temp,rearing_GDD5,rearing_prcp,
annual_snowpack,spawning_temp,spawning_prcp,SST_summer,
SST_winter,Pink_comp,Chum_comp,First_winter_FW,along=3)



print(covars)



# Write model code --------------------------------------------------------


jags_model = function(){
  
  #Priors
  
    for(c in 1:n.covars) {
      mu.coef[c] ~ dnorm(0, pow(25,-2))
      sigma.coef[c] ~ dnorm(0, pow(5,-2));T(1e-3,100) #original
      #sigma.coef[c] ~ dunif(0,100) #no difference
      dist.coef[c] ~ dnorm(mu.coef[c], pow(sigma.coef[c],-2))
    }#next c
  
  #population-specific
  
  for(p in 1:n.pops) {
    exp.alpha[p] ~ dunif(0,50)
    alpha[p] <- log(exp.alpha[p])
    #beta[p] ~ dnorm(0,pow(0.1,-2));T(0,100) original from Cook Inlet
     beta[p] ~ dnorm(0,pow(1000,-2));T(0,100) #beta from Fleischman et al to approximate uniform dist
    sigma.oe[p] ~ dnorm(0, pow(1,-2));T(1e-3,100)
    #sigma.oe[p] ~ dunif(0,100) #no difference
    
    #Covariate Effects
    for(c in 1:n.covars) {
      coef[p,c] ~ dnorm(mu.coef[c],pow(sigma.coef[c],-2))
    }
    # AR-1 Coeff
     phi[p] ~ dunif(-0.99, 0.99)
  }#next p
  

  
  #PREDICTIONS
  for(p in 1:n.pops) {
    # First Year
    for(c in 1:n.covars) {
        cov.eff[p,1,c] <- coef[p,c]*covars[p,1,c]
    }
    pred.rec.1[p,1] <- Sobs[p,1]*exp(alpha[p] - Sobs[p,1]*beta[p] + sum(cov.eff[p,1,1:n.covars]))
    
    log.pred.rec.1[p,1] <- log(pred.rec.1[p,1])
    log.pred.rec.2[p,1] <- log.pred.rec.1[p,1]
    log.resid[p,1] <- lnRobs[p,1] - log.pred.rec.1[p,1]
    
   #Subsequent Years
     for(y in 2:n.years[p]) {
      
       for(c in 1:n.covars) {
        cov.eff[p,y,c] <- coef[p,c]*covars[p,y,c]
      }
      
      pred.rec.1[p,y] <- Sobs[p,y]*exp(alpha[p] - Sobs[p,y]*beta[p] + sum(cov.eff[p,y,1:n.covars]))
      
     log.pred.rec.1[p,y] <- log(pred.rec.1[p,y])
     log.resid[p,y] <- lnRobs[p,y] - log.pred.rec.1[p,y]
     log.pred.rec.2[p,y] <- log.pred.rec.1[p,y]+log.resid[p,y-1]*phi[p]
     log.resid.2[p,y] <- lnRobs[p,y] - log.pred.rec.2[p,y]
      
        }#next y
    }#next p
  
  #LIKELIHOODS
  for(p in 1:n.pops) {
    
    #First Year
    lnRobs[p,1] ~ dnorm(log.pred.rec.2[p,1], pow(sigma.oe[p],-2))
    
    for(y in 2:n.years[p]) {
    lnRobs[p,y] ~ dnorm(log.pred.rec.2[p,y], pow(sigma.oe[p],-2))
      
    }#next y
  }#next p
}
# 


#### Specify initial values ####

jags_inits = function() {
  exp.alpha <- runif(n.pops,1,2)
  #beta <- runif(n.pops,0,0.001)
  beta <- runif(n.pops,0,0.3) #try increasing this based on results from popn diversity paper?
  sigma.oe <- runif(n.pops,0.1,1) #seems good based on prior data
  mu.coef <- rnorm(n.covars, 0, 1)
  sigma.coef <- rgamma(n.covars, 1, 1)
  
  #Return <- list(exp.alpha=exp.alpha, beta=beta, sigma.oe=sigma.oe,
                # mu.coef=mu.coef, sigma.coef=sigma.coef)
  Return <- list(exp.alpha=exp.alpha, beta=beta, sigma.oe=sigma.oe)
  return(Return)
}


#### Set nodes to monitor ####

jags_params = c("alpha","exp.alpha", "beta", "sigma.oe","mu.coef","sigma.coef",
                "coef","cov.eff","dist.coef","log.pred.rec.1","log.pred.rec.2","log.resid",
                "log.resid.2","phi")

##### Run the model with Jags #####


#bring in randomly sampled spawner-recruit data

load("Data/spawner_array_sep.Rdata")
load("Data/recruit_array_sep.Rdata")


#popn order is Carmacks, Lower Main, Middle Main, Pelly, Stewart, Teslin, Upper Lakes, White-Donjek



# Make a list for Jags ----------------------------------------------------

spawner_array_num <- spawner_array_num[,,1:100]
recruit_array_num <- recruit_array_num[,,1:100]

##Model 1##

Spawners <- spawner_array_num[,,1]
Recruits <- recruit_array_num[,,1]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out1 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out1, file=file.path(dir.output,"random sampling models 2 2/out1.rds"))
rm(out1)

##Model 2##
Spawners <- spawner_array_num[,,2]
Recruits <- recruit_array_num[,,2]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out2 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out2, file=file.path(dir.output,"random sampling models 2/out2.rds"))

rm(out2)

##Model 3##

Spawners <- spawner_array_num[,,3]
Recruits <- recruit_array_num[,,3]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out3 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out3, file=file.path(dir.output,"random sampling models 2/out3.rds"))
rm(out3)

##Model 4##
Spawners <- spawner_array_num[,,4]
Recruits <- recruit_array_num[,,4]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out4 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out4, file=file.path(dir.output,"random sampling models 2/out4.rds"))
rm(out4)

##Model 5##

Spawners <- spawner_array_num[,,5]
Recruits <- recruit_array_num[,,5]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out5 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out5, file=file.path(dir.output,"random sampling models 2/out5.rds"))
rm(out5)

##Model 6##
Spawners <- spawner_array_num[,,6]
Recruits <- recruit_array_num[,,6]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out6 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out6, file=file.path(dir.output,"random sampling models 2/out6.rds"))
rm(out6)


##Model 7##

Spawners <- spawner_array_num[,,7]
Recruits <- recruit_array_num[,,7]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out7 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out7, file=file.path(dir.output,"random sampling models 2/out7.rds"))
rm(out7)

##Model 4##
Spawners <- spawner_array_num[,,8]
Recruits <- recruit_array_num[,,8]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out8 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out8, file=file.path(dir.output,"random sampling models 2/out8.rds"))
rm(out8)

##Model 9##

Spawners <- spawner_array_num[,,9]
Recruits <- recruit_array_num[,,9]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out9 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out9, file=file.path(dir.output,"random sampling models 2/out9.rds"))
rm(out9)

##Model 10##
Spawners <- spawner_array_num[,,10]
Recruits <- recruit_array_num[,,10]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out10 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out10, file=file.path(dir.output,"random sampling models 2/out10.rds"))
rm(out10)


##Model 11##

Spawners <- spawner_array_num[,,11]
Recruits <- recruit_array_num[,,11]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out11 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out11, file=file.path(dir.output,"random sampling models 2/out11.rds"))
rm(out11)

##Model 12##
Spawners <- spawner_array_num[,,12]
Recruits <- recruit_array_num[,,12]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out12 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out12, file=file.path(dir.output,"random sampling models 2/out12.rds"))
rm(out12)


##Model 1##

Spawners <- spawner_array_num[,,13]
Recruits <- recruit_array_num[,,13]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out13 <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=50000, 
  n.burnin=10000)

#Save
saveRDS(out13, file=file.path(dir.output,"random sampling models 2/out13.rds"))
rm(out13)

##Model 2##
Spawners <- spawner_array_num[,,14]
Recruits <- recruit_array_num[,,14]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out14 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)


#Save
saveRDS(out14, file=file.path(dir.output,"random sampling models 2/out14.rds"))
rm(out14)


##Model 3##

Spawners <- spawner_array_num[,,15]
Recruits <- recruit_array_num[,,15]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out15 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)

#Save
saveRDS(out15, file=file.path(dir.output,"random sampling models 2/out15.rds"))
rm(out15)

##Model 4##
Spawners <- spawner_array_num[,,16]
Recruits <- recruit_array_num[,,16]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out16 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)


#Save
saveRDS(out16, file=file.path(dir.output,"random sampling models 2/out16.rds"))
rm(out16)

##Model 5##

Spawners <- spawner_array_num[,,17]
Recruits <- recruit_array_num[,,17]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out17 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)

#Save
saveRDS(out17, file=file.path(dir.output,"random sampling models 2/out17.rds"))
rm(out17)

##Model 6##
Spawners <- spawner_array_num[,,18]
Recruits <- recruit_array_num[,,18]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out18 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)


#Save
saveRDS(out18, file=file.path(dir.output,"random sampling models 2/out18.rds"))
rm(out18)


##Model 7##

Spawners <- spawner_array_num[,,19]
Recruits <- recruit_array_num[,,19]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out19 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)

#Save
saveRDS(out19, file=file.path(dir.output,"random sampling models 2/out19.rds"))
rm(out19)

##Model 4##
Spawners <- spawner_array_num[,,20]
Recruits <- recruit_array_num[,,20]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out20 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)


#Save
saveRDS(out20, file=file.path(dir.output,"random sampling models 2/out20.rds"))
rm(out20)

##Model 9##

Spawners <- spawner_array_num[,,21]
Recruits <- recruit_array_num[,,21]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out21 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)

#Save
saveRDS(out21, file=file.path(dir.output,"random sampling models 2/out21.rds"))
rm(out21)

##Model 10##
Spawners <- spawner_array_num[,,22]
Recruits <- recruit_array_num[,,22]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out22 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)


#Save
saveRDS(out22, file=file.path(dir.output,"random sampling models 2/out22.rds"))
rm(out22)


##Model 11##

Spawners <- spawner_array_num[,,23]
Recruits <- recruit_array_num[,,23]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out23 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)

#Save
saveRDS(out23, file=file.path(dir.output,"random sampling models 2/out23.rds"))
rm(out23)

##Model 12##
Spawners <- spawner_array_num[,,24]
Recruits <- recruit_array_num[,,24]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out24 <- jags.parallel(data=jags_data,
                      model.file=jags_model,
                      inits=jags_inits,
                      parameters.to.save=jags_params,
                      n.chains=3, 
                      n.thin=20, 
                      n.iter=50000, 
                      n.burnin=10000)


#Save
saveRDS(out24, file=file.path(dir.output,"random sampling models 2/out24.rds"))
rm(out24)




##Model 1##

Spawners <- spawner_array_num[,,25]
Recruits <- recruit_array_num[,,25]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out25 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out25, file=file.path(dir.output,"random sampling models 2/out25.rds"))
rm(out25)

##Model 2##
Spawners <- spawner_array_num[,,26]
Recruits <- recruit_array_num[,,26]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out26 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out26, file=file.path(dir.output,"random sampling models 2/out26.rds"))
rm(out26)


##Model 3##

Spawners <- spawner_array_num[,,27]
Recruits <- recruit_array_num[,,27]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out27 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out27, file=file.path(dir.output,"random sampling models 2/out27.rds"))
rm(out27)

##Model 4##
Spawners <- spawner_array_num[,,28]
Recruits <- recruit_array_num[,,28]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out28 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out28, file=file.path(dir.output,"random sampling models 2/out28.rds"))
rm(out28)

##Model 5##

Spawners <- spawner_array_num[,,29]
Recruits <- recruit_array_num[,,29]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out29 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out29, file=file.path(dir.output,"random sampling models 2/out29.rds"))
rm(out29)

##Model 6##
Spawners <- spawner_array_num[,,30]
Recruits <- recruit_array_num[,,30]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out30 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out30, file=file.path(dir.output,"random sampling models 2/out30.rds"))
rm(out30)


##Model 7##

Spawners <- spawner_array_num[,,31]
Recruits <- recruit_array_num[,,31]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out31 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out31, file=file.path(dir.output,"random sampling models 2/out31.rds"))
rm(out31)

##Model 4##
Spawners <- spawner_array_num[,,32]
Recruits <- recruit_array_num[,,32]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out32 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out32, file=file.path(dir.output,"random sampling models 2/out32.rds"))
rm(out32)

##Model 9##

Spawners <- spawner_array_num[,,33]
Recruits <- recruit_array_num[,,33]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out33 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out33, file=file.path(dir.output,"random sampling models 2/out33.rds"))
rm(out33)

##Model 10##
Spawners <- spawner_array_num[,,34]
Recruits <- recruit_array_num[,,34]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out34 <- jags.parallel(data=jags_data,
                        model.file=jags_model,
                        inits=jags_inits,
                        parameters.to.save=jags_params,
                        n.chains=3, 
                        n.thin=20, 
                        n.iter=50000, 
                        n.burnin=10000)


#Save
saveRDS(out34, file=file.path(dir.output,"random sampling models 2/out34.rds"))
rm(out34)


##Model 11##

Spawners <- spawner_array_num[,,35]
Recruits <- recruit_array_num[,,35]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out35 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out35, file=file.path(dir.output,"random sampling models 2/out35.rds"))
rm(out35)

##Model 12##
Spawners <- spawner_array_num[,,36]
Recruits <- recruit_array_num[,,36]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out36 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out36, file=file.path(dir.output,"random sampling models 2/out36.rds"))
rm(out36)





##Model 1##

Spawners <- spawner_array_num[,,37]
Recruits <- recruit_array_num[,,37]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out37 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out37, file=file.path(dir.output,"random sampling models 2/out37.rds"))
rm(out37)

##Model 2##
Spawners <- spawner_array_num[,,38]
Recruits <- recruit_array_num[,,38]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out38 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out38, file=file.path(dir.output,"random sampling models 2/out38.rds"))
rm(out38)


##Model 3##

Spawners <- spawner_array_num[,,39]
Recruits <- recruit_array_num[,,39]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out39 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out39, file=file.path(dir.output,"random sampling models 2/out39.rds"))
rm(out39)

##Model 4##
Spawners <- spawner_array_num[,,40]
Recruits <- recruit_array_num[,,40]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out40 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out40, file=file.path(dir.output,"random sampling models 2/out40.rds"))
rm(out40)

##Model 5##

Spawners <- spawner_array_num[,,41]
Recruits <- recruit_array_num[,,41]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out41 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out41, file=file.path(dir.output,"random sampling models 2/out41.rds"))

rm(out41)

##Model 6##
Spawners <- spawner_array_num[,,42]
Recruits <- recruit_array_num[,,42]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out42 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out42, file=file.path(dir.output,"random sampling models 2/out42.rds"))

rm(out42)

##Model 7##

Spawners <- spawner_array_num[,,43]
Recruits <- recruit_array_num[,,43]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out43 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out43, file=file.path(dir.output,"random sampling models 2/out43.rds"))

rm(out43)
##Model 4##
Spawners <- spawner_array_num[,,44]
Recruits <- recruit_array_num[,,44]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out44 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out44, file=file.path(dir.output,"random sampling models 2/out44.rds"))

rm(out44)
##Model 9##

Spawners <- spawner_array_num[,,45]
Recruits <- recruit_array_num[,,45]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out45 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out45, file=file.path(dir.output,"random sampling models 2/out45.rds"))
rm(out45)
##Model 10##
Spawners <- spawner_array_num[,,46]
Recruits <- recruit_array_num[,,46]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out46 <- jags.parallel(data=jags_data,
                        model.file=jags_model,
                        inits=jags_inits,
                        parameters.to.save=jags_params,
                        n.chains=3, 
                        n.thin=20, 
                        n.iter=50000, 
                        n.burnin=10000)


#Save
saveRDS(out46, file=file.path(dir.output,"random sampling models 2/out46.rds"))
rm(out46)


##Model 11##

Spawners <- spawner_array_num[,,47]
Recruits <- recruit_array_num[,,47]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out47 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out47, file=file.path(dir.output,"random sampling models 2/out47.rds"))

##Model 12##
Spawners <- spawner_array_num[,,48]
Recruits <- recruit_array_num[,,48]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out48 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out48, file=file.path(dir.output,"random sampling models 2/out48.rds"))
rm(out48)



##Model 1##

Spawners <- spawner_array_num[,,49]
Recruits <- recruit_array_num[,,49]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out49 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out49, file=file.path(dir.output,"random sampling models 2/out49.rds"))
rm(out49)
##Model 2##
Spawners <- spawner_array_num[,,50]
Recruits <- recruit_array_num[,,50]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out50 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out50, file=file.path(dir.output,"random sampling models 2/out50.rds"))
rm(out50)


##Model 3##

Spawners <- spawner_array_num[,,51]
Recruits <- recruit_array_num[,,51]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out51 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out51, file=file.path(dir.output,"random sampling models 2/out51.rds"))
rm(out51)
##Model 4##
Spawners <- spawner_array_num[,,52]
Recruits <- recruit_array_num[,,52]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out52 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out52, file=file.path(dir.output,"random sampling models 2/out52.rds"))
rm(out52)

##Model 5##

Spawners <- spawner_array_num[,,53]
Recruits <- recruit_array_num[,,53]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out53 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out53, file=file.path(dir.output,"random sampling models 2/out53.rds"))
rm(out53)
##Model 6##
Spawners <- spawner_array_num[,,54]
Recruits <- recruit_array_num[,,54]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out54 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out54, file=file.path(dir.output,"random sampling models 2/out54.rds"))
rm(out54)


##Model 7##

Spawners <- spawner_array_num[,,55]
Recruits <- recruit_array_num[,,55]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out55 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out55, file=file.path(dir.output,"random sampling models 2/out55.rds"))
rm(out55)
##Model 4##
Spawners <- spawner_array_num[,,56]
Recruits <- recruit_array_num[,,56]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out56 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out56, file=file.path(dir.output,"random sampling models 2/out56.rds"))
rm(out56)
##Model 9##

Spawners <- spawner_array_num[,,57]
Recruits <- recruit_array_num[,,57]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out57 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out57, file=file.path(dir.output,"random sampling models 2/out57.rds"))

##Model 10##
Spawners <- spawner_array_num[,,58]
Recruits <- recruit_array_num[,,58]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out58 <- jags.parallel(data=jags_data,
                        model.file=jags_model,
                        inits=jags_inits,
                        parameters.to.save=jags_params,
                        n.chains=3, 
                        n.thin=20, 
                        n.iter=50000, 
                        n.burnin=10000)


#Save
saveRDS(out58, file=file.path(dir.output,"random sampling models 2/out58.rds"))
rm(out58)


##Model 11##

Spawners <- spawner_array_num[,,59]
Recruits <- recruit_array_num[,,59]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out59 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out59, file=file.path(dir.output,"random sampling models 2/out59.rds"))
rm(out59)
##Model 12##
Spawners <- spawner_array_num[,,60]
Recruits <- recruit_array_num[,,60]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out60 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out60, file=file.path(dir.output,"random sampling models 2/out60.rds"))
rm(out60)




##Model 1##

Spawners <- spawner_array_num[,,61]
Recruits <- recruit_array_num[,,61]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out61 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out61, file=file.path(dir.output,"random sampling models 2/out61.rds"))
rm(out61)
##Model 2##
Spawners <- spawner_array_num[,,62]
Recruits <- recruit_array_num[,,62]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out62 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out62, file=file.path(dir.output,"random sampling models 2/out62.rds"))
rm(out62)


##Model 3##

Spawners <- spawner_array_num[,,63]
Recruits <- recruit_array_num[,,63]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out63 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out63, file=file.path(dir.output,"random sampling models 2/out63.rds"))
rm(out63)
##Model 4##
Spawners <- spawner_array_num[,,64]
Recruits <- recruit_array_num[,,64]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out64 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out64, file=file.path(dir.output,"random sampling models 2/out64.rds"))
rm(out64)

##Model 5##

Spawners <- spawner_array_num[,,65]
Recruits <- recruit_array_num[,,65]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out65 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out65, file=file.path(dir.output,"random sampling models 2/out65.rds"))
rm(out65)
##Model 6##
Spawners <- spawner_array_num[,,66]
Recruits <- recruit_array_num[,,66]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out66 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out66, file=file.path(dir.output,"random sampling models 2/out66.rds"))
rm(out66)


##Model 7##

Spawners <- spawner_array_num[,,67]
Recruits <- recruit_array_num[,,67]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out67 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out67, file=file.path(dir.output,"random sampling models 2/out67.rds"))
rm(out67)
##Model 4##
Spawners <- spawner_array_num[,,68]
Recruits <- recruit_array_num[,,68]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out68 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out68, file=file.path(dir.output,"random sampling models 2/out68.rds"))
rm(out68)
##Model 9##

Spawners <- spawner_array_num[,,69]
Recruits <- recruit_array_num[,,69]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out69 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out69, file=file.path(dir.output,"random sampling models 2/out69.rds"))
rm(out69)
##Model 10##
Spawners <- spawner_array_num[,,70]
Recruits <- recruit_array_num[,,70]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out70 <- jags.parallel(data=jags_data,
                        model.file=jags_model,
                        inits=jags_inits,
                        parameters.to.save=jags_params,
                        n.chains=3, 
                        n.thin=20, 
                        n.iter=50000, 
                        n.burnin=10000)


#Save
saveRDS(out70, file=file.path(dir.output,"random sampling models 2/out70.rds"))
rm(out70)


##Model 11##

Spawners <- spawner_array_num[,,71]
Recruits <- recruit_array_num[,,71]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out71 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out71, file=file.path(dir.output,"random sampling models 2/out71.rds"))
rm(out71)
##Model 12##
Spawners <- spawner_array_num[,,72]
Recruits <- recruit_array_num[,,72]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out72 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out72, file=file.path(dir.output,"random sampling models 2/out72.rds"))
rm(out72)





##Model 1##

Spawners <- spawner_array_num[,,73]
Recruits <- recruit_array_num[,,73]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out73 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out73, file=file.path(dir.output,"random sampling models 2/out73.rds"))
rm(out73)
##Model 2##
Spawners <- spawner_array_num[,,74]
Recruits <- recruit_array_num[,,74]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out74 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out74, file=file.path(dir.output,"random sampling models 2/out74.rds"))
rm(out74)


##Model 3##

Spawners <- spawner_array_num[,,75]
Recruits <- recruit_array_num[,,75]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out75 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out75, file=file.path(dir.output,"random sampling models 2/out75.rds"))
rm(out75)
##Model 4##
Spawners <- spawner_array_num[,,76]
Recruits <- recruit_array_num[,,76]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out76 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out76, file=file.path(dir.output,"random sampling models 2/out76.rds"))
rm(out76)

##Model 5##

Spawners <- spawner_array_num[,,77]
Recruits <- recruit_array_num[,,77]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out77 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out77, file=file.path(dir.output,"random sampling models 2/out77.rds"))
rm(out77)
##Model 6##
Spawners <- spawner_array_num[,,78]
Recruits <- recruit_array_num[,,78]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out78 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out78, file=file.path(dir.output,"random sampling models 2/out78.rds"))
rm(out78)


##Model 7##

Spawners <- spawner_array_num[,,79]
Recruits <- recruit_array_num[,,79]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out79 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out79, file=file.path(dir.output,"random sampling models 2/out79.rds"))
rm(out79)
##Model 4##
Spawners <- spawner_array_num[,,80]
Recruits <- recruit_array_num[,,80]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out80 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out80, file=file.path(dir.output,"random sampling models 2/out80.rds"))
rm(out80)
##Model 9##

Spawners <- spawner_array_num[,,81]
Recruits <- recruit_array_num[,,81]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out81 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out81, file=file.path(dir.output,"random sampling models 2/out81.rds"))
rm(out81)
##Model 10##
Spawners <- spawner_array_num[,,82]
Recruits <- recruit_array_num[,,82]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out82 <- jags.parallel(data=jags_data,
                        model.file=jags_model,
                        inits=jags_inits,
                        parameters.to.save=jags_params,
                        n.chains=3, 
                        n.thin=20, 
                        n.iter=50000, 
                        n.burnin=10000)


#Save
saveRDS(out82, file=file.path(dir.output,"random sampling models 2/out82.rds"))
rm(out82)


##Model 11##

Spawners <- spawner_array_num[,,83]
Recruits <- recruit_array_num[,,83]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out83 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out83, file=file.path(dir.output,"random sampling models 2/out83.rds"))
rm(out83)
##Model 12##
Spawners <- spawner_array_num[,,84]
Recruits <- recruit_array_num[,,84]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out84 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out84, file=file.path(dir.output,"random sampling models 2/out84.rds"))
rm(out84)

##Model 1##

Spawners <- spawner_array_num[,,85]
Recruits <- recruit_array_num[,,85]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out85 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out85, file=file.path(dir.output,"random sampling models 2/out85.rds"))
rm(out85)
##Model 2##
Spawners <- spawner_array_num[,,86]
Recruits <- recruit_array_num[,,86]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out86 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out86, file=file.path(dir.output,"random sampling models 2/out86.rds"))
rm(out86)


##Model 3##

Spawners <- spawner_array_num[,,87]
Recruits <- recruit_array_num[,,87]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out87 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out87, file=file.path(dir.output,"random sampling models 2/out87.rds"))
rm(out87)
##Model 4##
Spawners <- spawner_array_num[,,88]
Recruits <- recruit_array_num[,,88]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out88 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out88, file=file.path(dir.output,"random sampling models 2/out88.rds"))
rm(out88)

##Model 5##

Spawners <- spawner_array_num[,,89]
Recruits <- recruit_array_num[,,89]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out89 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out89, file=file.path(dir.output,"random sampling models 2/out89.rds"))
rm(out89)
##Model 6##
Spawners <- spawner_array_num[,,90]
Recruits <- recruit_array_num[,,90]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out90 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out90, file=file.path(dir.output,"random sampling models 2/out90.rds"))
rm(out90)


##Model 7##

Spawners <- spawner_array_num[,,91]
Recruits <- recruit_array_num[,,91]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out91 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out91, file=file.path(dir.output,"random sampling models 2/out91.rds"))

##Model 4##
Spawners <- spawner_array_num[,,92]
Recruits <- recruit_array_num[,,92]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out92 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out92, file=file.path(dir.output,"random sampling models 2/out92.rds"))

##Model 9##

Spawners <- spawner_array_num[,,93]
Recruits <- recruit_array_num[,,93]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out93 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out93, file=file.path(dir.output,"random sampling models 2/out93.rds"))

##Model 10##
Spawners <- spawner_array_num[,,94]
Recruits <- recruit_array_num[,,94]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out94 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out94, file=file.path(dir.output,"random sampling models 2/out94.rds"))



##Model 11##

Spawners <- spawner_array_num[,,95]
Recruits <- recruit_array_num[,,95]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out95 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out95, file=file.path(dir.output,"random sampling models 2/out95.rds"))

##Model 12##
Spawners <- spawner_array_num[,,96]
Recruits <- recruit_array_num[,,96]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out96 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out96, file=file.path(dir.output,"random sampling models 2/out96.rds"))





##Model 1##

Spawners <- spawner_array_num[,,97]
Recruits <- recruit_array_num[,,97]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out97 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out97, file=file.path(dir.output,"random sampling models 2/out97.rds"))

##Model 2##
Spawners <- spawner_array_num[,,98]
Recruits <- recruit_array_num[,,98]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out98 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out98, file=file.path(dir.output,"random sampling models 2/out98.rds"))




##Model 1##

Spawners <- spawner_array_num[,,99]
Recruits <- recruit_array_num[,,99]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out99 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)

#Save
saveRDS(out99, file=file.path(dir.output,"random sampling models 2/out99.rds"))

##Model 2##
Spawners <- spawner_array_num[,,100]
Recruits <- recruit_array_num[,,100]

ln.Recruits <- log(Recruits)

Sobs = Spawners
lnRobs = ln.Recruits

jags_data = list(
  "n.pops", "n.years", "n.covars",
  "Sobs", 
  "lnRobs",
  "covars"
  
)

#fast
out100 <- jags.parallel(data=jags_data,
                       model.file=jags_model,
                       inits=jags_inits,
                       parameters.to.save=jags_params,
                       n.chains=3, 
                       n.thin=20, 
                       n.iter=50000, 
                       n.burnin=10000)


#Save
saveRDS(out100, file=file.path(dir.output,"random sampling models 2/out100.rds"))




######load and merge models#####
out1 <- readRDS("Output/random sampling models 2/out1.rds")
out2 <- readRDS("Output/random sampling models 2/out2.rds")
out3 <- readRDS("Output/random sampling models 2/out3.rds")
out4 <- readRDS("Output/random sampling models 2/out4.rds")
out5 <- readRDS("Output/random sampling models 2/out5.rds")
out6 <- readRDS("Output/random sampling models 2/out6.rds")
out7 <- readRDS("Output/random sampling models 2/out7.rds")
out8 <- readRDS("Output/random sampling models 2/out8.rds")
out9 <- readRDS("Output/random sampling models 2/out9.rds")
out10 <- readRDS("Output/random sampling models 2/out10.rds")

mu.coef1 <- as.data.frame(out1$BUGSoutput$sims.list$mu.coef)
mu.coef2 <- as.data.frame(out2$BUGSoutput$sims.list$mu.coef)
mu.coef3 <- as.data.frame(out3$BUGSoutput$sims.list$mu.coef)
mu.coef4 <- as.data.frame(out4$BUGSoutput$sims.list$mu.coef)
mu.coef5 <- as.data.frame(out5$BUGSoutput$sims.list$mu.coef)
mu.coef6 <- as.data.frame(out6$BUGSoutput$sims.list$mu.coef)
mu.coef7 <- as.data.frame(out7$BUGSoutput$sims.list$mu.coef)
mu.coef8 <- as.data.frame(out8$BUGSoutput$sims.list$mu.coef)
mu.coef9 <- as.data.frame(out9$BUGSoutput$sims.list$mu.coef)
mu.coef10 <- as.data.frame(out10$BUGSoutput$sims.list$mu.coef)

coef1 <- as.data.frame.table(out1$BUGSoutput$sims.list$coef)
coef2 <- as.data.frame.table(out2$BUGSoutput$sims.list$coef)
coef3 <- as.data.frame.table(out3$BUGSoutput$sims.list$coef)
coef4 <- as.data.frame.table(out4$BUGSoutput$sims.list$coef)
coef5 <- as.data.frame.table(out5$BUGSoutput$sims.list$coef)
coef6 <- as.data.frame.table(out6$BUGSoutput$sims.list$coef)
coef7 <- as.data.frame.table(out7$BUGSoutput$sims.list$coef)
coef8 <- as.data.frame.table(out8$BUGSoutput$sims.list$coef)
coef9 <- as.data.frame.table(out9$BUGSoutput$sims.list$coef)
coef10 <- as.data.frame.table(out10$BUGSoutput$sims.list$coef)

out11 <- readRDS("Output/random sampling models 2/out11.rds")
out12 <- readRDS("Output/random sampling models 2/out12.rds")
out13 <- readRDS("Output/random sampling models 2/out13.rds")
out14 <- readRDS("Output/random sampling models 2/out14.rds")
out15 <- readRDS("Output/random sampling models 2/out15.rds")
out16 <- readRDS("Output/random sampling models 2/out16.rds")
out17 <- readRDS("Output/random sampling models 2/out17.rds")
out18 <- readRDS("Output/random sampling models 2/out18.rds")
out19 <- readRDS("Output/random sampling models 2/out19.rds")
out20 <- readRDS("Output/random sampling models 2/out20.rds")

mu.coef11 <- as.data.frame(out11$BUGSoutput$sims.list$mu.coef)
mu.coef12 <- as.data.frame(out12$BUGSoutput$sims.list$mu.coef)
mu.coef13 <- as.data.frame(out13$BUGSoutput$sims.list$mu.coef)
mu.coef14 <- as.data.frame(out14$BUGSoutput$sims.list$mu.coef)
mu.coef15 <- as.data.frame(out15$BUGSoutput$sims.list$mu.coef)
mu.coef16 <- as.data.frame(out16$BUGSoutput$sims.list$mu.coef)
mu.coef17 <- as.data.frame(out17$BUGSoutput$sims.list$mu.coef)
mu.coef18 <- as.data.frame(out18$BUGSoutput$sims.list$mu.coef)
mu.coef19 <- as.data.frame(out19$BUGSoutput$sims.list$mu.coef)
mu.coef20 <- as.data.frame(out20$BUGSoutput$sims.list$mu.coef)

coef11 <- as.data.frame.table(out11$BUGSoutput$sims.list$coef)
coef12 <- as.data.frame.table(out12$BUGSoutput$sims.list$coef)
coef13 <- as.data.frame.table(out13$BUGSoutput$sims.list$coef)
coef14 <- as.data.frame.table(out14$BUGSoutput$sims.list$coef)
coef15 <- as.data.frame.table(out15$BUGSoutput$sims.list$coef)
coef16 <- as.data.frame.table(out16$BUGSoutput$sims.list$coef)
coef17 <- as.data.frame.table(out17$BUGSoutput$sims.list$coef)
coef18 <- as.data.frame.table(out18$BUGSoutput$sims.list$coef)
coef19 <- as.data.frame.table(out19$BUGSoutput$sims.list$coef)
coef20 <- as.data.frame.table(out20$BUGSoutput$sims.list$coef)

mu.coef1_20 <- bind_rows(mu.coef1,mu.coef2,mu.coef3,mu.coef4,
                     mu.coef5,mu.coef6,mu.coef7,mu.coef8,
                     mu.coef9,mu.coef10,mu.coef11,mu.coef12,
                     mu.coef13,mu.coef14,mu.coef15,mu.coef16,
                     mu.coef17,mu.coef18,mu.coef19,mu.coef20)

write.csv(mu.coef1_20,"Output/random sampling models 2/mu.coef1_20.csv")

coef1_20 <- bind_rows(coef1,coef2,coef3,coef4,
                         coef5,coef6,coef7,coef8,
                         coef9,coef10,coef11,coef12,
                         coef13,coef14,coef15,coef16,
                         coef17,coef18,coef19,coef20)

write.csv(coef1_20,"Output/random sampling models 2/coef1_20.csv")

rm(list=ls())

out21 <- readRDS("Output/random sampling models 2/out21.rds")
out22 <- readRDS("Output/random sampling models 2/out22.rds")
out23 <- readRDS("Output/random sampling models 2/out23.rds")
out24 <- readRDS("Output/random sampling models 2/out24.rds")
out25 <- readRDS("Output/random sampling models 2/out25.rds")
out26 <- readRDS("Output/random sampling models 2/out26.rds")
out27 <- readRDS("Output/random sampling models 2/out27.rds")
out28 <- readRDS("Output/random sampling models 2/out28.rds")
out29 <- readRDS("Output/random sampling models 2/out29.rds")
out30 <- readRDS("Output/random sampling models 2/out30.rds")

mu.coef21 <- as.data.frame(out21$BUGSoutput$sims.list$mu.coef)
mu.coef22 <- as.data.frame(out22$BUGSoutput$sims.list$mu.coef)
mu.coef23 <- as.data.frame(out23$BUGSoutput$sims.list$mu.coef)
mu.coef24 <- as.data.frame(out24$BUGSoutput$sims.list$mu.coef)
mu.coef25 <- as.data.frame(out25$BUGSoutput$sims.list$mu.coef)
mu.coef26 <- as.data.frame(out26$BUGSoutput$sims.list$mu.coef)
mu.coef27 <- as.data.frame(out27$BUGSoutput$sims.list$mu.coef)
mu.coef28 <- as.data.frame(out28$BUGSoutput$sims.list$mu.coef)
mu.coef29 <- as.data.frame(out29$BUGSoutput$sims.list$mu.coef)
mu.coef30 <- as.data.frame(out30$BUGSoutput$sims.list$mu.coef)

coef21 <- as.data.frame.table(out21$BUGSoutput$sims.list$coef)
coef22 <- as.data.frame.table(out22$BUGSoutput$sims.list$coef)
coef23 <- as.data.frame.table(out23$BUGSoutput$sims.list$coef)
coef24 <- as.data.frame.table(out24$BUGSoutput$sims.list$coef)
coef25 <- as.data.frame.table(out25$BUGSoutput$sims.list$coef)
coef26 <- as.data.frame.table(out26$BUGSoutput$sims.list$coef)
coef27 <- as.data.frame.table(out27$BUGSoutput$sims.list$coef)
coef28 <- as.data.frame.table(out28$BUGSoutput$sims.list$coef)
coef29 <- as.data.frame.table(out29$BUGSoutput$sims.list$coef)
coef30 <- as.data.frame.table(out30$BUGSoutput$sims.list$coef)

out31 <- readRDS("Output/random sampling models 2/out31.rds")
out32 <- readRDS("Output/random sampling models 2/out32.rds")
out33 <- readRDS("Output/random sampling models 2/out33.rds")
out34 <- readRDS("Output/random sampling models 2/out34.rds")
out35 <- readRDS("Output/random sampling models 2/out35.rds")
out36 <- readRDS("Output/random sampling models 2/out36.rds")
out37 <- readRDS("Output/random sampling models 2/out37.rds")
out38 <- readRDS("Output/random sampling models 2/out38.rds")
out39 <- readRDS("Output/random sampling models 2/out39.rds")
out40 <- readRDS("Output/random sampling models 2/out40.rds")

mu.coef31 <- as.data.frame(out31$BUGSoutput$sims.list$mu.coef)
mu.coef32 <- as.data.frame(out32$BUGSoutput$sims.list$mu.coef)
mu.coef33 <- as.data.frame(out33$BUGSoutput$sims.list$mu.coef)
mu.coef34 <- as.data.frame(out34$BUGSoutput$sims.list$mu.coef)
mu.coef35 <- as.data.frame(out35$BUGSoutput$sims.list$mu.coef)
mu.coef36 <- as.data.frame(out36$BUGSoutput$sims.list$mu.coef)
mu.coef37 <- as.data.frame(out37$BUGSoutput$sims.list$mu.coef)
mu.coef38 <- as.data.frame(out38$BUGSoutput$sims.list$mu.coef)
mu.coef39 <- as.data.frame(out39$BUGSoutput$sims.list$mu.coef)
mu.coef40 <- as.data.frame(out40$BUGSoutput$sims.list$mu.coef)

coef31 <- as.data.frame.table(out31$BUGSoutput$sims.list$coef)
coef32 <- as.data.frame.table(out32$BUGSoutput$sims.list$coef)
coef33 <- as.data.frame.table(out33$BUGSoutput$sims.list$coef)
coef34 <- as.data.frame.table(out34$BUGSoutput$sims.list$coef)
coef35 <- as.data.frame.table(out35$BUGSoutput$sims.list$coef)
coef36 <- as.data.frame.table(out36$BUGSoutput$sims.list$coef)
coef37 <- as.data.frame.table(out37$BUGSoutput$sims.list$coef)
coef38 <- as.data.frame.table(out38$BUGSoutput$sims.list$coef)
coef39 <- as.data.frame.table(out39$BUGSoutput$sims.list$coef)
coef40 <- as.data.frame.table(out40$BUGSoutput$sims.list$coef)

library(tidyverse)
mu.coef21_40 <- bind_rows(mu.coef21,mu.coef22,mu.coef23,mu.coef24,
                          mu.coef25,mu.coef26,mu.coef27,mu.coef28,
                          mu.coef29,mu.coef30,mu.coef31,mu.coef32,
                          mu.coef33,mu.coef34,mu.coef35,mu.coef36,
                          mu.coef37,mu.coef38,mu.coef39,mu.coef40)

write.csv(mu.coef21_40,"Output/random sampling models 2/mu.coef21_40.csv")

coef21_40 <- bind_rows(coef21,coef22,coef23,coef24,
                          coef25,coef26,coef27,coef28,
                          coef29,coef30,coef31,coef32,
                          coef33,coef34,coef35,coef36,
                          coef37,coef38,coef39,coef40)

write.csv(coef21_40,"Output/random sampling models 2/coef21_40.csv")

rm(list=ls())

out41 <- readRDS("Output/random sampling models 2/out41.rds")
out42 <- readRDS("Output/random sampling models 2/out42.rds")
out43 <- readRDS("Output/random sampling models 2/out43.rds")
out44 <- readRDS("Output/random sampling models 2/out44.rds")
out45 <- readRDS("Output/random sampling models 2/out45.rds")
out46 <- readRDS("Output/random sampling models 2/out46.rds")
out47 <- readRDS("Output/random sampling models 2/out47.rds")
out48 <- readRDS("Output/random sampling models 2/out48.rds")
out49 <- readRDS("Output/random sampling models 2/out49.rds")
out50 <- readRDS("Output/random sampling models 2/out50.rds")

mu.coef41 <- as.data.frame(out41$BUGSoutput$sims.list$mu.coef)
mu.coef42 <- as.data.frame(out42$BUGSoutput$sims.list$mu.coef)
mu.coef43 <- as.data.frame(out43$BUGSoutput$sims.list$mu.coef)
mu.coef44 <- as.data.frame(out44$BUGSoutput$sims.list$mu.coef)
mu.coef45 <- as.data.frame(out45$BUGSoutput$sims.list$mu.coef)
mu.coef46 <- as.data.frame(out46$BUGSoutput$sims.list$mu.coef)
mu.coef47 <- as.data.frame(out47$BUGSoutput$sims.list$mu.coef)
mu.coef48 <- as.data.frame(out48$BUGSoutput$sims.list$mu.coef)
mu.coef49 <- as.data.frame(out49$BUGSoutput$sims.list$mu.coef)
mu.coef50 <- as.data.frame(out50$BUGSoutput$sims.list$mu.coef)

coef41 <- as.data.frame.table(out41$BUGSoutput$sims.list$coef)
coef42 <- as.data.frame.table(out42$BUGSoutput$sims.list$coef)
coef43 <- as.data.frame.table(out43$BUGSoutput$sims.list$coef)
coef44 <- as.data.frame.table(out44$BUGSoutput$sims.list$coef)
coef45 <- as.data.frame.table(out45$BUGSoutput$sims.list$coef)
coef46 <- as.data.frame.table(out46$BUGSoutput$sims.list$coef)
coef47 <- as.data.frame.table(out47$BUGSoutput$sims.list$coef)
coef48 <- as.data.frame.table(out48$BUGSoutput$sims.list$coef)
coef49 <- as.data.frame.table(out49$BUGSoutput$sims.list$coef)
coef50 <- as.data.frame.table(out50$BUGSoutput$sims.list$coef)

out51 <- readRDS("Output/random sampling models 2/out51.rds")
out52 <- readRDS("Output/random sampling models 2/out52.rds")
out53 <- readRDS("Output/random sampling models 2/out53.rds")
out54 <- readRDS("Output/random sampling models 2/out54.rds")
out55 <- readRDS("Output/random sampling models 2/out55.rds")
out56 <- readRDS("Output/random sampling models 2/out56.rds")
out57 <- readRDS("Output/random sampling models 2/out57.rds")
out58 <- readRDS("Output/random sampling models 2/out58.rds")
out59 <- readRDS("Output/random sampling models 2/out59.rds")
out60 <- readRDS("Output/random sampling models 2/out60.rds")

mu.coef51 <- as.data.frame(out51$BUGSoutput$sims.list$mu.coef)
mu.coef52 <- as.data.frame(out52$BUGSoutput$sims.list$mu.coef)
mu.coef53 <- as.data.frame(out53$BUGSoutput$sims.list$mu.coef)
mu.coef54 <- as.data.frame(out54$BUGSoutput$sims.list$mu.coef)
mu.coef55 <- as.data.frame(out55$BUGSoutput$sims.list$mu.coef)
mu.coef56 <- as.data.frame(out56$BUGSoutput$sims.list$mu.coef)
mu.coef57 <- as.data.frame(out57$BUGSoutput$sims.list$mu.coef)
mu.coef58 <- as.data.frame(out58$BUGSoutput$sims.list$mu.coef)
mu.coef59 <- as.data.frame(out59$BUGSoutput$sims.list$mu.coef)
mu.coef60 <- as.data.frame(out60$BUGSoutput$sims.list$mu.coef)

coef51 <- as.data.frame.table(out51$BUGSoutput$sims.list$coef)
coef52 <- as.data.frame.table(out52$BUGSoutput$sims.list$coef)
coef53 <- as.data.frame.table(out53$BUGSoutput$sims.list$coef)
coef54 <- as.data.frame.table(out54$BUGSoutput$sims.list$coef)
coef55 <- as.data.frame.table(out55$BUGSoutput$sims.list$coef)
coef56 <- as.data.frame.table(out56$BUGSoutput$sims.list$coef)
coef57 <- as.data.frame.table(out57$BUGSoutput$sims.list$coef)
coef58 <- as.data.frame.table(out58$BUGSoutput$sims.list$coef)
coef59 <- as.data.frame.table(out59$BUGSoutput$sims.list$coef)
coef60 <- as.data.frame.table(out60$BUGSoutput$sims.list$coef)

library(tidyverse)
mu.coef41_60 <- bind_rows(mu.coef41,mu.coef42,mu.coef43,mu.coef44,
                          mu.coef45,mu.coef46,mu.coef47,mu.coef48,
                          mu.coef49,mu.coef50,mu.coef51,mu.coef52,
                          mu.coef53,mu.coef54,mu.coef55,mu.coef56,
                          mu.coef57,mu.coef58,mu.coef59,mu.coef60)

write.csv(mu.coef41_60,"Output/random sampling models 2/mu.coef41_60.csv")

coef41_60 <- bind_rows(coef41,coef42,coef43,coef44,
                          coef45,coef46,coef47,coef48,
                          coef49,coef50,coef51,coef52,
                          coef53,coef54,coef55,coef56,
                          coef57,coef58,coef59,coef60)

write.csv(coef41_60,"Output/random sampling models 2/coef41_60.csv")

rm(list=ls())

out61 <- readRDS("Output/random sampling models 2/out61.rds")
out62 <- readRDS("Output/random sampling models 2/out62.rds")
out63 <- readRDS("Output/random sampling models 2/out63.rds")
out64 <- readRDS("Output/random sampling models 2/out64.rds")
out65 <- readRDS("Output/random sampling models 2/out65.rds")
out66 <- readRDS("Output/random sampling models 2/out66.rds")
out67 <- readRDS("Output/random sampling models 2/out67.rds")
out68 <- readRDS("Output/random sampling models 2/out68.rds")
out69 <- readRDS("Output/random sampling models 2/out69.rds")
out70 <- readRDS("Output/random sampling models 2/out70.rds")


mu.coef61 <- as.data.frame(out61$BUGSoutput$sims.list$mu.coef)
mu.coef62 <- as.data.frame(out62$BUGSoutput$sims.list$mu.coef)
mu.coef63 <- as.data.frame(out63$BUGSoutput$sims.list$mu.coef)
mu.coef64 <- as.data.frame(out64$BUGSoutput$sims.list$mu.coef)
mu.coef65 <- as.data.frame(out65$BUGSoutput$sims.list$mu.coef)
mu.coef66 <- as.data.frame(out66$BUGSoutput$sims.list$mu.coef)
mu.coef67 <- as.data.frame(out67$BUGSoutput$sims.list$mu.coef)
mu.coef68 <- as.data.frame(out68$BUGSoutput$sims.list$mu.coef)
mu.coef69 <- as.data.frame(out69$BUGSoutput$sims.list$mu.coef)
mu.coef70 <- as.data.frame(out70$BUGSoutput$sims.list$mu.coef)

coef61 <- as.data.frame.table(out61$BUGSoutput$sims.list$coef)
coef62 <- as.data.frame.table(out62$BUGSoutput$sims.list$coef)
coef63 <- as.data.frame.table(out63$BUGSoutput$sims.list$coef)
coef64 <- as.data.frame.table(out64$BUGSoutput$sims.list$coef)
coef65 <- as.data.frame.table(out65$BUGSoutput$sims.list$coef)
coef66 <- as.data.frame.table(out66$BUGSoutput$sims.list$coef)
coef67 <- as.data.frame.table(out67$BUGSoutput$sims.list$coef)
coef68 <- as.data.frame.table(out68$BUGSoutput$sims.list$coef)
coef69 <- as.data.frame.table(out69$BUGSoutput$sims.list$coef)
coef70 <- as.data.frame.table(out70$BUGSoutput$sims.list$coef)

out71 <- readRDS("Output/random sampling models 2/out71.rds")
out72 <- readRDS("Output/random sampling models 2/out72.rds")
out73 <- readRDS("Output/random sampling models 2/out73.rds")
out74 <- readRDS("Output/random sampling models 2/out74.rds")
out75 <- readRDS("Output/random sampling models 2/out75.rds")
out76 <- readRDS("Output/random sampling models 2/out76.rds")
out77 <- readRDS("Output/random sampling models 2/out77.rds")
out78 <- readRDS("Output/random sampling models 2/out78.rds")
out79 <- readRDS("Output/random sampling models 2/out79.rds")
out80 <- readRDS("Output/random sampling models 2/out80.rds")

mu.coef71 <- as.data.frame(out71$BUGSoutput$sims.list$mu.coef)
mu.coef72 <- as.data.frame(out72$BUGSoutput$sims.list$mu.coef)
mu.coef73 <- as.data.frame(out73$BUGSoutput$sims.list$mu.coef)
mu.coef74 <- as.data.frame(out74$BUGSoutput$sims.list$mu.coef)
mu.coef75 <- as.data.frame(out75$BUGSoutput$sims.list$mu.coef)
mu.coef76 <- as.data.frame(out76$BUGSoutput$sims.list$mu.coef)
mu.coef77 <- as.data.frame(out77$BUGSoutput$sims.list$mu.coef)
mu.coef78 <- as.data.frame(out78$BUGSoutput$sims.list$mu.coef)
mu.coef79 <- as.data.frame(out79$BUGSoutput$sims.list$mu.coef)
mu.coef80 <- as.data.frame(out80$BUGSoutput$sims.list$mu.coef)

coef71 <- as.data.frame.table(out71$BUGSoutput$sims.list$coef)
coef72 <- as.data.frame.table(out72$BUGSoutput$sims.list$coef)
coef73 <- as.data.frame.table(out73$BUGSoutput$sims.list$coef)
coef74 <- as.data.frame.table(out74$BUGSoutput$sims.list$coef)
coef75 <- as.data.frame.table(out75$BUGSoutput$sims.list$coef)
coef76 <- as.data.frame.table(out76$BUGSoutput$sims.list$coef)
coef77 <- as.data.frame.table(out77$BUGSoutput$sims.list$coef)
coef78 <- as.data.frame.table(out78$BUGSoutput$sims.list$coef)
coef79 <- as.data.frame.table(out79$BUGSoutput$sims.list$coef)
coef80 <- as.data.frame.table(out80$BUGSoutput$sims.list$coef)

library(tidyverse)
mu.coef61_80 <- bind_rows(mu.coef61,mu.coef62,mu.coef63,mu.coef64,
                          mu.coef65,mu.coef66,mu.coef67,mu.coef68,
                          mu.coef69,mu.coef70,mu.coef71,mu.coef72,
                          mu.coef73,mu.coef74,mu.coef75,mu.coef76,
                          mu.coef77,mu.coef78,mu.coef79,mu.coef80)

write.csv(mu.coef61_80,"Output/random sampling models 2/mu.coef61_80.csv")

coef61_80 <- bind_rows(coef61,coef62,coef63,coef64,
                          coef65,coef66,coef67,coef68,
                          coef69,coef70,coef71,coef72,
                          coef73,coef74,coef75,coef76,
                          coef77,coef78,coef79,coef80)

write.csv(coef61_80,"Output/random sampling models 2/coef61_80.csv")

rm(list=ls())

out81 <- readRDS("Output/random sampling models 2/out81.rds")
out82 <- readRDS("Output/random sampling models 2/out82.rds")
out83 <- readRDS("Output/random sampling models 2/out83.rds")
out84 <- readRDS("Output/random sampling models 2/out84.rds")
out85 <- readRDS("Output/random sampling models 2/out85.rds")
out86 <- readRDS("Output/random sampling models 2/out86.rds")
out87 <- readRDS("Output/random sampling models 2/out87.rds")
out88 <- readRDS("Output/random sampling models 2/out88.rds")
out89 <- readRDS("Output/random sampling models 2/out89.rds")
out90 <- readRDS("Output/random sampling models 2/out90.rds")

mu.coef81 <- as.data.frame(out81$BUGSoutput$sims.list$mu.coef)
mu.coef82 <- as.data.frame(out82$BUGSoutput$sims.list$mu.coef)
mu.coef83 <- as.data.frame(out83$BUGSoutput$sims.list$mu.coef)
mu.coef84 <- as.data.frame(out84$BUGSoutput$sims.list$mu.coef)
mu.coef85 <- as.data.frame(out85$BUGSoutput$sims.list$mu.coef)
mu.coef86 <- as.data.frame(out86$BUGSoutput$sims.list$mu.coef)
mu.coef87 <- as.data.frame(out87$BUGSoutput$sims.list$mu.coef)
mu.coef88 <- as.data.frame(out88$BUGSoutput$sims.list$mu.coef)
mu.coef89 <- as.data.frame(out89$BUGSoutput$sims.list$mu.coef)
mu.coef90 <- as.data.frame(out90$BUGSoutput$sims.list$mu.coef)

coef81 <- as.data.frame.table(out81$BUGSoutput$sims.list$coef)
coef82 <- as.data.frame.table(out82$BUGSoutput$sims.list$coef)
coef83 <- as.data.frame.table(out83$BUGSoutput$sims.list$coef)
coef84 <- as.data.frame.table(out84$BUGSoutput$sims.list$coef)
coef85 <- as.data.frame.table(out85$BUGSoutput$sims.list$coef)
coef86 <- as.data.frame.table(out86$BUGSoutput$sims.list$coef)
coef87 <- as.data.frame.table(out87$BUGSoutput$sims.list$coef)
coef88 <- as.data.frame.table(out88$BUGSoutput$sims.list$coef)
coef89 <- as.data.frame.table(out89$BUGSoutput$sims.list$coef)
coef90 <- as.data.frame.table(out90$BUGSoutput$sims.list$coef)

out91 <- readRDS("Output/random sampling models 2/out91.rds")
out92 <- readRDS("Output/random sampling models 2/out92.rds")
out93 <- readRDS("Output/random sampling models 2/out93.rds")
out94 <- readRDS("Output/random sampling models 2/out94.rds")
out95 <- readRDS("Output/random sampling models 2/out95.rds")
out96 <- readRDS("Output/random sampling models 2/out96.rds")
out97 <- readRDS("Output/random sampling models 2/out97.rds")
out98 <- readRDS("Output/random sampling models 2/out98.rds")
out99 <- readRDS("Output/random sampling models 2/out99.rds")
out100 <- readRDS("Output/random sampling models 2/out100.rds")

mu.coef91 <- as.data.frame(out91$BUGSoutput$sims.list$mu.coef)
mu.coef92 <- as.data.frame(out92$BUGSoutput$sims.list$mu.coef)
mu.coef93 <- as.data.frame(out93$BUGSoutput$sims.list$mu.coef)
mu.coef94 <- as.data.frame(out94$BUGSoutput$sims.list$mu.coef)
mu.coef95 <- as.data.frame(out95$BUGSoutput$sims.list$mu.coef)
mu.coef96 <- as.data.frame(out96$BUGSoutput$sims.list$mu.coef)
mu.coef97 <- as.data.frame(out97$BUGSoutput$sims.list$mu.coef)
mu.coef98 <- as.data.frame(out98$BUGSoutput$sims.list$mu.coef)
mu.coef99 <- as.data.frame(out99$BUGSoutput$sims.list$mu.coef)
mu.coef100 <- as.data.frame(out100$BUGSoutput$sims.list$mu.coef)

coef91 <- as.data.frame.table(out91$BUGSoutput$sims.list$coef)
coef92 <- as.data.frame.table(out92$BUGSoutput$sims.list$coef)
coef93 <- as.data.frame.table(out93$BUGSoutput$sims.list$coef)
coef94 <- as.data.frame.table(out94$BUGSoutput$sims.list$coef)
coef95 <- as.data.frame.table(out95$BUGSoutput$sims.list$coef)
coef96 <- as.data.frame.table(out96$BUGSoutput$sims.list$coef)
coef97 <- as.data.frame.table(out97$BUGSoutput$sims.list$coef)
coef98 <- as.data.frame.table(out98$BUGSoutput$sims.list$coef)
coef99 <- as.data.frame.table(out99$BUGSoutput$sims.list$coef)
coef100 <- as.data.frame.table(out100$BUGSoutput$sims.list$coef)

library(tidyverse)
mu.coef81_100 <- bind_rows(mu.coef81,mu.coef82,mu.coef83,mu.coef84,
                           mu.coef85,mu.coef86,mu.coef87,mu.coef88,
                           mu.coef89,mu.coef90,mu.coef91,mu.coef92,
                           mu.coef93,mu.coef94,mu.coef95,mu.coef96,
                           mu.coef97,mu.coef98,mu.coef99,mu.coef100)

write.csv(mu.coef81_100,"Output/random sampling models 2/mu.coef81_100.csv")

coef81_100 <- bind_rows(coef81,coef82,coef83,coef84,
                           coef85,coef86,coef87,coef88,
                           coef89,coef90,coef91,coef92,
                           coef93,coef94,coef95,coef96,
                           coef97,coef98,coef99,coef100)

write.csv(coef81_100,"Output/random sampling models 2/coef81_100.csv")

rm(list=ls())

library(tidyverse)
mu1 <- read.csv("Output/random sampling models 2/mu.coef1_20.csv")
mu2 <- read.csv("Output/random sampling models 2/mu.coef21_40.csv")
mu3 <- read.csv("Output/random sampling models 2/mu.coef41_60.csv")
mu4 <- read.csv("Output/random sampling models 2/mu.coef61_80.csv")
mu5 <- read.csv("Output/random sampling models 2/mu.coef81_100.csv")

mu_all <- bind_rows(mu1,mu2,mu3,mu4,mu5)
mu_all <- mu_all[,2:13]

write.csv(mu_all,"Output/mu_all_sep.csv")

coef1 <- read.csv("Output/random sampling models 2/coef1_20.csv")
coef2 <- read.csv("Output/random sampling models 2/coef21_40.csv")
coef3 <- read.csv("Output/random sampling models 2/coef41_60.csv")
coef4 <- read.csv("Output/random sampling models 2/coef61_80.csv")
coef5 <- read.csv("Output/random sampling models 2/coef81_100.csv")

coef_all <- bind_rows(coef1,coef2,coef3,coef4,coef5)
coef_all <- coef_all[,2:5]

write.csv(coef_all,"Output/coef_all_sep.csv")



# Visualize results -------------------------------------------------------


par(mfcol=c(1,1), mar=c(2,12,3,1), oma=c(2,2,1,1))
c <- 1
caterplot(mu_all,
          labels=names.covars,reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), 
          style='gray', col='blue',cex=1.1)
caterpoints(apply(mu_all,2,median), pch=21, col='red', bg='orange')
abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5)
mtext('Covariate', side=2, outer=TRUE, font=2, line=0.5)


pops2=c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
        "Upper Mainstem","White-Donjek")

#rearrange variables

mu.new <- mu_all[,c(2,6,7,5,12,3,4,1,8,9,10,11)]


mu.new <- as.data.frame(mu.new)

colnames(mu.new) <- c("Migration temperature","Spawning temperature","Spawning precipitation",
                                   "Snowpack","FW winter temperature","Rearing degree days",
                                   "Rearing precipitation","Ice out date","Summer SST",
                                   "Winter SST","Pink salmon","Chum salmon")

mu.coef.vector.long <- mu.new %>% gather("var","value",1:12)

mu.coef.vector.long <- mu.coef.vector.long %>% 
  mutate(stage=if_else(var=="Chum salmon"|var=="Pink salmon"|var=="Summer SST"|
                         var=="Winter SST","Marine",
                       if_else(var=="Ice out date","Outmigration",
                               if_else(var=="Rearing degree days"|var=="Rearing precipitation",
                                       "Freshwater juvenile",
                                       "Spawning and incubation"))))

# Source Useful Functions ==================================================
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.95 <- function(x) { return(quantile(x, probs=c(0.025,0.975))) }

library(tidybayes)

# Regional_effects


Regional_effects <- mu.coef.vector.long %>% 
  mutate(var=fct_relevel(var,"Chum salmon","Pink salmon","Winter SST",
                         "Summer SST","Ice out date","Rearing precipitation",
                         "Rearing degree days","FW winter temperature","Snowpack",
                         "Spawning precipitation","Spawning temperature","Migration temperature")) %>% 
  mutate(stage=fct_relevel(stage,"Spawning and incubation","Freshwater juvenile",
                           "Outmigration","Marine")) %>% 
  ggplot(aes(x=var, y=value,fill=stage)) +
  geom_hline(yintercept = 0, colour='red')+
  stat_eye()+
  theme_bw() +
  scale_fill_manual(values = c("#FE9000","#FFDD4A","#5ADBFF","#3C6997"))+
  coord_flip() +
  stat_summary(fun="q.95", colour="black", geom="line", lwd=0.5,show.legend=FALSE) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=1.25,show.legend=FALSE) +
  stat_summary(fun="median", colour="black", size=1.75, geom="point", pch=21, show.legend=FALSE) +
  stat_summary(fun="median", colour="gray", size=1.5, geom="point", pch=19,show.legend=FALSE) +
  theme(legend.position='top') +
  theme(legend.title=element_blank())+
  xlab("Environmental variable") + ylab("Change in log (recruits/spawner)") +
  theme(text = element_text(size=20),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=15,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))+
  theme(legend.text=element_text(size=15))+scale_y_continuous(limits = c(-0.35,0.4))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

Regional_effects

#785x620 saved dimensions svg


#data structure is 1:35250 (values), 1:8 (popns), 1:12 (covars)

names(coef_all) <- c("reps","population","var","value")

#rename factor levels
coef_all$population <- as.factor(coef_all$population)
levels(coef_all$population) <- c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
                                     "Upper Mainstem","White-Donjek")

coef_all$var <- as.factor(coef_all$var)
levels(coef_all$var) <- c("Ice out date","Migration temperature","Rearing degree days","Rearing precipitation",
                              "Snowpack","Spawning temperature","Spawning precipitation","Summer SST",
                              "Winter SST","Pink salmon","Chum salmon","FW winter temperature")

#add stage

coef.popn <- coef_all %>% 
  mutate(stage=if_else(var=="Chum salmon"|var=="Pink salmon"|var=="Summer SST"|
                         var=="Winter SST","Marine",
                       if_else(var=="Ice out date","Outmigration",
                               if_else(var=="Rearing degree days"|var=="Rearing precipitation",
                                       "Freshwater juvenile",
                                       "Spawning and incubation"))))


Popn_effects <- coef.popn %>% 
  mutate(var=fct_relevel(var,"Migration temperature","Spawning temperature","Spawning precipitation",
                         "Snowpack","FW winter temperature","Rearing degree days",
                         "Rearing precipitation","Ice out date","Summer SST",
                         "Winter SST","Pink salmon","Chum salmon")) %>% 
  mutate(stage=fct_relevel(stage,"Spawning and incubation","Freshwater juvenile",
                           "Outmigration","Marine")) %>% 
  mutate(population=fct_relevel(population,"White-Donjek","Upper Mainstem","Teslin",
                                "Stewart","Pelly","Middle Mainstem",
                                "Lower Mainstem","Carmacks")) %>% 
  ggplot(aes(x=population, y=value,fill=stage)) +
  geom_hline(yintercept = 0, colour='red')+
  stat_eye()+
  theme_bw() +
  scale_fill_manual(values = c("#FE9000","#FFDD4A","#5ADBFF","#3C6997"))+
  coord_flip() +
  stat_summary(fun="q.95", colour="black", geom="line", lwd=0.5,show.legend=FALSE) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=1.25,show.legend=FALSE) +
  stat_summary(fun="median", colour="black", size=1.75, geom="point", pch=21, show.legend=FALSE) +
  stat_summary(fun="median", colour="gray", size=1.5, geom="point", pch=19,show.legend=FALSE) +
  theme(legend.position='top') +
  theme(legend.title=element_blank())+
  facet_wrap(~var,ncol=3)+
  xlab("Population") + ylab("Change in log (recruits/spawner)") +
  theme(text = element_text(size=20),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=15,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=0)))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.text=element_text(size=20))+scale_y_continuous(limits = c(-0.45, 0.5),
                                                              labels=c(-0.4,-0.2,0,0.2,0.4),
                                                              breaks=c(-0.4,-0.2,0,0.2,0.4))

Popn_effects
#saved dimensions were 1200x1500
