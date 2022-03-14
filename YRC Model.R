
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

Spawners <- brood_table %>% select(population,BroodYear,S_med) %>% 
 spread(BroodYear,S_med) %>% select(-population)

ln.Recruits <- brood_table %>% select(population,BroodYear,R_med) %>% 
  mutate(lnR_med=log(R_med)) %>% select(-R_med) %>% 
  spread(BroodYear,lnR_med) %>% select(-population)

#popn order is Carmacks, Lower Main, Middle Main, Pelly, Stewart, Teslin, Upper Lakes, White-Donjek


# Merge covariate data ----------------------------------------------------

Ice_out <- read.csv(file.path(dir.data,"/Environmental data/Processed/Ice_out.csv"))
Migration_temp_t0 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Migration_temp_t0.csv"))
Migration_temp_returns <- read.csv(file.path(dir.data,"/Environmental data/Processed/Migration_temp_returns.csv"))
rearing_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_temp.csv"))
rearing_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_prcp.csv"))
annual_snowpack <- read.csv(file.path(dir.data,"/Environmental data/Processed/annual_snowpack.csv"))
spawning_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_temp.csv"))
spawning_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_prcp.csv"))

# Define covariate names
names.covars <- c("Ice_out","Migration_temp_t0","Migration_temp_returns","rearing_temp","rearing_prcp",
                  "annual_snowpack","spawning_temp","spawning_prcp")


n.covars <- length(names.covars)

covars <- array(data=NA,dim=c(n.pops, max(n.years), n.covars))


##put covariates into covars array for model

library(abind)

covars <- abind(Ice_out,Migration_temp_t0,Migration_temp_returns,rearing_temp,rearing_prcp,
                annual_snowpack,spawning_temp,spawning_prcp,along=3)

print(covars)


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
  
#Reminder here to try testing with and without autocorrelation term 

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
  sigma.oe <- runif(n.pops,0.1,1)
  mu.coef <- rnorm(n.covars, 0, 1)
  sigma.coef <- rgamma(n.covars, 1, 1)
  
  Return <- list(exp.alpha=exp.alpha, beta=beta, sigma.oe=sigma.oe,
                 mu.coef=mu.coef, sigma.coef=sigma.coef)
  return(Return)
}


#### Set nodes to monitor ####

jags_params = c("alpha","exp.alpha", "beta", "sigma.oe","mu.coef","sigma.coef",
                "coef","cov.eff","dist.coef","pred.rec")

##### Specify MCMC Dimensions #####

#Note to check convergence diagnostics and increase iterations and/or burn in if needed in the future

jags_dims = c(
  ni = 1e4,  # number of post-burn-in samples per chain
  nb = 1e4,  # number of burn-in samples
  nt = 200,     # thinning rate
  nc = 3      # number of chains
)


##### Run the model with Jags #####

out <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=50000, 
  n.burnin=20000)  


#when iterations lowered get this error: Error in checkForRemoteErrors(val) : 
#3 nodes produced errors; first error: n.iter must be a positive integer

#Save
saveRDS(out, file=file.path(dir.output,"out.rds"))

#Write Output File for Diagnostics
write.csv(out$BUGSoutput$summary, file=file.path(dir.figs,"out_summary.csv"))

##### STEP 7: CONVERGENCE DIAGNOSTICS #####
## Code adapted from Intro to Bayesian Stats with JAGS course by Ben Staton ##
diag_p = c("alpha", "beta", "sigma.coef", "sigma.oe","pred.rec")

# view convergence diagnostic summaries for nodes with priors
t(post_summ(as.mcmc(out), diag_p, Rhat = T, neff = T)[c("Rhat", "neff"),])

# view diagnostic plots
diag_plots(as.mcmc(out), diag_p, ext_device = T)

#Visualize results

pdf(file="Plots/Results.pdf",width=9,height=6)

#Hyper means for covariates
par(mfcol=c(2,3), mar=c(5,0,1,0), oma=c(1,1,3,1))
c <- 1
for(c in 1:n.covars) {
  plotPost(out$BUGSoutput$sims.list$mu.coef[,c], showCurve=TRUE, main='', xlab=names.covars[c],
           xlim=c(-0.5,0.5), rope=0)
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0,alpha=0.5))
}

#Population specific covariate "effects"
par(mfcol=c(2,3), mar=c(2,5,3,1), oma=c(2,2,1,1))
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

dev.off() 
