
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
Migration_temp_t0 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Weekly_max_migration_temp_t0.csv"))
Migration_temp_returns <- read.csv(file.path(dir.data,"/Environmental data/Processed/Migration_temp_returns.csv"))
rearing_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_temp.csv"))
rearing_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_prcp.csv"))
annual_snowpack <- read.csv(file.path(dir.data,"/Environmental data/Processed/annual_snowpack.csv"))
spawning_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_temp2.csv"))
spawning_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/spawning_prcp.csv"))
SST_summer <- read.csv(file.path(dir.data,"/Environmental data/Processed/SST_summer.csv"))
SST_winter <- read.csv(file.path(dir.data,"/Environmental data/Processed/SST_winter_new.csv"))
Threshold17 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Threshold17_runday.csv"))
Threshold17_count <- read.csv(file.path(dir.data,"/Environmental data/Processed/Threshold17_count.csv"))
Discharge <- read.csv(file.path(dir.data,"/Environmental data/Processed/Pilotdischarge.csv"))


# Define covariate names
names.covars <- c("Ice_out","Migration_temp_t0","rearing_temp","rearing_prcp",
                  "annual_snowpack","spawning_temp","spawning_prcp","SST_summer",
                  "SST_winter")



n.covars <- length(names.covars)

covars <- array(data=NA,dim=c(n.pops, max(n.years), n.covars))


##put covariates into covars array for model

library(abind)

covars <- abind(Ice_out,Migration_temp_t0,rearing_temp,rearing_prcp,
                annual_snowpack,spawning_temp,spawning_prcp,SST_summer,
                SST_winter,along=3)

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
    exp.alpha[p] ~ dunif(0,50)
    alpha[p] <- log(exp.alpha[p])
    beta[p] ~ dnorm(0,pow(0.1,-2));T(0,100)
    sigma.oe[p] ~ dnorm(0, pow(1,-2));T(1e-3,100)
    
    #Covariate Effects
    for(c in 1:n.covars) {
      coef[p,c] ~ dnorm(mu.coef[c],pow(sigma.coef[c],-2))
    }
    # AR-1 Coeff
    pre.phi[p] ~ dnorm(0,2)
    #scale to (1,1)
    phi[p] <-  2*exp(pre.phi[p])/(1+exp(pre.phi[p])) - 1
      # phi[p] ~ dunif(-0.99, 0.99)
  }#next p
  

  
  #PREDICTIONS
  for(p in 1:n.pops) {
    # First Year
    for(c in 1:n.covars) {
        cov.eff[p,1,c] <- coef[p,c]*covars[p,1,c]
    }
     pred.rec.1[p,1] <- Sobs[p,1]*exp(alpha[p] - Sobs[p,1]*beta[p] + sum(cov.eff[p,1,1:n.covars]))
    pred.rec.2[p,1] <- pred.rec.1[p,1]
     resid[p,1] <- exp(lnRobs[p,1]) - pred.rec.1[p,1]
     log.pred.rec.2[p,1] <- log(pred.rec.2[p,1]) 
    
   #Subsequent Years
     for(y in 2:n.years[p]) {
      
       for(c in 1:n.covars) {
        cov.eff[p,y,c] <- coef[p,c]*covars[p,y,c]
      }
      
      pred.rec.1[p,y] <- Sobs[p,y]*exp(alpha[p] - Sobs[p,y]*beta[p] + sum(cov.eff[p,y,1:n.covars]))
      resid[p,y] <- exp(lnRobs[p,y]) - pred.rec.1[p,y]
      
      pred.rec.2[p,y] <- pred.rec.1[p,y]+resid[p,y-1]*phi[p]
      log.pred.rec.2[p,y] <- log(pred.rec.2[p,y]) 
      
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

#without autocorrelation variable
# jags_model = function(){
# 
#   #Priors
# 
#     for(c in 1:n.covars) {
#       mu.coef[c] ~ dnorm(0, pow(25,-2))
#       sigma.coef[c] ~ dnorm(0, pow(5,-2));T(1e-3,100)
#       dist.coef[c] ~ dnorm(mu.coef[c], pow(sigma.coef[c],-2))
#     }#next c
# 
#   #population-specific
# 
#   for(p in 1:n.pops) {
#     exp.alpha[p] ~ dunif(0,25)
#     alpha[p] <- log(exp.alpha[p])
#     beta[p] ~ dnorm(0,pow(0.1,-2));T(0,100)
#     sigma.oe[p] ~ dnorm(0, pow(1,-2));T(1e-3,100)
# 
#     #Covariate Effects
#     for(c in 1:n.covars) {
#       coef[p,c] ~ dnorm(mu.coef[c],pow(sigma.coef[c],-2))
#     }
#   }#next p
# 
#     #PREDICTIONS
#   for(p in 1:n.pops) {
#     for(y in 1:n.years[p]) {
# 
#       for(c in 1:n.covars) {
#         cov.eff[p,y,c] <- coef[p,c]*covars[p,y,c]
#       }
# 
#       pred.rec[p,y] <- Sobs[p,y]*exp(alpha[p] - Sobs[p,y]*beta[p] + sum(cov.eff[p,y,1:n.covars]))
# 
#         }#next y
#     }#next p
# 
#   #LIKELIHOODS
#   for(p in 1:n.pops) {
#     for(y in 1:n.years[p]) {
#       lnRobs[p,y] ~ dnorm(log(pred.rec[p,y]), pow(sigma.oe[p],-2))
# 
#     }#next y
#   }#next p
# }


#### Specify initial values ####

jags_inits = function() {
  exp.alpha <- runif(n.pops,1,2)
  #beta <- runif(n.pops,0,0.001)
  beta <- runif(n.pops,0,0.3) #try increasing this based on results from popn diversity paper?
  sigma.oe <- runif(n.pops,0.1,1) #seems good based on prior data
  mu.coef <- rnorm(n.covars, 0, 1)
  sigma.coef <- rgamma(n.covars, 1, 1)
  
  Return <- list(exp.alpha=exp.alpha, beta=beta, sigma.oe=sigma.oe,
                 mu.coef=mu.coef, sigma.coef=sigma.coef)
  return(Return)
}


#### Set nodes to monitor ####

jags_params = c("alpha","exp.alpha", "beta", "sigma.oe","mu.coef","sigma.coef",
                "coef","cov.eff","dist.coef","pred.rec.1","pred.rec.2","log.resid","phi","pre.phi")


#without autocorrl
# jags_params = c("alpha","exp.alpha", "beta", "sigma.oe","mu.coef","sigma.coef",
#                 "coef","cov.eff","dist.coef","pred.rec")


##### Run the model with Jags #####

#fast
out <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=50000, 
  n.burnin=10000)  


out <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=500000, 
  n.burnin=20000)  

#Save
saveRDS(out, file=file.path(dir.output,"out.rds"))

out.mcmc <- as.mcmc(out)

#assess autocorrelation before correction

mean.pred.rec <- as.data.frame(out$BUGSoutput$mean$pred.rec) %>% 
  mutate(population=c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
        "Upper Lakes and Mainstem","White-Donjek")) %>% gather(year,pred.rec,1:28) %>% 
  mutate(year=rep(1985:2012,each=8))

mean.obs.rec <- brood_table %>% select(population,BroodYear,R_med) %>% 
  rename(year="BroodYear")

calc_resid <- left_join(mean.pred.rec,mean.obs.rec) %>% 
  mutate(resid=R_med-pred.rec,log.resid=log(R_med)-log(pred.rec)) %>% select(-R_med,-pred.rec) %>% 
  spread(population,resid) %>% select(-year)

par(mfrow=c(2,1))
acf(calc_resid[,4],main="Interpret the ARMA Order")
pacf(calc_resid[,4],main="")

#explore with arima models
mod1 <- arima(calc_resid[,2],order=c(1,0,0))
mod1

#Residuals Model 1
mod1.residuals=mod1$residuals
par(mfrow=c(2,1))
acf(mod1.residuals, main="ACF of Residuals")
pacf(mod1.residuals, main = "PACF of Residuals")


initial.acf <- acf(calc_resid[,2],plot=TRUE)
initial.acf$acf[2]

#populations 1, 3, 5, 6 have no autocorrelation
#population 2 has a sig lag at time step 4 for the pacf plot. An arima of (1,0,0) works but (2,0,0) works better
#pop 4, 7, and 8 have sig lag(s) and an arima of (1,0,0) works

#assess autocorrelation after adding phi

model_resids <- as.data.frame.table(out$BUGSoutput$sims.list$log.resid)
names(model_resids) <- c("reps","population","year","log.resid")

model_resids2 <- model_resids %>% group_by(population,year) %>% 
  summarise(mean_log.resid=mean(log.resid)) %>% spread(population,mean_log.resid) %>% select(-year)

par(mfrow=c(2,1))
acf(model_resids2[,4],main="Interpret the ARMA Order")
pacf(model_resids2[,4],main="")

test.acf <- acf(model_resids2[,8],plot=TRUE)
test.acf$acf[2]

#Write Output File for Diagnostics
write.csv(out$BUGSoutput$summary, file=file.path(dir.figs,"out_summary_Jun9_prephi.csv"))


##### STEP 7: CONVERGENCE DIAGNOSTICS #####
## Code adapted from Intro to Bayesian Stats with JAGS course by Ben Staton ##
diag_p = c("alpha", "beta", "mu.coef","sigma.coef", "coef","cov.eff","dist.coef","sigma.oe","pred.rec")
diag_p = c("cov.eff")


# view convergence diagnostic summaries for nodes with priors
diag <- t(post_summ(out.mcmc, diag_p, Rhat = T, neff = T)[c("Rhat", "neff"),])
diagb <- t(post_summ(out.mcmc, diag_p, Rhat = T, neff = T)[c("Rhat", "neff"),])

# view diagnostic plots
diag_plots(as.mcmc(out), diag_p, ext_device = T)



# Visualize results -------------------------------------------------------


pops2=c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
        "Upper Mainstem","White-Donjek")

labels2=c("Annual snowpack","Winter SST","Spawning temperature","Summer SST",
          "Ice out date","Migration temperature","Rearing temperature","Spawning precip","Rearing precip")


#rearrange variables

mu.coef.vector <- out$BUGSoutput$sims.list$mu.coef
mu.coef.vector.new <- mu.coef.vector[,c(5,9,6,8,1,2,3,7,4)]

#Hyper means summarized in order of median effect
par(mfcol=c(1,1), mar=c(2,12,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(out$BUGSoutput$sims.list$mu.coef,
            labels=labels2,reorder=TRUE, quantiles=list(0.025,0.25,0.75,0.975), 
            style='gray', col='blue',cex=1.1)
caterpoints(apply(mu.coef.vector.new,2,median), pch=21, col='red', bg='orange')
abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5,at=0.7)
mtext('Covariate', side=2, outer=TRUE, font=2, line=0.5)

#739x503 saved dimensions

pdf(file="Plots/Results_Jun8_prephi.pdf",width=9,height=6)

#Hyper means summarized
par(mfcol=c(1,1), mar=c(2,12,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(out$BUGSoutput$sims.list$mu.coef,
            labels=names.covars,reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), 
            style='plain', col='blue',cex=1.1)
caterpoints(apply(out$BUGSoutput$sims.list$mu.coef,2,median), pch=21, col='red', bg='orange')
abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5)
mtext('Covariate', side=2, outer=TRUE, font=2, line=0.5)


#Hyper means for covariates
par(mfcol=c(2,3), mar=c(5,0,1,0), oma=c(1,1,3,1))
c <- 1
for(c in 1:n.covars) {
  plotPost(out$BUGSoutput$sims.list$mu.coef[,c], showCurve=TRUE, main='', xlab=names.covars[c],
           xlim=c(-0.5,0.5), rope=0)
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0,alpha=0.5))
}

#Population specific covariate "effects"
par(mfcol=c(2,3), mar=c(2,7,3,1), oma=c(2,2,1,1))
c <- 1
for(c in 1:n.covars) {
  caterplot(out$BUGSoutput$sims.list$coef[,,c],
            labels=pops2, cex.labels=TRUE, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), style='plain', col='blue')
  mtext(names.covars[c], side=3, outer=FALSE, line=1)
  caterpoints(apply(out$BUGSoutput$sims.list$coef[,,c],2,median), reorder=FALSE, pch=21, col='red', bg='orange')
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
}
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5)
mtext('Population', side=2, outer=TRUE, font=2, line=0.5)
dev.off() 

#Population specific covariate "effects"
par(mfcol=c(1,1), mar=c(2,10,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(out$BUGSoutput$sims.list$coef[,,9],
            labels=pops2, cex=1.1, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), style='plain', col='blue')
  mtext("Winter SST", side=3, outer=FALSE, line=1,cex=1.1)
  caterpoints(apply(out$BUGSoutput$sims.list$coef[,,9],2,median), reorder=FALSE, pch=21, col='red', bg='orange')
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5,at=0.6)
mtext('Population', side=2, outer=TRUE, font=2, line=0.5)

#save dims 1000 x 730


#Population specific covariate "effects" saved by migration temp order

mig_temp <- out$BUGSoutput$sims.list$coef[,,2]
mig_temp_new <- mig_temp[,c(2,8,5,4,6,7,1,3)]
pops3=c("Lower Mainstem","White-Donjek","Stewart","Pelly","Teslin","Upper Mainstem","Carmacks","Middle Mainstem")
  
par(mfcol=c(1,1), mar=c(2,10,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(mig_temp_new,
            labels=pops3, cex=1.1, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), style='plain', col='blue')
  mtext("Migration temperature", side=3, outer=FALSE, line=1,cex=1.1)
  caterpoints(apply(mig_temp_new,2,median), reorder=FALSE, pch=21, col='red', bg='orange')
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Coefficient (Effect)', side=1, outer=TRUE, font=2, line=0.5,at=0.6)
#740x528



# Plot Fits =========

Rec <- exp(ln.Recruits)

pdf(file="Plots/Fitted_plot_noauto.pdf",width=9,height=6)
par(mfrow=c(2,2), mar=c(2,2,2,0), oma=c(2,2,1,1))
p <- 1 
for(p in 1:n.pops) {
  pred <- apply(out$BUGSoutput$sims.list$pred.rec[,p,1:n.years[p]],2, 
                quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  log.pred <- (log(pred))
  
  ylim <- c(0,max((ln.Recruits[p,]), log.pred, na.rm=TRUE))
  
  temp.yrs <- years[p,1:n.years[p]]
  
  plot(x=temp.yrs, y=(ln.Recruits[p,1:n.years[p]]), pch=21, bg=rgb(0,0,1,alpha=0.5), ylim=ylim)
  
  polygon(x=c(temp.yrs,rev(temp.yrs)), y=c(log.pred[1,],rev(log.pred[5,])), col=rgb(1,0,0, alpha=0.2),
          border=FALSE)
  polygon(x=c(temp.yrs,rev(temp.yrs)), y=c(log.pred[2,],rev(log.pred[4,])), col=rgb(1,0,0, alpha=0.2),
          border=FALSE)
  
  lines(x=temp.yrs, y=log.pred[3,], col='red')
  mtext(pops[p], side=3, line=0.25, font=2)
  if(p %in% c(1,5,9,13)) {
    mtext('Log Recruitment', side=2, font=2, outer=TRUE, line=0.5)
    mtext('Year', side=1, font=2, outer=TRUE, line=0.5)
  }
}
dev.off()

#quantitative inference

#(exp(median_coef)-1)*100 is percent change in recruitment for 1SD of that variable
#(exp(median_coef*2)-1)*100 is percent change in recruitment for 2SD of that variable
