
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
rearing_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_GDD.csv"))
rearing_prcp <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_prcp_max.csv"))
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
names.covars <- c("Ice_out","Migration_temp_t0","rearing_temp","rearing_prcp",
                  "annual_snowpack","spawning_temp","spawning_prcp","SST_summer",
                  "SST_winter","Pink_comp","Chum_comp","First_winter_FW")


# names.covars <- c("Migration_temp_t0","rearing_prcp",
#                   "annual_snowpack","spawning_temp",
#                   "SST_winter","Pink_comp")

n.covars <- length(names.covars)

covars <- array(data=NA,dim=c(n.pops, max(n.years), n.covars))



##put covariates into covars array for model

library(abind)

covars <- abind(Ice_out,Migration_temp_t0,rearing_temp,rearing_prcp,
annual_snowpack,spawning_temp,spawning_prcp,SST_summer,
SST_winter,Pink_comp,Chum_comp,First_winter_FW,along=3)


# covars <- abind(Migration_temp_t0,rearing_prcp,
#                 annual_snowpack,spawning_temp,
#                 SST_winter,Pink_comp,along=3)

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





# #model without covariates
# jags_model = function(){
# 
#   #Priors
# 
#   for(p in 1:n.pops) {
#     exp.alpha[p] ~ dunif(0,50)
#     alpha[p] <- log(exp.alpha[p])
#     beta[p] ~ dnorm(0,pow(0.1,-2));T(0,100)
#     sigma.oe[p] ~ dnorm(0, pow(1,-2));T(1e-3,100)
#      phi[p] ~ dunif(-0.99, 0.99)
# 
#   }#next p
# 
# 
# 
#   #PREDICTIONS
#   for(p in 1:n.pops) {
#     # First Year
# 
#     pred.rec.1[p,1] <- Sobs[p,1]*exp(alpha[p] - Sobs[p,1]*beta[p])
# 
#     log.pred.rec.1[p,1] <- log(pred.rec.1[p,1])
#     log.resid[p,1] <- lnRobs[p,1] - log.pred.rec.1[p,1]
#     log.pred.rec.2[p,1] <- log.pred.rec.1[p,1]
#     log.resid.2[p,1] <- log.resid[p,1]
# 
#    #Subsequent Years
#      for(y in 2:n.years[p]) {
# 
#     pred.rec.1[p,y] <- Sobs[p,y]*exp(alpha[p] - Sobs[p,y]*beta[p])
# 
#      log.pred.rec.1[p,y] <- log(pred.rec.1[p,y])
#      log.resid[p,y] <- lnRobs[p,y] - log.pred.rec.1[p,y]
#      log.pred.rec.2[p,y] <- log.pred.rec.1[p,y]+log.resid[p,y-1]*phi[p]
#      log.resid.2[p,y] <- lnRobs[p,y] - log.pred.rec.2[p,y]
# 
#         }#next y
#     }#next p
# 
#   #LIKELIHOODS
#   for(p in 1:n.pops) {
# 
#     #First Year
#     lnRobs[p,1] ~ dnorm(log.pred.rec.2[p,1], pow(sigma.oe[p],-2))
# 
#     for(y in 2:n.years[p]) {
#     lnRobs[p,y] ~ dnorm(log.pred.rec.2[p,y], pow(sigma.oe[p],-2))
# 
#     }#next y
#   }#next p
#   
# #for creating population-specific predicted SR relationships
#   
# for(p in 1:n.pops) {
#   for(y in 1:n.years[p]) {
#   for (i in 1:nSpred) {
#     # predict recruitment and sustained yield at each level of escapement
#     Rpred[p,i] <- Spred[p,i]*exp(alpha[p] - Spred[p,i]*beta[p])
#     
#   }
#   }
# }
# }


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

#fast
out <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=50000, 
  n.burnin=10000)

#increase number of simulations?

out <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=250000, 
  n.burnin=15000)  

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

#Write Output File for Diagnostics
write.csv(out_orig$BUGSoutput$summary, file=file.path(dir.figs,"out_Pink.csv"))
write.csv(out_nocomp$BUGSoutput$summary, file=file.path(dir.figs,"out_nocomp.csv"))
write.csv(out_spawnair$BUGSoutput$summary, file=file.path(dir.figs,"out_spawnair.csv"))

#assess autocorrelation

log.mean.pred.rec <- as.data.frame(out$BUGSoutput$mean$log.pred.rec.1) %>% 
  mutate(population=c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
        "Upper Lakes and Mainstem","White-Donjek")) %>% gather(year,log.pred.rec,1:28) %>% 
  mutate(year=rep(1985:2012,each=8))

log.mean.obs.rec <- brood_table %>% select(population,BroodYear,R_med) %>% 
  rename(year="BroodYear") %>% mutate(log_Robs=log(R_med))

calc_resid <- left_join(log.mean.pred.rec,log.mean.obs.rec) %>% 
  mutate(log.resid=log_Robs-log.pred.rec) %>% select(-log_Robs,-log.pred.rec,-R_med) %>% 
  spread(population,log.resid) %>% select(-year)

par(mfrow=c(2,1))
acf(calc_resid[,8],main="Interpret the ARMA Order")
pacf(calc_resid[,8],main="")

#explore with arima models
mod1 <- arima(calc_resid[,8],order=c(1,0,0))
mod1

#Residuals Model 1
mod1.residuals=mod1$residuals
par(mfrow=c(2,1))
acf(mod1.residuals, main="ACF of Residuals")
pacf(mod1.residuals, main = "PACF of Residuals")


initial.acf <- acf(calc_resid[,2],plot=TRUE)
initial.acf$acf[2]


log.resid <- as.data.frame(out$BUGSoutput$mean$log.resid) %>% 
  mutate(population=c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
        "Upper Lakes and Mainstem","White-Donjek")) %>% gather(year,log.resid,1:28) %>% 
  mutate(year=rep(1985:2012,each=8))



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

pdf(file="Plots/Fitted_plot_nocovars_July5.pdf",width=9,height=6)
par(mfrow=c(2,2), mar=c(2,2,2,0), oma=c(2,2,1,1))
p <- 1 
for(p in 1:n.pops) {
  log.pred <- apply(out_null$BUGSoutput$sims.list$log.pred.rec.2[,p,1:n.years[p]],2, 
                quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  
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


# Plot parameters ---------------------------------------------------------

library(bayesplot)

alpha.list <- out_null$BUGSoutput$sims.list$alpha
colnames(alpha.list) <- pops
alpha_plot <- mcmc_areas(alpha.list) +
  ggtitle('Ricker Alpha')
plot(alpha_plot)

#Beta

library(BayesTools)

#plot prior
p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1000))
plot(p1)

beta.list <- out$BUGSoutput$sims.list$beta
colnames(beta.list) <- pops
beta_plot <- mcmc_areas(beta.list) +
  ggtitle('Ricker Beta')
plot(beta_plot)


#plot prior
p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
plot(p1)

sigmaR.list <- out$BUGSoutput$sims.list$sigma.oe
colnames(sigmaR.list) <- pops
sigmaR_plot <- mcmc_areas(sigmaR.list) +
  ggtitle('Ricker Recruitment SD')
plot(sigmaR_plot)

#plot prior
p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 5))
plot(p1)

sigma_coef.list <- out$BUGSoutput$sims.list$sigma.coef
colnames(sigma_coef.list) <- names.covars
sigma_coef_plot <- mcmc_areas(sigma_coef.list) +
  ggtitle('Ricker Coef SD')
plot(sigma_coef_plot)

# Infer percent change in recruits ----------------------------------------

#(exp(median_coef)-1)*100 is percent change in recruitment for 1SD of that variable
#(exp(median_coef*2)-1)*100 is percent change in recruitment for 2SD of that variable



# Infer change in number of recruits --------------------------------------

#Processing
alpha <- as.data.frame(out$BUGSoutput$sims.list$alpha)

alpha_sample <- alpha %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,alpha,1:8)

beta <- as.data.frame(out$BUGSoutput$sims.list$beta)

beta_sample <- beta %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,beta,1:8)
  
mig_temp <- as.data.frame(out$BUGSoutput$sims.list$coef[,,2])

mig_temp_sample <- mig_temp %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,migtemp,1:8)

ice_out <- as.data.frame(out$BUGSoutput$sims.list$coef[,,1])

ice_out_sample <- ice_out %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,ice_out,1:8)

rear_temp <- as.data.frame(out$BUGSoutput$sims.list$coef[,,3])

rear_temp_sample <- rear_temp %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,reartemp,1:8)

rear_precip <- as.data.frame(out$BUGSoutput$sims.list$coef[,,4])

rear_precip_sample <- rear_precip %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,rear_precip,1:8)

snow <- as.data.frame(out$BUGSoutput$sims.list$coef[,,5])

snow_sample <- snow %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,snow,1:8)

spawn_temp <- as.data.frame(out$BUGSoutput$sims.list$coef[,,6])

spawn_temp_sample <- spawn_temp %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,spawntemp,1:8)

spawn_precip <- as.data.frame(out$BUGSoutput$sims.list$coef[,,7])

spawn_precip_sample <- spawn_precip %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,spawnprecip,1:8)

SST_summ <- as.data.frame(out$BUGSoutput$sims.list$coef[,,8])

SST_summ_sample <- SST_summ %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,SSTsummer,1:8)

SST_wint <- as.data.frame(out$BUGSoutput$sims.list$coef[,,9])

SST_wint_sample <- SST_wint %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,SSTwinter,1:8)

Pink <- as.data.frame(out$BUGSoutput$sims.list$coef[,,10])

Pink_sample <- Pink %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,Pink,1:8)

Chum <- as.data.frame(out$BUGSoutput$sims.list$coef[,,11])

Chum_sample <- Chum %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,Chum,1:8)

Winter <- as.data.frame(out$BUGSoutput$sims.list$coef[,,12])

Winter_sample <- Winter %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,WinterFW,1:8)

All_popns <- data.frame(alpha_sample,beta_sample,mig_temp_sample,
                     ice_out_sample,rear_temp_sample,rear_precip_sample,
                     snow_sample,spawn_temp_sample,spawn_precip_sample,
                     SST_summ_sample,SST_wint_sample,Pink_sample,
                     Chum_sample,Winter_sample) %>% select(1:2,4,6,8,10,12,14,16,18,20,22,24,
                                                           26,28)


All_popns$population <- as.factor(All_popns$population)
levels(All_popns$population)<- list("Lower Mainstem"="V2","White-Donjek"="V8","Middle Mainstem"="V3","Upper Lakes and Mainstem"="V7",
                                          Carmacks="V1",Teslin="V6",Stewart="V5",Pelly="V4")

All_popns_median <- All_popns

#create spawner estimates
Spred = seq(0, 32, 0.1)
Spred <- t(as.data.frame(Spred))
prefix <- "Spred"
suffix <- seq(1:321)
colnames(Spred) <- paste(prefix,suffix,sep="")
Spred <- Spred[rep(seq_len(nrow(Spred)),each=8000),]
rownames(Spred) <- NULL

All_popns <- data.frame(All_popns,Spred)


All_popns <- All_popns %>% gather(S_level,Spred,16:336)

All_popns <- All_popns %>% 
  mutate(Rpred_migtemp = Spred*exp(alpha-beta*Spred+migtemp*1), #change in 1SD for migration temp
         Rpred_iceout = Spred*exp(alpha-beta*Spred+ice_out*1),
         Rpred_reartemp = Spred*exp(alpha-beta*Spred+reartemp*1),
         Rpred_rearprcp = Spred*exp(alpha-beta*Spred+rear_precip*1),
         Rpred_snow = Spred*exp(alpha-beta*Spred+snow*1),
         Rpred_spawntemp = Spred*exp(alpha-beta*Spred+spawntemp*1),
         Rpred_spawnprcp = Spred*exp(alpha-beta*Spred+spawnprecip*1),
         Rpred_SSTsummer = Spred*exp(alpha-beta*Spred+SSTsummer*1),
         Rpred_SSTwinter = Spred*exp(alpha-beta*Spred+SSTwinter*1),
         Rpred_Pink = Spred*exp(alpha-beta*Spred+Pink*1),
         Rpred_Chum = Spred*exp(alpha-beta*Spred+Chum*1),
         Rpred_WinterFW = Spred*exp(alpha-beta*Spred+WinterFW*1),
        Rpred_null = Spred*exp(alpha-beta*Spred),
        Migtemp_diff=(Rpred_migtemp-Rpred_null)*1000, #multiplied by 1000 to give true numbers
        Iceout_diff=(Rpred_iceout-Rpred_null)*1000,
        Reartemp_diff=(Rpred_reartemp-Rpred_null)*1000,
        Rearprcp_diff=(Rpred_rearprcp-Rpred_null)*1000,
        Snow_diff=(Rpred_snow-Rpred_null)*1000,
        Spawntemp_diff=(Rpred_spawntemp-Rpred_null)*1000,
        Spawnprcp_diff=(Rpred_spawnprcp-Rpred_null)*1000,
        SSTsummer_diff=(Rpred_SSTsummer-Rpred_null)*1000,
        SSTwinter_diff=(Rpred_SSTwinter-Rpred_null)*1000,
        Pink_diff=(Rpred_Pink-Rpred_null)*1000,
        Chum_diff=(Rpred_Chum-Rpred_null)*1000,
        WinterFW_diff=(Rpred_WinterFW-Rpred_null)*1000,
        Spawners=Spred*1000)

All_popns_summary <- All_popns %>% select(1,31:43) %>% 
  gather(Covariate,Value,2:13) %>% 
  group_by(population,Spawners,Covariate) %>% 
  summarise(Med=median(Value),
            q5=quantile(Value,0.05),
            q95=quantile(Value,0.95))


#visualize
All_popns_plot <- All_popns_summary %>% filter(Covariate=="Migtemp_diff") %>% 
  ggplot()+
  geom_ribbon(aes(x=Spawners,ymin=q5,ymax=q95),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_line(aes(x = Spawners, y = Med), color="black", size = 0.75) +
  facet_wrap(~population)+
  xlab("Spawners") +
  ylab("Change in Recruits") +
  theme_bw() +
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#Maximum change in recruitment summary

Aggregate_change_neg <- All_popns_summary %>% filter(Med<0) %>% 
  group_by(population,Covariate) %>% 
  slice_min(Med) #using max change based on median values

Aggregate_change_pos <- All_popns_summary %>% filter(Med>0) %>% 
  group_by(population,Covariate) %>% 
  slice_max(Med) #using max change based on median values


Aggregate_change <- bind_rows(Aggregate_change_pos,Aggregate_change_neg)
  
Aggregate_change_summary <-  Aggregate_change %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") %>% #removing vars with no relationship
  group_by(Covariate) %>% 
  summarise(sum_spawners=sum(Spawners),
            sum_maxR_change=sum(Med))

Aggregate_change_summary$Covariate <- as.factor(Aggregate_change_summary$Covariate)

Max_change_plot <- Aggregate_change_summary %>%
  mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot(aes(y=sum_maxR_change,x=Covariate))+
  geom_col()+ylab("Max change in recruits")+xlab("Covariate")+theme_bw()+
  scale_x_discrete(labels=c("Iceout_diff"="Ice\nout","Migtemp_diff"="Migration\ntemp",
                            "Pink_diff"="Pink\nSalmon","Rearprcp_diff"="Rearing\nprecip",
                            "Reartemp_diff"="Rearing\ntemp",
                               "Snow_diff"="Snow","Spawnprcp_diff"="Spawning\nprecip",
                               "Spawntemp_diff"="Spawning\ntemp",
                            "SSTwinter_diff"="Winter\nSST",
                            "WinterFW_diff"="Winter\nfreshwater"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
  

#calculate change in recruitment at median popn-specific spawner values instead

summary_stats <- brood_table %>% 
  group_by(population) %>%
  summarise(Spawner_med=median(S_med),
            Recruit_med=median(R_med))

Spawner_med <- summary_stats %>% select(-Recruit_med)

Spawner_med$population <- as.factor(Spawner_med$population)
  
All_popns_median <- left_join(All_popns_median,Spawner_med)

All_popns_median2 <- All_popns_median %>% 
  mutate(Rpred_migtemp = Spawner_med*exp(alpha-beta*Spawner_med+migtemp*1), #change in 1SD for migration temp
         Rpred_iceout = Spawner_med*exp(alpha-beta*Spawner_med+ice_out*1),
         Rpred_reartemp = Spawner_med*exp(alpha-beta*Spawner_med+reartemp*1),
         Rpred_rearprcp = Spawner_med*exp(alpha-beta*Spawner_med+rear_precip*1),
         Rpred_snow = Spawner_med*exp(alpha-beta*Spawner_med+snow*1),
         Rpred_spawntemp = Spawner_med*exp(alpha-beta*Spawner_med+spawntemp*1),
         Rpred_spawnprcp = Spawner_med*exp(alpha-beta*Spawner_med+spawnprecip*1),
         Rpred_SSTsummer = Spawner_med*exp(alpha-beta*Spawner_med+SSTsummer*1),
         Rpred_SSTwinter = Spawner_med*exp(alpha-beta*Spawner_med+SSTwinter*1),
         Rpred_Pink = Spawner_med*exp(alpha-beta*Spawner_med+Pink*1),
         Rpred_Chum = Spawner_med*exp(alpha-beta*Spawner_med+Chum*1),
         Rpred_WinterFW = Spawner_med*exp(alpha-beta*Spawner_med+WinterFW*1),
        Rpred_null = Spawner_med*exp(alpha-beta*Spawner_med),
        Migtemp_diff=(Rpred_migtemp-Rpred_null)*1000, #multiplied by 1000 to give true numbers
        Iceout_diff=(Rpred_iceout-Rpred_null)*1000,
        Reartemp_diff=(Rpred_reartemp-Rpred_null)*1000,
        Rearprcp_diff=(Rpred_rearprcp-Rpred_null)*1000,
        Snow_diff=(Rpred_snow-Rpred_null)*1000,
        Spawntemp_diff=(Rpred_spawntemp-Rpred_null)*1000,
        Spawnprcp_diff=(Rpred_spawnprcp-Rpred_null)*1000,
        SSTsummer_diff=(Rpred_SSTsummer-Rpred_null)*1000,
        SSTwinter_diff=(Rpred_SSTwinter-Rpred_null)*1000,
        Pink_diff=(Rpred_Pink-Rpred_null)*1000,
        Chum_diff=(Rpred_Chum-Rpred_null)*1000,
        WinterFW_diff=(Rpred_WinterFW-Rpred_null)*1000,
        Spawners=Spawner_med*1000)

All_popns_median_summary <- All_popns_median2 %>% select(1,30:42) %>% 
  gather(Covariate,Value,2:13) %>% 
  group_by(population,Spawners,Covariate) %>% 
  summarise(Med=median(Value),
            q5=quantile(Value,0.05),
            q95=quantile(Value,0.95),
            sd=sd(Value))


All_popns_median_summary2 <-  All_popns_median_summary %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") #removing vars with no relationship

All_popns_median_summary2$Covariate <- as.factor(All_popns_median_summary2$Covariate)



#visualize by population
All_popns_median_plot <- All_popns_median_summary2 %>%
  mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot(aes(x=Covariate, y=Med))+
  geom_bar(stat="identity",width=0.75) +
  facet_wrap(~population)+
  xlab("Covariate") +
  ylab("Change in Recruits") +
  theme_bw() +
  geom_errorbar(aes(x=Covariate,ymin=Med-sd, ymax=Med+sd),width=.1) +
  scale_x_discrete(labels=c("Iceout_diff"="Ice\nout","Migtemp_diff"="Migration\ntemp",
                            "Pink_diff"="Pink\nSalmon","Rearprcp_diff"="Rearing\nprecip",
                            "Reartemp_diff"="Rearing\ntemp",
                               "Snow_diff"="Snow","Spawnprcp_diff"="Spawning\nprecip",
                               "Spawntemp_diff"="Spawning\ntemp",
                            "SSTwinter_diff"="Winter\nSST",
                            "WinterFW_diff"="Winter\nfreshwater"))+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#combined changed in recruitment based on median spawner values for each popn

Aggregate_change_median_neg <- All_popns_median_summary %>% filter(Med<0) 

Aggregate_change_median_pos <- All_popns_median_summary %>% filter(Med>0) 

Aggregate_change_median <- bind_rows(Aggregate_change_median_pos,Aggregate_change_median_neg)
  
Aggregate_change_median_summary <-  Aggregate_change_median %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") %>% 
  group_by(Covariate) %>% 
  summarise(sum_spawners=sum(Spawners),
            sum_maxR_change=sum(Med))

Aggregate_change_median_summary$Covariate <- as.factor(Aggregate_change_median_summary$Covariate)

#visualize combined
Recruit_change_median <- Aggregate_change_median_summary %>%
  mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot(aes(y=sum_maxR_change,x=Covariate))+
  geom_col()+ylab("Change in recruits")+xlab("Covariate")+theme_bw()+
  scale_x_discrete(labels=c("Iceout_diff"="Ice\nout","Migtemp_diff"="Migration\ntemp",
                            "Pink_diff"="Pink\nSalmon","Rearprcp_diff"="Rearing\nprecip",
                            "Reartemp_diff"="Rearing\ntemp",
                               "Snow_diff"="Snow","Spawnprcp_diff"="Spawning\nprecip",
                               "Spawntemp_diff"="Spawning\ntemp",
                            "SSTwinter_diff"="Winter\nSST",
                            "WinterFW_diff"="Winter\nfreshwater"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#is there a way to add some uncertainty to these combined plots?
  



  
