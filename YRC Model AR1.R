
####  Environmental drivers of YRC productivity Bayes hierarchical SR model ####

library(tidyverse)
library(postpack)
library(R2jags)
library(mcmcplots)
library(bayesboot)
library(lubridate)


#Code adapted from https://github.com/curryc2/Cook_Inlet_Chinook


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


#alternative analysis data
# 
# brood_table_alt <- read.csv(file.path(dir.data,"alt_brood_table.csv"))
# 
# brood_table_alt <- filter(brood_table_alt,year<2013)
# 
# Spawners <- brood_table_alt %>% filter(quantity=="S") %>%
#   select(pop,year,median) %>% mutate(median = median/1000) %>%
#     spread(year,median) %>%
#  select(-pop)
# 
# ln.Recruits <- brood_table_alt %>% filter(quantity=="logR") %>%
#   select(pop,year,median) %>% mutate(median=log(exp(median)/1000)) %>%
#   spread(year,median) %>% select(-pop)

# Merge covariate data ----------------------------------------------------

Ice_out <- read.csv(file.path(dir.data,"/Environmental data/Processed/Ice_out.csv"))
Migration_temp <- read.csv(file.path(dir.data,"/Environmental data/Processed/Weekly_max_migration_temp_t0.csv"))
rearing_GDD5 <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_GDD.csv"))
#rearing_GDD5 <- read.csv(file.path(dir.data,"/Environmental data/Processed/rearing_GDD_all.csv"))
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





#model without covariates
jags_model_null = function(){

  #Priors

  for(p in 1:n.pops) {
    exp.alpha[p] ~ dunif(0,50)
    alpha[p] <- log(exp.alpha[p])
    beta[p] ~ dnorm(0,pow(0.1,-2));T(0,100)
    sigma.oe[p] ~ dnorm(0, pow(1,-2));T(1e-3,100)
     phi[p] ~ dunif(-0.99, 0.99)

  }#next p



  #PREDICTIONS
  for(p in 1:n.pops) {
    # First Year

    pred.rec.1[p,1] <- Sobs[p,1]*exp(alpha[p] - Sobs[p,1]*beta[p])

    log.pred.rec.1[p,1] <- log(pred.rec.1[p,1])
    log.resid[p,1] <- lnRobs[p,1] - log.pred.rec.1[p,1]
    log.pred.rec.2[p,1] <- log.pred.rec.1[p,1]
    log.resid.2[p,1] <- log.resid[p,1]

   #Subsequent Years
     for(y in 2:n.years[p]) {

    pred.rec.1[p,y] <- Sobs[p,y]*exp(alpha[p] - Sobs[p,y]*beta[p])

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


#### Specify initial values ####

jags_inits = function() {
  exp.alpha <- runif(n.pops,1,2)
  #beta <- runif(n.pops,0,0.001)
  beta <- runif(n.pops,0,0.3) 
  sigma.oe <- runif(n.pops,0.1,1) 
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

out_null <- jags.parallel(data=jags_data,
  model.file=jags_model_null,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=50000, 
  n.burnin=10000)

#full iterations
out_250 <- jags.parallel(data=jags_data,
  model.file=jags_model,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=250000, 
  n.burnin=15000)  

out_null_250 <- jags.parallel(data=jags_data,
  model.file=jags_model_null,
  inits=jags_inits,
  parameters.to.save=jags_params,
  n.chains=3, 
  n.thin=20, 
  n.iter=250000, 
  n.burnin=15000)

#Save
saveRDS(out_250, file=file.path(dir.output,"out_250.rds"))

out.mcmc <- as.mcmc(out)

#Write Output File for Diagnostics
write.csv(out$BUGSoutput$summary, file=file.path(dir.figs,"out_Aug23.csv"))
write.csv(out_nocomp$BUGSoutput$summary, file=file.path(dir.figs,"out_nocomp.csv"))
write.csv(out_250$BUGSoutput$summary, file=file.path(dir.figs,"out_250.csv"))

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

labels3=c("Winter SST","Spawning temperature","FW winter temperature",
          "Annual snowpack","Summer SST","Chum salmon","Spawning precipitation",
          "Pink salmon","Rearing degree days","Ice out date",
          "Migration temperature","Rearing precipitation")

#ordered by life stage
labels4=c("Migration temperature","Spawning temperature","Spawning precipitation",
          "Snowpack","FW winter temperature","Rearing degree days",
          "Rearing precipitation","Ice out date","Summer SST",
          "Winter SST","Pink salmon","Chum salmon")

#rearrange variables

mu.coef.vector <- out$BUGSoutput$sims.list$mu.coef
mu.coef.vector.new <- mu.coef.vector[,c(5,9,6,8,1,2,3,7,4)]
mu.coef.vector.new2 <- mu.coef.vector[,c(9,6,12,5,8,11,7,10,1,3,2,4)]
mu.coef.vector.new3 <- mu.coef.vector[,c(2,6,7,5,12,3,4,1,8,9,10,11)]


#Hyper means summarized in order of median effect
par(mfcol=c(1,1), mar=c(2,12,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(out_250$BUGSoutput$sims.list$mu.coef,
            labels=labels3,reorder=TRUE, quantiles=list(0.025,0.25,0.75,0.975), 
            style='plain', col='blue',cex=1.1)
caterpoints(apply(mu.coef.vector.new2,2,median), pch=21, col='red', bg='orange')
abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Change in log (recruits/spawner)', side=1, outer=TRUE, font=2, line=0.5,at=0.7,cex=1.2)
mtext('Environmental variable', side=2, outer=TRUE, font=2, line=0.5,cex=1.2)

#Hyper means by life stage
par(mfcol=c(1,1), mar=c(2,12,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(mu.coef.vector.new3,
            labels=labels4,reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), 
            style='plain', col=c('blue'),cex=1.1)
caterpoints(apply(mu.coef.vector.new3,2,median), pch=21, col='red', bg='orange')
abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
mtext('Change in log (recruits/spawner)', side=1, outer=TRUE, font=2, line=0.5,at=0.7,cex=1.2)
mtext('Environmental variable', side=2, outer=TRUE, font=2, line=0.5,cex=1.2)

#739x503 saved dimensions

#try ggplot

mu.coef.vector.new3 <- as.data.frame(mu.coef.vector.new3)

colnames(mu.coef.vector.new3) <- c("Migration temperature","Spawning temperature","Spawning precipitation",
          "Snowpack","FW winter temperature","Rearing degree days",
          "Rearing precipitation","Ice out date","Summer SST",
          "Winter SST","Pink salmon","Chum salmon")

mu.coef.vector.long <- mu.coef.vector.new3 %>% gather("var","value",1:12)

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

# Regional_effects <- mu.coef.vector.long %>% 
#   mutate(var=fct_relevel(var,"Chum salmon","Pink salmon","Winter SST",
#                          "Summer SST","Ice out date","Rearing precipitation",
#                          "Rearing degree days","FW winter temperature","Snowpack",
#                          "Spawning precipitation","Spawning temperature","Migration temperature")) %>% 
#           ggplot(aes(x=var, y=value,fill=stage)) +
#   geom_hline(yintercept = 0, colour='red')+
#   stat_eye()+
#             theme_bw() +
#             stat_summary(fun="q.95", colour="black", geom="line", lwd=0.5) +
#             stat_summary(fun="q.50", colour="black", geom="line", lwd=1.25) +
#             stat_summary(fun="median", colour="white", size=1.75, geom="point", pch=19) +
#             stat_summary(fun="median", colour="black", size=1.75, geom="point", pch=21) +
#             scale_x_discrete(limits=rev(levels(mu.coef.vector.long$var))) +
#             coord_flip() +
#   theme(legend.position='top') +
#             xlab("Environmental variable") + ylab("Change in log (recruits/spawner)") +
#   theme(text = element_text(size=25),axis.text=element_text(size=15),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
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

#785x620 saved dimensions

Regional_effects_integrated <- mu.coef.vector.long %>% 
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
  theme(legend.text=element_text(size=15))+scale_y_continuous(limits = c(-0.55,0.6))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

Regional_effects_integrated

names.covars.updated <- c("Ice out date","Migration temperature","Rearing degree days","Rearing precipitation",
                  "Snowpack","Spawning temperature","Spawning precipitation","Summer SST",
                  "Winter SST","Pink salmon","Chum salmon","FW winter temperature")

#Population specific covariate "effects"
pdf(file="Plots/Popn_plot_Sep29_alternative.pdf", height=11, width=10)
par(mfcol=c(4,3), mar=c(2,7,4,1), oma=c(4,4,2,2))
c <- 1
for(c in 1:n.covars) {
  caterplot(out$BUGSoutput$sims.list$coef[,,c],
            labels=pops2, cex.labels=TRUE, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975),
            style='plain', col='blue',val.lim = c(-0.45, 0.5))
  mtext(names.covars.updated[c], side=3, outer=FALSE, line=1)
  caterpoints(apply(out$BUGSoutput$sims.list$coef[,,c],2,median), reorder=FALSE, pch=21, col='red', bg='orange')
  abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
}
mtext('Effect on salmon recruitment', side=1, outer=TRUE, font=2, cex=1.2,line=2)
mtext('Population', side=2, outer=TRUE, font=2, cex=1.2,line=2)

dev.off()

#re-do plot in ggplot

check <- out_250$BUGSoutput$sims.list$coef

#data structure is 1:35250 (values), 1:8 (popns), 1:12 (covars)

mu.coef.popn <- as.data.frame.table(out$BUGSoutput$sims.list$coef)
names(mu.coef.popn) <- c("reps","population","var","value")

#rename factor levels
levels(mu.coef.popn$population) <- c("Carmacks","Lower Mainstem","Middle Mainstem","Pelly","Stewart","Teslin",
        "Upper Mainstem","White-Donjek")

levels(mu.coef.popn$var) <- c("Ice out date","Migration temperature","Rearing degree days","Rearing precipitation",
                  "Snowpack","Spawning temperature","Spawning precipitation","Summer SST",
                  "Winter SST","Pink salmon","Chum salmon","FW winter temperature")

#add stage

mu.coef.popn <- mu.coef.popn %>% 
  mutate(stage=if_else(var=="Chum salmon"|var=="Pink salmon"|var=="Summer SST"|
                         var=="Winter SST","Marine",
                       if_else(var=="Ice out date","Outmigration",
                               if_else(var=="Rearing degree days"|var=="Rearing precipitation",
                                       "Freshwater juvenile",
                                         "Spawning and incubation"))))


Popn_effects <- mu.coef.popn %>% 
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

Popn_effects_alt <- mu.coef.popn %>% 
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
  theme(legend.text=element_text(size=20))+scale_y_continuous(limits = c(-0.6, 1),
                                                              labels=c(-0.5,0,0.5,1.0),
                                                              breaks=c(-0.5,0,0.5,1))

Popn_effects_alt
#saved dimensions were 1200x1500

# Popn_effects2 <- mu.coef.popn %>% 
#   mutate(var=fct_relevel(var,"Migration temperature","Spawning temperature","Spawning precipitation",
#                          "Snowpack","FW winter temperature","Rearing degree days",
#                          "Rearing precipitation","Ice out date","Summer SST",
#                          "Winter SST","Pink salmon","Chum salmon")) %>% 
#   mutate(stage=fct_relevel(stage,"Freshwater adult","Freshwater juvenile",
#                             "Outmigration","Marine")) %>% 
#   mutate(population=fct_relevel(population,"White-Donjek","Upper Mainstem","Teslin",
#                                 "Stewart","Pelly","Middle Mainstem",
#                                 "Lower Mainstem","Carmacks")) %>% 
#   ggplot(aes(x=population, y=value)) +
#   geom_hline(yintercept = 0, colour='red')+
#   theme_bw() +
#             coord_flip() +
#   stat_summary(fun="q.95", colour="black", geom="line", lwd=0.5,show.legend=FALSE) +
#   stat_summary(fun="q.50", colour="black", geom="line", lwd=1.25,show.legend=FALSE) +
#   stat_summary(fun="median", colour="black", size=2.5, geom="point", pch=21, show.legend=FALSE) +
#   stat_summary(fun="median", colour="gray", size=0.5, geom="point", pch=19,show.legend=FALSE) +
#            theme(legend.position='top') +
#   theme(legend.title=element_blank())+
#   facet_wrap(~var,ncol=3)+
#             xlab("Population") + ylab("Change in log (recruits/spawner)") +
#   theme(text = element_text(size=20),axis.text=element_text(size=10),
#         axis.title.x = element_text(margin=margin(t=15,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))+
#   theme(strip.text.x = element_text(size = 10))+
#   theme(legend.text=element_text(size=13))+scale_y_continuous(limits = c(-0.45, 0.5),
#                                                               labels=c(-0.4,-0.2,0,0.2,0.4),
#                                                               breaks=c(-0.4,-0.2,0,0.2,0.4))
# 
# Popn_effects2

# 
# tiff(file="Plots/Popn_plot_Oct25.tif", height=11, width=10,units="in",res=300)
# par(mfcol=c(4,3), mar=c(2,7,4,1), oma=c(4,4,2,2))
# c <- 1
# for(c in 1:n.covars) {
#   caterplot(out_250$BUGSoutput$sims.list$coef[,,c],
#             labels=pops2, cex.labels=TRUE, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975),
#             style='plain', col='blue',val.lim = c(-0.45, 0.5))
#   mtext(names.covars.updated[c], side=3, outer=FALSE, line=1)
#   caterpoints(apply(out_250$BUGSoutput$sims.list$coef[,,c],2,median), reorder=FALSE, pch=21, col='red', bg='orange')
#   abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
# }
# mtext('Change in log (recruits/spawner)', side=1, outer=TRUE, font=2, cex=1.2,line=2)
# mtext('Population', side=2, outer=TRUE, font=2, cex=1.2,line=2)
# 
# dev.off()
# 
# tiff(file="Plots/Popn_plot_Oct25_alt.tif", height=11, width=10,units="in",res=300)
# par(mfcol=c(4,3), mar=c(2,7,4,1), oma=c(4,4,2,2))
# c <- 1
# for(c in 1:n.covars) {
#   caterplot(out$BUGSoutput$sims.list$coef[,,c],
#             labels=pops2, cex.labels=TRUE, reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975),
#             style='plain', col='blue',val.lim = c(-0.6, 1))
#   mtext(names.covars.updated[c], side=3, outer=FALSE, line=1)
#   caterpoints(apply(out$BUGSoutput$sims.list$coef[,,c],2,median), reorder=FALSE, pch=21, col='red', bg='orange')
#   abline(v=0, lty=1, lwd=2, col=rgb(1,0,0, alpha=0.5))
# }
# mtext('Change in log (recruits/spawner)', side=1, outer=TRUE, font=2, cex=1.2,line=2)
# mtext('Population', side=2, outer=TRUE, font=2, cex=1.2,line=2)
# 
# dev.off()
# 
# pdf(file="Plots/Results_Aug5.pdf",width=9,height=6)

#Hyper means summarized
par(mfcol=c(1,1), mar=c(2,12,3,1), oma=c(2,2,1,1))
c <- 1
  caterplot(out$BUGSoutput$sims.list$mu.coef,
            labels=names.covars,reorder=FALSE, quantiles=list(0.025,0.25,0.75,0.975), 
            style='gray', col='blue',cex=1.1)
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
  log.pred <- apply(out$BUGSoutput$sims.list$log.pred.rec.2[,p,1:n.years[p]],2, 
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


#fast plot for Salmon Gathering

per_change <- read.csv(file.path(dir.data,"/Environmental data/Processed/Percent_change_Aug18.csv"))

per_change$Var <- as.factor(per_change$Var)

per_change <- per_change %>% filter(Var!="chum"&Var!="summer SST")

Per_change_median <- per_change %>%
  mutate(Var=fct_relevel(Var,"rearing prcp","migration temp",
                               "rearing temp","Ice out","pink",
                               "spawning prcp","snow","FW winter",
                               "spawning temp","winter SST")) %>%
  ggplot()+
  geom_point(aes(y=median,x=Var),size=10)+
  geom_hline(yintercept = 0, lty = "dotted",size=1) +
  ylab("Change in the number of returning salmon (%)")+xlab("Change in the environment")+theme_bw()+
  geom_errorbar(aes(x=Var,ymin=q5, ymax=q95),width=.1) +
  scale_x_discrete(labels=c("Ice out"="Ice\nout\n\n(+4 days)","migration temp"="Migration\ntemp\n\n(+1.2°C)",
                            "pink"="Pink\nSalmon\n\n(+61 million)","rearing prcp"="Rearing\nprecip\n\n(+29 mm)",
                            "rearing temp"="Rearing\nGDD\n\n(+165 GDD)",
                               "snow"="Snow\n\n\n(+61 mm)","spawning prcp"="Spawning\nprecip\n\n(+45 mm)",
                               "spawning temp"="Spawning\ntemp\n\n(+1.5°C)",
                            "winter SST"="Winter\nSST\n\n(+0.8°C)",
                            "FW winter"="Winter\nfreshwater\n\n(+2.7°C)"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

# Infer change in number of recruits --------------------------------------

#Processing
alpha <- as.data.frame(out_250$BUGSoutput$sims.list$alpha)

alpha_sample <- alpha %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,alpha,1:8)

beta <- as.data.frame(out_250$BUGSoutput$sims.list$beta)

beta_sample <- beta %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,beta,1:8)
  
mig_temp <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,2])

mig_temp_sample <- mig_temp %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,migtemp,1:8)

ice_out <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,1])

ice_out_sample <- ice_out %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,ice_out,1:8)

rear_temp <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,3])

rear_temp_sample <- rear_temp %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,reartemp,1:8)

rear_precip <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,4])

rear_precip_sample <- rear_precip %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,rear_precip,1:8)

snow <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,5])

snow_sample <- snow %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,snow,1:8)

spawn_temp <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,6])

spawn_temp_sample <- spawn_temp %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,spawntemp,1:8)

spawn_precip <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,7])

spawn_precip_sample <- spawn_precip %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,spawnprecip,1:8)

SST_summ <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,8])

SST_summ_sample <- SST_summ %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,SSTsummer,1:8)

SST_wint <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,9])

SST_wint_sample <- SST_wint %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,SSTwinter,1:8)

Pink <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,10])

Pink_sample <- Pink %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,Pink,1:8)

Chum <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,11])

Chum_sample <- Chum %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame() %>% gather(population,Chum,1:8)

Winter <- as.data.frame(out_250$BUGSoutput$sims.list$coef[,,12])

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
            Spawner25=quantile(S_med,0.25),
            Spawner75=quantile(S_med,0.75),
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

All_popns_median_summary <- All_popns_median2 %>% select(1,32:44) %>% 
  gather(Covariate,Value,2:13) %>% 
  group_by(population,Spawners,Covariate) %>% 
  summarise(Med=median(Value),
            q5=quantile(Value,0.05),
            q25=quantile(Value,0.25),
            q75=quantile(Value,0.75),
            q95=quantile(Value,0.95),
            q50=quantile(Value,0.5))


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
  

Aggregate_change_median_summary <-  All_popns_median_summary %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") %>% 
  group_by(Covariate) %>% 
  summarise(sum_spawners=sum(Spawners),
            sum_maxR_change=sum(Med),
            sum_q5=sum(q5),
            sum_q25=sum(q25),
            sum_q75=sum(q75),
            sum_q95=sum(q95),
            sum_q50=sum(q50))

#with additional variables
# Aggregate_change_median_summary <-  All_popns_median_summary %>% 
#   group_by(Covariate) %>% 
#   summarise(sum_spawners=sum(Spawners),
#             sum_maxR_change=sum(Med),
#             sum_q5=sum(q5),
#             sum_q95=sum(q95),
#             sum_q50=sum(q50))

Aggregate_change_median_summary$Covariate <- as.factor(Aggregate_change_median_summary$Covariate)

#visualize combined

Recruit_change_median <- Aggregate_change_median_summary %>%
  mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot()+
  geom_point(aes(y=sum_maxR_change,x=Covariate),size=10)+
  geom_hline(yintercept = 0, lty = "dotted",size=1) +
  ylab("Change in the number of recruits")+xlab("Change in the environment")+theme_bw()+
  geom_errorbar(aes(x=Covariate,ymin=sum_q5, ymax=sum_q95),width=.1) +
  scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
                            "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
                            "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
                               "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
                               "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
                            "SSTwinter_diff"="Winter\nSST\n\n(+0.8°C)",
                            "WinterFW_diff"="FW winter\ntemp\n\n(+2.7°C)"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Recruit_change_median_violin <- Aggregate_change_median_summary %>%
#   mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
#                                "Reartemp_diff","Iceout_diff","Pink_diff",
#                                "Spawnprcp_diff","Snow_diff","WinterFW_diff",
#                                "Spawntemp_diff","SSTwinter_diff")) %>%
#   ggplot()+
#   geom_boxplot(aes(y=c(ymin=sum_q5,lower=sum_q25,middle=sum_q50,upper=sum_q75,ymax=sum_q95),x=Covariate),size=10)+
#   geom_hline(yintercept = 0, lty = "dotted",size=1) +
#   ylab("Change in the number of recruits")+xlab("Change in the environment")+theme_bw()+
#     scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
#                             "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
#                             "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
#                                "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
#                                "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
#                             "SSTwinter_diff"="Winter\nSST\n\n(+0.8°C)",
#                             "WinterFW_diff"="FW winter\ntemp\n\n(+2.7°C)"))+
#   theme(text = element_text(size=25),axis.text=element_text(size=15),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


DF <- data.frame(x=Aggregate_change_median_summary$Covariate, 
                 min=Aggregate_change_median_summary$sum_q5, 
                 low=Aggregate_change_median_summary$sum_q25, 
                 mid=Aggregate_change_median_summary$sum_q50, 
                 top=Aggregate_change_median_summary$sum_q75, 
                 max=Aggregate_change_median_summary$sum_q95)

ggplot(DF, aes(x=x, ymin = min, lower = low, middle = mid, upper = top, ymax = max)) +
  geom_boxplot(stat = "identity",width=12)

Recruit_change_median_boxplot <- DF %>%
  mutate(x=fct_relevel(x,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot(aes(x=x, ymin = min, lower = low, middle = mid, upper = top, ymax = max)) +
  geom_boxplot(stat = "identity",width=0.5, position = position_dodge(width=5))+geom_hline(yintercept = 0, lty = "dotted",size=1) +
  ylab("Change in the number of recruits")+xlab("Change in the environment")+theme_bw()+
    scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
                            "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
                            "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
                               "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
                               "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
                            "SSTwinter_diff"="Winter\nSST\n\n(+0.8°C)",
                            "WinterFW_diff"="FW winter\ntemp\n\n(+2.7°C)"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#reordered by life stage


DF <- DF %>% 
  mutate(stage=if_else(x=="Pink_diff"|x=="SSTwinter_diff","Marine",
                       if_else(x=="Iceout_diff","Outmigration",
                               if_else(x=="Reartemp_diff"|x=="Rearprcp_diff",
                                       "Freshwater juvenile",
                                         "Spawning and incubation"))))

Stage_boxplot <- DF %>%
  mutate(x=fct_relevel(x,"Migtemp_diff","Spawntemp_diff","Spawnprcp_diff",
                       "Snow_diff","WinterFW_diff","Reartemp_diff","Rearprcp_diff",
                               "Iceout_diff","SSTwinter_diff","Pink_diff")) %>%
  mutate(stage=fct_relevel(stage,"Spawning and incubation","Freshwater juvenile",
                            "Outmigration","Marine")) %>% 
  ggplot(aes(x=x, ymin = min, lower = low, middle = mid, upper = top, ymax = max,fill=stage)) +
  geom_boxplot(stat = "identity",width=0.5, position = position_dodge(width=5))+geom_hline(yintercept = 0, lty = "dotted",size=1) +
  ylab("Change in the number of recruits")+xlab("Change in the environment")+theme_bw()+
  
  scale_fill_manual(values = c("#FE9000","#FFDD4A","#5ADBFF","#3C6997"))+
    scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
                            "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
                            "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
                               "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
                               "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
                            "SSTwinter_diff"="Winter\nSST\n\n(+0.8°C)",
                            "WinterFW_diff"="FW winter\ntemp\n\n(+2.7°C)"))+
  theme(legend.position='top') +
  theme(legend.title=element_blank())+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Stage_boxplot


#saved dimensions are: 

#same thing at 25th percentile

All_popns_median2 <- All_popns_median %>% 
  mutate(Rpred_migtemp = Spawner25*exp(alpha-beta*Spawner25+migtemp*1), #change in 1SD for migration temp
         Rpred_iceout = Spawner25*exp(alpha-beta*Spawner25+ice_out*1),
         Rpred_reartemp = Spawner25*exp(alpha-beta*Spawner25+reartemp*1),
         Rpred_rearprcp = Spawner25*exp(alpha-beta*Spawner25+rear_precip*1),
         Rpred_snow = Spawner25*exp(alpha-beta*Spawner25+snow*1),
         Rpred_spawntemp = Spawner25*exp(alpha-beta*Spawner25+spawntemp*1),
         Rpred_spawnprcp = Spawner25*exp(alpha-beta*Spawner25+spawnprecip*1),
         Rpred_SSTsummer = Spawner25*exp(alpha-beta*Spawner25+SSTsummer*1),
         Rpred_SSTwinter = Spawner25*exp(alpha-beta*Spawner25+SSTwinter*1),
         Rpred_Pink = Spawner25*exp(alpha-beta*Spawner25+Pink*1),
         Rpred_Chum = Spawner25*exp(alpha-beta*Spawner25+Chum*1),
         Rpred_WinterFW = Spawner25*exp(alpha-beta*Spawner25+WinterFW*1),
        Rpred_null = Spawner25*exp(alpha-beta*Spawner25),
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
        Spawners=Spawner25*1000)

All_popns_median_summary <- All_popns_median2 %>% select(1,32:44) %>% 
  gather(Covariate,Value,2:13) %>% 
  group_by(population,Spawners,Covariate) %>% 
  summarise(Med=median(Value),
            q5=quantile(Value,0.05),
            q95=quantile(Value,0.95),
            q50=quantile(Value,0.5))


All_popns_median_summary2 <-  All_popns_median_summary %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") #removing vars with no relationship

All_popns_median_summary2$Covariate <- as.factor(All_popns_median_summary2$Covariate)



#combined changed in recruitment based on 25th percentile values for each popn
  

Aggregate_change_median_summary25 <-  All_popns_median_summary %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") %>% 
  group_by(Covariate) %>% 
  summarise(sum_spawners=sum(Spawners),
            sum_maxR_change=sum(Med),
            sum_q5=sum(q5),
            sum_q95=sum(q95),
            sum_q50=sum(q50))

Aggregate_change_median_summary$Covariate <- as.factor(Aggregate_change_median_summary$Covariate)

#visualize combined

Recruit_change_25 <- Aggregate_change_median_summary25 %>%
  mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot()+
  geom_point(aes(y=sum_maxR_change,x=Covariate),size=10)+
  geom_hline(yintercept = 0, lty = "dotted",size=1) +
  ylab("Change in the number of recruits")+xlab("Change in the environment")+theme_bw()+
  geom_errorbar(aes(x=Covariate,ymin=sum_q5, ymax=sum_q95),width=.1) +
  scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
                            "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
                            "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
                               "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
                               "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
                            "SSTwinter_diff"="Winter\nmarine\n\n(+0.8°C)",
                            "WinterFW_diff"="Winter\nfreshwater\n\n(+2.7°C)"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#same thing at 75th percentile

All_popns_median2 <- All_popns_median %>% 
  mutate(Rpred_migtemp = Spawner75*exp(alpha-beta*Spawner75+migtemp*1), #change in 1SD for migration temp
         Rpred_iceout = Spawner75*exp(alpha-beta*Spawner75+ice_out*1),
         Rpred_reartemp = Spawner75*exp(alpha-beta*Spawner75+reartemp*1),
         Rpred_rearprcp = Spawner75*exp(alpha-beta*Spawner75+rear_precip*1),
         Rpred_snow = Spawner75*exp(alpha-beta*Spawner75+snow*1),
         Rpred_spawntemp = Spawner75*exp(alpha-beta*Spawner75+spawntemp*1),
         Rpred_spawnprcp = Spawner75*exp(alpha-beta*Spawner75+spawnprecip*1),
         Rpred_SSTsummer = Spawner75*exp(alpha-beta*Spawner75+SSTsummer*1),
         Rpred_SSTwinter = Spawner75*exp(alpha-beta*Spawner75+SSTwinter*1),
         Rpred_Pink = Spawner75*exp(alpha-beta*Spawner75+Pink*1),
         Rpred_Chum = Spawner75*exp(alpha-beta*Spawner75+Chum*1),
         Rpred_WinterFW = Spawner75*exp(alpha-beta*Spawner75+WinterFW*1),
        Rpred_null = Spawner75*exp(alpha-beta*Spawner75),
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
        Spawners=Spawner75*1000)

All_popns_median_summary <- All_popns_median2 %>% select(1,32:44) %>% 
  gather(Covariate,Value,2:13) %>% 
  group_by(population,Spawners,Covariate) %>% 
  summarise(Med=median(Value),
            q5=quantile(Value,0.05),
            q95=quantile(Value,0.95),
            q50=quantile(Value,0.5))


All_popns_median_summary2 <-  All_popns_median_summary %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") #removing vars with no relationship

All_popns_median_summary2$Covariate <- as.factor(All_popns_median_summary2$Covariate)



#combined changed in recruitment based on 25th percentile values for each popn
  

Aggregate_change_median_summary75 <-  All_popns_median_summary %>% 
  filter(Covariate!="Chum_diff"&Covariate!="SSTsummer_diff") %>% 
  group_by(Covariate) %>% 
  summarise(sum_spawners=sum(Spawners),
            sum_maxR_change=sum(Med),
            sum_q5=sum(q5),
            sum_q95=sum(q95),
            sum_q50=sum(q50))

Aggregate_change_median_summary$Covariate <- as.factor(Aggregate_change_median_summary$Covariate)

#visualize combined

Recruit_change_75 <- Aggregate_change_median_summary75 %>%
  mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
                               "Reartemp_diff","Iceout_diff","Pink_diff",
                               "Spawnprcp_diff","Snow_diff","WinterFW_diff",
                               "Spawntemp_diff","SSTwinter_diff")) %>%
  ggplot()+
  geom_point(aes(y=sum_maxR_change,x=Covariate),size=10)+
  geom_hline(yintercept = 0, lty = "dotted",size=1) +
  ylab("Change in the number of recruits")+xlab("Change in the environment")+theme_bw()+
  geom_errorbar(aes(x=Covariate,ymin=sum_q5, ymax=sum_q95),width=.1) +
  scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
                            "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
                            "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
                               "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
                               "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
                            "SSTwinter_diff"="Winter\nmarine\n\n(+0.8°C)",
                            "WinterFW_diff"="Winter\nfreshwater\n\n(+2.7°C)"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

# #original plot for gathering
# Recruit_change_median <- Aggregate_change_median_summary %>%
#   mutate(Covariate=fct_relevel(Covariate,"Rearprcp_diff","Migtemp_diff",
#                                "Reartemp_diff","Iceout_diff","Pink_diff",
#                                "Spawnprcp_diff","Snow_diff","WinterFW_diff",
#                                "Spawntemp_diff","SSTwinter_diff")) %>%
#   ggplot()+
#   geom_point(aes(y=sum_maxR_change,x=Covariate),size=10)+
#   geom_hline(yintercept = 0, lty = "dotted",size=1) +
#   ylab("Change in the number of returning salmon")+xlab("Change in the environment")+theme_bw()+
#   geom_errorbar(aes(x=Covariate,ymin=sum_q5, ymax=sum_q95),width=.1) +
#   scale_x_discrete(labels=c("Iceout_diff"="Ice\nout\n\n(+4 days)","Migtemp_diff"="Migration\ntemp\n\n(+1.2°C)",
#                             "Pink_diff"="Pink\nSalmon\n\n(+61 million)","Rearprcp_diff"="Rearing\nprecip\n\n(+29 mm)",
#                             "Reartemp_diff"="Degree\nDays\n\n(+165 DD)",
#                                "Snow_diff"="Snow\n\n\n(+61 mm)","Spawnprcp_diff"="Spawning\nprecip\n\n(+45 mm)",
#                                "Spawntemp_diff"="Spawning\ntemp\n\n(+1.5°C)",
#                             "SSTwinter_diff"="Winter\nmarine\n\n(+0.8°C)",
#                             "WinterFW_diff"="Winter\nfreshwater\n\n(+2.7°C)"))+
#   theme(text = element_text(size=25),axis.text=element_text(size=15),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))+
#   theme(plot.background = element_rect(fill="#F9F4ED"),
#         panel.background = element_rect(fill="#F9F4ED"))



#saved svg 1300x900

#organize vars by colour. 
#try adding low and high years also
#consider removing vars that weren't significant at the main level - spawn temp and prcp


pred_rec1 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,1,])
pred_rec2 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,2,])
pred_rec3 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,3,])
pred_rec4 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,4,])
pred_rec5 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,5,])
pred_rec6 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,6,])
pred_rec7 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,7,])
pred_rec8 <- as.data.frame(out$BUGSoutput$sims.list$log.pred.rec.2[,8,])

pred_rec1 <- pred_rec1 %>% mutate(population="Carmacks") %>% 
  gather(year,recruitment,1:28)
  
pred_rec2 <- pred_rec2 %>% mutate(population="Lower Mainstem") %>% 
  gather(year,recruitment,1:28)

pred_rec3 <- pred_rec3 %>% mutate(population="Middle Mainstem") %>% 
  gather(year,recruitment,1:28)

pred_rec4 <- pred_rec4 %>% mutate(population="Pelly") %>% 
  gather(year,recruitment,1:28)

pred_rec5 <- pred_rec5 %>% mutate(population="Stewart") %>% 
  gather(year,recruitment,1:28)

pred_rec6 <- pred_rec6 %>% mutate(population="Teslin") %>% 
  gather(year,recruitment,1:28)

pred_rec7 <- pred_rec7 %>% mutate(population="Upper Lakes and Mainstem") %>% 
  gather(year,recruitment,1:28)

pred_rec8 <- pred_rec8 %>% mutate(population="White-Donjek") %>% 
  gather(year,recruitment,1:28)

pred_rec_all <- bind_rows(pred_rec1,pred_rec2,pred_rec3,pred_rec4,pred_rec5,pred_rec6,pred_rec7,pred_rec8)

pred_rec_all_summary <- pred_rec_all %>% mutate(recruitment=exp(recruitment)) %>% 
  group_by(population,year) %>% 
  summarise(mean_recruitment=mean(recruitment),
            sd=sd(recruitment),
            q5=quantile(recruitment,0.05),
            q95=quantile(recruitment,0.95))
  
pred_rec_all_summary2 <- pred_rec_all_summary %>% 
  group_by(year) %>% 
  summarise(total_recruitment=sum(mean_recruitment),
            total_sd=sum(sd),
            total_5=sum(q5),
            total_95=sum(q95))

pred_rec_all_summary2 <- pred_rec_all_summary2 %>% mutate(year=as.numeric(substr(year,2,4))) %>% 
  arrange(year)

pred_rec_all_summary2$year <- 1985:2012

#need to check that this makes sense with adding the quantiles
  
Recruit_trend <- pred_rec_all_summary2 %>% ggplot()+
  geom_point(aes(y=total_recruitment,x=year),size=5)+
  geom_line(aes(y=total_recruitment,x=year)) + 
  geom_hline(yintercept = 110.57, lty = "dotted",size=1) +
  geom_ribbon(aes(x=year,ymin = total_5, ymax = total_95),alpha=0.4) +
                  ylab("Recruitment (000s)")+xlab("Year")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


# Productivity index ------------------------------------------------------

res1 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,1,])
res2 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,2,])
res3 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,3,])
res4 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,4,])
res5 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,5,])
res6 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,6,])
res7 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,7,])
res8 <- as.data.frame(out_null$BUGSoutput$sims.list$log.resid[,8,])

res1 <- res1 %>% mutate(population="Carmacks") %>% 
  gather(year,productivity,1:28)
  
res2 <- res2 %>% mutate(population="Lower Mainstem") %>% 
  gather(year,productivity,1:28)

res3 <- res3 %>% mutate(population="Middle Mainstem") %>% 
  gather(year,productivity,1:28)

res4 <- res4 %>% mutate(population="Pelly") %>% 
  gather(year,productivity,1:28)

res5 <- res5 %>% mutate(population="Stewart") %>% 
  gather(year,productivity,1:28)

res6 <- res6 %>% mutate(population="Teslin") %>% 
  gather(year,productivity,1:28)

res7 <- res7 %>% mutate(population="Upper Lakes and Mainstem") %>% 
  gather(year,productivity,1:28)

res8 <- res8 %>% mutate(population="White-Donjek") %>% 
  gather(year,productivity,1:28)

res_all <- bind_rows(res1,res2,res3,res4,res5,res6,res7,res8)

res_all_summary <- res_all %>% group_by(population) %>% 
  mutate(productivity=scale(productivity)) %>% 
  ungroup() %>% 
  group_by(population,year) %>% 
  summarise(med=quantile(productivity,0.5),
            sd=sd(productivity),
            q5=quantile(productivity,0.05),
            q95=quantile(productivity,0.95))

res_all_summary <- res_all_summary %>% mutate(year=as.numeric(substr(year,2,4))) %>%
  arrange(year)

res_all_summary$year <- rep(1985:2012,each=8)

ggplot(res_all_summary,aes(y=med,x=year))+
  geom_point()+geom_smooth(method="lm")

#write.csv(res_all_summary,file=file.path(dir.data,"/Environmental data/Processed/productivity_res_noAR1.csv"),row.names = FALSE)

# Productivity_trend_popn <- res_all_summary %>% ggplot()+
#   geom_point(aes(y=med,x=year),size=5)+
#   geom_line(aes(y=med,x=year)) + geom_smooth(aes(y=med,x=year),method="lm") + 
#   geom_hline(yintercept = 0, lty = "dotted",size=1) +
#   facet_wrap(~population, ncol=2) + 
#   geom_ribbon(aes(x=year,ymin = q5, ymax = q95),alpha=0.4) +
#                   ylab("Productivity index")+xlab("Year")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Productivity_trend_popn <- res_all_summary %>% ggplot()+
  geom_path(aes(y=med,x=year,group=population)) + geom_smooth(aes(y=med,x=year),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +
                  ylab("Productivity index")+xlab("Year")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

summary(lm(med~year,data=res_all_summary))
#1.6 unit decrease over the 28 year period

#add mean popn line

res_all_summary2 <- res_all_summary %>% 
  group_by(year) %>% 
  summarise(med=mean(med)) %>% 
  mutate(population="all")

res_all_summary3 <- bind_rows(res_all_summary2,res_all_summary)

data_breaks <- data.frame(start = c(1985, 1991.5, 2001.5, 2010.5), 
                          end = c(1991.5, 2001.5, 2010.5, 2012),
                          colors = c("Above-Average","Average","Below-Average","Average"))

Productivity_trend_popn2 <- res_all_summary3 %>% ggplot()+
  geom_rect(data=data_breaks,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf,fill=colors),alpha=0.5)+
  scale_fill_manual(values=c("#7a82ab","#C6D4FF","#307473","#C6D4FF"))+
  geom_line(aes(y=med,x=year,group=population))+
  labs(fill="")+
  geom_path(data= res_all_summary3[res_all_summary3$population=="all",],aes(y=med,x=year),size=2) +
    geom_hline(yintercept = 0, lty = "dotted",size=1) +
                  ylab("Productivity index")+xlab("Brood Year")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))



res_all_summary4 <- res_all_summary %>% 
  mutate(time_period=if_else(year<2000,"pre","post")) %>% 
  group_by(time_period) %>% 
  summarise(mean_prod=mean(med))

mean_trend <- res_all_summary3 %>% filter(population=="all")