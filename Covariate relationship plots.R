####Explore covariate relationships####

library(tidyverse)
library(mgcv)
library(visreg)

#Define Workflow Paths ====================================================
wd <- file.path(getwd())
dir.figs <- file.path(wd,"Plots")
dir.output <- file.path(wd,"Output")
dir.data <- file.path(wd,"Data")


brood_table <- read.csv(file.path(dir.data,"brood_table.csv"))

brood_table <- brood_table %>% mutate(lnRS=log(R_med/S_med)) %>% 
  filter(BroodYear<2013) %>% rename(Year="BroodYear")

# Fit simple linear Ricker models for each population
prod <- lm(lnRS ~ population*S_med, data = brood_table)
brood_table$lnRS_pred <- predict(prod)
brood_table$res <- brood_table$lnRS - brood_table$lnRS_pred

#add env data for plotting

Env_var <- read.csv(file.path(dir.data,"/Environmental data/Processed/all_env_unstd.csv"))

Env_var <- Env_var %>% rename(population="Population")

Env_var$population <- as.factor(Env_var$population)
levels(Env_var$population)<- list("Lower Mainstem"="Lower Mainstem","White-Donjek"="White Donjek","Middle Mainstem"="Middle Mainstem","Upper Lakes and Mainstem"="Upper Lakes and Mainstem",
   Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")

Env_prod <- left_join(brood_table,Env_var)


ggplot(Env_prod,aes(y=res,x=Pdischarge))+
  geom_point()+geom_smooth()

ggplot(Env_prod,aes(y=WeeklyMigtemp_t0,x=Pdischarge))+
  geom_point()+geom_smooth()

mod1 <- lm(res~Pdischarge*WeeklyMigtemp_t0+I(Pdischarge^2),data=Env_prod)
summary(mod1)

gam1 <- gam(res~s(migtemp_t0),data=Env_prod)
summary(gam1)
plot(gam1)

mod1 <- lm(res~migtemp_t0,data=Env_prod)
summary(mod1)

ggplot(Env_prod,aes(y=res,x=migtemp_t0))+
  geom_point(size=2)+geom_smooth()+xlab("Migration temperature (°C)")+ylab("SR residuals")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Env_prod,aes(y=res,x=migtemp_t0))+
  geom_point()+geom_smooth()

gam1 <- gam(res~s(rearing_prcp),data=Env_prod)
summary(gam1)

mod1 <- lm(res~rearing_prcp,data=Env_prod)
summary(mod1)
r.squaredGLMM(mod1)
visreg(mod1,gg=TRUE,points=list(size=2, pch=19),
       line=list(col="darkgray"),fill=list(fill="lightgray"))




# BRT exploration ---------------------------------------------------------

library(dismo)
set.seed(1000)
Prod.brt <- gbm.step(data=Env_prod,
                   gbm.x=c(13:25),gbm.y=12,
                   family="gaussian",tree.complexity = 5,
                   learning.rate = 0.001,bag.fraction = 0.7)

gbm.plot(Prod.brt,n.plots = 13,write.title=FALSE)


#try simplifying the model
set.seed(1000)
Prod.brt.simp <- gbm.simplify(Prod.brt,n.drops = "auto")

set.seed(1000)
Prod.brt.simplified <- gbm.step(data=Env_prod,gbm.x=Prod.brt.simp$pred.list[[3]],
                               gbm.y= 12,
                               family="gaussian",
                               tree.complexity = 5,learning.rate = 0.005,bag.fraction = 0.7)

gbm.plot(Prod.brt.simplified,n.plots = 6,write.title=FALSE)


Prod.int <- gbm.interactions (Prod.brt.simplified)
Prod.int$rank.list

gbm.perspec(Prod.brt.simplified,3,6,z.range=c(-2,2))
#interaction between break up day and Winter SST

Prod.int <- gbm.interactions (Prod.brt)
Prod.int$rank.list

gbm.perspec(Prod.brt.simplified,3,6,z.range=c(-2,2))




# Mixed model -------------------------------------------------------------

library(lme4)
library(lmerTest)

Env_prod$population <- as.factor(Env_prod$population)

Env_prod$Winter_SST_scale <- scale(Env_prod$Winter_SST)
Env_prod$Summer_SST_scale <- scale(Env_prod$Summer_SST)
Env_prod$rearing_prcp_scale <- scale(Env_prod$rearing_prcp)


mod1 <- lmer(res~rearing_prcp_scale+
               (1|population),data=Env_prod)
summary(mod1)
coef(mod1)

mod1 <- lm(res~rearing_prcp_scale+Winter_SST_scale+
             Summer_SST_scale,data=Env_prod)
summary(mod1)
coef(mod1)

library(sjPlot)
library(lattice)

dotplot(ranef(mod1, condVar=T))

library(visreg)

visreg(mod1)

Env_prod_C <- filter(Env_prod,population=="Carmacks")

prod_model=lm(lnRS ~ S_med,data=Env_prod_C)
summary(prod_model)

prod_res=as.data.frame(prod_model$residuals)

par(mfrow=c(2,1))
acf(prod_res, main="ACF of Residuals")
pacf(prod_res, main = "PACF of Residuals")

Mod3 <- arima(prod_res[,1],order=c(1,0,0))
Mod3
#Residuals Model 3
Mod3.residuals=Mod3$residuals
par(mfrow=c(2,1))
acf(Mod3.residuals, main="ACF of Residuals")
pacf(Mod3.residuals, main = "PACF of Residuals")

Env_prod_C <- Env_prod_C %>% mutate(lnRS1=lag(lnRS,1),resid1=lag(resid,1))
Env_prod_C <- Env_prod_C %>% mutate(resid2=lnRS1+resid1)
Env_prod_C$resid <- prod_model$residuals

prod_model_auto <- lm(lnRS~S_med+resid1,data=Env_prod_C)
summary(prod_model_auto)

auto.residuals=prod_model_auto$residuals
par(mfrow=c(2,1))
acf(auto.residuals, main="ACF of Residuals")
pacf(auto.residuals, main = "PACF of Residuals")

lnRS <- lm(lnRS~lnRS1,data=Env_prod_C)
summary(lnRS)

lnRS_resid <- lnRS$residuals

par(mfrow=c(2,1))
acf(lnRS_resid, main="ACF of Residuals")
pacf(lnRS_resid, main = "PACF of Residuals")

Mod1 <- arima(Env_prod_C[,10],order=c(1,0,0))
Mod1
#Residuals Model 1
Mod1.residuals=Mod1$residuals
par(mfrow=c(2,1))
acf(Mod1.residuals, main="ACF of Residuals")
pacf(Mod1.residuals, main = "PACF of Residuals")
