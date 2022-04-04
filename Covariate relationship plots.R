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


ggplot(Env_prod,aes(y=res,x=spawning_temp))+
  geom_point()+geom_smooth()

ggplot(Env_prod,aes(y=res,x=migtemp_t0))+
  geom_point()+geom_smooth()

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
