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


Prod_noAR1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/productivity_res_noAR1.csv")) %>% 
  rename(Year="year") %>% select(Year,population,med)

Env_prod <- left_join(Env_prod,Prod_noAR1)

Prod_noAR1_atall <- read.csv(file.path(dir.data,"/Environmental data/Processed/productivity_res_noAR1_atall.csv")) %>%
  rename(Year="year",med2="med") %>% select(Year,population,med2)

Env_prod <- left_join(Env_prod,Prod_noAR1_atall)


ggplot(Env_prod,aes(y=med2,x=med))+
  geom_point()+geom_smooth()


ggplot(Env_prod,aes(y=res,x=med))+
  geom_point()+geom_smooth()

ggplot(Env_prod,aes(y=med,x=WeeklyMigtemp_t0))+
  geom_point()+geom_smooth()



#add productivity year colors
Env_prod$Year

Env_prod <- Env_prod %>% 
  mutate(Productivity=if_else(Year<1992,"Above",
                        if_else(Year<2002,"Average",
                                if_else(Year<2011,"Below","Average"))))

# Env_prod_summary <- Env_prod %>% 
#   group_by(Productivity) %>% 
# summarise(migtemp=mean(WeeklyMigtemp_t0),
#           migtemp_sd=sd(WeeklyMigtemp_t0),
#           rearprcp=mean(rearing_prcp),
#           rearprcp_sd=sd(rearing_prcp))

Env_prod$population <- as.factor(Env_prod$population)

levels(Env_prod$population) <- c("Carmacks","Lower Mainstem","Middle Mainstem",
                                "Pelly","Stewart",
                                 "Teslin","Upper Mainstem","White-Donjek")

Productivity_trend_popn_check <- Env_prod %>% ggplot()+
  geom_line(aes(y=med,x=Year,group=population)) + geom_smooth(aes(y=med,x=Year),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +
                  ylab("Productivity index")+xlab("Year")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend1 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=WeeklyMigtemp_t0,color=Productivity)) + 
  geom_smooth(aes(y=med,x=WeeklyMigtemp_t0),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Migration temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend2 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=rearing_prcp_max,color=Productivity)) + 
  geom_smooth(aes(y=med,x=rearing_prcp_max),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Rearing precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend3 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=GDD5,color=Productivity)) + 
  geom_smooth(aes(y=med,x=GDD5),method="lm") + scale_x_continuous(breaks=c(1000,1300,1600))+
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Juvenile growing degree days")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend4 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Breakup_day,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Breakup_day),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Ice break up day of year")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend5 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Snowpack_mean,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Snowpack_mean),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Snowpack (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend5b <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Snowpack_mean)) + 
  geom_smooth(aes(y=med,x=Snowpack_mean,color=Productivity),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Snowpack")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend6 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=spawning_temp,color=Productivity)) + 
  geom_smooth(aes(y=med,x=spawning_temp),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Spawning temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend7 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=max_spawning_prcp,color=Productivity)) + 
  geom_smooth(aes(y=med,x=max_spawning_prcp),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Spawning precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend8 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Winter_SST,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Winter_SST),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Winter SST (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend9 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Summer_SST,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Summer_SST),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Summer SST (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend10 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Pink_total,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Pink_total),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Pink salmon abundance (millions)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend11 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Chum_total,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Chum_total),method="lm") + 
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Chum salmon abundance (millions)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Env_trend12 <- Env_prod %>% ggplot()+
  geom_point(aes(y=med,x=Wintering_temp,color=Productivity)) + 
  geom_smooth(aes(y=med,x=Wintering_temp),method="lm") + scale_x_continuous(breaks=c(-18,-13,-8))+
  geom_hline(yintercept = 0, lty = "dotted",size=1) +facet_wrap(~population)+
                  ylab("Productivity index")+xlab("Freshwater winter air temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


#looking at differences in covariates for the three time periods

Covar_trend1 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=WeeklyMigtemp_t0,x=Productivity))

Migtemp_by_tperiod <- Env_prod %>% 
  group_by(Productivity) %>% 
  summarise(mean_migtemp=mean(WeeklyMigtemp_t0))

summary(lm(WeeklyMigtemp_t0~Productivity,data=Env_prod))
Migtemp_mod1 <- aov(WeeklyMigtemp_t0~Productivity,data=Env_prod)
TukeyHSD(Migtemp_mod1)

#sig difference between below/above + below/ave where mig temps are sig higher during period of below ave productivity
#yes to higher mig temps lowering productivity during below ave years

Covar_trend2 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=rearing_prcp_max,x=Productivity))

summary(lm(rearing_prcp_max~Productivity,data=Env_prod))
rearprcp_mod1 <- aov(rearing_prcp_max~Productivity,data=Env_prod)
TukeyHSD(rearprcp_mod1)
#NS


Covar_trend3 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Breakup_day,x=Productivity))

summary(lm(Breakup_day~Productivity,data=Env_prod))
Break_mod1 <- aov(Breakup_day~Productivity,data=Env_prod)
TukeyHSD(Break_mod1)

#sig differences bw ave-above and ave-below, but average productivity years have latest break up + highest prod years have earliest
#not a clear trend


Covar_trend4 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Snowpack_mean,x=Productivity))

summary(lm(Snowpack_mean~Productivity,data=Env_prod))
Snow_mod1 <- aov(Snowpack_mean~Productivity,data=Env_prod)
TukeyHSD(Snow_mod1)

Snow_by_tperiod <- Env_prod %>% 
  group_by(Productivity) %>% 
  summarise(mean_snow=mean(Snowpack_mean))

#above average prod years have more snow compared to average and below. is snowpack decreasing then??
#same amount of snow in average and low productivity years

Covar_trend5 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Winter_SST,x=Productivity))

summary(lm(Winter_SST~Productivity,data=Env_prod))
Winter_mod1 <- aov(Winter_SST~Productivity,data=Env_prod)
TukeyHSD(Winter_mod1)

SST_by_tperiod <- Env_prod %>% 
  group_by(Productivity) %>% 
  summarise(mean_SST=mean(Winter_SST))

#all sig different but trend is nonlinear, with low prod years being cold and average prod years having the highest SST
#yes to low prod years being associated with colder temps

Covar_trend6 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=GDD5,x=Productivity))

GDD_by_tperiod <- Env_prod %>% 
  group_by(Productivity) %>% 
  summarise(mean_GDD=mean(GDD5))

summary(lm(GDD5~Productivity,data=Env_prod))
GDD_mod1 <- aov(GDD5~Productivity,data=Env_prod)
TukeyHSD(GDD_mod1)

#below is sig diff from average and above
#yes to low prod years being associated with more degree days

Covar_trend7 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Pink_total,x=Productivity))

summary(lm(Pink_total~Productivity,data=Env_prod))
Pink_mod1 <- aov(Pink_total~Productivity,data=Env_prod)
TukeyHSD(Pink_mod1)

Pink_by_tperiod <- Env_prod %>% 
  group_by(Productivity) %>% 
  summarise(mean_Pink=mean(Pink_total))

#above is sig diff from ave and below (less Pinks)
#below average years had same amount of pinks as average years

Covar_trend8 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=max_spawning_prcp,x=Productivity))

summary(lm(max_spawning_prcp~Productivity,data=Env_prod))

#NS

Covar_trend9 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=spawning_temp,x=Productivity))

summary(lm(spawning_temp~Productivity,data=Env_prod))

#NS


Covar_trend10 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Chum_total,x=Productivity))

summary(lm(Chum_total~Productivity,data=Env_prod))
Chum_mod1 <- aov(Chum_total~Productivity,data=Env_prod)
TukeyHSD(Chum_mod1)

Chum_by_tperiod <- Env_prod %>% 
  group_by(Productivity) %>% 
  summarise(mean_Chum=mean(Chum_total))

#high prod years had sig less chum than ave/above

Covar_trend11 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Summer_SST,x=Productivity))

summary(lm(Summer_SST~Productivity,data=Env_prod))
SummerSST_mod1 <- aov(Summer_SST~Productivity,data=Env_prod)
TukeyHSD(SummerSST_mod1)

#NS

Covar_trend12 <- Env_prod %>% ggplot()+
  geom_boxplot(aes(y=Wintering_temp,x=Productivity))

summary(lm(Wintering_temp~Productivity,data=Env_prod))

#NS


pdf("test.pdf")
grid.arrange(Env_trend1,Env_trend2,Env_trend3,
  nrow=5, ncol=1,
  layout_matrix = rbind(c(1, 1),
                        c(2, 2))
)

dev.off()


library(ggpubr)
g <- ggarrange(arrangeGrob(Env_trend1,Env_trend2,Env_trend3,Env_trend4,Env_trend5,ncol=1),
               newpage=TRUE)

g <- ggarrange(ggarrange(a,b,nrow =2, labels = c("a", "b")),c,
               ncol = 2, labels = c("","c"), hjust=0.7)

library(gridExtra)
g <- grid.arrange(arrangeGrob(Env_trend1,Env_trend2,Env_trend3,Env_trend4,Env_trend5,
                                     ncol=1))



g <- ggarrange(a,ggarrange(b,c,nrow =2, labels = c("b", "c")),
               ncol = 2, labels = c("a"))

jpeg("./04_figures/figures/figure4.jpeg", width = 6, height = 6, units = "in", res = 600)
print(g)
dev.off()

jpeg("./04_figures/figures/figure4.jpeg", width = 6, height = 6, units = "in", res = 600)
print(g)
dev.off()

#correlations 

Env_prod_corrs <- Env_prod %>% select(3,15:18,20:21,24:26,28:31,35)

cor.table <- cor(Env_prod_corrs)

testRes = cor.mtest(Env_prod_corrs, conf.level = 0.95)

library(corrplot)


plot <- corrplot(cor.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="p-value",sig.level=0.01)

corrplot(cor.table$p, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="p-value",sig.level=0.01)

corrplot(cor.table, p.mat = testRes$p, sig.level = 0.05, order = 'hclust', addrect = 2)

ggplot(Env_prod,aes(y=res,x=GDD5))+
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

ggplot(Env_prod,aes(y=res,x=population))+
  geom_boxplot()

gam1 <- gam(res~s(rearing_prcp),data=Env_prod)
summary(gam1)

mod1 <- lm(res~rearing_prcp,data=Env_prod)
summary(mod1)
r.squaredGLMM(mod1)
visreg(mod1,gg=TRUE,points=list(size=2, pch=19),
       line=list(col="darkgray"),fill=list(fill="lightgray"))




# BRT exploration ---------------------------------------------------------

Env_prod$population <- as.factor(Env_prod$population)

library(dismo)
set.seed(1000)
Prod.brt <- gbm.step(data=Env_prod,
                   gbm.x=c(2,15:18,20:21,23:26,28:30),gbm.y=12,
                   family="gaussian",tree.complexity = 5,
                   learning.rate = 0.004,bag.fraction = 0.7)

gbm.plot(Prod.brt,n.plots = 14,write.title=FALSE)


#try simplifying the model
set.seed(1000)
Prod.brt.simp <- gbm.simplify(Prod.brt,n.drops = "auto")

set.seed(1000)
Prod.brt.simplified <- gbm.step(data=Env_prod,gbm.x=Prod.brt.simp$pred.list[[4]],
                               gbm.y= 12,
                               family="gaussian",
                               tree.complexity = 5,learning.rate = 0.004,bag.fraction = 0.7)

gbm.plot(Prod.brt.simplified,n.plots = 12,write.title=FALSE)


Prod.int <- gbm.interactions (Prod.brt.simplified)
Prod.int$rank.list

gbm.perspec(Prod.brt.simplified,3,6,z.range=c(-2,2))
#interaction between break up day and Winter SST

Prod.int <- gbm.interactions (Prod.brt)
Prod.int$rank.list

gbm.perspec(Prod.brt,3,12,z.range=c(-2,2))



#simple regression tree


library(rpart)
library(rpart.plot)

prodtree <- rpart(res~rearing_prcp_max+Breakup_day+Chum_total+Pink_total+rearing_temp+
              Snowpack_mean+WeeklyMigtemp_t0+spawning_temp+max_spawning_prcp+
              Winter_SST+Summer_SST+Wintering_temp+Pdischarge, method="anova", xval =100,
                  data=Env_prod)

printcp(prodtree) # display the results
plotcp(prodtree) # visualize cross-validation results
#min is at 6, whereas 1-SE is 4


summary(prodtree) # detailed summary of splits
# plot tree
rpart.plot(prodtree)

pprodtree<- prune(prodtree, cp= prodtree$cptable[which.min(prodtree$cptable[,"xerror"]),"CP"])
pprodtree<- prune(prodtree, cp= 0.020123)

rpart.plot(pprodtree)
plot(pprodtree)
text(pprodtree)
summary(pprodtree)


# Mixed model -------------------------------------------------------------

library(lme4)
library(lmerTest)
library(MuMIn)
library(nlme)
library(arm)
library(car)

Env_prod$population <- as.factor(Env_prod$population)

Env_prod_std <- Env_prod %>% dplyr::select(13:30) %>% 
  scale()

#used recruitment res from model with no AR1
#productivity <- Env_prod$res
productivity <- Env_prod$med2
population <- Env_prod$population
year <- Env_prod$Year

Env_prod_std <- data.frame(population,year,productivity,Env_prod_std)


#with random intercept
mod1_lme_ar1 <- lme(productivity~rearing_prcp_max+Breakup_day+GDD5+Pink_total+Chum_total+
              Snowpack_mean+WeeklyMigtemp_t0+spawning_temp+max_spawning_prcp+
              Winter_SST+Summer_SST+Wintering_temp, 
               correlation=corAR1(form=~year|population),random= ~ 1|population,
            data=Env_prod_std)
summary(mod1_lme_ar1)
r.squaredGLMM(mod1_lme_ar1)
AIC(mod1_lme_ar1)
hist(mod1_lme_ar1$residuals)
shapiro.test(mod1_lme_ar1$residuals)
qqnorm(mod1_lme_ar1$residuals)

#adding pinks reduces variance explained quite a bit - not sure why because no strong corrls

# mod1_lme_ar1 <- lme(productivity~rearing_prcp+WeeklyMigtemp_t0+Breakup_day+
#               Snowpack_mean+spawning_temp+total_spawning_prcp+
#               Winter_SST+Summer_SST+rearing_temp+Pink_total+
#                Chum_total, random= ~ 1+spawning_temp+rearing_prcp+
#                 total_spawning_prcp+Breakup_day+Winter_SST|population,
#               correlation=corAR1(form=~year|population),
#             data=Env_prod_std)
# summary(mod1_lme_ar1)
# r.squaredGLMM(mod1_lme_ar1)
# AIC(mod1_lme_ar1)
# hist(mod1_lme_ar1$residuals)
# shapiro.test(mod1_lme_ar1$residuals)
# qqnorm(mod1_lme_ar1$residuals)

#without
mod1_lm_ar1 <- gls(productivity~rearing_prcp_max+Breakup_day+GDD5+Pink_total+Chum_total+
              Snowpack_mean+WeeklyMigtemp_t0+spawning_temp+max_spawning_prcp+
              Winter_SST+Summer_SST+Wintering_temp, correlation=corAR1(form=~year|population),
            data=Env_prod_std)
summary(mod1_lm_ar1)
AIC(mod1_lm_ar1)
hist(mod1_lm_ar1$residuals)
shapiro.test(mod1_lm_ar1$residuals)

anova(mod1_lme_ar1,mod1_lm_ar1)

#model without random effects not significantly different and slightly lower AIC


mod1_lme_ar1 <- lme(productivity~rescale(rearing_prcp)+rescale(WeeklyMigtemp_t0)+rescale(Breakup_day)+
              rescale(Snowpack_mean)+rescale(spawning_temp)+rescale(total_spawning_prcp)+
              rescale(Winter_SST)+rescale(Summer_SST)+rescale(rearing_temp)+rescale(Pink_total)+
               rescale(Chum_total), random= ~ 1+WeeklyMigtemp_t0|population,correlation=corAR1(form=~year|population),
            data=Env_prod)
summary(mod1_lme_ar1)
r.squaredGLMM(mod1_lme_ar1)
AIC(mod1_lme_ar1)
hist(mod1_lme_ar1$residuals)
shapiro.test(mod1_lme_ar1$residuals)


mod1_lme_ar1 <- lmer(productivity~rescale(rearing_prcp)+rescale(WeeklyMigtemp_t0)+rescale(Breakup_day)+
              rescale(Snowpack_mean)+rescale(spawning_temp)+rescale(total_spawning_prcp)+
              rescale(Winter_SST)+rescale(Summer_SST)+rescale(rearing_temp)+rescale(Pink_total)+
               rescale(Chum_total)+(1|population),
            data=Env_prod)

mod1_lme_ar1 <- lme(productivity~rescale(rearing_prcp)+rescale(WeeklyMigtemp_t0)+rescale(Breakup_day)+
              rescale(Snowpack_mean)+rescale(spawning_temp)+rescale(total_spawning_prcp)+
              rescale(Winter_SST)+rescale(Summer_SST)+rescale(rearing_temp)+rescale(Pink_total)+
               rescale(Chum_total), random= ~ 1|population,
            data=Env_prod)

ggplot(fortify(mod1_lme_ar1), aes(WeeklyMigtemp_t0, productivity, color=population)) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")


mod1_lm_ar1 <- gls(productivity~rescale(rearing_prcp)+rescale(WeeklyMigtemp_t0)+rescale(Breakup_day)+
              rescale(Snowpack_mean)+rescale(spawning_temp)+rescale(total_spawning_prcp)+
              rescale(Winter_SST)+rescale(Summer_SST)+rescale(rearing_temp)+rescale(Pink_total)+
               rescale(Chum_total), correlation=corAR1(form=~year|population),
            data=Env_prod)
summary(mod1_lm_ar1)
r.squaredGLMM(mod1_lm_ar1)
AIC(mod1_lm_ar1)
hist(mod1_lm_ar1$residuals)
shapiro.test(mod1_lm_ar1$residuals)

anova(mod1_lm_ar1,mod1_lme_ar1)
anova(mod1_lm_ar1,mod1_lm)

Env_prod_discharge <- Env_prod %>% filter(!is.na(Pdischarge))

mod1_lm_disc <- lm(res~population+Winter_SST+Breakup_day+Chum_total+Pdischarge+
                WeeklyMigtemp_t0+Snowpack_mean+Pink_total+Winter_SST+total_spawning_prcp+
                Summer_SST+rearing_temp+rearing_prcp+spawning_temp+I(Pdischarge^2),
            data=Env_prod_discharge)
summary(mod1_lm_disc)
visreg(mod1_lm_disc)
visreg(mod1_lm_disc,"Winter_SST","Breakup_day")
visreg(mod1_lm_disc,"Breakup_day","Winter_SST")

mod1_lm_disc_std <- standardize(mod1_lm_disc,standardize.y=TRUE)
summary(mod1_lm_disc_std)

options(na.action = "na.fail")
stdz.mod1_d <- dredge(mod1_lm_disc_std,rank="AIC")

Top_model <- get.models(stdz.mod1_d,subset=1) [[1]]
summary(Top_model)
vif(Top_model)
hist(Top_model$residuals)
visreg(Top_model)
visreg(Top_model,"z.Breakup_day","z.Winter_SST")
visreg(Top_model,"z.Pink_total","z.Winter_SST")


mod1_lm <- lm(res~Winter_SST+Breakup_day+Chum_total+
                WeeklyMigtemp_t0+Snowpack_mean+Pink_total+Winter_SST+total_spawning_prcp+
                Summer_SST+rearing_temp+rearing_prcp+spawning_temp+Wintering_temp,
            data=Env_prod)
summary(mod1_lm)

mod1_lm_std <- standardize(mod1_lm,standardize.y=TRUE)
summary(mod1_lm_std)

options(na.action = "na.fail")
stdz.mod1_d <- dredge(mod1_lm_std,rank="AIC")

Top_model <- get.models(stdz.mod1_d,subset=1) [[1]]
summary(Top_model)
vif(Top_model)
hist(Top_model$residuals)
visreg(Top_model)
visreg(Top_model,"z.Breakup_day","z.Winter_SST")
visreg(Top_model,"z.Pink_total","z.Winter_SST")
visreg(Top_model,"z.WeeklyMigtemp_t0","z.Snowpack_mean")

r.squaredGLMM(Top_model)
ranef(Top_model)

mod1 <- lme(productivity~rescale(rearing_prcp)+rescale(WeeklyMigtemp_t0)+rescale(Breakup_day)+
              rescale(Snowpack_mean)+rescale(spawning_temp)+rescale(total_spawning_prcp)+
              rescale(Winter_SST)+rescale(Summer_SST)+rescale(rearing_temp)+rescale(Pink_total)+
               rescale(Chum_total), random= ~ 1|population,correlation=corAR1(form=~year|population),
            data=Env_prod)
summary(mod1)

mod1 <- lme(productivity~rearing_prcp+WeeklyMigtemp_t0+Breakup_day+Snowpack_mean+
               spawning_temp+total_spawning_prcp+Winter_SST+Summer_SST+rearing_temp+Pink_total+
               Chum_total, random= ~ rearing_prcp+WeeklyMigtemp_t0+Breakup_day+Snowpack_mean+
               spawning_temp+total_spawning_prcp+Winter_SST+Summer_SST+rearing_temp+Pink_total+
               Chum_total|population,data=Env_prod_std)

mod1 <- lmer(productivity~(rearing_prcp|population)+WeeklyMigtemp_t0+Breakup_day+Snowpack_mean+
               spawning_temp+total_spawning_prcp+Winter_SST+Summer_SST+rearing_temp+Pink_total+
               Chum_total+(1|population),correlation=corAR1,data=Env_prod_std)
summary(mod1)
coef(mod1)
r.squaredGLMM(mod1)

mod1 <- lm(productivity~rearing_prcp+WeeklyMigtemp_t0+Breakup_day+Snowpack_mean+
               spawning_temp+total_spawning_prcp+Winter_SST+Summer_SST+rearing_temp+Pink_total+
               Chum_total,data=Env_prod_std)
summary(mod1)
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
