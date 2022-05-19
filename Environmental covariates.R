#### Processing Environmental Covariates ####

library(tidyverse)
library(lubridate)
library(mgcv)

#Define Workflow Paths ====================================================
wd <- file.path(getwd())
dir.figs <- file.path(wd,"Plots")
dir.output <- file.path(wd,"Output")
dir.data <- file.path(wd,"Data")

##### Ice break up date at Dawson in outmigration year (t = +2) ----------
#data from https://yukonriverbreakup.com/statistics

Ice_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Dawson Ice break up.csv"))

#view trend
ggplot(Ice_int1,aes(y=Julian_day,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Ice breakup date at Dawson")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


#index to outmigration year

Ice_int1 <- Ice_int1 %>% 
  mutate(Year2 = Year-2) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Ice_int2 <- filter(Ice_int1,Year>1984&Year<2013)

#standardize
Ice_int3 <- t(scale(Ice_int2$Julian_day))

#copy rows down
Ice_out <- Ice_int3[rep(seq_len(nrow(Ice_int3)),each=8),]

write.csv(Ice_out,file.path(dir.data,"/Environmental data/Processed/Ice_out.csv"),row.names=FALSE)

####### Migration temperature using water temp data from Emmonak (t = 0)##########################################################################

#Migration temp processing in Calculating Migration Temperature.R


#Max weekly temp
Weekly_mig_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Max_weekly_migration_temp_unstd.csv"))

Weekly_mig_temp_int2 <- filter(Weekly_mig_temp_int1,year<2013)

weekly_migration_temp_plot <- Weekly_mig_temp_int2 %>% 
  mutate(population=fct_relevel(population,"LowerMainstem","White-Donjek","Stewart","Pelly",
                                "Teslin","UpperMainstem","Carmacks","MiddleMainstem")) %>% 
ggplot(aes(y=water_temp,x=population))+
  geom_boxplot()+ylab("Weekly max mig temp (°C)")+xlab("Population")+theme_bw()+
  scale_x_discrete(labels=c("Carmacks"="Carmacks","UpperMainstem"="Upper","LowerMainstem"="Lower",
                     "Pelly"="Pelly","Stewart"="Stewart","White-Donjek"="White",
                     "Teslin"="Teslin","MiddleMainstem"="Middle"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


#standardize
Weekly_mig_temp_int3 <- Weekly_mig_temp_int2
Weekly_mig_temp_int3$water_temp <- scale(Weekly_mig_temp_int3$water_temp)

Weekly_migration_temp_t0 <- Weekly_mig_temp_int3 %>%  
  spread(year,water_temp) %>% select(-population)

write.csv(Weekly_migration_temp_t0,file.path(dir.data,"/Environmental data/Processed/Weekly_max_migration_temp_t0.csv"),row.names=FALSE)




#Mean daily temp
Daily_mig_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Mean_daily_migration_temp_unstd.csv"))

Daily_mig_temp_int2 <- filter(Daily_mig_temp_int1,year<2013)

Daily_migration_temp_plot <- Daily_mig_temp_int2 %>% 
  mutate(population=fct_relevel(population,"LowerMainstem","White-Donjek","Stewart","Pelly",
                                "Teslin","UpperMainstem","Carmacks","MiddleMainstem")) %>% 
ggplot(aes(y=water_temp,x=population))+
  geom_boxplot()+ylab("Daily mig temp (°C)")+xlab("Population")+theme_bw()+
  scale_x_discrete(labels=c("Carmacks"="Carmacks","UpperMainstem"="Upper","LowerMainstem"="Lower",
                     "Pelly"="Pelly","Stewart"="Stewart","White-Donjek"="White",
                     "Teslin"="Teslin","MiddleMainstem"="Middle"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


#standardize
Daily_mig_temp_int3 <- Daily_mig_temp_int2
Daily_mig_temp_int3$water_temp <- scale(Daily_mig_temp_int3$water_temp)

Daily_migration_temp_t0 <- Daily_mig_temp_int3 %>%  
  spread(year,water_temp) %>% select(-population)

write.csv(Daily_migration_temp_t0,file.path(dir.data,"/Environmental data/Processed/Daily_migration_temp_t0.csv"),row.names=FALSE)


#day of run that exceeds 17 degrees


Thres17_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Threshold17.csv"))

Thres17_int2 <- filter(Thres17_int1,year<2013)


#standardize
Thres17_int3 <- Thres17_int2
Thres17_int3$run_day <- scale(Thres17_int3$run_day)

Threshold17 <- Thres17_int3 %>%  select(year,population,run_day) %>% 
  spread(year,run_day) %>% select(-population)

write.csv(Threshold17,file.path(dir.data,"/Environmental data/Processed/Threshold17.csv"),row.names=FALSE)


##### Migration temperature (RETURN INDEX) ------------------------------------

# Migration_temp_returns_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/MigTemp_returns_unformatted.csv"))
# Migration_temp_returns_int1 <- Migration_temp_returns_int1 %>% 
#   select(-water_temp,-prop_4,-prop_5,-prop_6,-prop_7) %>% 
#   rename(water_temp="Mig_temp_returns")
# 
# #standardize
# Migration_temp_returns_int2 <- Migration_temp_returns_int1
# Migration_temp_returns_int2$water_temp <- scale(Migration_temp_returns_int2$water_temp)
# 
# ggplot(Migration_temp_returns_int1,aes(y=water_temp,x=Population))+
#  geom_boxplot()
# 
# Migration_temp_returns <- Migration_temp_returns_int2 %>% 
#   spread(year,water_temp) %>% select(-Population)
# 
# write.csv(Migration_temp_returns,file.path(dir.data,"/Environmental data/Processed/Migration_temp_returns.csv"),row.names=FALSE)
# 


##### Juvenile rearing temperature ------------------------------------------------

#Using Daymet monthly air temperature data processed in ArcGIS Pro to watershed-level means

Temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmonthlymean.csv"))

Temp_int2 <- select(Temp_int1,SUB_DRAINA,StdTime,MEAN)
Temp_int3 <- rename(Temp_int2,Population="SUB_DRAINA",AirTemp_mnth="MEAN")

#pull out year and month columns

Temp_int3$Year <- as.numeric(format(as.Date(Temp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Temp_int3$Month <- as.numeric(format(as.Date(Temp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Temp_int4 <- select(Temp_int3,-StdTime)

#index to t = +1

Temp_int5 <- Temp_int4 %>% 
  mutate(Year2 = Year-1) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Temp_int6 <- filter(Temp_int5,Year>1984&Year<2013)

#calculate annual rearing temperature from June to Aug

Temp_int7 <- filter(Temp_int6,Month==6|Month==7|Month==8) 

Temp_int8 <- Temp_int7 %>% 
  group_by(Population,Year) %>% 
  summarise(rearing_temp=mean(AirTemp_mnth))

#view trend
ggplot(Temp_int8,aes(y=rearing_temp,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Rearing air temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Temp_int8,aes(y=rearing_temp,x=Population))+
  geom_boxplot()+ylab("Rearing temperature")+theme_bw()+
  scale_x_discrete(labels=c("Carmacks"="Carmacks","Upper Lakes and Mainstem"="Upper","Lower Mainstem"="Lower",
                     "Pelly"="Pelly","Stewart"="Stewart","White Donjek"="White",
                     "Teslin"="Teslin","Middle Mainstem"="Middle"))+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))



Temp_int9 <- Temp_int8

Temp_int9$rearing_temp <- scale(Temp_int9$rearing_temp)

#organize into matrix for analyses
Temp_int10 <- Temp_int9 %>% 
  spread(key=Year,value=rearing_temp)
#checked popn order before slicing
rearing_temp <- as.matrix(Temp_int10[2:29])

write.csv(rearing_temp,file.path(dir.data,"/Environmental data/Processed/rearing_temp.csv"),row.names=FALSE)


##### Juvenile rearing precipitation ------------------------------------------------

#Using Daymet monthly precipitation data processed in ArcGIS Pro to watershed-level means

Prcp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pmonthlymean.csv"))

Prcp_int2 <- select(Prcp_int1,SUB_DRAINA,StdTime,MEAN)
Prcp_int3 <- rename(Prcp_int2,Population="SUB_DRAINA",Prcp_mnth="MEAN")

#pull out year and month columns

Prcp_int3$Year <- as.numeric(format(as.Date(Prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Prcp_int3$Month <- as.numeric(format(as.Date(Prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Prcp_int4 <- select(Prcp_int3,-StdTime)

#index to t = +1

Prcp_int5 <- Prcp_int4 %>% 
  mutate(Year2 = Year-1) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Prcp_int6 <- filter(Prcp_int5,Year>1984&Year<2013)

#calculate annual rearing precipitation from June to Aug

Prcp_int7 <- filter(Prcp_int6,Month==6|Month==7|Month==8) 

Prcp_int8 <- Prcp_int7 %>% 
  group_by(Population,Year) %>% 
  summarise(rearing_prcp=sum(Prcp_mnth))

Prcp_int9 <- Prcp_int8

#view trend
ggplot(Prcp_int8,aes(y=rearing_prcp,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Rearing precipitation (mm)")+theme_bw()+
  facet_wrap(~Population)+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Prcp_int8,aes(y=rearing_prcp,x=Population))+
  geom_boxplot()+ylab("Rearing precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Prcp_int8_2000 <- filter(Prcp_int8,Year>1999)
mod1 <- gam(rearing_prcp~s(Year),data=Prcp_int8_2000)
summary(mod1)
plot(mod1)

Prcp_int9$rearing_prcp <- scale(Prcp_int9$rearing_prcp)

#organize into matrix for analyses
Prcp_int10 <- Prcp_int9 %>% 
  spread(key=Year,value=rearing_prcp)
#checked popn order before slicing
rearing_prcp <- as.matrix(Prcp_int10[2:29])

write.csv(rearing_prcp,file.path(dir.data,"/Environmental data/Processed/rearing_prcp.csv"),row.names=FALSE)



##### Spawning and early incubation temperature -------------------------------

#Using Daymet monthly air temperature data processed in ArcGIS Pro to watershed-level means

Spawn_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmonthlymean.csv"))

Spawn_temp_int2 <- select(Spawn_temp_int1,SUB_DRAINA,StdTime,MEAN)
Spawn_temp_int3 <- rename(Spawn_temp_int2,Population="SUB_DRAINA",AirTemp_mnth="MEAN")

#pull out year and month columns

Spawn_temp_int3$Year <- as.numeric(format(as.Date(Spawn_temp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Spawn_temp_int3$Month <- as.numeric(format(as.Date(Spawn_temp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Spawn_temp_int4 <- select(Spawn_temp_int3,-StdTime)

#remove extra years i.e. <1985 and >2012
Spawn_temp_int5 <- filter(Spawn_temp_int4,Year>1984&Year<2013)

#calculate annual spawning and early incubation temperature 
Spawn_temp_int6 <- filter(Spawn_temp_int5,Month==8|Month==9) 

Spawn_temp_int7 <- Spawn_temp_int6 %>% 
  group_by(Population,Year) %>% 
  summarise(spawning_temp=mean(AirTemp_mnth))

#view trend
ggplot(Spawn_temp_int7,aes(y=spawning_temp,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Spawning air temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

mod1 <- gam(spawning_temp~s(Year),data=Spawn_temp_int7)
summary(mod1)

cor.test(Spawn_temp_int7$Year,Spawn_temp_int7$spawning_temp)

Spawn_temp_int8 <- Spawn_temp_int7

Spawn_temp_int8$spawning_temp <- scale(Spawn_temp_int8$spawning_temp)

#organize into matrix for analyses
Spawn_temp_int9 <- Spawn_temp_int8 %>% 
  spread(key=Year,value=spawning_temp)
#checked popn order before slicing
spawning_temp <- as.matrix(Spawn_temp_int9[2:29])

ggplot(Spawn_temp_int7,aes(y=spawning_temp,x=Year))+
  geom_point()+geom_smooth()+facet_wrap(~Population)

write.csv(spawning_temp,file.path(dir.data,"/Environmental data/Processed/spawning_temp.csv"),row.names=FALSE)


##### Spawning and early incubation precipitation -------------------------------

#Using Daymet monthly precipitation data processed in ArcGIS Pro to watershed-level means

Spawn_prcp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pmonthlymean.csv"))

Spawn_prcp_int2 <- select(Spawn_prcp_int1,SUB_DRAINA,StdTime,MEAN)
Spawn_prcp_int3 <- rename(Spawn_prcp_int2,Population="SUB_DRAINA",Prcp_mnth="MEAN")

#pull out year and month columns

Spawn_prcp_int3$Year <- as.numeric(format(as.Date(Spawn_prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Spawn_prcp_int3$Month <- as.numeric(format(as.Date(Spawn_prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Spawn_prcp_int4 <- select(Spawn_prcp_int3,-StdTime)

#remove extra years i.e. <1985 and >2012
Spawn_prcp_int5 <- filter(Spawn_prcp_int4,Year>1984&Year<2013)

#calculate annual spawning and early incubation precipitation

Spawn_prcp_int6 <- filter(Spawn_prcp_int5,Month==8|Month==9|Month==10) 

Spawn_prcp_int7 <- Spawn_prcp_int6 %>% 
  group_by(Population,Year) %>% 
  summarise(spawning_prcp=sum(Prcp_mnth))

#view trend
ggplot(Spawn_prcp_int7,aes(y=spawning_prcp,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Spawning precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Spawn_prcp_int7,aes(y=spawning_prcp,x=Population))+
  geom_boxplot()+ylab("Spawning precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Spawn_prcp_int8 <- Spawn_prcp_int7

Spawn_prcp_int8$spawning_prcp <- scale(Spawn_prcp_int8$spawning_prcp)

#organize into matrix for analyses
Spawn_prcp_int9 <- Spawn_prcp_int8 %>% 
  spread(key=Year,value=spawning_prcp)
#checked popn order before slicing
spawning_prcp <- as.matrix(Spawn_prcp_int9[2:29])

ggplot(Spawn_prcp_int7,aes(y=spawning_prcp,x=Year))+
  geom_point()+geom_smooth()+facet_wrap(~Population)

write.csv(spawning_prcp,file.path(dir.data,"/Environmental data/Processed/spawning_prcp.csv"),row.names=FALSE)


##### Maximum annual snow -----------------------------------------------------

#Influence on following freshet and summer water temperatures
#Try impact on outmigration smolts +2
#rearing fish +1
#spawning t = 0

#Using Daymet snow water equivalent data. Maximum values are averaged over each watershed by year

Swe_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Snow_max.csv"))


Swe_int2 <- select(Swe_int1,SUB_DRAINA,StdTime,MEAN)
Swe_int3 <- rename(Swe_int2,Population="SUB_DRAINA",Swe="MEAN")

#pull out year and month columns

Swe_int3$Year <- as.numeric(format(as.Date(Swe_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))

Swe_int4 <- select(Swe_int3,-StdTime)

ggplot(Swe_int4,aes(y=Swe,x=Population))+
  geom_boxplot()+ylab("Snowpack")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Swe_int4,aes(y=Swe,x=Year))+
  geom_point()+geom_smooth()+ylab("Snowpack")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

mod1 <- lm(Swe~Year,data=Swe_int4)
summary(mod1)

#index to t = 0

Swe_int5 <- Swe_int4 %>% 
  mutate(Year2 = Year-0) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Swe_int6 <- filter(Swe_int5,Year>1984&Year<2013)

Swe_int7 <- Swe_int6

Swe_int7$Swe <- scale(Swe_int7$Swe)

#organize into matrix for analyses
Swe_int8 <- Swe_int7 %>% 
  spread(key=Year,value=Swe)
#checked popn order before slicing
annual_snowpack <- as.matrix(Swe_int8[2:29])

write.csv(annual_snowpack,file.path(dir.data,"/Environmental data/Processed/annual_snowpack.csv"),row.names=FALSE)




##### SST ---------------------------------------------------------------------

# #Bering Sea data - time series not long enough
# MaySST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/MaySST.csv")) %>% select(-"ï..OID_")
# M2SST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/M2SST.csv")) %>% select(-lat,-lon,-depth,-"ï..OID_")
# PISST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/PISST.csv")) %>% select(-lat,-lon,-depth,-"ï..OID_")
# NEBtrawl_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Northeastern Bering Sea Trawl Temps.csv"))
# 
# SST_all <- full_join(MaySST_int1,M2SST_int1)
# SST_all <- full_join(SST_all,PISST_int1)
# 
# SST_all$Year <- as.numeric(format(as.Date(SST_all$time, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# 
# SST_all <- SST_all %>% select(-time) %>% 
#   left_join(NEBtrawl_int1)
# 
# #compare NOAA SST to Pribilof data
# Prib_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pribilof_SST.csv")) %>% 
#   select(StdTime,sst)
# 
# Prib_int1$Year <- as.numeric(format(as.Date(Prib_int1$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# Prib_int1$Month <- as.numeric(format(as.Date(Prib_int1$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))
# 
# Prib_winter_int1 <- Prib_int1 %>% group_by(Year) %>% 
#   filter(Month=="1"|Month=="2"|Month=="3") %>%
#   summarise(Winter_SST=mean(sst))
# 
# SST_all <- left_join(SST_all,Prib_winter_int1)
# 
# ggplot(SST_all,aes(y=oc_PISST_SSTanom,x=Winter_SST))+
#   geom_point(size=2)+geom_smooth(method="lm")
# 
# summary(lm(oc_PISST_SSTanom~Winter_SST,data=SST_all))
# 
# #index to first winter in marine environment
# 
# Prib_winter_int2 <- Prib_winter_int1 %>%
#   mutate(Year2 = Year-3) %>%
#   select(-Year) %>%
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Prib_winter_int3 <- filter(Prib_winter_int2,Year>1984&Year<2013)
# 
# #standardize
# Prib_winter_int4 <- t(scale(Prib_winter_int3$Winter_SST))
# 
# #copy rows down
# Prib_winter <- Prib_winter_int4[rep(seq_len(nrow(Prib_winter_int4)),each=8),]
# 
# write.csv(Prib_winter,file.path(dir.data,"/Environmental data/Processed/SST_Prib_winter.csv"),row.names=FALSE)
# 
# #estimate Pribilof data with NOAA extra years
# 
# Prib_new_winter_int1 <- left_join(Prib_winter_int1,SST_all) %>% select(Year,Winter_SST,oc_PISST_SSTanom)
# 
# Prib_new_winter_int2 <- Prib_new_winter_int1 %>% 
#   mutate(New_Winter_SST=if_else(is.na(oc_PISST_SSTanom),(-8.9870+2.7645*Winter_SST),oc_PISST_SSTanom))
# 
# #index to first winter in marine environment
# 
# Prib_new_winter_int3 <- Prib_new_winter_int2 %>%
#   mutate(Year2 = Year-3) %>%
#   select(-Year) %>%
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Prib_new_winter_int4 <- filter(Prib_new_winter_int3,Year>1984&Year<2013)
# 
# #standardize
# Prib_new_winter_int5 <- t(scale(Prib_new_winter_int4$New_Winter_SST))
# 
# #copy rows down
# Prib_new_winter <- Prib_new_winter_int5[rep(seq_len(nrow(Prib_new_winter_int5)),each=8),]
# 
# write.csv(Prib_new_winter,file.path(dir.data,"/Environmental data/Processed/SST_Prib_new_winter.csv"),row.names=FALSE)


#NOAA SST derived for Northern Bering Sea less than 50m
SST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/SST_monthly.csv")) %>% 
  select(StdTime,MEAN)

SST_int1$Year <- as.numeric(format(as.Date(SST_int1$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
SST_int1$Month <- as.numeric(format(as.Date(SST_int1$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

#derive summer SST

SST_summer_int1 <- SST_int1 %>% group_by(Year) %>% 
  filter(Month=="6"|Month=="7"|Month=="8") %>% 
  summarise(Summer_SST=mean(MEAN))

#view trend
ggplot(SST_summer_int1,aes(y=Summer_SST,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Summer SST")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

mod1 <- gam(Summer_SST~s(Year),data=SST_summer_int1)
summary(mod1)
cor.test(SST_summer_int1$Summer_SST,SST_summer_int1$Year)

#index to first summer in marine environment

SST_summer_int2 <- SST_summer_int1 %>% 
  mutate(Year2 = Year-2) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
SST_summer_int3 <- filter(SST_summer_int2,Year>1984&Year<2013)

#standardize
SST_summer_int4 <- t(scale(SST_summer_int3$Summer_SST))

#copy rows down
SST_summer <- SST_summer_int4[rep(seq_len(nrow(SST_summer_int4)),each=8),]

write.csv(SST_summer,file.path(dir.data,"/Environmental data/Processed/SST_summer.csv"),row.names=FALSE)

#winter SST - use other version below based on southern bering sea
# 
# SST_winter_int1 <- SST_int1 %>% group_by(Year) %>% 
#   filter(Month=="1"|Month=="2"|Month=="3") %>% 
#   summarise(Winter_SST=mean(MEAN))
# 
# #index to first winter in marine environment
# 
# SST_winter_int2 <- SST_winter_int1 %>% 
#   mutate(Year2 = Year-3) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# SST_winter_int3 <- filter(SST_winter_int2,Year>1984&Year<2013)
# 
# #standardize
# SST_winter_int4 <- t(scale(SST_winter_int3$Winter_SST))
# 
# #copy rows down
# SST_winter <- SST_winter_int4[rep(seq_len(nrow(SST_winter_int4)),each=8),]
# 
# write.csv(SST_winter,file.path(dir.data,"/Environmental data/Processed/SST_winter.csv"),row.names=FALSE)

#NOAA data precalculated for northern and southern bering sea
#https://www.fisheries.noaa.gov/resource/data/current-sea-surface-temperatures-eastern-bering-sea-gulf-alaska-and-aleutian-islands

#SSTnew_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/BS-SST-2022-03-25.csv"))

SEBS_SST_int1<- read.csv(file.path(dir.data,"/Environmental data/Raw/Bering_Shelf_SST.csv"))%>% 
  select(StdTime,MEAN)

#winter SST
# SEBS_SST_int1 <- SSTnew_int1 %>% filter(Ecosystem_sub=="Southeastern Bering Sea") %>% 
#   select(date,meansst)

SEBS_SST_int1$Year <- as.numeric(format(as.Date(SEBS_SST_int1$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
SEBS_SST_int1$Month <- as.numeric(format(as.Date(SEBS_SST_int1$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

SEBS_SST_int2 <- SEBS_SST_int1 %>% group_by(Year) %>% 
  filter(Month=="1"|Month=="2"|Month=="3") %>% 
  summarise(Winter_SST=mean(MEAN))

#view trend
ggplot(SEBS_SST_int2,aes(y=Winter_SST,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Winter SST")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


#index to first winter in marine environment

SEBS_SST_int3 <- SEBS_SST_int2 %>% 
  mutate(Year2 = Year-3) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
SEBS_SST_int4 <- filter(SEBS_SST_int3,Year>1984&Year<2013)

#standardize
SEBS_SST_int5 <- t(scale(SEBS_SST_int4$Winter_SST))

#copy rows down
SST_winter_new <- SEBS_SST_int5[rep(seq_len(nrow(SEBS_SST_int5)),each=8),]

write.csv(SST_winter_new,file.path(dir.data,"/Environmental data/Processed/SST_winter_new.csv"),row.names=FALSE)



#summer SST same dataset but not cropped to juvenile salmon distribution area (includes deeper areas)
# 
# NEBS_SST_int1 <- SSTnew_int1 %>% filter(Ecosystem_sub=="Northern Bering Sea") %>% 
#   select(date,meansst)
# 
# NEBS_SST_int1$Year <- as.numeric(format(as.Date(NEBS_SST_int1$date, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# NEBS_SST_int1$Month <- as.numeric(format(as.Date(NEBS_SST_int1$date, format="%Y-%m-%d", "%H:%M:%S"),"%m"))
# 
# 
# NEBS_SST_int2 <- NEBS_SST_int1 %>% group_by(Year) %>% 
#   filter(Month=="6"|Month=="7"|Month=="8") %>% 
#   summarise(summer_SST=mean(meansst))
# 
# #index to first summer in marine environment
# 
# NEBS_SST_int3 <- NEBS_SST_int2 %>% 
#   mutate(Year2 = Year-2) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# NEBS_SST_int4 <- filter(NEBS_SST_int3,Year>1984&Year<2013)
# 
# #standardize
# NEBS_SST_int5 <- t(scale(NEBS_SST_int4$summer_SST))
# 
# #copy rows down
# SST_summer_new <- NEBS_SST_int5[rep(seq_len(nrow(NEBS_SST_int5)),each=8),]
# 
# write.csv(SST_summer_new,file.path(dir.data,"/Environmental data/Processed/SST_summer_new.csv"),row.names=FALSE)


##### cold pool ---------------------------------------------------------------------

#Bering Sea data - time series not long enough
Coldpool_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/ColdPool.csv"))

#index to first summer in marine environment

Coldpool_int2 <- Coldpool_int1 %>% 
  mutate(Year2 = Year-2) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Coldpool_int3 <- filter(Coldpool_int2,Year>1984&Year<2013)

#standardize
Coldpool_int4 <- t(scale(as.numeric(Coldpool_int3$Value)))

#copy rows down
Coldpool <- Coldpool_int4[rep(seq_len(nrow(Coldpool_int4)),each=8),]

write.csv(Coldpool,file.path(dir.data,"/Environmental data/Processed/Coldpool.csv"),row.names=FALSE)

# Checking correlations between vars -----------------------------------------------

#wrangling to long data
Migtemp_t0_for_corrl <- rename(Daily_mig_temp_int2,Year="year",migtemp_t0="water_temp",Population="population")
Migtemp_t0_for_corrl$Population <- as.factor(Migtemp_t0_for_corrl$Population)
levels(Migtemp_t0_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
                                          Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")



# Migtemp_returns_for_corrl <- rename(Migration_temp_returns_int1,Year="year",migtemp_returns="water_temp")
# Migtemp_returns_for_corrl$Population <- as.factor(Migtemp_returns_for_corrl$Population)
# levels(Migtemp_returns_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
#                                           Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")

Ice_for_corrl_int1 <- data.frame(rename(Ice_int2,Breakup_day="Julian_day"))
Ice_for_corrl_int2 <- select(Ice_for_corrl_int1,-Date)
Ice_for_corrl_int2$Year <- as.numeric(Ice_for_corrl_int2$Year)
Ice_for_corrl <- Ice_for_corrl_int2

rearing_precip_for_corrl <- Prcp_int8
rearing_temp_for_corrl <-Temp_int8
annual_snowpack_for_corrl <- Swe_int6
spawning_temp_for_corrl <- Spawn_temp_int7
spawning_prcp_for_corrl <-Spawn_prcp_int7

SST_winter_for_corrl <- SEBS_SST_int4
SST_summer_for_corrl <- SST_summer_int3

master_corrl <- full_join(rearing_temp_for_corrl,rearing_precip_for_corrl,by=c("Year","Population"))
master_corrl <- left_join(master_corrl,Migtemp_t0_for_corrl)
#master_corrl <- left_join(master_corrl,Migtemp_returns_for_corrl)
master_corrl <- left_join(master_corrl,Ice_for_corrl)
master_corrl <- left_join(master_corrl,annual_snowpack_for_corrl)
master_corrl <- left_join(master_corrl,spawning_temp_for_corrl)
master_corrl <- left_join(master_corrl,spawning_prcp_for_corrl)
master_corrl <- left_join(master_corrl,SST_winter_for_corrl)
master_corrl <- left_join(master_corrl,SST_summer_for_corrl)
master_corrl_all <- master_corrl[,3:11]

write.csv(master_corrl,file.path(dir.data,"/Environmental data/Processed/all_env_unstd.csv"),row.names=FALSE)

cor.table <- cor(master_corrl_all)

library(corrplot)
corrplot(cor.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

ggplot(master_corrl,aes(y=rearing_prcp,x=rearing_temp))+
  geom_point()+geom_smooth()+facet_wrap(~Population)

ggplot(master_corrl,aes(y=rearing_prcp,x=spawning_prcp))+
  geom_point()+geom_smooth()+facet_wrap(~Population)

#correlations by population

Carmacks_master <- filter(master_corrl,Population=="Carmacks")
Carmacks_master <- Carmacks_master[,3:12]
Carmacks.table <- cor(Carmacks_master)
corrplot(Carmacks.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Lower_Mainstem_master <- filter(master_corrl,Population=="Lower Mainstem")
Lower_Mainstem_master <- Lower_Mainstem_master[,3:12]
Lower_Mainstem.table <- cor(Lower_Mainstem_master)
corrplot(Lower_Mainstem.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Middle_Mainstem_master <- filter(master_corrl,Population=="Middle Mainstem")
Middle_Mainstem_master <- Middle_Mainstem_master[,3:12]
Middle_Mainstem.table <- cor(Middle_Mainstem_master)
corrplot(Middle_Mainstem.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Pelly_master <- filter(master_corrl,Population=="Pelly")
Pelly_master <- Pelly_master[,3:12]
Pelly.table <- cor(Pelly_master)
corrplot(Pelly.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Stewart_master <- filter(master_corrl,Population=="Stewart")
Stewart_master <- Stewart_master[,3:12]
Stewart.table <- cor(Stewart_master)
corrplot(Stewart.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Teslin_master <- filter(master_corrl,Population=="Teslin")
Teslin_master <- Teslin_master[,3:12]
Teslin.table <- cor(Teslin_master)
corrplot(Teslin.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Upper_Lakes_master <- filter(master_corrl,Population=="Upper Lakes and Mainstem")
Upper_Lakes_master <- Upper_Lakes_master[,3:12]
Upper_Lakes.table <- cor(Upper_Lakes_master)
corrplot(Upper_Lakes.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

White_Donjek_master <- filter(master_corrl,Population=="White Donjek")
White_Donjek_master <- White_Donjek_master[,3:12]
White_Donjek.table <- cor(White_Donjek_master)
corrplot(White_Donjek.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

