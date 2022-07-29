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


#view trend
ggplot(Weekly_mig_temp_int2 ,aes(y=water_temp,x=year))+
  geom_point(size=2)+geom_smooth()+ylab("Migration temperature (°C)")+theme_bw()+facet_wrap(~population)+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#standardize
Weekly_mig_temp_int3 <- Weekly_mig_temp_int2
Weekly_mig_temp_int3$water_temp <- scale(Weekly_mig_temp_int3$water_temp)

Weekly_migration_temp_t0 <- Weekly_mig_temp_int3 %>%  
  spread(year,water_temp) %>% select(-population)

write.csv(Weekly_migration_temp_t0,file.path(dir.data,"/Environmental data/Processed/Weekly_max_migration_temp_t0.csv"),row.names=FALSE)

#weighted weekly temp using mean distance to spawning areas in each subbasin

dist_mig_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Yukon_rkm.csv")) %>% 
  dplyr::select(SUB_DRAINA,Site,Distance) %>% rename(population=SUB_DRAINA)

dist_mig_temp_int2 <- dist_mig_temp_int1 %>% group_by(population) %>% 
  summarise(mean_dist=mean(Distance))

summary(dist_mig_temp_int2$mean_dist)

dist_mig_temp_int3 <- dist_mig_temp_int2 %>% mutate(dist_scale=mean_dist/(2158))
dist_mig_temp_int3$population <- as.factor(dist_mig_temp_int3$population)

levels(dist_mig_temp_int3$population)<- list("LowerMainstem"="Lower Mainstem","White-Donjek"="White Donjek","MiddleMainstem"="Middle Mainstem","UpperMainstem"="Upper Lakes and Mainstem",
                                          Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")

dist_mig_temp_int4 <- left_join(Weekly_mig_temp_int2,dist_mig_temp_int3) %>% 
  mutate(water_temp_adj=water_temp*dist_scale)

dist_mig_temp_int4$water_temp_adj <- scale(dist_mig_temp_int4$water_temp_adj)

Weekly_migration_temp_weighted <- dist_mig_temp_int4 %>% dplyr::select(-mean_dist,-dist_scale,-water_temp) %>%   
  spread(year,water_temp_adj) %>% dplyr::select(-population)

write.csv(Weekly_migration_temp_weighted,file.path(dir.data,"/Environmental data/Processed/Weekly_max_migration_temp_weighted.csv"),row.names=FALSE)


# #Mean daily temp
# Daily_mig_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Mean_daily_migration_temp_unstd.csv"))
# 
# Daily_mig_temp_int2 <- filter(Daily_mig_temp_int1,year<2013)
# 
# Daily_migration_temp_plot <- Daily_mig_temp_int2 %>% 
#   mutate(population=fct_relevel(population,"LowerMainstem","White-Donjek","Stewart","Pelly",
#                                 "Teslin","UpperMainstem","Carmacks","MiddleMainstem")) %>% 
# ggplot(aes(y=water_temp,x=population))+
#   geom_boxplot()+ylab("Daily mig temp (°C)")+xlab("Population")+theme_bw()+
#   scale_x_discrete(labels=c("Carmacks"="Carmacks","UpperMainstem"="Upper","LowerMainstem"="Lower",
#                      "Pelly"="Pelly","Stewart"="Stewart","White-Donjek"="White",
#                      "Teslin"="Teslin","MiddleMainstem"="Middle"))+
#   theme(text = element_text(size=25),axis.text=element_text(size=15),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# 
# #standardize
# Daily_mig_temp_int3 <- Daily_mig_temp_int2
# Daily_mig_temp_int3$water_temp <- scale(Daily_mig_temp_int3$water_temp)
# 
# Daily_migration_temp_t0 <- Daily_mig_temp_int3 %>%  
#   spread(year,water_temp) %>% select(-population)
# 
# write.csv(Daily_migration_temp_t0,file.path(dir.data,"/Environmental data/Processed/Daily_migration_temp_t0.csv"),row.names=FALSE)
# 
# 
# #day of year that exceeds 17 degrees
# 
# 
# Thres17_yday_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Threshold17_unstd.csv"))
# 
# Thres17_yday_int2 <- filter(Thres17_yday_int1,year<2013)
# 
# 
# #standardize
# Thres17_yday_int3 <- Thres17_yday_int2
# Thres17_yday_int3$yday17 <- scale(Thres17_yday_int3$yday17)
# 
# Threshold17_yday <- Thres17_yday_int3 %>%  select(year,population,yday17) %>% 
#   spread(year,yday17) %>% select(-population)
# 
# write.csv(Threshold17_yday,file.path(dir.data,"/Environmental data/Processed/Threshold17_yday.csv"),row.names=FALSE)
# 
# 
# 
# #day of run that exceeds 17 degrees
# 
# 
# Thres17_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Threshold17_unstd.csv"))
# 
# Thres17_int2 <- filter(Thres17_int1,year<2013)
# 
# 
# #standardize
# Thres17_int3 <- Thres17_int2
# Thres17_int3$run_day <- scale(Thres17_int3$run_day)
# 
# Threshold17 <- Thres17_int3 %>%  select(year,population,run_day) %>% 
#   spread(year,run_day) %>% select(-population)
# 
# write.csv(Threshold17,file.path(dir.data,"/Environmental data/Processed/Threshold17_runday.csv"),row.names=FALSE)
# 
# #number of days exceeding 17 degrees
# 
# Thres17_count_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Threshold17_numberdays_unstd.csv"))
# 
# Thres17_count_int2 <- filter(Thres17_count_int1,year<2013)
# 
# #standardize
# Thres17_count_int3 <- Thres17_count_int2
# Thres17_count_int3$count <- scale(Thres17_count_int3$count)
# 
# Threshold17_count <- Thres17_count_int3 %>%  dplyr::select(year,population,count) %>% 
#   spread(year,count) %>% dplyr::select(-population)
# 
# write.csv(Threshold17_count,file.path(dir.data,"/Environmental data/Processed/Threshold17_count.csv"),row.names=FALSE)


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



##### Juvenile GDD ------------------------------------------------------------

#growing degree days above five degrees using mean temperatures from daily Daymet data

GDD_int1a <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmax_daily1.csv"))
GDD_int1b <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmax_daily2.csv"))
GDD_int1c <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmin_daily1.csv"))
GDD_int1d <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmin_daily2.csv"))
GDD_int1e <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmin_daily2b.csv"))

GDD_int2a <- bind_rows(GDD_int1a,GDD_int1b) %>% rename(population="VALUE",tmax="MEAN") %>% select(population,StdTime,tmax)

GDD_int2b <- bind_rows(GDD_int1c,GDD_int1d,GDD_int1e) %>% rename(population="VALUE",tmin="MEAN") %>% select(population,StdTime,tmin)

GDD_int3 <- left_join(GDD_int2a,GDD_int2b)

GDD_int3$tmean=rowMeans(GDD_int3[,c("tmin","tmax")])
  
GDD_int4 <- GDD_int3 %>% 
  mutate(Year=as.numeric(format(as.Date(StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y")),
         Month=as.numeric(format(as.Date(StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m")),
         yday=yday(StdTime)) %>% 
  select(population,Year,Month,yday,tmean)
  
GDD_int5 <- GDD_int4 %>% filter(tmean>5) %>% 
  group_by(population,Year) %>% 
  mutate(GDD5=sum(tmean)) %>% 
  select(population,Year,GDD5) %>% distinct()

#view trend
ggplot(GDD_int5,aes(y=GDD5,x=Year))+facet_wrap(~population)+
  geom_point(size=2)+geom_smooth()+ylab("Growing degree days")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#index to t = +1

GDD_int5$Year <- as.numeric(GDD_int5$Year)
GDD_int5 <- as.data.frame(GDD_int5)

GDD_int6 <- GDD_int5 %>% 
  mutate(Year2 = Year-1) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
GDD_int7 <- filter(GDD_int6,Year>1984&Year<2013)

GDD_int8 <- GDD_int7

GDD_int8$GDD5 <- scale(GDD_int8$GDD5)

#organize into matrix for analyses
GDD_int9 <- GDD_int8 %>% 
  spread(key=Year,value=GDD5)

#checked popn order before slicing
rearing_GDD <- as.matrix(GDD_int9[2:29])

write.csv(rearing_GDD,file.path(dir.data,"/Environmental data/Processed/rearing_GDD.csv"),row.names=FALSE)



##### Juvenile rearing precipitation ------------------------------------------------

#Using Daymet monthly precipitation data processed in ArcGIS Pro to watershed-level means

Prcp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pmonthlymean.csv"))

Prcp_int2 <- select(Prcp_int1,SUB_DRAINA,StdTime,MEAN)
Prcp_int3 <- rename(Prcp_int2,Population="SUB_DRAINA",Prcp_mnth="MEAN")

#pull out year and month columns

Prcp_int3$Year <- as.numeric(format(as.Date(Prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Prcp_int3$Month <- as.numeric(format(as.Date(Prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Prcp_int4 <- select(Prcp_int3,-StdTime)

ggplot(Prcp_int4,aes(y=Prcp_mnth,x=as.factor(Month)))+
  geom_boxplot()+ylab("Monthly precipitation (mm)")+theme_bw()+
  facet_wrap(~Population)+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

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



#Using Daymet monthly precipitation data processed in ArcGIS Pro to watershed-level max

Prcp_max_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pmonthlymean.csv"))

Prcp_max_int2 <- select(Prcp_max_int1,SUB_DRAINA,StdTime,MEAN)
Prcp_max_int3 <- rename(Prcp_max_int2,Population="SUB_DRAINA",Prcp_max_mnth="MEAN")

#pull out year and month columns

Prcp_max_int3$Year <- as.numeric(format(as.Date(Prcp_max_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Prcp_max_int3$Month <- as.numeric(format(as.Date(Prcp_max_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Prcp_max_int4 <- select(Prcp_max_int3,-StdTime)

#index to t = +1

Prcp_max_int5 <- Prcp_max_int4 %>% 
  mutate(Year2 = Year-1) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Prcp_max_int6 <- filter(Prcp_max_int5,Year>1984&Year<2013)

#calculate annual rearing precipitation from June to Aug

Prcp_max_int7 <- filter(Prcp_max_int6,Month==6|Month==7|Month==8) 

Prcp_max_int8 <- Prcp_max_int7 %>% 
  group_by(Population,Year) %>% 
  summarise(rearing_prcp_max=max(Prcp_max_mnth))

Prcp_max_int9 <- Prcp_max_int8

#view trend
ggplot(Prcp_max_int8,aes(y=rearing_prcp_max,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Rearing precipitation (mm)")+theme_bw()+
  facet_wrap(~Population)+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Prcp_max_int8,aes(y=rearing_prcp_max,x=Population))+
  geom_boxplot()+ylab("Rearing precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Prcp_max_int8_2000 <- filter(Prcp_max_int8,Year>1999)
mod1 <- gam(rearing_Prcp_max~s(Year),data=Prcp_max_int8_2000)
summary(mod1)
plot(mod1)

Prcp_max_int9$rearing_prcp_max <- scale(Prcp_max_int9$rearing_prcp_max)

#organize into matrix for analyses
Prcp_max_int10 <- Prcp_max_int9 %>% 
  spread(key=Year,value=rearing_prcp_max)
#checked popn order before slicing
rearing_prcp_max <- as.matrix(Prcp_max_int10[2:29])

write.csv(rearing_prcp_max,file.path(dir.data,"/Environmental data/Processed/rearing_prcp_max.csv"),row.names=FALSE)


##### Spawning and early incubation temperature -------------------------------

#Using Daymet monthly air temperature data processed in ArcGIS Pro to watershed-level means
# 
# Spawn_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Tmonthlymean.csv"))
# 
# Spawn_temp_int2 <- select(Spawn_temp_int1,SUB_DRAINA,StdTime,MEAN)
# Spawn_temp_int3 <- rename(Spawn_temp_int2,Population="SUB_DRAINA",AirTemp_mnth="MEAN")
# 
# #pull out year and month columns
# 
# Spawn_temp_int3$Year <- as.numeric(format(as.Date(Spawn_temp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# Spawn_temp_int3$Month <- as.numeric(format(as.Date(Spawn_temp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))
# 
# Spawn_temp_int4 <- select(Spawn_temp_int3,-StdTime)
# 
# #remove extra years i.e. <1985 and >2012
# Spawn_temp_int5 <- filter(Spawn_temp_int4,Year>1984&Year<2013)
# 
# #calculate annual spawning and early incubation temperature 
# Spawn_temp_int6 <- filter(Spawn_temp_int5,Month==8|Month==9) 
# 
# Spawn_temp_int7 <- Spawn_temp_int6 %>% 
#   group_by(Population,Year) %>% 
#   summarise(spawning_temp=mean(AirTemp_mnth))
# 
# #view trend
# ggplot(Spawn_temp_int7,aes(y=spawning_temp,x=Year))+
#   geom_point(size=2)+geom_smooth()+ylab("Spawning air temperature (°C)")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# mod1 <- gam(spawning_temp~s(Year),data=Spawn_temp_int7)
# summary(mod1)
# 
# cor.test(Spawn_temp_int7$Year,Spawn_temp_int7$spawning_temp)
# 
# Spawn_temp_int8 <- Spawn_temp_int7
# 
# Spawn_temp_int8$spawning_temp <- scale(Spawn_temp_int8$spawning_temp)
# 
# #organize into matrix for analyses
# Spawn_temp_int9 <- Spawn_temp_int8 %>% 
#   spread(key=Year,value=spawning_temp)
# #checked popn order before slicing
# spawning_temp <- as.matrix(Spawn_temp_int9[2:29])
# 
# ggplot(Spawn_temp_int7,aes(y=spawning_temp,x=Year))+
#   geom_point()+geom_smooth()+facet_wrap(~Population)
# 
# write.csv(spawning_temp,file.path(dir.data,"/Environmental data/Processed/spawning_temp.csv"),row.names=FALSE)
# 


#using air temperatures at Brown et al spawning sites
#prior data processing done in YRC Water Temperatures project, Climate spawning sites.R
#this is a non-weighted spawning air temperature proxy using the mean of all sites within each subbasin
#daily air temps were a running average of the past 14 days, and then summarized for dates in Aug - Sep (tmean14)

Brown_air_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/spawning_temp_non_weighted.csv"))


#view trend
ggplot(Brown_air_int1,aes(y=tmean14,x=year))+
  geom_point(size=2)+geom_smooth()+ylab("Spawning air temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Brown_air_int2 <- Brown_air_int1 %>% select(subbasin,year,tmean14) %>% filter(year<2013) 

Brown_air_int3 <- Brown_air_int2

Brown_air_int3$tmean14 <- scale(Brown_air_int3$tmean14)

#organize into matrix for analyses
Brown_air_int4 <- Brown_air_int3 %>% 
  spread(key=year,value=tmean14)
#checked popn order before slicing
spawning_brown_air_temp <- as.matrix(Brown_air_int4[2:29])

write.csv(spawning_brown_air_temp,file.path(dir.data,"/Environmental data/Processed/spawning_temp_air.csv"),row.names=FALSE)



#using estimated historical water temperatures at Brown et al spawning sites

Brown_water_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Spawn_sites_predict.csv"))


ggplot(Brown_water_int1,aes(y=WTpredicted,x=yday))+facet_wrap(~subbasin_orig)+
  geom_point(size=2)+geom_smooth()+ylab("Spawning est water temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

#try a few options
# 1 # non-weighted average with all sites
# 2 # major producers only

# 1 # non-weighted spawning average

Brown_water_int2 <- Brown_water_int1 %>% group_by(subbasin_orig,year) %>% 
  summarise(wtmean=mean(WTpredicted))

#view trend
ggplot(Brown_water_int2,aes(y=wtmean,x=year))+
  geom_point(size=2)+geom_smooth()+ylab("Spawning est water temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Brown_water_int3 <- Brown_water_int2 %>% select(subbasin_orig,year,wtmean)

Brown_water_int4 <- Brown_water_int3

Brown_water_int4$wtmean <- scale(Brown_water_int4$wtmean)

#organize into matrix for analyses
Brown_water_int5 <- Brown_water_int4 %>% 
  spread(key=year,value=wtmean)
#checked popn order before slicing
spawning_brown_water_temp <- as.matrix(Brown_water_int5[2:29])

write.csv(spawning_brown_water_temp,file.path(dir.data,"/Environmental data/Processed/spawning_est_water_temp1.csv"),row.names=FALSE)


# 2 # major producers only
# 
# Brown_water_major_int2 <- Brown_water_int1 %>% filter(Prod_level=="major") %>% 
#   group_by(subbasin_orig,year) %>% 
#   summarise(wtmean=mean(WTpredicted))
# 
# #view trend
# ggplot(Brown_water_major_int2,aes(y=wtmean,x=year))+
#   geom_point(size=2)+geom_smooth()+ylab("Spawning est water temperature (°C)")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Brown_water_major_int3 <- Brown_water_major_int2 %>% select(subbasin_orig,year,wtmean)
# 
# Brown_water_major_int4 <- Brown_water_major_int3
# 
# Brown_water_major_int4$wtmean <- scale(Brown_water_major_int4$wtmean)
# 
# #organize into matrix for analyses
# Brown_water_major_int5 <- Brown_water_major_int4 %>% 
#   spread(key=year,value=wtmean)
# #checked popn order before slicing
# spawning_brown_water_temp <- as.matrix(Brown_water_major_int5[2:29])
# 
# write.csv(spawning_brown_water_temp,file.path(dir.data,"/Environmental data/Processed/spawning_temp_major.csv"),row.names=FALSE)
# 

#compare air temp proxy to estimated historical water temps
# 
# compare_spawning_temp_int1 <- Brown_water_int3 %>% rename(subbasin="subbasin_orig")
# 
# compare_spawning_temp <- left_join(compare_spawning_temp_int1,Brown_air_int2)
# 
# compare_spawning_temp$subbasin <- as.factor(compare_spawning_temp$subbasin)
# ggplot(compare_spawning_temp,aes(y=tmean14,x=wtmean))+
#   geom_point(size=2)+geom_smooth()+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# cor.test(compare_spawning_temp$tmean14,compare_spawning_temp$wtmean)

##### Spawning and early incubation precipitation -------------------------------

#Using Daymet monthly precipitation data processed in ArcGIS Pro to watershed-level means
# 
# Spawn_prcp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pmonthlymean.csv"))
# 
# Spawn_prcp_int2 <- select(Spawn_prcp_int1,SUB_DRAINA,StdTime,MEAN)
# Spawn_prcp_int3 <- rename(Spawn_prcp_int2,Population="SUB_DRAINA",Prcp_mnth="MEAN")
# 
# #pull out year and month columns
# 
# Spawn_prcp_int3$Year <- as.numeric(format(as.Date(Spawn_prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# Spawn_prcp_int3$Month <- as.numeric(format(as.Date(Spawn_prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))
# 
# Spawn_prcp_int4 <- select(Spawn_prcp_int3,-StdTime)
# 
# #remove extra years i.e. <1985 and >2012
# Spawn_prcp_int5 <- filter(Spawn_prcp_int4,Year>1984&Year<2013)
# 
# #calculate annual spawning and early incubation precipitation
# 
# Spawn_prcp_int6 <- filter(Spawn_prcp_int5,Month==8|Month==9|Month==10) 
# 
# Spawn_prcp_int7 <- Spawn_prcp_int6 %>% 
#   group_by(Population,Year) %>% 
#   summarise(spawning_prcp=sum(Prcp_mnth))
# 
# #view trend
# ggplot(Spawn_prcp_int7,aes(y=spawning_prcp,x=Year))+
#   geom_point(size=2)+geom_smooth()+ylab("Spawning precipitation (mm)")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# ggplot(Spawn_prcp_int7,aes(y=spawning_prcp,x=Population))+
#   geom_boxplot()+ylab("Spawning precipitation (mm)")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Spawn_prcp_int8 <- Spawn_prcp_int7
# 
# Spawn_prcp_int8$spawning_prcp <- scale(Spawn_prcp_int8$spawning_prcp)
# 
# #organize into matrix for analyses
# Spawn_prcp_int9 <- Spawn_prcp_int8 %>% 
#   spread(key=Year,value=spawning_prcp)
# #checked popn order before slicing
# spawning_prcp <- as.matrix(Spawn_prcp_int9[2:29])
# 
# ggplot(Spawn_prcp_int7,aes(y=spawning_prcp,x=Year))+
#   geom_point()+geom_smooth()+facet_wrap(~Population)
# 
# write.csv(spawning_prcp,file.path(dir.data,"/Environmental data/Processed/spawning_prcp.csv"),row.names=FALSE)
# 



#Using Daymet monthly precipitation data processed in ArcGIS Pro to spawning watershed-level means

Spawn_prcp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Prcpmonthly_spawning.csv"))

Spawn_prcp_int2 <- select(Spawn_prcp_int1,Population,Site_1,StdTime,MEAN)
Spawn_prcp_int3 <- rename(Spawn_prcp_int2,Site="Site_1",Prcp_mnth="MEAN")

#fix Wolf R empty population to Teslin

Spawn_prcp_int3 <- Spawn_prcp_int3 %>% mutate(Population=if_else(Site!="WolfR",
                                                                 Population,"Teslin"))

#pull out year and month columns

Spawn_prcp_int3$Year <- as.numeric(format(as.Date(Spawn_prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Spawn_prcp_int3$Month <- as.numeric(format(as.Date(Spawn_prcp_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Spawn_prcp_int4 <- select(Spawn_prcp_int3,-StdTime)

#remove extra years i.e. <1985 and >2012
Spawn_prcp_int5 <- filter(Spawn_prcp_int4,Year>1984&Year<2013)

#calculate annual spawning and early incubation precipitation

Spawn_prcp_int6 <- filter(Spawn_prcp_int5,Month==8|Month==9|Month==10) 

Spawn_prcp_int7 <- Spawn_prcp_int6 %>% 
  group_by(Population,Year,Month) %>% 
  summarise(mean_monthly_spawning_prcp=mean(Prcp_mnth))

Spawn_prcp_int7b <- Spawn_prcp_int7 %>% 
group_by(Population,Year) %>% 
  summarise(total_spawning_prcp=sum(mean_monthly_spawning_prcp))
  
#view trend


ggplot(Spawn_prcp_int7b,aes(y=spawning_prcp,x=Population))+
  geom_boxplot()+ylab("Spawning precipitation total (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Spawn_prcp_int8 <- Spawn_prcp_int7b

Spawn_prcp_int8$total_spawning_prcp <- scale(Spawn_prcp_int8$total_spawning_prcp)

#organize into matrix for analyses
Spawn_prcp_int9 <- Spawn_prcp_int8 %>% 
  spread(key=Year,value=total_spawning_prcp)
#checked popn order before slicing
spawning_prcp <- as.matrix(Spawn_prcp_int9[2:29])

ggplot(Spawn_prcp_int7b,aes(y=spawning_prcp,x=Year))+
  geom_point()+geom_smooth()+facet_wrap(~Population)


write.csv(spawning_prcp,file.path(dir.data,"/Environmental data/Processed/spawning_prcp_total.csv"),row.names=FALSE)




#Using Daymet monthly precipitation data processed in ArcGIS Pro to spawning watershed-level max monthly prcp

Spawn_prcp_max_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Prcpmonthly_spawning.csv"))

Spawn_prcp_max_int2 <- select(Spawn_prcp_max_int1,Population,Site_1,StdTime,MEAN)
Spawn_prcp_max_int3 <- rename(Spawn_prcp_max_int2,Site="Site_1",Prcp_mnth="MEAN")

#fix Wolf R empty population to Teslin

Spawn_prcp_max_int3 <- Spawn_prcp_max_int3 %>% mutate(Population=if_else(Site!="WolfR",
                                                                 Population,"Teslin"))

#pull out year and month columns

Spawn_prcp_max_int3$Year <- as.numeric(format(as.Date(Spawn_prcp_max_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Spawn_prcp_max_int3$Month <- as.numeric(format(as.Date(Spawn_prcp_max_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Spawn_prcp_max_int4 <- select(Spawn_prcp_max_int3,-StdTime)

#remove extra years i.e. <1985 and >2012
Spawn_prcp_max_int5 <- filter(Spawn_prcp_max_int4,Year>1984&Year<2013)

#calculate annual spawning and early incubation precipitation

Spawn_prcp_max_int6 <- filter(Spawn_prcp_max_int5,Month==8|Month==9|Month==10) 

Spawn_prcp_max_int7 <- Spawn_prcp_max_int6 %>% 
  group_by(Population,Year,Month) %>% 
  summarise(mean_monthly_spawning_prcp=mean(Prcp_mnth))

ggplot(Spawn_prcp_max_int7,aes(y=mean_monthly_spawning_prcp,x=Month))+
  geom_boxplot()+ylab("Spawning precipitation (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Spawn_prcp_max_int7b <- Spawn_prcp_max_int7 %>% 
group_by(Population,Year) %>% 
  summarise(max_spawning_prcp=max(mean_monthly_spawning_prcp))
  
#view trend
ggplot(Spawn_prcp_max_int7b,aes(y=max_spawning_prcp,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Spawning precipitation max (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Spawn_prcp_max_int7b,aes(y=max_spawning_prcp,x=Population))+
  geom_boxplot()+ylab("Spawning precipitation max (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Spawn_prcp_max_int8 <- Spawn_prcp_max_int7b

Spawn_prcp_max_int8$max_spawning_prcp <- scale(Spawn_prcp_max_int8$max_spawning_prcp)

#organize into matrix for analyses
Spawn_prcp_max_int9 <- Spawn_prcp_max_int8 %>% 
  spread(key=Year,value=max_spawning_prcp)
#checked popn order before slicing
spawning_prcp_max <- as.matrix(Spawn_prcp_max_int9[2:29])

write.csv(spawning_prcp_max,file.path(dir.data,"/Environmental data/Processed/spawning_prcp_max.csv"),row.names=FALSE)





#Using Daymet monthly precipitation data processed in ArcGIS Pro to spawning watershed-level means

Spawn_prcp_mean_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Prcpmonthly_spawning.csv"))

Spawn_prcp_mean_int2 <- select(Spawn_prcp_mean_int1,Population,Site_1,StdTime,MEAN)
Spawn_prcp_mean_int3 <- rename(Spawn_prcp_mean_int2,Site="Site_1",Prcp_mnth="MEAN")

#fix Wolf R empty population to Teslin

Spawn_prcp_mean_int3 <- Spawn_prcp_mean_int3 %>% mutate(Population=if_else(Site!="WolfR",
                                                                 Population,"Teslin"))

#pull out year and month columns

Spawn_prcp_mean_int3$Year <- as.numeric(format(as.Date(Spawn_prcp_mean_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
Spawn_prcp_mean_int3$Month <- as.numeric(format(as.Date(Spawn_prcp_mean_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

Spawn_prcp_mean_int4 <- select(Spawn_prcp_mean_int3,-StdTime)

#remove extra years i.e. <1985 and >2012
Spawn_prcp_mean_int5 <- filter(Spawn_prcp_mean_int4,Year>1984&Year<2013)

#calculate annual spawning and early incubation precipitation

Spawn_prcp_mean_int6 <- filter(Spawn_prcp_mean_int5,Month==8|Month==9|Month==10) 

Spawn_prcp_mean_int7 <- Spawn_prcp_mean_int6 %>% 
  group_by(Population,Year,Month) %>% 
  summarise(mean_monthly_spawning_prcp=mean(Prcp_mnth))

Spawn_prcp_mean_int7b <- Spawn_prcp_mean_int7 %>% 
group_by(Population,Year) %>% 
  summarise(mean_spawning_prcp=mean(mean_monthly_spawning_prcp))
  
#view trend
ggplot(Spawn_prcp_mean_int7b,aes(y=mean_spawning_prcp,x=Year))+
  geom_point(size=2)+geom_smooth()+ylab("Spawning precipitation mean (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

ggplot(Spawn_prcp_mean_int7b,aes(y=spawning_prcp,x=Population))+
  geom_boxplot()+ylab("Spawning precipitation mean (mm)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Spawn_prcp_mean_int8 <- Spawn_prcp_mean_int7b

Spawn_prcp_mean_int8$mean_spawning_prcp <- scale(Spawn_prcp_mean_int8$mean_spawning_prcp)

#organize into matrix for analyses
Spawn_prcp_mean_int9 <- Spawn_prcp_mean_int8 %>% 
  spread(key=Year,value=mean_spawning_prcp)
#checked popn order before slicing
spawning_prcp <- as.matrix(Spawn_prcp_mean_int9[2:29])


write.csv(spawning_prcp,file.path(dir.data,"/Environmental data/Processed/spawning_prcp_mean.csv"),row.names=FALSE)



##### Maximum annual snow -----------------------------------------------------

#Influence on following freshet and summer water temperatures
#Try impact on outmigration smolts +2
#rearing fish +1
#spawning t = 0

#Using Daymet snow water equivalent data. Maximum values are averaged over each watershed by year
# 
# Swe_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Snow_max.csv"))
# 
# 
# Swe_int2 <- select(Swe_int1,SUB_DRAINA,StdTime,MEAN)
# Swe_int3 <- rename(Swe_int2,Population="SUB_DRAINA",Swe="MEAN")
# 
# #pull out year and month columns
# 
# Swe_int3$Year <- as.numeric(format(as.Date(Swe_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# 
# Swe_int4 <- select(Swe_int3,-StdTime)
# 
# ggplot(Swe_int4,aes(y=Swe,x=Population))+
#   geom_boxplot()+ylab("Snowpack")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# ggplot(Swe_int4,aes(y=Swe,x=Year))+
#   geom_point()+geom_smooth()+ylab("Snowpack")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# mod1 <- lm(Swe~Year,data=Swe_int4)
# summary(mod1)
# 
# #index to t = 0
# 
# Swe_int5 <- Swe_int4 %>% 
#   mutate(Year2 = Year-0) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Swe_int6 <- filter(Swe_int5,Year>1984&Year<2013)
# 
# Swe_int7 <- Swe_int6
# 
# Swe_int7$Swe <- scale(Swe_int7$Swe)
# 
# #organize into matrix for analyses
# Swe_int8 <- Swe_int7 %>% 
#   spread(key=Year,value=Swe)
# #checked popn order before slicing
# annual_snowpack <- as.matrix(Swe_int8[2:29])
# 
# write.csv(annual_snowpack,file.path(dir.data,"/Environmental data/Processed/annual_snowpack.csv"),row.names=FALSE)
# 
# 
# #Max snowpack by year averaged over spawning watersheds
# 
# Swe_spawn_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Spawning_snow_max.csv"))
# 
# 
# Swe_spawn_int2 <- select(Swe_spawn_int1,Population,Site_1,StdTime,MEAN)
# Swe_spawn_int3 <- rename(Swe_spawn_int2,Site="Site_1",Swe="MEAN")
# 
# #pull out year and month columns
# 
# Swe_spawn_int3$Year <- as.numeric(format(as.Date(Swe_spawn_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# 
# Swe_spawn_int4 <- select(Swe_spawn_int3,-StdTime)
# 
# #summarise to watershed
# 
# Swe_spawn_int4b <- Swe_spawn_int4 %>% group_by(Population,Year) %>% 
#   summarise(Snowpack_mean=mean(Swe))
# 
# 
# #index to t = 0
# 
# Swe_spawn_int5 <- Swe_spawn_int4b %>% 
#   mutate(Year2 = Year-0) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Swe_spawn_int6 <- filter(Swe_spawn_int5,Year>1984&Year<2013)
# 
# Swe_spawn_int7 <- Swe_spawn_int6
# 
# Swe_spawn_int7$Snowpack_mean <- scale(Swe_spawn_int7$Snowpack_mean)
# 
# #organize into matrix for analyses
# Swe_spawn_int8 <- Swe_spawn_int7 %>% 
#   spread(key=Year,value=Snowpack_mean)
# #checked popn order before slicing
# annual_snowpack <- as.matrix(Swe_spawn_int8[2:29])
# 
# write.csv(annual_snowpack,file.path(dir.data,"/Environmental data/Processed/annual_snowpack_spawning.csv"),row.names=FALSE)
# 
# 
# 
# #Max snowpack by year averaged over spawning watersheds
# 
# Swe_spawn_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Spawning_snow_max_sum.csv"))
# 
# 
# Swe_spawn_int2 <- select(Swe_spawn_int1,Population,Site_1,StdTime,SUM)
# Swe_spawn_int3 <- rename(Swe_spawn_int2,Site="Site_1",Swe="SUM")
# 
# #pull out year and month columns
# 
# Swe_spawn_int3$Year <- as.numeric(format(as.Date(Swe_spawn_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# 
# #remove MSYukon site - very high outlier skews data too much
# 
# Swe_spawn_int3 <- Swe_spawn_int3 %>% filter(Site!="MSYukonR")
# 
# Swe_spawn_int4 <- select(Swe_spawn_int3,-StdTime)
# 
# 
# 
# #summarise to watershed
# 
# Swe_spawn_int4b <- Swe_spawn_int4 %>% group_by(Population,Year) %>% 
#   summarise(Snowpack_mean=mean(Swe))
# 
# 
# #index to t = 0
# 
# Swe_spawn_int5 <- Swe_spawn_int4b %>% 
#   mutate(Year2 = Year-0) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Swe_spawn_int6 <- filter(Swe_spawn_int5,Year>1984&Year<2013)
# 
# ggplot(Swe_spawn_int6,aes(y=Snowpack_mean,x=Population))+
#   geom_boxplot()+ylab("Snowpack")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Swe_spawn_int7 <- Swe_spawn_int6
# 
# Swe_spawn_int7$Snowpack_mean <- scale(Swe_spawn_int7$Snowpack_mean)
# 
# #organize into matrix for analyses
# Swe_spawn_int8 <- Swe_spawn_int7 %>% 
#   spread(key=Year,value=Snowpack_mean)
# #checked popn order before slicing
# annual_snowpack <- as.matrix(Swe_spawn_int8[2:29])
# 
# write.csv(annual_snowpack,file.path(dir.data,"/Environmental data/Processed/annual_snowpack_spawning_sum.csv"),row.names=FALSE)
# 
# 

#Snow on Apr1 by year averaged over spawning watersheds

Swe_spawn_Apr_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Spawning_snow_Apr1_mean.csv"))


Swe_spawn_Apr_int2 <- select(Swe_spawn_Apr_int1,Population,Site_1,StdTime,MEAN)
Swe_spawn_Apr_int3 <- rename(Swe_spawn_Apr_int2,Site="Site_1",Swe="MEAN")

#pull out year and month columns

Swe_spawn_Apr_int3$Year <- as.numeric(format(as.Date(Swe_spawn_Apr_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))

#remove MSYukon site - very high outlier skews data too much
# 
# Swe_spawn_Apr_int3 <- Swe_spawn_Apr_int3 %>% filter(Site!="MSYukonR")

Swe_spawn_Apr_int4 <- select(Swe_spawn_Apr_int3,-StdTime)



#summarise to watershed

Swe_spawn_Apr_int4b <- Swe_spawn_Apr_int4 %>% group_by(Population,Year) %>% 
  summarise(Snowpack_mean=mean(Swe))


#index to t = 0

Swe_spawn_Apr_int5 <- Swe_spawn_Apr_int4b %>% 
  mutate(Year2 = Year-0) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Swe_spawn_Apr_int6 <- filter(Swe_spawn_Apr_int5,Year>1984&Year<2013)


Swe_spawn_Apr_int7 <- Swe_spawn_Apr_int6

ggplot(Swe_spawn_Apr_int6,aes(y=Snowpack_mean,x=Year))+facet_wrap(~Population)+
  geom_point()+geom_smooth()+ylab("Snowpack")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Swe_spawn_Apr_int7$Snowpack_mean <- scale(Swe_spawn_Apr_int7$Snowpack_mean)

#organize into matrix for analyses
Swe_spawn_Apr_int8 <- Swe_spawn_Apr_int7 %>% 
  spread(key=Year,value=Snowpack_mean)
#checked popn order before slicing
annual_snowpack <- as.matrix(Swe_spawn_Apr_int8[2:29])

write.csv(annual_snowpack,file.path(dir.data,"/Environmental data/Processed/snowpack_spawning_mean_Apr1.csv"),row.names=FALSE)




#Snowpack on Apr1 by year averaged over spawning watersheds (Sum snowpack)
# 
# Swe_spawn_Apr_sum_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Spawning_snow_Apr1_sum.csv"))
# 
# 
# Swe_spawn_Apr_sum_int2 <- select(Swe_spawn_Apr_sum_int1,Population,Site_1,StdTime,SUM)
# Swe_spawn_Apr_sum_int3 <- rename(Swe_spawn_Apr_sum_int2,Site="Site_1",Swe="SUM")
# 
# #pull out year and month columns
# 
# Swe_spawn_Apr_sum_int3$Year <- as.numeric(format(as.Date(Swe_spawn_Apr_sum_int3$StdTime, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
# 
# #remove MSYukon site - very high outlier skews data too much
# 
# Swe_spawn_Apr_sum_int3 <- Swe_spawn_Apr_sum_int3 %>% filter(Site!="MSYukonR")
# 
# Swe_spawn_Apr_sum_int4 <- select(Swe_spawn_Apr_sum_int3,-StdTime)
# 
# 
# 
# #summarise to watershed
# 
# Swe_spawn_Apr_sum_int4b <- Swe_spawn_Apr_sum_int4 %>% group_by(Population,Year) %>% 
#   summarise(Snowpack_mean=mean(Swe))
# 
# 
# #index to t = 0
# 
# Swe_spawn_Apr_sum_int5 <- Swe_spawn_Apr_sum_int4b %>% 
#   mutate(Year2 = Year-0) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Swe_spawn_Apr_sum_int6 <- filter(Swe_spawn_Apr_sum_int5,Year>1984&Year<2013)
# 
# ggplot(Swe_spawn_Apr_sum_int6,aes(y=Snowpack_mean,x=Population))+
#   geom_boxplot()+ylab("Snowpack")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Swe_spawn_Apr_sum_int7 <- Swe_spawn_Apr_sum_int6
# 
# Swe_spawn_Apr_sum_int7$Snowpack_mean <- scale(Swe_spawn_Apr_sum_int7$Snowpack_mean)
# 
# #organize into matrix for analyses
# Swe_spawn_Apr_sum_int8 <- Swe_spawn_Apr_sum_int7 %>% 
#   spread(key=Year,value=Snowpack_mean)
# #checked popn order before slicing
# annual_snowpack <- as.matrix(Swe_spawn_Apr_sum_int8[2:29])
# 
# write.csv(annual_snowpack,file.path(dir.data,"/Environmental data/Processed/snowpack_spawning_sum_Apr1.csv"),row.names=FALSE)



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


##### Marine competition ------------------------------------------------------

#try indices from Cunningham paper 

#Data from Ruggerone and Irvine 2018 DOI: 10.1002/mcf2.10023

#japanese hatchery chum
# 
# Jap_chum_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/JapChumHatchery.csv"))
# 
# #index to second summer in marine environment
# 
# Jap_chum_int2 <- Jap_chum_int1 %>% 
#   mutate(Year2 = Year-3) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Jap_chum_int3 <- filter(Jap_chum_int2,Year>1984&Year<2013)
# 
# #standardize
# Jap_chum_int4 <- t(scale(as.numeric(Jap_chum_int3$Returns)))
# 
# #copy rows down
# Jap_chum <- Jap_chum_int4[rep(seq_len(nrow(Jap_chum_int4)),each=8),]
# 
# write.csv(Jap_chum,file.path(dir.data,"/Environmental data/Processed/Jap_chum.csv"),row.names=FALSE)

#Pink salmon all competition

Pink_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pink_all.csv"))

#index to second summer in marine environment

Pink_int2 <- Pink_int1 %>% 
  mutate(Year2 = Year-3) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Pink_int3 <- filter(Pink_int2,Year>1984&Year<2013)

Pink_int4 <- Pink_int3 %>% select(Year,Total)

#standardize
Pink_int5 <- t(scale(as.numeric(Pink_int4$Total)))

#copy rows down
Pink <- Pink_int5[rep(seq_len(nrow(Pink_int5)),each=8),]

write.csv(Pink,file.path(dir.data,"/Environmental data/Processed/Pink.csv"),row.names=FALSE)

#Chum salmon all competition

Chum_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Chum_all.csv"))

#index to second summer in marine environment

Chum_int2 <- Chum_int1 %>% 
  mutate(Year2 = Year-3) %>% 
  select(-Year) %>% 
  rename(Year=Year2)

#remove extra years i.e. <1985 and >2012
Chum_int3 <- filter(Chum_int2,Year>1984&Year<2013)

Chum_int4 <- Chum_int3 %>% select(Year,Total)

#standardize
Chum_int5 <- t(scale(as.numeric(Chum_int4$Total)))

#copy rows down
Chum <- Chum_int5[rep(seq_len(nrow(Chum_int5)),each=8),]

write.csv(Chum,file.path(dir.data,"/Environmental data/Processed/Chum.csv"),row.names=FALSE)


#Sockeye salmon all competition
# 
# Sockeye_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Sockeye_all.csv"))
# 
# #index to second summer in marine environment
# 
# Sockeye_int2 <- Sockeye_int1 %>% 
#   mutate(Year2 = Year-3) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Sockeye_int3 <- filter(Sockeye_int2,Year>1984&Year<2013)
# 
# Sockeye_int4 <- Sockeye_int3 %>% select(Year,Total)
# 
# #standardize
# Sockeye_int5 <- t(scale(as.numeric(Sockeye_int4$Total)))
# 
# #copy rows down
# Sockeye <- Sockeye_int5[rep(seq_len(nrow(Sockeye_int5)),each=8),]
# 
# write.csv(Sockeye,file.path(dir.data,"/Environmental data/Processed/Sockeye.csv"),row.names=FALSE)
# 
# 
# #Chum salmon hatchery competition
# 
# Chum_hatch_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Chum_hatchery_all.csv"))
# 
# #index to second summer in marine environment
# 
# Chum_hatch_int2 <- Chum_hatch_int1 %>% 
#   mutate(Year2 = Year-3) %>% 
#   select(-Year) %>% 
#   rename(Year=Year2)
# 
# #remove extra years i.e. <1985 and >2012
# Chum_hatch_int3 <- filter(Chum_hatch_int2,Year>1984&Year<2013)
# 
# Chum_hatch_int4 <- Chum_hatch_int3 %>% select(Year,Total)
# 
# #standardize
# Chum_hatch_int5 <- t(scale(as.numeric(Chum_hatch_int4$Total)))
# 
# #copy rows down
# Chum_hatch <- Chum_hatch_int5[rep(seq_len(nrow(Chum_hatch_int5)),each=8),]
# 
# write.csv(Chum_hatch,file.path(dir.data,"/Environmental data/Processed/Chum_hatch_all.csv"),row.names=FALSE)



##### YR discharge at Eagle ---------------------------------------------------

YRdischarge_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Eagle_YR_discharge_monthly.csv"))


YRdischarge_int2 <- YRdischarge_int1 %>% mutate(year=as.numeric(substr(MM..YYYY,5,9)),
                                                month=as.numeric(substr(MM..YYYY,1,2))) %>% 
  filter(month==7) %>% dplyr::select(Value,year) %>% rename(discharge=Value) %>% 
  filter(year>1984&year<2013)


#standardize
YRdischarge_int3 <- t(scale(as.numeric(YRdischarge_int2$discharge)))

#copy rows down
YRdischarge <- YRdischarge_int3[rep(seq_len(nrow(YRdischarge_int3)),each=8),]

write.csv(YRdischarge,file.path(dir.data,"/Environmental data/Processed/YRDischarge.csv"),row.names=FALSE)



##### YR discharge at Pilot ---------------------------------------------------

Pilotdischarge_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pilot_YR_discharge_monthly.csv"))

Pilotdischarge_int2 <- Pilotdischarge_int1 %>% filter(month_nu==7)%>%
  dplyr::select(year_nu,mean_va) %>% rename(Pdischarge=mean_va,Year=year_nu) %>% 
  filter(Year>1984&Year<2013)

#need to add in missing years as NA

Discharge_years <- as.data.frame(1985:2012)

colnames(Discharge_years) <- "Year"

Pilotdischarge_int4 <- left_join(Discharge_years,Pilotdischarge_int2)

#standardize
Pilotdischarge_int5 <- as.data.frame(scale(as.numeric(Pilotdischarge_int4$Pdischarge)))

#add mean values for missing years

Pilotdischarge_int6 <- Pilotdischarge_int5 %>% mutate(V1=if_else(!is.na(V1),V1,0)) 

Pilotdischarge_int7 <- t(Pilotdischarge_int6)

#copy rows down
Pilotdischarge <- Pilotdischarge_int7[rep(seq_len(nrow(Pilotdischarge_int7)),each=8),]

write.csv(Pilotdischarge,file.path(dir.data,"/Environmental data/Processed/PilotDischarge.csv"),row.names=FALSE)


#Daily

Pilotdischarge_daily_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pilot_YR_discharge_daily.csv")) %>% 
  select(Date,Discharge)

Pilotdischarge_daily_int2 <- Pilotdischarge_daily_int1 %>% 
 mutate(year = as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%Y")),
        month = as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%m")),
        week = as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%W")),
        yday=yday(Date))

Timing_int1 <- read.csv(file.path(dir.data,"MigTiming_intervals.csv")) 

Pilotdischarge_daily_int3 <- left_join(Timing_int1,Pilotdischarge_daily_int2)

ggplot(Pilotdischarge_daily_int3,aes(y=Discharge,x=yday))+
  geom_point()+geom_smooth()+facet_wrap(~year)

#missing years will need estimating: 1997 - 2000

Discharge_stats <- Pilotdischarge_daily_int3 %>% group_by(yday) %>% 
  summarise(mean_discharge=mean(Discharge,na.rm=TRUE))

ggplot(Discharge_stats,aes(y=mean_discharge,x=yday))+
  geom_point()+geom_smooth()

#fill in missing years

Pilotdischarge_daily_int4 <- left_join(Pilotdischarge_daily_int3,Discharge_stats)
Pilotdischarge_daily_int4$Discharge <- as.numeric(Pilotdischarge_daily_int4$Discharge)

Pilotdischarge_daily_int5 <- Pilotdischarge_daily_int4 %>% 
  mutate(Discharge=if_else(year<1997|year>2000,Discharge,mean_discharge))

ggplot(Pilotdischarge_daily_int5,aes(y=Discharge,x=yday))+
  geom_point()+geom_smooth()+facet_wrap(~year)

#interpolate - end of season in 1998

# library(imputeTS)
# 
# Pilotdischarge_daily_int6 <- Pilotdischarge_daily_int5 %>% na.interpolation()
# 
# ggplot(Pilotdischarge_daily_int6,aes(y=Discharge,x=yday))+
#   geom_point()+geom_smooth()+facet_wrap(~year)


Pilotdischarge_daily_int6 <- Pilotdischarge_daily_int5 %>% group_by(year,population) %>% 
  summarise(Discharge=mean(Discharge,na.rm=TRUE))

#remove extra years i.e. <1985 and >2012
Pilotdischarge_daily_int7 <- filter(Pilotdischarge_daily_int6,year>1984&year<2013)


Pilotdischarge_daily_int8 <- Pilotdischarge_daily_int7

Pilotdischarge_daily_int8$Discharge <- scale(Pilotdischarge_daily_int8$Discharge)

#organize into matrix for analyses
Pilotdischarge_daily_int9 <- Pilotdischarge_daily_int8 %>% 
  spread(key=year,value=Discharge)

#checked popn order before slicing
daily_discharge <- as.matrix(Pilotdischarge_daily_int9[2:29])

write.csv(daily_discharge,file.path(dir.data,"/Environmental data/Processed/PilotDischarge_daily.csv"),row.names=FALSE)



# Freshwater winter -------------------------------------------------------

Winter_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/monthly_winter_temps_spawningsites.csv"))

#calculate annual winter temperature for incubating eggs
#Nov and Dec t=0
#Jan - Apr t=1

Winter_temp_int5a <- filter(Winter_temp_int1,month==11|month==12)
Winter_temp_int5b <- filter(Winter_temp_int1,month==1|month==2|month==3|month==4) %>% 
  mutate(year=year-1)

Winter_temp_int5c <- bind_rows(Winter_temp_int5a,Winter_temp_int5b)

Winter_temp_int6 <- Winter_temp_int5c %>%
  group_by(subbasin,year) %>%
  summarise(Wintering_temp=mean(tmean))

#remove extra years i.e. <1985 and >2012
Winter_temp_int7 <- filter(Winter_temp_int6,year>1984&year<2013)

#view trend
ggplot(Winter_temp_int7,aes(y=Wintering_temp,x=year))+facet_wrap(~subbasin)+
  geom_point(size=2)+geom_smooth()+ylab("Wintering air temperature (°C)")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Winter_temp_int8 <- Winter_temp_int7

Winter_temp_int8$Wintering_temp <- scale(Winter_temp_int8$Wintering_temp)

#organize into matrix for analyses
Winter_temp_int9 <- Winter_temp_int8 %>%
  spread(key=year,value=Wintering_temp)
#checked popn order before slicing
Wintering_temp <- as.matrix(Winter_temp_int9[2:29])


write.csv(Wintering_temp,file.path(dir.data,"/Environmental data/Processed/Winter_temp.csv"),row.names=FALSE)

#Temp variation over full winter period
# 
# Winter_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/monthly_winter_temps_spawningsites.csv"))
# 
# #calculate annual winter temperature for incubating eggs
# #Nov and Dec t=0
# #Jan - Apr t=1
# 
# Winter_temp_int5a <- filter(Winter_temp_int1,month==11|month==12)
# Winter_temp_int5b <- filter(Winter_temp_int1,month==1|month==2|month==3|month==4) %>% 
#   mutate(year=year-1)
# 
# Winter_temp_int5c <- bind_rows(Winter_temp_int5a,Winter_temp_int5b)
# 
# Winter_temp_int6 <- Winter_temp_int5c %>%
#   group_by(subbasin,year) %>%
#   summarise(Winter_CV=(sd(tmean)/mean(tmean*-1)*100))
# 
# #remove extra years i.e. <1985 and >2012
# Winter_temp_int7 <- filter(Winter_temp_int6,year>1984&year<2013)
# 
# #view trend
# ggplot(Winter_temp_int7,aes(y=Winter_CV,x=year))+facet_wrap(~subbasin)+
#   geom_point(size=2)+geom_smooth()+ylab("Wintering air temperature (°C)")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Winter_temp_int8 <- Winter_temp_int7
# 
# Winter_temp_int8$Winter_CV <- scale(Winter_temp_int8$Winter_CV)
# 
# #organize into matrix for analyses
# Winter_temp_int9 <- Winter_temp_int8 %>%
#   spread(key=year,value=Winter_CV)
# #checked popn order before slicing
# Winter_CV <- as.matrix(Winter_temp_int9[2:29])
# 
# 
# write.csv(Winter_CV,file.path(dir.data,"/Environmental data/Processed/Winter_temp_CV.csv"),row.names=FALSE)
# 
# 
# 
# 
# #Temp variation over winter shoulder season months
# 
# Winter_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/monthly_winter_temps_spawningsites.csv"))
# 
# #calculate annual winter temperature for incubating eggs
# #Nov and Dec t=0
# #Jan - Apr t=1
# 
# Winter_temp_int5a <- filter(Winter_temp_int1,month==10|month==11)
# Winter_temp_int5b <- filter(Winter_temp_int1,month==3|month==4) %>% 
#   mutate(year=year-1)
# Winter_temp_int5c <- bind_rows(Winter_temp_int5a,Winter_temp_int5b)
# 
# Winter_temp_int6 <- Winter_temp_int5c %>%
#   group_by(subbasin,year) %>%
#   summarise(Winter_CV=(sd(tmean)/mean(tmean*-1)*100))
# 
# #remove extra years i.e. <1985 and >2012
# Winter_temp_int7 <- filter(Winter_temp_int6,year>1984&year<2013)
# 
# #view trend
# ggplot(Winter_temp_int7,aes(y=Winter_CV,x=year))+facet_wrap(~subbasin)+
#   geom_point(size=2)+geom_smooth()+ylab("Wintering air temperature (°C)")+theme_bw()+
#   theme(text = element_text(size=25),axis.text=element_text(size=20),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
#         axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))
# 
# Winter_temp_int8 <- Winter_temp_int7
# 
# Winter_temp_int8$Winter_CV <- scale(Winter_temp_int8$Winter_CV)
# 
# #organize into matrix for analyses
# Winter_temp_int9 <- Winter_temp_int8 %>%
#   spread(key=year,value=Winter_CV)
# #checked popn order before slicing
# Winter_CV <- as.matrix(Winter_temp_int9[2:29])
# 
# write.csv(Winter_CV,file.path(dir.data,"/Environmental data/Processed/Winter_temp_shoulder_CV.csv"),row.names=FALSE)
# 

#winter length

#read in data

Winter_length_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/winter_length.csv"))

#remove extra years i.e. <1985 and >2012
Winter_length_int2 <- filter(Winter_length_int1,year>1984&year<2013) %>% select(-X)

#view trend
ggplot(Winter_length_int2,aes(y=winter_length,x=year))+facet_wrap(~subbasin)+
  geom_point(size=2)+geom_smooth()+ylab("Winter length")+theme_bw()+
  theme(text = element_text(size=25),axis.text=element_text(size=20),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))

Winter_length_int3 <- Winter_length_int2

Winter_length_int3$winter_length <- scale(Winter_length_int3$winter_length)

#organize into matrix for analyses
Winter_length_int4 <- Winter_length_int3 %>%
  spread(key=year,value=winter_length)
#checked popn order before slicing
Winter_length <- as.matrix(Winter_length_int4[2:29])


write.csv(Winter_length,file.path(dir.data,"/Environmental data/Processed/Winter_length.csv"),row.names=FALSE)

# Checking correlations between vars -----------------------------------------------

#wrangling to long data
# DailyMigtemp_t0_for_corrl <- rename(Daily_mig_temp_int2,Year="year",DailyMigtemp_t0="water_temp",Population="population")
# DailyMigtemp_t0_for_corrl$Population <- as.factor(DailyMigtemp_t0_for_corrl$Population)
# levels(DailyMigtemp_t0_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
#                                           Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")


WeeklyMigtemp_t0_for_corrl <- rename(Weekly_mig_temp_int2,Year="year",WeeklyMigtemp_t0="water_temp",Population="population")
WeeklyMigtemp_t0_for_corrl$Population <- as.factor(WeeklyMigtemp_t0_for_corrl$Population)
levels(WeeklyMigtemp_t0_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
                                          Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")


# Thres_yday_for_corrl <- rename(Thres17_yday_int2,Year="year",Population="population")
# Thres_yday_for_corrl$Population <- as.factor(Thres_yday_for_corrl$Population)
# levels(Thres_yday_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
#                                           Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")
# 
# Thres_count_for_corrl <- rename(Thres17_count_int2,Year="year",Population="population")
# Thres_count_for_corrl$Population <- as.factor(Thres_count_for_corrl$Population)
# levels(Thres_count_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
#                                           Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")




# Migtemp_returns_for_corrl <- rename(Migration_temp_returns_int1,Year="year",migtemp_returns="water_temp")
# Migtemp_returns_for_corrl$Population <- as.factor(Migtemp_returns_for_corrl$Population)
# levels(Migtemp_returns_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
#                                           Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")

Ice_for_corrl_int1 <- data.frame(rename(Ice_int2,Breakup_day="Julian_day"))
Ice_for_corrl_int2 <- select(Ice_for_corrl_int1,-Date)
Ice_for_corrl_int2$Year <- as.numeric(Ice_for_corrl_int2$Year)
Ice_for_corrl <- Ice_for_corrl_int2

rearing_precip_for_corrl <- Prcp_int8
rearing_max_precip_for_corrl <- Prcp_max_int8
rearing_temp_for_corrl <-Temp_int8
rearing_GDD_for_corrl <- GDD_int7 %>% rename(Population=population)
annual_snowpack_for_corrl <- Swe_spawn_Apr_int6
spawning_temp_for_corrl <- Brown_water_int3 %>% rename(Population=subbasin_orig,
                                                       Year=year,spawning_temp=wtmean)
spawning_prcp_for_corrl <-Spawn_prcp_int7b
spawning_max_prcp_for_corrl <-Spawn_prcp_max_int7b

SST_winter_for_corrl <- SEBS_SST_int4
SST_summer_for_corrl <- SST_summer_int3
Discharge_for_corrl <- YRdischarge_int2 %>% rename(Year=year)
Pilot_Discharge_for_corrl <- Pilotdischarge_daily_int7 %>% 
  rename(Year=year,Population=population,Pdischarge=Discharge)
Pilot_Discharge_for_corrl$Population <-  as.factor(Pilot_Discharge_for_corrl$Population)
levels(Pilot_Discharge_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
                                          Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")


PinkComp_for_corrl <- Pink_int4 %>% rename(Pink_total=Total)
ChumComp_for_corrl <- Chum_int4 %>% rename(Chum_total=Total)
FW_winter_for_corrl <- Winter_temp_int7%>% rename(Year=year,Population=subbasin)
FW_winter_length_for_corrl <- Winter_length_int2%>% rename(Year=year,Population=subbasin)

master_corrl <- full_join(rearing_temp_for_corrl,rearing_precip_for_corrl,by=c("Year","Population"))
#master_corrl <- left_join(master_corrl,DailyMigtemp_t0_for_corrl)
master_corrl <- left_join(master_corrl,WeeklyMigtemp_t0_for_corrl)
#master_corrl <- left_join(master_corrl,Thres_yday_for_corrl)
#master_corrl <- left_join(master_corrl,Thres_count_for_corrl)
#master_corrl <- left_join(master_corrl,Migtemp_returns_for_corrl)
master_corrl <- left_join(master_corrl,Ice_for_corrl)
master_corrl <- left_join(master_corrl,annual_snowpack_for_corrl)
master_corrl <- left_join(master_corrl,spawning_temp_for_corrl)
master_corrl <- left_join(master_corrl,spawning_prcp_for_corrl)
master_corrl <- left_join(master_corrl,SST_winter_for_corrl)
master_corrl <- left_join(master_corrl,SST_summer_for_corrl)
master_corrl <- left_join(master_corrl,Discharge_for_corrl)
master_corrl <- left_join(master_corrl,Pilot_Discharge_for_corrl)
master_corrl <- left_join(master_corrl,PinkComp_for_corrl)
master_corrl <- left_join(master_corrl,ChumComp_for_corrl)
master_corrl <- left_join(master_corrl,FW_winter_for_corrl)
master_corrl <- left_join(master_corrl,FW_winter_length_for_corrl)
master_corrl <- left_join(master_corrl,spawning_max_prcp_for_corrl)
master_corrl <- left_join(master_corrl,rearing_max_precip_for_corrl)
master_corrl <- left_join(master_corrl,rearing_GDD_for_corrl)
master_corrl_all <- master_corrl[,c(3:20)]

write.csv(master_corrl,file.path(dir.data,"/Environmental data/Processed/all_env_unstd.csv"),row.names=FALSE)

cor.table <- cor(master_corrl_all)

library(corrplot)
corrplot(cor.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

ggplot(master_corrl,aes(y=Wintering_temp,x=winter_length))+
  geom_point()+geom_smooth()+facet_wrap(~Population)

ggplot(master_corrl,aes(y=rearing_prcp,x=spawning_prcp))+
  geom_point()+geom_smooth()+facet_wrap(~Population)

#correlations by population

Carmacks_master <- filter(master_corrl,Population=="Carmacks")
Carmacks_master <- Carmacks_master[,c(3:12,14:16)]
Carmacks.table <- cor(Carmacks_master)
corrplot(Carmacks.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Lower_Mainstem_master <- filter(master_corrl,Population=="Lower Mainstem")
Lower_Mainstem_master <- Lower_Mainstem_master[,c(3:12,14:16)]
Lower_Mainstem.table <- cor(Lower_Mainstem_master)
corrplot(Lower_Mainstem.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Middle_Mainstem_master <- filter(master_corrl,Population=="Middle Mainstem")
Middle_Mainstem_master <- Middle_Mainstem_master[,c(3:12,14:16)]
Middle_Mainstem.table <- cor(Middle_Mainstem_master)
corrplot(Middle_Mainstem.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Pelly_master <- filter(master_corrl,Population=="Pelly")
Pelly_master <- Pelly_master[,c(3:12,14:16)]
Pelly.table <- cor(Pelly_master)
corrplot(Pelly.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Stewart_master <- filter(master_corrl,Population=="Stewart")
Stewart_master <- Stewart_master[,c(3:12,14:16)]
Stewart.table <- cor(Stewart_master)
corrplot(Stewart.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Teslin_master <- filter(master_corrl,Population=="Teslin")
Teslin_master <- Teslin_master[,c(3:12,14:16)]
Teslin.table <- cor(Teslin_master)
corrplot(Teslin.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Upper_Lakes_master <- filter(master_corrl,Population=="Upper Lakes and Mainstem")
Upper_Lakes_master <- Upper_Lakes_master[,c(3:12,14:16)]
Upper_Lakes.table <- cor(Upper_Lakes_master)
corrplot(Upper_Lakes.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

White_Donjek_master <- filter(master_corrl,Population=="White Donjek")
White_Donjek_master <- White_Donjek_master[,c(3:12,14:16)]
White_Donjek.table <- cor(White_Donjek_master)
corrplot(White_Donjek.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

