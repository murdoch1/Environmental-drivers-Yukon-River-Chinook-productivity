#### Processing Environmental Covariates ####

library(tidyverse)
library(lubridate)

#Define Workflow Paths ====================================================
wd <- file.path(getwd())
dir.figs <- file.path(wd,"Plots")
dir.output <- file.path(wd,"Output")
dir.data <- file.path(wd,"Data")

##### Ice break up date at Dawson in outmigration year (t = +2) ----------
#data from https://yukonriverbreakup.com/statistics

Ice_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Dawson Ice break up.csv"))

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

#data from https://www.adfg.alaska.gov/CF_R3/external/sites/aykdbms_website/DataSelection.aspx
#search for Lower Yukon Test Fishing data

Emmonak_water_int1a <- read.csv(file.path(dir.data,"/Environmental data/Raw/Emmonak Water Temp Data1.csv"))
Emmonak_water_int1b <- read.csv(file.path(dir.data,"/Environmental data/Raw/Emmonak Water Temp Data2.csv"))
Emmonak_water_int1a$Instrument.Site <- as.character(Emmonak_water_int1a$Instrument.Site)
Emmonak_water_int1 <- bind_rows(Emmonak_water_int1a,Emmonak_water_int1b)

#sampling locations across year change, look at summary

Emmonak_water_summary1  <- Emmonak_water_int1 %>% 
  group_by(Project.Year,Location) %>% count(Location)

#two years missing middle mouth and one missing big eddy

ggplot(Emmonak_water_summary1,aes(fill=Location,y=n,x=Project.Year))+
  geom_bar(stat="identity",position="dodge")+scale_y_continuous(trans="log2")

#sample timing differs too. before 2004 temperatures were only taken at 8 am and 8 pm
#how many have 8 am and 8 pm?

Emmonak_water_int2 <- Emmonak_water_int1 %>% 
  select(-Time) %>% 
  rename(Date_Time=Date) %>% 
  separate(Date_Time, into = c("Date", "Time"), sep = " ", remove = FALSE) %>%
  mutate(Date = lubridate::as_date(Date, format = "%Y-%m-%d"),
         Time = hms::as_hms(str_c(Time, ":00")))
#note missing date and time data just for one entry in 2014 gives warning message

Emmonak_water_summary2  <- Emmonak_water_int2 %>% 
  group_by(Project.Year,Location) %>% count(Time)

#examining data before 2004
Emmonak_water_summary3 <- filter(Emmonak_water_summary2,Project.Year<2004)

ggplot(Emmonak_water_summary3,aes(fill=Location,y=n,x=Time))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#most consistent data collection is for Big Eddy data at 8 AM. Will be missing for three years:
#1990 was only middle mouth, 1995 not sampled at 8 AM, 2003 also only sampled at MM

#detour to examine differences for later replacement

#1. Look at significant differences between AM and PM samples for Big Eddy and adjust for missing morning samples in 1995
Emmonak_water_int2$Time <- as.character(Emmonak_water_int2$Time)
AM_PM_int1 <- Emmonak_water_int2 %>%
  filter((Time=="08:00:00"|Time=="08:01:00"|Time=="20:00:00")&Location=="Big Eddy") %>% 
  select(Date,Time,Temperature..avg.) %>% 
  mutate(Time = as.factor(case_when(Time=="20:00:00"  ~ "Evening",
                             Time %in% c("08:00:00","08:01:00") ~ "Morning"))) 

AM_PM_int2 <- AM_PM_int1 %>% 
  group_by(Date,Time) %>% 
  summarize(water_temp = mean(Temperature..avg.))

AM_PM_compare <- AM_PM_int2 %>% spread(Time,water_temp)

#visualize differences
ggplot(AM_PM_compare,aes(y=Evening,x=Morning))+
  geom_point()+geom_smooth()

ggplot(AM_PM_int2,aes(y=water_temp,x=Time))+
  geom_boxplot()

#model differences
temp_time <- summary(lm(Morning~Evening,data=AM_PM_compare))
temp_time

#morning temps can be calculated using the linear model y = 0.980783*x + 0.021345

#2. Look at significant differences between Big Eddy and Middle Mouth samples for years with missing Big Eddy

BE_MM_compare_int1 <- Emmonak_water_int2 %>% 
  filter(Time=="08:00:00"|Time=="08:01:00") %>% 
  select(Location,Date,Temperature..avg.) 

BE_MM_compare_int2 <- BE_MM_compare_int1 %>% 
  group_by(Location,Date) %>% 
  summarize(water_temp = mean(Temperature..avg.)) 

BE_MM_compare <- BE_MM_compare_int2 %>% spread(Location,water_temp) %>% 
  rename(Middle_Mouth="Middle Mouth",Big_Eddy="Big Eddy")

#visualize differences
ggplot(BE_MM_compare,aes(y=Middle_Mouth,x=Big_Eddy))+
  geom_point()+geom_smooth(method = "lm")

ggplot(BE_MM_compare_int2,aes(y=water_temp,x=Location))+
  geom_boxplot()

#model differences
temp_location <- summary(lm(Big_Eddy~Middle_Mouth,data=BE_MM_compare))
temp_location

#Big Eddy missing temps can be calculated using the linear model y = 0.80138*x + 3.02411

#Continue data processing
Emmonak_water_int3 <- filter(Emmonak_water_int2,Time=="08:00:00"|Time=="08:01:00")

Emmonak_water_summary4  <- Emmonak_water_int3 %>% 
  group_by(Project.Year,Location) %>% count(Time)

ggplot(Emmonak_water_summary4,aes(fill=Location,y=n,x=Location))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#extract Big Eddy
Emmonak_water_int4 <- filter(Emmonak_water_int3,Location=="Big Eddy")

#extract MM for 1990 and 2003 and adjust based on linear model
Emmonak_water_int5 <- Emmonak_water_int3 %>% filter(Location=="Middle Mouth"&Project.Year=="1990"|
    Location=="Middle Mouth"&Project.Year=="2003") %>% 
    mutate(Temperature..avg.=Temperature..avg.*0.80138+3.02411)

#extract 8 PM for Big Eddy in 1995 and adjust based on linear model
Emmonak_water_int6 <- Emmonak_water_int2 %>% 
  filter(Time=="20:00:00"&Location=="Big Eddy"&Project.Year=="1995") %>% 
  mutate(Temperature..avg.=Temperature..avg.*0.980783+0.021345)


#bind dataframes together
Emmonak_water_int7 <- bind_rows(Emmonak_water_int4,Emmonak_water_int5,Emmonak_water_int6)

#add julian date
Emmonak_water_int8 <- mutate(Emmonak_water_int7,julian_date=yday(Date))

ggplot(Emmonak_water_int8,aes(y=Temperature..avg.,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

Emmonak_water_int9 <- filter(Emmonak_water_int8,Temperature..avg.!="NA")

Emmonak_water_summary5  <- Emmonak_water_int9 %>% 
  group_by(Project.Year) %>% summarise(
    count = n(),
    min = min(julian_date),
    max=max(julian_date))




#adjust migration interval by population and year using derived passage timing from population diversity paper

Timing_int1 <- read.csv(file.path(dir.data,"MigTiming_intervals.csv"))

#join to temperature data
Emmonak_water_int10 <- Emmonak_water_int9 %>% 
  select(Project.Year,Temperature..avg.,julian_date) %>% 
  rename(year=Project.Year,yday=julian_date,water_temp=Temperature..avg.)

Emmonak_water_int11 <- left_join(Timing_int1,Emmonak_water_int10)

#look at years with missing data
Mig_NAs_summary  <- Emmonak_water_int11 %>% 
  group_by(year,population) %>% summarise(
    count = n(),
    min = min(yday),
    max=max(yday),
    minWT=min(water_temp,na.rm=TRUE),
    meanWT=mean(water_temp,na.rm=TRUE),
    maxWT=max(water_temp,na.rm=TRUE),
    NA_col=(sum(is.na(water_temp)))) %>% 
  mutate(prop_NA=NA_col/count*100)

Mig_NAs_summary_int1 <- Mig_NAs_summary %>% select(year,population,prop_NA)

Emmonak_water_int12 <- left_join(Emmonak_water_int11,Mig_NAs_summary_int1) %>% 
  group_by(year,population) %>% summarise(count = n(),
                                     water_temp=mean(water_temp,na.rm=TRUE),
                                     prop_NA=mean(prop_NA))

#for estimating population means - only use years with greater than 75% data availability
Mig_temp_bypopn <- left_join(Emmonak_water_int11,Mig_NAs_summary_int1) %>% 
  filter(prop_NA<25) %>% 
  group_by(population) %>% summarise(meanWT=mean(water_temp,na.rm=TRUE))

#replace years with >50% NA's with population level means

Emmonak_water_int13 <- left_join(Emmonak_water_int12,Mig_temp_bypopn) %>% 
  mutate(water_temp=if_else(prop_NA>50,meanWT,water_temp)) %>% 
  select(-meanWT,-prop_NA,-count)


#data for calculating return index
write.csv(Emmonak_water_int13,file.path(dir.data,"/Environmental data/Processed/Migration_temp_unstd.csv"),row.names=FALSE)

ggplot(Emmonak_water_int13,aes(y=water_temp,x=population,))+
 geom_boxplot()

#standardize
Emmonak_water_int14 <- Emmonak_water_int13
Emmonak_water_int14$water_temp <- scale(Emmonak_water_int14$water_temp)

Migration_temp_t0 <- Emmonak_water_int14 %>% filter(year<2013) %>%   
  spread(year,water_temp) %>% select(-population)

write.csv(Migration_temp_t0,file.path(dir.data,"/Environmental data/Processed/Migration_temp_t0.csv"),row.names=FALSE)



##### Migration temperature (RETURN INDEX) ------------------------------------

Migration_temp_returns_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/MigTemp_returns_unformatted.csv"))
Migration_temp_returns_int1 <- Migration_temp_returns_int1 %>% 
  select(-water_temp,-prop_4,-prop_5,-prop_6,-prop_7) %>% 
  rename(water_temp="Mig_temp_returns")

#standardize
Migration_temp_returns_int2 <- Migration_temp_returns_int1
Migration_temp_returns_int2$water_temp <- scale(Migration_temp_returns_int2$water_temp)

ggplot(Migration_temp_returns_int1,aes(y=water_temp,x=Population,))+
 geom_boxplot()

Migration_temp_returns <- Migration_temp_returns_int2 %>% 
  spread(year,water_temp) %>% select(-Population)

write.csv(Migration_temp_returns,file.path(dir.data,"/Environmental data/Processed/Migration_temp_returns.csv"),row.names=FALSE)



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

Spawn_prcp_int6 <- filter(Spawn_prcp_int5,Month==8|Month==9) 

Spawn_prcp_int7 <- Spawn_prcp_int6 %>% 
  group_by(Population,Year) %>% 
  summarise(spawning_prcp=sum(Prcp_mnth))

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

#index to t = +2

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


# Checking correlations between vars -----------------------------------------------

#wrangling to long data
Migtemp_t0_for_corrl <- rename(Emmonak_water_int13,Year="year",migtemp_t0="water_temp",Population="population")
Migtemp_t0_for_corrl$Population <- as.factor(Migtemp_t0_for_corrl$Population)
levels(Migtemp_t0_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
                                          Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")



Migtemp_returns_for_corrl <- rename(Migration_temp_returns_int1,Year="year",migtemp_returns="water_temp")
Migtemp_returns_for_corrl$Population <- as.factor(Migtemp_returns_for_corrl$Population)
levels(Migtemp_returns_for_corrl$Population)<- list("Lower Mainstem"="LowerMainstem","White Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
                                          Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")

Ice_for_corrl_int1 <- data.frame(rename(Ice_int2,Breakup_day="Julian_day"))
Ice_for_corrl_int2 <- select(Ice_for_corrl_int1,-Date)
Ice_for_corrl_int2$Year <- as.numeric(Ice_for_corrl_int2$Year)
Ice_for_corrl <- Ice_for_corrl_int2

rearing_precip_for_corrl <- Prcp_int8
rearing_temp_for_corrl <-Temp_int8
annual_snowpack_for_corrl <- Swe_int6
spawning_temp_for_corrl <- Spawn_temp_int7
spawning_prcp_for_corrl <-Spawn_prcp_int7

master_corrl <- full_join(rearing_temp_for_corrl,rearing_precip_for_corrl,by=c("Year","Population"))
master_corrl <- left_join(master_corrl,Migtemp_t0_for_corrl)
master_corrl <- left_join(master_corrl,Migtemp_returns_for_corrl)
master_corrl <- left_join(master_corrl,Ice_for_corrl)
master_corrl <- left_join(master_corrl,annual_snowpack_for_corrl)
master_corrl <- left_join(master_corrl,spawning_temp_for_corrl)
master_corrl <- left_join(master_corrl,spawning_prcp_for_corrl)
master_corrl_all <- master_corrl[,3:10]

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
Carmacks_master <- Carmacks_master[,3:10]
Carmacks.table <- cor(Carmacks_master)
corrplot(Carmacks.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Lower_Mainstem_master <- filter(master_corrl,Population=="Lower Mainstem")
Lower_Mainstem_master <- Lower_Mainstem_master[,3:10]
Lower_Mainstem.table <- cor(Lower_Mainstem_master)
corrplot(Lower_Mainstem.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Middle_Mainstem_master <- filter(master_corrl,Population=="Middle Mainstem")
Middle_Mainstem_master <- Middle_Mainstem_master[,3:10]
Middle_Mainstem.table <- cor(Middle_Mainstem_master)
corrplot(Middle_Mainstem.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Pelly_master <- filter(master_corrl,Population=="Pelly")
Pelly_master <- Pelly_master[,3:10]
Pelly.table <- cor(Pelly_master)
corrplot(Pelly.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Stewart_master <- filter(master_corrl,Population=="Stewart")
Stewart_master <- Stewart_master[,3:10]
Stewart.table <- cor(Stewart_master)
corrplot(Stewart.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Teslin_master <- filter(master_corrl,Population=="Teslin")
Teslin_master <- Teslin_master[,3:10]
Teslin.table <- cor(Teslin_master)
corrplot(Teslin.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

Upper_Lakes_master <- filter(master_corrl,Population=="Upper Lakes and Mainstem")
Upper_Lakes_master <- Upper_Lakes_master[,3:10]
Upper_Lakes.table <- cor(Upper_Lakes_master)
corrplot(Upper_Lakes.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

White_Donjek_master <- filter(master_corrl,Population=="White Donjek")
White_Donjek_master <- White_Donjek_master[,3:10]
White_Donjek.table <- cor(White_Donjek_master)
corrplot(White_Donjek.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)

#other possible variables to include:
#Sea surface temp (t = +2, +3)
#climate drivers (t = +2, +3)
#hatchery fish (t = +2, +3)


