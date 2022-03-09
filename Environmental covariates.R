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
#note missing date and time data just for one entry in 2014

Emmonak_water_summary2  <- Emmonak_water_int2 %>% 
  group_by(Project.Year,Location) %>% count(Time)

#examining data before 2004
Emmonak_water_summary3 <- filter(Emmonak_water_summary2,Project.Year<2004)

ggplot(Emmonak_water_summary3,aes(fill=Location,y=n,x=Time))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#missing a few PM samples and only one AM sample in 1995. 
#Use AM for now and sub in PM for 1995 until this can be adjusted

Emmonak_water_int2$Time <- as.character(Emmonak_water_int2$Time)
Emmonak_water_int3 <- filter(Emmonak_water_int2,Time=="08:00:00"|Time=="08:01:00")

ggplot(Emmonak_water_int3,aes(x=Project.Year,y=Temperature..avg.))+
  geom_point()

Emmonak_water_summary4  <- Emmonak_water_int3 %>% 
  group_by(Project.Year,Location) %>% count(Time)

ggplot(Emmonak_water_summary4,aes(fill=Location,y=n,x=Location))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#best bet for now is to use Big Eddy data at 8 AM. Will be missing for three years:
#1990 was only middle mouth, 1995 not sampled at 8 AM, 2003 also only sampled at MM
#For now sub in MM for BE in 1990 and 2003 and use 8 PM for 1995 from Big Eddy** fix or remove later

#extract Big Eddy
Emmonak_water_int4 <- filter(Emmonak_water_int3,Location=="Big Eddy")

#extract MM for 1990 and 2003
Emmonak_water_int5 <- filter(Emmonak_water_int3,Location=="Middle Mouth"&Project.Year=="1990"|
                               Location=="Middle Mouth"&Project.Year=="2003")

#extract 8 PM for Big Eddy in 1995
Emmonak_water_int6 <- filter(Emmonak_water_int2,Time=="20:00:00"&Location=="Big Eddy"&Project.Year=="1995")

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

Timing_int1 <- rename(Timing_int1,yday=yday_adj)

#join to temperature data
#note some years have missing data (see Temp_data_availability for comparison in main folder)
#estimate later ?

Emmonak_water_int10 <- Emmonak_water_int9 %>% 
  select(Project.Year,Temperature..avg.,julian_date) %>% 
  rename(year=Project.Year,yday=julian_date,water_temp=Temperature..avg.)

Emmonak_water_int11 <- left_join(Timing_int1,Emmonak_water_int10)

ggplot(Emmonak_water_int11,aes(y=water_temp,x=yday,))+
  geom_point()+geom_smooth()+facet_wrap(~year)

Emmonak_water_summary7  <- Emmonak_water_int11 %>% 
  group_by(year,population) %>% summarise(
    count = n(),
    min = min(yday),
    max=max(yday),
    minWT=min(water_temp),
    meanWT=mean(water_temp),
    maxZWT=max(water_temp))


#calculate mean daily water temp during migration

Emmonak_water_int12 <- Emmonak_water_int11 %>% 
  group_by(year,population) %>% summarise(water_temp=mean(water_temp,na.rm=TRUE))

#remove extra years
Emmonak_water_int13 <- filter(Emmonak_water_int12,year>1984&year<2020)

Emmonak_water_summary8  <- Emmonak_water_int13 %>% 
  group_by(population) %>% summarise(
    count = n(),
    meanWT=mean(water_temp,na.rm=TRUE))

#add missing data
Emmonak_water_int13[32,3] <- 12 #WD 1988
Emmonak_water_int13[38,3] <- 14 #Teslin 1989
Emmonak_water_int13[84,3] <- 14.5 #Pelly 1995
Emmonak_water_int13[90,3] <- 8.2 #LM 1996
Emmonak_water_int13[109,3] <- 14.7 #Stewart 1998
Emmonak_water_int13[121,3] <- 15.1 #Carmacks 2000
Emmonak_water_int13[122,3] <- 12.7 #LM 2000
Emmonak_water_int13[123,3] <- 15.6 #MM 2000
Emmonak_water_int13[124,3] <- 14.2 #Pelly 2000
Emmonak_water_int13[125,3] <- 13.9 #Stewart 2000
Emmonak_water_int13[126,3] <- 14.7 #Teslin 2000
Emmonak_water_int13[127,3] <- 16.4 #UM 2000
Emmonak_water_int13[128,3] <- 13.8 #WJ 2000
Emmonak_water_int13[135,3] <- 14.4 #UM 2001
Emmonak_water_int13[173,3] <- 11.5 #Stewart 2006

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

Migration_temp_returns_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/MigTemp_returns.csv"))
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

#Using Daymet monthly air precipitation data processed in ArcGIS Pro to watershed-level means

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

master_corrl <- full_join(rearing_temp_for_corrl,rearing_precip_for_corrl,by=c("Year","Population"))
master_corrl <- left_join(master_corrl,Migtemp_t0_for_corrl)
master_corrl <- left_join(master_corrl,Migtemp_returns_for_corrl)
master_corrl <- left_join(master_corrl,Ice_for_corrl)
master_corrl <- left_join(master_corrl,annual_snowpack_for_corrl)
master_corrl <- master_corrl[,3:8]

cor.table <- cor(master_corrl)

library(corrplot)
corrplot(cor.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)


#To do; check for missing data within each year and determine if needs estimation
#TO DO: derive weekly maximums?

#other possible variables to include:
#Mean precipitation during spawning and early incubation (Aug - Nov) in basin (t = 0)
#Temperature for spawning/incubation in basin (t = 0)
#Sea surface temp (t = +2, +3)
#climate drivers (t = +2, +3)
#hatchery fish (t = +2, +3)


