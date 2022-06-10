### Calculating Migration Temperature###

library(tidyverse)
library(lubridate)
library(imputeTS)

#Define Workflow Paths ====================================================
wd <- file.path(getwd())
dir.figs <- file.path(wd,"Plots")
dir.output <- file.path(wd,"Output")
dir.data <- file.path(wd,"Data")

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

#Use Instrument Site 1 for consistency
#Not indicated for year 2015 where all measurements used

Emmonak_water_int2 <- Emmonak_water_int1 %>% 
  select(-Time) %>% 
  rename(Date_Time=Date) %>% 
  separate(Date_Time, into = c("Date", "Time"), sep = " ", remove = FALSE) %>%
  mutate(Date = lubridate::as_date(Date, format = "%Y-%m-%d"),
         Time = hms::as_hms(str_c(Time, ":00"))) %>% filter(if_else(Project.Year!=2015,
                                Instrument.Site=="1",Instrument.Site=="1"|Instrument.Site!="1")) 
#note missing date and time data just for one entry in 2014 gives warning message

#sample timing differs too. before 2004 temperatures were only taken at 8 am and 8 pm
#how many have 8 am and 8 pm?

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
  filter((Time=="08:00:00"|Time=="08:01:00"|Time=="08:01:06"|Time=="20:00:00")&Location=="Big Eddy") %>% 
  select(Date,Time,Temperature..avg.) %>% 
  mutate(Time = as.factor(case_when(Time=="20:00:00"  ~ "Evening",
                             Time %in% c("08:00:00","08:01:00","08:01:06") ~ "Morning"))) 

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

#morning temps can be calculated using the linear model y = 0.978142*x + 0.043432

#2. Look at significant differences between Big Eddy and Middle Mouth morning samples

BE_MM_compare_int1 <- Emmonak_water_int2 %>% 
  filter(Time=="08:00:00"|Time=="08:01:00"|Time=="08:01:06") %>% 
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

#Big Eddy missing temps can be calculated using the linear model y = 0.78387*x + 3.25694

#3. Look at significant differences between Big Eddy (morning) and Middle Mouth (evening) samples

BE_MM_timing_compare_int1 <- Emmonak_water_int2 %>% 
  filter(((Time=="08:00:00"|Time=="08:01:00"|Time=="08:01:06")&Location=="Big Eddy")|(Time=="20:00:00"&Location=="Middle Mouth")) %>% 
  select(Location,Date,Temperature..avg.) 

BE_MM_timing_compare_int2 <- BE_MM_timing_compare_int1 %>% 
  group_by(Location,Date) %>% 
  summarize(water_temp = mean(Temperature..avg.)) 

BE_MM_timing_compare <- BE_MM_timing_compare_int2 %>% spread(Location,water_temp) %>% 
  rename(Middle_Mouth="Middle Mouth",Big_Eddy="Big Eddy")

#visualize differences
ggplot(BE_MM_timing_compare,aes(y=Middle_Mouth,x=Big_Eddy))+
  geom_point()+geom_smooth(method = "lm")

ggplot(BE_MM_timing_compare_int2,aes(y=water_temp,x=Location))+
  geom_boxplot()

#model differences
temp_time_location <- summary(lm(Big_Eddy~Middle_Mouth,data=BE_MM_timing_compare))
temp_time_location

#Big Eddy missing temps can be calculated using the linear model y = 0.73623*x + 3.69341

#4. Look at differences between Big Eddy (morning) and Pilot Station

#use instrument 2 for more consistent coverage in years needed
#broadly including times from 8 am to 10 am

Pilot_water_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Pilot Station water temp data.csv"))
Pilot_water_int1$Instrument.Site <- as.character(Pilot_water_int1$Instrument.Site)

Pilot_water_int2 <- Pilot_water_int1 %>% 
  select(-Time) %>% 
  rename(Date_Time=Date) %>% 
  separate(Date_Time, into = c("Date", "Time"), sep = " ", remove = FALSE) %>%
  mutate(Date = lubridate::as_date(Date, format = "%Y-%m-%d"),
         Time = hms::as_hms(str_c(Time, ":00"))) %>% filter(Instrument.Site=="2")

Pilot_water_int2$Time <- as.character(Pilot_water_int2$Time)
Pilot_water_int3 <- Pilot_water_int2 %>% filter(str_detect(Time,"^08|^09")) %>% 
       select(Location,Date,Temperature..avg.)

Big_eddy_morning <- Emmonak_water_int2 %>% 
  filter((Time=="08:00:00"|Time=="08:01:00"|Time=="08:01:06")&Location=="Big Eddy") %>% 
  select(Location,Date,Temperature..avg.) 

BE_PS_compare_int1 <- full_join(Pilot_water_int3,Big_eddy_morning) %>% 
group_by(Location,Date) %>% 
  summarize(water_temp = mean(Temperature..avg.)) 

BE_PS_compare <- BE_PS_compare_int1 %>% spread(Location,water_temp) %>% 
  rename(Pilot="Pilot Station Sonar",Big_Eddy="Big Eddy")

#visualize differences
ggplot(BE_PS_compare,aes(y=Big_Eddy,x=Pilot))+
  geom_point()+geom_smooth(method = "lm")

ggplot(BE_PS_compare_int1,aes(y=water_temp,x=Location))+
  geom_boxplot()

#model differences
Pilot_Emmonak <- summary(lm(Big_Eddy~Pilot,data=BE_PS_compare))
Pilot_Emmonak

#Big Eddy missing temps can be calculated using the linear model y = 0.86042*x + 1.74941



#Continue data processing
#use broader time between 8-9 Am for years 2009,2010,2014
Emmonak_water_int3 <- Emmonak_water_int2 %>% filter(if_else(Project.Year!=2009&Project.Year!=2010&Project.Year!=2014,
 (Time=="08:00:00"|Time=="08:01:06"|Time=="08:01:00"),str_detect(Time,"^08")))

Emmonak_water_summary4  <- Emmonak_water_int3 %>% 
  group_by(Project.Year,Location) %>% count(Time)

ggplot(Emmonak_water_summary4,aes(fill=Location,y=n,x=Location))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#extract Big Eddy and add julian date
Emmonak_water_int4 <- Emmonak_water_int3 %>% filter(Location=="Big Eddy") %>% 
  mutate(julian_date=yday(Date),
         week = as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%W")))

#some years have manual and logger readings. Average these for all years except 2004 due to very different measurements
Emmonak_water_int5a <- filter(Emmonak_water_int4,Project.Year=="2004"&Recording.Method=="Logger")%>% 
  select(Project.Year,week,Temperature..avg.,julian_date)

Emmonak_water_int5b <- Emmonak_water_int4 %>% filter(Project.Year!="2004") %>% 
  group_by(Project.Year,week,julian_date) %>% summarise(new_temp=mean(Temperature..avg.)) %>% 
  select(Project.Year,new_temp,julian_date) %>% rename(Temperature..avg.=new_temp)


Emmonak_water_int5c <- full_join(Emmonak_water_int5a,Emmonak_water_int5b)

#adjust migration interval by population and year using derived passage timing from population diversity paper

Timing_int1 <- read.csv(file.path(dir.data,"MigTiming_intervals.csv"))

#join to temperature data
Emmonak_water_int6 <- Emmonak_water_int5c %>% 
  rename(year=Project.Year,yday=julian_date,water_temp=Temperature..avg.)

Emmonak_water_int7 <- left_join(Timing_int1,Emmonak_water_int6)


#fill in missing data using estimated relationships above

#Eq1 = 0.978142*x + 0.043432
#Eq2 = 0.78387*x + 3.25694
#Eq3 = 0.73623*x + 3.69341
#Eq4 = 0.86042*x + 1.74941

#add other temp data for estimation

AM_PM_compare_est <- AM_PM_compare %>% mutate(yday=yday(Date)) %>% 
  mutate(year=as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))) %>% 
  select(year,yday,Evening) %>% rename(BE_evening="Evening")

AM_PM_compare_est <- as.data.frame(AM_PM_compare_est)
AM_PM_compare_est <- select(AM_PM_compare_est,-Date)

BE_MM_compare_est <- BE_MM_compare %>% mutate(yday=yday(Date)) %>% 
  mutate(year=as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))) %>% 
  select(year,yday,Middle_Mouth) %>% rename(MM_morn=Middle_Mouth) 

BE_MM_timing_compare_est <- BE_MM_timing_compare %>% mutate(yday=yday(Date)) %>% 
  mutate(year=as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))) %>% 
  select(year,yday,Middle_Mouth) %>% rename(MM_eve=Middle_Mouth) 

BE_PS_compare_est <- BE_PS_compare %>% mutate(yday=yday(Date)) %>% 
  mutate(year=as.numeric(format(as.Date(Date, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))) %>% 
  select(year,yday,Pilot)


Emmonak_water_int8a <- left_join(Emmonak_water_int7,AM_PM_compare_est)
Emmonak_water_int8b <- left_join(Emmonak_water_int8a,BE_MM_compare_est)
Emmonak_water_int8c <- left_join(Emmonak_water_int8b,BE_MM_timing_compare_est)
Emmonak_water_int8d <- left_join(Emmonak_water_int8c,BE_PS_compare_est)

Emmonak_water_int9 <- Emmonak_water_int8d %>%
  mutate(water_temp=if_else(!is.na(water_temp),water_temp,
                          if_else(!is.na(BE_evening),(0.978142*BE_evening + 0.043432),
                          if_else(!is.na(MM_morn),(0.78387*MM_morn + 3.25694),
                          if_else(!is.na(MM_eve),(0.73623*MM_eve + 3.69341),(0.86042*Pilot + 1.74941))))))
                                  
ggplot(Emmonak_water_int9,aes(y=water_temp,x=yday))+
  geom_point()+geom_smooth()+facet_wrap(~year)

#interpolate missing data for all years

Emmonak_water_int10 <- Emmonak_water_int9 %>% select(year,yday,water_temp) %>% distinct()

Emmonak_water_int11 <- Emmonak_water_int10 %>% spread(year,water_temp) %>%  
  na.interpolation()

Emmonak_water_int12 <- Emmonak_water_int11 %>% gather(year,water_temp2,2:36)
Emmonak_water_int12$year <- as.integer(Emmonak_water_int12$year)
Emmonak_water_int13 <- left_join(Emmonak_water_int9,Emmonak_water_int12)

Emmonak_water_int14 <- Emmonak_water_int13 %>% select(year,week,population,yday,water_temp2) %>% 
  rename(water_temp=water_temp2)

#fix 2000 due to very long missing time series at beginning of the season, days 165:187

#for estimating population means

Mig_temp_bypopn <- Emmonak_water_int14 %>% 
  filter(year!=2000) %>% 
  group_by(population,yday) %>% summarise(meanWT=mean(water_temp,na.rm=TRUE))

Emmonak_water_int15 <- left_join(Emmonak_water_int14,Mig_temp_bypopn) 

Emmonak_water_int16 <- Emmonak_water_int15 %>% 
  mutate(water_temp=if_else(year!=2000|(year=2000&yday>187),water_temp,meanWT)) %>% select(-meanWT)
    
ggplot(Emmonak_water_int16,aes(y=water_temp,x=yday))+
  geom_point()+geom_smooth()+facet_wrap(~year)

ggplot(Emmonak_water_int16,aes(y=water_temp,x=yday,color=population))+
  geom_point()+geom_smooth()+facet_wrap(~year)




#missing temp notes by year
#1985 - good time series just impute few missing values
#1986 - same
#1987 - same
#1988 - same
#1989 - a week of missing data
#1990 - good needs imputing
#1991 - same
#1992 - same
#1993 - same
#1994 - same
#1995 - same
#1996 - a chunk of missing data 
#1997 - a chunk of missing data 
#1998 - a chunk of missing data 
#1999 - a chunk of missing data 
#2000 - lots of data missing; will need extra step for estimation
#2001 - can use middle mouth to extend time series by two weeks. still doesn't cover upper mainstem dates
#2002 - needs imputing with Pilot for final missing chunk
#2003 - NO NA's
#2004 - NO NA's
#2005 - NO NA's
#2006 - impute missing data
#2007 - impute missing data
#2008 - NO NA's
#2009 - check this one - looks like big eddy data close to 8 am but varying numbers to add in - fixed
#2010 - same as 2009 - fixed
#2011 - NO NA's
#2012 - NO NA's
#2013 - NO NA's
#2014 - same as 2009 - fixed
#2015 - data appears to be there - check why this was screened out - fixed
#2016 - NO NA's
#2017 - data appears to be there - check why this was screened out - fixed
#2018 - NO NA's
#2019 - NO NA's







#add Eagle data and create weighted migration temperature

# Eagle_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Eagle_water_temps.csv")) %>% 
#   select(Project.Year,julian_date,water_temp) %>% rename(year=Project.Year,yday=julian_date)
# 
# Eagle_temp_int2 <- left_join(Timing_int1,Eagle_temp_int1) %>% rename(Eagle_water_temp=water_temp)
# 
# Mig_temp_int1 <- full_join(Emmonak_water_int16,Eagle_temp_int2) %>% rename(Emmonak_water_temp=water_temp)
  

#compare estimated data with actual for historical period

# Eagle_temp_old_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Eagle historical biweekly water temp data.csv")) %>% 
#   select(Sample.time,Temperature..water...Lab..80...VMV..2061.) %>% 
#   rename(water_temp_actual=Temperature..water...Lab..80...VMV..2061.) %>% 
#   mutate(year=as.numeric(format(as.Date(Sample.time, format="%Y-%m-%d", "%H:%M:%S"),"%Y")),
#          yday=yday(Sample.time)) %>% select(-Sample.time)
# 
# Eagle_temp_predicted <- Eagle_temp_int1 %>% rename(Eagle_water_temp=water_temp)
# 
# Eagle_temp_compare <- left_join(Eagle_temp_old_int1,Eagle_temp_predicted)
# 
# #visualize differences
# ggplot(Eagle_temp_compare,aes(y=Eagle_water_temp,x=water_temp_actual))+
#   geom_point()+geom_smooth(method = "lm")
# 
# #model differences
# Eagle_predictvsactual <- summary(lm(Eagle_water_temp~water_temp_actual,data=Eagle_temp_compare))
# Eagle_predictvsactual
# 
# 
# library(Metrics)
# Eagle_temp_compare_NA <- Eagle_temp_compare %>% filter(!is.na(Eagle_water_temp)) %>% 
#   filter(!is.na(water_temp_actual))
# rmse(Eagle_temp_compare_NA$water_temp_actual,Eagle_temp_compare_NA$Eagle_water_temp)

#model overestimates temperatures
#consider modifying by time to compare morning temps only
#could try adding historical data to dataset to model to see if it improves fit ? although only afternoon temps available

#create migration timing weighting

# Mig_temp_summary <- Mig_temp_int1 %>% group_by(year,population) %>% 
#   summarise(count = n(),
#     min = min(yday),
#     max=max(yday))
# 
# Mig_30days <- as.numeric(c(seq(from = 0, to = 1, length.out = 30),rep("NA",times=4)))
# #Mig_33days <- seq(from = 0, to = 1, length.out = 33)
# #Mig_34days <- seq(from = 0, to = 1, length.out = 34)
# 
# Mig_30days <- bind_cols(Mig_30days,(1:34))
# 
# #Mig_33days <- bind_cols(Mig_33days,(1:33),rep(33, times=33))
# #Mig_34days <- bind_cols(Mig_34days,(1:34),rep(34, times=34))
# 
# colnames(Mig_30days) <- c("weight","day")
# #colnames(Mig_33days) <- c("weight","day","count")
# #colnames(Mig_34days) <- c("weight","day","count")
# 
# Mig_temp_summary_count <- select(Mig_temp_summary,year,population,count)
# Mig_temp_int2 <- left_join(Mig_temp_int1,Mig_temp_summary_count) 
# 
# 
#          
# 
# Mig_30days_adj34 <- Mig_30days %>% 
#   mutate(weight_adj1=lag(weight,n=1),
#          weight_adj2=lag(weight,n=2),
#         weight_adj3=lag(weight,n=3),
#          weight_adj4=lag(weight,n=4))
# 
# Mig_30days_adj34$weighted_33 <- rowMeans(Mig_30days_adj34[,c(1,3:5)],na.rm=TRUE)
# Mig_30days_adj34$weighted_34 <- rowMeans(Mig_30days_adj34[,c(1,3:6)],na.rm=TRUE)
# 
# #split into two dataframes with different number of migration timing days
# Mig_temp_int2a <- Mig_temp_int2 %>% filter(count==33)
# Mig_temp_int2b <- Mig_temp_int2 %>% filter(count==34)
# 
# #process 33 day migrations
# 
# Mig_temp_int2c <- Mig_temp_int2a %>% mutate(yearpopn=paste(year,population,sep="")) %>% 
#                                 select(yearpopn,yday,count) %>% spread(yearpopn,count)
#   
# Mig_temp_int2d <- select(Mig_temp_int2c,2:31)
# Mig_temp_int2d_pos <- which(Mig_temp_int2d>0,arr.ind = T)
# weighted33 <- Mig_30days_adj34 %>% select(weighted_33) %>% filter(!is.na(weighted_33))
# Mig_temp_int2d[Mig_temp_int2d_pos] <- weighted33$weighted_33
# 
# #bind back
# 
# Mig_temp_int2e <- bind_cols(Mig_temp_int2c$yday,Mig_temp_int2d) %>% rename(yday="...1")
# 
# Mig_temp_int2f <- Mig_temp_int2e %>% gather(yearpopn,weight,2:31) %>% 
#   mutate(year=substr(yearpopn,1, 4),population=substr(yearpopn,5,20))  %>% select(-yearpopn)
#   
# #process 34 day migrations
# 
# Mig_temp_int3a <- Mig_temp_int2b %>% mutate(yearpopn=paste(year,population,sep="")) %>% 
#                                 select(yearpopn,yday,count) %>% spread(yearpopn,count)
#   
# Mig_temp_int3b <- select(Mig_temp_int3a,2:251)
# Mig_temp_int3b_pos <- which(Mig_temp_int3b>0,arr.ind = T)
# weighted34 <- Mig_30days_adj34 %>% select(weighted_34) %>% filter(!is.na(weighted_34))
# Mig_temp_int3b[Mig_temp_int3b_pos] <- weighted34$weighted_34
# 
# #bind back
# 
# Mig_temp_int3c <- bind_cols(Mig_temp_int3a$yday,Mig_temp_int3b) %>% rename(yday="...1")
# 
# Mig_temp_int3d <- Mig_temp_int3c %>% gather(yearpopn,weight,2:251) %>% 
#   mutate(year=substr(yearpopn,1, 4),population=substr(yearpopn,5,20))  %>% select(-yearpopn)
# 
# #combine both groups
# Mig_temp_int2f$weight <- as.numeric(Mig_temp_int2f$weight)
# Mig_temp_int3d$weight <- as.numeric(Mig_temp_int3d$weight)
# Mig_temp_int2f$year <- as.numeric(Mig_temp_int2f$year)
# Mig_temp_int3d$year <- as.numeric(Mig_temp_int3d$year)
# 
# Mig_temp_int4 <- bind_rows(Mig_temp_int2f,Mig_temp_int3d)
# 
# #join weights back to original dataset and create weighted temp
# 
# Mig_temp_int5 <- left_join(Mig_temp_int1,Mig_temp_int4) %>% rename(weight_Eagle=weight) %>% 
#   mutate(weight_Emmonak=1-weight_Eagle,
#          migtemp_weighted=(Emmonak_water_temp*weight_Emmonak)+(Eagle_water_temp*weight_Eagle))

#Mig_temp_int5b <- left_join(Mig_temp_int5,Mig_temp_summary)

#daily migration temp

Mean_Mig_temp_int1 <- Emmonak_water_int16 %>% group_by(year,population) %>% 
  summarise(water_temp=mean(water_temp))

write.csv(Mean_Mig_temp_int1,file.path(dir.data,"/Environmental data/Processed/Mean_daily_migration_temp_unstd.csv"),row.names=FALSE)


#max weekly migration temp

Emmonak_water_int16_summary <- Emmonak_water_int16 %>% 
  group_by(year,population) %>% summarise(max_temp=max(water_temp))

Mig_temp_int1 <- Emmonak_water_int16 %>% group_by(year,population,week) %>% 
  summarise(Max_weekly_temp=max(water_temp))

Mig_temp_int2 <- Mig_temp_int1  %>% group_by(year,population) %>%  
  summarise(water_temp=mean(Max_weekly_temp))

write.csv(Mig_temp_int2,file.path(dir.data,"/Environmental data/Processed/Max_weekly_migration_temp_unstd.csv"),row.names=FALSE)

#first day when river temperature exceeded a specific threshold

Threshold17_int1 <- Emmonak_water_int16 %>% arrange(year,population,yday) %>%
  group_by(year,population) %>% mutate(run_day=1:n(),count=n()) 
  
#add final day of run for those that don't cross 17

 # Threshold17_int2 <- Threshold17_int1 %>% filter(count==33) %>%
 #   mutate(new_temp=if_else(run_day<33,water_temp,17))
 # 
 # Threshold17_int3 <- Threshold17_int1 %>% filter(count==34) %>%
 #   mutate(new_temp=if_else(run_day<34,water_temp,17))
 # 
 # Threshold17_int4 <- bind_rows(Threshold17_int2,Threshold17_int3) %>% arrange(year,population,yday) %>%
 #   group_by(year,population) %>% filter(new_temp>=17) %>%
 #     slice(1) %>% rename(yday17=yday) %>% dplyr::select(year,population,yday17,run_day)

Threshold17_int4 <- Threshold17_int1 %>% arrange(year,population,yday) %>%
  group_by(year,population) %>% filter(water_temp>=17) %>% 
    slice(1) %>% rename(yday17=yday) %>% select(year,population,yday17,run_day)

ggplot(Threshold17_int4,aes(y=run_day,x=year))+
  geom_point()+geom_smooth() + facet_wrap(~population)

write.csv(Threshold17_int4,file.path(dir.data,"/Environmental data/Processed/Threshold17_unstd.csv"),row.names=FALSE)

#Number of days exceeding threshold values

DaysThreshold_int1 <- Emmonak_water_int16 %>% filter(water_temp>=17) %>% 
  group_by(year,population) %>% summarise(count=n())

DaysThreshold_int2 <- Emmonak_water_int16 %>% group_by(year,population) %>%
  summarise(total_count=n())

DaysThreshold_int3 <- left_join(DaysThreshold_int2,DaysThreshold_int1) %>%
  mutate(count=if_else(is.na(count),0,as.numeric(count))) %>% dplyr::select(-total_count)

write.csv(DaysThreshold_int3,file.path(dir.data,"/Environmental data/Processed/Threshold17_numberdays_unstd.csv"),row.names=FALSE)
