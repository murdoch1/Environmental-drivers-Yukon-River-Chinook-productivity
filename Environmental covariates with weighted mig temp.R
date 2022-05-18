#### Processing Environmental Covariates ####

library(tidyverse)
library(lubridate)
library(mgcv)
library(imputeTS)

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
  mutate(julian_date=yday(Date))

#some years have manual and logger readings. Average these for all years except 2004 due to very different measurements
Emmonak_water_int5a <- filter(Emmonak_water_int4,Project.Year=="2004"&Recording.Method=="Logger")%>% 
  select(Project.Year,Temperature..avg.,julian_date)

Emmonak_water_int5b <- Emmonak_water_int4 %>% filter(Project.Year!="2004") %>% 
  group_by(Project.Year,julian_date) %>% summarise(new_temp=mean(Temperature..avg.)) %>% 
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
                                  
ggplot(Emmonak_water_int9,aes(y=new_temp,x=yday))+
  geom_point()+geom_smooth()+facet_wrap(~year)

#interpolate missing data for all years

Emmonak_water_int10 <- Emmonak_water_int9 %>% select(year,yday,water_temp) %>% distinct()

Emmonak_water_int11 <- Emmonak_water_int10 %>% spread(year,water_temp) %>%  
  na.interpolation()

Emmonak_water_int12 <- Emmonak_water_int11 %>% gather(year,water_temp2,2:36)
Emmonak_water_int12$year <- as.integer(Emmonak_water_int12$year)
Emmonak_water_int13 <- left_join(Emmonak_water_int9,Emmonak_water_int12)

Emmonak_water_int14 <- Emmonak_water_int13 %>% select(year,population,yday,water_temp2) %>% 
  rename(water_temp=water_temp2)

#fix 2000 due to very long missing time series at beginning of the season, days 165:187

#for estimating population means
#need to break this down and see data by yday also
Mig_temp_bypopn <- Emmonak_water_int14 %>% 
  filter(year!=2000) %>% 
  group_by(population,yday) %>% summarise(meanWT=mean(water_temp,na.rm=TRUE))

Emmonak_water_int15 <- left_join(Emmonak_water_int14,Mig_temp_bypopn) 

Emmonak_water_int16 <- Emmonak_water_int15 %>% 
  mutate(water_temp=if_else(year!=2000|(year=2000&yday>187),water_temp,meanWT)) %>% select(-meanWT)
    
ggplot(Emmonak_water_int16,aes(y=water_temp,x=yday))+
  geom_point()+geom_smooth()+facet_wrap(~year)




#temp notes by year
#1985 - good time series just impute few missing values
#1986 - same
#1987 - same
#1988 - same
#1989 - a week of missing data could try middle mouth there vs imputing - fixed
#1990 - good needs imputing
#1991 - same
#1992 - same
#1993 - same
#1994 - same
#1995 - same
#1996 - a chunk of missing data could try middle mouth there vs imputing - fixed
#1997 - a chunk of missing data could try middle mouth there vs imputing - fixed - remove outlier?
#1998 - a chunk of missing data could try middle mouth there vs imputing - left off: NEED TO ADD PILOT FOR LAST FEW
#1999 - a chunk of missing data could try middle mouth there vs imputing
#2000 - lots of data missing. use means from other years or try a mix of pilot and middle mouth data
#2001 - can use middle mouth to extend time series by two weeks. still doesn't cover upper mainstem dates
#2002 - needs imputing. maybe look at relationship with pilot for final missing chunk?
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

Eagle_temp_int1 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Eagle_water_temps.csv")) %>% 
  select(Project.Year,julian_date,water_temp) %>% rename(year=Project.Year,yday=julian_date)

Eagle_temp_int2 <- left_join(Timing_int1,Eagle_temp_int1) %>% rename(Eagle_water_temp=water_temp)

Mig_temp_int1 <- full_join(Emmonak_water_int16,Eagle_temp_int2) %>% rename(Emmonak_water_temp=water_temp)
  

#compare estimated data with actual for historical period

Eagle_temp_old_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Eagle historical biweekly water temp data.csv")) %>% 
  select(Sample.time,Temperature..water...Lab..80...VMV..2061.) %>% 
  rename(water_temp_actual=Temperature..water...Lab..80...VMV..2061.) %>% 
  mutate(year=as.numeric(format(as.Date(Sample.time, format="%Y-%m-%d", "%H:%M:%S"),"%Y")),
         yday=yday(Sample.time)) %>% select(-Sample.time)

Eagle_temp_predicted <- Eagle_temp_int1 %>% rename(Eagle_water_temp=water_temp)

Eagle_temp_compare <- left_join(Eagle_temp_old_int1,Eagle_temp_predicted)

#visualize differences
ggplot(Eagle_temp_compare,aes(y=Eagle_water_temp,x=water_temp_actual))+
  geom_point()+geom_smooth(method = "lm")

#model differences
Eagle_predictvsactual <- summary(lm(Eagle_water_temp~water_temp_actual,data=Eagle_temp_compare))
Eagle_predictvsactual


library(Metrics)
Eagle_temp_compare_NA <- Eagle_temp_compare %>% filter(!is.na(Eagle_water_temp)) %>% 
  filter(!is.na(water_temp_actual))
rmse(Eagle_temp_compare_NA$water_temp_actual,Eagle_temp_compare_NA$Eagle_water_temp)

#model overestimates temperatures
#consider modifying by time to compare morning temps only
#could try adding historical data to dataset to model to see if it improves fit ? although only afternoon temps available

#create migration timing weighting

Mig_temp_summary <- Mig_temp_int1 %>% group_by(year,population) %>% 
  summarise(count = n(),
    min = min(yday),
    max=max(yday))

Mig_30days <- as.numeric(c(seq(from = 0, to = 1, length.out = 30),rep("NA",times=4)))
#Mig_33days <- seq(from = 0, to = 1, length.out = 33)
#Mig_34days <- seq(from = 0, to = 1, length.out = 34)

Mig_30days <- bind_cols(Mig_30days,(1:34))

#Mig_33days <- bind_cols(Mig_33days,(1:33),rep(33, times=33))
#Mig_34days <- bind_cols(Mig_34days,(1:34),rep(34, times=34))

colnames(Mig_30days) <- c("weight","day")
#colnames(Mig_33days) <- c("weight","day","count")
#colnames(Mig_34days) <- c("weight","day","count")

Mig_temp_summary_count <- select(Mig_temp_summary,year,population,count)
Mig_temp_int2 <- left_join(Mig_temp_int1,Mig_temp_summary_count) 


         

Mig_30days_adj34 <- Mig_30days %>% 
  mutate(weight_adj1=lag(weight,n=1),
         weight_adj2=lag(weight,n=2),
        weight_adj3=lag(weight,n=3),
         weight_adj4=lag(weight,n=4))

Mig_30days_adj34$weighted_33 <- rowMeans(Mig_30days_adj34[,c(1,3:5)],na.rm=TRUE)
Mig_30days_adj34$weighted_34 <- rowMeans(Mig_30days_adj34[,c(1,3:6)],na.rm=TRUE)

#split into two dataframes with different number of migration timing days
Mig_temp_int2a <- Mig_temp_int2 %>% filter(count==33)
Mig_temp_int2b <- Mig_temp_int2 %>% filter(count==34)

#process 33 day migrations

Mig_temp_int2c <- Mig_temp_int2a %>% mutate(yearpopn=paste(year,population,sep="")) %>% 
                                select(yearpopn,yday,count) %>% spread(yearpopn,count)
  
Mig_temp_int2d <- select(Mig_temp_int2c,2:31)
Mig_temp_int2d_pos <- which(Mig_temp_int2d>0,arr.ind = T)
weighted33 <- Mig_30days_adj34 %>% select(weighted_33) %>% filter(!is.na(weighted_33))
Mig_temp_int2d[Mig_temp_int2d_pos] <- weighted33$weighted_33

#bind back

Mig_temp_int2e <- bind_cols(Mig_temp_int2c$yday,Mig_temp_int2d) %>% rename(yday="...1")

Mig_temp_int2f <- Mig_temp_int2e %>% gather(yearpopn,weight,2:31) %>% 
  mutate(year=substr(yearpopn,1, 4),population=substr(yearpopn,5,20))  %>% select(-yearpopn)
  
#process 34 day migrations

Mig_temp_int3a <- Mig_temp_int2b %>% mutate(yearpopn=paste(year,population,sep="")) %>% 
                                select(yearpopn,yday,count) %>% spread(yearpopn,count)
  
Mig_temp_int3b <- select(Mig_temp_int3a,2:251)
Mig_temp_int3b_pos <- which(Mig_temp_int3b>0,arr.ind = T)
weighted34 <- Mig_30days_adj34 %>% select(weighted_34) %>% filter(!is.na(weighted_34))
Mig_temp_int3b[Mig_temp_int3b_pos] <- weighted34$weighted_34

#bind back

Mig_temp_int3c <- bind_cols(Mig_temp_int3a$yday,Mig_temp_int3b) %>% rename(yday="...1")

Mig_temp_int3d <- Mig_temp_int3c %>% gather(yearpopn,weight,2:251) %>% 
  mutate(year=substr(yearpopn,1, 4),population=substr(yearpopn,5,20))  %>% select(-yearpopn)

#combine both groups
Mig_temp_int2f$weight <- as.numeric(Mig_temp_int2f$weight)
Mig_temp_int3d$weight <- as.numeric(Mig_temp_int3d$weight)
Mig_temp_int2f$year <- as.numeric(Mig_temp_int2f$year)
Mig_temp_int3d$year <- as.numeric(Mig_temp_int3d$year)

Mig_temp_int4 <- bind_rows(Mig_temp_int2f,Mig_temp_int3d)

#join weights back to original dataset and create weighted temp

Mig_temp_int5 <- left_join(Mig_temp_int1,Mig_temp_int4) %>% rename(weight_Eagle=weight) %>% 
  mutate(weight_Emmonak=1-weight_Eagle,
         migtemp_weighted=(Emmonak_water_temp*weight_Emmonak)+(Eagle_water_temp*weight_Eagle))

#try max weekly temps
Mig_temp_int5b <- left_join(Mig_temp_int5,Mig_temp_summary)
Mig_temp_int5c <- Mig_temp_int5b %>% filter(count==33) %>% arrange(year,population,yday) %>% 
  mutate(week=rep(c((rep(1:4,each=7)),(rep(5,times=5))),times=30))
Mig_temp_int5d <- left_join(Mig_temp_int5,Mig_temp_summary)
Mig_temp_int5e <- Mig_temp_int5d %>% filter(count==34) %>% arrange(year,population,yday) %>% 
  mutate(week=rep(c((rep(1:4,each=7)),(rep(5,times=6))),times=250))
  
Mig_temp_int5f  <- bind_rows(Mig_temp_int5c,Mig_temp_int5e)

Mig_temp_int5g <- Mig_temp_int5f %>% group_by(year,population,week) %>% 
  summarise(Max_weekly_temp=max(migtemp_weighted))

Mig_temp_int6 <- Mig_temp_int5g %>% filter(week<3) %>% group_by(year,population) %>%  
  summarise(water_temp=mean(Max_weekly_temp))
  
#try max temp in first two or three weeks
#try number of days over 18
#try day of year where exceeds 17 degrees

#try top two weekly max temps
#Mig_temp_int6 <- Mig_temp_int5g %>% arrange(desc(Max_weekly_temp)) %>% 
 # group_by(year,population) %>% slice(1:2) %>% summarise(water_temp=mean(Max_weekly_temp))


#data for calculating return index
#write.csv(Emmonak_water_int13,file.path(dir.data,"/Environmental data/Processed/Migration_temp_unstd.csv"),row.names=FALSE)


Mig_temp_int3 <- read.csv(file.path(dir.data,"/Environmental data/Processed/Max_weekly_migration_temp_unstd.csv"))


Mig_temp_int6 <- filter(Mig_temp_int3,year<2013)

migration_temp_plot <- Mig_temp_int6 %>% 
  mutate(population=fct_relevel(population,"LowerMainstem","White-Donjek","Stewart","Pelly",
                                "Teslin","UpperMainstem","Carmacks","MiddleMainstem")) %>% 
ggplot(aes(y=water_temp,x=population))+
  geom_boxplot()+ylab("Migration temperature (°C)")+xlab("Population")+theme_bw()+
  scale_x_discrete(labels=c("Carmacks"="Carmacks","UpperMainstem"="Upper","LowerMainstem"="Lower",
                     "Pelly"="Pelly","Stewart"="Stewart","White-Donjek"="White",
                     "Teslin"="Teslin","MiddleMainstem"="Middle"))+
  theme(text = element_text(size=25),axis.text=element_text(size=15),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y = element_text(margin=margin(t=20,r=20,b=20,l=20)))


#standardize
Mig_temp_int7 <- Mig_temp_int6
Mig_temp_int7$water_temp <- scale(Mig_temp_int7$water_temp)

Migration_temp_t0 <- Mig_temp_int7 %>%  
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

ggplot(Migration_temp_returns_int1,aes(y=water_temp,x=Population))+
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

#Bering Sea data - time series not long enough
MaySST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/MaySST.csv")) %>% select(-"ï..OID_")
M2SST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/M2SST.csv")) %>% select(-lat,-lon,-depth,-"ï..OID_")
PISST_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/PISST.csv")) %>% select(-lat,-lon,-depth,-"ï..OID_")
NEBtrawl_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Northeastern Bering Sea Trawl Temps.csv"))

SST_all <- full_join(MaySST_int1,M2SST_int1)
SST_all <- full_join(SST_all,PISST_int1)

SST_all$Year <- as.numeric(format(as.Date(SST_all$time, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))

SST_all <- SST_all %>% select(-time) %>% 
  left_join(NEBtrawl_int1)

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

SSTnew_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/BS-SST-2022-03-25.csv"))


#winter SST
SEBS_SST_int1 <- SSTnew_int1 %>% filter(Ecosystem_sub=="Southeastern Bering Sea") %>% 
  select(date,meansst)

SEBS_SST_int1$Year <- as.numeric(format(as.Date(SEBS_SST_int1$date, format="%Y-%m-%d", "%H:%M:%S"),"%Y"))
SEBS_SST_int1$Month <- as.numeric(format(as.Date(SEBS_SST_int1$date, format="%Y-%m-%d", "%H:%M:%S"),"%m"))

SEBS_SST_int2 <- SEBS_SST_int1 %>% group_by(Year) %>% 
  filter(Month=="1"|Month=="2"|Month=="3") %>% 
  summarise(Winter_SST=mean(meansst))

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

SST_winter_for_corrl <- SEBS_SST_int4
SST_summer_for_corrl <- SST_summer_int3

master_corrl <- full_join(rearing_temp_for_corrl,rearing_precip_for_corrl,by=c("Year","Population"))
master_corrl <- left_join(master_corrl,Migtemp_t0_for_corrl)
master_corrl <- left_join(master_corrl,Migtemp_returns_for_corrl)
master_corrl <- left_join(master_corrl,Ice_for_corrl)
master_corrl <- left_join(master_corrl,annual_snowpack_for_corrl)
master_corrl <- left_join(master_corrl,spawning_temp_for_corrl)
master_corrl <- left_join(master_corrl,spawning_prcp_for_corrl)
master_corrl <- left_join(master_corrl,SST_winter_for_corrl)
master_corrl <- left_join(master_corrl,SST_summer_for_corrl)
master_corrl_all <- master_corrl[,3:12]

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

