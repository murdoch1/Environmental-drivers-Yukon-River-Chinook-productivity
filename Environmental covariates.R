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

Emmonak_water_int1 <- read.csv(file.path(dir.data,"/Environmental data/Raw/Emmonak Water Temp Data.csv"))

#sampling locations across year change, look at summary

Emmonak_water_summary1  <- Emmonak_water_int1 %>% 
  group_by(Project.Year,Location) %>% count(Location)

#two years missing middle mouth and one missing big eddy

#ggplot(Emmonak_water_int1_summary,aes(fill=Location,y=n,x=Project.Year))+
#geom_bar(stat="identity",position="dodge")+scale_y_continuous(trans="log2")

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

Emmonak_water_summary2  <- Emmonak_water_int2 %>% 
  group_by(Project.Year,Location) %>% count(Time)


Emmonak_water_summary3 <- filter(Emmonak_water_summary2,Project.Year<2004)

ggplot(Emmonak_water_summary3,aes(fill=Location,y=n,x=Time))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#missing a few PM samples and only one AM sample in 1995. 
#Use AM for now and sub in PM for 1995 until this can be adjusted

Emmonak_water_int2$Time <- as.character(Emmonak_water_int2$Time)
Emmonak_water_int3 <- filter(Emmonak_water_int2,Time=="08:00:00")

ggplot(Emmonak_water_int3,aes(x=Project.Year,y=Temperature..avg.))+
  geom_point()

Emmonak_water_summary4  <- Emmonak_water_int3 %>% 
  group_by(Project.Year,Location) %>% count(Time)

ggplot(Emmonak_water_summary4,aes(fill=Location,y=n,x=Location))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Project.Year)

#best bet for now is to use Big Eddy data at 8 AM. Will be missing for three years:
#1990 was only middle mouth, 1995 not sampled at 8 AM, 2003 also only sampled at MM
#For now sub in MM for BE in 1990 and 2003 and use 8 PM for 1995 from Big Eddy

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

#figure out general timing window for passing Emmonak for Can-origin 
#travel time is approx one month between river mouth and the border
#estimate late May to early Aug but pick window based on comparability of data to start

#looks like data missing before July 15 in 2000
#1999 started a bit late - June 15
#can maybe pull data from Pilot for this year if needed?
#check Jones paper - does the model just infer absences? or is better to estimate inputs

Emmonak_water_int9 <- filter(Emmonak_water_int8,Temperature..avg.!="NA")

Emmonak_water_summary5  <- Emmonak_water_int9 %>% 
  group_by(Project.Year) %>% summarise(
    count = n(),
    min = min(julian_date),
    max=max(julian_date))

#most years have data from June 9 - Aug 1. use this window for now

Emmonak_water_int10 <- filter(Emmonak_water_int9,julian_date>159&julian_date<214)

ggplot(Emmonak_water_int10,aes(y=Temperature..avg.,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

Emmonak_water_summary6  <- Emmonak_water_int10 %>% 
  group_by(Project.Year) %>% summarise(
    count = n(),
    min = min(julian_date),
    max=max(julian_date))


#calculate mean daily water temp during migration as temporary covariate

Emmonak_water_int11 <- Emmonak_water_int10 %>% 
  group_by(Project.Year) %>% summarise(water_temp=mean(Temperature..avg.))

#remove extra years i.e. <1985 and >2012
Emmonak_water_int12 <- filter(Emmonak_water_int11,Project.Year>1984&Project.Year<2013)

ggplot(Emmonak_water_int12,aes(fill=water_temp,y=water_temp,x=Project.Year))+
  geom_bar(stat="identity")

Migration_temp_int1 <- select(Emmonak_water_int12,water_temp)

#standardize
Migration_temp_int2 <- t(scale(Migration_temp_int1$water_temp))

#copy rows down
Migration_temp <- as.matrix(Migration_temp_int2[rep(seq_len(nrow(Migration_temp_int2)),each=8),])

write.csv(Migration_temp,file.path(dir.data,"/Environmental data/Processed/Migration_temp.csv"),row.names=FALSE)

#To do; check for missing data within each year and determine if needs estimation
#To do: check for obvious outliers/ measurement errors
#TO DO: derive weekly maximums?
#TO DO: estimate run timing past Emmonak and adjust migration temp by population and maybe even by year
#To do: confirm that Emmonak temps are outside of delta influence
#To do: figure out if there is a way to use migration temps for returning years also - this will matter if returns are too stressed to hit the border count

#other possible variables to include:
#Mean summer air temperature (June - Aug) in basin for rearing juveniles (t = +1)
#Mean summer precipitation (June - Aug) in basin for rearing juveniles (t = +1)
#Mean precipitation during spawning and early incubation (Aug - Nov) in basin (t = 0)
#Temperature for spawning/incubation in basin (t = 0)
#Sea surface temp (t = +2, +3)
#climate drivers (t = +2, +3)
#hatchery fish (t = +2, +3)


