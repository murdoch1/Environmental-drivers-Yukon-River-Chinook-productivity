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



# Juvenile rearing temperature ------------------------------------------------

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


# Juvenile rearing precipitation ------------------------------------------------

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



# Maximum annual snow -----------------------------------------------------

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
Migration_for_corrl <- rename(Emmonak_water_int12,Year="Project.Year",migration_temp="water_temp")

Ice_for_corrl_int1 <- data.frame(rename(Ice_int2,Breakup_day="Julian_day"))
Ice_for_corrl_int2 <- select(Ice_for_corrl_int1,-Date)
Ice_for_corrl_int2$Year <- as.numeric(Ice_for_corrl_int2$Year)
Ice_for_corrl <- Ice_for_corrl_int2

rearing_precip_for_corrl <- Prcp_int8
rearing_temp_for_corrl <-Temp_int8
annual_snowpack_for_corrl <- Swe_int6

master_corrl <- full_join(rearing_temp_for_corrl,rearing_precip_for_corrl,by=c("Year","Population"))
master_corrl <- left_join(master_corrl,Migration_for_corrl)
master_corrl <- left_join(master_corrl,Ice_for_corrl)
master_corrl <- left_join(master_corrl,annual_snowpack_for_corrl)
master_corrl <- master_corrl[,3:7]

cor.table <- cor(master_corrl)
write.csv(cor.table,file.path(dir.data,"/Environmental data/Processed/cor_table.csv"),row.names=FALSE)

library(corrplot)
corrplot(cor.table, method="number",type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, insig="blank",sig.level=0.01)


#To do; check for missing data within each year and determine if needs estimation
#To do: check for obvious outliers/ measurement errors
#TO DO: derive weekly maximums?
#TO DO: estimate run timing past Emmonak and adjust migration temp by population and maybe even by year
#To do: confirm that Emmonak temps are outside of delta influence
#To do: figure out if there is a way to use migration temps for returning years also - this will matter if returns are too stressed to hit the border count

#other possible variables to include:
#Mean precipitation during spawning and early incubation (Aug - Nov) in basin (t = 0)
#Temperature for spawning/incubation in basin (t = 0)
#Sea surface temp (t = +2, +3)
#climate drivers (t = +2, +3)
#hatchery fish (t = +2, +3)


