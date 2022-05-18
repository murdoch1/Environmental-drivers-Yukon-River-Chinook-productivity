#Emmonak vs Eagle data

library(tidyverse)
library(lubridate)
library(zoo)

Eagle_water_int1 <- read.csv("Data/Environmental data/Raw/Eagle station water temperature data.csv")

#hand held station only taken at 6 PM - will need to standardize to left bank station 1 and also HH-1 at 6 PM


Eagle_water_int1$Instrument.Site <- as.character(Eagle_water_int1$Instrument.Site)

Eagle_water_int2 <- Eagle_water_int1 %>% filter(Instrument.Site=="1"|Instrument.Site=="HH_1") %>% 
  select(-Time) %>% 
  rename(Date_Time=Date) %>% 
  separate(Date_Time, into = c("Date", "Time"), sep = " ", remove = FALSE) %>%
  mutate(Date = lubridate::as_date(Date, format = "%Y-%m-%d"),
         Time = hms::as_hms(str_c(Time, ":00"))) 

Eagle_water_int2$Time <- as.character(Eagle_water_int2$Time)

#Eagle_water_int3 <- Eagle_water_int2 %>% 
 # filter(Time=="18:00:00")

Eagle_water_int3 <- Eagle_water_int2 %>% 
  filter(Time=="08:00:00"|Time=="08:08:59")

#add julian date
Eagle_water_int4 <- mutate(Eagle_water_int3,julian_date=yday(Date))

ggplot(Eagle_water_int4,aes(y=Temperature..avg.,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

Eagle_water_temp <- Eagle_water_int4 %>% filter(Temperature..avg.!="NA")%>% 
  select(Project.Year,julian_date,Temperature..avg.) %>% 
  rename(Eagle_temp="Temperature..avg.")


#bring in Emmonak data and compare

Emmonak_water_int1a <- read.csv("Data/Environmental data/Raw/Emmonak Water Temp Data1.csv")
Emmonak_water_int1b <- read.csv("Data/Environmental data/Raw/Emmonak Water Temp Data2.csv")
Emmonak_water_int1a$Instrument.Site <- as.character(Emmonak_water_int1a$Instrument.Site)
Emmonak_water_int1 <- bind_rows(Emmonak_water_int1a,Emmonak_water_int1b)

Emmonak_water_int2 <- Emmonak_water_int1 %>% 
  select(-Time) %>% 
  rename(Date_Time=Date) %>% 
  separate(Date_Time, into = c("Date", "Time"), sep = " ", remove = FALSE) %>%
  mutate(Date = lubridate::as_date(Date, format = "%Y-%m-%d"),
         Time = hms::as_hms(str_c(Time, ":00")))
#note missing date and time data just for one entry in 2014 gives warning message

Emmonak_water_int2$Time <- as.character(Emmonak_water_int2$Time)

#Emmonak_water_int3 <- Emmonak_water_int2 %>%
 # filter((Time=="18:00:00"|Time=="18:15:15"|Time=="18:01:06"|Time=="17:28:43"|Time=="18:28:43"|
  #        Time=="17:41:03")&Location=="Big Eddy")

Emmonak_water_int3 <- Emmonak_water_int2 %>%
  filter((Time=="08:00:00"|Time=="08:01:00")&Location=="Big Eddy")

Emmonak_summary1 <- Emmonak_water_int3 %>% group_by(Project.Year,Time) %>% 
  summarise(count = n())

#add julian date
Emmonak_water_int4 <- mutate(Emmonak_water_int3,julian_date=yday(Date))

ggplot(Emmonak_water_int4,aes(y=Temperature..avg.,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

Emmonak_water_temp <- Emmonak_water_int4 %>% filter(Temperature..avg.!="NA") %>% 
  select(Project.Year,julian_date,Temperature..avg.) %>% 
  rename(Emmonak_temp="Temperature..avg.")

Compare <- left_join(Emmonak_water_temp,Eagle_water_temp)

Compare <- Compare %>% filter(!is.na(Eagle_temp))

ggplot(Compare,aes(y=Emmonak_temp,x=Eagle_temp))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

mod1 <- lm(Emmonak_temp~Eagle_temp,data=Compare)
summary(mod1)  

#R2 = 0.5



#try vs air temp
#Emmonak

Emmonak_air_int1 <- read.csv("Data/Environmental data/Raw/Emmonak_airtemp.csv") 

Emmonak_air <- Emmonak_air_int1 %>% filter(measurement=="tmax..deg.c."|measurement=="tmin..deg.c.") %>% 
  group_by(year,yday) %>% 
  mutate(tmean_air=mean(value)) %>% select(year,yday,tmean_air) %>% distinct() %>% 
  rename(julian_date="yday",Project.Year="year")

Emmonak_air_water <- left_join(Emmonak_water_temp,Emmonak_air)

ggplot(Emmonak_air_water,aes(y=Emmonak_temp,x=tmean_air))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)


Emmonak_air_water_long <- gather(Emmonak_air_water,type,temperature,3:4)

ggplot(Emmonak_air_water_long,aes(y=temperature,x=julian_date,color=type))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

mod1 <- lm(Emmonak_temp~tmean_air+Project.Year,data=Emmonak_air_water)
summary(mod1)

#not very close relationship - Emmonak related to upstream discharge and heat balance

#try for Eagle

Eagle_air_int1 <- read.csv("Data/Environmental data/Raw/Eagle_airtemp.csv") 

Eagle_air <- Eagle_air_int1 %>% filter(measurement=="tmax..deg.c."|measurement=="tmin..deg.c.") %>% 
  group_by(year,yday) %>% 
  mutate(tmean_air=mean(value)) %>% select(year,yday,tmean_air) %>% distinct() %>% 
  rename(julian_date="yday",Project.Year="year")

Eagle_air_water <- left_join(Eagle_water_temp,Eagle_air)

ggplot(Eagle_air_water,aes(y=Eagle_temp,x=julian_date))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)


Eagle_air_water_long <- gather(Eagle_air_water,type,temperature,3:4)

ggplot(Eagle_air_water_long,aes(y=temperature,x=julian_date,color=type))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

mod1 <- lm(Eagle_temp~tmean_air,data=Eagle_air_water)
summary(mod1)

#strong relationship


#try with lags

Eagle_air_water = Eagle_air_water %>% 
  arrange(Project.Year,julian_date) %>%
  group_by(Project.Year) %>% 
                                    mutate(airtemp.1=lag(tmean_air),
                                            airtemp.2=lag(tmean_air,n=2),
                                            airtemp.3=lag(tmean_air,n=3),
                                            airtemp.4=lag(tmean_air,n=4),
                                            airtemp.5=lag(tmean_air,n=5),
                                            airtemp.6=lag(tmean_air,n=6),
                                            airtemp.7=lag(tmean_air,n=7),
                                            airtemp.8=lag(tmean_air,n=8),
                                            airtemp.9=lag(tmean_air,n=9),
                                            airtemp.10=lag(tmean_air,n=10),
                                            airtemp.11=lag(tmean_air,n=11),
                                            airtemp.12=lag(tmean_air,n=12),
                                            airtemp.13=lag(tmean_air,n=13),
                                            airtemp.14=lag(tmean_air,n=14),
                                            airtemp.15=lag(tmean_air,n=15),
                                            airtemp.16=lag(tmean_air,n=16))


mod1 <- lm(Eagle_temp~airtemp.1,data=Eagle_air_water)
summary(mod1)

#strongest lagged 1-2 days

#try rolling means

Eagle_air_water = Eagle_air_water %>%
  group_by(Project.Year) %>%
  arrange(Project.Year,julian_date) %>%
  mutate(airtemp.lag1 = rollmean(tmean_air,1,align="right",fill=NA)) %>% 
  mutate(airtemp.lag2 = rollmean(tmean_air,2,align="right",fill=NA)) %>% 
  mutate(airtemp.lag3=rollmean(tmean_air,3,align="right",fill=NA)) %>% 
  mutate(airtemp.lag4=rollmean(tmean_air,4,align="right",fill=NA)) %>% 
  mutate(airtemp.lag5=rollmean(tmean_air,5,align="right",fill=NA)) %>% 
  mutate(airtemp.lag6=rollmean(tmean_air,6,align="right",fill=NA)) %>% 
  mutate(airtemp.lag7=rollmean(tmean_air,7,align="right",fill=NA)) %>% 
  mutate(airtemp.lag8=rollmean(tmean_air,8,align="right",fill=NA)) %>% 
  mutate(airtemp.lag9=rollmean(tmean_air,9,align="right",fill=NA)) %>% 
  mutate(airtemp.lag10=rollmean(tmean_air,10,align="right",fill=NA)) %>% 
  mutate(airtemp.lag11=rollmean(tmean_air,11,align="right",fill=NA)) %>% 
  mutate(airtemp.lag12=rollmean(tmean_air,12,align="right",fill=NA)) %>% 
  mutate(airtemp.lag13=rollmean(tmean_air,13,align="right",fill=NA)) %>% 
  mutate(airtemp.lag14=rollmean(tmean_air,14,align="right",fill=NA)) 

hist(sqrt(Eagle_air_water$airtemp.lag9))

mod1 <- lm(Eagle_temp~airtemp.lag8,data=Eagle_air_water)
summary(mod1)
plot(mod1$residuals)
shapiro.test(mod1$residuals)
qqnorm(mod1$residuals)



#r2 94% is max around 8-9 days lag including day of water temp

ggplot(Eagle_air_water,aes(y=Eagle_temp,x=airtemp.lag8))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

Eagle_air_water_long_lag <- Eagle_air_water %>% select(Project.Year,julian_date,Eagle_temp,airtemp.lag8) %>% 
                                                       gather(type,temperature,3:4)

ggplot(Eagle_air_water_long_lag,aes(y=temperature,x=julian_date,color=type))+
  geom_point()+geom_smooth()+facet_wrap(~Project.Year)

#add predicted water temperatures for Eagle

Eagle_air_water <- Eagle_air_water %>% mutate(watertemp_pred=3.559898+0.771535*airtemp.lag8)

ggplot(Eagle_air_water,aes(y=watertemp_pred,x=Eagle_temp))+
  geom_point()+geom_smooth()

mod1 <- summary(lm(watertemp_pred~Eagle_temp,data=Eagle_air_water))
mod1

Eagle_air <- Eagle_air %>% group_by(Project.Year) %>%
  arrange(Project.Year,julian_date) %>%
  mutate(airtemp.lag8 = rollmean(tmean_air,8,align="right",fill=NA)) %>% 
mutate(watertemp_pred=3.559898+0.771535*airtemp.lag8)

Eagle_air_water_all <- full_join(Eagle_air,Eagle_water_temp) %>% rename(watertemp_actual=Eagle_temp) %>% 
  filter(julian_date>150&julian_date<280) %>% 
  mutate(water_temp=if_else(is.na(watertemp_actual),watertemp_pred,watertemp_actual))

#model based on dates roughly 180-190 (June 29 - July 9) to 280 (Oct 7) but need temps back to 150 based on mig timing

write.csv(Eagle_air_water_all,"Data/Environmental Data/Processed/Eagle_water_temps.csv",row.names = FALSE)

#to do
#merge eagle data and make weighted mig temp between emmonak and eagle
#try accumulated degree days (maybe extend further upstream)
#try maximum weekly instead