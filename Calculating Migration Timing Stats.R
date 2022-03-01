#Extract run timing by population

load("Data/MigTiming.Rdata")
library(tidyverse)

#Data structure is 1:126 (day of year), 1:8 (eight populations), 1:35 (number of years)

arrayAsDF <- as.data.frame.table(Timing)
names(arrayAsDF) <- c("reps","population","year","prop")
arrayAsDF[,"Percent"] <- arrayAsDF$prop*100

#Add day of year column
doy=rep(160:285,times=(8*35))
arrayAsDF[,"doy"] <- doy

summarizedDOY <- arrayAsDF %>%
  group_by(year, population) %>%
  summarize(reps = rep(doy, Percent)) %>%
  summarize(min= min(reps),
            q_25 = quantile(reps, 0.25),
            q_50 = quantile(reps, 0.5),
            q_75 = quantile(reps, 0.75),
            max=max(reps)) %>%
  data.frame()
summarizedDOY

#sub in year and population names
#need to verify order
popn=rep(c("LowerMainstem","White-Donjek","Pelly","Stewart","Carmacks","Teslin","MiddleMainstem","UpperMainstem"),times=(35))
summarizedDOY[,"population"] <- popn

years=rep(1985:2019,each=8)
summarizedDOY[,"year"] <- years

write.csv(summarizedDOY,"Data/MigTiming_stats.csv",row.names=FALSE)

#calculate interval needed by year

Timing_interval_25_75 <- summarizedDOY %>%
  group_by(year) %>%
  summarize(min= min(q_25),
            max=max(q_75)) %>%
  data.frame()
Timing_interval_25_75

#list of day of year that fish are passing

summarizedDOY_long <- arrayAsDF %>%
  group_by(year, population) %>%
  summarize(yday = rep(doy, Percent))

#add stats

summarizedDOY_long_int1 <- left_join(summarizedDOY_long,summarizedDOY)

MigTiming_intervals_int1 <- summarizedDOY_long_int1 %>% 
  filter(yday>=q_25&yday<=q_75) %>% 
  mutate(yday_adj=yday-30) %>% 
  select(1:2,9)

#sub in year and population names
levels(MigTiming_intervals_int1$year) <- 1985:2019
levels(MigTiming_intervals_int1$population) <- (c("LowerMainstem","White-Donjek","Pelly","Stewart","Carmacks","Teslin","MiddleMainstem","UpperMainstem"))

MigTiming_intervals <- distinct(MigTiming_intervals_int1)



write.csv(summarizedDOY,"Data/MigTiming_stats.csv",row.names=FALSE)