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




#full list of day of year that fish are passing

summarizedDOY_long <- arrayAsDF %>%
  group_by(year, population) %>%
  summarize(yday = rep(doy, Percent))

#calculate interval needed by year using dates between 25th and 75th percentiles of the run
#adjust migration timing by one month earlier when fish pass lower Yukon where temperature data is taken

summarizedDOY_long_int1 <- left_join(summarizedDOY_long,summarizedDOY)

MigTiming_intervals_int1 <- summarizedDOY_long_int1 %>% 
  filter(yday>=q_25&yday<=q_75) %>% 
  mutate(yday_adj=yday-30) %>% 
  select(1:2,9)

#sub in year and population names
levels(MigTiming_intervals_int1$year) <- 1985:2019
levels(MigTiming_intervals_int1$population) <- (c("LowerMainstem","White-Donjek","Pelly","Stewart","Carmacks","Teslin","MiddleMainstem","UpperMainstem"))

#remove duplicates
MigTiming_intervals_int2 <- distinct(MigTiming_intervals_int1)

#add thirty day period from maximum
MigTiming_intervals_int3 <- MigTiming_intervals_int2 %>% mutate(yday_30=yday_adj+30) %>% 
  gather(stat,yday,3:4) %>% group_by(year,population) %>% summarize(yday_max=max(yday),yday_min=min(yday))

library(purrr)

MigTiming_intervals <- MigTiming_intervals_int3 %>%
  mutate(yday= map2(yday_min, yday_max, seq)) %>% unnest(cols=yday) %>% select(-yday_min,-yday_max)

write.csv(MigTiming_intervals,"Data/MigTiming_intervals.csv",row.names=FALSE)
