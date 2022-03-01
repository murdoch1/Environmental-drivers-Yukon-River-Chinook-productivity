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

years=rep(1985:2019,times=8)
summarizedDOY[,"year"] <- years
