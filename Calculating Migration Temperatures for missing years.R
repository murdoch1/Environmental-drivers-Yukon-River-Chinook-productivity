### Estimating migration temperatures for years with missing data ###


Mig_NAs_int1 <- read.csv("Data/Environmental data/Processed/Mig_temp_NAs.csv")

ggplot(Mig_NAs_int1,aes(y=water_temp,x=yday,color=population))+
  geom_point()+geom_smooth()+facet_wrap(~year)

#summarise migration temps for years and populations with at least 75% of dates covered

Mig_NAs_summary  <- Mig_NAs_int1 %>% 
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

#for exploring by population and year
Mig_NAs_int2 <- left_join(Mig_NAs_int1,Mig_NAs_summary_int1) %>% 
  filter(prop_NA<25) %>% 
  group_by(year,population) %>% summarise(count = n(),
                                     meanWT=mean(water_temp,na.rm=TRUE))

ggplot(Mig_NAs_int2,aes(y=meanWT,x=year))+
  geom_point()+geom_smooth()+facet_wrap(~population)

mod1 <- lm(meanWT~population+year,data=Mig_NAs_int2)
summary(mod1)

#for estimating population means
Mig_NAs <- left_join(Mig_NAs_int1,Mig_NAs_summary_int1) %>% 
  filter(prop_NA<25) %>% 
  group_by(population) %>% summarise(count = n(),
                                     meanWT=mean(water_temp,na.rm=TRUE))

write.csv(Mig_NAs,"Data/Environmental data/Processed/Mig_temp_by_popn.csv",row.names = FALSE)
