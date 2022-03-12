##### Calculating migration temperature index for returning cohorts ####

#method for deriving returning props needs qaqc#

library(tidyverse)

Age_prop_int1 <- read.csv("Data/Age_prop.csv")

#remove duplicate rows
Age_prop_int2 <- distinct(Age_prop_int1)
  
Age_prop_int3 <- select(Age_prop_int2,year,Fish.Age,age_prop,Population,count)

#use median spawners from brood table to estimate fish escapement by age, year, popn

brood_table_int1 <- read.csv("Data/brood_table.csv")

brood_table_int2 <- brood_table_int1 %>% select(BroodYear,population,S_med) %>% 
  rename(year=BroodYear,Population=population) 

brood_table_int3 <- brood_table_int1 %>% 
  group_by(BroodYear) %>% summarize(S_med=sum(S_med)) %>% rename(year=BroodYear) %>% 
  mutate(Population="Aggregate")

brood_table_int4 <- bind_rows(brood_table_int2,brood_table_int3)

Age_prop_int6 <- left_join(Age_prop_int3,brood_table_int4)

Age_prop_int7 <- mutate(Age_prop_int6,return_est=age_prop*S_med*1000)


#calculate proportion of returns by brood year

Age_prop_int8 <- Age_prop_int7 %>%  
  select(year,Population,return_est,Fish.Age,count) 

Age_prop_int9 <- Age_prop_int8 %>% 
  spread(key=Fish.Age,value=return_est) %>% 
  rename(age_4="4",age_5="5",age_6="6",age_7="7")

Age_prop_int9$year <- as.numeric(Age_prop_int9$year)


#arrange and add zeroes to missing ages
Age_prop_int10 <- Age_prop_int9 %>% arrange(Population,year) 

Age_prop_int10[is.na(Age_prop_int10)] <- 0

#dataframe with all years

year=rep(1985:2019,times=9)
Population=rep(c("Aggregate","Carmacks","Lower Mainstem","Middle Mainstem","Pelly",
             "Stewart","Teslin","Upper Lakes and Mainstem","White-Donjek"),each=35)
Age_prop_int11 <- data.frame(year,Population)
  
Age_prop_int12 <- left_join(Age_prop_int11,Age_prop_int10)

#Sub in years with population specific proportions using aggregate data proportions for that year
#missing years are 1986, 1988, 1989, 1990, 1996, 1997, 1998, 2001, 2005, 2006, 2007, 2019

#add in props of spawners by popn and year
population_props_year_int1 <- brood_table_int2 %>% 
  group_by(year) %>% 
  summarize(S_med_total=sum(S_med))

population_props_year_int2 <- left_join(brood_table_int2,population_props_year_int1)
population_props_year_int3 <- mutate(population_props_year_int2,popn_prop=S_med/S_med_total)

Age_prop_int13 <- left_join(Age_prop_int12,population_props_year_int3)

#replacing data with less than ten samples with aggregate proportions as well

Age_prop_int13b <- Age_prop_int13 %>% filter(Population=="Aggregate") %>% 
  rename(age_4_aggr=age_4,age_5_aggr=age_5,age_6_aggr=age_6,age_7_aggr=age_7) %>% 
  select(year,age_4_aggr,age_5_aggr,age_6_aggr,age_7_aggr)

Age_prop_int13c <- left_join(Age_prop_int13,Age_prop_int13b)

Age_prop_int14 <- Age_prop_int13c %>% 
  mutate(age_4=if_else(is.na(count)|count<10&Population!="Aggregate",age_4_aggr*popn_prop,age_4),
         age_5=if_else(is.na(count)|count<10&Population!="Aggregate",age_5_aggr*popn_prop,age_5),
         age_6=if_else(is.na(count)|count<10&Population!="Aggregate",age_6_aggr*popn_prop,age_6),
         age_7=if_else(is.na(count)|count<10&Population!="Aggregate",age_7_aggr*popn_prop,age_7))


#calculate age proportions based on returns by broodyear

n_offset_4=rep(4:318)
n_offset_5=rep(5:319)
n_offset_6=rep(6:320)
n_offset_7=rep(7:321)

Age_prop_int15 <- Age_prop_int14 %>% 
  mutate(age_4_adj=age_4[row_number(1)+n_offset_4],
         age_5_adj=age_5[row_number(1)+n_offset_5],
         age_6_adj=age_6[row_number(1)+n_offset_6],
         age_7_adj=age_7[row_number(1)+n_offset_7])
#this creates some junk years when rows above start drawing from the wrong popn which are removed below

Age_prop_int16 <- Age_prop_int15 %>% 
  mutate(total_returns=age_4_adj+age_5_adj+age_6_adj+age_7_adj) %>% 
  mutate(prop_4=age_4_adj/total_returns,prop_5=age_5_adj/total_returns,
         prop_6=age_6_adj/total_returns,prop_7=age_7_adj/total_returns) %>% 
  filter(year<2013)

Age_prop <- select(Age_prop_int16,year,Population,prop_4,prop_5,prop_6,prop_7)


#load unstandardized water temperature data from env covariate sheet

MigTemp_int1 <- read.csv("Data/Environmental data/Processed/Migration_temp_unstd.csv")

MigTemp_int2 <- rename(MigTemp_int1,Population=population)
MigTemp_int2$Population <- as.factor(MigTemp_int2$Population)

levels(MigTemp_int2$Population)<- list("Lower Mainstem"="LowerMainstem","White-Donjek"="White-Donjek","Middle Mainstem"="MiddleMainstem","Upper Lakes and Mainstem"="UpperMainstem",
   Carmacks="Carmacks",Teslin="Teslin",Stewart="Stewart",Pelly="Pelly")


MigTemp_int3 <- left_join(MigTemp_int2,Age_prop)

MigTemp_int3 <- arrange(MigTemp_int3,Population,year)

n_offset_4_MT=rep(4:283)
n_offset_5_MT=rep(5:284)
n_offset_6_MT=rep(6:285)
n_offset_7_MT=rep(7:286)

MigTemp <- MigTemp_int3 %>% 
  mutate(Mig_temp_returns=prop_4*water_temp[row_number(1)+n_offset_4_MT]+
  prop_5*water_temp[row_number(1)+n_offset_5_MT]+
  prop_6*water_temp[row_number(1)+n_offset_6_MT]+
  prop_7*water_temp[row_number(1)+n_offset_7_MT]) %>% 
  filter(year<2013)

write.csv(MigTemp,"Data/Environmental Data/Processed/MigTemp_returns_unformatted.csv",row.names = FALSE)
