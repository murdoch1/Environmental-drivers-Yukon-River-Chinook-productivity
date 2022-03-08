##### Calculating migration temperature index for returning cohorts ####

#method for deriving returning props needs qaqc#

library(tidyverse)

Age_prop_int1 <- read.csv("Data/Age_prop.csv")

#re-analyze to have proportions add to one by popn and year
#Age_prop_int2 <-Age_prop_int1 %>% 
  #group_by(year,Population) %>% 
 #summarize(sum_age=sum(age_prop)) 

#remove duplicate rows
Age_prop_int2 <- distinct(Age_prop_int1)
  
#Age_prop_int3 <- left_join(Age_prop_int1,Age_prop_int2)

#Age_prop_int4 <- Age_prop_int3 %>% mutate(age_prop=(age_prop/sum_age))%>% 
  #select(year,Fish.Age,age_prop,Population)

Age_prop_int3 <- select(Age_prop_int2,year,Fish.Age,age_prop,Population)

#Age_prop_int5 <- Age_prop_int4 %>% 
  #group_by(year,Population,Fish.Age) %>% 
  #summarize(age_prop=sum(age_prop))

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

#Age_prop_int8 <- Age_prop_int7 %>% filter(Population=="Aggregate") %>% 
  #select(year,return_est,Fish.Age) 

Age_prop_int8 <- Age_prop_int7 %>%  
  select(year,Population,return_est,Fish.Age) 

#Age_prop_int9 <- Age_prop_int8 %>% mutate(year_pop=paste(Age_prop_int8$Population,Age_prop_int8$year,sep="")) %>% 
  #select(-year,-Population) 

Age_prop_int9 <- Age_prop_int8 %>% 
  spread(key=Fish.Age,value=return_est) %>% 
  rename(age_4="4",age_5="5",age_6="6",age_7="7")

#for visualization of offset years
Age_prop_int10 <- Age_prop_int9 %>% mutate(year_4=year-4,year_5=year-5,year_6=year-6,year_7=year-7)

Age_prop_int10$year <- as.numeric(Age_prop_int10$year)


#arrange and add zeroes to missing ages
Age_prop_int11 <- Age_prop_int10 %>% arrange(Population,year) 

Age_prop_int11[is.na(Age_prop_int11)] <- 0

#dataframe with all years

year=rep(1985:2019,times=9)
Population=rep(c("Aggregate","Carmacks","Lower Mainstem","Middle Mainstem","Pelly",
             "Stewart","Teslin","Upper Lakes and Mainstem","White-Donjek"),each=35)
Age_prop_int12 <- data.frame(year,Population)
  
Age_prop_int13 <- left_join(Age_prop_int12,Age_prop_int11)

#Sub in years with population specific proportions using aggregate data proportions for that year
#missing years are 1986, 1988, 1989, 1990, 1996, 1997, 1998, 2001, 2005, 2006, 2007, 2019

#add in props of spawners by popn and year


population_props_year_int1 <- brood_table_int2 %>% 
  group_by(year) %>% 
  summarize(S_med_total=sum(S_med))

population_props_year_int2 <- left_join(brood_table_int2,population_props_year_int1)
population_props_year_int3 <- mutate(population_props_year_int2,popn_prop=S_med/S_med_total)

Age_prop_int13 <- left_join(Age_prop_int13,population_props_year_int3)

Age_prop_int14 <- Age_prop_int13 %>% 
  mutate(age_4=if_else(year=="1986"&is.na(age_4),Age_prop_int13[2,3]*popn_prop,
                    if_else(year=="1988"&is.na(age_4),Age_prop_int13[4,3]*popn_prop,
                    if_else(year=="1989"&is.na(age_4),Age_prop_int13[5,3]*popn_prop,
                    if_else(year=="1990"&is.na(age_4),Age_prop_int13[6,3]*popn_prop,
                    if_else(year=="1996"&is.na(age_4),Age_prop_int13[12,3]*popn_prop,
                    if_else(year=="1997"&is.na(age_4),Age_prop_int13[13,3]*popn_prop,
                    if_else(year=="1998"&is.na(age_4),Age_prop_int13[14,3]*popn_prop,
                    if_else(year=="2001"&is.na(age_4),Age_prop_int13[17,3]*popn_prop,
                    if_else(year=="2005"&is.na(age_4),Age_prop_int13[21,3]*popn_prop,
                    if_else(year=="2006"&is.na(age_4),Age_prop_int13[22,3]*popn_prop,
                    if_else(year=="2007"&is.na(age_4),Age_prop_int13[23,3]*popn_prop,       
                    if_else(year=="2019"&is.na(age_4),Age_prop_int13[35,3]*popn_prop,age_4)))))))))))),
                   
         age_5=if_else(year=="1986"&is.na(age_5),Age_prop_int13[2,4]*popn_prop,
                    if_else(year=="1988"&is.na(age_5),Age_prop_int13[4,4]*popn_prop,
                    if_else(year=="1989"&is.na(age_5),Age_prop_int13[5,4]*popn_prop,
                    if_else(year=="1990"&is.na(age_5),Age_prop_int13[6,4]*popn_prop,
                    if_else(year=="1996"&is.na(age_5),Age_prop_int13[12,4]*popn_prop,
                    if_else(year=="1997"&is.na(age_5),Age_prop_int13[13,4]*popn_prop,
                    if_else(year=="1998"&is.na(age_5),Age_prop_int13[14,4]*popn_prop,
                    if_else(year=="2001"&is.na(age_5),Age_prop_int13[17,4]*popn_prop,
                    if_else(year=="2005"&is.na(age_5),Age_prop_int13[21,4]*popn_prop,
                    if_else(year=="2006"&is.na(age_5),Age_prop_int13[22,4]*popn_prop,
                    if_else(year=="2007"&is.na(age_5),Age_prop_int13[23,4]*popn_prop,       
                    if_else(year=="2019"&is.na(age_5),Age_prop_int13[35,4]*popn_prop,age_5)))))))))))),
         
         age_6=if_else(year=="1986"&is.na(age_6),Age_prop_int13[2,5]*popn_prop,
                    if_else(year=="1988"&is.na(age_6),Age_prop_int13[4,5]*popn_prop,
                    if_else(year=="1989"&is.na(age_6),Age_prop_int13[5,5]*popn_prop,
                    if_else(year=="1990"&is.na(age_6),Age_prop_int13[6,5]*popn_prop,
                    if_else(year=="1996"&is.na(age_6),Age_prop_int13[12,5]*popn_prop,
                    if_else(year=="1997"&is.na(age_6),Age_prop_int13[13,5]*popn_prop,
                    if_else(year=="1998"&is.na(age_6),Age_prop_int13[14,5]*popn_prop,
                    if_else(year=="2001"&is.na(age_6),Age_prop_int13[17,5]*popn_prop,
                    if_else(year=="2005"&is.na(age_6),Age_prop_int13[21,5]*popn_prop,
                    if_else(year=="2006"&is.na(age_6),Age_prop_int13[22,5]*popn_prop,
                    if_else(year=="2007"&is.na(age_6),Age_prop_int13[23,5]*popn_prop,       
                    if_else(year=="2019"&is.na(age_6),Age_prop_int13[35,5]*popn_prop,age_6)))))))))))),
         
         age_7=if_else(year=="1986"&is.na(age_7),Age_prop_int13[2,6]*popn_prop,
                    if_else(year=="1988"&is.na(age_7),Age_prop_int13[4,6]*popn_prop,
                    if_else(year=="1989"&is.na(age_7),Age_prop_int13[5,6]*popn_prop,
                    if_else(year=="1990"&is.na(age_7),Age_prop_int13[6,6]*popn_prop,
                    if_else(year=="1996"&is.na(age_7),Age_prop_int13[12,6]*popn_prop,
                    if_else(year=="1997"&is.na(age_7),Age_prop_int13[13,6]*popn_prop,
                    if_else(year=="1998"&is.na(age_7),Age_prop_int13[14,6]*popn_prop,
                    if_else(year=="2001"&is.na(age_7),Age_prop_int13[17,6]*popn_prop,
                    if_else(year=="2005"&is.na(age_7),Age_prop_int13[21,6]*popn_prop,
                    if_else(year=="2006"&is.na(age_7),Age_prop_int13[22,6]*popn_prop,
                    if_else(year=="2007"&is.na(age_7),Age_prop_int13[23,6]*popn_prop,       
                    if_else(year=="2019"&is.na(age_7),Age_prop_int13[35,6]*popn_prop,age_7)))))))))))))


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

write.csv(MigTemp,"Data/Environmental Data/Processed/MigTemp_returns.csv",row.names = FALSE)
