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

left <- 1 #one as left bound because aggregate samples have zero in the count column and these need to be treated differently 
right <- 9  

Age_prop_int14 <- Age_prop_int13 %>% 
  mutate(age_4=if_else(year=="1985"&(is.na(count)|between(count,left,right)),Age_prop_int13[1,4]*popn_prop,
                    if_else(year=="1986"&(is.na(count)|between(count,left,right)),Age_prop_int13[2,4]*popn_prop,
                    if_else(year=="1987"&(is.na(count)|between(count,left,right)),Age_prop_int13[3,4]*popn_prop,
                    if_else(year=="1988"&(is.na(count)|between(count,left,right)),Age_prop_int13[4,4]*popn_prop,
                    if_else(year=="1989"&(is.na(count)|between(count,left,right)),Age_prop_int13[5,4]*popn_prop,
                    if_else(year=="1990"&(is.na(count)|between(count,left,right)),Age_prop_int13[6,4]*popn_prop,
                    if_else(year=="1991"&(is.na(count)|between(count,left,right)),Age_prop_int13[7,4]*popn_prop,
                    if_else(year=="1992"&(is.na(count)|between(count,left,right)),Age_prop_int13[8,4]*popn_prop,
                    if_else(year=="1993"&(is.na(count)|between(count,left,right)),Age_prop_int13[9,4]*popn_prop,
                    if_else(year=="1994"&(is.na(count)|between(count,left,right)),Age_prop_int13[10,4]*popn_prop,  
                    if_else(year=="1995"&(is.na(count)|between(count,left,right)),Age_prop_int13[11,4]*popn_prop,
                    if_else(year=="1996"&(is.na(count)|between(count,left,right)),Age_prop_int13[12,4]*popn_prop,
                    if_else(year=="1997"&(is.na(count)|between(count,left,right)),Age_prop_int13[13,4]*popn_prop,
                    if_else(year=="1998"&(is.na(count)|between(count,left,right)),Age_prop_int13[14,4]*popn_prop,
                    if_else(year=="1999"&(is.na(count)|between(count,left,right)),Age_prop_int13[15,4]*popn_prop,
                    if_else(year=="2000"&(is.na(count)|between(count,left,right)),Age_prop_int13[16,4]*popn_prop,
                    if_else(year=="2001"&(is.na(count)|between(count,left,right)),Age_prop_int13[17,4]*popn_prop,
                    if_else(year=="2002"&(is.na(count)|between(count,left,right)),Age_prop_int13[18,4]*popn_prop,
                    if_else(year=="2003"&(is.na(count)|between(count,left,right)),Age_prop_int13[19,4]*popn_prop,
                    if_else(year=="2004"&(is.na(count)|between(count,left,right)),Age_prop_int13[20,4]*popn_prop,
                    if_else(year=="2005"&(is.na(count)|between(count,left,right)),Age_prop_int13[21,4]*popn_prop,
                    if_else(year=="2006"&(is.na(count)|between(count,left,right)),Age_prop_int13[22,4]*popn_prop,
                    if_else(year=="2007"&(is.na(count)|between(count,left,right)),Age_prop_int13[23,4]*popn_prop,
                    if_else(year=="2008"&(is.na(count)|between(count,left,right)),Age_prop_int13[24,4]*popn_prop,
                    if_else(year=="2009"&(is.na(count)|between(count,left,right)),Age_prop_int13[25,4]*popn_prop,
                    if_else(year=="2010"&(is.na(count)|between(count,left,right)),Age_prop_int13[26,4]*popn_prop,
                    if_else(year=="2011"&(is.na(count)|between(count,left,right)),Age_prop_int13[27,4]*popn_prop,
                    if_else(year=="2012"&(is.na(count)|between(count,left,right)),Age_prop_int13[28,4]*popn_prop,
                    if_else(year=="2013"&(is.na(count)|between(count,left,right)),Age_prop_int13[29,4]*popn_prop,
                    if_else(year=="2014"&(is.na(count)|between(count,left,right)),Age_prop_int13[30,4]*popn_prop,
                    if_else(year=="2015"&(is.na(count)|between(count,left,right)),Age_prop_int13[31,4]*popn_prop,
                    if_else(year=="2016"&(is.na(count)|between(count,left,right)),Age_prop_int13[32,4]*popn_prop,
                    if_else(year=="2017"&(is.na(count)|between(count,left,right)),Age_prop_int13[33,4]*popn_prop,
                    if_else(year=="2018"&(is.na(count)|between(count,left,right)),Age_prop_int13[34,4]*popn_prop,
                    if_else(year=="2019"&(is.na(count)|between(count,left,right)),Age_prop_int13[35,4]*popn_prop,age_4))))))))))))))))))))))))))))))))))),
                   
         age_5=if_else(year=="1985"&(is.na(count)|between(count,left,right)),Age_prop_int13[1,5]*popn_prop,
                    if_else(year=="1986"&(is.na(count)|between(count,left,right)),Age_prop_int13[2,5]*popn_prop,
                    if_else(year=="1987"&(is.na(count)|between(count,left,right)),Age_prop_int13[3,5]*popn_prop,
                    if_else(year=="1988"&(is.na(count)|between(count,left,right)),Age_prop_int13[4,5]*popn_prop,
                    if_else(year=="1989"&(is.na(count)|between(count,left,right)),Age_prop_int13[5,5]*popn_prop,
                    if_else(year=="1990"&(is.na(count)|between(count,left,right)),Age_prop_int13[6,5]*popn_prop,
                    if_else(year=="1991"&(is.na(count)|between(count,left,right)),Age_prop_int13[7,5]*popn_prop,
                    if_else(year=="1992"&(is.na(count)|between(count,left,right)),Age_prop_int13[8,5]*popn_prop,
                    if_else(year=="1993"&(is.na(count)|between(count,left,right)),Age_prop_int13[9,5]*popn_prop,
                    if_else(year=="1994"&(is.na(count)|between(count,left,right)),Age_prop_int13[10,5]*popn_prop,  
                    if_else(year=="1995"&(is.na(count)|between(count,left,right)),Age_prop_int13[11,5]*popn_prop,
                    if_else(year=="1996"&(is.na(count)|between(count,left,right)),Age_prop_int13[12,5]*popn_prop,
                    if_else(year=="1997"&(is.na(count)|between(count,left,right)),Age_prop_int13[13,5]*popn_prop,
                    if_else(year=="1998"&(is.na(count)|between(count,left,right)),Age_prop_int13[14,5]*popn_prop,
                    if_else(year=="1999"&(is.na(count)|between(count,left,right)),Age_prop_int13[15,5]*popn_prop,
                    if_else(year=="2000"&(is.na(count)|between(count,left,right)),Age_prop_int13[16,5]*popn_prop,
                    if_else(year=="2001"&(is.na(count)|between(count,left,right)),Age_prop_int13[17,5]*popn_prop,
                    if_else(year=="2002"&(is.na(count)|between(count,left,right)),Age_prop_int13[18,5]*popn_prop,
                    if_else(year=="2003"&(is.na(count)|between(count,left,right)),Age_prop_int13[19,5]*popn_prop,
                    if_else(year=="2004"&(is.na(count)|between(count,left,right)),Age_prop_int13[20,5]*popn_prop,
                    if_else(year=="2005"&(is.na(count)|between(count,left,right)),Age_prop_int13[21,5]*popn_prop,
                    if_else(year=="2006"&(is.na(count)|between(count,left,right)),Age_prop_int13[22,5]*popn_prop,
                    if_else(year=="2007"&(is.na(count)|between(count,left,right)),Age_prop_int13[23,5]*popn_prop,
                    if_else(year=="2008"&(is.na(count)|between(count,left,right)),Age_prop_int13[24,5]*popn_prop,
                    if_else(year=="2009"&(is.na(count)|between(count,left,right)),Age_prop_int13[25,5]*popn_prop,
                    if_else(year=="2010"&(is.na(count)|between(count,left,right)),Age_prop_int13[26,5]*popn_prop,
                    if_else(year=="2011"&(is.na(count)|between(count,left,right)),Age_prop_int13[27,5]*popn_prop,
                    if_else(year=="2012"&(is.na(count)|between(count,left,right)),Age_prop_int13[28,5]*popn_prop,
                    if_else(year=="2013"&(is.na(count)|between(count,left,right)),Age_prop_int13[29,5]*popn_prop,
                    if_else(year=="2014"&(is.na(count)|between(count,left,right)),Age_prop_int13[30,5]*popn_prop,
                    if_else(year=="2015"&(is.na(count)|between(count,left,right)),Age_prop_int13[31,5]*popn_prop,
                    if_else(year=="2016"&(is.na(count)|between(count,left,right)),Age_prop_int13[32,5]*popn_prop,
                    if_else(year=="2017"&(is.na(count)|between(count,left,right)),Age_prop_int13[33,5]*popn_prop,
                    if_else(year=="2018"&(is.na(count)|between(count,left,right)),Age_prop_int13[34,5]*popn_prop,
                    if_else(year=="2019"&(is.na(count)|between(count,left,right)),Age_prop_int13[35,5]*popn_prop,age_5))))))))))))))))))))))))))))))))))),
        
           age_6=if_else(year=="1985"&(is.na(count)|between(count,left,right)),Age_prop_int13[1,6]*popn_prop,
                    if_else(year=="1986"&(is.na(count)|between(count,left,right)),Age_prop_int13[2,6]*popn_prop,
                    if_else(year=="1987"&(is.na(count)|between(count,left,right)),Age_prop_int13[3,6]*popn_prop,
                    if_else(year=="1988"&(is.na(count)|between(count,left,right)),Age_prop_int13[4,6]*popn_prop,
                    if_else(year=="1989"&(is.na(count)|between(count,left,right)),Age_prop_int13[5,6]*popn_prop,
                    if_else(year=="1990"&(is.na(count)|between(count,left,right)),Age_prop_int13[6,6]*popn_prop,
                    if_else(year=="1991"&(is.na(count)|between(count,left,right)),Age_prop_int13[7,6]*popn_prop,
                    if_else(year=="1992"&(is.na(count)|between(count,left,right)),Age_prop_int13[8,6]*popn_prop,
                    if_else(year=="1993"&(is.na(count)|between(count,left,right)),Age_prop_int13[9,6]*popn_prop,
                    if_else(year=="1994"&(is.na(count)|between(count,left,right)),Age_prop_int13[10,6]*popn_prop,  
                    if_else(year=="1995"&(is.na(count)|between(count,left,right)),Age_prop_int13[11,6]*popn_prop,
                    if_else(year=="1996"&(is.na(count)|between(count,left,right)),Age_prop_int13[12,6]*popn_prop,
                    if_else(year=="1997"&(is.na(count)|between(count,left,right)),Age_prop_int13[13,6]*popn_prop,
                    if_else(year=="1998"&(is.na(count)|between(count,left,right)),Age_prop_int13[14,6]*popn_prop,
                    if_else(year=="1999"&(is.na(count)|between(count,left,right)),Age_prop_int13[15,6]*popn_prop,
                    if_else(year=="2000"&(is.na(count)|between(count,left,right)),Age_prop_int13[16,6]*popn_prop,
                    if_else(year=="2001"&(is.na(count)|between(count,left,right)),Age_prop_int13[17,6]*popn_prop,
                    if_else(year=="2002"&(is.na(count)|between(count,left,right)),Age_prop_int13[18,6]*popn_prop,
                    if_else(year=="2003"&(is.na(count)|between(count,left,right)),Age_prop_int13[19,6]*popn_prop,
                    if_else(year=="2004"&(is.na(count)|between(count,left,right)),Age_prop_int13[20,6]*popn_prop,
                    if_else(year=="2005"&(is.na(count)|between(count,left,right)),Age_prop_int13[21,6]*popn_prop,
                    if_else(year=="2006"&(is.na(count)|between(count,left,right)),Age_prop_int13[22,6]*popn_prop,
                    if_else(year=="2007"&(is.na(count)|between(count,left,right)),Age_prop_int13[23,6]*popn_prop,
                    if_else(year=="2008"&(is.na(count)|between(count,left,right)),Age_prop_int13[24,6]*popn_prop,
                    if_else(year=="2009"&(is.na(count)|between(count,left,right)),Age_prop_int13[25,6]*popn_prop,
                    if_else(year=="2010"&(is.na(count)|between(count,left,right)),Age_prop_int13[26,6]*popn_prop,
                    if_else(year=="2011"&(is.na(count)|between(count,left,right)),Age_prop_int13[27,6]*popn_prop,
                    if_else(year=="2012"&(is.na(count)|between(count,left,right)),Age_prop_int13[28,6]*popn_prop,
                    if_else(year=="2013"&(is.na(count)|between(count,left,right)),Age_prop_int13[29,6]*popn_prop,
                    if_else(year=="2014"&(is.na(count)|between(count,left,right)),Age_prop_int13[30,6]*popn_prop,
                    if_else(year=="2015"&(is.na(count)|between(count,left,right)),Age_prop_int13[31,6]*popn_prop,
                    if_else(year=="2016"&(is.na(count)|between(count,left,right)),Age_prop_int13[32,6]*popn_prop,
                    if_else(year=="2017"&(is.na(count)|between(count,left,right)),Age_prop_int13[33,6]*popn_prop,
                    if_else(year=="2018"&(is.na(count)|between(count,left,right)),Age_prop_int13[34,6]*popn_prop,
                    if_else(year=="2019"&(is.na(count)|between(count,left,right)),Age_prop_int13[35,6]*popn_prop,age_6))))))))))))))))))))))))))))))))))),
        
         age_7=if_else(year=="1985"&(is.na(count)|between(count,left,right)),Age_prop_int13[1,7]*popn_prop,
                    if_else(year=="1986"&(is.na(count)|between(count,left,right)),Age_prop_int13[2,7]*popn_prop,
                    if_else(year=="1987"&(is.na(count)|between(count,left,right)),Age_prop_int13[3,7]*popn_prop,
                    if_else(year=="1988"&(is.na(count)|between(count,left,right)),Age_prop_int13[4,7]*popn_prop,
                    if_else(year=="1989"&(is.na(count)|between(count,left,right)),Age_prop_int13[5,7]*popn_prop,
                    if_else(year=="1990"&(is.na(count)|between(count,left,right)),Age_prop_int13[6,7]*popn_prop,
                    if_else(year=="1991"&(is.na(count)|between(count,left,right)),Age_prop_int13[7,7]*popn_prop,
                    if_else(year=="1992"&(is.na(count)|between(count,left,right)),Age_prop_int13[8,7]*popn_prop,
                    if_else(year=="1993"&(is.na(count)|between(count,left,right)),Age_prop_int13[9,7]*popn_prop,
                    if_else(year=="1994"&(is.na(count)|between(count,left,right)),Age_prop_int13[10,7]*popn_prop,  
                    if_else(year=="1995"&(is.na(count)|between(count,left,right)),Age_prop_int13[11,7]*popn_prop,
                    if_else(year=="1996"&(is.na(count)|between(count,left,right)),Age_prop_int13[12,7]*popn_prop,
                    if_else(year=="1997"&(is.na(count)|between(count,left,right)),Age_prop_int13[13,7]*popn_prop,
                    if_else(year=="1998"&(is.na(count)|between(count,left,right)),Age_prop_int13[14,7]*popn_prop,
                    if_else(year=="1999"&(is.na(count)|between(count,left,right)),Age_prop_int13[15,7]*popn_prop,
                    if_else(year=="2000"&(is.na(count)|between(count,left,right)),Age_prop_int13[16,7]*popn_prop,
                    if_else(year=="2001"&(is.na(count)|between(count,left,right)),Age_prop_int13[17,7]*popn_prop,
                    if_else(year=="2002"&(is.na(count)|between(count,left,right)),Age_prop_int13[18,7]*popn_prop,
                    if_else(year=="2003"&(is.na(count)|between(count,left,right)),Age_prop_int13[19,7]*popn_prop,
                    if_else(year=="2004"&(is.na(count)|between(count,left,right)),Age_prop_int13[20,7]*popn_prop,
                    if_else(year=="2005"&(is.na(count)|between(count,left,right)),Age_prop_int13[21,7]*popn_prop,
                    if_else(year=="2006"&(is.na(count)|between(count,left,right)),Age_prop_int13[22,7]*popn_prop,
                    if_else(year=="2007"&(is.na(count)|between(count,left,right)),Age_prop_int13[23,7]*popn_prop,
                    if_else(year=="2008"&(is.na(count)|between(count,left,right)),Age_prop_int13[24,7]*popn_prop,
                    if_else(year=="2009"&(is.na(count)|between(count,left,right)),Age_prop_int13[25,7]*popn_prop,
                    if_else(year=="2010"&(is.na(count)|between(count,left,right)),Age_prop_int13[26,7]*popn_prop,
                    if_else(year=="2011"&(is.na(count)|between(count,left,right)),Age_prop_int13[27,7]*popn_prop,
                    if_else(year=="2012"&(is.na(count)|between(count,left,right)),Age_prop_int13[28,7]*popn_prop,
                    if_else(year=="2013"&(is.na(count)|between(count,left,right)),Age_prop_int13[29,7]*popn_prop,
                    if_else(year=="2014"&(is.na(count)|between(count,left,right)),Age_prop_int13[30,7]*popn_prop,
                    if_else(year=="2015"&(is.na(count)|between(count,left,right)),Age_prop_int13[31,7]*popn_prop,
                    if_else(year=="2016"&(is.na(count)|between(count,left,right)),Age_prop_int13[32,7]*popn_prop,
                    if_else(year=="2017"&(is.na(count)|between(count,left,right)),Age_prop_int13[33,7]*popn_prop,
                    if_else(year=="2018"&(is.na(count)|between(count,left,right)),Age_prop_int13[34,7]*popn_prop,
                    if_else(year=="2019"&(is.na(count)|between(count,left,right)),Age_prop_int13[35,7]*popn_prop,age_7))))))))))))))))))))))))))))))))))))
        
        


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
