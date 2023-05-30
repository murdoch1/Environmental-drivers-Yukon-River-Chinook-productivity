## Sensitivity analysis processing ###

library(tidyverse)

integrated_pos_int1 <- as.data.frame(readRDS("Data/seperated-posterior.rds"))

integrated_pos_int2 <- integrated_pos_int1 %>% 
  select(starts_with(c("S[","R[")))

#summaries

integrated_pos_summary <- summary(integrated_pos_int2)

integrated_pos_int3 <- integrated_pos_int2 %>% 
  sample_n(100) %>% as.data.frame()

#select spawners

integrated_pos_int4 <- integrated_pos_int3 %>% select(1:280)

integrated_pos_int4b <- integrated_pos_int4[,c(1,12,23,30:35,2:11,13:22,24:29,
                                               36,47,58,65:70,37:46,48:57,59:64,
                                               71,82,93,100:105,72:81,83:92,94:99,
                                               106,117,128,135:140,107:116,118:127,129:134,
                                               141,152,163,170:175,142:151,153:162,164:169,
                                               176,187,198,205:210,177:186,188:197,199:204,
                                               211,222,233,240:245,212:221,223:232,234:239,
                                               246,257,268,275:280,247:256,258:267,269:274)]

  
integrated_pos_int5 <- as.data.frame(t(integrated_pos_int4b)) %>% 
  mutate(year=rep(1985:2019,times=8)) %>% 
  mutate(population=rep(c("Lower Mainstem","White-Donjek","Pelly","Stewart","Carmacks","Teslin",
                      "Middle Mainstem","Upper Lakes and Mainstem"),each=35)) %>%
  gather(sample,spawners,1:100) %>% mutate(sample=rep(1:100,each=280)) %>% 
  filter(year<2013)

#process recruits and then match


integrated_pos_int6 <- integrated_pos_int3 %>% select(281:584)

integrated_pos_int6b <- integrated_pos_int6[,c(1,12,23,33:38,2:11,13:22,24:32,
                                               39,50,61,71:76,40:49,51:60,62:70,
                                               77,88,99,109:114,78:87,89:98,100:108,
                                               115,126,137,147:152,116:125,127:136,138:146,
                                               153,164,175,185:190,154:163,165:174,176:184,
                                               191,202,213,223:228,192:201,203:212,214:222,
                                               229,240,251,261:266,230:239,241:250,252:260,
                                               267,278,289,299:304,268:277,279:288,290:298)]

integrated_pos_int7 <- as.data.frame(t(integrated_pos_int6b)) %>% 
  mutate(year=rep(1978:2015,times=8)) %>% 
  mutate(population=rep(c("Lower Mainstem","White-Donjek","Pelly","Stewart","Carmacks","Teslin",
                          "Middle Mainstem","Upper Lakes and Mainstem"),each=38)) %>% 
  gather(sample,recruits,1:100) %>% mutate(sample=rep(1:100,each=304)) %>% 
  filter(year>1984&year<2013)

#join

SR_data <- left_join(integrated_pos_int5,integrated_pos_int7)
SR_data$sample <- as.factor(SR_data$sample)

#Spawner array
Spawner_int1 <- SR_data %>% select(1:4) %>% spread(year,spawners)

Spawner_int2 <- split(Spawner_int1,Spawner_int1$sample,drop=FALSE)

library(abind)
Spawner_int3 <- abind(Spawner_int2,along=3)

spawner_array <- Spawner_int3[1:8,3:30,1:100]

spawner_array <- as.numeric(spawner_array)

spawner_array_num <- array(as.numeric(unlist(spawner_array)),dim=c(8,28,100))
spawner_array_num

save(spawner_array_num,file="Data/spawner_array_sep.Rdata")

#Recruit array

Recruit_int1 <- SR_data %>% select(1:3,5) %>% spread(year,recruits)

Recruit_int2 <- split(Recruit_int1,Recruit_int1$sample,drop=FALSE)

Recruit_int3 <- abind(Recruit_int2,along=3)
recruit_array <- Recruit_int3[1:8,3:30,1:100]

recruit_array <- as.numeric(recruit_array)

recruit_array_num <- array(as.numeric(unlist(recruit_array)),dim=c(8,28,100))
recruit_array_num

save(recruit_array_num,file="Data/recruit_array_sep.Rdata")

#popn order index is Carmacks, Lower Main, Middle Main, Pelly, Stewart, Teslin, Upper Lakes, White-Donjek




# 
# 
# x2 <- array(unlist(x),dim=c(224,5,100),dimnames=list(c.names,r.names,m.names))
# 
# r.names <- c("year","population","sample","S","R")
# c.names <- 1:224
# m.names <- 1:100
# 
# x3 <- x2
# 
# x3 <- apply(x2,1:3,as.numeric)
# 
# 
# spawners_long_int1 <- x2[1:224,c(1:2,4),1:100]
# library(reshape2)
# spawners_long_int2 <- acast(spawners_long_int1,population~year,value.var="S")
# 
# 
# 
# z <- array(1:24, dim = 2:4)
# zseq <- apply(z, 1:2, function(x) seq_len(max(x)))
# zseq         ## a 2 x 3 matrix

# #need to figure out how to loop all of this
# 
# pops = sort(unique(brood_table$population))
# n.pops = length(pops)
# n.years = vector(length=n.pops)
# years <- matrix(nrow=n.pops,ncol=28)
# 
# for(p in 1:n.pops) {
#   n.years[p] <- length(unique(brood_table$BroodYear[brood_table$population==pops[p]]))
#   years[p,1:n.years[p]] <- sort(unique(brood_table$BroodYear[brood_table$population==pops[p]]))
# }
# 
# 
# n.datasets=100
# 
# spawner_array <- array(data=NA,dim=c(8,28,100))
# 
# for (i in 1:n.datasets){
#   spawner_array <- integrated_pos_int2 %>% 
#     sample_n(100) %>% select(305:584) %>% 
#    t %>% 
#     mutate(year=rep(1985:2019,times=8)) %>% 
#     mutate(population=rep(c("Lower Mainstem","White-Donjek","Pelly","Stewart","Carmacks","Teslin",
#                             "Middle Mainstem","Upper Lakes and Mainstem"),each=35)) %>% 
#     filter(year<2013) %>% spread(year,V1) %>% select(2:29)
#   
# }
# 
