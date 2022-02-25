#Extract run timing by population

Timing <- rr_mod $rho_dst

Timing_pop1 <- Timing[,1,]
dim(Timing_pop1)
rownames(Timing_pop1) <- NULL

write.csv(Timing_pop2,"MigTiming_Popn1.csv")

Yday=(160:285)

Timing_pop2 <- bind_cols(Yday,Timing_pop1)

#Calculate migration timing by population and year

Timing_select <- select(Timing_pop2,1:2)
Reps1 <- rep(Timing_select$...1, (Timing_select$...2)*100)
Reps1 <- as.numeric(Reps1)
P1Y1_stats <- data.frame(min= min(Reps1),
                            q_25 = quantile(Reps1, 0.25),
                   q_50 = quantile(Reps1, 0.5),
                   q_75 = quantile(Reps1, 0.75),
                   max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,3)]
Reps1 <- rep(Timing_select$...1, Timing_select$...3*100)
Reps1 <- as.numeric(Reps1)
P1Y2_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,4)]
Reps1 <- rep(Timing_select$...1, Timing_select$...4*100)
Reps1 <- as.numeric(Reps1)
P1Y3_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,5)]
Reps1 <- rep(Timing_select$...1, Timing_select$...5*100)
Reps1 <- as.numeric(Reps1)
P1Y4_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,6)]
Reps1 <- rep(Timing_select$...1, Timing_select$...6*100)
Reps1 <- as.numeric(Reps1)
P1Y5_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,7)]
Reps1 <- rep(Timing_select$...1, Timing_select$...7*100)
Reps1 <- as.numeric(Reps1)
P1Y6_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,8)]
Reps1 <- rep(Timing_select$...1, Timing_select$...8*100)
Reps1 <- as.numeric(Reps1)
P1Y7_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,9)]
Reps1 <- rep(Timing_select$...1, Timing_select$...9*100)
Reps1 <- as.numeric(Reps1)
P1Y8_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,10)]
Reps1 <- rep(Timing_select$...1, Timing_select$...10*100)
Reps1 <- as.numeric(Reps1)
P1Y9_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,11)]
Reps1 <- rep(Timing_select$...1, Timing_select$...11*100)
Reps1 <- as.numeric(Reps1)
P1Y10_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,12)]
Reps1 <- rep(Timing_select$...1, Timing_select$...12*100)
Reps1 <- as.numeric(Reps1)
P1Y11_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,13)]
Reps1 <- rep(Timing_select$...1, Timing_select$...13*100)
Reps1 <- as.numeric(Reps1)
P1Y12_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,14)]
Reps1 <- rep(Timing_select$...1, Timing_select$...14*100)
Reps1 <- as.numeric(Reps1)
P1Y13_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,15)]
Reps1 <- rep(Timing_select$...1, Timing_select$...15*100)
Reps1 <- as.numeric(Reps1)
P1Y14_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,16)]
Reps1 <- rep(Timing_select$...1, Timing_select$...16*100)
Reps1 <- as.numeric(Reps1)
P1Y15_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,17)]
Reps1 <- rep(Timing_select$...1, Timing_select$...17*100)
Reps1 <- as.numeric(Reps1)
P1Y16_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,18)]
Reps1 <- rep(Timing_select$...1, Timing_select$...18*100)
Reps1 <- as.numeric(Reps1)
P1Y17_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,19)]
Reps1 <- rep(Timing_select$...1, Timing_select$...19*100)
Reps1 <- as.numeric(Reps1)
P1Y18_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,20)]
Reps1 <- rep(Timing_select$...1, Timing_select$...20*100)
Reps1 <- as.numeric(Reps1)
P1Y19_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,21)]
Reps1 <- rep(Timing_select$...1, Timing_select$...21*100)
Reps1 <- as.numeric(Reps1)
P1Y20_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,22)]
Reps1 <- rep(Timing_select$...1, Timing_select$...22*100)
Reps1 <- as.numeric(Reps1)
P1Y21_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,23)]
Reps1 <- rep(Timing_select$...1, Timing_select$...23*100)
Reps1 <- as.numeric(Reps1)
P1Y22_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,24)]
Reps1 <- rep(Timing_select$...1, Timing_select$...24*100)
Reps1 <- as.numeric(Reps1)
P1Y23_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,25)]
Reps1 <- rep(Timing_select$...1, Timing_select$...25*100)
Reps1 <- as.numeric(Reps1)
P1Y24_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,26)]
Reps1 <- rep(Timing_select$...1, Timing_select$...26*100)
Reps1 <- as.numeric(Reps1)
P1Y25_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,27)]
Reps1 <- rep(Timing_select$...1, Timing_select$...27*100)
Reps1 <- as.numeric(Reps1)
P1Y26_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,28)]
Reps1 <- rep(Timing_select$...1, Timing_select$...28*100)
Reps1 <- as.numeric(Reps1)
P1Y27_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,29)]
Reps1 <- rep(Timing_select$...1, Timing_select$...29*100)
Reps1 <- as.numeric(Reps1)
P1Y28_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,30)]
Reps1 <- rep(Timing_select$...1, Timing_select$...30*100)
Reps1 <- as.numeric(Reps1)
P1Y29_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,31)]
Reps1 <- rep(Timing_select$...1, Timing_select$...31*100)
Reps1 <- as.numeric(Reps1)
P1Y30_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,32)]
Reps1 <- rep(Timing_select$...1, Timing_select$...32*100)
Reps1 <- as.numeric(Reps1)
P1Y31_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,33)]
Reps1 <- rep(Timing_select$...1, Timing_select$...33*100)
Reps1 <- as.numeric(Reps1)
P1Y32_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,34)]
Reps1 <- rep(Timing_select$...1, Timing_select$...34*100)
Reps1 <- as.numeric(Reps1)
P1Y33_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,35)]
Reps1 <- rep(Timing_select$...1, Timing_select$...35*100)
Reps1 <- as.numeric(Reps1)
P1Y34_stats <- data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1)) 

Timing_select <- Timing_pop2[,c(1,36)]
Reps1 <- rep(Timing_select$...1, Timing_select$...36*100)
Reps1 <- as.numeric(Reps1)
P1Y35_stats <- (data.frame(min= min(Reps1),
                         q_25 = quantile(Reps1, 0.25),
                         q_50 = quantile(Reps1, 0.5),
                         q_75 = quantile(Reps1, 0.75),
                         max=max(Reps1))) 

Stat <- c("min","q_25","q_50","q_75","max")

P1_stats <- (rbind(P1Y1_stats,P1Y2_stats,P1Y3_stats,P1Y4_stats,P1Y5_stats,
                  P1Y6_stats,P1Y7_stats,P1Y8_stats,P1Y9_stats,P1Y10_stats,
                  P1Y11_stats,P1Y12_stats,P1Y13_stats,P1Y14_stats,P1Y15_stats,
                  P1Y16_stats,P1Y17_stats,P1Y18_stats,P1Y19_stats,P1Y20_stats,
                  P1Y21_stats,P1Y22_stats,P1Y23_stats,P1Y24_stats,P1Y25_stats,
                  P1Y26_stats,P1Y27_stats,P1Y28_stats,P1Y29_stats,P1Y30_stats,
                  P1Y31_stats,P1Y32_stats,P1Y33_stats,P1Y34_stats,P1Y35_stats))
rownames(P1_stats) <- NULL
P1_stats <- mutate(P1_stats,Year=rep(1985:2019),Popn="1")
P1_stats2 <- gather(P1_stats,"Stat","Yday",1:5)
colnames(P1Y35_stats) <- c("Stat","Yday")



#try loop over all years

Days <- select(Timing_pop2,1)
Years <- select(Timing_pop2,2)

Reps1 <- rep(Test1$...1, Test1$...2) 

Timing_pop2 <- as.matrix(Timing_pop2) 

Timing_pop3 <- split(Timing_pop2,rep(1:ncol(Timing_pop2),each=nrow(Timing_pop2)))

for(i in Timing_pop3[2:36]){
  print(i)
}

for(i in Timing_pop3[2:36]){
  Reps[i] <- rep(Timing_pop3[1], (Timing_pop3[i+1]))
  print(Reps)
}

for(i in Timing_pop2[,2:36]){
  Reps[i] <- as.numeric(unlist(rep(Timing_pop2[1], (Timing_pop2[i]))))
  print(Reps)
}


MT <- matrix(nrow=n.days,ncol=n.years)

Reps <- matrix(nrow=n.days,ncol=n.years)

p <- 1
for(p in 1:n.days) {
  y <- 1
  for(y in 1:n.years[p]) {
    year <- years[p,y]
    
    #Spawners
    Spawners[p,y] <- brood_table$S_med[brood_table$population==pops[p] & brood_table$BroodYear==year]
    
    Reps[p,y] <- rep(Timing_pop1[,1], d$Frequency)
    
    }#next y
}#next p


#try loop over years and populations

p <- 1
for(p in 1:n.days) {
  y <- 1
  for(y in 1:n.years[p]) {
    year <- years[p,y]
    
    #Spawners
    Spawners[p,y] <- brood_table$S_med[brood_table$population==pops[p] & brood_table$BroodYear==year]
    
    Reps[p,y] <- rep(Timing_pop1[,1], d$Frequency)
    
  }#next y
}#next p

d <- read.table(header = TRUE, text = "Score     Frequency
 100         10
 200         30
 300         40")

d2 <- rep(d$Score, d$Frequency)  ## expands the data by frequency of score

multi.fun <- function(x) {
  c(mean = mean(x), median = median(x), var = var(x), sd = sd(x))
}

multi.fun(d2)




alphabeta_df3_perc <- Timing_pop1 %>% 
  group_by(population, abund) %>%
  dplyr::summarise(q_05 = quantile(y, probs=0.05),
                   q_95 = quantile(y, probs=0.95),
                   q_50 = quantile(y, probs=0.5)) %>%
  as.data.frame()