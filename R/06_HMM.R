
library(momentuHMM)
library(tidyverse)
library(bayesmove)

source('R/Iterate HMMs.R')


######################################################
### Analyze simulations with unusual distributions ###
######################################################

### Prep Data ###

d<- read.csv("data/CRW_MM_sim_weird.csv", as.is = T)
names(d)[1]<- "ID"
d.list<- df_to_list(d, ind = "ID")

#format data
move.d<- map(d.list, ~prepData(.[,1:3], type = "UTM", coordNames = c("x","y")))

#Plot tracks
par(mfrow=c(2,2), ask=T)
for (i in 1:length(move.d)) {
  plot(y~x, data = move.d[[i]], main = paste("ID",names(d.list)[i]), type="o", col="gray65",
       asp=1, cex=0, pch=19, xlab="Longitude",ylab="Latitude")
}
par(mfrow=c(1,1), ask=F)



### Run HMMs ###

### Test with 2-4 behavioral states and then perform order selection via AIC/BIC
set.seed(3)
hmm.res_1k<- run.HMMs(move.d, 1:5)
names(hmm.res_1k)<- names(move.d)[1:5]

set.seed(3)
hmm.res_5k<- run.HMMs(move.d, 6:10)
hmm.res_5k<- hmm.res_5k[6:10]
names(hmm.res_5k)<- names(move.d)[6:10]

set.seed(1)
hmm.res_10k<- run.HMMs(move.d, 11:15)
hmm.res_10k<- hmm.res_10k[11:15]
names(hmm.res_10k)<- names(move.d)[11:15]

set.seed(1)
hmm.res_50k<- run.HMMs(move.d, 16:20)
hmm.res_50k<- hmm.res_50k[16:20]
names(hmm.res_50k)<- names(move.d)[16:20]


## Combine all model results
hmm.res<- c(hmm.res_1k, hmm.res_5k, hmm.res_10k, hmm.res_50k)



### Extract param values from fitted HMMs

hmm.params<- map(hmm.res, ~{map(.$models[2], getPar0) %>% 
    map(., ~pluck(., "Par")) %>% 
    flatten(.)
})

#step length params
hmm.step.params<- map(hmm.params, ~pluck(., "step")[1:6]) %>% 
  bind_rows() %>% 
  as.data.frame()

#turning angle params
hmm.angle.params<- map(hmm.params, ~pluck(., "angle")) %>% 
  bind_rows() %>% 
  as.data.frame()

#save model params
# write.csv(hmm.step.params, "data/HMM result step params_weird.csv")
# write.csv(hmm.angle.params, "data/HMM result angle params_weird.csv")




#Compare elapsed times
time<- map(hmm.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))





## View output from each model and inspect pseudo-residuals

for (i in 1:length(hmm.res)) {
  plot(hmm.res[[i]]$models[[1]])
  plot(hmm.res[[i]]$models[[2]])
  plot(hmm.res[[i]]$models[[3]])
}

par(ask=T)
for (i in 1:length(hmm.res)) {
  plotPR(hmm.res[[i]]$models[[1]])
  plotPR(hmm.res[[i]]$models[[2]])
  plotPR(hmm.res[[i]]$models[[3]])
}
par(ask=F)


## Make inference via AIC
for (i in 1:length(hmm.res)) {
  k.2<- hmm.res[[i]]$models[[1]]  # 2 states
  k.3<- hmm.res[[i]]$models[[2]]  # 3 states  
  k.4<- hmm.res[[i]]$models[[3]]  # 4 states
  
  print(names(hmm.res)[i])
  print(AIC(k.2, k.3, k.4))
  print(AICweights(k.2, k.3, k.4))  #4 states is far and away the best model
}

# Identify K per AIC
k.optim_AIC<- c(3, 3, 3, 3, 3,
                3, 3, 3, 3, 3,
                3, 3, 3, 3, 4,
                4, 4, 3, 3, 3)

## Make inference via BIC
# BIC = -2*logL + p*log(T)

for (i in 1:length(hmm.res)) {
  
  k.2<- hmm.res[[i]]$models[[1]]  # 2 states
  k.3<- hmm.res[[i]]$models[[2]]  # 3 states  
  k.4<- hmm.res[[i]]$models[[3]]  # 4 states
  
  print(names(hmm.res)[i])
  #K=2
  print(2*k.2$mod$minimum) + (length(k.2$mod$estimate)*log(nrow(d.list[[i]])))  
  #K=3
  print(2*k.3$mod$minimum) + (length(k.3$mod$estimate)*log(nrow(d.list[[i]])))  
  #K=4
  print(2*k.4$mod$minimum) + (length(k.4$mod$estimate)*log(nrow(d.list[[i]])))  
}

# Identify K per BIC
k.optim_BIC<- c(3, 3, 3, 3, 3,
                3, 3, 3, 3, 3,
                4, 3, 3, 3, 4,
                4, 4, 3, 3, 3)


# AIC suggested 3 sims w/ 4 states
# BIC suggested 4 sims w/ 4 states


### Sticking w/ K=3 for direct comparison

hmm.states<- hmm.res %>% 
  map(., ~pluck(., 1, 2)) %>% 
  map(viterbi) %>% 
  map2(d.list, ., ~cbind(.x, hmm.state = c(NA, .y[-length(.y)]))) %>% 
  bind_rows()




### Export Results ###

# write.csv(hmm.states, "data/HMM results_weird.csv", row.names = F)
# write.csv(time, "data/HMM_elapsed_time_weird.csv", row.names = F)  #units = min










#########################################################
### Analyze simulations with parametric distributions ###
#########################################################

### Prep Data ###

d<- read.csv("data/CRW_MM_sim_parametric.csv", as.is = T)
names(d)[1]<- "ID"
d.list<- df_to_list(d, ind = "ID")

#format data
move.d<- map(d.list, ~prepData(.[,1:3], type = "UTM", coordNames = c("x","y")))

#Plot tracks
par(mfrow=c(2,2), ask=T)
for (i in 1:length(move.d)) {
  plot(y~x, data = move.d[[i]], main = paste("ID",names(d.list)[i]), type="o", col="gray65",
       asp=1, cex=0, pch=19, xlab="Longitude",ylab="Latitude")
}
par(mfrow=c(1,1), ask=F)



### Run HMMs ###

### Test with 2-4 behavioral states and then perform order selection via AIC/BIC
set.seed(3)
hmm.res_1k<- run.HMMs(move.d, 1:5)
names(hmm.res_1k)<- names(move.d)[1:5]

set.seed(3)
hmm.res_5k<- run.HMMs(move.d, 6:10)
hmm.res_5k<- hmm.res_5k[6:10]
names(hmm.res_5k)<- names(move.d)[6:10]

set.seed(1)
hmm.res_10k<- run.HMMs(move.d, 11:15)
hmm.res_10k<- hmm.res_10k[11:15]
names(hmm.res_10k)<- names(move.d)[11:15]

set.seed(1)
hmm.res_50k<- run.HMMs(move.d, 16:20)
hmm.res_50k<- hmm.res_50k[16:20]
names(hmm.res_50k)<- names(move.d)[16:20]


## Combine all model results
hmm.res<- c(hmm.res_1k, hmm.res_5k, hmm.res_10k, hmm.res_50k)



### Extract param values from fitted HMMs

hmm.params<- map(hmm.res, ~{map(.$models[2], getPar0) %>% 
    map(., ~pluck(., "Par")) %>% 
    flatten(.)
})

#step length params
hmm.step.params<- map(hmm.params, ~pluck(., "step")[1:6]) %>% 
  bind_rows() %>% 
  as.data.frame()

#turning angle params
hmm.angle.params<- map(hmm.params, ~pluck(., "angle")) %>% 
  bind_rows() %>% 
  as.data.frame()

#save model params
# write.csv(hmm.step.params, "data/HMM result step params_parametric.csv")
# write.csv(hmm.angle.params, "data/HMM result angle params_parametric.csv")




#Compare elapsed times
time<- map(hmm.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))





## View output from each model and inspect pseudo-residuals

for (i in 1:length(hmm.res)) {
  plot(hmm.res[[i]]$models[[1]])
  plot(hmm.res[[i]]$models[[2]])
  plot(hmm.res[[i]]$models[[3]])
}

par(ask=T)
for (i in 1:length(hmm.res)) {
  plotPR(hmm.res[[i]]$models[[1]])
  plotPR(hmm.res[[i]]$models[[2]])
  plotPR(hmm.res[[i]]$models[[3]])
}
par(ask=F)


## Make inference via AIC
for (i in 1:length(hmm.res)) {
  k.2<- hmm.res[[i]]$models[[1]]  # 2 states
  k.3<- hmm.res[[i]]$models[[2]]  # 3 states  
  k.4<- hmm.res[[i]]$models[[3]]  # 4 states
  
  print(names(hmm.res)[i])
  print(AIC(k.2, k.3, k.4))
  print(AICweights(k.2, k.3, k.4))  #4 states is far and away the best model
}

# Identify K per AIC
k.optim_AIC<- c(3, 3, 3, 3, 3,
                3, 3, 3, 3, 3,
                3, 3, 4, 4, 4,
                3, 3, 3, 4, 3)

## Make inference via BIC
# BIC = -2*logL + p*log(T)

for (i in 1:length(hmm.res)) {
  
  k.2<- hmm.res[[i]]$models[[1]]  # 2 states
  k.3<- hmm.res[[i]]$models[[2]]  # 3 states  
  k.4<- hmm.res[[i]]$models[[3]]  # 4 states
  
  print(names(hmm.res)[i])
  #K=2
  print(2*k.2$mod$minimum) + (length(k.2$mod$estimate)*log(nrow(d.list[[i]])))  
  #K=3
  print(2*k.3$mod$minimum) + (length(k.3$mod$estimate)*log(nrow(d.list[[i]])))  
  #K=4
  print(2*k.4$mod$minimum) + (length(k.4$mod$estimate)*log(nrow(d.list[[i]])))  
}

# Identify K per BIC
k.optim_BIC<- c(3, 3, 3, 3, 3,
                3, 3, 3, 3, 3,
                3, 3, 4, 4, 4,
                3, 3, 3, 4, 3)

# AIC suggested 4 sims w/ 4 states
# BIC suggested 4 sims w/ 4 states


### Sticking w/ K=3 for direct comparison

hmm.states<- hmm.res %>% 
  map(., ~pluck(., 1, 2)) %>% 
  map(viterbi) %>% 
  map2(d.list, ., ~cbind(.x, hmm.state = c(NA, .y[-length(.y)]))) %>% 
  bind_rows()




### Export Results ###

# write.csv(hmm.states, "data/HMM results_parametric.csv", row.names = F)
# write.csv(time, "data/HMM_elapsed_time_parametric.csv", row.names = F)  #units = min
