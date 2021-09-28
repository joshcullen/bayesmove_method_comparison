

########################
#### Run EMbC Model ####
########################

library(tidyverse)
library(bayesmove)
library(future)
library(furrr)
library(progressr)
library(EMbC)

source('R/Iterate EMbC.R')


######################################################
### Analyze simulations with unusual distributions ###
######################################################

### Prep Data ###

dat<- read.csv("data/CRW_MM_sim_weird.csv", as.is = T)
dat$TA.abs<- abs(dat$TA)  #for proper analysis by EMbC

dat.list<- df_to_list(dat, ind = "id") %>% 
  map(., ~dplyr::select(., SL, TA.abs)) %>%   #filter columns for EMbC model
  map(., ~slice(., -1))  #remove first row w/ NAs


### Run EMbC ###

handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

with_progress({
  p <- progressor(steps = length(dat.list[1:5]))
  embc.res_1k<- future_map(dat.list[1:5], run.embc,
                               .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[6:10]))
  embc.res_5k<- future_map(dat.list[6:10], run.embc,
                           .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[11:15]))
  embc.res_10k<- future_map(dat.list[11:15], run.embc,
                           .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[16:20]))
  embc.res_50k<- future_map(dat.list[16:20], run.embc,
                            .options = furrr_options(seed = TRUE))
})

plan(sequential)


## Combine all model results
embc.res<- c(embc.res_1k, embc.res_5k, embc.res_10k, embc.res_50k)

embc.mod<- map(embc.res, ~pluck(., "model"))


## Plot results
par(ask = T)
for (i in 1:length(embc.mod)) {
  lkhp(embc.mod[[i]])
}

for (i in 1:length(embc.mod)) {
  sctr(embc.mod[[i]])
}
par(ask = F)



### Extract param values from fitted models

embc.state.params<- map(embc.mod, extract.embc.params) %>% 
  bind_rows(.id = "id")

#save model params
# write.csv(embc.state.params, "data/EMbC result state params_weird.csv")



### Add state estimates to observations

dat2<- dat %>% 
  drop_na(TA)
dat2$state<- map(embc.mod, pluck, "A") %>% 
  unlist() %>% 
  recode(., `1` = "LL",
         `2` = "LH",
         `3` = "HL",
         `4` = "HH",
         `5` = "NC")



#Compare elapsed times
time<- map(embc.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))



######################
### Export Results ###
######################

# write.csv(dat2, "data/EMbC results_weird.csv", row.names = F)
# write.csv(time, "data/EMbC_elapsed_time_weird.csv", row.names = F)  #units = min







#########################################################
### Analyze simulations with parametric distributions ###
#########################################################

### Prep Data ###

dat<- read.csv("data/CRW_MM_sim_parametric.csv", as.is = T)
dat$TA.abs<- abs(dat$TA)  #for proper analysis by EMbC

dat.list<- df_to_list(dat, ind = "id") %>% 
  map(., ~dplyr::select(., SL, TA.abs)) %>%   #filter columns for EMbC model
  map(., ~slice(., -1))  #remove first row w/ NAs


### Run EMbC ###


handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

with_progress({
  p <- progressor(steps = length(dat.list[1:5]))
  embc.res_1k<- future_map(dat.list[1:5], run.embc,
                           .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[6:10]))
  embc.res_5k<- future_map(dat.list[6:10], run.embc,
                           .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[11:15]))
  embc.res_10k<- future_map(dat.list[11:15], run.embc,
                            .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[16:20]))
  embc.res_50k<- future_map(dat.list[16:20], run.embc,
                            .options = furrr_options(seed = TRUE))
})

plan(sequential)


## Combine all model results
embc.res<- c(embc.res_1k, embc.res_5k, embc.res_10k, embc.res_50k)

embc.mod<- map(embc.res, ~pluck(., "model"))


## Plot results
par(ask = T)
for (i in 1:length(embc.mod)) {
  lkhp(embc.mod[[i]])
}

for (i in 1:length(embc.mod)) {
  sctr(embc.mod[[i]])
}
par(ask = F)



### Extract param values from fitted models

embc.state.params<- map(embc.mod, extract.embc.params) %>% 
  bind_rows(.id = "id")

#save model params
# write.csv(embc.state.params, "data/EMbC result state params_parametric.csv")



### Add state estimates to observations

dat2<- dat %>% 
  drop_na(TA)
dat2$state<- map(embc.mod, pluck, "A") %>% 
  unlist() %>% 
  recode(., `1` = "LL",
         `2` = "LH",
         `3` = "HL",
         `4` = "HH",
         `5` = "NC")



#Compare elapsed times
time<- map(embc.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))



######################
### Export Results ###
######################

# write.csv(dat2, "data/EMbC results_parametric.csv", row.names = F)
# write.csv(time, "data/EMbC_elapsed_time_parametric.csv", row.names = F)  #units = min











##########################################################
### Analyze HMM simulations with unusual distributions ###
##########################################################

### Prep Data ###

dat<- read.csv("data/HMM_sim_weird.csv", as.is = T)
dat$TA.abs<- abs(dat$TA)  #for proper analysis by EMbC

dat.list<- df_to_list(dat, ind = "id") %>% 
  map(., ~dplyr::select(., SL, TA.abs)) %>%   #filter columns for EMbC model
  map(., ~slice(., -1))  #remove first row w/ NAs


### Run EMbC ###

handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

with_progress({
  p <- progressor(steps = length(dat.list[1:5]))
  embc.res_1k<- future_map(dat.list[1:5], run.embc,
                           .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[6:10]))
  embc.res_5k<- future_map(dat.list[6:10], run.embc,
                           .options = furrr_options(seed = TRUE))
})

plan(sequential)


## Combine all model results
embc.res<- c(embc.res_1k, embc.res_5k)

embc.mod<- map(embc.res, ~pluck(., "model"))


## Plot results
par(ask = T)
for (i in 1:length(embc.mod)) {
  lkhp(embc.mod[[i]])
}

for (i in 1:length(embc.mod)) {
  sctr(embc.mod[[i]])
}
par(ask = F)



### Extract param values from fitted models

embc.state.params<- map(embc.mod, extract.embc.params) %>% 
  bind_rows(.id = "id")

#save model params
# write.csv(embc.state.params, "data/EMbC(HMM) result state params_weird.csv")



### Add state estimates to observations

dat2<- dat %>% 
  drop_na(TA)
dat2$state<- map(embc.mod, pluck, "A") %>% 
  unlist() %>% 
  recode(., `1` = "LL",
         `2` = "LH",
         `3` = "HL",
         `4` = "HH",
         `5` = "NC")



#Compare elapsed times
time<- map(embc.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k'), each = 5) %>% 
  factor(., levels = c('1k','5k'))



######################
### Export Results ###
######################

# write.csv(dat2, "data/EMbC(HMM) results_weird.csv", row.names = F)
# write.csv(time, "data/EMbC(HMM)_elapsed_time_weird.csv", row.names = F)  #units = min











#########################################################
### Analyze HMM simulations with common distributions ###
#########################################################

### Prep Data ###

dat<- read.csv("data/HMM_sim.csv", as.is = T)
dat$TA.abs<- abs(dat$TA)  #for proper analysis by EMbC

dat.list<- df_to_list(dat, ind = "id") %>% 
  map(., ~dplyr::select(., SL, TA.abs)) %>%   #filter columns for EMbC model
  map(., ~slice(., -1))  #remove first row w/ NAs


### Run EMbC ###

handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

with_progress({
  p <- progressor(steps = length(dat.list[1:5]))
  embc.res_1k<- future_map(dat.list[1:5], run.embc,
                           .options = furrr_options(seed = TRUE))
})

with_progress({
  p <- progressor(steps = length(dat.list[6:10]))
  embc.res_5k<- future_map(dat.list[6:10], run.embc,
                           .options = furrr_options(seed = TRUE))
})

plan(sequential)


## Combine all model results
embc.res<- c(embc.res_1k, embc.res_5k)

embc.mod<- map(embc.res, ~pluck(., "model"))


## Plot results
par(ask = T)
for (i in 1:length(embc.mod)) {
  lkhp(embc.mod[[i]])
}

for (i in 1:length(embc.mod)) {
  sctr(embc.mod[[i]])
}
par(ask = F)



### Extract param values from fitted models

embc.state.params<- map(embc.mod, extract.embc.params) %>% 
  bind_rows(.id = "id")

#save model params
# write.csv(embc.state.params, "data/EMbC(HMM) result state params_parametric.csv")



### Add state estimates to observations

dat2<- dat %>% 
  drop_na(TA)
dat2$state<- map(embc.mod, pluck, "A") %>% 
  unlist() %>% 
  recode(., `1` = "LL",
         `2` = "LH",
         `3` = "HL",
         `4` = "HH",
         `5` = "NC")



#Compare elapsed times
time<- map(embc.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k'), each = 5) %>% 
  factor(., levels = c('1k','5k'))



######################
### Export Results ###
######################

# write.csv(dat2, "data/EMbC(HMM) results_parametric.csv", row.names = F)
# write.csv(time, "data/EMbC(HMM)_elapsed_time_parametric.csv", row.names = F)  #units = min
