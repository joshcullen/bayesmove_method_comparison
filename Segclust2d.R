
library(segclust2d)
library(tidyverse)
library(bayesmove)
library(furrr)
library(future)
library(progressr)

source('helper functions.R')
source('Iterate segclust2d.R')


######################################################
### Analyze simulations with unusual distributions ###
######################################################

### Prep Data ###

dat<- read.csv("CRW_MM_sim_weird.csv", as.is = T)
true.brkpts<- read.csv("CRW_MM_sim_brkpts_weird.csv", as.is = T)

dat<- dat %>% 
  drop_na(TA) %>%   #remove any missing values
  mutate(abs.TA = abs(TA))  #create new var for absolute value of TA
dat.list<- df_to_list(dat, ind = "id")



### Run Segclust2d ###

### Test with 2-4 behavioral states (Lavielle method performs order selection via BIC)

plan(multisession, workers = 5)
handlers(handler_progress(format = "[:bar] :percent in :elapsed",
                          enable = TRUE,
                          clear = FALSE))

# takes 1 min
with_progress({
  p <- progressor(steps = length(dat.list[1:5]))
  segclust.res_1k<- future_map(dat.list[1:5], run.segclust2d_map,
                               .options = furrr_options(seed = TRUE))
})

# takes 1 hr
with_progress({
  p <- progressor(steps = length(dat.list[6:10]))
  segclust.res_5k<- future_map(dat.list[6:10], run.segclust2d_map,
                               .options = furrr_options(seed = TRUE))
})

# takes x hrs
with_progress({
  p <- progressor(steps = length(dat.list[11:15]))
  segclust.res_10k<- future_map(dat.list[11:15], run.segclust2d_map,
                               .options = furrr_options(seed = TRUE))
})

# with_progress({
#   p <- progressor(steps = length(dat.list[16:20]))
#   segclust.res_50k<- future_map(dat.list[16:20], run.segclust2d_map,
#                                .options = furrr_options(seed = TRUE))
# })

plan(sequential)

## Combine all model results
segclust.res<- c(segclust.res_1k, segclust.res_5k, segclust.res_10k)



segclust.mod<- map(segclust.res, ~pluck(., "models"))
# segclust.time<- map(segclust.res, ~pluck(., "elapsed.time")) %>% 
#   unlist()


plot(segclust.mod$`1_1`)
# plot(segclust.mod[[i]], ncluster = 3)
# plot(segclust.mod[[i]], ncluster = 3, nseg = 10)
# plot_likelihood(segclust.mod[[i]])
plot_BIC(segclust.mod$`1_1`)
# segmap(segclust.mod[[i]])
# stateplot(segclust.mod[[i]])


boop$`1_1`$outputs$`3 class - 10 segments`
augment(boop$`1_1`, ncluster = 3, nseg = 10)





# Find ncluster and nseg that maximize BIC per ID/nclust
segclust.BIC<- map(segclust.mod, ~pluck(., "BIC")) %>% 
  bind_rows(.id = "id") %>% 
  group_by(id, ncluster) %>% 
  summarize(max.BIC = max(BIC)) %>% 
  group_by(id) %>% 
  mutate(d.BIC = max(max.BIC) - max.BIC)

# Identify nclust per BIC
nclust.optim_BIC<- c(3, 3, 3, 3, 4,
                     4, 4, 4, 4, 4,
                     4, 4, 4, 3, 3)



#Compare elapsed times
time<- map(segclust.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k'))




# View output from each model; doesn't currently plot output
par(ask=T)
for (i in 1:length(segclust.mod)) {
  plot(segclust.mod[[i]])
  # plot(segclust.mod[[i]], ncluster = 3)
  # plot(segclust.mod[[i]], ncluster = 3, nseg = 10)
  # plot_likelihood(segclust.mod[[i]])
  plot_BIC(segclust.mod[[i]])
  # segmap(segclust.mod[[i]])
  # stateplot(segclust.mod[[i]])
}
par(ask=F)



### Extract param values from fitted models

segclust.state.params<- map(segclust.mod, extract.segclust2d.behav.params, ncluster = 3) %>% 
  bind_rows(.id = "id")

#save model params
# write.csv(segclust.state.params, "Segclust result state params_weird.csv")




### Extract breakpoints from best model

# Augment state estimates
segclust.aug<- segclust.mod %>% 
  map2(., as.list(nclust.optim_BIC), .f = ~augment(.x, ncluster = .y))

# Extract breakpoints
brks<- list()
for (i in 1:length(segclust.aug)) {
  brks[[i]]<- which(diff(segclust.aug[[i]]$state, lag = 1) != 0)
}


# Compare true vs modeled breakpoints
all.brkpts<- list()
for (i in 1:length(brks)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = brks[[i]],
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}

names(all.brkpts)<- names(dat.list)[1:15]
all.brkpts<- bind_rows(all.brkpts, .id = 'id')




### Sticking w/ ncluster=3 for all simulations for direct comparison

# Augment state estimates
segclust.aug2<- segclust.mod %>% 
  map(., augment, ncluster = 3) %>% 
  map(., 
      ~rbind(c(unique(.$id), 0, 0, rep(NA, 4), unique(.$track_length), rep(NA, 10)),
             .x)) %>% 
  bind_rows()



######################
### Export Results ###
######################

# write.csv(all.brkpts, "Segclust2d allbrkpts_weird.csv", row.names = F)
# write.csv(segclust.aug2, "Segclust2d results_weird.csv", row.names = F)
# write.csv(time, "Segclust2d_elapsed_time_weird.csv", row.names = F)  #units = min








#########################################################
### Analyze simulations with parametric distributions ###
#########################################################

### Prep Data ###

dat<- read.csv("CRW_MM_sim_parametric.csv", as.is = T)
true.brkpts<- read.csv("CRW_MM_sim_brkpts_parametric.csv", as.is = T)

dat<- dat %>% 
  drop_na(TA) %>%   #remove any missing values
  mutate(abs.TA = abs(TA))  #create new var for absolute value of TA
dat.list<- df_to_list(dat, ind = "id")



### Run Segclust2d ###

### Test with 2-4 behavioral states (Lavielle method performs order selection via BIC)

plan(multisession, workers = 5)
handlers(handler_progress(format = "[:bar] :percent in :elapsed",
                          enable = TRUE,
                          clear = FALSE))

# takes 1 min to run
with_progress({
  p <- progressor(steps = length(dat.list[1:5]))
  segclust.res_1k<- future_map(dat.list[1:5], run.segclust2d_map,
                               .options = furrr_options(seed = TRUE))
})

# takes 1 hr to run
with_progress({
  p <- progressor(steps = length(dat.list[6:10]))
  segclust.res_5k<- future_map(dat.list[6:10], run.segclust2d_map,
                               .options = furrr_options(seed = TRUE))
})

# takes 7 hrs to run
with_progress({
  p <- progressor(steps = length(dat.list[11:15]))
  segclust.res_10k<- future_map(dat.list[11:15], run.segclust2d_map,
                                .options = furrr_options(seed = TRUE))
})

# with_progress({
#   p <- progressor(steps = length(dat.list[16:20]))
#   segclust.res_50k<- future_map(dat.list[16:20], run.segclust2d_map,
#                                .options = furrr_options(seed = TRUE))
# })

plan(sequential)

## Combine all model results
segclust.res<- c(segclust.res_1k, segclust.res_5k, segclust.res_10k)

segclust.mod<- map(segclust.res, ~pluck(., "models"))



plot(segclust.mod$`5_3`)
# plot(segclust.mod[[i]], ncluster = 3)
# plot(segclust.mod[[i]], ncluster = 3, nseg = 10)
# plot_likelihood(segclust.mod[[i]])
plot_BIC(segclust.mod$`5_3`)
# segmap(segclust.mod[[i]])
# stateplot(segclust.mod[[i]])





# Find ncluster and nseg that maximize BIC per ID/nclust
segclust.BIC<- map(segclust.mod, ~pluck(., "BIC")) %>% 
  bind_rows(.id = "id") %>% 
  group_by(id, ncluster) %>% 
  summarize(max.BIC = max(BIC)) %>% 
  group_by(id) %>% 
  mutate(d.BIC = max(max.BIC) - max.BIC)

# Identify nclust per BIC
nclust.optim_BIC<- c(3, 3, 3, 4, 3,
                     4, 4, 4, 3, 4,
                     4, 4, 4, 4, 4)



#Compare elapsed times
time<- map(segclust.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k'))





### Extract param values from fitted models

segclust.state.params<- map(segclust.mod, extract.segclust2d.behav.params, ncluster = 3) %>% 
  bind_rows(.id = "id")

#save model params
# write.csv(segclust.state.params, "Segclust result state params_parametric.csv")




### Extract breakpoints from best model

# Augment state estimates
segclust.aug<- segclust.mod %>% 
  map2(., as.list(nclust.optim_BIC), .f = ~augment(.x, ncluster = .y))

# Extract breakpoints
brks<- list()
for (i in 1:length(segclust.aug)) {
  brks[[i]]<- which(diff(segclust.aug[[i]]$state, lag = 1) != 0)
}


# Compare true vs modeled breakpoints
all.brkpts<- list()
for (i in 1:length(brks)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = brks[[i]],
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}

names(all.brkpts)<- names(dat.list)[1:15]
all.brkpts<- bind_rows(all.brkpts, .id = 'id')




### Sticking w/ ncluster=3 for all simulations for direct comparison

# Augment state estimates
segclust.aug2<- segclust.mod %>% 
  map(., augment, ncluster = 3) %>% 
  map(., 
      ~rbind(c(unique(.$id), 0, 0, rep(NA, 4), unique(.$track_length), rep(NA, 10)),
             .x)) %>% 
  bind_rows()



######################
### Export Results ###
######################

# write.csv(all.brkpts, "Segclust2d allbrkpts_parametric.csv", row.names = F)
# write.csv(segclust.aug2, "Segclust2d results_parametric.csv", row.names = F)
# write.csv(time, "Segclust2d_elapsed_time_parametric.csv", row.names = F)  #units = min