library(bcpa)
library(tidyverse)

source('R/helper functions.R')

set.seed(1)



###################################################################
### Analyze segment-based simulations w/ uncommon distributions ###
###################################################################


### Prep Data ###

d<- read.csv("data/CRW_MM_sim_weird.csv", as.is = T)
true.brkpts<- read.csv("data/CRW_MM_sim_brkpts_weird.csv", as.is = T)

d.list<- bayesmove::df_to_list(d, ind = "id")
d.list<- map(d.list, ~mutate(., time = seq(c(ISOdate(2020,4,29)), by = "hour",
                                           length.out = n())))  #create fake hourly times 

mytrack<- map(d.list, ~MakeTrack(.$x, .$y, .$time))

#Plot tracks
par(mfrow=c(2,2), ask=T)
for (i in 1:length(mytrack)) {
  plot(mytrack[[i]], main = paste("ID",names(d.list)[i]))
}
par(mfrow=c(1,1), ask=F)



### Run BCPA ###

Simp1<- map(mytrack, GetVT)

Simp.ws<- list()
elapsed.time<- list()
cp<- list()
cw=30
#clusterwidth = the number of times a changepoint must be identified in the moving window for it to count. 1=all changepoints
#windowsize = # points to consider each time the model is run, k=2: default sensitivity (lower=less sensitive)

for (i in 1:length(Simp1)) {
  start.time<- Sys.time()
  
  Simp.ws[[i]]<- WindowSweep(Simp1[[i]], "V*cos(Theta)", windowsize=80, progress=T, K=2)
  cp[[i]]<-ChangePointSummary(Simp.ws[[i]], clusterwidth=cw)
  
  end.time<- Sys.time()
  elapsed.time[[i]]<- difftime(end.time, start.time, units = "min")
}

names(elapsed.time)<- names(d.list)

#Plot results
par(mfrow=c(2,2), ask=T)
for (i in 1:length(Simp.ws)) {
  plot(Simp.ws[[i]], type="flat", clusterwidth=cw, main = paste("ID",names(d.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Plot diagnostics
par(ask=T)
for (i in 1:length(Simp.ws)) {
  DiagPlot(Simp.ws[[i]])
}
par(ask=F)


# Store model results in DFs
y<- list()
for (i in 1:length(d.list)) {
  y[[i]]<- cbind(d.list[[i]][3:nrow(d.list[[i]]),],
                 t = Simp.ws[[i]]$t)  #number of lines in Simp.ws
  
  y[[i]]$group<-0
  len<- length(cp[[i]]$phases$t0)
  
  for (j in 1:len){ #at least the number of changepoints
    ind1 = cp[[i]]$phases$t0[j] 
    ind2 = cp[[i]]$phases$t1[j] 
    test1 = (y[[i]]$t>=ind1 & y[[i]]$t<=ind2)
    y[[i]]$group[test1]=j  #assign time segments to obs
  }
}



## Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:length(cp)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = round(cp[[i]]$breaks$middle, 0),
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}



#Compare elapsed time
time<- map_dfr(elapsed.time, `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
         time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))



### Export Results ###

names(all.brkpts)<- names(d.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')


# write.csv(all.brkpts, "data/BCPA_allbrkpts_weird.csv", row.names = F)
# write.csv(time, "data/BCPA_elapsed_time_weird.csv", row.names = F)  #units = min







###################################################################
### Analyze segment-based simulations w/ common distributions ###
###################################################################

set.seed(1)

### Prep Data ###

d<- read.csv("data/CRW_MM_sim_parametric.csv", as.is = T)
true.brkpts<- read.csv("data/CRW_MM_sim_brkpts_parametric.csv", as.is = T)

d.list<- bayesmove::df_to_list(d, ind = "id")
d.list<- map(d.list, ~mutate(., time = seq(c(ISOdate(2020,4,29)), by = "hour",
                                           length.out = n())))  #create fake hourly times 

mytrack<- map(d.list, ~MakeTrack(.$x, .$y, .$time))

#Plot tracks
par(mfrow=c(2,2), ask=T)
for (i in 1:length(mytrack)) {
  plot(mytrack[[i]], main = paste("ID",names(d.list)[i]))
}
par(mfrow=c(1,1), ask=F)



### Run BCPA ###

Simp1<- map(mytrack, GetVT)

Simp.ws<- list()
elapsed.time<- list()
cp<- list()
cw=30
#clusterwidth = the number of times a changepoint must be identified in the moving window for it to count. 1=all changepoints
#windowsize = # points to consider each time the model is run, k=2: default sensitivity (lower=less sensitive)

for (i in 1:length(Simp1)) {
  start.time<- Sys.time()
  
  Simp.ws[[i]]<- WindowSweep(Simp1[[i]], "V*cos(Theta)", windowsize=80, progress=T, K=2)
  cp[[i]]<-ChangePointSummary(Simp.ws[[i]], clusterwidth=cw)
  
  end.time<- Sys.time()
  elapsed.time[[i]]<- difftime(end.time, start.time, units = "min")
}

names(elapsed.time)<- names(d.list)

#Plot results
par(mfrow=c(2,2), ask=T)
for (i in 1:length(Simp.ws)) {
  plot(Simp.ws[[i]], type="flat", clusterwidth=cw, main = paste("ID",names(d.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Plot diagnostics
par(ask=T)
for (i in 1:length(Simp.ws)) {
  DiagPlot(Simp.ws[[i]])
}
par(ask=F)


# Store model results in DFs
y<- list()
for (i in 1:length(d.list)) {
  y[[i]]<- cbind(d.list[[i]][3:nrow(d.list[[i]]),],
                 t = Simp.ws[[i]]$t)  #number of lines in Simp.ws
  
  y[[i]]$group<-0
  len<- length(cp[[i]]$phases$t0)
  
  for (j in 1:len){ #at least the number of changepoints
    ind1 = cp[[i]]$phases$t0[j] 
    ind2 = cp[[i]]$phases$t1[j] 
    test1 = (y[[i]]$t>=ind1 & y[[i]]$t<=ind2)
    y[[i]]$group[test1]=j  #assign time segments to obs
  }
}



## Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:length(cp)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = round(cp[[i]]$breaks$middle, 0),
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}



#Compare elapsed time
time<- map_dfr(elapsed.time, `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))



### Export Results ###

names(all.brkpts)<- names(d.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')


# write.csv(all.brkpts, "data/BCPA_allbrkpts_parametric.csv", row.names = F)
# write.csv(time, "data/BCPA_elapsed_time_parametric.csv", row.names = F)  #units = min


