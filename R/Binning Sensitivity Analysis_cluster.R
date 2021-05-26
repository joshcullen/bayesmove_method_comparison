
##################################
#### Run Sensitivity Analysis ####
##################################

library(bayesmove)
library(tidyverse)
library(lubridate)
library(future)
library(furrr)
library(progressr)

source('R/helper functions.R')



####################################
### Assess 5 bins of equal width ###
####################################

set.seed(1)

#get data
dat<- read.csv("data/Sensitivity_results_5bins_equal.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=5000
nburn=ngibbs/2
nmaxclust=7
ndata.types=length(nbins)

#prior
gamma1=0.1
alpha=0.1

#run Gibbs sampler
obs.list<- df_to_list(dat = obs, ind = "id")
res<- list()
elapsed.time<- vector()

for (i in 1:length(obs.list)) {
  start.time<- Sys.time()
  res[[i]] = cluster_segments(dat=obs.list[[i]], gamma1=gamma1, alpha=alpha,
                              ngibbs=ngibbs, nmaxclust=nmaxclust,
                              nburn=nburn, ndata.types=ndata.types)
  end.time<- Sys.time()
  elapsed.time[i]<- difftime(end.time, start.time, units = "min")
}



#Check traceplot of log likelihood
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#Determine proportion of behaviors (across all time segments)
#top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  # 10/10 > 0.90


## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 2,])  #only select the top 2 behaviors
names(behav.res)<- names(dat.list)

#Plot histograms of proportion data; order color scale from slow to fast
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    ggplot(behav.res[[i]], aes(x = bin, y = prop, fill = as.factor(behav))) +
      geom_bar(stat = 'identity') +
      labs(x = "\nBin", y = "Proportion\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
            axis.text.x.bottom = element_text(size = 12),
            strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
      scale_fill_viridis_d(guide = F) +
      facet_grid(behav ~ var, scales = "free_x")
  )
}
par(ask=F)

### Only 2 states distinguishable for all tracks






#####################################
### Assess 10 bins of equal width ###
#####################################

set.seed(1)

#get data
dat<- read.csv("data/Sensitivity_results_10bins_equal.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(10,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=5000
nburn=ngibbs/2
nmaxclust=7  
ndata.types=length(nbins)

#prior
gamma1=0.1
alpha=0.1

#run Gibbs sampler
obs.list<- df_to_list(dat = obs, ind = "id")
res<- list()
elapsed.time<- vector()

for (i in 1:length(obs.list)) {
  start.time<- Sys.time()
  res[[i]] = cluster_segments(dat=obs.list[[i]], gamma1=gamma1, alpha=alpha,
                              ngibbs=ngibbs, nmaxclust=nmaxclust,
                              nburn=nburn, ndata.types=ndata.types)
  end.time<- Sys.time()
  elapsed.time[i]<- difftime(end.time, start.time, units = "min")
}



#Check traceplot of log likelihood
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#Determine proportion of behaviors (across all time segments)
#top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  # 9/10 > 0.90
#top 3 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3]))  # 10/10 > 0.90


## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)

#Plot histograms of proportion data; order color scale from slow to fast
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    ggplot(behav.res[[i]], aes(x = bin, y = prop, fill = as.factor(behav))) +
      geom_bar(stat = 'identity') +
      labs(x = "\nBin", y = "Proportion\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
            axis.text.x.bottom = element_text(size = 12),
            strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
      scale_fill_viridis_d(guide = F) +
      facet_grid(behav ~ var, scales = "free_x")
  )
}
par(ask=F)

### 3 states distinguishable for 3/5 10k tracks; otherwise 2 states










######################################
### Assess 10 bins using quantiles ###
######################################

set.seed(1)

#get data
dat<- read.csv("data/Sensitivity_results_10bins_quantile.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(10,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=5000
nburn=ngibbs/2
nmaxclust=7  
ndata.types=length(nbins)

#prior
gamma1=0.1
alpha=0.1

#run Gibbs sampler
obs.list<- df_to_list(dat = obs, ind = "id")
res<- list()
elapsed.time<- vector()

for (i in 1:length(obs.list)) {
  start.time<- Sys.time()
  res[[i]] = cluster_segments(dat=obs.list[[i]], gamma1=gamma1, alpha=alpha,
                              ngibbs=ngibbs, nmaxclust=nmaxclust,
                              nburn=nburn, ndata.types=ndata.types)
  end.time<- Sys.time()
  elapsed.time[i]<- difftime(end.time, start.time, units = "min")
}



#Check traceplot of log likelihood
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#Determine proportion of behaviors (across all time segments)
#top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  # 0/10 > 0.90
#top 3 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3]))  # 10/10 > 0.90


## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)

#Plot histograms of proportion data; order color scale from slow to fast
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    ggplot(behav.res[[i]], aes(x = bin, y = prop, fill = as.factor(behav))) +
      geom_bar(stat = 'identity') +
      labs(x = "\nBin", y = "Proportion\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
            axis.text.x.bottom = element_text(size = 12),
            strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
      scale_fill_viridis_d(guide = F) +
      facet_grid(behav ~ var, scales = "free_x")
  )
}
par(ask=F)

### 3 states distinguishable for all tracks