#######################
#### Run LDA Model ####
#######################

set.seed(123)

library(bayesmove)
library(tidyverse)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(ggnewscale)
library(circular)
library(gridExtra)


source('R/helper functions.R')


######################################################
### Analyze simulations with unusual distributions ###
######################################################

#get data
dat<- read.csv("data/CRW_MM_tsegs_weird.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
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
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]),
       xlab = "Iteration", ylab = "Log Likelihood")
  abline(v = ngibbs/2, col = "red", lwd = 1.5)
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#export theta.estim if need to re-analyze later
theta.estim_export<- theta.estim %>%
  map(., as.data.frame) %>%
  bind_rows(., .id = 'id')
# write.csv(theta.estim_export, "data/theta_estim_weird.csv", row.names = F)

#read-in data
# theta.estim<- read.csv("data/theta_estim_weird.csv", as.is=T)
# theta.estim<- theta.estim %>% 
#   split(., .$id) %>% 
#   map(., dplyr::select, -id)
# theta.estim<- theta.estim[c(seq(1, 20, by = 5),  #reorder to match original list
#                             seq(2, 20, by = 5),
#                             seq(3, 20, by = 5),
#                             seq(4, 20, by = 5))]



#boxplots
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
  boxplot(theta.estim[[i]], xlab="Behavior", ylab="Proportion of Total Behavior",
          main = paste("ID",names(theta.estim)[i]))
}
par(mfrow=c(1,1), ask=F)


#Determine proportion of behaviors (across all time segments)
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  #for top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3]))  #for top 3 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>%
  purrr::map(., ~sum(.[1:4]))  #for top 4 states
  

## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)
behav.res_exp<- bind_rows(behav.res, .id = "id")

#export behav.res values
# write.csv(behav.res_exp, "data/CRW MM LDA Phi values_weird.csv", row.names = F)


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





#Need to manually inspect all histograms and assign proper order (slowest to fastest)
behav.order<- list(c(1,3,2), c(1,3,2), c(2,3,1), c(2,1,3), c(1,3,2),
                   c(2,1,3), c(1,2,3), c(2,1,3), c(1,2,3), c(2,1,3),
                   c(1,2,3), c(2,1,3), c(1,2,3), c(1,3,2), c(1,3,2),
                   c(2,3,1), c(3,1,2), c(3,1,2), c(3,2,1), c(2,3,1))


names(behav.order)<- names(theta.estim)

behav.order_exp<- bind_rows(behav.order) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(id = names(behav.order))
# write.csv(behav.order_exp, "data/CRW MM LDA behavior order_weird.csv", row.names = F)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim.long<- list()
for (i in 1:length(dat.list)) {
  theta.estim.long[[i]]<- dat.list[[i]] %>% 
    drop_na(behav_fine) %>% 
    mutate(date=time1-1) %>% 
    expand_behavior(dat = .,
                    theta.estim = theta.estim[[i]],
                    obs = obs.list[[i]],
                    nbehav = 3,
                    behav.names = c("X1","X2","X3"),
                    behav.order = behav.order[[i]]) %>% 
    mutate_at("behavior", as.factor) %>% 
    mutate_at("behavior", ~recode(., 'X1' = behav.order[[i]][1],
                                  'X2' = behav.order[[i]][2], 'X3' = behav.order[[i]][3]))
}








##True proportions by simulation ID
true.behavior.long<- list()
for (i in 1:length(dat.list)) {
  true.behavior.long[[i]]<- data.frame(true.tseg = rep(1:(dat.list[[i]]$track_length[1]/100),
                                                       each = 300),
                                       behav.coarse = rep(dat.list[[i]]$behav_coarse[-1],
                                                          each = 3),
                                       behav_fine = rep(dat.list[[i]]$behav_fine[-1],
                                                        each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(dat.list[[i]]$track_length[1]),each = 3))
  
  true.behavior.long[[i]]$prop<- 0.1
  
  cond<- true.behavior.long[[i]][,"behav.coarse"]
  ind<- which(true.behavior.long[[i]][,"behavior"] == cond)
  
  true.behavior.long[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long[[i]][cond2, "prop"]<- true.behavior.long[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long)<- names(dat.list)




## Plot traces of true and modeled behavior proportions across time segments
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data = true.behavior.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = theta.estim.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)



#assign behavior from sim to data
dat2<- list()
for (i in 1:length(theta.estim.long)) {  #assign behaviors to all obs
  dat2[[i]]<- assign_behavior(dat.orig = dat %>% 
                                filter(id == unique(dat$id)[i]),
                              dat.seg.list = map(dat.list[i], ~slice(., -1)),
                              theta.estim.long=theta.estim.long[[i]],
                              behav.names = c("1","2","3"))
}

dat2<- bind_rows(dat2)


## Plot tracks with modeled behaviors

ggplot(data = dat2 %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




# Elapsed time
time<- elapsed.time %>% 
  unlist() %>% 
  data.frame(time = .)

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))




#export results
# write.csv(dat2, "data/Modeled MM Sim Tracks w Behav_weird.csv", row.names = F)
# write.csv(time, "data/LDA_elapsed_time_weird.csv", row.names = F)  #units = min






##################
#### Figure 5 ####
##################

# Load breakpoints
bayes.brkpts<- read.csv("data/Bayesian_allbreakpts_weird.csv")
true.brkpts_weird<- read.csv("data/CRW_MM_sim_brkpts_weird.csv")


## Part a: segmentation heatmap
data<- dat.list$`2_3`  #ID 2_3

# Reformat data into long form using preferred time variable
dat.long<- data[,c("id", "time1", "dist", "rel.angle")] %>%
  tidyr::pivot_longer(cols = -c(1:2), names_to = "var", values_to = "value")

# Store number of variables/data streams
var.len<- 2
var_names<- c("dist","rel.angle")
var_labels<- c("Step Length (units)", "Turning Angle (rad)")

#Relabel variables is var_labels specified
if (!is.null(var_labels)) {
  for (i in 1:var.len) {
    dat.long$var<- gsub(x = dat.long$var, pattern = var_names[i], replacement = var_labels[i])
  }
}


breakpt<- bayes.brkpts %>% 
  filter(id == unique(data$id) & type == "Model") %>% 
  dplyr::select(brks) %>% 
  rename(breaks = brks)

true.brks<- bayes.brkpts %>% 
  filter(id == unique(data$id) & type == "True") %>% 
  dplyr::select(brks) %>% 
  rename(breaks = brks)

ticks.bottom<- data.frame(x=rep(true.brks$breaks, 2),
                          y=rep(c(-5,-3.5), each = nrow(true.brks)),
                          xend=rep(true.brks$breaks, 2),
                          yend=rep(c(10,-2.5), each = nrow(true.brks)),
                          var=rep(c("Step Length (units)","Turning Angle (rad)"),
                                  each = nrow(true.brks)))
ticks.top<- data.frame(x=rep(true.brks$breaks, 2),
                       y=rep(c(65,2.5), each = nrow(true.brks)),
                       xend=rep(true.brks$breaks, 2),
                       yend=rep(c(80,3.5), each = nrow(true.brks)),
                       var=rep(c("Step Length (units)","Turning Angle (rad)"),
                               each = nrow(true.brks)))





p.seg<- ggplot(data = dat.long, aes(x=time1, y=value, color=var)) +
  geom_line(na.rm = TRUE, size = 0.25) +
  facet_wrap(~var, scales = 'free', nrow = var.len, strip.position = "left") +
  scale_color_brewer("", palette = "Dark2", guide = FALSE) +
  geom_vline(data = breakpt, aes(xintercept = breaks - 0.5),
             color = "black", size = 1.5, alpha = 1) +
  geom_segment(data = ticks.bottom, aes(x = x, xend = xend, y = y, yend = yend),
               color = "purple", size = 1, alpha = 1) +
  geom_segment(data = ticks.top, aes(x = x, xend = xend, y = y, yend = yend),
               color = "purple", size = 1, alpha = 1) +
  labs(x = "\nTime", y = "") +
  scale_y_continuous(expand = c(0.005,0.005)) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.title = element_text(size = 20),
        panel.grid.minor = element_blank())





## Part b: behavior histogram

#define bin number and limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

dist.bin.lims=quantile(dat.list[[12]][dat.list[[12]]$dt == 3600,]$dist,
                       c(0,0.25,0.50,0.75,0.90,1), na.rm=T) #5 bins


#calculate true proportions of SL and TA by behavior for all bins for ID 2_2
SL.params_weird<- data.frame(par1 = c(0.25, 2, exp(2)), par2 = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0), par2 = c(0.5, pi, 1))

true.b_weird<- extract.behav.props_weird(params = list(SL.params_weird, TA.params),
                                         lims = list(dist.bin.lims, angle.bin.lims),
                                         behav.names = c("Encamped","ARS","Transit"))

bayes.b<- behav.res[[12]]
bayes.b$behav<- bayes.b$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "ARS") %>% 
  str_replace_all(., "2", "Encamped") %>%
  str_replace_all(., "3", "Transit") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.hist<- ggplot(true.b_weird, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = bayes.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2, fill="white",
             color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank()) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")




##Part c: generate time series plots comparing behavior proportions
dat1<- theta.estim.long[[12]]
dat1$behavior<- dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

bayes.list<- df_to_list(dat2, "id")
true.behavior.long<- list()
for (i in 1:length(bayes.list)) {
  true.behavior.long[[i]]<- data.frame(true.tseg = rep(1:(bayes.list[[i]]$track_length[1]/100),
                                                       each = 300),
                                       behav_coarse = rep(bayes.list[[i]]$behav_coarse,
                                                          each = 3),
                                       behav_fine = rep(bayes.list[[i]]$behav_fine,
                                                        each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(bayes.list[[i]]$track_length[1]),each = 3))
  
  true.behavior.long[[i]]$prop<- 0.1
  
  cond<- true.behavior.long[[i]][,"behav_coarse"]
  ind<- which(true.behavior.long[[i]][,"behavior"] == cond)
  
  true.behavior.long[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long[[i]][cond2, "prop"]<- true.behavior.long[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long)<- names(bayes.list)

true.dat1<- true.behavior.long[[12]]
true.dat1$behavior<- true.dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

p.prop<- ggplot() +
  geom_path(data = true.dat1,
            aes(x=time1, y=prop, color = behavior),
            size = 1.5, linejoin = "round", lineend = "round") +
  scale_color_manual(values = c(viridis(n=20)[c(1,9)], "gold3"), guide=F) +
  new_scale_color() +
  geom_path(data = dat1,
            aes(x=time1-1, y=prop, color = behavior),
            size = 0.75, linejoin = "round", lineend = "round") +
  scale_color_manual(values = c(viridis(n=20)[c(7,13)], "gold2"), guide=F) +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  facet_wrap(~behavior, nrow = 3)




## Make composite

# png("Figure 5 (results from sim).png", width = 14, height = 10, units = "in", res = 330)
grid.arrange(p.seg, p.hist, p.prop, heights = c(0.2, 1, 0.1, 1),
             layout_matrix = rbind(c(NA, NA),
                                   c(1, 2),
                                   c(NA, NA),
                                   c(3, 3)))
# dev.off()














#########################################################
### Analyze simulations with parametric distributions ###
#########################################################

set.seed(123)

#get data
dat<- read.csv("data/CRW_MM_tsegs_parametric.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
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
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]),
       xlab = "Iteration", ylab = "Log Likelihood")
  abline(v = ngibbs/2, col = "red", lwd = 1.5)
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#export theta.estim if need to re-analyze later
theta.estim_export<- theta.estim %>%
  map(., as.data.frame) %>%
  bind_rows(., .id = 'id')
# write.csv(theta.estim_export, "data/theta_estim_parametric.csv", row.names = F)

#read-in data
# theta.estim<- read.csv("data/theta_estim_parametric.csv", as.is=T)
# theta.estim<- theta.estim %>% 
#   split(., .$id) %>% 
#   map(., dplyr::select, -id)
# theta.estim<- theta.estim[c(seq(1, 20, by = 5),  #reorder to match original list
#                             seq(2, 20, by = 5),
#                             seq(3, 20, by = 5),
#                             seq(4, 20, by = 5))]



#boxplots
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
  boxplot(theta.estim[[i]], xlab="Behavior", ylab="Proportion of Total Behavior",
          main = paste("ID",names(theta.estim)[i]))
}
par(mfrow=c(1,1), ask=F)


#Determine proportion of behaviors (across all time segments)
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  #for top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3]))  #for top 3 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>%
  purrr::map(., ~sum(.[1:4]))  #for top 4 states


## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)
behav.res_exp<- bind_rows(behav.res, .id = "id")

#export behav.res values
# write.csv(behav.res_exp, "data/CRW MM LDA Phi values_parametric.csv", row.names = F)


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





#Need to manually inspect all histograms and assign proper order (slowest to fastest)
behav.order<- list(c(1,3,2), c(1,2,3), c(3,1,2), c(2,1,3), c(3,2,1),
                   c(2,1,3), c(2,1,3), c(2,3,1), c(1,3,2), c(2,3,1),
                   c(2,1,3), c(2,1,3), c(1,2,3), c(1,3,2), c(2,1,3),
                   c(2,1,3), c(1,2,3), c(2,1,3), c(2,3,1), c(3,1,2))


names(behav.order)<- names(theta.estim)

behav.order_exp<- bind_rows(behav.order) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(id = names(behav.order))
# write.csv(behav.order_exp, "data/CRW MM LDA behavior order_parametric.csv", row.names = F)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim.long<- list()
for (i in 1:length(dat.list)) {
  theta.estim.long[[i]]<- dat.list[[i]] %>% 
    drop_na(behav_fine) %>% 
    mutate(date=time1-1) %>% 
    expand_behavior(dat = .,
                    theta.estim = theta.estim[[i]],
                    obs = obs.list[[i]],
                    nbehav = 3,
                    behav.names = c("X1","X2","X3"),
                    behav.order = behav.order[[i]]) %>% 
    mutate_at("behavior", as.factor) %>% 
    mutate_at("behavior", ~recode(., 'X1' = behav.order[[i]][1],
                                  'X2' = behav.order[[i]][2], 'X3' = behav.order[[i]][3]))
}








##True proportions by simulation ID
true.behavior.long<- list()
for (i in 1:length(dat.list)) {
  true.behavior.long[[i]]<- data.frame(true.tseg = rep(1:(dat.list[[i]]$track_length[1]/100),
                                                       each = 300),
                                       behav.coarse = rep(dat.list[[i]]$behav_coarse[-1],
                                                          each = 3),
                                       behav_fine = rep(dat.list[[i]]$behav_fine[-1],
                                                        each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(dat.list[[i]]$track_length[1]),each = 3))
  
  true.behavior.long[[i]]$prop<- 0.1
  
  cond<- true.behavior.long[[i]][,"behav.coarse"]
  ind<- which(true.behavior.long[[i]][,"behavior"] == cond)
  
  true.behavior.long[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long[[i]][cond2, "prop"]<- true.behavior.long[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long)<- names(dat.list)




## Plot traces of true and modeled behavior proportions across time segments
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data = true.behavior.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = theta.estim.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)



#assign behavior from sim to data
dat2<- list()
for (i in 1:length(theta.estim.long)) {  #assign behaviors to all obs
  dat2[[i]]<- assign_behavior(dat.orig = dat %>% 
                                filter(id == unique(dat$id)[i]),
                              dat.seg.list = map(dat.list[i], ~slice(., -1)),
                              theta.estim.long=theta.estim.long[[i]],
                              behav.names = c("1","2","3"))
}

dat2<- bind_rows(dat2)


## Plot tracks with modeled behaviors

ggplot(data = dat2 %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




# Elapsed time
time<- elapsed.time %>% 
  unlist() %>% 
  data.frame(time = .)

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))




#export results
# write.csv(dat2, "data/Modeled MM Sim Tracks w Behav_parametric.csv", row.names = F)
# write.csv(time, "data/LDA_elapsed_time_parametric.csv", row.names = F)  #units = min








#########################################################
### Analyze HMM simulations with common distributions ###
#########################################################

set.seed(123)

#get data
dat<- read.csv("data/HMM_tsegs_parametric.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=10000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
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
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]),
       xlab = "Iteration", ylab = "Log Likelihood")
  abline(v = ngibbs/2, col = "red", lwd = 1.5)
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#export theta.estim if need to re-analyze later
theta.estim_export<- theta.estim %>%
  map(., as.data.frame) %>%
  bind_rows(., .id = 'id')
# write.csv(theta.estim_export, "data/theta(HMM)_estim_parametric.csv", row.names = F)

#read-in data
# theta.estim<- read.csv("data/theta(HMM)_estim_parametric.csv", as.is=T)
# theta.estim<- theta.estim %>% 
#   split(., .$id) %>% 
#   map(., dplyr::select, -id)
# theta.estim<- theta.estim[c(1,3,5,7,9,2,4,6,8,10)]



#boxplots
par(mfrow=c(2,2), ask = T)
for (i in 1:length(theta.estim)) {
  boxplot(theta.estim[[i]], xlab="Behavior", ylab="Proportion of Total Behavior",
          main = paste("ID",names(theta.estim)[i]))
}
par(mfrow=c(1,1), ask=F)


#Determine proportion of behaviors (across all time segments)
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  #for top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3]))  #for top 3 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>%
  purrr::map(., ~sum(.[1:4]))  #for top 4 states


## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)
behav.res_exp<- bind_rows(behav.res, .id = "id")

#export behav.res values
# write.csv(behav.res_exp, "data/HMM LDA Phi values_parametric.csv", row.names = F)


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





#Need to manually inspect all histograms and assign proper order (slowest to fastest)
behav.order<- list(c(1,3,2), c(3,1,2), c(3,1,2), c(2,3,1), c(2,3,1),
                   c(3,2,1), c(3,1,2), c(2,3,1), c(2,1,3), c(1,3,2))


names(behav.order)<- names(theta.estim)

behav.order_exp<- bind_rows(behav.order) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(id = names(behav.order))
# write.csv(behav.order_exp, "data/HMM LDA behavior order_parametric.csv", row.names = F)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim.long<- list()
for (i in 1:length(dat.list)) {
  theta.estim.long[[i]]<- dat.list[[i]] %>% 
    drop_na(state) %>% 
    mutate(date=time1-1) %>% 
    expand_behavior(dat = .,
                    theta.estim = theta.estim[[i]],
                    obs = obs.list[[i]],
                    nbehav = 3,
                    behav.names = c("X1","X2","X3"),
                    behav.order = behav.order[[i]]) %>% 
    mutate_at("behavior", as.factor) %>% 
    mutate_at("behavior", ~recode(., 'X1' = behav.order[[i]][1],
                                  'X2' = behav.order[[i]][2], 'X3' = behav.order[[i]][3]))
}
names(theta.estim.long)<- names(dat.list)




##True proportions by simulation ID

# Calculate "true" breakpoints as observations where state changes
true.brks<- list()
for (i in 1:length(dat.list)) {
  true.brks[[i]]<- which(diff(dat.list[[i]]$state) != 0)
}
names(true.brks)<- names(dat.list)

max.length<- max(sapply(true.brks, length))
true.brks2<- lapply(true.brks, function(x) { c(x, rep(NA, max.length-length(x)))}) %>%
  dplyr::bind_rows() %>%
  t() %>%
  data.frame()

true.brks2<- cbind(id = names(dat.list), true.brks2)
     

# assign true track segments to all observations                  
dat.list2<- assign_tseg(dat.list, true.brks2)


# calculate the true proportion of each state exhibited per time segment (for comparison)
theta.estim.true<- list()
for (j in 1:n_distinct(dat.list2$id)) {
  ind1<- dat.list2 %>%
    filter(id == unique(.$id)[j]) %>% 
    drop_na(state) %>%
    mutate(across(c("tseg", "state"), factor)) %>% 
    group_by(tseg, state, .drop = FALSE) %>%
    count() %>%
    group_by(tseg) %>%
    mutate(id = unique(dat.list2$id)[j],
           prop = n/sum(n)) %>% 
    ungroup()
  
  nobs<- ind1 %>% 
    filter(prop == 1) %>% 
    dplyr::pull(n)
  
  
  for (i in 1:length(nobs)) {
    
    if (i == 1) {
      ind2<- ind1 %>% 
        slice(rep(which(.$tseg == i), each = nobs[i]))
    } else {
      tmp<- ind1 %>% 
        slice(rep(which(.$tseg == i), each = nobs[i]))
      
      ind2<- rbind(ind2, tmp)
    }
  }
  
  
  ind2<- ind2[order(ind2$state),]
  ind2$time1<- rep(2:length(which(dat.list2$id == unique(dat.list2$id)[j])), 3)
  
  theta.estim.true[[j]]<- ind2 
}




## Plot traces of true and modeled behavior proportions across time segments
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data = theta.estim.true[[i]] %>% 
                  drop_na(state) %>% 
                  rename(behavior = state),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = theta.estim.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)



#assign behavior from sim to data
dat2<- list()
for (i in 1:length(theta.estim.long)) {  #assign behaviors to all obs
  dat2[[i]]<- assign_behavior(dat.orig = dat %>% 
                                filter(id == unique(dat$id)[i]),
                              dat.seg.list = map(dat.list[i], ~slice(., -1)),
                              theta.estim.long=theta.estim.long[[i]],
                              behav.names = c("1","2","3"))
}

dat2<- bind_rows(dat2)


## Plot tracks with modeled behaviors

ggplot(data = dat2 %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




# Elapsed time
time<- elapsed.time %>% 
  unlist() %>% 
  data.frame(time = .)

time$track_length<- rep(c('1k','5k'), each = 5) %>% 
  factor(., levels = c('1k','5k'))




#export results
# write.csv(dat2, "data/Modeled HMM Sim Tracks w Behav_parametric.csv", row.names = F)
# write.csv(time, "data/LDA(HMM)_elapsed_time_parametric.csv", row.names = F)  #units = min











#########################################################
### Analyze HMM simulations with uncommon distributions ###
#########################################################

set.seed(123)

#get data
dat<- read.csv("data/HMM_tsegs_weird.csv", as.is = T)
dat.list<- df_to_list(dat = dat, ind = "id")  #for later behavioral assignment
nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% 
  dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- summarize_tsegs(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=10000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
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
  plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]),
       xlab = "Iteration", ylab = "Log Likelihood")
  abline(v = ngibbs/2, col = "red", lwd = 1.5)
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.estim<- map(res, extract_prop, ngibbs = ngibbs, nburn = ngibbs/2, nmaxclust = nmaxclust)
names(theta.estim)<- names(dat.list)

#export theta.estim if need to re-analyze later
theta.estim_export<- theta.estim %>%
  map(., as.data.frame) %>%
  bind_rows(., .id = 'id')
# write.csv(theta.estim_export, "data/theta(HMM)_estim_weird.csv", row.names = F)

#read-in data
# theta.estim<- read.csv("data/theta(HMM)_estim_weird.csv", as.is=T)
# theta.estim<- theta.estim %>% 
#   split(., .$id) %>% 
#   map(., dplyr::select, -id)
# theta.estim<- theta.estim[c(1,3,5,7,9,2,4,6,8,10)]



#boxplots
par(mfrow=c(2,2), ask = T)
for (i in 1:length(theta.estim)) {
  boxplot(theta.estim[[i]], xlab="Behavior", ylab="Proportion of Total Behavior",
          main = paste("ID",names(theta.estim)[i]))
}
par(mfrow=c(1,1), ask=F)


#Determine proportion of behaviors (across all time segments)
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:2]))  #for top 2 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3]))  #for top 3 states
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>%
  purrr::map(., ~sum(.[1:4]))  #for top 4 states


## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)
behav.res_exp<- bind_rows(behav.res, .id = "id")

#export behav.res values
# write.csv(behav.res_exp, "data/HMM LDA Phi values_weird.csv", row.names = F)


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





#Need to manually inspect all histograms and assign proper order (slowest to fastest)
behav.order<- list(c(1,3,2), c(1,3,2), c(3,1,2), c(1,3,2), c(3,1,2),
                   c(3,2,1), c(2,3,1), c(2,3,1), c(2,1,3), c(2,3,1))


names(behav.order)<- names(theta.estim)

behav.order_exp<- bind_rows(behav.order) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(id = names(behav.order))
# write.csv(behav.order_exp, "data/HMM LDA behavior order_weird.csv", row.names = F)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim.long<- list()
for (i in 1:length(dat.list)) {
  theta.estim.long[[i]]<- dat.list[[i]] %>% 
    drop_na(state) %>% 
    mutate(date=time1-1) %>% 
    expand_behavior(dat = .,
                    theta.estim = theta.estim[[i]],
                    obs = obs.list[[i]],
                    nbehav = 3,
                    behav.names = c("X1","X2","X3"),
                    behav.order = behav.order[[i]]) %>% 
    mutate_at("behavior", as.factor) %>% 
    mutate_at("behavior", ~recode(., 'X1' = behav.order[[i]][1],
                                  'X2' = behav.order[[i]][2], 'X3' = behav.order[[i]][3]))
}
names(theta.estim.long)<- names(dat.list)




##True proportions by simulation ID

# Calculate "true" breakpoints as observations where state changes
true.brks<- list()
for (i in 1:length(dat.list)) {
  true.brks[[i]]<- which(diff(dat.list[[i]]$state) != 0)
}
names(true.brks)<- names(dat.list)

max.length<- max(sapply(true.brks, length))
true.brks2<- lapply(true.brks, function(x) { c(x, rep(NA, max.length-length(x)))}) %>%
  dplyr::bind_rows() %>%
  t() %>%
  data.frame()

true.brks2<- cbind(id = names(dat.list), true.brks2)


# assign true track segments to all observations                  
dat.list2<- assign_tseg(dat.list, true.brks2)


# calculate the true proportion of each state exhibited per time segment (for comparison)
theta.estim.true<- list()
for (j in 1:n_distinct(dat.list2$id)) {
  ind1<- dat.list2 %>%
    filter(id == unique(.$id)[j]) %>% 
    drop_na(state) %>%
    mutate(across(c("tseg", "state"), factor)) %>% 
    group_by(tseg, state, .drop = FALSE) %>%
    count() %>%
    group_by(tseg) %>%
    mutate(id = unique(dat.list2$id)[j],
           prop = n/sum(n)) %>% 
    ungroup()
  
  nobs<- ind1 %>% 
    filter(prop == 1) %>% 
    dplyr::pull(n)
  
  
  for (i in 1:length(nobs)) {
    
    if (i == 1) {
      ind2<- ind1 %>% 
        slice(rep(which(.$tseg == i), each = nobs[i]))
    } else {
      tmp<- ind1 %>% 
        slice(rep(which(.$tseg == i), each = nobs[i]))
      
      ind2<- rbind(ind2, tmp)
    }
  }
  
  
  ind2<- ind2[order(ind2$state),]
  ind2$time1<- rep(2:length(which(dat.list2$id == unique(dat.list2$id)[j])), 3)
  
  theta.estim.true[[j]]<- ind2 
}




## Plot traces of true and modeled behavior proportions across time segments
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data = theta.estim.true[[i]] %>% 
                  drop_na(state) %>% 
                  rename(behavior = state),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = theta.estim.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)



#assign behavior from sim to data
dat2<- list()
for (i in 1:length(theta.estim.long)) {  #assign behaviors to all obs
  dat2[[i]]<- assign_behavior(dat.orig = dat %>% 
                                filter(id == unique(dat$id)[i]),
                              dat.seg.list = map(dat.list[i], ~slice(., -1)),
                              theta.estim.long=theta.estim.long[[i]],
                              behav.names = c("1","2","3"))
}

dat2<- bind_rows(dat2)


## Plot tracks with modeled behaviors

ggplot(data = dat2 %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




# Elapsed time
time<- elapsed.time %>% 
  unlist() %>% 
  data.frame(time = .)

time$track_length<- rep(c('1k','5k'), each = 5) %>% 
  factor(., levels = c('1k','5k'))




#export results
# write.csv(dat2, "data/Modeled HMM Sim Tracks w Behav_weird.csv", row.names = F)
# write.csv(time, "data/LDA(HMM)_elapsed_time_weird.csv", row.names = F)  #units = min
