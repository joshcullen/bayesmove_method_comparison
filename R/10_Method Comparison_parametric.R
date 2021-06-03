#########################
### Method Comparison ###
#########################

library(tidyverse)
library(wesanderson)
library(lubridate)
library(cowplot)
library(viridis)
library(ggnewscale)

source('R/helper functions.R')


# Load elapsed time
seg.time<- read.csv("data/Bayesian_elapsed_time_parametric.csv")
lda.time<- read.csv("data/LDA_elapsed_time_parametric.csv")
bcpa.time<- read.csv("data/BCPA_elapsed_time_parametric.csv")
hmm.time<- read.csv("data/HMM_elapsed_time_parametric.csv")
segclust.time<- read.csv("data/Segclust2d_elapsed_time_parametric.csv")
embc.time<- read.csv("data/EMbC_elapsed_time_parametric.csv")

# Load breakpoints
bayes.brkpts<- read.csv("data/Bayesian_allbreakpts_parametric.csv")
bcpa.brkpts<- read.csv("data/BCPA_allbrkpts_parametric.csv")
segclust.brkpts<- read.csv("data/Segclust2d allbrkpts_parametric.csv")

# Load results
bayes.res_parametric<- read.csv("data/Modeled MM Sim Tracks w Behav_parametric.csv")
hmm.res_parametric<- read.csv("data/HMM results_parametric.csv")
segclust.res_parametric<- read.csv("data/Segclust2d results_parametric.csv")
embc.res_parametric<- read.csv("data/EMbC results_parametric.csv")

# Load true breakpoints
true.brkpts_parametric<- read.csv("data/CRW_MM_sim_brkpts_parametric.csv")


############################
### Compare Elapsed Time ###
############################

#Add times together for segmentation and LDA model
bayes.time<- data.frame(time = seg.time$time + lda.time$time,
                        track_length = seg.time$track_length)

time<- rbind(bayes.time, bcpa.time, hmm.time, segclust.time, embc.time)
time$method<- c(rep(c("M4", "BCPA", "HMM"), each = 20),
                rep("Segclust2d", 15),
                rep("EMbC", 20)) %>% 
  factor(., levels = c('M4','HMM','Segclust2d','EMbC','BCPA'))
time$track_length<- time$track_length %>% 
  factor(., levels = c('1k','5k','10k','50k'))

pal1<- c(wes_palette("Darjeeling1", 5)[-4], "mediumorchid")


p.time<- ggplot(time, aes(track_length, time, fill = method, color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "point", shape = 15, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "Elapsed Time (min)") +
  # labs(x = '', y = expression(log[10]~'(Elapsed Time) (min)')) +
  scale_x_discrete(labels = c(1000,5000,10000,50000)) +
  scale_y_log10(breaks = c(10^(-2:2), 200, 400),
                labels = c(0.01, 0.1, 1, 10, 100, 200, 400)) +
  annotation_logticks(sides = "lr") +
  scale_fill_manual("", values = pal1) +
  scale_color_manual("", values = pal1) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = c(0.15,0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 12))


###########################
### Compare Breakpoints ###
###########################

bayes.brkpts$method<- rep("M4", nrow(bayes.brkpts))
bcpa.brkpts$method<- rep("BCPA", nrow(bcpa.brkpts))
segclust.brkpts$method<- rep("Segclust2d", nrow(segclust.brkpts))
all.brkpts<- rbind(bayes.brkpts, bcpa.brkpts, segclust.brkpts)


brkpt.acc<- all.brkpts %>% 
  group_by(method, id, acc) %>% 
  filter(type == "Model") %>% 
  tally() %>% 
  mutate(freq = n/sum(n)) %>% 
  mutate(track_length = case_when(str_detect(id, "_1") ~ "1000",
                                  str_detect(id, "_2") ~ "5000",
                                  str_detect(id, "_3") ~ "10000",
                                  str_detect(id, "_4") ~ "50000")) %>% 
  group_by(method, track_length, id) %>% 
  filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
  summarise(freq = sum(freq)) %>%  #and calculate accuracy across all sims combined
  ungroup()

brkpt.acc$track_length<- brkpt.acc$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))
brkpt.acc$method<- factor(brkpt.acc$method, levels = c('M4','Segclust2d','BCPA'))


#Accuracy (includes 'accurate' and 'accurate duplicate' classifications)
p.brk<- ggplot(brkpt.acc, aes(track_length, freq, fill = method, color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "point", shape = 15, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "Accuracy of Breakpoints\n") +
  scale_fill_manual("", values = pal1[c(1,3,5)], guide = F) +
  scale_color_manual("", values = pal1[c(1,3,5)], guide = F) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  ylim(0,1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank())


#calc mean accuracy per track_length and method
brkpt.acc %>% 
  group_by(method, track_length) %>% 
  summarise(mean=mean(freq))

# number of breakpoints missed
all.brkpts %>% 
  filter(acc == "Missing") %>% 
  group_by(method) %>% 
  tally()

# total number of true breakpoints 
all.brkpts %>% 
  filter(type == "True") %>% 
  group_by(method) %>% 
  tally()






##############################################
### Compare Accuracy of Behavior Estimates ###
##############################################

# Add "blank" first row to embc.res to match nrow of other objects
embc.res_parametric<- embc.res_parametric %>% 
  bayesmove::df_to_list("id") %>% 
  map(., 
      ~rbind(c(unique(.$id), 0, 0, rep(NA, 4), unique(.$track_length), rep(NA, 2)),
             .x)) %>% 
  bind_rows()

# Assign identifiers by method and make consistent behavior colname
bayes.res_parametric$method<- rep("M4", nrow(bayes.res_parametric))
hmm.res_parametric$method<- rep("HMM", nrow(hmm.res_parametric))
segclust.res_parametric$method<- rep("Segclust2d", nrow(segclust.res_parametric))
embc.res_parametric$method<- rep("EMbC", nrow(embc.res_parametric))

bayes.res_parametric<- bayes.res_parametric %>% 
  rename(state = behav) %>% 
  mutate_at("state", ~factor(., levels = 1:3))
hmm.res_parametric<- hmm.res_parametric %>% 
  rename(state = hmm.state, id = ID)
embc.res_parametric.agg<- embc.res_parametric %>%  #combine HL and HH states into one
  mutate_at("state", ~case_when(state == "LH" ~ 1,
                                state == "LL" ~ 2,
                                state == "HL" | state == "HH" ~ 3))


# Modify hmm.res_parametric, segclust.res_parametric, and embc.res_parametric to be same as bayes.res_parametric format
# calc proportions of behaviors by true time segment and then identify dominant behavior
hmm.res_parametric2<- hmm.res_parametric %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_parametric) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  dplyr::select(-c(time2, n)) %>% 
  pivot_wider(names_from = behavior, values_from = prop) %>% 
  mutate(track_length = max(time1)) %>% 
  ungroup() %>% 
  mutate(state = apply(.[,4:6], 1, which.max)) %>% 
  mutate_at("id", as.character) %>% 
  arrange(track_length, id) %>% 
  dplyr::select(-track_length)

hmm.res_parametric3<- hmm.res_parametric %>% 
  dplyr::select(-state) %>% 
  bayesmove::df_to_list("id") %>% 
  map2(.,
       bayesmove::df_to_list(hmm.res_parametric2, "id") %>% 
         map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()




segclust.res_parametric2<- segclust.res_parametric %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_segclust) %>% 
  bind_rows() %>% 
  drop_na("state") %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  dplyr::select(-c(time2, n)) %>% 
  pivot_wider(names_from = behavior, values_from = prop) %>% 
  mutate(track_length = max(time1)) %>% 
  ungroup() %>% 
  mutate(state = apply(.[,4:6], 1, which.max)) %>% 
  mutate_at("id", as.character) %>% 
  arrange(track_length, id) %>%
  dplyr::select(-track_length)


segclust.res_parametric3<- segclust.res_parametric %>% 
  dplyr::select(-state) %>% 
  bayesmove::df_to_list("id") %>% 
  map2(.,
       bayesmove::df_to_list(segclust.res_parametric2, "id") %>% 
         map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()






embc.res_parametric2<- embc.res_parametric.agg %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_parametric) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  dplyr::select(-c(time2, n)) %>% 
  pivot_wider(names_from = behavior, values_from = prop) %>% 
  mutate(track_length = max(time1)) %>% 
  ungroup() %>% 
  mutate(state = apply(.[,4:6], 1, which.max)) %>% 
  mutate_at("id", as.character) %>% 
  arrange(track_length, id) %>% 
  dplyr::select(-track_length)

embc.res_parametric3<- embc.res_parametric.agg %>% 
  dplyr::select(-state) %>% 
  bayesmove::df_to_list("id") %>% 
  map2(.,
       bayesmove::df_to_list(embc.res_parametric2, "id") %>% 
         map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()



# Combine all datasets
res_parametric<- rbind(bayes.res_parametric[,c("id","behav_fine","behav_coarse","track_length",
                                               "state", "method")],
                  hmm.res_parametric3[,c("id","behav_fine","behav_coarse","track_length",
                                         "state","method")],
                  segclust.res_parametric3[,c("id","behav_fine","behav_coarse","track_length",
                                              "state", "method")],
                  embc.res_parametric3[,c("id","behav_fine","behav_coarse","track_length"
                                          ,"state","method")])
res_parametric$method<- factor(res_parametric$method, levels = c('M4','HMM','Segclust2d',
                                                                 'EMbC'))


## Overall

#Coarse-scale behavior
res_parametric %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate_at("track_length", as.numeric) %>%
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))

summ.stats_coarse_parametric<- res_parametric %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate_at("track_length", as.numeric) %>%
  mutate(acc = n/track_length) %>% 
  ungroup()

summ.stats_coarse_parametric$track_length<- summ.stats_coarse_parametric$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))

p.coarse_parametric<- ggplot(summ.stats_coarse_parametric, aes(track_length,acc,fill = method,
                                                     color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "point", shape = 15, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  ylim(0,1) +
  labs(x="\nTrack Length (observations)", y = "Accuracy of Behavior Estimates\n") +
  scale_fill_manual("", values = pal1[c(1:4)]) +
  scale_color_manual("", values = pal1[c(1:4)]) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))




#### Compare Accuracy of Bayesian, HMM, and Segclust2d Proportion Estimates ####

##True proportions by simulation ID
bayes.list_parametric<- bayesmove::df_to_list(bayes.res_parametric, "id")

true.behavior.long_parametric<- list()
for (i in 1:length(bayes.list_parametric)) {
  true.behavior.long_parametric[[i]]<- 
    data.frame(true.tseg = rep(1:(bayes.list_parametric[[i]]$track_length[1]/100), each = 300),
               behav_coarse = rep(bayes.list_parametric[[i]]$behav_coarse[-1],
                                  each = 3),
               behav_fine = rep(bayes.list_parametric[[i]]$behav_fine[-1],
                                each = 3),
               behavior = rep(1:3, 1000),
               time1 = rep(1:(bayes.list_parametric[[i]]$track_length[1]),
                           each = 3))
  
  true.behavior.long_parametric[[i]]$prop<- 0.1
  
  cond<- true.behavior.long_parametric[[i]][,"behav_coarse"]
  ind<- which(true.behavior.long_parametric[[i]][,"behavior"] == cond)
  
  true.behavior.long_parametric[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long_parametric[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long_parametric[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long_parametric[[i]][cond2, "prop"]<- true.behavior.long_parametric[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long_parametric)<- names(bayes.list_parametric)

## True proportions for HMMs (from time segments using true breakpoints)
hmm.props_parametric<- hmm.res_parametric %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_parametric) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  ungroup()

par(ask=T)
for (i in 1:length(unique(as.character(hmm.props_parametric$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=true.behavior.long_parametric[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = hmm.props_parametric %>% 
                  filter(id == unique(as.character(hmm.res_parametric$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(hmm.res_parametric$id))[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)





## True proportions for Segclust2d (from best model w/ 3 states)
segclust.props_parametric<- segclust.res_parametric %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_segclust) %>% 
  bind_rows() %>% 
  drop_na("state") %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  ungroup()

par(ask=T)
for (i in 1:length(unique(as.character(segclust.props_parametric$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=true.behavior.long_parametric[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = segclust.props_parametric %>% 
                  filter(id == unique(as.character(segclust.res_parametric$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(segclust.res_parametric$id))[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)






## True proportions for EMbC (from time segments using true breakpoints)
embc.props_parametric<- embc.res_parametric.agg %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_parametric) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  ungroup()

par(ask=T)
for (i in 1:length(unique(as.character(embc.props_parametric$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=true.behavior.long_parametric[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = embc.props_parametric %>% 
                  filter(id == unique(as.character(embc.res_parametric$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(embc.res_parametric$id))[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)






## True proportions for Bayesian model (from modeled time segments)
bayes.props_parametric<- bayes.res_parametric %>% 
  rename(Encamped = X1, ARS = X2, Transit = X3) %>% 
  drop_na() %>% 
  pivot_longer(., cols = c(Encamped, ARS, Transit), names_to = "behavior",
               values_to = "prop") %>% 
  dplyr::select(id, tseg.y, behavior, prop, time1) %>% 
  mutate_at("behavior", ~recode(., 'Encamped' = 1, 'ARS' = 2, 'Transit' = 3))
bayes.props_parametric$time1<- bayes.props_parametric$time1 - 1


## Calculate RMSE
true.behavior_parametric<- true.behavior.long_parametric %>% 
  bind_rows(.id = "id")

hmm.rmse_parametric<- vector()
for (i in 1:length(unique(hmm.res_parametric$id))) {
  ind<- unique(as.character(hmm.res_parametric$id))[i]
  
  hmm.rmse_parametric[i]<- sqrt(sum((hmm.props_parametric[hmm.props_parametric$id == ind, "prop"] - 
                                  true.behavior_parametric[true.behavior_parametric$id == ind, "prop"])^2) / nrow(hmm.props_parametric[hmm.props_parametric$id == ind,]))
}



segclust.rmse_parametric<- vector()
for (i in 1:length(unique(segclust.res_parametric$id))) {
  ind<- unique(as.character(segclust.res_parametric$id))[i]
  
  segclust.rmse_parametric[i]<- sqrt(sum((segclust.props_parametric[segclust.props_parametric$id == ind, "prop"] - 
                                       true.behavior_parametric[true.behavior_parametric$id == ind, "prop"])^2) / nrow(segclust.props_parametric[segclust.props_parametric$id == ind,]))
}



embc.rmse_parametric<- vector()
for (i in 1:length(unique(embc.res_parametric.agg$id))) {
  ind<- unique(as.character(embc.res_parametric.agg$id))[i]
  
  embc.rmse_parametric[i]<- sqrt(sum((embc.props_parametric[embc.props_parametric$id == ind, "prop"] - 
                                   true.behavior_parametric[true.behavior_parametric$id == ind, "prop"])^2) / nrow(embc.props_parametric[embc.props_parametric$id == ind,]))
}



bayes.rmse_parametric<- vector()
for (i in 1:length(unique(bayes.res_parametric$id))) {
  ind<- unique(as.character(bayes.res_parametric$id))[i]
  
  bayes.rmse_parametric[i]<- sqrt(sum((bayes.props_parametric[bayes.props_parametric$id == ind, "prop"] - 
                                    true.behavior_parametric[true.behavior_parametric$id == ind, "prop"])^2) /
                               nrow(bayes.props_parametric[bayes.props_parametric$id == ind,]))
}


segclust.rmse.df<- data.frame(id = unique(as.character(segclust.res_parametric$id)),
                              track_length = factor(rep(c(1000,5000,10000), each = 5),
                                                    levels = c("1000","5000","10000","50000")),
                              rmse = segclust.rmse_parametric,
                              method = rep("Segclust2d", 15))
rmse.df_parametric<- data.frame(id = rep(unique(as.character(hmm.res_parametric$id)), 3),
                           track_length = factor(rep(rep(c(1000,5000,10000,50000), each = 5), 3),
                                                 levels = c("1000","5000","10000","50000")),
                           rmse = c(bayes.rmse_parametric, hmm.rmse_parametric,
                                    embc.rmse_parametric),
                           method = rep(c("M4","HMM","EMbC"), each = 20))
rmse.df_parametric<- rbind(rmse.df_parametric, segclust.rmse.df)
rmse.df_parametric$method<- factor(rmse.df_parametric$method, levels = c('M4','HMM',
                                                                         'Segclust2d','EMbC'))


#Plot results

p.rmse_parametric<- ggplot(rmse.df_parametric, aes(track_length, rmse, fill = method,
                                                   color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "point", shape = 15, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="\nTrack Length (observations)", y = "RMSE\n") +
  scale_fill_manual("", values = pal1[c(1:4)]) +
  scale_color_manual("", values = pal1[c(1:4)]) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), limits = c(0, 0.65)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))




res.comp<- plot_grid(NULL, NULL, NULL,
          p.time + theme(legend.position = "none"), NULL, p.brk,
          # NULL, NULL, NULL,
          # p.coarse, NULL, p.rmse,
          NULL, NULL, NULL,
          p.coarse_parametric, NULL, p.rmse_parametric,
          align = "hv", nrow = 4, rel_widths = c(1,0.1,1), rel_heights = c(0.2,1,0.1,1))

# extract the legend from one of the plots
legend.comp_res<- get_legend(p.time + theme(legend.position="top",
                                            legend.text = element_text(size = 14)))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(legend.comp_res, res.comp, ncol = 1, rel_heights = c(0.1, 1))

# ggsave("Figure 3 (method comparison)_updatedParametric.png", width = 12, height = 9, units = "in", dpi = 330)









###############################################################################
### Compare Characterization of Step Length and Turning Angle Distributions ###
###############################################################################

library(circular)


## Figure S1 ##

### Define bin limits
dat_parametric<- read.csv("data/CRW_MM_sim_parametric.csv", as.is = T)
dat_parametric$dt<- 3600
names(dat_parametric)[4:5]<- c("dist","rel.angle")

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

dist.bin.lims=quantile(dat_parametric[dat_parametric$dt == 3600,]$dist,
                       c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins

hmm.SL.params_parametric<- read.csv("data/HMM result step params_parametric.csv", as.is = T)  
hmm.TA.params_parametric<- read.csv("data/HMM result angle params_parametric.csv", as.is = T)  

#Manipulate param dfs to reformat for extract.behav.props()
hmm.SL.params2_parametric<- list()
for (i in 1:nrow(hmm.SL.params_parametric)) {
  hmm.SL.params2_parametric[[i]]<- data.frame(par1 = as.numeric(hmm.SL.params_parametric[i,2:4]),
                                         par2 = as.numeric(hmm.SL.params_parametric[i,5:7]))
}
names(hmm.SL.params2_parametric)<- hmm.SL.params_parametric$X

hmm.TA.params2_parametric<- list()
for (i in 1:nrow(hmm.TA.params_parametric)) {
  hmm.TA.params2_parametric[[i]]<- data.frame(par1 = as.numeric(hmm.TA.params_parametric[i,2:4]),
                                         par2 = as.numeric(hmm.TA.params_parametric[i,5:7]))
}
names(hmm.TA.params2_parametric)<- hmm.TA.params_parametric$X



#convert mean and sd to shape and rate params for gamma dist
for (j in 1:length(hmm.SL.params2_parametric)) {
  for (i in 1:nrow(hmm.SL.params2_parametric[[j]])) {
    shape<- (hmm.SL.params2_parametric[[j]][i,1]^2) / (hmm.SL.params2_parametric[[j]][i,2]^2)
    rate<- hmm.SL.params2_parametric[[j]][i,1] / (hmm.SL.params2_parametric[[j]][i,2]^2)
    
    params<- c(shape, rate)
    
    hmm.SL.params2_parametric[[j]][i,]<- params
  }
}



segclust.params_parametric<- read.csv("data/Segclust result state params_parametric.csv", as.is = T) %>%
  dplyr::select(-1)

#Manipulate param dfs to reformat for extract.behav.props()
segclust.params2_parametric<- list()
for (i in 1:length(unique(segclust.params_parametric$id))) {
  tmp<- segclust.params_parametric %>% 
    filter(id == unique(segclust.params_parametric$id)[i])
  segclust.params2_parametric[[i]]<- tmp[,4:7]
  names(segclust.params2_parametric[[i]])<- paste0("par", 1:4)
}
names(segclust.params2_parametric)<- unique(segclust.params_parametric$id)






samp.size<- embc.res_parametric %>%  #for use in merging Gaussians of HL and HH
  mutate_at("state", factor, levels = c("LL","LH","HL","HH")) %>% 
  group_by(id, state) %>% 
  tally() %>% 
  drop_na()
embc.params_parametric<- read.csv("data/EMbC result state params_parametric.csv", as.is = T) %>% 
  dplyr::select(-1) %>% 
  mutate(n = samp.size$n)

#Manipulate param dfs to reformat for extract.behav.props()
embc.params2_parametric<- list()
for (i in 1:length(unique(embc.params_parametric$id))) {
  tmp<- embc.params_parametric %>% 
    filter(id == unique(embc.params_parametric$id)[i])
  
  # Update SL for merged HL and HH states
  SL.mu.hat<- ((tmp$n[3]*tmp$SL.mu[3]) + (tmp$n[4]*tmp$SL.mu[4])) / sum(tmp$n[3:4])
  
  SL.sd.num<- ((tmp$SL.sd[3]^2 + tmp$SL.mu[3]^2)*tmp$n[3]) + ((tmp$SL.sd[4]^2 + 
                                                                 tmp$SL.mu[4]^2)*tmp$n[4])
  SL.sd.den<- sum(tmp$n[3:4])
  SL.sd.hat<- sqrt((SL.sd.num/SL.sd.den) - SL.mu.hat^2)
  
  
  # Update absTA for merged HL and HH states
  absTA.mu.hat<- ((tmp$n[3]*tmp$absTA.mu[3]) + (tmp$n[4]*tmp$absTA.mu[4])) / sum(tmp$n[3:4])
  
  absTA.sd.num<- ((tmp$absTA.sd[3]^2 + tmp$absTA.mu[3]^2)*tmp$n[3]) + ((tmp$absTA.sd[4]^2 + 
                                                                          tmp$absTA.mu[4]^2)*tmp$n[4])
  absTA.sd.den<- sum(tmp$n[3:4])
  absTA.sd.hat<- sqrt((absTA.sd.num/absTA.sd.den) - absTA.mu.hat^2)
  
  
  # Merge all SL and TA params together
  H<- c(SL.mu.hat, absTA.mu.hat, SL.sd.hat, absTA.sd.hat)
  
  embc.params2_parametric[[i]]<- rbind(tmp[2:1, 2:5], H)
  names(embc.params2_parametric[[i]])<- paste0("par", 1:4)
}
names(embc.params2_parametric)<- unique(embc.params_parametric$id)





SL.params<- data.frame(par1 = c(0.25, 2, 10), par2 = c(1, 1, 1))
TA.params<- data.frame(par1 = c(pi, pi, 0), par2 = c(0.8, 0, 0.8))

true.b_parametric<- extract.behav.props(params = list(SL.params, TA.params),
                                         lims = list(dist.bin.lims, angle.bin.lims),
                                         behav.names = c("Encamped","ARS","Transit"))


hmm.b_parametric<- map2(hmm.SL.params2_parametric, hmm.TA.params2_parametric,
                   ~extract.behav.props(params = list(.x, .y),
                                        lims = list(dist.bin.lims, angle.bin.lims),
                                        behav.names = c("Encamped","ARS","Transit"))
)


segclust.b_parametric<- map(segclust.params2_parametric,
                            ~extract.behav.props_norm(params = .x,
                                                      lims = list(dist.bin.lims,angle.bin.lims),
                                                      behav.names = c("Encamped","ARS","Transit"))
)


embc.b_parametric<- map(embc.params2_parametric,
                   ~extract.behav.props_norm(params = .x,
                                                 lims = list(dist.bin.lims, angle.bin.lims),
                                                 behav.names = c("Encamped","ARS","Transit"))
)




## Bayesian
behav.res_parametric<-  read.csv("data/CRW MM LDA Phi values_parametric.csv", as.is = T)
behav.res_parametric<- bayesmove::df_to_list(behav.res_parametric, "id")

behav.order_parametric<- read.csv("data/CRW MM LDA behavior order_parametric.csv", as.is = T)
behav.order_parametric<- bayesmove::df_to_list(behav.order_parametric, "id")
behav.order_parametric<- map(behav.order_parametric, as.numeric)


#calculate true proportions of SL and TA by behavior for all bins for ID 2_2 for Bayesian model
bayes.b_parametric<- behav.res_parametric[[15]]
bayes.b_parametric$behav<- bayes.b_parametric$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "ARS") %>% 
  str_replace_all(., "2", "Encamped") %>%
  str_replace_all(., "3", "Transit") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.distfit<- ggplot(true.b_parametric, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = hmm.b_parametric[[15]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[2], color="black", stroke=1) +
  geom_point(data = segclust.b_parametric[[15]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[3], color="black", stroke=1) +
  geom_point(data = embc.b_parametric[[15]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[4], color="black", stroke=1) +
  geom_point(data = bayes.b_parametric, aes(x=bin, y=prop, group = behav), pch=21, size = 2,
             fill=pal1[1], color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")



### Calculate root mean square error (RMSE)

#Bayesian
bayes.rmse_parametric<- list()
for (i in 1:length(behav.res_parametric)) {
  bayes.b_parametric<- behav.res_parametric[[i]]
  
  behavs<- c("Encamped","ARS","Transit")
  
  bayes.b_parametric$behav<- bayes.b_parametric$behav %>% 
    factor(levels = behav.order_parametric[[i]]) %>% 
    str_replace_all(., "1", behavs[behav.order_parametric[[i]][1]]) %>% 
    str_replace_all(., "2", behavs[behav.order_parametric[[i]][2]]) %>%
    str_replace_all(., "3", behavs[behav.order_parametric[[i]][3]]) %>%
    factor(., levels = c("Encamped","ARS","Transit"))
  bayes.b_parametric<- bayes.b_parametric %>% 
    group_by(var) %>% 
    arrange(behav) %>% 
    ungroup()
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(bayes.b_parametric$var))) {
    tmp<- bayes.b_parametric %>% 
      filter(var == unique(bayes.b_parametric$var)[j])
    true.tmp<- true.b_parametric %>% 
      filter(var == unique(bayes.b_parametric$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  bayes.rmse_parametric[[i]]<- data.frame(vec)
  names(bayes.rmse_parametric)[i]<- names(behav.res_parametric)[i]
}
bayes.rmse_parametric<- bind_rows(bayes.rmse_parametric) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))




#HMM
hmm.rmse_parametric<- list()
for (i in 1:length(behav.res_parametric)) {
  hmm.b_parametric<- extract.behav.props(params = list(hmm.SL.params2_parametric[[i]],
                                                  hmm.TA.params2_parametric[[i]]),
                                    lims = list(dist.bin.lims, angle.bin.lims),
                                    behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(hmm.b_parametric$var))) {
    tmp<- hmm.b_parametric %>% 
      filter(var == unique(hmm.b_parametric$var)[j])
    true.tmp<- true.b_parametric %>% 
      filter(var == unique(hmm.b_parametric$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  hmm.rmse_parametric[[i]]<- data.frame(vec)
  names(hmm.rmse_parametric)[i]<- names(behav.res_parametric)[i]
}
hmm.rmse_parametric<- bind_rows(hmm.rmse_parametric) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))





#Segclust2d
segclust.rmse_parametric<- list()
for (i in 1:length(segclust.params2_parametric)) {
  segclust.b_parametric<- extract.behav.props_norm(params = segclust.params2_parametric[[i]],
                                                  lims = list(dist.bin.lims, angle.bin.lims),
                                                  behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(segclust.b_parametric$var))) {
    tmp<- segclust.b_parametric %>% 
      filter(var == unique(segclust.b_parametric$var)[j])
    true.tmp<- true.b_parametric %>% 
      filter(var == unique(segclust.b_parametric$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  segclust.rmse_parametric[[i]]<- data.frame(vec)
  names(segclust.rmse_parametric)[i]<- names(behav.res_parametric)[i]
}
segclust.rmse_parametric<- bind_rows(segclust.rmse_parametric) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 15))






#EMbC
embc.rmse_parametric<- list()
for (i in 1:length(embc.params2_parametric)) {
  embc.b_parametric<- extract.behav.props_norm(params = embc.params2_parametric[[i]],
                                              lims = list(dist.bin.lims, angle.bin.lims),
                                              behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(embc.b_parametric$var))) {
    tmp<- embc.b_parametric %>% 
      filter(var == unique(embc.b_parametric$var)[j])
    true.tmp<- true.b_parametric %>% 
      filter(var == unique(embc.b_parametric$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  embc.rmse_parametric[[i]]<- data.frame(vec)
  names(embc.rmse_parametric)[i]<- names(behav.res_parametric)[i]
}
embc.rmse_parametric<- bind_rows(embc.rmse_parametric) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))




segclust.rmse.df<- data.frame(id = rep(unique(as.character(segclust.res_parametric$id)), each = 2),
                              track_length = factor(rep(rep(c(1000,5000,10000), each = 5),
                                                        each = 2),
                                                    levels = c("1000","5000","10000","50000")),
                              rmse = segclust.rmse_parametric,
                              method = rep("Segclust2d", 30))
rmse.df_parametric<- data.frame(id = rep(rep(names(behav.res_parametric), 3), each = 2),
                           track_length = factor(rep(rep(rep(c(1000,5000,10000,50000),
                                                             each = 5), 3), each = 2),
                                                 levels = c("1000","5000","10000","50000")),
                           rmse = rbind(bayes.rmse_parametric, hmm.rmse_parametric,
                                        embc.rmse_parametric),
                           method = rep(c("M4","HMM","EMbC"), each = 40))
rmse.df_parametric<- rbind(rmse.df_parametric, segclust.rmse.df)
rmse.df_parametric$method<- factor(rmse.df_parametric$method,
                                   levels = c('M4','HMM','Segclust2d','EMbC'))


p.rmse_parametric<- ggplot(rmse.df_parametric, aes(track_length, rmse.value, fill = method,
                                         color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="\nTrack Length (observations)", y = "RMSE\n") +
  scale_fill_manual("", values = pal1[c(1:4)]) +
  scale_color_manual("", values = pal1[c(1:4)]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = c(0.9,0.88),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold")) +
  facet_wrap(~rmse.var, ncol = 2)




## Make composite
library(gridExtra)
# png("Figure S1 (rmse from sim)_parametric.png", width = 14.5, height = 5.5, units = "in", res = 330)

grid.arrange(p.distfit, p.rmse_parametric, heights = c(0.2, 1),
             widths = c(1, 0.2, 1, 0.5),
             layout_matrix = rbind(c(NA, NA, NA, NA),
                                   c(3, NA, 4, 4)))
# dev.off()










## Figure 5 ##

### Make plot comparing HMM gamma distribs and Segclust2d normal distribs against true distributions
hmm.SL.df_parametric<- hmm.SL.params2_parametric %>% 
  bind_rows(., .id = "id") %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

segclust.SL.df_parametric<- segclust.params2_parametric %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par2, par4)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 15),
         track_length = factor(rep(c(1000,5000,10000), each = 15),
                               levels = c("1000","5000","10000","50000")))

embc.SL.df_parametric<- embc.params2_parametric %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par2, par4)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

true.SL.params_parametric<- data.frame(id = 1:3,
                                  SL.params,
                                  behavior = c("Encamped","ARS","Transit"),
                                  track_length = "True")


plot_data_SL_hmm<- 
  pmap_df(hmm.SL.df_parametric,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dgamma(x, shape = par1, rate = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_SL_hmm$behavior<- factor(plot_data_SL_hmm$behavior, levels = c("Encamped","ARS","Transit"))


plot_data_SL_segclust<- 
  pmap_df(segclust.SL.df_parametric,
          function(id, par1, par3, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dnorm(x, mean = par1, sd = par3),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_SL_segclust$behavior<- factor(plot_data_SL_segclust$behavior,
                                        levels = c("Encamped","ARS","Transit"))


plot_data_SL_embc<- 
  pmap_df(embc.SL.df_parametric,
          function(id, par1, par3, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dnorm(x, mean = par1, sd = par3),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_SL_embc$behavior<- factor(plot_data_SL_embc$behavior,
                                    levels = c("Encamped","ARS","Transit"))


true_plot_data_SL<- 
  pmap_df(true.SL.params_parametric,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dgamma(x, shape = par1, rate = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
true_plot_data_SL$behavior<- factor(true_plot_data_SL$behavior,
                                    levels = c("Encamped","ARS","Transit"))


# Plot overlapping densities per behavior (step lengths)
p.SL_hmm<- ggplot(data = plot_data_SL_hmm, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_data_SL, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nStep Length (units)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3,  scales = "free")


p.SL_segclust<- ggplot(data = plot_data_SL_segclust, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_data_SL, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nStep Length (units)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3,  scales = "free")


p.SL_embc<- ggplot(data = plot_data_SL_embc, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_data_SL, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nStep Length (units)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3,  scales = "free")





### Make plot comparing HMM wrapped Cauchy distribs and Segclust2d normal distribs against true generating distribs
hmm.TA.df_parametric<- hmm.TA.params2_parametric %>% 
  bind_rows(., .id = "id") %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

segclust.TA.df_parametric<- segclust.params2_parametric %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par1, par3)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 15),
         track_length = factor(rep(c(1000,5000,10000), each = 15),
                               levels = c("1000","5000","10000","50000")))

embc.TA.df_parametric<- embc.params2_parametric %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par1, par3)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

true.TA.params_parametric<- data.frame(id = 1:3,
                                  TA.params,
                                  behavior = c("Encamped","ARS","Transit"),
                                  track_length = "True")


plot_data_TA_hmm<- 
  pmap_df(hmm.TA.df_parametric,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(-pi, pi, by = (pi/512)),
                   y = dwrappedcauchy(x, mu = circular(par1), rho = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_TA_hmm$behavior<- factor(plot_data_TA_hmm$behavior,
                                   levels = c("Encamped","ARS","Transit"))


plot_data_TA_segclust<- 
  pmap_df(segclust.TA.df_parametric,
          function(id, par2, par4, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, pi, by = (pi/256)),
                   y = dnorm(x, mean = par2, sd = par4),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_TA_segclust$behavior<- factor(plot_data_TA_segclust$behavior,
                                    levels = c("Encamped","ARS","Transit"))



plot_data_TA_embc<- 
  pmap_df(embc.TA.df_parametric,
          function(id, par2, par4, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, pi, by = (pi/256)),
                   y = dnorm(x, mean = par2, sd = par4),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_TA_embc$behavior<- factor(plot_data_TA_embc$behavior,
                                    levels = c("Encamped","ARS","Transit"))



true_plot_data_TA<- 
  pmap_df(true.TA.params_parametric,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(-pi, pi, by = (pi/512)),
                   y = dwrappedcauchy(x, mu = circular(par1), rho = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
true_plot_data_TA$behavior<- factor(true_plot_data_TA$behavior,
                                    levels = c("Encamped","ARS","Transit"))


# Plot overlapping densities per behavior (turning angles)
p.TA_hmm<- ggplot(data = plot_data_TA_hmm, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_data_TA, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nTurning Angle (radians)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3, scales = "free")


true_plot_TA_abs<- true_plot_data_TA %>% 
  filter(x >= 0) %>% 
  mutate_at("y", ~{.*2})
p.TA_segclust<- ggplot(data = plot_data_TA_segclust, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_TA_abs, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nTurning Angle (radians)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3, scales = "free")


p.TA_embc<- ggplot(data = plot_data_TA_embc, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_TA_abs, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nTurning Angle (radians)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3, scales = "free")



## Create composite plot w/ shared legend for HMM
p.comp_hmm<- plot_grid(p.SL_hmm + theme(legend.position="none") + 
                         ggtitle("HMM", subtitle = "Common distributions"),
                       NA,
                       p.TA_hmm + theme(legend.position="none") + 
                         ggtitle("", subtitle = ""),
                       rel_widths = c(1, 0.1, 1),
                       hjust = -1,
                       nrow = 1)

# extract the legend from one of the plots
legend.comp_hmm<- get_legend(p.SL_hmm + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(p.comp_hmm, legend.comp_hmm, ncol = 1, rel_heights = c(1, 0.1))

# ggsave("HMM Distribs_parametric.png", width = 7, height = 6, units = "in", dpi = 330)



## Create composite plot w/ shared legend for segclust
p.comp_segclust<- plot_grid(p.SL_segclust + theme(legend.position="none") + 
                              ggtitle("Segclust2d", subtitle = "Common distributions"),
                            NA,
                            p.TA_segclust + theme(legend.position="none") + 
                              ggtitle("", subtitle = ""),
                            rel_widths = c(1, 0.1, 1),
                            hjust = -1,
                            nrow = 1)

# extract the legend from one of the plots
legend.comp_segclust<- get_legend(p.SL_segclust + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(p.comp_segclust, legend.comp_segclust, ncol = 1, rel_heights = c(1, 0.1))

# ggsave("Segclust2d distribs_parametric.png", width = 7, height = 6, units = "in", dpi = 330)




## Create composite plot w/ shared legend for EMbC
p.comp_embc<- plot_grid(p.SL_embc + theme(legend.position="none") + 
                          ggtitle("EMbC", subtitle = "Common distributions"),
                            NA,
                            p.TA_embc + theme(legend.position="none") + 
                          ggtitle("", subtitle = ""),
                            rel_widths = c(1, 0.1, 1),
                            hjust = -1,
                            nrow = 1)

# extract the legend from one of the plots
legend.comp_embc<- get_legend(p.SL_embc + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(p.comp_embc, legend.comp_embc, ncol = 1, rel_heights = c(1, 0.1))

# ggsave("EMbC distribs_parametric.png", width = 7, height = 6, units = "in", dpi = 330)
