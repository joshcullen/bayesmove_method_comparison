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
seg.time<- read.csv("data/Bayesian_elapsed_time_weird.csv")
lda.time<- read.csv("data/LDA_elapsed_time_weird.csv")
bcpa.time<- read.csv("data/BCPA_elapsed_time_weird.csv")
hmm.time<- read.csv("data/HMM_elapsed_time_weird.csv")
segclust.time<- read.csv("data/Segclust2d_elapsed_time_weird.csv")
embc.time<- read.csv("data/EMbC_elapsed_time_weird.csv")

# Load breakpoints
bayes.brkpts<- read.csv("data/Bayesian_allbreakpts_weird.csv")
bcpa.brkpts<- read.csv("data/BCPA_allbrkpts_weird.csv")
segclust.brkpts<- read.csv("data/Segclust2d allbrkpts_weird.csv")

# Load results
bayes.res_weird<- read.csv("data/Modeled MM Sim Tracks w Behav_weird.csv")
hmm.res_weird<- read.csv("data/HMM results_weird.csv")
segclust.res_weird<- read.csv("data/Segclust2d results_weird.csv")
embc.res_weird<- read.csv("data/EMbC results_weird.csv")

# Load true breakpoints
true.brkpts_weird<- read.csv("data/CRW_MM_sim_brkpts_weird.csv")


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
  scale_y_log10(breaks = c(10^(-2:2), 300),
                labels = c(0.01, 0.1, 1, 10, 100, 300)) +
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
embc.res_weird<- embc.res_weird %>% 
  bayesmove::df_to_list("id") %>% 
  map(., 
      ~rbind(c(unique(.$id), 0, 0, rep(NA, 4), unique(.$track_length), rep(NA, 2)),
             .x)) %>% 
  bind_rows()

# Assign identifiers by method and make consistent behavior colname
bayes.res_weird$method<- rep("M4", nrow(bayes.res_weird))
hmm.res_weird$method<- rep("HMM", nrow(hmm.res_weird))
segclust.res_weird$method<- rep("Segclust2d", nrow(segclust.res_weird))
embc.res_weird$method<- rep("EMbC", nrow(embc.res_weird))

bayes.res_weird<- bayes.res_weird %>% 
  rename(state = behav) %>% 
  mutate_at("state", ~factor(., levels = 1:3))
hmm.res_weird<- hmm.res_weird %>% 
  rename(state = hmm.state, id = ID)
embc.res_weird.agg<- embc.res_weird %>%  #combine HL and HH states into one
  mutate_at("state", ~case_when(state == "LH" ~ 1,
                                state == "LL" ~ 2,
                                state == "HL" | state == "HH" ~ 3))


# Modify hmm.res_weird, segclust.res_weird, and embc.res_weird to be same as bayes.res_weird format
# calc proportions of behaviors by true time segment and then identify dominant behavior
hmm.res_weird2<- hmm.res_weird %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_weird) %>% 
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

hmm.res_weird3<- hmm.res_weird %>% 
  dplyr::select(-state) %>% 
  bayesmove::df_to_list("id") %>% 
  map2(.,
       bayesmove::df_to_list(hmm.res_weird2, "id") %>% 
         map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()




segclust.res_weird2<- segclust.res_weird %>% 
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


segclust.res_weird3<- segclust.res_weird %>% 
  dplyr::select(-state) %>% 
  bayesmove::df_to_list("id") %>% 
  map2(.,
       bayesmove::df_to_list(segclust.res_weird2, "id") %>% 
         map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()







embc.res_weird2<- embc.res_weird.agg %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_weird) %>% 
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

embc.res_weird3<- embc.res_weird.agg %>% 
  dplyr::select(-state) %>% 
  bayesmove::df_to_list("id") %>% 
  map2(.,
       bayesmove::df_to_list(embc.res_weird2, "id") %>% 
         map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()



# Combine all datasets
res_weird<- rbind(bayes.res_weird[,c("id","behav_fine","behav_coarse","track_length","state",
                                     "method")],
            hmm.res_weird3[,c("id","behav_fine","behav_coarse","track_length","state","method")],
            segclust.res_weird3[,c("id","behav_fine","behav_coarse","track_length","state",
                                   "method")],
            embc.res_weird3[,c("id","behav_fine","behav_coarse","track_length","state",
                                   "method")])
res_weird$method<- factor(res_weird$method, levels = c('M4','HMM','Segclust2d','EMbC'))



## Overall

#Coarse-scale behavior
res_weird %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate_at("track_length", as.numeric) %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))

summ.stats_coarse_weird<- res_weird %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate_at("track_length", as.numeric) %>%
  mutate(acc = n/track_length) %>% 
  ungroup()

summ.stats_coarse_weird$track_length<- summ.stats_coarse_weird$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))

p.coarse_weird<- ggplot(summ.stats_coarse_weird, aes(track_length,acc,fill = method,
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
bayes.list_weird<- bayesmove::df_to_list(bayes.res_weird, "id")

true.behavior.long_weird<- list()
for (i in 1:length(bayes.list_weird)) {
  true.behavior.long_weird[[i]]<- 
    data.frame(true.tseg = rep(1:(bayes.list_weird[[i]]$track_length[1]/100), each = 300),
                                       behav_coarse = rep(bayes.list_weird[[i]]$behav_coarse[-1],
                                                          each = 3),
                                       behav_fine = rep(bayes.list_weird[[i]]$behav_fine[-1],
                                                        each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(bayes.list_weird[[i]]$track_length[1]),
                                                   each = 3))
  
  true.behavior.long_weird[[i]]$prop<- 0.1
  
  cond<- true.behavior.long_weird[[i]][,"behav_coarse"]
  ind<- which(true.behavior.long_weird[[i]][,"behavior"] == cond)
  
  true.behavior.long_weird[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long_weird[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long_weird[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long_weird[[i]][cond2, "prop"]<- true.behavior.long_weird[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long_weird)<- names(bayes.list_weird)

## True proportions for HMMs (from time segments using true breakpoints)
hmm.props_weird<- hmm.res_weird %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_weird) %>% 
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
for (i in 1:length(unique(as.character(hmm.props_weird$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=true.behavior.long_weird[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = hmm.props_weird %>% 
                  filter(id == unique(as.character(hmm.res_weird$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(hmm.res_weird$id))[i]) +
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
segclust.props_weird<- segclust.res_weird %>% 
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
for (i in 1:length(unique(as.character(segclust.props_weird$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=true.behavior.long_weird[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = segclust.props_weird %>% 
                  filter(id == unique(as.character(segclust.res_weird$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(segclust.res_weird$id))[i]) +
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
embc.props_weird<- embc.res_weird.agg %>% 
  bayesmove::df_to_list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_weird) %>% 
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
for (i in 1:length(unique(as.character(embc.props_weird$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=true.behavior.long_weird[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = embc.props_weird %>% 
                  filter(id == unique(as.character(embc.res_weird$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(embc.res_weird$id))[i]) +
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
bayes.props_weird<- bayes.res_weird %>% 
  rename(Encamped = X1, ARS = X2, Transit = X3) %>% 
  drop_na() %>% 
  pivot_longer(., cols = c(Encamped, ARS, Transit), names_to = "behavior",
               values_to = "prop") %>% 
  dplyr::select(id, tseg.y, behavior, prop, time1) %>% 
  mutate_at("behavior", ~recode(., 'Encamped' = 1, 'ARS' = 2, 'Transit' = 3))
bayes.props_weird$time1<- bayes.props_weird$time1 - 1


## Calculate RMSE
true.behavior_weird<- true.behavior.long_weird %>% 
  bind_rows(.id = "id")

hmm.rmse_weird<- vector()
for (i in 1:length(unique(hmm.res_weird$id))) {
  ind<- unique(as.character(hmm.res_weird$id))[i]
  
  hmm.rmse_weird[i]<- sqrt(sum((hmm.props_weird[hmm.props_weird$id == ind, "prop"] - 
                           true.behavior_weird[true.behavior_weird$id == ind, "prop"])^2) / nrow(hmm.props_weird[hmm.props_weird$id == ind,]))
}



segclust.rmse_weird<- vector()
for (i in 1:length(unique(segclust.res_weird$id))) {
  ind<- unique(as.character(segclust.res_weird$id))[i]
  
  segclust.rmse_weird[i]<- sqrt(sum((segclust.props_weird[segclust.props_weird$id == ind, "prop"] - 
                                  true.behavior_weird[true.behavior_weird$id == ind, "prop"])^2) / nrow(segclust.props_weird[segclust.props_weird$id == ind,]))
}



embc.rmse_weird<- vector()
for (i in 1:length(unique(embc.res_weird.agg$id))) {
  ind<- unique(as.character(embc.res_weird.agg$id))[i]
  
  embc.rmse_weird[i]<- sqrt(sum((embc.props_weird[embc.props_weird$id == ind, "prop"] - 
                                       true.behavior_weird[true.behavior_weird$id == ind, "prop"])^2) / nrow(embc.props_weird[embc.props_weird$id == ind,]))
}



bayes.rmse_weird<- vector()
for (i in 1:length(unique(bayes.res_weird$id))) {
  ind<- unique(as.character(bayes.res_weird$id))[i]
  
  bayes.rmse_weird[i]<- sqrt(sum((bayes.props_weird[bayes.props_weird$id == ind, "prop"] - 
                             true.behavior_weird[true.behavior_weird$id == ind, "prop"])^2) /
                        nrow(bayes.props_weird[bayes.props_weird$id == ind,]))
}


segclust.rmse.df<- data.frame(id = unique(as.character(segclust.res_weird$id)),
                              track_length = factor(rep(c(1000,5000,10000), each = 5),
                                                    levels = c("1000","5000","10000","50000")),
                              rmse = segclust.rmse_weird,
                              method = rep("Segclust2d", 15))
rmse.df_weird<- data.frame(id = rep(unique(as.character(hmm.res_weird$id)), 3),
                    track_length = factor(rep(rep(c(1000,5000,10000,50000), each = 5), 3),
                                          levels = c("1000","5000","10000","50000")),
                    rmse = c(bayes.rmse_weird, hmm.rmse_weird, embc.rmse_weird),
                    method = rep(c("M4","HMM","EMbC"), each = 20))
rmse.df_weird<- rbind(rmse.df_weird, segclust.rmse.df)
rmse.df_weird$method<- factor(rmse.df_weird$method, levels = c('M4','HMM','Segclust2d','EMbC'))


#Plot results

p.rmse_weird<- ggplot(rmse.df_weird, aes(track_length, rmse, fill = method, color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "point", shape = 15, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="\nTrack Length (observations)", y = "RMSE\n") +
  scale_fill_manual("", values = pal1[c(1:4)]) +
  scale_color_manual("", values = pal1[c(1:4)]) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.55)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))



# create composite
res.comp<- plot_grid(NULL, NULL, NULL,
          p.time + theme(legend.position = "none"), NULL, p.brk,
          # NULL, NULL, NULL,
          # p.coarse, NULL, p.rmse,
          NULL, NULL, NULL,
          p.coarse_weird, NULL, p.rmse_weird,
          align = "hv", nrow = 4, rel_widths = c(1,0.1,1), rel_heights = c(0.2,1,0.1,1))

# extract the legend from one of the plots
legend.comp_res<- get_legend(p.time + theme(legend.position="top",
                                            legend.text = element_text(size = 14)))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(legend.comp_res, res.comp, ncol = 1, rel_heights = c(0.1, 1))

# ggsave("Figure 3 (method comparison)_updatedWeird.png", width = 12, height = 9, units = "in", dpi = 330)









###############################################################################
### Compare Characterization of Step Length and Turning Angle Distributions ###
###############################################################################

library(circular)


## Figure S1 ##

### Define bin limits
dat_weird<- read.csv("data/CRW_MM_sim_weird.csv", as.is = T)
dat_weird$dt<- 3600
names(dat_weird)[4:5]<- c("dist","rel.angle")

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

dist.bin.lims=quantile(dat_weird[dat_weird$dt == 3600,]$dist,
                       c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins

hmm.SL.params_weird<- read.csv("data/HMM result step params_weird.csv", as.is = T)  
hmm.TA.params_weird<- read.csv("data/HMM result angle params_weird.csv", as.is = T)  

#Manipulate param dfs to reformat for extract.behav.props()
hmm.SL.params2_weird<- list()
for (i in 1:nrow(hmm.SL.params_weird)) {
  hmm.SL.params2_weird[[i]]<- data.frame(par1 = as.numeric(hmm.SL.params_weird[i,2:4]),
                                         par2 = as.numeric(hmm.SL.params_weird[i,5:7]))
}
names(hmm.SL.params2_weird)<- hmm.SL.params_weird$X

hmm.TA.params2_weird<- list()
for (i in 1:nrow(hmm.TA.params_weird)) {
  hmm.TA.params2_weird[[i]]<- data.frame(par1 = as.numeric(hmm.TA.params_weird[i,2:4]),
                                         par2 = as.numeric(hmm.TA.params_weird[i,5:7]))
}
names(hmm.TA.params2_weird)<- hmm.TA.params_weird$X



#convert mean and sd to shape and rate params for gamma dist
for (j in 1:length(hmm.SL.params2_weird)) {
  for (i in 1:nrow(hmm.SL.params2_weird[[j]])) {
    shape<- (hmm.SL.params2_weird[[j]][i,1]^2) / (hmm.SL.params2_weird[[j]][i,2]^2)
    rate<- hmm.SL.params2_weird[[j]][i,1] / (hmm.SL.params2_weird[[j]][i,2]^2)
    
    params<- c(shape, rate)
    
    hmm.SL.params2_weird[[j]][i,]<- params
  }
}



segclust.params_weird<- read.csv("data/Segclust result state params_weird.csv", as.is = T) %>% 
  dplyr::select(-1)

#Manipulate param dfs to reformat for extract.behav.props()
segclust.params2_weird<- list()
for (i in 1:length(unique(segclust.params_weird$id))) {
  tmp<- segclust.params_weird %>% 
    filter(id == unique(segclust.params_weird$id)[i])
  segclust.params2_weird[[i]]<- tmp[,4:7]
  names(segclust.params2_weird[[i]])<- paste0("par", 1:4)
}
names(segclust.params2_weird)<- unique(segclust.params_weird$id)





samp.size<- embc.res_weird %>%  #for use in merging Gaussians of HL and HH
  mutate_at("state", factor, levels = c("LL","LH","HL","HH")) %>% 
  group_by(id, state) %>% 
  tally() %>% 
  drop_na()
embc.params_weird<- read.csv("data/EMbC result state params_weird.csv", as.is = T) %>% 
  dplyr::select(-1) %>% 
  mutate(n = samp.size$n)

#Manipulate param dfs to reformat for extract.behav.props()
embc.params2_weird<- list()
for (i in 1:length(unique(embc.params_weird$id))) {
  tmp<- embc.params_weird %>% 
    filter(id == unique(embc.params_weird$id)[i])
  
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
  
  embc.params2_weird[[i]]<- rbind(tmp[2:1, 2:5], H)
  names(embc.params2_weird[[i]])<- paste0("par", 1:4)
}
names(embc.params2_weird)<- unique(embc.params_weird$id)





SL.params_weird<- data.frame(par1 = c(0.25, 2, exp(2)), par2 = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0), par2 = c(0.5, pi, 1))

true.b_weird<- extract.behav.props_weird(params = list(SL.params_weird, TA.params),
                                         lims = list(dist.bin.lims, angle.bin.lims),
                                         behav.names = c("Encamped","ARS","Transit"))


hmm.b_weird<- map2(hmm.SL.params2_weird, hmm.TA.params2_weird,
                   ~extract.behav.props(params = list(.x, .y),
                                        lims = list(dist.bin.lims, angle.bin.lims),
                                        behav.names = c("Encamped","ARS","Transit"))
)


segclust.b_weird<- map(segclust.params2_weird,
                       ~extract.behav.props_norm(params = .x,
                                                 lims = list(dist.bin.lims, angle.bin.lims),
                                                 behav.names = c("Encamped","ARS","Transit"))
)


embc.b_weird<- map(embc.params2_weird,
                   ~extract.behav.props_norm(params = .x,
                                             lims = list(dist.bin.lims, angle.bin.lims),
                                             behav.names = c("Encamped","ARS","Transit"))
)




## Bayesian
behav.res_weird<-  read.csv("data/CRW MM LDA Phi values_weird.csv", as.is = T)
behav.res_weird<- bayesmove::df_to_list(behav.res_weird, "id")

behav.order_weird<- read.csv("data/CRW MM LDA behavior order_weird.csv", as.is = T)
behav.order_weird<- bayesmove::df_to_list(behav.order_weird, "id")
behav.order_weird<- map(behav.order_weird, as.numeric)


#calculate true proportions of SL and TA by behavior for all bins for ID 2_2 for Bayesian model
bayes.b_weird<- behav.res_weird[[15]] %>% rename(., var = param)
bayes.b_weird$behav<- bayes.b_weird$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "Transit") %>%
  str_replace_all(., "3", "ARS") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.distfit<- ggplot(true.b_weird, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = hmm.b_weird[[15]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[2], color="black", stroke=1) +
  geom_point(data = segclust.b_weird[[15]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[3], color="black", stroke=1) +
  geom_point(data = embc.b_weird[[15]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[4], color="black", stroke=1) +
  geom_point(data = bayes.b_weird, aes(x=bin, y=prop, group = behav), pch=21, size = 2,
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
bayes.rmse_weird<- list()
for (i in 1:length(behav.res_weird)) {
  bayes.b_weird<- behav.res_weird[[i]] %>% 
    rename(., var = param)
  
  behavs<- c("Encamped","ARS","Transit")
  
  bayes.b_weird$behav<- bayes.b_weird$behav %>% 
    factor(levels = behav.order_weird[[i]]) %>% 
    str_replace_all(., "1", behavs[behav.order_weird[[i]][1]]) %>% 
    str_replace_all(., "2", behavs[behav.order_weird[[i]][2]]) %>%
    str_replace_all(., "3", behavs[behav.order_weird[[i]][3]]) %>%
    factor(., levels = c("Encamped","ARS","Transit"))
  bayes.b_weird<- bayes.b_weird %>% 
    group_by(var) %>% 
    arrange(behav) %>% 
    ungroup()
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(bayes.b_weird$var))) {
    tmp<- bayes.b_weird %>% 
      filter(var == unique(bayes.b_weird$var)[j])
    true.tmp<- true.b_weird %>% 
      filter(var == unique(bayes.b_weird$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  bayes.rmse_weird[[i]]<- data.frame(vec)
  names(bayes.rmse_weird)[i]<- names(behav.res_weird)[i]
}
bayes.rmse_weird<- bind_rows(bayes.rmse_weird) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))




#HMM
hmm.rmse_weird<- list()
for (i in 1:length(behav.res_weird)) {
  hmm.b_weird<- extract.behav.props(params = list(hmm.SL.params2_weird[[i]],
                                                  hmm.TA.params2_weird[[i]]),
                                    lims = list(dist.bin.lims, angle.bin.lims),
                                    behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(hmm.b_weird$var))) {
    tmp<- hmm.b_weird %>% 
      filter(var == unique(hmm.b_weird$var)[j])
    true.tmp<- true.b_weird %>% 
      filter(var == unique(hmm.b_weird$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  hmm.rmse_weird[[i]]<- data.frame(vec)
  names(hmm.rmse_weird)[i]<- names(behav.res_weird)[i]
}
hmm.rmse_weird<- bind_rows(hmm.rmse_weird) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))





#Segclust2d
segclust.rmse_weird<- list()
for (i in 1:length(segclust.params2_weird)) {
  segclust.b_weird<- extract.behav.props_norm(params = segclust.params2_weird[[i]],
                                                  lims = list(dist.bin.lims, angle.bin.lims),
                                                  behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(segclust.b_weird$var))) {
    tmp<- segclust.b_weird %>% 
      filter(var == unique(segclust.b_weird$var)[j])
    true.tmp<- true.b_weird %>% 
      filter(var == unique(segclust.b_weird$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  segclust.rmse_weird[[i]]<- data.frame(vec)
  names(segclust.rmse_weird)[i]<- names(behav.res_weird)[i]
}
segclust.rmse_weird<- bind_rows(segclust.rmse_weird) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 15))





#EMbC
embc.rmse_weird<- list()
for (i in 1:length(embc.params2_weird)) {
  embc.b_weird<- extract.behav.props_norm(params = embc.params2_weird[[i]],
                                                  lims = list(dist.bin.lims, angle.bin.lims),
                                                  behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(embc.b_weird$var))) {
    tmp<- embc.b_weird %>% 
      filter(var == unique(embc.b_weird$var)[j])
    true.tmp<- true.b_weird %>% 
      filter(var == unique(embc.b_weird$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  embc.rmse_weird[[i]]<- data.frame(vec)
  names(embc.rmse_weird)[i]<- names(behav.res_weird)[i]
}
embc.rmse_weird<- bind_rows(embc.rmse_weird) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))




segclust.rmse.df<- data.frame(id = rep(unique(as.character(segclust.res_weird$id)), each = 2),
                              track_length = factor(rep(rep(c(1000,5000,10000), each = 5),
                                                        each = 2),
                                                    levels = c("1000","5000","10000","50000")),
                              rmse = segclust.rmse_weird,
                              method = rep("Segclust2d", 30))
rmse.df_weird<- data.frame(id = rep(rep(names(behav.res_weird), 3), each = 2),
                           track_length = factor(rep(rep(rep(c(1000,5000,10000,50000),
                                                             each = 5), 3), each = 2),
                                                 levels = c("1000","5000","10000","50000")),
                           rmse = rbind(bayes.rmse_weird, hmm.rmse_weird, embc.rmse_weird),
                           method = rep(c("M4","HMM","EMbC"), each = 40))
rmse.df_weird<- rbind(rmse.df_weird, segclust.rmse.df)
rmse.df_weird$method<- factor(rmse.df_weird$method, levels = c('M4','HMM','Segclust2d','EMbC'))


p.rmse_weird<- ggplot(rmse.df_weird, aes(track_length, rmse.value, fill = method,
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
png("Figure S1 (rmse from sim)_updated.png", width = 14.5, height = 5.5, units = "in", res = 330)

grid.arrange(p.distfit, p.rmse_weird, heights = c(0.2, 1),
             widths = c(1, 0.2, 1, 0.5),
             layout_matrix = rbind(c(NA, NA, NA, NA),
                                   c(3, NA, 4, 4)))
dev.off()










## Figure 5 ##

### Make plot comparing HMM gamma distribs and Segclust2d/EMbC normal distribs against true distributions
hmm.SL.df_weird<- hmm.SL.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

segclust.SL.df_weird<- segclust.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par2, par4)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 15),
         track_length = factor(rep(c(1000,5000,10000), each = 15),
                               levels = c("1000","5000","10000","50000")))

embc.SL.df_weird<- embc.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par2, par4)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

true.SL.params_weird<- data.frame(id = 1:3,
                                  SL.params_weird,
                                  behavior = c("Encamped","ARS","Transit"),
                                  track_length = "True")


plot_data_SL_hmm<- 
  pmap_df(hmm.SL.df_weird,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dgamma(x, shape = par1, rate = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_SL_hmm$behavior<- factor(plot_data_SL_hmm$behavior, levels = c("Encamped","ARS","Transit"))


plot_data_SL_segclust<- 
  pmap_df(segclust.SL.df_weird,
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
  pmap_df(embc.SL.df_weird,
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
  pmap_df(true.SL.params_weird,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dtnorm(x, mean1 = par1, sd1 = par2, lo = 0, hi = Inf),
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





### Make plot comparing HMM wrapped Cauchy distribs and Segclust2d/EMbC normal distribs against true generating distribs
hmm.TA.df_weird<- hmm.TA.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

segclust.TA.df_weird<- segclust.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par1, par3)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 15),
         track_length = factor(rep(c(1000,5000,10000), each = 15),
                               levels = c("1000","5000","10000","50000")))

embc.TA.df_weird<- embc.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  dplyr::select(-c(par1, par3)) %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

true.TA.params_weird<- data.frame(id = 1:3,
                                  TA.params,
                                  behavior = c("Encamped","ARS","Transit"),
                                  track_length = "True")


plot_data_TA_hmm<- 
  pmap_df(hmm.TA.df_weird,
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
  pmap_df(segclust.TA.df_weird,
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
  pmap_df(embc.TA.df_weird,
          function(id, par2, par4, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, pi, by = (pi/256)),
                   y = dnorm(x, mean = par2, sd = par4),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_TA_embc$behavior<- factor(plot_data_TA_embc$behavior,
                                        levels = c("Encamped","ARS","Transit"))



#True Encamped TA
true_plot_data_TA1<- 
  pmap_df(true.TA.params_weird[1,],
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 1, length.out = 1025),
                   y = dbeta(x, shape1 = par1, shape2 = par2),
                   behavior = behavior,
                   track_length = track_length)
          }) %>% 
  mutate_at("x", ~{.*2*pi-pi}) %>% 
  mutate_at("y", ~{. / (2*pi)})

#True ARS TA
true_plot_data_TA2<- 
  pmap_df(true.TA.params_weird[2,],
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(-pi, pi, by = (2*pi/512)),
                   y = dunif(x, min = par1, max = par2),
                   behavior = behavior,
                   track_length = track_length)
          })

#True Transit TA
true_plot_data_TA3<- 
  pmap_df(true.TA.params_weird[3,],
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(-pi, pi, by = (2*pi/512)),
                   y = dtnorm(x, mean1 = par1, sd1 = par2, lo = -pi, hi = pi),
                   behavior = behavior,
                   track_length = track_length)
          })

#Merge all TA estimates
true_plot_data_TA<- rbind(true_plot_data_TA1, true_plot_data_TA2, true_plot_data_TA3)

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
                         ggtitle("HMM", subtitle = "Uncommon distributions"),
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

# ggsave("HMM Distribs_weird.png", width = 7, height = 6, units = "in", dpi = 330)



## Create composite plot w/ shared legend for segclust
p.comp_segclust<- plot_grid(p.SL_segclust + theme(legend.position="none") + 
                              ggtitle("Segclust2d", subtitle = "Uncommon distributions"),
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

# ggsave("Segclust2d distribs_weird.png", width = 7, height = 6, units = "in", dpi = 330)



## Create composite plot w/ shared legend for EMbC
p.comp_embc<- plot_grid(p.SL_embc + theme(legend.position="none") + 
                          ggtitle("EMbC", subtitle = "Uncommon distributions"),
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

# ggsave("EMbC distribs_weird.png", width = 7, height = 6, units = "in", dpi = 330)
