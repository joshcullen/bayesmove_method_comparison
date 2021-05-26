
###################################################################
### Comparison of thresholds used to assess breakpoint accuracy ###
###################################################################

library(bayesmove)
library(tidyverse)
library(wesanderson)

source('R/helper functions.R')


### Load breakpoint estimates
bayes.brkpts<- read.csv("data/Bayesian_allbreakpts_weird.csv")
bcpa.brkpts<- read.csv("data/BCPA_allbrkpts_weird.csv")
segclust.brkpts<- read.csv("data/Segclust2d allbrkpts_weird.csv")

# Load true breakpoints
true.brkpts<- read.csv("data/CRW_MM_sim_brkpts_weird.csv")





### Accurate w/in 5 observations, missing beyond 15 observations ###

#Bayesian
bayes.list<- bayes.brkpts %>%  #convert to list of only a breakpoint vector per ID
  filter(type == "Model") %>% 
  df_to_list(., "id") %>% 
  map(., ~{.$brks})

# Compare true vs modeled breakpoints
bayes.list2<- list()
for (i in 1:length(bayes.list)) {
  bayes.list2[[i]]<- brkpt.accuracy(model.brkpts = bayes.list[[i]],
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 5, dup.tol = 1, miss.tol = 15)
}

names(bayes.list2)<- names(bayes.list)
bayes.5_15<- bind_rows(bayes.list2, .id = 'id') %>% 
  mutate(method = "M4", scenario = 1)


#BCPA
bcpa.list<- bcpa.brkpts %>%  #convert to list of only a breakpoint vector per ID
  filter(type == "Model") %>% 
  df_to_list(., "id") %>% 
  map(., ~{.$brks})

# Compare true vs modeled breakpoints
bcpa.list2<- list()
for (i in 1:length(bcpa.list)) {
  bcpa.list2[[i]]<- brkpt.accuracy(model.brkpts = bcpa.list[[i]],
                                    true.brkpts = true.brkpts[i,-1],
                                    acc.tol = 5, dup.tol = 1, miss.tol = 15)
}

names(bcpa.list2)<- names(bcpa.list)
bcpa.5_15<- bind_rows(bcpa.list2, .id = 'id') %>% 
  mutate(method = "BCPA", scenario = 1)


#Segclust2d
segclust.list<- segclust.brkpts %>%  #convert to list of only a breakpoint vector per ID
  filter(type == "Model") %>% 
  df_to_list(., "id") %>% 
  map(., ~{.$brks})

# Compare true vs modeled breakpoints
segclust.list2<- list()
for (i in 1:length(segclust.list)) {
  segclust.list2[[i]]<- brkpt.accuracy(model.brkpts = segclust.list[[i]],
                                    true.brkpts = true.brkpts[i,-1],
                                    acc.tol = 5, dup.tol = 1, miss.tol = 15)
}

names(segclust.list2)<- names(segclust.list)[1:15]
segclust.5_15<- bind_rows(segclust.list2, .id = 'id') %>% 
  mutate(method = "Segclust2d", scenario = 1)






### Accurate w/in 20 observations, missing beyond 50 observations ###

#Bayesian

# Compare true vs modeled breakpoints
bayes.list2<- list()
for (i in 1:length(bayes.list)) {
  bayes.list2[[i]]<- brkpt.accuracy(model.brkpts = bayes.list[[i]],
                                    true.brkpts = true.brkpts[i,-1],
                                    acc.tol = 20, dup.tol = 1, miss.tol = 50)
}

names(bayes.list2)<- names(bayes.list)
bayes.20_50<- bind_rows(bayes.list2, .id = 'id') %>% 
  mutate(method = "M4", scenario = 3)


#BCPA

# Compare true vs modeled breakpoints
bcpa.list2<- list()
for (i in 1:length(bcpa.list)) {
  bcpa.list2[[i]]<- brkpt.accuracy(model.brkpts = bcpa.list[[i]],
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 20, dup.tol = 1, miss.tol = 50)
}

names(bcpa.list2)<- names(bcpa.list)
bcpa.20_50<- bind_rows(bcpa.list2, .id = 'id') %>% 
  mutate(method = "BCPA", scenario = 3)


#Segclust2d

# Compare true vs modeled breakpoints
segclust.list2<- list()
for (i in 1:length(segclust.list)) {
  segclust.list2[[i]]<- brkpt.accuracy(model.brkpts = segclust.list[[i]],
                                       true.brkpts = true.brkpts[i,-1],
                                       acc.tol = 20, dup.tol = 1, miss.tol = 50)
}

names(segclust.list2)<- names(segclust.list)[1:15]
segclust.20_50<- bind_rows(segclust.list2, .id = 'id') %>% 
  mutate(method = "Segclust2d", scenario = 3)







### Wrangle data and plot for comparison

bayes.brkpts2<- bayes.brkpts %>% 
  mutate(method = "M4", scenario = 2)
bcpa.brkpts2<- bcpa.brkpts %>% 
  mutate(method = "BCPA", scenario = 2)
segclust.brkpts2<- segclust.brkpts %>% 
  mutate(method = "Segclust2d", scenario = 2)


all.brkpts<- rbind(bayes.5_15, bayes.brkpts2, bayes.20_50,
               bcpa.5_15, bcpa.brkpts2, bcpa.20_50,
               segclust.5_15, segclust.brkpts2, segclust.20_50)


### Compare accuracy
brkpt.acc<- all.brkpts %>% 
  group_by(scenario, method, id, acc) %>% 
  filter(type == "Model") %>% 
  tally() %>% 
  mutate(freq = n/sum(n)) %>% 
  mutate(track_length = case_when(str_detect(id, "_1") ~ "1000",
                                  str_detect(id, "_2") ~ "5000",
                                  str_detect(id, "_3") ~ "10000",
                                  str_detect(id, "_4") ~ "50000")) %>% 
  group_by(scenario, method, track_length, id) %>% 
  filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
  summarise(freq = sum(freq)) %>%  #and calculate accuracy across all sims combined
  ungroup()

brkpt.acc$track_length<- brkpt.acc$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))
brkpt.acc$method<- factor(brkpt.acc$method, levels = c('M4','Segclust2d','BCPA'))

pal1<- c(wes_palette("Darjeeling1", 5)[-4], "mediumorchid")

#Accuracy (includes 'accurate' and 'accurate duplicate' classifications)
ggplot(brkpt.acc, aes(track_length, freq, fill = method, color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x = "\nTrack Length", y = "Accuracy of Breakpoints\n") +
  scale_fill_manual("", values = pal1[c(1,3,5)]) +
  scale_color_manual("", values = pal1[c(1,3,5)]) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  ylim(0,1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = "top") +
  facet_wrap(~ scenario)

# ggsave("Breakpoint threshold comparison_accuracy.png", width = 9, height = 5, units = "in", dpi = 330)



#Compare missingness
n.true.brks<- all.brkpts %>% 
  group_by(scenario, method, id) %>% 
  filter(type == "True") %>% 
  count(.drop = FALSE) %>% 
  ungroup() %>% 
  select(n)

brkpt.miss<- all.brkpts %>% 
  filter(type == "True") %>% 
  mutate_at("acc", factor, levels = c("True","Missing")) %>% 
  group_by(scenario, method, id, acc, .drop = FALSE) %>% 
  count() %>% 
  filter(acc == "Missing") %>% 
  ungroup() %>% 
  mutate(tot.n = n.true.brks$n, freq = n/tot.n) %>% 
  mutate(track_length = case_when(str_detect(id, "_1") ~ "1000",
                                  str_detect(id, "_2") ~ "5000",
                                  str_detect(id, "_3") ~ "10000",
                                  str_detect(id, "_4") ~ "50000")) %>% 
  group_by(scenario, method, track_length, id) %>% 
  # filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
  summarise(freq = sum(freq)) %>%  #and calculate accuracy across all sims combined
  ungroup()

brkpt.miss$track_length<- brkpt.miss$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))
brkpt.miss$method<- factor(brkpt.miss$method, levels = c('M4','Segclust2d','BCPA'))


#Accuracy (includes 'accurate' and 'accurate duplicate' classifications)
ggplot(brkpt.miss, aes(track_length, freq, fill = method, color = method)) +
  geom_boxplot(width = 0.75) +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x = "\nTrack Length", y = "Proportion of Missed Breakpoints\n") +
  scale_fill_manual("", values = pal1[c(1,3,5)]) +
  scale_color_manual("", values = pal1[c(1,3,5)]) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  ylim(0,1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = "top") +
  facet_wrap(~ scenario)

# ggsave("Breakpoint threshold comparison_miss.png", width = 9, height = 5, units = "in", dpi = 330)
