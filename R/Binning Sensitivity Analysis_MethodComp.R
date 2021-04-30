
##################################################
### Method Comparison for Sensitivity Analysis ###
##################################################

library(tidyverse)
library(wesanderson)
library(lubridate)
library(cowplot)
library(viridis)
library(ggnewscale)

source('helper functions.R')


# Load breakpoints
equal.5bins.brkpts<- read.csv("Sensitivity_allbreakpts_5bins_equal.csv")
equal.10bins.brkpts<- read.csv("Sensitivity_allbreakpts_10bins_equal.csv")
quant.5bins.brkpts<- read.csv("Bayesian_allbreakpts_weird.csv")  #from original analysis
quant.10bins.brkpts<- read.csv("Sensitivity_allbreakpts_10bins_quantile.csv")

# Load results
# bayes.res_weird<- read.csv("Modeled MM Sim Tracks w Behav_weird.csv")
# hmm.res_weird<- read.csv("HMM results_weird.csv")
# segclust.res_weird<- read.csv("Segclust2d results_weird.csv")

# Load true breakpoints
true.brkpts_weird<- read.csv("CRW_MM_sim_brkpts_weird.csv")



# Filter 'quant.5bins.brkpts' and 'true.brkpts_weird' to only include matching sim. tracks
quant.5bins.brkpts<- quant.5bins.brkpts %>% 
  filter(id %in% unique(equal.5bins.brkpts$id))
true.brkpts_weird<- true.brkpts_weird %>% 
  filter(id %in% unique(equal.5bins.brkpts$id))



###########################
### Compare Breakpoints ###
###########################

equal.5bins.brkpts$method<- rep("equal_5bins", nrow(equal.5bins.brkpts))
equal.10bins.brkpts$method<- rep("equal_10bins", nrow(equal.10bins.brkpts))
quant.5bins.brkpts$method<- rep("quant_5bins", nrow(quant.5bins.brkpts))
quant.10bins.brkpts$method<- rep("quant_10bins", nrow(quant.10bins.brkpts))
all.brkpts<- rbind(equal.5bins.brkpts, equal.10bins.brkpts, quant.5bins.brkpts, quant.10bins.brkpts)


brkpt.acc<- all.brkpts %>% 
  group_by(method, id, acc) %>% 
  filter(type == "Model") %>% 
  tally() %>% 
  mutate(freq = n/sum(n)) %>% 
  mutate(track_length = case_when(str_detect(id, "_2") ~ "5000",
                                  str_detect(id, "_3") ~ "10000")) %>% 
  group_by(method, track_length, id) %>% 
  filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
  summarise(freq = sum(freq)) %>%  #and calculate accuracy across all sims combined
  ungroup()

brkpt.acc$track_length<- brkpt.acc$track_length %>% 
  factor(., levels = c('5000','10000'))


#Accuracy (includes 'accurate' and 'accurate duplicate' classifications)
ggplot(brkpt.acc, aes(track_length, freq, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "Proportion of Accurate Breakpoints\n") +
  scale_fill_viridis_d("") +
  scale_color_viridis_d("") +
  ylim(0,1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank())


#calc mean accuracy per track_length and method
brkpt.acc %>% 
  group_by(method, track_length) %>% 
  summarise(mean=mean(freq), lo = min(freq), hi = max(freq))


# number of breakpoints missed
tot<- all.brkpts %>% 
  group_by(method, id) %>% 
  filter(type == "True") %>% 
  tally() %>% 
  ungroup()

all.brkpts %>% 
  group_by(method, id) %>% 
  filter(acc == "Missing") %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(n.true = tot$n, prop = n/n.true) %>% 
  mutate(track_length = case_when(str_detect(id, "_2") ~ "5000",
                                  str_detect(id, "_3") ~ "10000")) %>% 
  group_by(method, track_length) %>% 
  summarise(freq.mu = mean(prop), lo = min(prop), hi = max(prop)) %>%  #and calculate accuracy across all sims combined
  ungroup()

