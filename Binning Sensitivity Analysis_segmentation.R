
##################################
#### Run Sensitivity Analysis ####
##################################

### Running sensitivity analysis on the method used for binning step lengths (i.e, # of bins, method)

### Evaluating simulated tracks (w/ weird distributions) that have 5k or 10k observations

### Since we use 5 bins w/ quantiles to discretize step lengths in the body of the manuscript, we will test 5 bins of equal widths, 10 bins using quantiles, and 10 bins of equal widths

### Results are provided as summary statistics on the estimated breakpoints as well as accuracy of state assignments



library(bayesmove)
library(tidyverse)
library(lubridate)
library(future)
library(furrr)
library(progressr)

source('helper functions.R')


#load and manipulate data
dat<- read.csv("CRW_MM_sim_weird.csv", as.is = T)
true.brkpts<- read.csv("CRW_MM_sim_brkpts_weird.csv", as.is = T)
dat$dt<- 3600  #set time step
names(dat)[4:5]<- c("dist","rel.angle")  #change names for step length and turning angle
dat.list<- df_to_list(dat=dat, ind = "id")
dat.list<- dat.list[6:15]  #only keep tracks w/ 5k or 10k obs

#filter data for tstep of interest
behav.list<- filter_time(dat.list = dat.list, int = 3600)  #add move params and filter by 3600 s interval

#define bin number and limits for turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins



####################################
### Assess 5 bins of equal width ###
####################################

set.seed(1)

#define bin number and limits for step lengths
hist(dat$dist, 100)
dist.bin.lims = c(seq(0, 40, length.out = 5), max(dat$dist, na.rm = T))  #5 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=dist)) +
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = dist.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nStep Length", y = "Density\n")

ggplot(behav.df, aes(x=rel.angle)) +
  geom_density(fill = "indianred") +
  geom_vline(xintercept = angle.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nTurning Angle (rad)", y = "Density\n")


#assign bins to obs
behav.list1<- map(behav.list, discrete_move_var, lims = list(dist.bin.lims, angle.bin.lims),
                 varIn = c("dist", "rel.angle"), varOut = c("SL", "TA"))
behav.list2<- lapply(behav.list1, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment




## Run RJMCMC

#prior
alpha = 1

ngibbs = 40000

handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

# track_length == 5000; takes 7.5 min for 40000 iterations
dat.res_5k<- segment_behavior(data = behav.list2[1:5], ngibbs = ngibbs, nbins = c(5,8),
                                alpha = alpha)

# track_length == 10000; takes 16 min for 40000 iterations
dat.res_10k<- segment_behavior(data = behav.list2[6:10], ngibbs = ngibbs, nbins = c(5,8),
                                 alpha = alpha)

plan(sequential)  #closes background workers


## If sims analyzed separately, merge all runs together in single list
dat.res<- mapply(rbind, dat.res_5k, dat.res_10k, SIMPLIFY=FALSE)


## Reclassify and restructure data as needed
dat.res$nbrks[,2:ncol(dat.res$nbrks)]<- apply(dat.res$nbrks[,2:ncol(dat.res$nbrks)], 2,
                                              function(x) as.numeric(as.character(x)))
dat.res$LML[,2:ncol(dat.res$LML)]<- apply(dat.res$LML[,2:ncol(dat.res$LML)], 2,
                                          function(x) as.numeric(as.character(x)))
dat.res$brkpts<- dat.res$brkpts %>% 
  split(., row(.)) %>% 
  flatten()  #deals with dimensional list
names(dat.res$brkpts)<- dat.res$nbrks[,1]


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")


##Determine maximum a posteriori (MAP) estimate for selecting breakpoints
MAP.est<- get_MAP(dat.res$LML, nburn = ngibbs/2)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

## Visualize breakpoints over data streams
plot_breakpoints(data = behav.list1, as_date = FALSE, var_names = c("dist","rel.angle"),
                 var_labels = c("Step Length","Turning Angle"), brkpts = brkpts)


# Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:nrow(brkpts)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = brkpts[i,-1], true.brkpts = true.brkpts[i+5,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}


#assign time seg and make as DF
dat_out<- assign_tseg(dat = behav.list1, brkpts = brkpts)

names(all.brkpts)<- names(dat.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')


# Export results
# write.csv(dat_out, "Sensitivity_results_5bins_equal.csv", row.names = F)
# write.csv(all.brkpts, "Sensitivity_allbreakpts_5bins_equal.csv", row.names = F)







#####################################
### Assess 10 bins of equal width ###
#####################################

set.seed(1)

#define bin number and limits for step lengths
hist(dat$dist, 100)
dist.bin.lims = c(seq(0, 40, length.out = 10), max(dat$dist, na.rm = T))  #10 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=dist)) +
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = dist.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nStep Length", y = "Density\n")



#assign bins to obs
behav.list1<- map(behav.list, discrete_move_var, lims = list(dist.bin.lims, angle.bin.lims),
                 varIn = c("dist", "rel.angle"), varOut = c("SL", "TA"))
behav.list2<- lapply(behav.list1, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment




## Run RJMCMC

#prior
alpha = 1

ngibbs = 60000

# handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

# track_length == 5000; takes 10 min for 60000 iterations
dat.res_5k<- segment_behavior(data = behav.list2[1:5], ngibbs = ngibbs, nbins = c(10,8),
                              alpha = alpha)

# track_length == 10000; takes 24 min for 60000 iterations
dat.res_10k<- segment_behavior(data = behav.list2[6:10], ngibbs = ngibbs, nbins = c(10,8),
                               alpha = alpha)

plan(sequential)  #closes background workers


## If sims analyzed separately, merge all runs together in single list
dat.res<- mapply(rbind, dat.res_5k, dat.res_10k, SIMPLIFY=FALSE)


## Reclassify and restructure data as needed
dat.res$nbrks[,2:ncol(dat.res$nbrks)]<- apply(dat.res$nbrks[,2:ncol(dat.res$nbrks)], 2,
                                              function(x) as.numeric(as.character(x)))
dat.res$LML[,2:ncol(dat.res$LML)]<- apply(dat.res$LML[,2:ncol(dat.res$LML)], 2,
                                          function(x) as.numeric(as.character(x)))
dat.res$brkpts<- dat.res$brkpts %>% 
  split(., row(.)) %>% 
  flatten()  #deals with dimensional list
names(dat.res$brkpts)<- dat.res$nbrks[,1]


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")


##Determine maximum a posteriori (MAP) estimate for selecting breakpoints
MAP.est<- get_MAP(dat.res$LML, nburn = ngibbs/2)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

## Visualize breakpoints over data streams
plot_breakpoints(data = behav.list1, as_date = FALSE, var_names = c("dist","rel.angle"),
                 var_labels = c("Step Length","Turning Angle"), brkpts = brkpts)


# Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:nrow(brkpts)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = brkpts[i,-1], true.brkpts = true.brkpts[i+5,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}


#assign time seg and make as DF
dat_out<- assign_tseg(dat = behav.list1, brkpts = brkpts)

names(all.brkpts)<- names(dat.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')


# Export results
# write.csv(dat_out, "Sensitivity_results_10bins_equal.csv", row.names = F)
# write.csv(all.brkpts, "Sensitivity_allbreakpts_10bins_equal.csv", row.names = F)






######################################
### Assess 10 bins using quantiles ###
######################################

set.seed(1)

#define bin number and limits for step lengths
hist(dat$dist, 100)
dist.bin.lims = quantile(dat$dist, seq(0, 1, length.out = 11), na.rm = T)  #10 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=dist)) +
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = dist.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nStep Length", y = "Density\n")



#assign bins to obs
behav.list1<- map(behav.list, discrete_move_var, lims = list(dist.bin.lims, angle.bin.lims),
                  varIn = c("dist", "rel.angle"), varOut = c("SL", "TA"))
behav.list2<- lapply(behav.list1, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment




## Run RJMCMC

#prior
alpha = 1

ngibbs = 40000

# handlers(handler_progress(clear = FALSE))
plan(multisession, workers = 5)

# track_length == 5000; takes 13.5 min for 40000 iterations
dat.res_5k<- segment_behavior(data = behav.list2[1:5], ngibbs = ngibbs, nbins = c(10,8),
                              alpha = alpha)

# track_length == 10000; takes 26 min for 40000 iterations
dat.res_10k<- segment_behavior(data = behav.list2[6:10], ngibbs = ngibbs, nbins = c(10,8),
                               alpha = alpha)

plan(sequential)  #closes background workers


## If sims analyzed separately, merge all runs together in single list
dat.res<- mapply(rbind, dat.res_5k, dat.res_10k, SIMPLIFY=FALSE)


## Reclassify and restructure data as needed
dat.res$nbrks[,2:ncol(dat.res$nbrks)]<- apply(dat.res$nbrks[,2:ncol(dat.res$nbrks)], 2,
                                              function(x) as.numeric(as.character(x)))
dat.res$LML[,2:ncol(dat.res$LML)]<- apply(dat.res$LML[,2:ncol(dat.res$LML)], 2,
                                          function(x) as.numeric(as.character(x)))
dat.res$brkpts<- dat.res$brkpts %>% 
  split(., row(.)) %>% 
  flatten()  #deals with dimensional list
names(dat.res$brkpts)<- dat.res$nbrks[,1]


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")


##Determine maximum a posteriori (MAP) estimate for selecting breakpoints
MAP.est<- get_MAP(dat.res$LML, nburn = ngibbs/2)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

## Visualize breakpoints over data streams
plot_breakpoints(data = behav.list1, as_date = FALSE, var_names = c("dist","rel.angle"),
                 var_labels = c("Step Length","Turning Angle"), brkpts = brkpts)


# Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:nrow(brkpts)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = brkpts[i,-1], true.brkpts = true.brkpts[i+5,-1],
                                   acc.tol = 10, dup.tol = 1, miss.tol = 30)
}


#assign time seg and make as DF
dat_out<- assign_tseg(dat = behav.list1, brkpts = brkpts)

names(all.brkpts)<- names(dat.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')


# Export results
# write.csv(dat_out, "Sensitivity_results_10bins_quantile.csv", row.names = F)
# write.csv(all.brkpts, "Sensitivity_allbreakpts_10bins_quantile.csv", row.names = F)




#########################################################
### Create Figure Showing Comparison of Bin Locations ###
#########################################################

equal.5bins<- c(seq(0, 40, length.out = 5), max(dat$dist, na.rm = T)) %>%   #5 bins
  data.frame(lims = .) %>% 
  mutate(method = "Equal width - 5 bins")
equal.10bins<- c(seq(0, 40, length.out = 10), max(dat$dist, na.rm = T)) %>%   #10 bins
  data.frame(lims = .) %>% 
  mutate(method = "Equal width - 10 bins")
quant.5bins<- quantile(dat$dist, c(0, 0.25, 0.50, 0.75, 0.90, 1), na.rm = T) %>%   #5 bins
  data.frame(lims = .) %>% 
  mutate(method = "Quantiles - 5 bins")
quant.10bins<- quantile(dat$dist, seq(0, 1, length.out = 11), na.rm = T) %>%   #10 bins
  data.frame(lims = .) %>% 
  mutate(method = "Quantiles - 10 bins")

all.lims<- rbind(equal.5bins, equal.10bins, quant.5bins, quant.10bins)
all.lims$method<- factor(all.lims$method, levels = unique(all.lims$method))



ggplot() +
  geom_density(data = behav.df, aes(x=dist), fill = "lightblue") +
  geom_vline(data = all.lims, aes(xintercept = lims), linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold")) +
  labs(x = "\nStep Length", y = "Density\n") +
  facet_wrap(~ method)

#Export figure
setwd("~/Documents/Manuscripts/Bayesian Behavior Estimation Model/Figures")

# ggsave('Bin Comparison for Step Lengths.png', width = 9, height = 5, units = "in", dpi = 330)
