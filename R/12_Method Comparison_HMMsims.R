#########################
### Method Comparison ###
#########################

library(tidyverse)
library(wesanderson)
library(lubridate)
library(cowplot)
library(viridis)
library(ggnewscale)
library(circular)
library(gridExtra)

source('R/helper functions.R')



###########################################################################
### Compare HMM-based simulations generated from uncommon distributions ###
###########################################################################

# Load results
bayes.res_weird<- read.csv("data/Modeled HMM Sim Tracks w Behav_weird.csv")
hmm.res_weird<- read.csv("data/HMM(HMM) results_weird.csv")
segclust.res_weird<- read.csv("data/Segclust2d(HMM) results_weird.csv")
embc.res_weird<- read.csv("data/EMbC(HMM) results_weird.csv")


## Define color palette
pal1<- c(wes_palette("Darjeeling1", 5)[-4], "mediumorchid")


# Assign identifiers by method and make consistent behavior colname
bayes.res_weird$method<- "M4"
hmm.res_weird$method<- "HMM"
segclust.res_weird$method<- "Segclust2d"
embc.res_weird$method<- "EMbC"

bayes.res_weird<- bayes.res_weird %>% 
  mutate_at("behav", ~factor(., levels = 1:3))
hmm.res_weird<- hmm.res_weird %>% 
  rename(id = ID)




### Compare Step Length and Turning Angle Distributions ###


## Figure S11.1 ##

### Define bin limits
dat_weird<- read.csv("data/HMM_sim_weird.csv", as.is = T)
dat_weird$dt<- 3600
names(dat_weird)[5:6]<- c("dist","rel.angle")

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

dist.bin.lims=quantile(dat_weird[dat_weird$dt == 3600,]$dist,
                       c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins




hmm.SL.params_weird<- read.csv("data/HMM(HMM) result step params_weird.csv", as.is = T)  
hmm.TA.params_weird<- read.csv("data/HMM(HMM) result angle params_weird.csv", as.is = T)  

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



segclust.params_weird<- read.csv("data/Segclust(HMM) result state params_weird.csv", as.is = T) %>% 
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
embc.params_weird<- read.csv("data/EMbC(HMM) result state params_weird.csv", as.is = T) %>% 
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
behav.res_weird<-  read.csv("data/HMM LDA Phi values_weird.csv", as.is = T)
behav.res_weird<- bayesmove::df_to_list(behav.res_weird, "id")

behav.order_weird<- read.csv("data/HMM LDA behavior order_weird.csv", as.is = T)
behav.order_weird<- bayesmove::df_to_list(behav.order_weird, "id")
behav.order_weird<- map(behav.order_weird, {. %>% 
    dplyr::select(-id) %>% 
    as.numeric()}
    )


#calculate true proportions of SL and TA by behavior for all bins for ID 5_2 for Bayesian model
bayes.b_weird<- behav.res_weird[[10]]
bayes.b_weird$behav<- bayes.b_weird$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "ARS") %>% 
  str_replace_all(., "2", "Transit") %>%
  str_replace_all(., "3", "Encamped") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.distfit<- ggplot(true.b_weird, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = hmm.b_weird[[10]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[2], color="black", stroke=1) +
  geom_point(data = segclust.b_weird[[10]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[3], color="black", stroke=1) +
  geom_point(data = embc.b_weird[[10]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
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
  bayes.b_weird<- behav.res_weird[[i]]
  
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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))




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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))





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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))





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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))




rmse.df_weird<- data.frame(id = rep(rep(names(behav.res_weird), 4), each = 2),
                           track_length = factor(rep(rep(rep(c(1000,5000),
                                                             each = 5), 4), each = 2),
                                                 levels = c("1000","5000")),
                           rmse = rbind(bayes.rmse_weird, hmm.rmse_weird, segclust.rmse_weird,
                                        embc.rmse_weird),
                           method = rep(c("M4","HMM","Segclust2d","EMbC"), each = 20))
# rmse.df_weird<- rbind(rmse.df_weird, segclust.rmse.df)
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
# png("Figure S11.1 (rmse from HMM sim).png", width = 14.5, height = 5.5, units = "in", res = 330)

grid.arrange(p.distfit, p.rmse_weird, heights = c(0.2, 1),
             widths = c(1, 0.2, 1, 0.5),
             layout_matrix = rbind(c(NA, NA, NA, NA),
                                   c(3, NA, 4, 4)))
# dev.off()











#########################################################################
### Compare HMM-based simulations generated from common distributions ###
#########################################################################

# Load results
bayes.res_parametric<- read.csv("data/Modeled HMM Sim Tracks w Behav_parametric.csv")
hmm.res_parametric<- read.csv("data/HMM(HMM) results_parametric.csv")
segclust.res_parametric<- read.csv("data/Segclust2d(HMM) results_parametric.csv")
embc.res_parametric<- read.csv("data/EMbC(HMM) results_parametric.csv")


## Define color palette
pal1<- c(wes_palette("Darjeeling1", 5)[-4], "mediumorchid")


# Assign identifiers by method and make consistent behavior colname
bayes.res_parametric$method<- "M4"
hmm.res_parametric$method<- "HMM"
segclust.res_parametric$method<- "Segclust2d"
embc.res_parametric$method<- "EMbC"

bayes.res_parametric<- bayes.res_parametric %>% 
  mutate_at("behav", ~factor(., levels = 1:3))
hmm.res_parametric<- hmm.res_parametric %>% 
  rename(id = ID)




### Compare Step Length and Turning Angle Distributions ###


## Figure S11.2 ##

### Define bin limits
dat_parametric<- read.csv("data/HMM_sim.csv", as.is = T)
dat_parametric$dt<- 3600
names(dat_parametric)[5:6]<- c("dist","rel.angle")

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

dist.bin.lims=quantile(dat_parametric[dat_parametric$dt == 3600,]$dist,
                       c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins



hmm.SL.params_parametric<- read.csv("data/HMM(HMM) result step params_parametric.csv", as.is = T)  
hmm.TA.params_parametric<- read.csv("data/HMM(HMM) result angle params_parametric.csv", as.is = T)  

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



segclust.params_parametric<- read.csv("data/Segclust(HMM) result state params_parametric.csv", as.is = T) %>% 
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
embc.params_parametric<- read.csv("data/EMbC(HMM) result state params_parametric.csv", as.is = T) %>% 
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





SL.params_parametric<- data.frame(par1 = c(0.25, 2, 10), par2 = c(1, 1, 1))
TA.params<- data.frame(par1 = c(pi, pi, 0), par2 = c(0.8, 0, 0.8))

true.b_parametric<- extract.behav.props(params = list(SL.params_parametric, TA.params),
                                         lims = list(dist.bin.lims, angle.bin.lims),
                                         behav.names = c("Encamped","ARS","Transit"))


hmm.b_parametric<- map2(hmm.SL.params2_parametric, hmm.TA.params2_parametric,
                   ~extract.behav.props(params = list(.x, .y),
                                        lims = list(dist.bin.lims, angle.bin.lims),
                                        behav.names = c("Encamped","ARS","Transit"))
)


segclust.b_parametric<- map(segclust.params2_parametric,
                       ~extract.behav.props_norm(params = .x,
                                                 lims = list(dist.bin.lims, angle.bin.lims),
                                                 behav.names = c("Encamped","ARS","Transit"))
)


embc.b_parametric<- map(embc.params2_parametric,
                   ~extract.behav.props_norm(params = .x,
                                             lims = list(dist.bin.lims, angle.bin.lims),
                                             behav.names = c("Encamped","ARS","Transit"))
)




## Bayesian
behav.res_parametric<-  read.csv("data/HMM LDA Phi values_parametric.csv", as.is = T)
behav.res_parametric<- bayesmove::df_to_list(behav.res_parametric, "id")

behav.order_parametric<- read.csv("data/HMM LDA behavior order_parametric.csv", as.is = T)
behav.order_parametric<- bayesmove::df_to_list(behav.order_parametric, "id")
behav.order_parametric<- map(behav.order_parametric, {. %>% 
    dplyr::select(-id) %>% 
    as.numeric()}
)


#calculate true proportions of SL and TA by behavior for all bins for ID 5_2 for Bayesian model
bayes.b_parametric<- behav.res_parametric[[10]]
bayes.b_parametric$behav<- bayes.b_parametric$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "Transit") %>%
  str_replace_all(., "3", "ARS") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.distfit<- ggplot(true.b_parametric, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = hmm.b_parametric[[10]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[2], color="black", stroke=1) +
  geom_point(data = segclust.b_parametric[[10]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill=pal1[3], color="black", stroke=1) +
  geom_point(data = embc.b_parametric[[10]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))




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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))





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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))





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
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 10))




rmse.df_parametric<- data.frame(id = rep(rep(names(behav.res_parametric), 4), each = 2),
                           track_length = factor(rep(rep(rep(c(1000,5000),
                                                             each = 5), 4), each = 2),
                                                 levels = c("1000","5000")),
                           rmse = rbind(bayes.rmse_parametric, hmm.rmse_parametric, segclust.rmse_parametric,
                                        embc.rmse_parametric),
                           method = rep(c("M4","HMM","Segclust2d","EMbC"), each = 20))
# rmse.df_parametric<- rbind(rmse.df_parametric, segclust.rmse.df)
rmse.df_parametric$method<- factor(rmse.df_parametric$method, levels = c('M4','HMM','Segclust2d','EMbC'))


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
# png("Figure S11.2 (rmse from HMM sim).png", width = 14.5, height = 5.5, units = "in", res = 330)

grid.arrange(p.distfit, p.rmse_parametric, heights = c(0.2, 1),
             widths = c(1, 0.2, 1, 0.5),
             layout_matrix = rbind(c(NA, NA, NA, NA),
                                   c(3, NA, 4, 4)))
# dev.off()

