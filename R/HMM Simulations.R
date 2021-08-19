#################################################
#### HMM Simulations w/ Common Distributions ####
#################################################

library(tidyverse)
library(circular)
library(tictoc)

source('R/Simulation Functions.R')

### Simulate full track w/ different lengths ###
## Create tracks w/ 1000, 5000, 10000, and 50000 observations
## Generate 5 different versions of each of these tracks (all else being equal)

#step lengths are drawn from a gamma dist; turning angles are drawn from a wrapped Cauchy dist
#sample latent states based on transition probability matrix (each state has 0.9 probability of remaining within state from time t to t+1)

set.seed(2)


#simulate track
nobs<- rep(c(1000,5000,10000,50000), each = 5)
SL.params<- data.frame(shape=c(0.25, 2, 10), scale = c(1, 1, 1))
TA.params<- data.frame(mu=c(pi, pi, 0), rho = c(0.8, 1e-12, 0.8))  #error if rho == 0


tic()
tracks<- HMM.sim(nsim = 20,
                 nobs = nobs,
                 SL.params = SL.params,
                 TA.params = TA.params,
                 Z0 = c(0,0))
toc()
#takes ~48 sec to run


tracks.df<- bind_rows(tracks) %>% 
  mutate(rep = str_sub(id, start = 1, end = 1))


# Plot time series of step length
ggplot(tracks.df, aes(time1, SL)) +
  geom_line() +
  geom_point(aes(color = factor(state))) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_grid(rep ~ track_length, scales = "free")

# Plot time series of turning angle
ggplot(tracks.df, aes(time1, TA)) +
  geom_line() +
  geom_point(aes(color = factor(state))) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_grid(rep ~ track_length, scales = "free")

# Plot tracks
ggplot(tracks.df, aes(x, y)) +
  geom_path() +
  geom_point(aes(color = factor(state))) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_wrap(~ id, scales = "free")



## Compare empirical state-dependent distributions of SL and TA ##

ggplot(data = tracks.df[tracks.df$track_length == 1000,], aes(SL, color=factor(state))) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks.df[tracks.df$track_length == 5000,], aes(SL, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks.df[tracks.df$track_length == 10000,], aes(SL, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks.df[tracks.df$track_length == 50000,], aes(SL, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(data = tracks.df[tracks.df$track_length == 1000,], aes(TA, color=factor(state))) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks.df[tracks.df$track_length == 5000,], aes(TA, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks.df[tracks.df$track_length == 10000,], aes(TA, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks.df[tracks.df$track_length == 50000,], aes(TA, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


# Export simulations

# write.csv(tracks.df, "data/HMM_sim.csv", row.names = F)
