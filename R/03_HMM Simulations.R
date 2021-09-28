##############################################################
#### HMM Simulations w/ Common and Uncommon Distributions ####
##############################################################

library(tidyverse)
library(circular)
library(tictoc)

source('R/Simulation Functions.R')
source('R/helper functions.R')


### Simulate full track w/ different lengths ###
## Create tracks w/ 1000 and 5000 observations (which is similar in length to many datasets currently being analyzed)
## Generate 5 different versions of each of these tracks (all else being equal)

## Common distributions
#step lengths are drawn from a gamma dist; turning angles are drawn from a wrapped Cauchy dist

## Uncommon distributions
#step lengths are drawn from a truncated normal dist; turning angles are drawn from either a beta, uniform, or truncated normal dist

## Sample latent states based on transition probability matrix (each state has 0.9 probability of remaining within state from time t to t+1)
# This is within range of what has been found in previous studies (e.g., Hance et al 2021; Jonsen et al 2005; McClintock et al 2020; McClintock and Michelot, 2018)



################################################
### Create simulations using common distribs ###
################################################

set.seed(2)

#simulate track
nobs<- rep(c(1000,5000), each = 5)
SL.params<- data.frame(shape=c(0.25, 2, 10), scale = c(1, 1, 1))
TA.params<- data.frame(mu=c(pi, pi, 0), rho = c(0.8, 1e-12, 0.8))  #error if rho == 0


tic()
tracks<- HMM.sim(nsim = 10,
                 nobs = nobs,
                 SL.params = SL.params,
                 TA.params = TA.params,
                 Z0 = c(0,0))
toc()
#takes ~5 sec to run


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
  facet_wrap(track_length ~ id, scales = "free")



## Compare empirical state-dependent distributions of SL and TA ##

ggplot(data = tracks.df[tracks.df$track_length == 1000,], aes(SL, color=factor(state))) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks.df[tracks.df$track_length == 5000,], aes(SL, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(data = tracks.df[tracks.df$track_length == 1000,], aes(TA, color=factor(state))) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks.df[tracks.df$track_length == 5000,], aes(TA, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()






##################################################
### Create simulations using uncommon distribs ###
##################################################

set.seed(2)

#simulate track
nobs<- rep(c(1000,5000), each = 5)
SL.params<- data.frame(lo = c(0, 0, 0),
                       hi = c(Inf, Inf, Inf),
                       mu = c(0.25, 2, exp(2)),
                       sig = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0),
                       par2 = c(0.5, pi, 1))


tic()
tracks2<- HMM.sim2(nsim = 10,
                 nobs = nobs,
                 SL.params = SL.params,
                 TA.params = TA.params,
                 Z0 = c(0,0))
toc()
#takes ~3 sec to run


tracks.df2<- bind_rows(tracks2) %>% 
  mutate(rep = str_sub(id, start = 1, end = 1))


# Plot time series of step length
ggplot(tracks.df2, aes(time1, SL)) +
  geom_line() +
  geom_point(aes(color = factor(state))) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_grid(rep ~ track_length, scales = "free")

# Plot time series of turning angle
ggplot(tracks.df2, aes(time1, TA)) +
  geom_line() +
  geom_point(aes(color = factor(state))) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_grid(rep ~ track_length, scales = "free")

# Plot tracks
ggplot(tracks.df2, aes(x, y)) +
  geom_path() +
  geom_point(aes(color = factor(state))) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_wrap(track_length ~ id, scales = "free")



## Compare empirical state-dependent distributions of SL and TA ##

ggplot(data = tracks.df2[tracks.df2$track_length == 1000,], aes(SL, color=factor(state))) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks.df2[tracks.df2$track_length == 5000,], aes(SL, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(data = tracks.df2[tracks.df2$track_length == 1000,], aes(TA, color=factor(state))) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks.df2[tracks.df2$track_length == 5000,], aes(TA, color=factor(state)),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()





# Export simulations

# write.csv(tracks.df, "data/HMM_sim.csv", row.names = F)
# write.csv(tracks.df2, "data/HMM_sim_weird.csv", row.names = F)
