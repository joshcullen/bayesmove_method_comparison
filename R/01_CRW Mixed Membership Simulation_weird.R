############################################
#### Simulations w/ Weird Distributions ####
############################################

library(tidyverse)
library(tictoc)
library(cowplot)

source('R/Simulation Functions.R')
source('R/helper functions.R')

### Simulate full track w/ different numbers of time segments ###
## Create tracks w/ 10, 50, 100, and 500 segments while keeping 100 obs. per tseg.
## Generate 5 different versions of each of these tracks (all else being equal)

#define behaviors and sample them from Categorical distribution
#weight probabilities so that behavior 1 (Resting) occurs 50%, behavior 2
#(Area-restricted search) occurs 30%, and behavior 3 (Transit) occurs 20%

set.seed(2)


#simulate track
ntseg<- c(10, 50, 100, 500)
nstep<- 100
SL.params<- data.frame(lo = c(0, 0, 0),
                       hi = c(Inf, Inf, Inf),
                       mu=c(0.25, 2, exp(2)),
                       sig = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0),
                       par2 = c(0.5, pi, 1))


tic()
track.sim<- CRW.sim2(nsim=5, ntseg = ntseg, nstep = nstep, SL.params = SL.params,
                    TA.params = TA.params, Z0=c(0,0))
toc()
#takes ~36 s to run


#extract tracks
tracks<- track.sim$tracks %>%
  modify_depth(2, ~modify_at(., "id", as.character)) %>%
  modify_depth(1, ~map_dfr(., `[`)) %>%
  map_dfr(`[`)

#extract breakpoints
max.length<- track.sim$brkpts %>%  #to set max number of columns for DF
  modify_depth(2, ~length(.)) %>%
  unlist() %>%
  max()

brkpts<- track.sim$brkpts %>%
  modify_depth(2, function(x) {c(x, rep(NA, max.length-length(x)))}) %>%
  modify_depth(1, ~map_dfr(., `[`)) %>%
  map(t) %>%
  map(as.data.frame) %>% 
  map_dfr(., `[`)
brkpts<- cbind(id = unique(tracks$id), brkpts)
names(brkpts)<- c('id', paste0("Brk_",1:(ncol(brkpts)-1)))


## Plot tracks ##

ggplot(data = tracks %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_coarse), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




## Comparison of distributions based on observation-level behavior classification
ggplot(tracks[tracks$track_length == 1000,], aes(SL, color=behav_fine)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(SL, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(SL, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(SL, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(tracks[tracks$track_length == 1000,], aes(TA, color=behav_fine)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(TA, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(TA, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(TA, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()



## Comparison of distributions based on segment-level behavior classification
ggplot(tracks[tracks$track_length == 1000,], aes(SL, color=behav_coarse)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(SL, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(SL, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(SL, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(tracks[tracks$track_length == 1000,], aes(TA, color=behav_coarse)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(TA, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(TA, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(TA, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()




##############
## Figure 2 ##
##############

SL.params<- data.frame(par1 = c(0.25, 2, exp(2)), par2 = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0), par2 = c(0.5, pi, 1))




### Step Lengths

SL.params2<- data.frame(SL.params, behavior = c("Encamped","ARS","Transit"))


true_plot_data_SL<- 
  pmap_df(SL.params2,
          function(par1, par2, behavior) {
            tibble(x = seq(0, 40, by = 0.025),
                   y = dtnorm(x, mean1 = par1, sd1 = par2, lo = 0, hi = Inf),
                   behavior = behavior)
          })
true_plot_data_SL$behavior<- factor(true_plot_data_SL$behavior,
                                    levels = c("Encamped","ARS","Transit"))


SL.plot<- ggplot(data = true_plot_data_SL, aes(x = x, y = y)) +
  geom_area(aes(color = behavior, fill = behavior), size = 0.9, alpha  = 0.5,
            position = "dodge") +
  labs(x = "\nStep Length (units)", y = "Density\n") +
  scale_color_viridis_d("") +
  scale_fill_viridis_d("") +
  scale_y_continuous(breaks = c(0, 0.5, 1.00)) +
  annotate(geom = "rect", xmin = 3, xmax = 19, ymin = 0.53, ymax = 0.57, color = "#440154FF",
           fill = "#440154FF", alpha = 0.1) +
  annotate(geom = "text", x = 11, y = 0.55, label = 'TN(0.25, 1, 0, \U221E)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = 7, xmax = 21, ymin = 0.18, ymax = 0.22, color = "#21908CFF",
           fill = "#21908CFF", alpha = 0.1) +
  annotate(geom = "text", x = 14, y = 0.2, label = 'TN(2, 2, 0, \U221E)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = 18, xmax = 32, ymin = 0.08, ymax = 0.12, color = "#FDE725FF",
           fill = "#FDE725FF", alpha = 0.1) +
  annotate(geom = "text", x = 25, y = 0.1, label = 'TN(7, 20, 0, \U221E)',
           fontface = "italic") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75)))








### Turning Angles

TA.params2<- data.frame(TA.params, behavior = c("Encamped","ARS","Transit"))


true_plot_data_TA<- list()

true_plot_data_TA[[1]]<- 
  pmap_df(TA.params2[1,],
          function(par1, par2, behavior) {
            tibble(x = seq(0.001, 0.999, length.out = 513),
                   y = dbeta(x, shape1 = par1, shape2 = par2),
                   behavior = behavior)
          })
true_plot_data_TA[[1]]$x<- (true_plot_data_TA[[1]]$x * (2*pi)) - pi
true_plot_data_TA[[1]]$y<- true_plot_data_TA[[1]]$y / max(true_plot_data_TA[[1]]$y)

true_plot_data_TA[[2]]<- 
  pmap_df(TA.params2[2,],
          function(par1, par2, behavior) {
            tibble(x = seq(-pi, pi, by = (2*pi/512)),
                   y = dunif(x, min = par1, max = par2),
                   behavior = behavior)
          })
true_plot_data_TA[[3]]<- 
  pmap_df(TA.params2[3,],
          function(par1, par2, behavior) {
            tibble(x = seq(-pi, pi, by = (2*pi/512)),
                   y = dtnorm(x, mean1 = par1, sd1 = par2, lo = -pi, hi = pi),
                   behavior = behavior)
          })

true_plot_data_TA<- bind_rows(true_plot_data_TA)
true_plot_data_TA$behavior<- factor(true_plot_data_TA$behavior,
                                    levels = c("Encamped","ARS","Transit"))


TA.plot<- ggplot(data = true_plot_data_TA, aes(x = x, y = y)) +
  geom_area(aes(color = behavior, fill = behavior), size = 0.9, alpha  = 0.5,
            position = "dodge") +
  labs(x = "\nTurning Angle (radians)", y = "Density\n") +
  scale_color_viridis_d("") +
  scale_fill_viridis_d("") +
  scale_y_continuous(breaks = c(0, 0.5, 1.00)) +
  annotate(geom = "rect", xmin = -2.9, xmax = 0.5, ymin = 0.82, ymax = 0.88,
           color = "#440154FF", fill = "#440154FF", alpha = 0.1) +
  annotate(geom = "text", x = -1.2, y = 0.85, label = 'Beta(0.5, 0.5) \U00D7 2\U03C0 - \U03C0',
           fontface = "italic") +
  annotate(geom = "rect", xmin = -2.6, xmax = -1, ymin = 0.34, ymax = 0.4,
           color = "#21908CFF", fill = "#21908CFF", alpha = 0.1) +
  annotate(geom = "text", x = -1.8, y = 0.37, label = 'Unif(-\U03C0, \U03C0)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = 0.2, xmax = 2.6, ymin = 0.47, ymax = 0.53, color = "#FDE725FF",
           fill = "#FDE725FF", alpha = 0.1) +
  annotate(geom = "text", x = 1.4, y = 0.5, label = 'TN(7, 20, -\U03C0, \U03C0)',
           fontface = "italic") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75)))





## Create composite plot w/ shared legend
p.comp<- plot_grid(SL.plot + theme(legend.position="none"),
                   NA,
                   TA.plot + theme(legend.position="none"),
                   align = 'v',
                   rel_widths = c(1, 0.1, 1),
                   hjust = -1,
                   nrow = 1)

# extract the legend from one of the plots
legend.comp<- get_legend(SL.plot + theme(legend.position="top"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(legend.comp, p.comp, ncol = 1, rel_heights = c(0.1, 1))

# ggsave("Figure 2_weird.png", width = 9, height = 5, units = "in", dpi = 330)







# write.csv(tracks, "data/CRW_MM_sim_weird.csv", row.names = F)
# write.csv(brkpts, "data/CRW_MM_sim_brkpts_weird.csv", row.names = F)

