#####################################
#####################################
#### Mixed-Membership Simulation ####
#####################################
#####################################

library(tidyverse)
library(circular)
library(tictoc)

source('Simulation Functions.R')

### Simulate full track w/ different numbers of time segments ###
## Create tracks w/ 10, 50, 100, and 500 segments while keeping 100 obs per tseg
## Generate 5 different versions of each of these tracks (all else being equal)

#define behaviors and sample them from Categorical distribution
#weight probs so that behavior 1 (Resting) occurs 50%, behavior 2 (Area-restricted search) occurs 30%, and behavior 3 (Transit) occurs 20%

set.seed(2)


#simulate track
ntseg<- c(10, 50, 100, 500)
nstep<- 100
SL.params<- data.frame(shape=c(0.25, 2, 10), scale = c(1, 1, 1))
TA.params<- data.frame(mu=c(pi, pi, 0), rho = c(0.8, 0, 0.8))


tic()
track.sim<- CRW.sim(nsim=5, ntseg = ntseg, nstep = nstep, SL.params = SL.params,
                    TA.params = TA.params, Z0=c(0,0))
toc()
#takes ~1.5 min to run


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
  # coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




## Compare distributions of SL and TA among fine-scale behaviors ##

# behav_fine
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


# behav_coarse
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






# write.csv(tracks, "CRW_MM_sim_parametric.csv", row.names = F)
# write.csv(brkpts, "CRW_MM_sim_brkpts_parametric.csv", row.names = F)





##############
## Figure 2 ##
##############

SL.params<- data.frame(par1 = c(0.25, 2, 10), par2 = c(1, 1, 1))
TA.params<- data.frame(par1 = c(pi, pi, 0), par2 = c(0.8, 0, 0.8))




### Step Lengths

SL.params2<- data.frame(SL.params, behavior = c("Encamped","ARS","Transit"))


true_plot_data_SL<- 
  pmap_df(SL.params2,
          function(par1, par2, behavior) {
            tibble(x = seq(0, 40, by = 0.025),
                   y = dgamma(x, shape = par1, scale = par2),
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
  # scale_y_continuous(breaks = c(0, 0.5, 1.00)) +
  annotate(geom = "rect", xmin = 0.5, xmax = 15.5, ymin = 3.35, ymax = 3.65, color = "#440154FF",
           fill = "#440154FF", alpha = 0.1) +
  annotate(geom = "text", x = 8, y = 3.5, label = 'Gamma(0.25, 1)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = 4, xmax = 16, ymin = 1.15, ymax = 1.45, color = "#21908CFF",
           fill = "#21908CFF", alpha = 0.1) +
  annotate(geom = "text", x = 10, y = 1.3, label = 'Gamma(2, 1)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = 8, xmax = 22, ymin = 0.35, ymax = 0.65, color = "#FDE725FF",
           fill = "#FDE725FF", alpha = 0.1) +
  annotate(geom = "text", x = 15, y = 0.5, label = 'Gamma(10, 1)',
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

true_plot_data_TA<- 
  pmap_df(TA.params2,
          function(par1, par2, behavior) {
            tibble(x = seq(-pi, pi, by = (pi/512)),
                   y = dwrappedcauchy(x, mu = circular(par1), rho = par2),
                   behavior = behavior)
          })
true_plot_data_TA$behavior<- factor(true_plot_data_TA$behavior,
                                    levels = c("Encamped","ARS","Transit"))


TA.plot<- ggplot(data = true_plot_data_TA, aes(x = x, y = y)) +
  geom_area(aes(color = behavior, fill = behavior), size = 0.9, alpha  = 0.5,
            position = "dodge") +
  labs(x = "\nTurning Angle (radians)", y = "Density\n") +
  scale_color_viridis_d("") +
  scale_fill_viridis_d("") +
  # scale_y_continuous(breaks = c(0, 0.5, 1.00)) +
  annotate(geom = "rect", xmin = -2.9, xmax = -1.3, ymin = 1.25, ymax = 1.35,
           color = "#440154FF", fill = "#440154FF", alpha = 0.1) +
  annotate(geom = "text", x = -2.1, y = 1.3, label = 'WC(\U03C0, 0.8)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = 0.4, xmax = 2, ymin = 0.75, ymax = 0.85,
           color = "#FDE725FF", fill = "#FDE725FF", alpha = 0.1) +
  annotate(geom = "text", x = 1.2, y = 0.8, label = 'WC(0, 0.8)',
           fontface = "italic") +
  annotate(geom = "rect", xmin = -2.4, xmax = -1, ymin = 0.3, ymax = 0.4, color = "#21908CFF",
           fill = "#21908CFF", alpha = 0.1) +
  annotate(geom = "text", x = -1.7, y = 0.35, label = 'WC(\U03C0, 0)',
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
plot_grid(SL.plot + theme(legend.position="none"),
                   NA,
                   TA.plot + theme(legend.position="none"),
                   rel_widths = c(1, 0.1, 1),
                   hjust = -1,
                   nrow = 1)

# ggsave("Figure 2_parametric.png", width = 9, height = 5, units = "in", dpi = 330)
