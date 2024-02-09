#load packages
library(tidyverse)
library(vegan)
library(cowplot)
library(brms)
library(car)
library(jtools)
library(tidybayes)
library(emmeans)
library(posterior)

### SET WORKING DIRECTORY ###
setwd("")

# load data
NSD <- read.csv("damselfish_NSD_full.csv", header = T, stringsAsFactors = F) # collated behavioural outputs from BORIS v.8.6.2 (Friard & Gamba, 2016) 
NSD$type <- as.factor(NSD$type)
NSD$size_bin <- as.factor(NSD$size_bin) # -1-smaller, 1-larger
NSD <- NSD %>%
  mutate(total_aggression = bite + display)

# Bayes analysis
subset <- NSD %>%
  filter(type != "control")

# aggression
# setting priors
aggression_priors <- brms::set_prior("normal(0,10)", class = "b")

# run model
aggression_mod <- brm(total_aggression ~ type * size_bin + (1|damselID), 
                        family = negbinomial(),
                        data = subset,
                        iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                        control = list(adapt_delta = 0.95), 
                        prior = aggression_priors,
                        sample_prior = "yes")
print(aggression_mod)

# check plots
pp_check(aggression_mod) # looks good
plot(aggression_mod, ask = FALSE)
bayes_R2(aggression_mod)

# extract estimates
# to report median estimates for each group, on response scale: 
aggression_mod %>% 
  emmeans(~ type * size_bin,
          epred = TRUE,
          re_formula = NA)

#type       size_bin emmean lower.HPD upper.HPD
#neighbour  -1        10.17     1.704      25.6
#neighbour  1          3.94     0.493      12.2
#stranger   -1         7.03     1.170      19.3
#stranger   1          5.27     0.832      13.9

# contrasts: full output for interaction effect
aggression_mod %>% 
  emmeans(c("type", "size_bin"),
          type = "response")  %>%
  contrast(method = "pairwise")

# averaged over size bin
aggression_mod %>% 
  emmeans(~ type,
          by = c("size_bin"),
          epred = TRUE,
          re_formula = NA)

# averaged over type
aggression_mod %>% 
  emmeans(~ size_bin,
          by = c("type"),
          epred = TRUE,
          type = "response")

# make plot
plot_theme <-
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(), 
        legend.position = c(.85,.85),
        rect = element_rect(colour = "white"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))

pred_aggression <-
  aggression_mod %>% 
  epred_draws(newdata = expand_grid(type = c("neighbour", "stranger"), 
                                    size_bin = c("-1", "1"),
                                    reps = 1000), 
              re_formula=NA)

aggression_plot <-
  ggplot(data = pred_aggression, aes(x = type, y = .epred, color = size_bin, fill = size_bin, shape = size_bin)) +
  geom_point(data = subset, aes(x = type, y = total_aggression, color = size_bin, fill = size_bin), color = "black", shape = 18,
             position = position_dodge(width = .3)) +
  stat_pointinterval(point_interval = median_hdi, .width = c(.9, .7), fatten_point = 2, 
                     interval_alpha = 0.5, position = position_dodge(width = .3)) +
  scale_color_manual(values = c("#0067A5", "#F99379"),
                     breaks = c("-1",  "1"),
                     labels = c("smaller",  "larger"),
                     name = "" ) +
  scale_fill_manual(values = c("#0067A5", "#F99379"),
                    breaks=c("-1",  "1"),
                    labels=c("smaller",  "larger"),
                    name = "" ) + 
  scale_shape_manual(values = c(19, 15),
                     breaks=c("-1",  "1"),
                     labels=c("smaller",  "larger"),
                     name = "" ) +
  scale_x_discrete(name = "Stimulus type", labels = c("Neighbour", "Stranger")) +
  scale_y_continuous(name = "Total aggressive displays") +
  plot_theme
aggression_plot # Figure S1

# posterior distribution
plot_theme <-
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = 'none',
        legend.background = element_rect(colour = "white"),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))

aggression_post <-
  pred_aggression %>% 
  mutate(Type = str_c(type, size_bin, sep =".")) %>%
  mutate(Type = as.factor(Type)) %>%
  mutate(Type = fct_relevel(Type, c("neighbour.1", "neighbour.-1", "stranger.1", "stranger.-1")))%>%
  ggplot(aes(x = .epred, y = Type, color = Type, fill = Type)) +
  stat_slab(aes(x = .epred, color = Type, fill = Type), point_interval = median_hdi, slab_alpha = .4) +
  stat_pointinterval(aes(x = .epred, color = Type, fill = Type),
                       .width=c(.9, .7), position = position_dodge(width = .3, preserve = "single"), shape = 24) +
  scale_color_manual(values = c("#0067A5", "#A1CAF1", "#F99379", "#BE0032"),
                     name = "",
                     breaks=c("neighbour.1", "neighbour.-1", "stranger.1", "stranger.-1"),
                     labels=c("Larger neighbour", "Smaller neighbour", "Larger stranger", "Smaller stranger"))+
  scale_fill_manual(values = c("#0067A5", "#A1CAF1", "#F99379", "#BE0032"),
                    name = "",
                    breaks=c("neighbour.1", "neighbour.-1", "stranger.1", "stranger.-1"),
                    labels=c("Larger neighbour", "Smaller neighbour", "Larger stranger", "Smaller stranger")) + #light blue = #A1CAF1, light red/pink = #F99379, orange = #F38400, purple = #875692
  scale_y_discrete(name = "Posterior density", labels = c("Larger neighbour", "Smaller neighbour", "Larger stranger", "Smaller stranger")) +
  scale_x_continuous(name = "Total aggressive displays", lim = c(0, 30)) +
  plot_theme
aggression_post # Figure 2


# vicinity
# setting priors
vicinity_priors <- set_prior("normal(0,10)", class = "b")

# run model
vicinity_mod <- brm(vicinity ~ type * size_bin + (1|damselID), 
                      family = gaussian(),
                      data = subset,
                      iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                      control = list(adapt_delta = 0.95), 
                      prior = vicinity_priors,
                      sample_prior = "yes")
print(vicinity_mod)

# check plots
pp_check(vicinity_mod) # looks good
plot(vicinity_mod, ask = FALSE)
bayes_R2(vicinity_mod)

# extract estimates
# to report median estimates for each group, on response scale:
emm_options(opt.digits = FALSE)
vicinity_mod %>% 
  emmeans(~ type * size_bin,
          epred = TRUE,
          re_formula = NA) 

#type       size_bin emmean lower.HPD upper.HPD
#neighbour  -1         97.0      65.9       125.1
#neighbour  1          94.4      63.9       124.5
#stranger   -1         92.6      63.4       123.8
#stranger   1          81.7      49.8       112.3

# make plot
plot_theme <-
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(.85, .85),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))

pred_vicinity <-
  vicinity_mod %>% 
  epred_draws(newdata = expand_grid(type = c("neighbour", "stranger"), 
                                    size_bin = c("-1", "1"),
                                    reps = 1000), 
              re_formula=NA)

vicinity_plot <-
  pred_vicinity %>%
  ggplot(aes(x = type, y = .epred, color = size_bin, fill = size_bin, shape = size_bin)) +
  geom_point(data = NSD, aes(x = factor(type, level = c("neighbour", "stranger", "control")), y = vicinity, color = size_bin, fill = size_bin), color = "black", shape = 18,
             position = position_dodge(width = .3)) +
  stat_pointinterval(point_interval = median_hdi, .width = c(.9, .7), fatten_point = 2, 
                     interval_alpha = 0.5, position = position_dodge(width = .3)) +
  scale_color_manual(values = c("#0067A5", "#F99379"),
                     breaks = c("-1",  "1"),
                     labels = c("smaller",  "larger"),
                     name = "" ) +
  scale_fill_manual(values = c("#0067A5", "#F99379"),
                    breaks=c("-1",  "1"),
                    labels=c("smaller",  "larger"),
                    name = "" ) + 
  scale_shape_manual(values = c(19, 15),
                     breaks=c("-1",  "1"),
                     labels=c("smaller",  "larger"),
                     name = "" ) +
  scale_x_discrete(name = "Stimulus type", breaks = c("neighbour", "stranger", "control"), 
                   labels = c("Neighbour", "Stranger", "Control")) +
  scale_y_continuous(name = "Time spent in vicinity") +
  plot_theme
vicinity_plot # Figure S2

# posterior distribution 
vicinity_post <-
  pred_vicinity %>% 
  mutate(Type = str_c(type, size_bin, sep =".")) %>%
  mutate(Type = as.factor(Type)) %>%
  mutate(Type = fct_relevel(Type, c("neighbour.1", "neighbour.-1", "stranger.1", "stranger.-1")))%>%
  ggplot(aes(x = .epred, y = Type, color = Type, fill = Type)) +
  stat_slab(aes(x = .epred, color = Type, fill = Type), point_interval = median_hdi, slab_alpha = .4) +
  stat_pointinterval(aes(x = .epred, color = Type, fill = Type),
                     .width=c(.9, .7), position = position_dodge(width = .3, preserve = "single"), shape = 24) +
  scale_color_manual(values = c("#0067A5", "#A1CAF1", "#F99379", "#BE0032"),
                     name = "",
                     breaks=c("neighbour.1", "neighbour.-1", "stranger.1", "stranger.-1"),
                     labels=c("Larger neighbour", "Smaller neighbour", "Larger stranger", "Smaller stranger"))+
  scale_fill_manual(values = c("#0067A5", "#A1CAF1", "#F99379", "#BE0032"),
                    name = "",
                    breaks=c("neighbour.1", "neighbour.-1", "stranger.1", "stranger.-1"),
                    labels=c("Larger neighbour", "Smaller neighbour", "Larger stranger", "Smaller stranger")) + #light blue = #A1CAF1, light red/pink = #F99379, orange = #F38400, purple = #875692
  scale_y_discrete(name = "Posterior density", labels = c("Larger neighbour", "Smaller neighbour", "Larger stranger", "Smaller stranger")) +
  scale_x_continuous(name = "Time spent within close proximity (< 15 cm) of stimulus (seconds)", lim = c(10, 180)) +
  plot_theme
vicinity_post # Figure 3

