#load packages
library(tidyverse)
library(vegan)
library(cowplot)
library(brms)
library(car)
library(jtools)
library(tidybayes)
library(posterior)

### SET WORKING DIRECTORY ###
setwd("")

# load and prepare data
dear_enemy <- read.csv("damselfish_NSD_full.csv", header = T, stringsAsFactors = F) # collated behavioural outputs from BORIS v.8.6.2 (Friard & Gamba, 2016) 
dear_enemy$type <- as.factor(dear_enemy$type)
dear_enemy$order <- as.factor(dear_enemy$order)
dear_enemy$size_bin <- as.factor(dear_enemy$size_bin) # -1-smaller, 1-larger
dear_enemy <- dear_enemy %>%
  mutate(total_aggression = bite + display)

# setting priors for bayes models
aggression_priors <- brms::set_prior("normal(0,10)", class = "b")

# Hypothesis 1a: territory holders are less aggressive towards neighbours than non-neighbours
subset <- dear_enemy %>%
  filter(type != "control")

# bayes model
aggression_mod <- brm(total_aggression ~ 0 + type * size_bin + (1|damselID), 
                      family = poisson(),
                      data = subset,
                      iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                      control = list(adapt_delta = 0.95), 
                      prior = aggression_priors,
                      sample_prior = "yes")
summary(aggression_mod)
get_variables(aggression_mod)

# check plots
plot(aggression_mod, ask = FALSE)
pp_check(aggression_mod, resp = 'type', ndraws = 100) # looks good
pp_check(aggression_mod, resp = 'size_bin', ndraws = 100) # looks good
bayes_R2(aggression_mod)

# hypothesis testing 
aggression_mod.hyp <-
  hypothesis(aggression_mod, c("typeneighbour > typestranger"))
print(aggression_mod.hyp)

aggression_mod.hyp_draws<-(aggression_mod.hyp[["samples"]])
as.data.frame(aggression_mod.hyp_draws)
is.data.frame(aggression_mod.hyp_draws)

# hypothesis 1a plot
# set theme
plot_theme <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "Black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=10, colour="Black", face=2),
        legend.position = "none",
        axis.text.y=element_text(size= 10, colour="Black"),
        axis.text.x = element_text( vjust=0.5, size=10, colour="Black"),
        strip.background = element_rect( fill="White"),  
        plot.title=element_text(size=10, face="bold", hjust = 0.5),
        strip.text.x = element_blank())

aggression_mod.hyp_plot <-
  aggression_mod.hyp_draws %>%
  ggplot(aes(x=H1, fill = after_stat(x > 0))) +
  scale_x_continuous(name = NULL, lim = c(-1.2, 1.4)) +
  stat_slab(alpha=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Hypothesis: Territory holders are more aggressive
          towards neighbours than strangers") +
  scale_fill_manual(values=c("grey","#F99379")) +
  annotate("text", x = -1.2, y = 0.99, label = "PP = 1.00\nER = 431.43", hjust = 0, vjust = 1) +
  plot_theme
aggression_mod.hyp_plot

# Hypothesis 1b: territory holders are less aggressive towards larger neighbours than larger non-neighbours
larger <- dear_enemy %>%
  filter(type != "control" & size_bin == 1)

# bayes model
aggression_mod_large <- brm(total_aggression ~ 0 + type + (1|damselID), 
                            family = poisson(),
                            data = larger,
                            iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                            control = list(adapt_delta = 0.99), 
                            prior = aggression_priors,
                            sample_prior = "yes")
summary(aggression_mod_large)
get_variables(aggression_mod_large)

# check plots
plot(aggression_mod_large, ask = FALSE)
pp_check(aggression_mod_large, resp = 'type', ndraws = 100) # looks good
bayes_R2(aggression_mod_large)

# hypothesis testing 
aggression_mod_large.hyp<-
  hypothesis(aggression_mod_large, c("typeneighbour > typestranger"))
print(aggression_mod_large.hyp)

aggression_mod_large.hyp_draws<-(aggression_mod_large.hyp[["samples"]])
as.data.frame(aggression_mod_large.hyp_draws)
is.data.frame(aggression_mod_large.hyp_draws)

# hypothesis 1b plot (larger conspecifics only)
aggression_mod_large.hyp_plot <-
  aggression_mod_large.hyp_draws %>%
  ggplot(aes(x=H1, fill = after_stat(x > 0))) +
  scale_x_continuous(name = NULL, lim = c(-1.2, 1.4)) +
  stat_slab(alpha=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Hypothesis: Territory holders are more aggressive towards
          larger neighbours than larger strangers")+
  scale_fill_manual(values=c("grey","#F99379")) +
  annotate("text", x = -1.2, y = 0.99, label = "PP = 0.66\nER = 1.91", hjust = 0, vjust = 1) +
  plot_theme
aggression_mod_large.hyp_plot


# Hypothesis 1c: territory holders are less aggressive towards smaller neighbours than non-neighbours
smaller <- dear_enemy %>%
  filter(type != "control" & size_bin == -1)

aggression_mod_small <- brm(total_aggression ~ 0 + type + (1|damselID), 
                            family = poisson(),
                            data = smaller,
                            iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                            control = list(adapt_delta = 0.99), 
                            prior = aggression_priors,
                            sample_prior = "yes")
summary(aggression_mod_small)
get_variables(aggression_mod_small)

# check plots
plot(aggression_mod_small, ask = FALSE)
pp_check(aggression_mod_small, resp = 'type', ndraws = 100) # looks good
bayes_R2(aggression_mod_small)

#Hypothesis testing 
aggression_mod_small.hyp <-
  hypothesis(aggression_mod_small, c("typeneighbour > typestranger"))
print(aggression_mod_small.hyp)

aggression_mod_small.hyp_draws<-(aggression_mod_small.hyp[["samples"]])
as.data.frame(aggression_mod_small.hyp_draws)
is.data.frame(aggression_mod_small.hyp_draws)

# smaller-only hypothesis plot
aggression_mod_small.hyp_plot <-
  aggression_mod_small.hyp_draws %>%
  ggplot(aes(x=H1, fill = after_stat(x > 0))) +
  scale_x_continuous(name = "Estimated difference in aggression", lim = c(-1.2, 1.4)) +
  stat_slab(alpha=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Hypothesis: Territory holders are more aggressive towards
          smaller neighbours than smaller strangers")+
  scale_fill_manual(values=c("grey","#F99379")) +
  annotate("text", x = -1.2, y = 0.99, label = "PP = 1.00\nER = 1332.33", hjust = 0, vjust = 1) +
  plot_theme
aggression_mod_small.hyp_plot

cowplot::plot_grid(aggression_mod.hyp_plot, 
                   aggression_mod_large.hyp_plot, 
                   aggression_mod_small.hyp_plot, 
                   ncol=1, align = "vh",
                   labels = c('a', 'b', 'c'))


# Hypothesis 2a: territory holders are more aggressive towards smaller conspecifics than larger conspecifics

# bayes model
aggression_mod2 <- brm(total_aggression ~ 0 + size_bin * type + (1|damselID), 
                      family = poisson(),
                      data = subset,
                      iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                      control = list(adapt_delta = 0.95), 
                      prior = aggression_priors,
                      sample_prior = "yes")
summary(aggression_mod2)
get_variables(aggression_mod2)

# hypothesis testing 
aggression_mod2.hyp <-
  hypothesis(aggression_mod2, c("size_binM1 > size_bin1"))
print(aggression_mod2.hyp)

# check plots
plot(aggression_mod2, ask = FALSE)
pp_check(aggression_mod2, resp = 'type', ndraws = 100) # looks good
pp_check(aggression_mod2, resp = 'size_bin', ndraws = 100) # looks good
bayes_R2(aggression_mod2)

aggression_mod2.hyp_draws<-(aggression_mod2.hyp[["samples"]])
as.data.frame(aggression_mod2.hyp_draws)
is.data.frame(aggression_mod2.hyp_draws)

# hypothesis 2a plot
aggression_mod2.hyp_plot <-
  aggression_mod2.hyp_draws %>%
  ggplot(aes(x=H1, fill = after_stat(x > 0))) +
  scale_x_continuous(name = NULL, lim = c(-5, 8), 
                     breaks = c(-5.0, -2.5, 0.0, 2.5, 5.0, 7.5), labels = c(-5.0, -2.5, 0.0, 2.5, 5.0, 7.5)) +
  stat_slab(alpha=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Hypothesis: Territory holders are more aggressive
          towards smaller conspecifics than larger conspecifics") +
  scale_fill_manual(values=c("grey","#F99379")) +
  annotate("text", x = -5, y = 0.99, label = "PP = 1.00\nER = 7999", hjust = 0, vjust = 1) +
  plot_theme
aggression_mod2.hyp_plot

# Hypothesis 2b: territory holders are more aggressive towards smaller neighbours than larger neighbours
neighbours <- dear_enemy %>%
  filter(type == "neighbour")

# bayes model (neighbour only)
aggression_mod_neigh <- brm(total_aggression ~ 0 + size_bin + (1|damselID), 
                            family = poisson(),
                            data = neighbours,
                            iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                            control = list(adapt_delta = 0.99), 
                            prior = aggression_priors,
                            sample_prior = "yes")
summary(aggression_mod_neigh)
get_variables(aggression_mod_neigh)

# check plots
plot(aggression_mod_neigh, ask = FALSE)
pp_check(aggression_mod_neigh, resp = 'type', ndraws = 100) # looks good
bayes_R2(aggression_mod_neigh)

# hypothesis testing 
aggression_mod_neigh.hyp<-
  hypothesis(aggression_mod_neigh, c("size_binM1 > size_bin1"))
print(aggression_mod_neigh.hyp)

aggression_mod_neigh.hyp_draws<-(aggression_mod_neigh.hyp[["samples"]])
as.data.frame(aggression_mod_neigh.hyp_draws)
is.data.frame(aggression_mod_neigh.hyp_draws)

# hypothesis 2b plot (neighbour only)
aggression_mod_neigh.hyp_plot <-
  aggression_mod_neigh.hyp_draws %>%
  ggplot(aes(x=H1, fill = after_stat(x > 0))) +
  scale_x_continuous(name = NULL, lim = c(-5, 8), 
                     breaks = c(-5.0, -2.5, 0.0, 2.5, 5.0, 7.5), labels = c(-5.0, -2.5, 0.0, 2.5, 5.0, 7.5)) +
  stat_slab(alpha=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Hypothesis: Territory holders are more aggressive towards
          smaller neighbours than larger neighbours")+
  scale_fill_manual(values=c("grey","#F99379")) +
  annotate("text", x = -5, y = 0.99, label = "PP = 0.96\nER = 22.77", hjust = 0, vjust = 1) +
  plot_theme
aggression_mod_neigh.hyp_plot


# Hypothesis 2c: territory holders are more aggressive towards smaller non-neighbours than larger non-neighbours
non_neighbour <- dear_enemy %>%
  filter(type == "stranger")

# bayes model (non-neighbour only)
aggression_mod_nn <- brm(total_aggression ~ 0 + size_bin + (1|damselID), 
                            family = poisson(),
                            data = non_neighbour,
                            iter = 5000, warmup = 1000, cores = 4, chains = 4, 
                            control = list(adapt_delta = 0.99), 
                            prior = aggression_priors,
                            sample_prior = "yes")
summary(aggression_mod_nn)
get_variables(aggression_mod_nn)

# check plots
plot(aggression_mod_nn, ask = FALSE)
pp_check(aggression_mod_nn, resp = 'type', ndraws = 100) # looks good
bayes_R2(aggression_mod_nn)

# hypothesis testing 
aggression_mod_nn.hyp <-
  hypothesis(aggression_mod_nn, c("size_binM1 > size_bin1"))
print(aggression_mod_nn.hyp)

aggression_mod_nn.hyp_draws<-(aggression_mod_nn.hyp[["samples"]])
as.data.frame(aggression_mod_nn.hyp_draws)
is.data.frame(aggression_mod_nn.hyp_draws)

# hypothesis 2c plot
aggression_mod_nn.hyp_plot <-
  aggression_mod_nn.hyp_draws %>%
  ggplot(aes(x=H1, fill = after_stat(x > 0))) +
  scale_x_continuous(name = "Estimated difference in aggression", lim = c(-5, 8), 
                     breaks = c(-5.0, -2.5, 0.0, 2.5, 5.0, 7.5), labels = c(-5.0, -2.5, 0.0, 2.5, 5.0, 7.5)) +
  stat_slab(alpha=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Hypothesis: Territory holders are more aggressive towards
          smaller non-neighbours than larger non-neighbours")+
  scale_fill_manual(values=c("grey","#F99379")) +
  annotate("text", x = -5, y = 0.99, label = "PP = 0.82\nER = 4.63", hjust = 0, vjust = 1) +
  plot_theme
aggression_mod_nn.hyp_plot

# Figure 3
cowplot::plot_grid(aggression_mod2.hyp_plot, 
                   aggression_mod_neigh.hyp_plot, 
                   aggression_mod_nn.hyp_plot, 
                   ncol=1, align = "vh",
                   labels = c('a', 'b', 'c'))


# Full size data
sizes <- read.csv("damselfish_NSD_sizes.csv", header = T, stringsAsFactors = F)
sizes$size_bin <- as.factor(sizes$size_bin)

plot_theme <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "Black"), 
        axis.title.y = element_text(size=10, colour="Black", face=2),
        axis.title.x = element_text(size=10, colour="Black", face=2),
        legend.position = "none",
        axis.text.y=element_text(size= 10, colour="Black"),
        axis.text.x = element_text( vjust=0.5, size=10, colour="Black"),
        strip.background = element_rect( fill="White"),  
        plot.title=element_text(size=10, face="bold", hjust = 0.5),
        strip.text.x = element_blank())

sizes %>%
  group_by(type) %>%
  summarise(mean = mean(size_dif),
            min = min(size_dif),
            max = max(size_dif))

S1a <- sizes %>%
  ggplot(aes(type, size_dif)) +
  geom_boxplot() +
  scale_y_continuous(name = "Body size difference (mm)") +
  scale_x_discrete(name = "Stimulus type", labels = c("Neighbour", "Non-neighbour")) +
  plot_theme

sizes %>%
  group_by(size_bin) %>%
  summarise(mean = mean(size_dif),
            min = min(size_dif),
            max = max(size_dif))

S1b <- sizes %>%
  mutate(size_dif = abs(size_dif)) %>%
  ggplot(aes(size_bin, size_dif)) +
  geom_boxplot() +
  scale_y_continuous(name = "Absolute body size difference (mm)") +
  scale_x_discrete(name = "Conspecific size relative to territory holder", labels = c("Smaller", "Larger")) +
  plot_theme

cowplot::plot_grid(S1a, S1b, ncol=1, align = "vh",labels = c('a', 'b'))


# Posterior distribution plots (Figure S2 and S3)
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
             position = position_jitterdodge(dodge.width = .8)) +
  stat_pointinterval(point_interval = median_hdi, .width = c(.9, .7), fatten_point = 2, 
                     interval_alpha = 0.5, position = position_dodge(width = .8)) +
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
  scale_x_discrete(name = "Stimulus type", labels = c("Neighbour", "Non-neighbour")) +
  scale_y_continuous(name = "Total aggressive displays") +
  plot_theme

# Figure S2
aggression_plot

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

# Figure S3
aggression_post

# Median estimate and HPDI (Table S2)
pred_aggression %>%
  dplyr::group_by(size_bin, type) %>%
  summarise(median = median(.epred))

plot_info <- ggplot_build(aggression_plot)$`data`[[2]]
plot_info['ymax']
plot_info['ymin']