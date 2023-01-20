library(tidyr)
library(dplyr)
library(patchwork)
library(ggforce)
library(spatstat)
library(vroom)
library(purrr)

setwd("/N/project/endemic_covid/data/raw/VE/30/total/")

vax_list <- list.files(pattern=".csv")

title = "Model dynamics for vaccine effectiveness = 30%"

stoch_map_1 <- map_dfr(vax_list, time_series)

saveRDS(stoch_map_1, "cum_age_map.RDS")

#time series====
time_plot <- time_series_graph_func(title)

ggsave("/N/project/endemic_covid/data/raw/VE/30/viz/time_series_30.png", 
       type="cairo", 
       height=7,
       width=10,
       units = "in",
       plot=time_plot)

#rt====

rt_plot <- rt_plot_func(title)

ggsave("/N/project/endemic_covid/data/raw/VE/30/viz/rt_ve_30.png", 
       type="cairo", 
       height=7,
       width=9,
       units = "in",
       plot=rt_plot)

#age of infection====

age_map <- age_heat_map(vax_list)

saveRDS(age_map, "strat_age_map.RDS")

age_inf_map <- age_inf_plot_func(title)

ggsave("/N/project/endemic_covid/data/raw/VE/30/viz/age_diff_30.png", 
       type="cairo", 
       height=7,
       width=10,
       units = "in",
       plot=age_inf_map)

#death heat map====
heat1 <- heat_death_func(title)

ggsave("/N/project/endemic_covid/data/raw/VE/30/viz/heat_maps_deaths_30.png", 
       type="cairo", 
       height=8,
       width=10,
       units = "in",
       plot=heat1)

#infection heat map====
heat2 <- heat_inf_func(title)

ggsave("/N/project/endemic_covid/data/raw/VE/30/viz/heat_maps_infections_30.png", 
       type="cairo", 
       height=8,
       width=10,
       units = "in",
       plot=heat2)
