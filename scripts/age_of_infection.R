library(data.table)
library(tidyverse)
library(patchwork)
library(ggforce)
library(future)
library(furrr)
library(spatstat)


vax_list <- c(
  "data/raw/30_run/agetotal_sim_vax_level_00.csv",
  "data/raw/30_run/agetotal_sim_vax_level_10.csv",
  "data/raw/30_run/agetotal_sim_vax_level_20.csv",
  "data/raw/30_run/agetotal_sim_vax_level_30.csv",
  "data/raw/30_run/agetotal_sim_vax_level_40.csv",
  "data/raw/30_run/agetotal_sim_vax_level_50.csv",
  "data/raw/30_run/agetotal_sim_vax_level_60.csv",
  "data/raw/30_run/agetotal_sim_vax_level_70.csv",
  "data/raw/30_run/agetotal_sim_vax_level_80.csv",
  "data/raw/30_run/agetotal_sim_vax_level_90.csv",
  "data/raw/30_run/agetotal_sim_vax_level_100.csv"
)

mean_vax <- function(x) {
  test <- fread(x) %>%
  filter(vartype %in% c("Xsi1","N")) %>%
  group_by(age,vv, vartype) %>% 
  mutate(quarter = ceiling(t / (365/4)),
         year = ceiling(t/365)) %>% 
  group_by(quarter,year,vv,vartype) %>% 
  summarise(hand_calc=sum(age*value)/sum(value),
            mean=weighted.mean(age,value),
            median=weighted.median(age,value)) %>% 
  mutate(vax = as.character(vv),
         vax = factor(
           case_when(
             vax == "0" ~ "0%",
             vax == "1.21" ~ "10%",
             vax == "2.42" ~ "20%",
             vax == "3.64" ~ "30%",
             vax == "4.84" ~ "40%",
             vax == "6.05" ~ "50%",
             vax == "7.26" ~ "60%",
             vax == "8.47" ~ "70%",
             vax == "9.68" ~ "80%",
             vax == "10.89" ~ "90%",
             vax == "12.1" ~ "100%",
             TRUE ~ vax
           ),
           levels = c(
             "0%",
             "10%",
             "20%",
             "30%",
             "40%",
             "50%",
             "60%",
             "70%",
             "80%",
             "90%",
             "100%"
           )
         ))
}
  
  age_inf_map <- map_dfr(vax_list, mean_vax)
  
  age_diff <- age_inf_map %>% 
    ungroup() %>% 
    pivot_wider(c(1,2,8,4,6),names_from = "vartype",values_from = mean) %>% 
    mutate(diff=N-Xsi1)

ggplot(age_diff) +
         geom_line(aes(quarter,diff,color=vax))+
  labs(y="Difference",
       title="Population normalized average age of infection",
       x="Quarter",
       color="Weighted \nvaccination \ncoverage")+
  scale_color_brewer(palette="RdBu",guide=guide_legend(reverse = TRUE)) +
  theme_minimal()

ggsave("diff_age_of_infection_plot.png", type="cairo")
