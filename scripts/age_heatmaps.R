library(data.table)
library(tidyverse)
library(patchwork)
library(ggforce)
library(spatstat)
library(furrr)
library(future)

age_group <- c(-0.5, 5, 18, 35, 50, 65, Inf)

vax_list <- list.files(path = "/N/project/endemic_covid/data/totals/", full.names = TRUE, pattern=".csv")

mean_vax <- function(x) {
  test <- fread(x) %>%
    filter(vartype %in% c("XD", "Xsi1", "N")) %>%
    group_by(t, age, vv, vartype) %>%
    summarise(value = mean(value)) %>% #summarise simulations
    ungroup() %>%
    mutate(age_group = cut(age, age_group)) %>%
    group_by(age_group,vv,vartype) %>%
    mutate(year = ceiling(t/365)) %>%
    filter(year>0) %>%
    select(-t,age) %>%
    group_by(age_group,year,vv,vartype) %>%
    summarise(value=sum(value)) %>%
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
  return(test)
}
plan("multisession")

age_map <- future_map_dfr(vax_list, mean_vax)

x1 <- age_map %>%
  filter(vartype == "XD") %>%
  ungroup() %>%
  mutate(value = big_round(value),
         cut_value = cut(
           value,
           breaks = c(0, 1e3, 5e3, 1e4, 1e5, 2.5e5, 5e5, 7.5e5, 1e6, 1.5e6,2e6, Inf),
           labels = c(
             "0-1k",
             "1k-5k",
             "5k-10k",
             "10k-100k",
             "100k-250k",
             "250k-500k",
             "500k-750k",
             "750k-1M",
             "1M-1.5M",
             "1.5M-2M",
             "2M+"
           )
         ))

ggplot(x1) +
  geom_tile(aes(year, vax, fill = cut_value)) +
  scale_fill_brewer(palette = "Reds")  +
  labs(fill = "Deaths",
       y = "Coverage",
       x = "Years",
       title = "Varying COVID-19 Vaccination Coverage (Deaths)") +
  scale_x_continuous(limits = c(1, 20)) +
  theme_minimal() +
  facet_wrap( ~ age_group)

#infection

x2 <- age_map %>%
  filter(vartype == "Xsi1") %>%
  ungroup() %>%
  mutate(value1 = big_round(value),
         cut_value = cut(
           value1,
           breaks = c(0, 1e7, 2e7, 4e7, 6e7, 8e7, 1e8, 1.5e8, 2e8),
           labels = c(
             "0-10M",
             "10M-20M",
             "20M-40M",
             "40M-60M",
             "60M-80M",
             "80M-100M",
             "100M-150M",
             "150M-200M"
           )
         ))

ggplot(x2) +
  geom_tile(aes(year, vax, fill = cut_value)) +
  scale_fill_brewer(palette = "Reds")  +
  labs(fill = "Infections",
       y = "Coverage",
       x = "Years",
       title = "Varying COVID-19 Vaccination Coverage (Infections)") +
  scale_x_continuous(limits = c(0, 19)) +
  theme_minimal() +
  facet_wrap( ~ age_group)
