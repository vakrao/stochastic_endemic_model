#Data prep R03====
R01 <- 3
variant_start_R02 <- 5

var_dat <-
  grep(
    paste0("_variant_", variant_start_R02, ".*"),
    list.files(path = "data/temp_model/varun/variant_only/", full.names = TRUE),
    perl = T,
    value = T
  ) %>%
  map_dfr(readRDS)

var_sum_dat <- var_dat %>%
  group_by(vv,Var1,age_group) %>% 
  filter(time > 365+230) %>% 
  mutate(season=1:6705 ) %>% 
  mutate(year=ceiling(season/365)) %>% 
  filter(year< 19) %>% 
  group_by(vv,Var1,year, age_group) %>% 
  summarise(value=sum(value)) %>% 
  mutate(vv = as.factor(vv))


#Deaths====

##Variant====
var_sum_infections <- var_sum_dat %>% 
  filter(Var1 %in% c("New Variant Infections")) %>%   
  ungroup() %>% 
  mutate(value = big_round(value)) %>% 
  mutate(cut_value = cut_width(value,8e6,dig.lab=10,boundary=0))

levels(var_sum_infections$cut_value) <- pretty_labels(levels(var_sum_infections$cut_value))
levels(var_sum_infections$vv) <- scales::percent(seq(0, 1, by = .1))

var_infections_heat <- ggplot(var_sum_infections) +
  geom_tile(aes(year, vv, fill=cut_value)) +
  scale_fill_brewer(palette = "Reds")  + 
  labs(title = "Population Weighted Average Vaccination Coverage",
       subtitle = "Wild Type + Variant Model",
       y = "Weighted Vaccination Coverage",
       fill = "Infections",
       x = "Years",
       tag = param_label()
  ) +
  scale_x_continuous(limits=c(0,19)) +
  theme_minimal() +  
  facet_wrap(~age_group) +
  theme(plot.tag = element_text(hjust=0),
        plot.tag.position = c(0.93, 0.25),
        plot.margin=unit(c(1,3,1,1),"cm"))

var_infections_heat

ggsave(
  paste0("pngs/var_infections_heat_",variant_start_R02,".png"),
  var_infections_heat,
  width = 16,
  height = 8
)

##Proportion====
###Variant====
age_prop <- var_sum_dat %>% 
  filter(Var1 %in%  c("New Variant Infections")) %>% 
  ungroup() %>% 
  group_by(vv,year) %>% 
  summarise(prop=value/sum(value),
            Coverage=vv,
            age_group)

levels(age_prop$Coverage) <- scales::percent(seq(0, 1, by = .1))  

ggplot(age_prop) +
  geom_col(aes(year,prop,fill=age_group)) +
  facet_wrap(~Coverage,labeller = label_both) +
  scale_y_continuous(labels = scales::percent) +
  labs(y="New infections proportion",
       x="Year",
       fill="Age Group",
       title = "Proportion of new infections by strain",
       subtitle="Wild Type + Variant model")


#Data prep R07=========================================================================
R01 <- 3
variant_start_R02 <- 7

var_dat <-
  grep(
    paste0("_variant_", variant_start_R02, ".*"),
    list.files(path = "data/temp_model/varun/variant_only/", full.names = TRUE),
    perl = T,
    value = T
  ) %>%
  map_dfr(readRDS)

var_sum_dat <- var_dat %>%
  group_by(vv,Var1,age_group) %>% 
  filter(time > 365+230) %>% 
  mutate(season=1:6705 ) %>% 
  mutate(year=ceiling(season/365)) %>% 
  filter(year< 19) %>% 
  group_by(vv,Var1,year, age_group) %>% 
  summarise(value=sum(value)) %>% 
  mutate(vv = as.factor(vv))


#Deaths====

##Variant====
var_sum_infections <- var_sum_dat %>% 
  filter(Var1 %in% c("New Variant Type Infections")) %>%   
  ungroup() %>% 
  mutate(value = big_round(value)) %>% 
  mutate(cut_value = cut_width(value,8e6,dig.lab=10,boundary=0))

levels(var_sum_infections$cut_value) <- pretty_labels(levels(var_sum_infections$cut_value))
levels(var_sum_infections$vv) <- scales::percent(seq(0, 1, by = .1))

var_infections_heat <- ggplot(var_sum_infections) +
  geom_tile(aes(year, vv, fill=cut_value)) +
  scale_fill_brewer(palette = "Reds")  + 
  labs(title = "Population Weighted Average Vaccination Coverage",
       subtitle = "Wild Type + Variant Model",
       y = "Weighted Vaccination Coverage",
       fill = "Infections",
       x = "Years",
       tag = param_label()
  ) +
  facet_wrap(~age_group) +
  scale_x_continuous(limits=c(0,19)) +
  theme_minimal() +
  theme(plot.tag = element_text(hjust=0),
        plot.tag.position = c(0.93, 0.25),
        plot.margin=unit(c(1,3,1,1),"cm"))

var_infections_heat

ggsave(
  paste0("pngs/var_infections_heat_",variant_start_R02,".png"),
  var_infections_heat,
  width = 16,
  height = 8
)


#Data prep R09====
R01 <- 3
variant_start_R02 <- 9

var_dat <-
  grep(
    paste0("_variant_", variant_start_R02, ".*"),
    list.files(path = "data/temp_model/varun/variant_only/", full.names = TRUE),
    perl = T,
    value = T
  ) %>%
  map_dfr(readRDS)

var_sum_dat <- var_dat %>%
  group_by(vv,Var1,age_group) %>% 
  filter(time > 365+230) %>% 
  mutate(season=1:6705 ) %>% 
  mutate(year=ceiling(season/365)) %>% 
  filter(year< 19) %>% 
  group_by(vv,Var1,year, age_group) %>% 
  summarise(value=sum(value)) %>% 
  mutate(vv = as.factor(vv))


#Deaths====

##Variant====
var_sum_infections <- var_sum_dat %>% 
  filter(Var1 %in% c("New Variant Infections")) %>%   
  ungroup() %>% 
  mutate(value = big_round(value)) %>% 
  mutate(cut_value = cut_width(value,8e6,dig.lab=10,boundary=0))

levels(var_sum_infections$cut_value) <- pretty_labels(levels(var_sum_infections$cut_value))
levels(var_sum_infections$vv) <- scales::percent(seq(0, 1, by = .1))

var_infections_heat <- ggplot(var_sum_infections) +
  geom_tile(aes(year, vv, fill=cut_value)) +
  scale_fill_brewer(palette = "Reds")  + 
  labs(title = "Population Weighted Average Vaccination Coverage",
       subtitle = "Wild Type + Variant Model",
       y = "Weighted Vaccination Coverage",
       fill = "Infections",
       x = "Years",
       tag = param_label()
  ) +
  scale_x_continuous(limits=c(0,19)) +
  theme_minimal() +
  facet_wrap(~age_group) +
  theme(plot.tag = element_text(hjust=0),
        plot.tag.position = c(0.93, 0.25),
        plot.margin=unit(c(1,3,1,1),"cm"))

var_infections_heat

ggsave(
  paste0("pngs/var_infections_heat_",variant_start_R02,".png"),
  var_infections_heat,
  width = 16,
  height = 8
)