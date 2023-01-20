library(data.table)
library(tidyverse)
library(patchwork)
library(ggforce)
library(future)
library(furrr)
library(spatstat)

vax_list <- c(
  "data/raw/30_run/agetotal_sim_vax_level_00.csv",
  "data/raw/30_run/agetotal_sim_vax_level_50.csv",
  "data/raw/30_run/agetotal_sim_vax_level_100.csv"
)

mean_vax <- function(x) {
  test <- fread(x) %>%
  filter(vartype %in% c("XD", "Xsi1")) %>%
  group_by(t,vv,vartype,age) %>% 
  summarise(value=mean(value)) %>% 
    ungroup() %>% 
    mutate(year = ceiling(t/365)) %>%
    group_by(year,vv,vartype) %>% 
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
}

stoch_map <- map_dfr(vax_list, mean_vax)

# p1 <- ggplot(test) +
#   geom_smooth(aes(y=mean,x=quarter,color=vax)) +
#   labs(x="Quarter",
#        y="Mean Age of Infection",
#        color="Weighted Vaccination Coverage") +
#   theme_minimal()
# 
# p2 <- ggplot(test) +
#   geom_point(aes(y=median,x=quarter,color=vax)) +
#   geom_smooth(aes(y=median,x=quarter,color=vax),se=FALSE,span=.3) +
#   labs(x="Quarter",
#        y="Median Age of Infection",
#        color="Weighted Vaccination Coverage") +
#   theme_minimal()
# 
# p1 / p2 + plot_layout(guides = 'collect')
# 
# ggsave("age_of_infection_plot.png", type="cairo")

# mean_vax_test <- mean_vax %>% 
#   #mutate(age= cut(age, age_group)) %>% 
#            group_by(t,age,vv,vartype) %>% 
#            summarise(value=mean(value))

# 
# quant_vax <- list.files(path = "data/quantile/", full.names = TRUE, pattern=".csv") %>%
#   map_dfr(fread)

#no ages data
#mean_vax <- fread("data/raw/mean_exp1.csv")
# 
# quant_vax <- fread("data/all_vax_quantile_values_exp1.csv")
# 

mean_vax1 <- mean_vax_test %>% filter(vartype %in% c("Xsi1")) %>%
  mutate(vv = round(vv, 2)) %>%
  # mutate(year = ceiling(t / 30)) %>%
  # group_by(year, vv) %>%
  # summarise(value = sum(value)) %>%
  mutate(vax = as.character(vv),
         vax = factor(
           case_when(
             vax == "0" ~ "0%",
             vax == "1.21" ~ "10%",
             vax == "2.42" ~ "20%",
             vax == "3.64" ~ "30%",
             vax == "4.84" ~ "40%",
             vax == "6.05" ~ "50%",
             vax == "7.62" ~ "60%",
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

mean_vax2 <- mean_vax_test %>% filter(vartype == "XD") %>%
  mutate(vv = round(vv, 2)) %>%
  # mutate(month = ceiling(t / 30)) %>%
  # group_by(vartype, t, vv) %>%
  # summarise(value = sum(value)) %>%
  mutate(vax = as.character(vv),
         vax = factor(
           case_when(
             vax == "0" ~ "0%",
             vax == "1.21" ~ "10%",
             vax == "2.42" ~ "20%",
             vax == "3.64" ~ "30%",
             vax == "4.84" ~ "40%",
             vax == "6.05" ~ "50%",
             vax == "7.62" ~ "60%",
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

#time series facet (black)

ggplot(stoch_map %>% filter(vartype=="Xsi1")) +
  geom_line(aes(year,value,color=vax),size=1) +
  scale_y_continuous(labels = scales::label_number_si(accuracy = .05)) +
  scale_color_brewer(palette="Set1") +
  labs(title="Stochastic Mean Infections by Weighted Vaccination Coverage",
       color="Weighted \nVaccination \nCoverage",
       y="Infections") +
  facet_wrap(~age_group)

ggplot(test_map %>% filter(vartype=="XD",year>1)) +
  geom_line(aes(year,value,color=vax),size=1) +
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_color_brewer(palette="Set1") +
  labs(title="Stochastic Mean Deaths by Weighted Vaccination Coverage",
       color="Weighted \nVaccination \nCoverage",
       y="Deaths") +
  facet_wrap(~age_group)

#time series by vax (color)
vaccineindex <- data.frame(xmin = c(0,seq(365,7299,by=365)-120), xmax=c(0,345,seq(730,7299,by=365)-60),ymin=0,ymax=Inf)

schoolindex <- data.frame(xmax = c(0,seq(365,7299,by=365)-120), xmin=c(0,seq(365,7299,by=365)-120-95),ymin=0,ymax=Inf)

year_label <- paste("Year:", 1:20)


ggplot(stoch_map %>% filter(vartype=="Xsi1")) +
  geom_line(aes(year,value,color=vax),size=1.5) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Infections",
       x="Days") +
  geom_vline(xintercept = seq(0,7299,by=365),color="darkblue",alpha=.4) +
  annotate("text",x=seq(0,7299,by=365)+30,y=12000,label=paste("Year:", 1:20),hjust=0,angle=90,size=3,alpha=.5) +
  geom_rect(data=vaccineindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="Vaccine"),alpha=.2) +
  geom_rect(data=schoolindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill="Summer"),alpha=.2) +
  scale_fill_manual('Season',
                    values = c('darkred', "darkblue"),
                    guide = guide_legend(override.aes = list(alpha = 0.2))) +
  theme_minimal() +
  scale_y_continuous(labels = scales::label_number_si(accuracy = 0.5),limits = c(0,2e6)) +
  xlim(300,7000)

plot2 <- ggplot(mean_vax2) +
  geom_line(aes(t,value,color=vax),size=1.5) +
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(labels = seq(365,7299,by=365),
                     breaks = seq(365,7299,by=365)) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Deaths",
       x="Days") +
  xlim(300,7000) +
  ylim(0,7.5e3) +
  geom_vline(xintercept = seq(0,7299,by=365),color="darkblue",alpha=.4) +
  geom_rect(data=vaccineindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="Vaccine"),alpha=.2) +
  geom_rect(data=schoolindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill="Summer"),alpha=.2) +
  scale_fill_manual('Season',
                    values = c('darkred', "darkblue"),
                    guide = guide_legend(override.aes = list(alpha = 0.2))) +
  theme_minimal()

plot1 / plot2 + plot_layout(guides = 'collect')

#quantiles
# 
# quant_vax1 <- quant_vax %>% filter(vartype == "Xsi1") %>%
#   mutate(vv = round(vv, 2)) %>%
#   # mutate(month = ceiling(t / 30)) %>%
#   # group_by(vartype, t, vv, quantile = `Unnamed: 2`) %>%
#   # summarise(value = sum(value)) %>%
#   mutate(vax = as.character(vv),
#          vax = factor(
#            case_when(
#              vax == "0" ~ "0%",
#              vax == "1.21" ~ "10%",
#              vax == "2.42" ~ "20%",
#              vax == "3.64" ~ "30%",
#              vax == "4.84" ~ "40%",
#              vax == "6.05" ~ "50%",
#              vax == "7.26" ~ "60%",
#              vax == "8.47" ~ "70%",
#              vax == "9.68" ~ "80%",
#              vax == "10.89" ~ "90%",
#              vax == "12.1" ~ "100%",
#              TRUE ~ vax
#            ),
#            levels = c(
#              "0%",
#              "10%",
#              "20%",
#              "30%",
#              "40%",
#              "50%",
#              "60%",
#              "70%",
#              "80%",
#              "90%",
#              "100%"
#            )
#          ))
# 
# 
# quant_comp <-
#   quant_vax1 %>% mutate(
#     quant_names = case_when(
#       V3 == .25 ~ "Lower",
#       V3 == .5 ~ "Median",
#       V3 == .75 ~ "Upper"
#     )
#   ) %>% 
#   pivot_wider(c(1,2,5), names_from = quant_names, values_from = value)

# ggplot(quant_comp %>% filter(t>365, t<3000)) +
#   geom_line(aes(t, Median)) +
#   geom_point(aes(t, Lower), color="red",alpha=.3) +
#   geom_ribbon(aes(t, ymin = Lower, ymax = Upper), fill = "grey20",alpha=.3) +
#   scale_y_continuous(labels = scales::label_number_si(accuracy = 1)) +
#   labs(title = "Stochastic Quantile Infections by Weighted Vaccination Coverage",
#        color = "Weighted \nVaccination \nCoverage",
#        y = "Infections") +
#   facet_wrap( ~ vax) +
#   theme_minimal() +
#   ylim(1e6,1.2e6)
# 
# ggplot(quant_comp) +
#   geom_line(aes(t, Median)) +
#   geom_ribbon(aes(t, ymin = Lower, ymax = Upper), fill = "grey20",alpha=.3) +
#   scale_y_continuous(labels = scales::label_number_si()) +
#   labs(title = "Stochastic Quantile Cases by Weighted Vaccination Coverage",
#        color = "Weighted \nVaccination \nCoverage",
#        y = "Cases") +
#   facet_wrap( ~ vax) +
#   theme_minimal() +
#   ylim(0,25000) +
#   xlim(0,365*3)
  

ggplot(na.omit(mean_vax1) %>% filter(t>300)) +
  geom_line(aes(t,value,color=vax),size=1) +
  scale_y_continuous(labels = scales::label_number_si(accuracy = .05)) +
  scale_color_brewer(palette="RdBu") +
  labs(title="Stochastic Mean Infections by Weighted Vaccination Coverage",
       color="Weighted \nVaccination \nCoverage",
       y="Infections") +
  facet_wrap(~age)



ggplot(na.omit(mean_vax2) %>% filter(t>300)) +
  geom_line(aes(t,value,color=vax),size=1) +
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_color_brewer(palette="RdBu") +
  labs(title="Stochastic Mean Deaths by Weighted Vaccination Coverage",
       color="Weighted \nVaccination \nCoverage",
       y="Deaths") +
  facet_wrap(~age)
