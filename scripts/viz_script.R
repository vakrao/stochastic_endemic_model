library(tidyr)
library(dplyr)
library(patchwork)
library(ggforce)
library(spatstat)
library(vroom)
library(purrr)

big_round <<- function(x) {
  case_when(
    x < 1 ~ 0,
    x <= 1e2 ~ round(x, 1),
    x <= 1e3 & x > 1e2 ~ round(x, -1),
    x <= 1e4 & x > 1e3 ~ round(x, -2),
    x <= 1e5 & x > 1e4 ~ round(x, -3),
    x <= 1e6 & x > 1e5 ~ round(x, -4),
    x <= 1e7 & x > 1e6 ~ round(x, -5),
    x <= 1e8 & x > 1e7 ~ round(x, -6),
    x <= 1e9 & x > 1e8 ~ round(x, -7),
    x <= 1e10 & x > 1e9 ~ round(x, -8),
    x <= 1e11 & x > 1e10 ~ round(x, -9)) 
  
}

SI_format <<- function(x) {
  dplyr::case_when(
    x < 1e3 ~ as.character(x),
    x < 1e6 ~ paste0(as.character(x/1e3), "K"),
    x < 1e9 ~ paste0(as.character(x/1e6), "M"),
    x < 1e12 ~ paste0(as.character(x/1e9), "B"),
    x < 1e15 ~ paste0(as.character(x/1e12), "T"),
    TRUE ~ "..."
  )
}

pretty_labels <<- function(x) {
  
  first_bracket <- str_extract(x,"^.") 
  last_bracket <- str_extract(x,".$")
  first_number <- SI_format(as.numeric(str_extract(x,"(\\d+)")))
  last_number <- SI_format(as.numeric(str_extract(x,"(\\d+)(?!.*\\d)")))
  
  new_label <- paste0(first_bracket,first_number,",",last_number,last_bracket)
  return(new_label)
}

time_series <- function(x) {
  test <- vroom(x,
                delim = ",",
                col_types = c(
                  t = "d",
                  vv = "d",
                  age = "d",
                  value = "d",
                  sim_number = "d",
                  vartype = "c"
                )
  ) %>%
    filter(vartype %in% c("XD", "Xsi1", "Rt")) %>% 
    group_by(t,vv,vartype,age) %>% #mean of simulations
    summarise(value=mean(value)) %>% 
    ungroup() %>% 
    group_by(t,vv,vartype) %>% 
    summarise(value=sum(value)) %>% #sum of ages
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

age_group <- c(-0.5, 5, 18, 35, 50, 65, Inf)

age_inf <- function(x) {
  test <- vroom(x,
                delim = ",",
                col_types = c(
                  t = "d",
                  vv = "d",
                  age = "d",
                  value = "d",
                  sim_number = "d",
                  vartype = "c"
                )
  ) %>%
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

heat_maps <- function(x) {
  test <- vroom(x,
                delim = ",",
                col_types = c(
                  t = "d",
                  vv = "d",
                  age = "d",
                  value = "d",
                  sim_number = "d",
                  vartype = "c"
                )
  ) %>%
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
  
comb_func <- function(x) {
  dat <- vroom(x,
               delim = ",",
               col_types = c(
                 t = "d",
                 vv = "d",
                 age = "d",
                 value = "d",
                 sim_number = "d",
                 vartype = "c"
               )
  ) 
  
  age_inf <- dat %>%
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
  
  heat_maps <- dat %>% 
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
  df_list <- list(age_inf, heat_maps) # we create the list
  
  
  return(df_list)
  
}

age_heat_map <- function(x) {

test_dat <- map(x, comb_func)

return(map(transpose(test_dat), bind_rows))

}

#plots

time_series_graph_func <- function(title) {
  
vaccineindex <- data.frame(xmin = c(0,seq(365,7299,by=365)-120), xmax=c(0,345,seq(730,7299,by=365)-60),ymin=0,ymax=Inf)

schoolindex <- data.frame(xmax = c(0,seq(365,7299,by=365)-120), xmin=c(0,seq(365,7299,by=365)-120-95),ymin=0,ymax=Inf)

year_label <- paste("Year:", 1:20)

plot1 <- ggplot(stoch_map_1 %>% filter(vartype=="Xsi1")) +
  geom_line(aes(t,value,color=vax),size=1.5) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Infections",
       x="Days",
       title=title) +
  geom_vline(xintercept = seq(0,7299,by=365),color="darkblue",alpha=.4) +
  annotate("text",x=seq(0,7299,by=365)+30,y=12000,label=paste("Year:", 1:20),hjust=0,angle=90,size=3,alpha=.5) +
  geom_rect(data=vaccineindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="Vaccine"),alpha=.2) +
  geom_rect(data=schoolindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill="Summer"),alpha=.2) +
  scale_fill_manual('Season',
                    values = c('darkred', "darkblue"),
                    guide = guide_legend(override.aes = list(alpha = 0.2))) +
  theme_light() +
  scale_y_continuous(labels = scales::label_number_si(accuracy = 0.5),limits = c(0,2e6)) +
  xlim(300,7000)

plot2 <- ggplot(stoch_map_1 %>% filter(vartype=="XD"))+
  geom_line(aes(t,value,color=vax),size=1.5) +
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(labels = seq(365,7299,by=365),
                     breaks = seq(365,7299,by=365)) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Deaths",
       x="Days") +
  geom_vline(xintercept = seq(0,7299,by=365),color="darkblue",alpha=.4) +
  geom_rect(data=vaccineindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="Vaccine"),alpha=.2) +
  geom_rect(data=schoolindex, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill="Summer"),alpha=.2) +
  scale_fill_manual('Season',
                    values = c('darkred', "darkblue"),
                    guide = guide_legend(override.aes = list(alpha = 0.2))) +
  theme_light() +
  xlim(300,7000) +
  ylim(0,2.5e3)

time_plot <- plot1 / plot2 + plot_layout(guides = 'collect')

return(time_plot)

}

rt_plot_func <- function(title) {

rt_plot <- ggplot(stoch_map_1 %>% filter(vartype=="Rt")) +
  geom_line(aes(t,value, color=vax)) +
  geom_hline(aes(yintercept=1),color="red") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Rt",
       x="Days",
       title=title)+
  scale_color_brewer(palette="RdBu") +
  ylim(0.5,1.5) +
  xlim(300,7300) +
  theme_light()

return(rt_plot)

}

age_inf_plot_func <- function(title) {

age_diff <- age_map[[1]] %>% 
  ungroup() %>% 
  pivot_wider(c(1,2,6),names_from = "vartype",values_from = "mean") %>% 
  mutate(diff=Xsi1-N)

age_inf_map <- ggplot(age_diff) +
  geom_line(aes(quarter,diff,color=vax))+
  labs(y="Difference",
       title=paste("Population normalized average age of infection,", title),
       x="Quarter",
       color="Weighted \nvaccination \ncoverage")+
  scale_color_brewer(palette="RdBu",guide=guide_legend(reverse = TRUE)) +
  theme_light()

return(age_inf_map)
}

heat_death_func <- function(title) {
  
x1 <-  age_map[[2]]%>%
  filter(vartype == "XD") %>%
  ungroup() %>%
  mutate(value = big_round(value),
         cut_value = cut(
           value,
           breaks = c(0, 1e3, 1e4, 1e5, 2.5e5, 5e5, 7.5e5, 1e6, 1.5e6, Inf),
           labels = c(
             "0-1k",
             "1k-10k",
             "10k-100k",
             "100k-250k",
             "250k-500k",
             "500k-750k",
             "750k-1M",
             "1M-1.5M",
             "1.5M+"
           )
         ))

heat1 <- ggplot(x1) +
  geom_tile(aes(year, vax, fill = cut_value)) +
  scale_fill_brewer(palette = "Reds")  +
  labs(fill = "Deaths",
       y = "Coverage",
       x = "Years",
       title = paste("Varying COVID-19 Vaccination Coverage (Deaths),", title)) +
  scale_x_continuous(limits = c(1, 20)) +
  theme_light() +
  facet_wrap( ~ age_group)

return(heat1)

}

heat_inf_func <- function(title) {

x2 <- age_map[[2]]%>%
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

heat2 <- ggplot(x2) +
  geom_tile(aes(year, vax, fill = cut_value)) +
  scale_fill_brewer(palette = "Reds")  +
  labs(fill = "Infections",
       y = "Coverage",
       x = "Years",
       title = paste("Varying COVID-19 Vaccination Coverage (Infections),", title)) +
  scale_x_continuous(limits = c(0, 19)) +
  theme_light() +
  facet_wrap( ~ age_group)

return(heat2) 

}

hist_inf_func <- function(x) {
  test <- vroom(x,
                delim = ",",
                col_types = c(
                  t = "d",
                  vv = "d",
                  age = "d",
                  value = "d",
                  sim_number = "d",
                  vartype = "c"
                )
  ) %>%
    filter(vartype %in% c("Xsi1")) %>% 
    group_by(t,vv,vartype,age) %>% #mean of simulations
    summarise(value=mean(value)) %>% 
    ungroup() %>% 
    group_by(t,vv,vartype) %>% 
    summarise(value=sum(value)) %>% #sum of ages
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
           )) %>% 
           group_by(vax,age) %>% 
           mutate(year = t %/% 365) %>%
           ungroup() %>%
           group_by(vax,age,year) %>%
           summarise(infections=sum(value))
   hist <- ggplot(test %>% filter(vax="50%")) + 
              geom_col(aes(age,infections)) + 
              facet_wrap(~year,labeller=label_both)
   return(hist)

}


