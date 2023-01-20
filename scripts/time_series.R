library(data.table)
library(tidyverse)
library(patchwork)
library(ggforce)
library(future)
library(furrr)
library(spatstat)
library(vroom)

setwd("/N/project/endemic_covid/data/raw/R0/R0_5/total/")

vax_list <- list.files()

mean_vax <- function(x) {
  test <- vroom(x) %>%
    filter(vartype %in% c("XD", "Xsi1")) %>%
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

stoch_map_1 <- map_dfr(vax_list, mean_vax)

#time series by vax (color)
vaccineindex <- data.frame(xmin = c(0,seq(365,7299,by=365)-120), xmax=c(0,345,seq(730,7299,by=365)-60),ymin=0,ymax=Inf)

schoolindex <- data.frame(xmax = c(0,seq(365,7299,by=365)-120), xmin=c(0,seq(365,7299,by=365)-120-95),ymin=0,ymax=Inf)

year_label <- paste("Year:", 1:20)

plot1 <- ggplot(stoch_map_1 %>% filter(vartype=="Xsi1")) +
  geom_line(aes(t,value,color=vax),size=1.5) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Infections",
       x="Days",
       title="Model dynamics for R0=5") +
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
  theme_minimal() +
  xlim(300,7000) +
  ylim(0,4e4)

plot1 / plot2 + plot_layout(guides = 'collect')

ggsave("/N/project/endemic_covid/data/raw/R0/R0_5/viz/time_series_R0_5.png", 
       type="cairo", 
       height=7,
       width=10,
       units = "in")

###R0=7

setwd("/N/project/endemic_covid/data/raw/R0/R0_7/total/")

vax_list <- list.files()

mean_vax <- function(x) {
  test <- vroom(x) %>%
    filter(vartype %in% c("XD", "Xsi1")) %>%
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

stoch_map_2 <- map_dfr(vax_list, mean_vax)

#time series by vax (color)
vaccineindex <- data.frame(xmin = c(0,seq(365,7299,by=365)-120), xmax=c(0,345,seq(730,7299,by=365)-60),ymin=0,ymax=Inf)

schoolindex <- data.frame(xmax = c(0,seq(365,7299,by=365)-120), xmin=c(0,seq(365,7299,by=365)-120-95),ymin=0,ymax=Inf)

year_label <- paste("Year:", 1:20)

plot1 <- ggplot(stoch_map_2 %>% filter(vartype=="Xsi1")) +
  geom_line(aes(t,value,color=vax),size=1.5) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Infections",
       x="Days",
       title="Model dynamics for R0=7") +
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

plot2 <- ggplot(stoch_map_2 %>% filter(vartype=="XD"))+
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
  theme_minimal() +
  xlim(300,7000) +
  ylim(0,4e4)

plot1 / plot2 + plot_layout(guides = 'collect')

ggsave("/N/project/endemic_covid/data/raw/R0/R0_7/viz/time_series_R0_7.png", 
       type="cairo", 
       height=7,
       width=10,
       units = "in")

###R0=9

vax_list <- list.files("/N/project/endemic_covid/data/raw/R0/R0_9/total")

mean_vax <- function(x) {
  test <- vroom(x) %>%
    filter(vartype %in% c("XD", "Xsi1")) %>%
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

stoch_map_3 <- map_dfr(vax_list, mean_vax)

#time series by vax (color)
vaccineindex <- data.frame(xmin = c(0,seq(365,7299,by=365)-120), xmax=c(0,345,seq(730,7299,by=365)-60),ymin=0,ymax=Inf)

schoolindex <- data.frame(xmax = c(0,seq(365,7299,by=365)-120), xmin=c(0,seq(365,7299,by=365)-120-95),ymin=0,ymax=Inf)

year_label <- paste("Year:", 1:20)


plot1 <- ggplot(stoch_map %>% filter(vartype=="Xsi1")) +
  geom_line(aes(t,value,color=vax),size=1.5) +
  scale_color_brewer(palette="RdBu") +
  labs(color="Weighted \nVaccination \nCoverage",
       y="Infections",
       x="Days",
       title="Model dynamics for R0=9") +
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

plot2 <- ggplot(stoch_map %>% filter(vartype=="XD"))+
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
  theme_minimal() +
  xlim(300,7000) +
  ylim(0,4e4)

plot1 / plot2 + plot_layout(guides = 'collect')

ggsave("/N/project/endemic_covid/data/raw/R0/R0_9/viz/time_series_R0_9.png", 
       type="cairo", 
       height=7,
       width=10,
       units = "in")
