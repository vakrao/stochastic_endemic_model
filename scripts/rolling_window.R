library(tidyverse)
library(ggplot2)
epi_data = read.csv("dat.csv")
# given age-agnostic data, perform 
# rolling window analysis 
# RETURNS: post-processed csv data
window = 365
epi_data %>% 
  filter(vartype == "Xsi1") %>% 
  mutate(winsize = int(time/window)) %>% 
  group_by(winsize,vv) %>% 
  summarise(win_infec = mean(value,na.rm=TRUE))


ggplot(epi_data) +
  geom_point(aes(x=winsize,y=value, color=vv))  + 
  labs(color="Vaccination Coverage Rate")

write.csv(epi_data,"rolling_window_data.csv")

  
  
