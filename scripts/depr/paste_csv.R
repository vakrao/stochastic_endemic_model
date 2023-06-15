#Load libraries====
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggforce)
library(data.table)
library(viridis)
library(furrr)
library(future.apply)
library(progressr)
library(scales)
library(schoolmath)
library(extraDistr)
library(ggpubr)
library(patchwork)

data_all <- list.files(path = "C:\\Users\\varun\\endemic_covid\\stoch_model\\src",     # Identify all csv files in folder
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                            # Store all files in list
  bind_rows                                                       # Combine data sets into one data set 
data_all        