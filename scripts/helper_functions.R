aging <<- function(X) {
  Y = rep(0, length(X))
  Y[2:length(X)] = unlist(lapply(c(2:length(Y)), function(i)
    Y[i] = X[i - 1]))
  return(Y)
}
#Finds vv, modifier for age-based coverage, for the weighted vaccine coverage
#after VS days. G refers to the population level in each age group. 
#ABR-> age-based vaccine coverage
#Returns: list of all vv values
find_vv <<- function(ABR,G,VS) {
  G = G / sum(G)
  WVC = seq(0,1,.1)
  vv_vals = rep(0,11)
  counter = 1
  for(k in WVC){
    for(i in 1:13000){
      test_vv = i/1000
      test_WVC = sum(ABR*G)*test_vv*VS
      if(k - test_WVC < .0001){
        #print(test_vv)
        vv_vals[counter] = test_vv
        counter = counter+ 1
        break
      }
    }
  }
}

vv = find_individual_vv(ABR, 30000000, 60, total_vax_percent)

find_individual_vv <<- function(ABR,G,VS,WVC) {
  G = G / sum(G)
  counter = 1
  vv_val = vv_vals
  for(i in 1:13000){
    test_vv = i/1000
    test_WVC = sum(ABR*G)*test_vv*VS
    if(test_vv - WVC < .0001){
      #print(test_vv)
      vv_val = test_vv
      counter = counter+ 1
      break
    }
  }
}

#Finds area under the curve for some vector X, based on number of years Y 
# and based on some vector Z
finding_area_under_the_curve <<- function(X,Y,Z){
  Y= Y/365.0
  years = rep(seq(1,Y,by=1),10)
  vax_percentage = c()
  all_vals = c()
  for (i in Z){
    relevant_values = filter(X,vv == as.numeric(i))
    print(relevant_values)
    yearly_auc_values = rep(0,Y)
    vax_level = rep(i,Y)
    
    for(j in (1:Y)){
      pruned_values = filter(relevant_values,(year == j))
      total_values = area_under_curve(x=pruned_values$time,y=pruned_values$value,method="step")
      print(pruned_values)
      yearly_auc_values[j] = total_values
      all_vals = c(all_vals,total_values)
    }
    vax_percentage = c(vax_percentage,vax_level)
    #print(vax_percentage)
  }
  complete_values = data.frame(all_vals,years,vax_percentage)
  return(complete_values)
}
#Funcion that finds area under the curve for some dataframe 
#over some period of time given some age groups
age_based_finding_area_under_the_curve <<- function(X,Y,AN){
  all_percent_vaccinated = seq(0.1,1,by=.1)
  years = rep(seq(1,Y,by=1),10)
  vax_percentage = c()
  all_vals = c()
  age_groups = c()
  for (i in all_percent_vaccinated){
    print(i)
    relevant_values = filter(X,percent_vaccinated == as.numeric(i))
    # print(relevant_values)
    yearly_auc_values = rep(0,Y)
    vax_level = rep(i,Y)
    # print(vax_level)
    for(k in AN){
      for(j in (1:Y)){
        pruned_values = filter(relevant_values,(year == j) & (age_group == k))
        # print(pruned_values)
        total_values = area_under_curve(x=pruned_values$time,y=pruned_values$value,method="trapezoid")
        yearly_auc_values[j] = total_values
        all_vals = c(all_vals,total_values)
        age_groups = c(age_groups,k)
        # vax_percentage = c(vax_percentage,vax_level)
      }
    }
    
    #all_vals = c(all_vals,yearly_auc_values)
    vax_percentage = c(vax_percentage,vax_level)
    #print(vax_percentage)
  }
  
  complete_values = data.frame(all_vals,years,vax_percentage,age_groups)
  return(complete_values)
}
param_label <<- function(...) {
  if(include_variant == 1){
  paste0(
    "Parameters: \nStrain 1 R0=",
    R01,
    "\nStrain 1 infection cross protection=",
    percent(C2),
    "\nStrain 2 R0=",
    variant_start_R02,
    "\nStrain 2 infection cross protection=",
    percent(C1),
    "\nVaccine immunity duration=",
    time_of_immunity / 365,
    "y \nInfection immunity duration=",
    time_of_waning_natural / 365,
    "y \nStrain 1 vaccine efficacy=",
    percent(ve1),
    "\nStrain 1 vaccine cross protection=",
    percent(variant_cp),
    "\nStrain 2 vaccine efficacy=",
    percent(variant_eff_2),
    "\nStrain 2 vaccine cross protection=",
    percent(VC2),
    "\nVaccine death reduction=",
    percent(1-sigmad),
    "\nVaccine hospitalization reduction=",
    percent(1-sigmah),
    ...
  ) 
  } else {
    paste0(
      "Parameters: \nStrain 1 R0=",
      R01,
      "\nVaccine immunity duration=",
      time_of_immunity / 365,
      "y \nInfection immunity duration=",
      time_of_waning_natural / 365,
      "y \nStrain 1 vaccine efficacy=",
      percent(ve1),
      "\nVaccine death reduction=",
      percent(sigmad),
      "\nVaccine hospitalization reduction=",
      percent(sigmah),
      ...
    )
  }
}

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


consolidate_files <<- function(filename) {
  
  dat <- list.files(path = "data/new_model/temp/", full.names = TRUE, pattern=".RDS") %>%
  map_dfr(readRDS)
  
  saveRDS(dat,paste0("data/new_model/",filename))
  
  rtdat <- list.files(path = "data/new_model/temp/RT/", full.names = TRUE, pattern=".RDS") %>%
    map_dfr(readRDS)
  
  saveRDS(rtdat,paste0("data/new_model/RT/",filename))
  
  file.remove(list.files("data/new_model/temp/",pattern=".RDS",full.names = TRUE))
  file.remove(list.files("data/new_model/temp/RT",pattern=".RDS",full.names = TRUE))
  
}


load_sim <<- function() {
  if(include_variant == 1){
    readRDS(paste0("data/new_model/wtr0_",R01,"_wtve_",ve1,"_varr0_",variant_start_R02,"_varve_",variant_eff_2,".RDS"))
} else {
    readRDS(paste0("data/new_model/no_var/wtr0_",R01,"_wtve_",ve1,".RDS"))
}
}

load_sim_rt <<- function() {
  if(include_variant == 1){
    readRDS(paste0("data/new_model/RT/wtr0_",R01,"_wtve_",ve1,"_varr0_",variant_start_R02,"_varve_",variant_eff_2,".RDS"))
  } else {
    readRDS(paste0("data/new_model/RT/no_var/wtr0_",R01,"_wtve_",ve1,".RDS"))
  }
}
