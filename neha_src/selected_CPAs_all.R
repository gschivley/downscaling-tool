setwd('/Users/npatankar/D/US_Policy/downscaling_Jan23/capacity_downscaling/result_data')
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(RSQLite)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(reshape)
library(arrow)
library(broom)
library(geojsonio)
library(mapproj)
library(rjson)
library(data.table)
library(ggmap)
library(gridExtra)
library(readr)

#################################################################################################################################
# Assumptions for the input data
#################################################################################################################################
#Total site boundary power densities:
#Solar: 45 MW/km2
#Wind: 2.7 MW/km2
#offWind floating: 8 MW/km2
#offWind fixed: 5 MW/km2
#Maximum solar development density per candidate project area/cluster: 20%
#Maximum onshore wind development density per candidate project area/cluster: 50%
#Maximum offshore wind development density per candidate project area/cluster: 50%
#Directly impacted lands (roads, panels/turbines, substation etc) as % of total site boundary:
#Solar: 91%
#Wind: 1%
#Derating based on population density
#RESOURCE_POP_DEN_RETAIN = {
#  "solar": {
#    "pop_den_(0-44]": 1,
#    "pop_den_(44-140]": 1 - 0.75,
#    "pop_den_(140-310]": 1 - 0.9,
#    "pop_den_(310-400]": 1 - 0.95,
#    "pop_den_400+": 1 - 1,
#  },
#  "onshorewind": {
#    "pop_den_(0-7]": 1,
#    "pop_den_(7-13]": 1 - 0.75,
#    "pop_den_(13-18]": 1 - 0.9,
#    "pop_den_(18-20]": 1 - 0.95,
#    "pop_den_20+": 1 - 1,
#  }
#Max MW = KM2 x site boundary power density x max development density per area

#################################################################################################################################
# Read the input data
#################################################################################################################################
All_CPAs <- read.csv("All_potential_cpa_data.csv", header=TRUE)

# read in rawcap file, filter to only keep relevant rows
rawcap_ <- read.csv("../data/ira_final_20230117/capacity_y.csv")
rawcap_ <- subset(rawcap_, value >= 0.05)
rawcap <- rawcap_%>% filter(tech..outputs_group_aggregate %in% c("onshore wind", "offshore wind","solar")) %>%
  filter(!tech%in%c("existing_onshore wind|1", "existing_offshore wind|1", "existing_transmission-sited solar pv power plant_1","existing_rooftop pv - com","existing_rooftop pv - res","existing_rooftop pv - pro","existing_solar thermal with energy storage|1","rooftop solar - res"))%>%
  filter(!zone %in% c("alaska agg.", "hawaii agg."))%>%
  mutate(tech = as.character(tech))%>%
  mutate(net_cluster = substr(tech, nchar(tech), nchar(tech)))%>%
  dplyr::rename(tech_ = tech)%>%
  mutate(value_MW = value*1000)%>%
  mutate(tech..outputs_group_detailed = ifelse(grepl('fixed_0', tech_),"offshore wind fixed_np - nonpreferred",
                                               ifelse(grepl('fixed_1', tech_), "offshore wind fixed - preferred",
                                               ifelse(grepl('floating_0',tech_),"offshore wind floating - nonpreferred",
                                               ifelse(grepl('floating_1', tech_), "offshore wind floating - preferred",tech_)))))

# rename columns and cells to allow easier merge
names(rawcap)[names(rawcap)=="zone"] <- "region"
rawcap$tech[rawcap$tech..outputs_group_detailed %in% c("utilitypv_class1_1","utilitypv_class1_2","utilitypv_class1_3")] <- "solar"
rawcap$tech[rawcap$tech..outputs_group_detailed %in% c("landbasedwind_class4_1","landbasedwind_class4_2","landbasedwind_class4_3","landbasedwind_class4_4","landbasedwind_class4_5")] <- "onshore_wind"
rawcap$tech[rawcap$tech..outputs_group_detailed == "offshore wind floating - preferred"] <- "offshore_wind_floating"
rawcap$tech[rawcap$tech..outputs_group_detailed == "offshore wind fixed - preferred"] <- "offshore_wind_fixed"
rawcap$tech[rawcap$tech..outputs_group_detailed == "offshore wind floating - nonpreferred"] <- "offshore_wind_floating_np"
rawcap$tech[rawcap$tech..outputs_group_detailed == "offshore wind fixed_np - nonpreferred"] <- "offshore_wind_fixed_np"

# find all unique combinations of MGA iteration and policy
policy_list <- unique(rawcap$run.name)

# repeat loop for each policy
for(m in policy_list)
{
  
  # subset rawcap for only policy, find number of unique MGA iterations for
  # this policy
  curr_policy <- subset(rawcap, run.name == m)

  merged_iteration <- merge(curr_policy, All_CPAs)
  merged_iteration <- merged_iteration %>% mutate(name = paste0(tech,region, net_cluster))
  subset_df <- merged_iteration[0,]
  
  # runs for each unique case of capacity
  for (y in unique(merged_iteration$year)){
    curr_capframe_ <- merged_iteration[merged_iteration$year==as.integer(y),]
    print(y)
  for(i in unique(curr_capframe_$name))
  {
    #print(i)
    # this will store the current iteration's unique capacity dataframe
    curr_capframe <- curr_capframe_[curr_capframe_$name == i,]
    cap_sum = 0
    rowcount = length(curr_capframe$lcoe)
    curr_capframe<- curr_capframe %>% arrange(lcoe)
    # finds number of CPAs that will fit in capacity
    for (j in 1:length(curr_capframe$lcoe))
    {
      cap_sum = cap_sum + curr_capframe$cpa_mw[j]
      if(cap_sum>= curr_capframe$value_MW[j])
      {
        rowcount <- j
        break
      }
    }
    
    # transfers appropriate number of CPA_id rows to subset_df
    sliced_capframe <- slice_head(curr_capframe, n= rowcount)
    subset_df <- rbind(subset_df, sliced_capframe)
  }
  }
  # all final parq variables from each iteration are stored here
  write.csv(subset_df, paste0("Selected_CPAs_all_", m, ".csv"))

  
}
