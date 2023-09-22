setwd('/Users/npatankar/D/US_Policy/downscaling_Jan23/capacity_downscaling/result_data')
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(readr)


# read in csv files
csv_solarpv <- read.csv("../metadata/rio-emm_utilitypv_(metro_county_)_False_site_cluster.csv", header=TRUE)
csv_solarpv$tech <- "solar" 
csv_offwind_fixed <- read.csv("../metadata/rio-emm_offshorewind_(metro_county_)_False_fixed_True_site_cluster.csv", header=TRUE)
csv_offwind_fixed$tech <- "offshore_wind_fixed"
csv_offwind_fixed_np <- read.csv("../metadata/rio-emm_offshorewind_(metro_county_)_False_fixed_False_site_cluster.csv", header=TRUE)
csv_offwind_fixed_np$tech <- "offshore_wind_fixed_np"
csv_offwind_floating <- read.csv("../metadata/rio-emm_offshorewind_(metro_county_)_False_floating_True_site_cluster.csv", header=TRUE)
csv_offwind_floating$tech <- "offshore_wind_floating"
csv_offwind_floating_np <- read.csv("../metadata/rio-emm_offshorewind_(metro_county_)_False_floating_False_site_cluster.csv", header=TRUE)
csv_offwind_floating_np$tech <- "offshore_wind_floating_np"
csv_wind <- read.csv("../metadata/rio-emm_landbasedwind_(metro_county_)_False_site_cluster.csv", header=TRUE)
csv_wind$tech <- "onshore_wind"

csv <- rbind(csv_solarpv, csv_offwind_fixed, csv_offwind_floating,csv_offwind_fixed_np, csv_offwind_floating_np, csv_wind)


# read in cluster ID/region linking csv
sample_sites <- read.csv("../data/TECH_CAPITAL_COST.csv", header=TRUE) %>%
  dplyr::select(region, name, cluster, zone, ids, resource_type) %>% 
  filter(resource_type %in% c("solar", "onshore_wind", "offshore_wind_floating", "offshore_wind_fixed","offshore_wind_floating_np", "offshore_wind_fixed_np")) %>%
  filter(ids !=0) %>%
  dplyr::rename(technology = resource_type)
sample_sites$technology <- as.character(sample_sites$technology)

techs <- c("solar", "onshore_wind", "offshore_wind_floating", "offshore_wind_fixed","offshore_wind_floating_np", "offshore_wind_fixed_np")

sample_scenarios <- data.frame()
# repeat loop for each policy
for(m in techs)
{
  
    # subsets Cluster_ids csv file based on the technology
    sample_sites_ <- subset(sample_sites, sample_sites$technology == m)
    
    # stores a list of cluster IDs in order
    cpa_id_list <- sample_sites_$ids
    cpa_id_list <- as.character(cpa_id_list)
    cpa_id_new<- gsub("[^0-9,-]", "", cpa_id_list)
    cpas <- strsplit(cpa_id_new, ",")
    
    # will sort through csv to choose only those rows that are linked to 
    # current policy and iteration
    relevant_csv_1 <- subset(csv, csv$tech==m)
    relevant_csv <- subset(csv, csv$tech==m)
    relevant_csv$region <- ""
    relevant_csv$net_cluster <- 0
    relevant_csv <- relevant_csv[0,]
    
    # iterates through each merged folder row, finding associated cpa_id
    # rows. These rows are extracted and built up into relevant_csv.
    for(i in 1:length(cpas))
    {
      temp <- relevant_csv_1[relevant_csv_1$id %in% cpas[[i]],]
      temp$region <- sample_sites_$region[i]
      temp$net_cluster <- sample_sites_$cluster[i]
      row.names(temp) <- NULL
      relevant_csv <- rbind(relevant_csv, temp)
    }
    
    # this will store the final csv files, after the CPA site selection is done
    sample_scenarios <- rbind(sample_scenarios, relevant_csv)

}
sample_scenarios <- sample_scenarios %>% dplyr::rename(CPA_ID = cpa_id)
write.csv(sample_scenarios, "All_potential_cpa_ids.csv")

####################################################################################################
# This section is for extracting area, location and capacity of the selected CPA
####################################################################################################

solar_lcoe <- read.csv("../lcoe/solar_lcoe_rio_emm_metro_county.csv")
solar_lcoe$tech <- "solar"
solar_lcoe <- solar_lcoe %>% dplyr::select(CPA_ID, lcoe, tech, Area, latitude, longitude,
              interconnect_annuity, exFacil, plFacil, metro_id, cpa_mw, m_popden, m_HMI, state)

wind_lcoe <- read.csv("../lcoe/onshorewind_lcoe_rio_emm_metro_county.csv")
wind_lcoe$tech <- "onshore_wind"
wind_lcoe <- wind_lcoe %>% dplyr::select(CPA_ID, lcoe, tech, Area, latitude, longitude,
            interconnect_annuity, exFacil, plFacil, metro_id, cpa_mw, m_popden, m_HMI, state)

offwind_lcoe_floating <- read.csv("../lcoe/offshorewind_lcoe_rio_emm_metro_county.csv") %>%
  filter(turbine_type == "floating", prefSite==1)
offwind_lcoe_floating$tech <- "offshore_wind_floating"
offwind_lcoe_floating <- offwind_lcoe_floating %>% dplyr::select(CPA_ID, lcoe, tech, Area, latitude, longitude,interconnect_annuity,  metro_id, cpa_mw, state)
offwind_lcoe_floating$exFacil <- 100
offwind_lcoe_floating$plFacil <- 100
offwind_lcoe_floating$m_popden <- 0
offwind_lcoe_floating$m_HMI <- 0

offwind_lcoe_floating_np <- read.csv("../lcoe/offshorewind_lcoe_rio_emm_metro_county.csv") %>%
  filter(turbine_type == "floating", prefSite==0)
offwind_lcoe_floating_np$tech <- "offshore_wind_floating_np"
offwind_lcoe_floating_np <- offwind_lcoe_floating_np %>% dplyr::select(CPA_ID, lcoe, tech, Area, latitude, longitude,interconnect_annuity,  metro_id, cpa_mw, state)
offwind_lcoe_floating_np$exFacil <- 100
offwind_lcoe_floating_np$plFacil <- 100
offwind_lcoe_floating_np$m_popden <- 0
offwind_lcoe_floating_np$m_HMI <- 0

offwind_lcoe_fixed <- read.csv("../lcoe/offshorewind_lcoe_rio_emm_metro_county.csv") %>%
  filter(turbine_type == "fixed", prefSite==1)
offwind_lcoe_fixed$tech <- "offshore_wind_fixed"
offwind_lcoe_fixed <- offwind_lcoe_fixed %>% dplyr::select(CPA_ID, lcoe, tech, Area, latitude, longitude,interconnect_annuity,  metro_id, cpa_mw, state)
offwind_lcoe_fixed$exFacil <- 100
offwind_lcoe_fixed$plFacil <- 100
offwind_lcoe_fixed$m_popden <- 0
offwind_lcoe_fixed$m_HMI <- 0

offwind_lcoe_fixed_np <- read.csv("../lcoe/offshorewind_lcoe_rio_emm_metro_county.csv") %>%
  filter(turbine_type == "fixed", prefSite==0)
offwind_lcoe_fixed_np$tech <- "offshore_wind_fixed_np"
offwind_lcoe_fixed_np <- offwind_lcoe_fixed_np %>% dplyr::select(CPA_ID, lcoe, tech, Area, latitude, longitude,interconnect_annuity, metro_id, cpa_mw, state)
offwind_lcoe_fixed_np$exFacil <- 100
offwind_lcoe_fixed_np$plFacil <- 100
offwind_lcoe_fixed_np$m_popden <- 0
offwind_lcoe_fixed_np$m_HMI <- 0

lcoe_all <- bind_rows(solar_lcoe, wind_lcoe, offwind_lcoe_floating, offwind_lcoe_fixed,offwind_lcoe_floating_np, offwind_lcoe_fixed_np)
#sample_scenarios <- sample_scenarios %>% dplyr::rename(CPA_ID = cpa_id)

final_cpa_1 <- merge(sample_scenarios, lcoe_all)

write.csv(final_cpa_1, paste0("All_potential_cpa_data.csv"))

