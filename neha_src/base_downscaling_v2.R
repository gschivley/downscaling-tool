setwd('/Users/npatankar/D/US_Policy/downscaling_Jan23/capacity_downscaling/plots/')
getwd()
#tools, install packages, type in the name of the packages 
library(ggplot2)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(RSQLite)
library(gridExtra)
library(dplyr)
library(reshape)
library(MASS)
library(reshape2)
library(GGally)
library(lattice)
library(data.table)
library(splitstackshape)
library(ggvis)
library(tidyr)
library(viridis)
library(rgdal)
library(sf)
library(arrow)
library(purrr)
library(png)
set.seed(4561)
library(extrafont)
loadfonts(device = "win")
rm(list = setdiff ( ls() , "") )

Raw_data_W <- read.csv(file="../lcoe/onshorewind_lcoe_rio_emm_metro_county.csv", header=TRUE, sep=",")
Raw_data_S <- read.csv(file="../lcoe/solar_lcoe_rio_emm_metro_county.csv", header=TRUE, sep=",")

Raw_data <- read.csv(file="../result_data/All_potential_cpa_data.csv", header=TRUE, sep=",")
#Raw_data_S <- subset(Raw_data, Raw_data$tech=="solar")
#Raw_data_W <- subset(Raw_data, Raw_data$tech=="onshore_wind")

Existing_S <- Raw_data_S %>% 
  dplyr::select(exFacil, plFacil, CPA_ID, latitude, longitude, Area, m_popden, m_HMI, incap, metro_region, lcoe,state) %>%
  mutate(ex = exFacil+plFacil)%>%
  filter(ex>=1) %>%
  mutate(tech = "solar")

# Location of existinf onshore wind sites
Existing_W <- Raw_data_W  %>% 
  dplyr::select(exFacil, plFacil, CPA_ID, latitude, longitude, Area, m_popden, m_HMI, incap,  metro_region, lcoe,state) %>%
  mutate(ex = exFacil+plFacil)%>%
  filter(ex>=1) %>%
  mutate(tech = "onshore_wind")
Ex_S <- data.frame(Existing_S %>% group_by(state) %>% summarize(Capacity_sum = sum(incap))%>% mutate(year = 2020, tech="solar"))
Ex_W <- data.frame(Existing_W %>% group_by(state) %>% summarize(Capacity_sum = sum(incap))%>% mutate(year = 2020, tech="onshore_wind"))

#read the shape files for the CPAs
solar_cpa_shp <- readOGR(dsn = file.path("../../Creating_transmission/raw_data/CandidateProjectAreas_WindAndSolar_20210623/CandidateProjectArea_SolarPV.shp"), stringsAsFactors = F)
solar_cpa_shp <- spTransform(solar_cpa_shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
solar_cpa_shp_existing <- solar_cpa_shp[(solar_cpa_shp$exFacil ==1) | (solar_cpa_shp$plFacil ==1),]

wind_cpa_shp <- readOGR(dsn = file.path("../../Creating_transmission/raw_data/CandidateProjectAreas_WindAndSolar_20210623/CandidateProjectArea_OnshoreWind.shp"), stringsAsFactors = F)
wind_cpa_shp <- spTransform(wind_cpa_shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
wind_cpa_shp_existing <- wind_cpa_shp[(wind_cpa_shp$exFacil ==1) | (wind_cpa_shp$plFacil==1),]

library(maps)
library(ggthemes)

gcounty <- map_data("county")
gusa <- map_data("state")
fipstab <- transmute(county.fips, fips, county = sub(":.*", "", polyname))
fipstab <- separate(fipstab, county, c("region", "subregion"), sep = ",")
gcounty <- left_join(gcounty, fipstab, c("region", "subregion"))



###########################################################################################
# Plot downscaling resutls for the Base case of all policies
###########################################################################################
library(geojsonio)
library(rgeos)
library(raster)
library(sp)
library(maptools)
library(broom)

inter_regional_routes <- read.csv("../../transmission_downscaling/transmission/interregional_network_connections_rio_emm_updated.csv")%>%
  mutate(route_name = paste(start_region, "_to_",dest_region))

inter_regional_results <- read.csv("../data/ira_final_20230117/tx_capacity_y.csv") %>% filter(blend=="electricity") %>%
  mutate(route_name = paste(zone, "_to_",zone_to))%>%
  filter(route_name %in% inter_regional_routes$route_name)

backbone <- read.csv(file="../../transmission_downscaling/transmission/intraregion_msa_connections_rio_emm.csv", header=TRUE, sep=",")

shp_msa <- readOGR(dsn = file.path("../../transmission_downscaling/transmission/intraregion_msa_connections_rio_emm/intraregion_msa_connections_rio_emm.shp"), stringsAsFactors = F)
shp_proj <- spTransform(shp_msa, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

infinite_sinks <- read.csv(file="../../Creating_transmission/derived_data/infinite_sink_msas.csv", header=TRUE, sep=",")
spdf <- geojson_read("../../Creating_transmission/derived_data/substations_urban_convex_20210603.geojson",  what = "sp")
spdf_ <- spTransform(spdf, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
sinks = spdf_[spdf_$CBSAFP %in% unique(infinite_sinks$msa_id),]

####################################################################################################################
# For plotting selected substation to msa
####################################################################################################################
shp_substation <- readOGR(dsn = file.path("../data/substation_linestrings/substation_linestrings.shp"), stringsAsFactors = F)
shp_substation <- spTransform(shp_substation, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

solar_path <- read.csv("../lcoe/solar_lcoe_rio_emm_metro_county.csv")%>% dplyr::select(tech, metro_id, end_msa_substation_1,	end_msa_substation_2,	end_msa_substation_3,	end_msa_substation_4)
wind_path <- read.csv("../lcoe/onshorewind_lcoe_rio_emm_metro_county.csv")%>% dplyr::select(tech, metro_id, end_msa_substation_1,	end_msa_substation_2,	end_msa_substation_3,	end_msa_substation_4)
offwind_path <- read.csv("../lcoe/offshorewind_lcoe_rio_emm_metro_county.csv") %>% dplyr::select(tech, metro_id, end_msa_substation_1,	end_msa_substation_2,	end_msa_substation_3,	end_msa_substation_4)

test <- rbind(solar_path, wind_path, offwind_path)
test_ <- subset(test, !is.na(test$end_msa_substation_2))
test_1 <- test_ %>% dplyr::select(metro_id, end_msa_substation_1, end_msa_substation_2) %>% dplyr::rename(start_substation = end_msa_substation_1, end_substation = end_msa_substation_2)
test_2 <- test_ %>% dplyr::select(metro_id, end_msa_substation_2, end_msa_substation_3) %>% dplyr::rename(start_substation = end_msa_substation_2, end_substation = end_msa_substation_3) %>% na.omit
test_3 <- test_ %>% dplyr::select(metro_id, end_msa_substation_3, end_msa_substation_4) %>% dplyr::rename(start_substation = end_msa_substation_3, end_substation = end_msa_substation_4) %>% na.omit

test_sub_sub <- unique(rbind(test_1, test_2, test_3))
####################################################################################################################
# For plotting selected inter and intra-regional transmission
####################################################################################################################
shp_inter_region <- readOGR(dsn = file.path("../../transmission_downscaling/transmission/interregional_connections_rio_emm/interregional_connections_rio_emm.shp"), stringsAsFactors = F)
shp_inter_region <- spTransform(shp_inter_region, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

####################################################################################################################
#plotting convex hull for infinite sink MSAs
####################################################################################################################
rawcap <- read.csv("../data/ira_final_20230117/capacity_y.csv")

data_df <- as.data.frame(spdf_)
data_sinks <- as.data.frame(sinks)
hull <- data_sinks %>%
  group_by(CBSAFP) %>%
  slice(chull(coords.x1, coords.x2))

centroids <- read.csv(file="../data/centroids.csv", header=TRUE, sep=",")
state_cost_distribution_renewable <- data.frame()

  
for (P in unique(rawcap$run.name)) {
  print (P)
  Raw_data_sel <- read.csv(file=paste("../result_data/Selected_CPAs_all_",P,".csv",sep = ""), header=TRUE, sep=",")
  inter_regional_results_temp <- inter_regional_results %>% filter(as.character(run.name)==as.character(P))

for (y in c(2022,2024,2026,2028,2030,2032, 2035,2040,2050)) {
  print(y)
  
  plot <- ggplot() +
    geom_polygon(aes(long, lat, group = group),fill = "white", data = gusa, color = "black", size=0.5)+
    geom_polygon(data = hull, aes(coords.x1, coords.x2 ,group=CBSAFP),fill = "white", color="black", alpha = 0.7, size=0.05)+
    with(centroids, annotate(geom="text", x = longitude, y=latitude, label = state, size = 5, color="grey"))+ 
    coord_map("albers",lat0=39, lat1=45, xlim = c(-124, -75),ylim = c(25, 60))+ 
    ggthemes::theme_map()+
    ggtitle(paste(y,": ",P," scenario",sep = ""))+
    theme(plot.title = element_text(size = 40, face = "bold",hjust = 0.5, vjust=0.1))+
    annotate("point",y = 55, x = -124,colour = "orange", size=5)+
    annotate("point",y = 57, x = -104,colour = "red", size=5)+
    annotate("point",y = 56.3, x = -84,colour = "violet", size=5)+
    annotate("point",y = 53, x = -123,colour = "lightblue", size=5)+
    annotate("point",y = 55, x = -103.8,colour = "royalblue4", size=5)+
    annotate("point",y = 54.3, x = -84.7,colour = "purple", size=5)+
    annotate("segment",y = 50.5, x = -123,yend=50.85,xend=-121,colour = "dark red", linewidth=2)+
    annotate("segment",y = 52.6, x = -104.5,yend=52.7,xend=-102.5,colour = "indianred", linewidth=2)+
    annotate("segment",y = 52.2, x = -86.3,yend=52,xend=-84.3,colour = "salmon", linewidth=2)+
    annotate(geom="text", x = -122, y =55.5,size=6, label = str_wrap("Existing Utility-scale Solar", width=35), hjust = "left")+
    annotate(geom="text", x = -102, y =57,size=6, label = str_wrap("New Utility-Scale Solar", width=35), hjust = "left")+
    annotate(geom="text", x = -82, y =56,size=6, label = str_wrap("New Offshore Wind - Preferred", width=35), hjust = "left")+
    annotate(geom="text", x = -121, y =53.5,size=6, label = str_wrap("Existing Onshore Wind", width=35), hjust = "left")+
    annotate(geom="text", x = -101.8, y =55,size=6, label = str_wrap("New Onshore Wind", width=35), hjust = "left")+
    annotate(geom="text", x = -82.7, y =54,size=6, label = str_wrap("New Offshore Wind - Non-preferred", width=35), hjust = "left")+
    annotate(geom="text", x = -119.8, y =51,size=6, label = str_wrap("Inter-regional transmission (>500 kV HVDC)", width=30), hjust = "left")+
    annotate(geom="text", x = -101.5, y =52.5,size=6, label = str_wrap("Intra-region transmission (<500 kV,>240 kV)", width=30), hjust = "left")+
    annotate(geom="text", x = -83.5, y =51.7,size=6, label = str_wrap("CPA-substation to metro-substation (240 kV)", width=35), hjust = "left")
  
  # creating a dataframe for state-specific cost distribution
  Raw_data_sel_S_ <- subset(Raw_data_sel, Raw_data_sel$year==y & Raw_data_sel$tech=="solar")
  Raw_data_sel_W_ <- subset(Raw_data_sel, Raw_data_sel$year==y & Raw_data_sel$tech=="onshore_wind")
  Raw_data_sel_OWFL_temp<- subset(Raw_data_sel, Raw_data_sel$year==y & Raw_data_sel$tech=="offshore_wind_floating")
  Raw_data_sel_OWFL <- Raw_data_sel_OWFL_temp[!(is.na(Raw_data_sel_OWFL_temp$state) | Raw_data_sel_OWFL_temp$state==""), ]
  Raw_data_sel_OWF<- subset(Raw_data_sel, Raw_data_sel$year==y & Raw_data_sel$tech=="offshore_wind_fixed")
  Raw_data_sel_OWFL_np <- subset(Raw_data_sel, Raw_data_sel$year==y & Raw_data_sel$tech=="offshore_wind_floating_np")
  Raw_data_sel_OWF_np <- subset(Raw_data_sel, Raw_data_sel$year==y & Raw_data_sel$tech=="offshore_wind_fixed_np")
  
  state_cost_distribution_renewable <- bind_rows(state_cost_distribution_renewable,
                                                 Raw_data_sel_S_ %>% group_by(state) %>% summarize(Capacity_sum = sum(cpa_mw)) %>% mutate(tech = "solar", year = y, Policy = P),
                                                 Raw_data_sel_W_ %>% group_by(state) %>% summarize(Capacity_sum = sum(cpa_mw)) %>% mutate(tech = "onshore_wind", year = y, Policy = P),
                                                 Raw_data_sel_OWFL %>% group_by(state) %>% summarize(Capacity_sum = sum(cpa_mw)) %>% mutate(tech = "offshore_wind_floating", year = y, Policy = P),
                                                 Raw_data_sel_OWF %>% group_by(state) %>% summarize(Capacity_sum = sum(cpa_mw)) %>% mutate(tech = "offshore_wind_fixed", year = y, Policy = P),
                                                 Raw_data_sel_OWFL_np %>% group_by(state) %>% summarize(Capacity_sum = sum(cpa_mw)) %>% mutate(tech = "offshore_wind_floating_np", year = y, Policy = P),
                                                 Raw_data_sel_OWF_np %>% group_by(state) %>% summarize(Capacity_sum = sum(cpa_mw)) %>% mutate(tech = "offshore_wind_fixed_np", year = y, Policy = P))
  
  #filtering the selected CPAs for a given year
  km <- kmeans(cbind(Raw_data_sel_S_$latitude, Raw_data_sel_S_$longitude), centers = length(Raw_data_sel_S_[,1])/1.1)
  Raw_data_sel_S_$cluster <- km$cluster
  Raw_data_sel_S_0 <- subset(Raw_data_sel_S_, Raw_data_sel_S_$m_popden <= 5)
  Raw_data_sel_S_0_c <- Raw_data_sel_S_0 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.2)) %>% mutate(unique_id=1:NROW(.))
  Raw_data_sel_S_0 <- Raw_data_sel_S_0 %>% group_by(cluster)  %>% right_join(Raw_data_sel_S_0_c)
  Raw_data_sel_S_1 <-subset(Raw_data_sel_S_, Raw_data_sel_S_$m_popden > 5 & Raw_data_sel_S_$m_popden <=10)
  Raw_data_sel_S_1_c <- Raw_data_sel_S_1 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.1)) %>% mutate(unique_id=1:NROW(.))
  Raw_data_sel_S_1 <- Raw_data_sel_S_1 %>% group_by(cluster)  %>% right_join(Raw_data_sel_S_1_c)
  Raw_data_sel_S_2 <-subset(Raw_data_sel_S_, Raw_data_sel_S_$m_popden > 10 & Raw_data_sel_S_$m_popden <= 30)
  if (y!=2020){
    Raw_data_sel_S_2_c <- Raw_data_sel_S_2 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.05)) %>% mutate(unique_id=1:NROW(.))
    Raw_data_sel_S_2 <- Raw_data_sel_S_2 %>% group_by(cluster)  %>% right_join(Raw_data_sel_S_2_c)
    Raw_data_sel_S_3 <-subset(Raw_data_sel_S_, Raw_data_sel_S_$m_popden > 30 & Raw_data_sel_S_$m_popden <= 40)
    Raw_data_sel_S_3_c <- Raw_data_sel_S_3 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.025)) %>% mutate(unique_id=1:NROW(.))
    Raw_data_sel_S_3 <- Raw_data_sel_S_3 %>% group_by(cluster)  %>% right_join(Raw_data_sel_S_3_c)
    Raw_data_sel_S <- rbind(Raw_data_sel_S_0, Raw_data_sel_S_1,Raw_data_sel_S_2,Raw_data_sel_S_3)
  }
  if (y==2020){
    Raw_data_sel_S <- rbind(Raw_data_sel_S_0, Raw_data_sel_S_1)
  }




  #Raw_data_sel_W_ <- Raw_data_W %>% filter(m_HMI<0.8, m_popden<=40)
  km <- kmeans(cbind(Raw_data_sel_W_$latitude, Raw_data_sel_W_$longitude), centers = length(Raw_data_sel_W_[,1])/6)
  Raw_data_sel_W_$cluster <- km$cluster
  Raw_data_sel_W_0 <- subset(Raw_data_sel_W_, Raw_data_sel_W_$m_popden <= 1)
  Raw_data_sel_W_0_c <- Raw_data_sel_W_0 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.6)) %>% mutate(unique_id=1:NROW(.))
  Raw_data_sel_W_0 <- Raw_data_sel_W_0 %>% group_by(cluster)  %>% right_join(Raw_data_sel_W_0_c)
  Raw_data_sel_W_1 <-subset(Raw_data_sel_W_, Raw_data_sel_W_$m_popden > 1 & Raw_data_sel_W_$m_popden <=3)
  Raw_data_sel_W_1_c <- Raw_data_sel_W_1 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.2)) %>% mutate(unique_id=1:NROW(.))
  Raw_data_sel_W_1 <- Raw_data_sel_W_1 %>% group_by(cluster)  %>% right_join(Raw_data_sel_W_1_c)
  if (y!=2020){
  Raw_data_sel_W_2 <-subset(Raw_data_sel_W_, Raw_data_sel_W_$m_popden > 3 & Raw_data_sel_W_$m_popden <=5)
  Raw_data_sel_W_2_c <- Raw_data_sel_W_2 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.1)) %>% mutate(unique_id=1:NROW(.))
  Raw_data_sel_W_2 <- Raw_data_sel_W_2 %>% group_by(cluster)  %>% right_join(Raw_data_sel_W_2_c)
  Raw_data_sel_W_3 <-subset(Raw_data_sel_W_, Raw_data_sel_W_$m_popden > 5 & Raw_data_sel_W_$m_popden <=10)
  Raw_data_sel_W_3_c <- Raw_data_sel_W_3 %>% group_by(cluster) %>% dplyr::summarise() %>% sample_n(ceiling(n()*0.05)) %>% mutate(unique_id=1:NROW(.))
  Raw_data_sel_W_3 <- Raw_data_sel_W_3 %>% group_by(cluster)  %>% right_join(Raw_data_sel_W_3_c)
  Raw_data_sel_W <- rbind(Raw_data_sel_W_0, Raw_data_sel_W_1,Raw_data_sel_W_2,Raw_data_sel_W_3)
  }
  if (y==2020){
    Raw_data_sel_W <- rbind(Raw_data_sel_W_0, Raw_data_sel_W_1)
  }

  #plot selected transmission-substation to MSA-substation lines
  Raw_data_sel_y <- bind_rows(Raw_data_sel_S, Raw_data_sel_W, Raw_data_sel_OWFL, Raw_data_sel_OWF,Raw_data_sel_OWFL_np, Raw_data_sel_OWF_np)
  test_msa <- subset(test_sub_sub, as.character(test_sub_sub$metro_id) %in% unique(Raw_data_sel_y$metro_id))
  shp_substation_ <- subset(shp_substation, shp_substation$start_id %in% test_msa$start_substation)
  shp_substation_ <- subset(shp_substation_, shp_substation_$dest_id %in% test_msa$end_substation)
  plot <- plot+geom_path(data=shp_substation_, aes(x=long, y=lat, group=group), color="salmon", size=0.06, inherit.aes = FALSE, alpha=0.7)

  #plotting inter-regional transmission line for a selected year and scenario
  inter_regional_results_ <- inter_regional_results_temp %>% filter(year==y) %>% right_join(inter_regional_routes)

  for (row in 1:nrow(inter_regional_results_)){
    #print(inter_regional_results_[row,]$route_name)
    shp_inter_region_ <- subset(shp_inter_region, (shp_inter_region$start_id == inter_regional_results_[row,]$start_id & shp_inter_region$dest_id == inter_regional_results_[row,]$dest_id))
    plot <- plot + geom_path(data=shp_inter_region_, aes(x=long, y=lat, group=group), color="Dark Red", size=inter_regional_results_[row,]$value/20, inherit.aes = FALSE)

    #ploting the corresponding backbone
    backbone_ <- backbone %>% filter(region %in% c(as.character(inter_regional_results_[row,]$zone),as.character(inter_regional_results_[row,]$zone_to)))%>%
      filter(mst_conn == 'True')%>%
      distinct(start_msa, dest_msa, interconnect_cost_mw ,region)
    if (dim(backbone_)[1]!=0){
      for (b in 1:length(backbone_[,1])){
        backbone_temp <- shp_proj[shp_proj$msa_1==as.character(backbone_[b,]$start_msa) & shp_proj$msa_2==as.character(backbone_[b,]$dest_msa),]
        plot <- plot+ geom_path(data = backbone_temp, aes(x = long, y = lat, group = group), colour = "indianred", size=inter_regional_results_[row,]$value/40, inherit.aes = FALSE)
      }
    }
  }

  #plotting the selected CPAs
  solar_cpa_shp_tmp <- solar_cpa_shp[as.integer(solar_cpa_shp$CPA_ID) %in% unique(Raw_data_sel_S$CPA_ID),]
  wind_cpa_shp_tmp <- wind_cpa_shp[as.integer(wind_cpa_shp$CPA_ID) %in% unique(Raw_data_sel_W$CPA_ID),]
  plot <- plot +
    geom_polygon(data = wind_cpa_shp_existing, aes(x=long, y=lat, group=group), fill="lightblue", inherit.aes = FALSE) +
    geom_polygon(data = solar_cpa_shp_existing, aes(x=long, y=lat, group=group), fill="orange", inherit.aes = FALSE) +
    geom_polygon(data = wind_cpa_shp_tmp, aes(x=long, y=lat, group=group), fill="royalblue4", inherit.aes = FALSE) +
    geom_polygon(data = solar_cpa_shp_tmp, aes(x=long, y=lat, group=group), fill="red", inherit.aes = FALSE) +
    geom_point(data = Raw_data_sel_OWFL, aes(longitude, latitude), size=Raw_data_sel_OWFL$Area/1e8, color="violet") +
    geom_point(data = Raw_data_sel_OWF, aes(longitude, latitude), size=Raw_data_sel_OWF$Area/1e8, color="violet") +
    geom_point(data = Raw_data_sel_OWFL_np, aes(longitude, latitude), size=Raw_data_sel_OWFL_np$Area/1e8, color="purple") +
    geom_point(data = Raw_data_sel_OWF_np, aes(longitude, latitude), size=Raw_data_sel_OWF_np$Area/1e8, color="purple") +
    theme(text = element_text(family = "Futura"),
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank())

  ggsave(paste(y,"/",P,".png",sep = ""), width=25, height=25, dpi=500)

}
}
state_cost_distribution_renewable_ <- state_cost_distribution_renewable %>% mutate(year = paste('y',year,sep=''))%>% spread(year, Capacity_sum) %>%
  replace(is.na(.), 0)%>%
  mutate(Delta_2024 = y2024 - y2022, Delta_2026 = y2026-y2024, Delta_2028 = y2028-y2026,
         Delta_2030 = y2030 - y2028,Delta_2032 = y2032 - y2030, Delta_2035 = y2035-y2032, Delta_2040 = y2040-y2035,
         Delta_2050 = y2050 - y2040)%>%
  dplyr::select(state, tech, Policy,y2022,Delta_2024,Delta_2026,Delta_2028,Delta_2030,Delta_2032,Delta_2035,Delta_2040,Delta_2050)%>%
  dplyr::rename('2022' = y2022, '2024' = Delta_2024,'2026' = Delta_2026,'2028'=Delta_2028,
                '2030' = Delta_2030,'2032' = Delta_2032,'2035' = Delta_2035,'2040'=Delta_2040,'2050'=Delta_2050)

state_cost_distribution_renewable_final <- melt(state_cost_distribution_renewable_)%>%dplyr::rename(year=variable, Capacity_sum=value)
#savingt the CSV file for cost distribution
write.csv(state_cost_distribution_renewable_final,"../result_data/renewable_cost_distribution_delta.csv")
write.csv(state_cost_distribution_renewable,"../result_data/renewable_cost_distribution.csv")
