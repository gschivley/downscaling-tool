from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Tuple, Union
import numpy as np
import rasterio
import rasterio.features
import geopandas as gpd
import pandas as pd
from tqdm import tqdm
import logging
from IPython import embed as IP
import os
import geopandas as gpd

CWD = Path.cwd()


transmission_data = pd.read_csv("transmission/interregional_line_cost_by_state_rio_emm.csv") 
RIO_results = pd.read_csv("../capacity_downscaling/data/ira_final_20230117/tx_capacity_y.csv") 
transmission_line_data = pd.read_csv("transmission/interregional_network_connections_rio_emm_updated.csv")

RIO_results_elec = RIO_results[(RIO_results['blend']=="electricity")].reset_index()
RIO_results_filtered = RIO_results_elec.drop(['zone', 'zone_to', 'unit', 'blend'], axis=1).groupby(['year', 'name', 'run name']).mean().reset_index()

transmission_line_data['path1'] = transmission_line_data['start_region'] + '||' + transmission_line_data['dest_region'] + '||electricity||1'
transmission_line_data['path2'] = transmission_line_data['dest_region'] + '||' + transmission_line_data['start_region'] + '||electricity||1'
transmission_data['path1'] = transmission_data['sink_to'] + '||' + transmission_data['sink_from'] + '||electricity||1'
transmission_data['path2'] = transmission_data['sink_from'] + '||' + transmission_data['sink_to'] + '||electricity||1'
#RIO_results_filtered = RIO_results_elec_avg[(RIO_results_elec_avg.name.isin(list(transmission_line_data['path1']) + list(transmission_line_data['path2'])))].reset_index()

transmission_results = pd.DataFrame()
#downscaling state-specific cost distribution
for index, row in RIO_results_filtered.iterrows():

	print(row['name'])

	transmission_results_ = transmission_data[(transmission_data['path1']==row['name']) | (transmission_data['path2']==row['name'])]
	transmission_results_['final_cost (overnight capex $)'] = transmission_results_['state_specific_cost'] * row['value']*1000
	transmission_results_['year'] = row['year']
	transmission_results_['run_name'] = row['run name']

	miles = transmission_data[(transmission_data['path1']==row['name']) | (transmission_data['path2']==row['name'])]['state_specific_dist']*0.621371

	transmission_results_['GW_mile'] = row['value'] * miles
	transmission_results_['name'] = row['name']
	
	transmission_results_['GW_state'] = row['value'] * transmission_results_['frac_total_capex']

	transmission_results = transmission_results.append(transmission_results_, ignore_index=True)

num = transmission_results._get_numeric_data()
num[num < 0] = 0
transmission_results.to_csv("transmission_cost_distribution_results_y.csv")

#downscaling the backbone and transmission network
backbone = pd.read_csv("transmission/intraregion_msa_connections_rio_emm.csv") 
backbone_results = pd.DataFrame()
for index, row in RIO_results_elec.iterrows(): 
	print(row['name'])

	backbone_from_ = backbone[backbone['region']==row['zone']]
	backbone_to_ = backbone[backbone['region']==row['zone_to']]
	backbone_results_ = backbone_from_.append(backbone_to_)
	backbone_results_['year'] = row['year']
	backbone_results_['run_name'] = row['run name']

	backbone_results = backbone_results.append(backbone_results_, ignore_index=True)

backbone_results.to_csv("backbone_selected_results_y.csv")

# #computing the length of the backbone lines
# intra_region_routes = pd.read_csv("../Creating_transmission/raw_data/all_intraregion_msa_connections_no_HVDC2/all_intraregion_msa_connections.csv") 

# shapefile = gpd.read_file("../Creating_transmission/raw_data/all_intraregion_msa_connections_no_HVDC2/all_intraregion_msa_connections_no_HVDC2.shp")

# shapefile['dist_km'] = 0
# for i in range(0,len(shapefile)):
# 	shapefile['dist_km'][i]=shapefile.iloc[i]['geometry'].length/1000 

# shapefile1 = shapefile.drop('geometry',1)
# shapefile1.to_csv("../Creating_transmission/raw_data/all_intraregion_msa_connections_no_HVDC2/all_intraregion_msa_connections_updated.csv")

#compute MW-miles for inter-regional transmission lines. Note: I am considering the total length include the regional backbone
 	
RIO_results_filtered['Total_GW_mile']=0
RIO_results_filtered['Total_mile']=0
RIO_results_filtered['Inter_regional_mile']=0
RIO_results_filtered['Inter_regional_GW_mile']=0
RIO_results_filtered['Backbone_GW_mile']=0
for i in range(0,len(RIO_results_filtered)): 
	# using Miles instead of km

	RIO_results_filtered['Inter_regional_mile'][i] = transmission_line_data[(transmission_line_data['path1']==RIO_results_filtered.iloc[i]['name']) | (transmission_line_data['path2']==RIO_results_filtered.iloc[i]['name'])]['distance_km']*0.621371
	RIO_results_filtered['Total_mile'][i] = transmission_line_data[(transmission_line_data['path1']==RIO_results_filtered.iloc[i]['name']) | (transmission_line_data['path2']==RIO_results_filtered.iloc[i]['name'])]['updated_dist']*0.621371
	
	RIO_results_filtered['Total_GW_mile'][i] = RIO_results_filtered['Total_mile'][i] * RIO_results_filtered['value'][i]
	RIO_results_filtered['Inter_regional_GW_mile'][i] = RIO_results_filtered['Inter_regional_mile'][i] * RIO_results_filtered['value'][i]
	RIO_results_filtered['Backbone_GW_mile'][i] = RIO_results_filtered['Total_GW_mile'][i] - RIO_results_filtered['Inter_regional_GW_mile'][i]

	RIO_results_filtered['Total_GW_mile'][RIO_results_filtered['Total_GW_mile'] < 0] = 0
	RIO_results_filtered['Inter_regional_GW_mile'][RIO_results_filtered['Inter_regional_GW_mile'] < 0] = 0
	RIO_results_filtered['Backbone_GW_mile'][RIO_results_filtered['Backbone_GW_mile'] < 0] = 0
	


RIO_results_filtered.to_csv("RIO_Transmission_GW_mile_y.csv")



