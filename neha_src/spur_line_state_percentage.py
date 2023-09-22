from shapely.strtree import STRtree
from pathlib import Path
from typing import Dict, List, Tuple, Union
import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from IPython import embed as IP
import pickle
import numbers

CWD = Path.cwd()
DATA = CWD / "data" / "cpa_linestrings"

def cpa_transmission_stats_by_state(
    start_ids: List[int],
    cpa_ls_dict,
    cpa_route_dict,
    assigned_df: pd.DataFrame,
    resource: str = None,
):

    v_classes = [
        "230kV_single",
        "230kV_double",
        "345kV_single",
        "345kV_double",
        "500kV_single",
        "500kV_double",
        "500kV_HVDC",
    ]


    if resource:
        assigned_df = assigned_df.query("tech == @resource")
    assigned_df = assigned_df.set_index("CPA_ID")

    state_gdf = gpd.read_file(DATA / "cb_2019_us_state_5m")
    state_gdf = state_gdf.to_crs("EPSG:5070")
    states = state_gdf.NAME.unique()

    # Store state name in gdf index as polygon object attributes for use in STRtree
    state_geoms = []
    state_geoms_temp = []
    for r in state_gdf.itertuples():
        poly = gpd.GeoDataFrame(geometry=[r.geometry])
        poly['state_name']=r.NAME
        poly['index']=r.Index
        state_geoms_temp.append(poly)
        state_geoms.append(poly.geometry.values[0])
    
    tree = STRtree(state_geoms)
    results_dict = {}
    summary_dict = {}
    for state in states:
        results_dict[state] = {}
        summary_dict[state] = {}
        for v in v_classes:
            summary_dict[state][v] = {"dist_km": 0, "capex": 0, "capacity_mw":0, "capex_km_mw":0, "route_mw_km":0}
    for start in tqdm(start_ids):
        IP()
        cpa_mw = np.array(assigned_df.loc[start, "cpa_mw"]).sum()
        tech = assigned_df['tech'].unique()[0]
        ls_list = cpa_ls_dict[start]
        route_list = cpa_route_dict[start]
        #print(start)
        #IP()
        for i, ls in enumerate(ls_list):
            if isinstance(route_list[i]['distance_km'], numbers.Number)==False:
                route_capex_km = 0
            elif route_list[i]["distance_km"] ==0 or ls=='error':
                route_capex_km = 0
            else:
                #print(i, start)
                capex_mw = np.array(route_list[i]['capex_mw'])
                distance_km = np.array(route_list[i]['distance_km'])
                route_capex_km = (capex_mw / distance_km).sum()
                route_capex_mw_km = (route_capex_km / cpa_mw).sum()
                #print(route_list[i])
                intersect_states = tree.query(ls)
                for s in intersect_states:
                    s_ls = tree.geometries[s].intersection(ls)
                    if i==0 and tech == 'offshorewind':
                        dist_km = distance_km
                    else:
                        dist_km = s_ls.length / 1000
                    cpa_mw_temp = cpa_mw * dist_km/distance_km
                    s_capex = dist_km * route_capex_km * cpa_mw_temp
                    v_class = route_list[i]["v_class"]
                    route_mw_km = (cpa_mw_temp * dist_km).sum()
                    if start not in results_dict[state_geoms_temp[s].state_name[0]]: #s.state_name]:
                        results_dict[state_geoms_temp[s].state_name[0]][start] = [
                            {
                                "v_class": v_class,
                                "capex": s_capex,
                                "capacity_mw": cpa_mw_temp,
                                "capex_km_mw":route_capex_km,
                                "dist_km": dist_km,
                                "linestring": s_ls,
                                "route_mw_km": route_mw_km,
                            }
                        ]
                    else:
                        results_dict[state_geoms_temp[s].state_name[0]][start].append(
                            {
                                "v_class": v_class,
                                "capex": s_capex,
                                "capacity_mw": cpa_mw_temp,
                                "capex_km_mw":route_capex_km,
                                "dist_km": dist_km,
                                "linestring": s_ls,
                                "route_mw_km": route_mw_km,
                            }
                        )
                    #IP()
                    if isinstance(v_class, str):
                        summary_dict[state_geoms_temp[s].state_name[0]][v_class]["dist_km"] += dist_km
                        summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capex"] += s_capex
                        summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capacity_mw"] += cpa_mw_temp
                        summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capex_km_mw"] = (summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capex_km_mw"] + route_capex_km)/2
                        summary_dict[state_geoms_temp[s].state_name[0]][v_class]["route_mw_km"] += route_mw_km
                    else:
                        for i in range(0,len(v_class)):
                            summary_dict[state_geoms_temp[s].state_name[0]][v_class[1]]["dist_km"] += dist_km
                            summary_dict[state_geoms_temp[s].state_name[0]][v_class[1]]["capex"] += s_capex
                            summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capacity_mw"] += cpa_mw_temp
                            summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capex_km_mw"] = (summary_dict[state_geoms_temp[s].state_name[0]][v_class]["capex_km_mw"] + route_capex_km)/2
                            summary_dict[state_geoms_temp[s].state_name[0]][v_class]["route_mw_km"] += route_mw_km
                IP()

    for k in results_dict.copy():
        if not results_dict[k]:
            results_dict.pop(k)
    for k in summary_dict.copy():
        for v_class in v_classes:
            if summary_dict[k][v_class]["capex"] == 0:
                summary_dict[k].pop(v_class)
        if not summary_dict[k]:
            summary_dict.pop(k)

    return results_dict, summary_dict

capacity_results = pd.read_csv("data/ira_final_20230117/capacity_delta_y.csv")
Policies = list(capacity_results['run name'].unique())

# concatenating selected CPAs
selected_cpas = pd.DataFrame()
for p in Policies:
    selected_cpas_temp = pd.read_csv("result_data/Selected_CPAs_all_" + p + ".csv")
    selected_cpas_temp['policy'] = p
    selected_cpas = selected_cpas.append(selected_cpas_temp)

# # selected_cpas = pd.read_csv("result_data/Selected_CPAs_all_core_copy.csv")
# # selected_cpas['policy'] = "core"
# # calculating solar spur line cost distribution
# selected_solar_cpas = selected_cpas[(selected_cpas['tech']=='solar') & (~selected_cpas['CPA_ID'].isin([395071, 405135, 377240, 395763]))]
# all_solar_cpas = pd.read_csv("lcoe/solar_lcoe_rio_emm_metro_county.csv")
# solar_ls_dict = pickle.load(open(DATA / "REPEAT_pkl_files/solar_ls_dict.pkl", "rb"))
# solar_route_dict = pd.read_pickle(DATA / "REPEAT_pkl_files/solar_route_dict.pkl")

# state_distribution_solar = pd.DataFrame()
# for year in selected_solar_cpas['year'].unique():
#     for p in Policies:
#         selected_solar_cpas_temp = selected_solar_cpas[(selected_solar_cpas['year']==year) & (selected_solar_cpas['policy']==p)]

#         results_dict, summary_dict = cpa_transmission_stats_by_state(list(selected_solar_cpas_temp['CPA_ID'].unique()), solar_ls_dict, solar_route_dict, all_solar_cpas)
#         A = pd.DataFrame.from_dict({(i,j): summary_dict[i][j] for i in summary_dict.keys() for j in summary_dict[i].keys()}, orient='index')
#         A['year'] = year
#         A['policy'] = p
#         state_distribution_solar = state_distribution_solar.append(A)
# state_distribution_solar['tech'] = 'solar'


# # calculating wind spur line cost distribution
# selected_wind_cpas = selected_cpas[(selected_cpas['tech']=='onshore_wind')]# & (~selected_cpas['CPA_ID'].isin([395071, 405135, 377240, 395763]))]
# all_wind_cpas = pd.read_csv("lcoe/onshorewind_lcoe_rio_emm_metro_county.csv")
# wind_ls_dict = pickle.load(open(DATA / "REPEAT_pkl_files/onshorewind_ls_dict.pkl", "rb"))
# wind_route_dict = pd.read_pickle(DATA / "REPEAT_pkl_files/onshorewind_route_dict.pkl")

# state_distribution_wind = pd.DataFrame()
# for year in selected_wind_cpas['year'].unique():
#     for p in Policies:
#         selected_wind_cpas_temp = selected_wind_cpas[(selected_wind_cpas['year']==year) & (selected_wind_cpas['policy']==p)]
#         results_dict, summary_dict = cpa_transmission_stats_by_state(list(selected_wind_cpas_temp['CPA_ID'].unique()), wind_ls_dict, wind_route_dict, all_wind_cpas)
#         A = pd.DataFrame.from_dict({(i,j): summary_dict[i][j] for i in summary_dict.keys() for j in summary_dict[i].keys()}, orient='index')
#         A['year'] = year
#         A['policy'] = p
#         IP()
#         state_distribution_wind = state_distribution_wind.append(A)
# state_distribution_wind['tech'] = 'wind'

# calculating offshore spur line cost distribution

selected_offshorewind_cpas = selected_cpas[(selected_cpas.tech.isin(['offshore_wind_fixed','offshore_wind_floating','offshore_wind_fixed_np','offshore_wind_floating_np']))]# & (~selected_cpas['CPA_ID'].isin([395071, 405135, 377240, 395763]))]
all_offshorewind_cpas = pd.read_csv("lcoe/offshorewind_lcoe_rio_emm_metro_county.csv")
offshorewind_ls_dict = pickle.load(open(DATA / "REPEAT_pkl_files/offshore_ls_dict_20230711.pkl", "rb"))
offshorewind_route_dict = pd.read_pickle(DATA / "REPEAT_pkl_files/offshore_route_dict_20230711.pkl")

state_distribution_offshorewind = pd.DataFrame()
for year in selected_offshorewind_cpas['year'].unique():
    for p in Policies:
        selected_offshorewind_cpas_temp = selected_offshorewind_cpas[(selected_offshorewind_cpas['year']==year) & (selected_offshorewind_cpas['policy']==p)]
        results_dict, summary_dict = cpa_transmission_stats_by_state(list(selected_offshorewind_cpas_temp['CPA_ID'].unique()), offshorewind_ls_dict, offshorewind_route_dict, all_offshorewind_cpas)
        A = pd.DataFrame.from_dict({(i,j): summary_dict[i][j] for i in summary_dict.keys() for j in summary_dict[i].keys()}, orient='index')
        A['year'] = year
        A['policy'] = p
        state_distribution_offshorewind = state_distribution_offshorewind.append(A)
state_distribution_offshorewind['tech'] = 'offshorewind'

# ##temp script
# state_distribution = pd.read_csv('/Users/npatankar/D/US_Policy/downscaling_Jan23/capacity_downscaling/result_data/spur_line/spur_line_state_distribution.csv')
# state_distribution = state_distribution.append(state_distribution_offshorewind)
# state_distribution.to_csv('result_data/spur_line_state_distribution.csv')
# ##

state_distribution = state_distribution_solar.append(state_distribution_wind)
state_distribution = state_distribution.append(state_distribution_offshorewind)
state_distribution.to_csv('result_data/spur_line_state_distribution.csv')

# finding incremental capex added by state
#IP()
state_distribution1 = pd.read_csv("result_data/spur_line_state_distribution_correct.csv")
state_distribution1 = state_distribution1.fillna(0)
state_distribution1.sort_values(['year'], inplace=True)
grouped_df = state_distribution1.groupby(['policy', 'tech', 'Unnamed: 0', 'Unnamed: 1'])
state_distribution1['gw-mile_Difference'] = grouped_df['gw-mile'].diff()
state_distribution1['capex_km_mw_Difference'] = grouped_df['capex_km_mw'].diff()
state_distribution1['dist_km_Difference'] = grouped_df['dist_km'].diff()
state_distribution1['capex_Difference'] = grouped_df['capex'].diff()
state_distribution1.to_csv('result_data/spur_line_state_distribution_delta_y.csv')

