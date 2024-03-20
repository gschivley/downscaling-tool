"Plot downscaled results"

# %load_ext autoreload
# %autoreload 2

import os

os.environ["USE_PYGEOS"] = "0"
import pandas as pd
import geopandas as gpd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import altair as alt
from joblib import Parallel, delayed

from downscale_results import downscale_cluster_cpas

sns.set_style(
    "white",
    {
        "axes.spines.left": False,
        "axes.spines.bottom": False,
        "axes.spines.right": False,
        "axes.spines.top": False,
    },
)
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300

results_folder = Path.cwd() / "v3"

geo_files = {
    "solar": "",
    "onshore_wind": "",
    "offshore_wind": "",
}

tech_color = {
    "solar_0": "#FF0000",  # "#E69800", "red",
    "solar_existing": "#FFA500",
    "onshore_wind_0": "#27408B",  # "#005CE6", "#27408B" (royalblue4)
    "onshore_wind_1": "#38A800",
    "onshore_wind_existing": "#ADD8E6",
    "offshore_wind_0": "#EE82EE",  # "#DF73FF",
}


state_gdf = gpd.read_file("cb_2022_us_state_5m.zip")
state_gdf = state_gdf.query("~STUSPS.isin(['HI', 'AK', 'PR', 'GU', 'VI', 'MP', 'AS'])")
state_lat_lon = state_gdf.to_crs("EPSG:4326")
state_gdf = state_gdf.to_crs("EPSG:5070")

msa_gdf = gpd.read_file("")
tx_gdf = gpd.read_file("")


def mpl_maps(
    results_folder: Path,
    save_folder: Path,
    n_jobs=1,
    format="png",
    annotate=False,
    **kwargs,
):
    tech_gdf_albers_dict = {}
    base_map = state_gdf.plot(color="white", ec="#6E6E6E", lw=0.3, figsize=(10, 7.5))
    base_map = msa_gdf.boundary.plot(ax=base_map, color="#CCCCCC", lw=0.2)  # "#CCCCCC"
    base_map = tx_gdf.plot(ax=base_map, color="#686868", lw=0.1)
    tech_gdf_albers_dict = {
        tech: gpd.read_parquet(tech_path).to_crs("EPSG:5070")
        for tech, tech_path in geo_files.items()
    }
    Parallel(n_jobs=n_jobs)(
        delayed(single_mpl_map)(
            save_folder, tech_gdf_albers_dict, base_map, f, format, annotate, **kwargs
        )
        for f in list(results_folder.rglob("*.csv"))
    )


def single_mpl_map(
    save_folder, tech_gdf_albers_dict, base_map, f, format, annotate, **kwargs
):
    scenario = f.name.split("_")[0]
    year = f.name.split("_")[1]
    df = pd.read_csv(f)
    final_map = None
    for (tech, ep), _df in df.groupby(["technology", "energy_park"]):
        if tech_gdf_albers_dict.get(tech) is None:
            tech_gdf_albers_dict[tech] = gpd.read_parquet(geo_files[tech]).to_crs(
                "EPSG:5070"
            )

            # if tech == "onshore_wind":
        final_map = (
            tech_gdf_albers_dict[tech]
            .loc[tech_gdf_albers_dict[tech]["CPA_ID"].isin(_df["CPA_ID"]), :]
            .plot(
                ax=(final_map or base_map),
                color=tech_color[f"{tech}_{ep}"],
                ec=tech_color[f"{tech}_{ep}"],
                lw=0.05,
            )
        )
        if "exFacil" in tech_gdf_albers_dict[tech].columns and ep == 0:
            final_map = (
                tech_gdf_albers_dict[tech]
                .loc[tech_gdf_albers_dict[tech]["exFacil"] == 1, :]
                .plot(
                    ax=final_map,
                    color=tech_color[f"{tech}_existing"],
                    ec=tech_color[f"{tech}_existing"],
                    lw=0.05,
                    # alpha=0.5,
                )
            )
    if annotate is True:
        base_map = add_mpl_annotations(base_map)
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    save_folder.mkdir(exist_ok=True)
    plt.savefig(
        save_folder / f"{scenario}_{year}_cpa_map.{format}", **kwargs
    )  # , dpi=300, bbox_inches="tight"
    # )


def add_mpl_annotations(base_map: plt.Axes):

    text = [
        "Existing Utility-scale Solar",
        "New Utility-Scale Solar",
        "Existing Onshore Wind",
        "New Onshore Wind",
        "New Onshore Wind Energy Park",
        "New Offshore Wind",
    ]
    x = [2.5e6, 2.5e6, 2.5e6, 2.5e6, 2.5e6, 2.5e6]
    y = [2.6e6, 2.4e6, 2.2e6, 2.0e6, 1.8e6, 1.6e6]
    color = ["#FFA500", "#FF0000", "#ADD8E6", "#27408B", "#38A800", "#EE82EE"]

    for _t, _x, _y, _c in zip(text, x, y, color):
        base_map.text(
            _x + 1e5,
            _y,
            _t,
            horizontalalignment="left",
            verticalalignment="center_baseline",
            fontsize=9,
        )
        base_map.scatter(_x, _y, c=_c, s=30)

    return base_map


####################
# Altair charts
####################

_msa_gdf = msa_gdf.to_crs("EPSG:4326")
_tx_gdf = tx_gdf.to_crs("EPSG:4326")

map_charts = []
base_map = (
    alt.Chart(state_lat_lon)
    .mark_geoshape(fill="white", stroke="#6E6E6E", strokeWidth=0.3)
    .properties(
        width=500,
        height=300,
    )
    .project("albersUsa")
)
map_charts.append(base_map)
msa_map = (
    alt.Chart(_msa_gdf)
    .mark_geoshape(fill=None, stroke="#CCCCCC", strokeWidth=0.2)
    .properties(
        width=500,
        height=300,
    )
    .project("albersUsa")
)
map_charts.append(msa_map)
tx_map = (
    alt.Chart(_tx_gdf)
    .mark_geoshape(fill=None, stroke="#686868", strokeWidth=0.1)
    .properties(
        width=500,
        height=300,
    )
    .project("albersUsa")
)
map_charts.append(tx_map)


def alt_maps(
    results_folder: Path,
    save_folder: Path,
    map_charts,
    n_jobs=1,
    format="png",
    **kwargs,
):
    tech_gdf_ll_dict = {
        tech: gpd.read_parquet(tech_path).to_crs("EPSG:4326")
        for tech, tech_path in geo_files.items()
    }

    for f in list(results_folder.glob("*.csv")):
        single_alt_map(save_folder, map_charts, tech_gdf_ll_dict, f, format, **kwargs)

    Parallel(n_jobs=n_jobs)(
        delayed(single_mpl_map)(save_folder, tech_gdf_ll_dict, f, format, **kwargs)
        for f in list(results_folder.rglob("*.csv"))
    )


def single_alt_map(save_folder, map_charts, tech_gdf_ll_dict, f, format, **kwargs):
    df = pd.read_csv(f)
    scenario = f.name.split("_")[0]
    year = f.name.split("_")[1]
    tech_charts = []
    for (tech, ep), _df in df.groupby(["technology", "energy_park"]):
        if tech_gdf_ll_dict.get(tech) is None:
            tech_gdf_ll_dict[tech] = gpd.read_parquet(geo_files[tech]).to_crs(
                "EPSG:4326"
            )

        data = tech_gdf_ll_dict[tech].loc[
            tech_gdf_ll_dict[tech]["CPA_ID"].isin(_df["CPA_ID"]), :
        ]
        tech_map = (
            alt.Chart(data)
            .mark_geoshape(
                fill=tech_color[f"{tech}_{ep}"],
                stroke=tech_color[f"{tech}_{ep}"],
                strokeWidth=0.05,
            )
            .properties(
                width=500,
                height=300,
            )
            .project("albersUsa")
        )
        tech_charts.append(tech_map)

        if "exFacil" in tech_gdf_ll_dict[tech].columns and ep == 0:
            ex_map = (
                alt.Chart(data)
                .mark_geoshape(
                    fill=tech_color[f"{tech}_existing"],
                    stroke=tech_color[f"{tech}_existing"],
                    strokeWidth=0.05,
                )
                .properties(
                    width=500,
                    height=300,
                )
                .project("albersUsa")
            )
            tech_charts.append(ex_map)

    chart = alt.layer(*map_charts, *tech_charts)
    save_folder.mkdir(exist_ok=True)
    chart.save(
        save_folder / f"{scenario}_{year}_cpa_map.{format}", **kwargs
    )  # ppi=300)


##############
# All downscaling/clustering results
##############

Parallel(n_jobs=-1)(
    delayed(downscale_cluster_cpas)(
        year=year,
        save_path=Path.cwd() / "v4",
        solar_lcoe_path="",
        onshorewind_lcoe_path="",
        offshorewind_lcoe_path="",
        offshore_distance=15,
    )
    for year in list(range(2025, 2051, 5))
)

##############
# All figs using mpl
##############

mpl_maps(
    results_folder=Path.cwd() / "v4",
    save_folder=Path.cwd() / "v4" / "maps_existing_png",
    n_jobs=-1,
    format="png",
    dpi=300,
    bbox_inches="tight",
)
