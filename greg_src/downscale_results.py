from typing import List, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn import cluster, preprocessing
from typer import run


# %%
def load_lcoe_file(path: str) -> pd.DataFrame:
    """Load a LCOE file from the `path` location, filter out existing facilities, and
    return with a subset of the columns.

    Parameters
    ----------
    path : str
        Path of file

    Returns
    -------
    pd.DataFrame
        Data from file
    """
    keep_cols = [
        "CPA_ID",
        "incap",
        "cpa_mw",
        "lcoe",
        "total_interconnect_km",
        "d_shore_m",
        "latitude",
        "longitude",
        "m_popden",
        "exFacil",
    ]
    df = pd.read_csv(
        path,
    )
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df.loc[:, keep_cols]
    if "exFacil" in df.columns:
        df = df.query("exFacil==0")
    df = df.sort_values("lcoe", ignore_index=True)
    return df


def select_cpas(
    lcoe_df: pd.DataFrame,
    cpa_ids: List[int],
    mw_build: float,
    offshore_distance: float = None,
    run: str = None,
    region: str = None,
) -> pd.DataFrame:
    """Select the subset of all CPAs that satisfy the `mw_built` requirement.

    Parameters
    ----------
    lcoe_df : pd.DataFrame
        LCOE file with columns "lcoe" and "cpa_mw"
    cpa_ids : List[int]
        All CPAs that were included in the original resource cluster
    mw_build : float
        The capacity (MW) of a resource that was built by the model
    offshore_distance : float, optional
        The distance to shore (km) to try limiting offshore wind, by default None
    run : str, optional
        Name of the run, by default None
    region : str, optional
        Name of the region, by default None

    Returns
    -------
    pd.DataFrame
        _description_
    """
    df = lcoe_df.loc[lcoe_df["CPA_ID"].isin(cpa_ids), :].reset_index(drop=True)
    if offshore_distance is not None and "d_shore_m" in lcoe_df.columns:
        total_cap = df["cpa_mw"].sum()
        cap_with_buffer = df.query("d_shore_m >= @offshore_distance")["cpa_mw"].sum()
        if df.query("d_shore_m >= @offshore_distance")["cpa_mw"].sum() >= mw_build:
            df = df.query("d_shore_m >= @offshore_distance")
        else:
            df = sort_mask(mw_build, df, "d_shore_m")
            min_d_shore = df["d_shore_m"].min()
            print(
                f"Run {run}, region {region} has {cap_with_buffer} out of {total_cap}, model built {mw_build} ({cap_with_buffer / mw_build}. Min distance is {min_d_shore.round(1)})"
            )
            return df

    _df = sort_mask(mw_build, df)

    return _df


def sort_mask(
    mw_build: float, df: pd.DataFrame, sort_col: str = "lcoe"
) -> pd.DataFrame:
    """Sort the df by "sort_col" and return the rows that satisfy the mw_build value

    Parameters
    ----------
    mw_build : float
        Capacity (MW) of a resource that is built
    df : pd.DataFrame
        DataFrame with individual candidate sites in each row. Should have a column "cpa_mw"
        and a column matching `sort_col`.
    sort_col : str, optional
        Name of the column used to sort candidate sites, by default "lcoe"

    Returns
    -------
    pd.DataFrame
        Dataframe with only the least cost candidate sites that are needed to match the
        mw_build
    """
    _df = df.sort_values(sort_col)
    mask = np.ones(len(_df), dtype=bool)
    # temp = (_df.loc[mask, "incap"].cumsum() < mw_build).values
    temp = (_df.loc[mask, "cpa_mw"].cumsum() < mw_build).values
    temp[temp.argmin()] = True
    mask[mask] = temp
    return _df.loc[mask, :]


def tech_capitalize(s: str) -> str:
    return (
        s.replace("class", "Class")
        .replace("offshorewind", "OffShoreWind")
        .replace("landbasedwind", "LandbasedWind")
        .replace("utilitypv", "UtilityPV")
        .replace("moderate", "Moderate")
        .replace("conservative", "Conservative")
        .replace("advanced", "Advanced")
    )


# %%
# Load data and correct tech names
def load_capacity_results() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load the EER files "binned_resource_df.csv" and "capacity_y.csv". Format the capacity
    column "tech" to match the binned resource names.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Binned resource and capacity dataframes.
    """
    binned_resource = pd.read_csv(
        "binned_resource_df.csv", usecols=["cpa_id", "us_eer", "RIO_tech_name"]
    )
    binned_resource = binned_resource.dropna(subset=["RIO_tech_name"])
    capacity = pd.read_csv("capacity_y.csv")
    capacity["tech"] = capacity["tech"].str.replace(
        "onshore wind|", "onshore_wind|", regex=False
    )
    capacity.loc[
        capacity["tech||outputs_group_detailed"].str.contains("offshore wind - f"),
        "tech",
    ] = (
        capacity.loc[
            capacity["tech||outputs_group_detailed"].str.contains("offshore wind - f"),
            "tech||outputs_group_detailed",
        ]
        .str.split(" - ")
        .str[-1]
        + "_"
        + capacity.loc[
            capacity["tech||outputs_group_detailed"].str.contains("offshore wind - f"),
            "tech",
        ]
    ).str.replace(
        " ", "_"
    )
    capacity.loc[capacity["tech"] == "utility-scale solar pv|1", "tech"] = "solar|1"
    capacity = capacity.loc[
        capacity["tech"].isin(binned_resource["RIO_tech_name"].unique()),
        :,
    ]
    capacity["value"] *= 1000
    capacity["unit"] = "megawatt"

    return capacity, binned_resource


def random_select_cluster_cpas(df: pd.DataFrame, tech: str) -> pd.DataFrame:
    """Randomly select clusters of CPAs within each population density bin and return
    all selected CPAs.

    CPAs within each population density bin are clustered according to lat/lon, with the
    number of clusters roughly accounting for how many CPAs would be in a typical project
    (my assumption of what the values represent). Cluster factors of 1.1 for solar and 6
    for onshore wind are based on original method by Neha Patankar.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of candidate sites that are available for random selection.
    tech : str
        Name of the technology (solar or onshore_wind)

    Returns
    -------
    pd.DataFrame
        Randomly selected CPAs
    """
    if df.empty:
        return pd.DataFrame()
    # Tuple values are the capacity retained and upper limit of population density bin.
    # Solar values are multiplied by 0.2 to reflect automatic 80% derating of CPAs
    pop_frac = {
        "solar": [
            (1 * 0.2, 5),
            (0.5 * 0.2, 10),
            (0.25 * 0.2, 30),
            (0.125 * 0.2, 40),
            (0.1 * 0.2, 60),
        ],
        "onshore_wind": [(0.6, 1), (0.2, 3), (0.1, 5), (0.05, 10)],
    }

    # Number of clusters is based on some assumption about how many CPAs should be adjacent
    # in a typical project.
    cluster_factor = {"solar": 1.1, "onshore_wind": 6}
    cluster_cols = ["latitude", "longitude"]
    n_clusters = max(1, int(len(df) // cluster_factor[tech]))
    # Use StandardScaler to normalize lat/lon values and avoid bias
    clusters = cluster.KMeans(n_clusters=n_clusters).fit(
        preprocessing.StandardScaler().fit_transform(df[cluster_cols])
    )
    df["cluster"] = clusters.labels_
    if len(df) == 1:
        return df
    df_list = []
    popden_low = 0
    for frac, popden_high in pop_frac[tech]:
        selected_df = df.query("@popden_low < m_popden <= @popden_high")
        selected_clust = (
            pd.DataFrame(selected_df["cluster"].unique())
            .sample(frac=frac, replace=False)
            .rename(columns={0: "cluster"})
        )
        final_df = selected_df.loc[
            selected_df["cluster"].isin(selected_clust["cluster"])
        ]
        df_list.append(final_df)
        popden_low = popden_high

    return pd.concat(df_list, ignore_index=True)


# def random_select_tech_cpas(cluster_cpas: List[int], )


def select_final_cpas(
    tech: str, row: pd.Series, cpa_map: pd.DataFrame, lcoe_df: pd.DataFrame
) -> pd.DataFrame:
    energy_park = 0
    if "energy_park" in row["tech"]:
        energy_park = 1
        # r = "onshore_wind|energy_park"
    cap_built = row["value"]
    cluster_cpas = cpa_map.loc[
        (cpa_map["RIO_tech_name"] == row["tech"]) & (cpa_map["us_eer"] == row["zone"]),
        "cpa_id",
    ].to_list()
    if cluster_cpas:
        cpas_built = select_cpas(
            lcoe_df.copy(),
            cluster_cpas,
            cap_built,
            offshore_distance=None,
            run=run,
            region=row["zone"],
        )
        if tech != "offshore_wind":
            cpas_built = random_select_cluster_cpas(cpas_built, tech)
        cpas_built["technology"] = tech
        cpas_built["energy_park"] = energy_park
    return cpas_built


# %%
def downscale_cluster_cpas(
    save_path: Path,
    year=2050,
    solar_lcoe_path: str = None,
    onshorewind_lcoe_path: str = None,
    offshorewind_lcoe_path: str = None,
    offshore_distance=None,
):
    cap_results, cpa_map = load_capacity_results()
    cap_results = cap_results.loc[cap_results["year"] == year, :]
    cap_results = cap_results.loc[
        ~cap_results["zone"].isin(["alaska_emm", "hawaii_emm"])
    ]

    for run, cap_result in cap_results.groupby("run name"):
        print(run)

        df_list = []
        for r, lcoe_path in zip(
            ["solar", "onshore_wind", "offshore_wind"],
            [solar_lcoe_path, onshorewind_lcoe_path, offshorewind_lcoe_path],
        ):
            tech_df_list = []
            if lcoe_path is None:
                continue
            print(r)
            lcoe_df = load_lcoe_file(lcoe_path)
            cap_built = cap_result.loc[
                cap_result["tech"].str.contains(r, case=False), "value"
            ]
            # tech_list = Parallel(n_jobs=-1, verbose=10)(
            #     delayed(select_final_cpas)(r, row, cpa_map, lcoe_df)
            #     for idx, row in cap_result.loc[
            #         cap_result["tech"].str.contains(r, case=False), :
            #     ].iterrows()
            # )
            # df_list.extend(tech_list)
            for idx, row in cap_result.loc[
                cap_result["tech"].str.contains(r, case=False), :
            ].iterrows():
                energy_park = 0
                if "energy_park" in row["tech"]:
                    energy_park = 1
                    # r = "onshore_wind|energy_park"
                cap_built = row["value"]
                cluster_cpas = cpa_map.loc[
                    (cpa_map["RIO_tech_name"] == row["tech"])
                    & (cpa_map["us_eer"] == row["zone"]),
                    "cpa_id",
                ].to_list()

                if cluster_cpas:
                    cpas_built = select_cpas(
                        lcoe_df.copy(),
                        cluster_cpas,
                        cap_built,
                        offshore_distance=offshore_distance,
                        run=run,
                        region=row["zone"],
                    )
                    cpas_built["zone"] = row["zone"]
                    cpas_built["technology"] = r
                    cpas_built["energy_park"] = energy_park
                    tech_df_list.append(cpas_built)

            tech_cpas_built = pd.concat(tech_df_list)
            if r != "offshore_wind":
                built_list = []
                for _, _df in tech_cpas_built.groupby(["zone", "energy_park"]):
                    _df = random_select_cluster_cpas(_df, r)
                    built_list.append(_df)
                # cpas_built = random_select_cpas(cpas_built, r)
                selected_cpas_built = pd.concat(built_list)

                df_list.append(selected_cpas_built)
            else:
                df_list.append(tech_cpas_built)
        # out_folder = Path(__file__).parent / "v4"
        save_path.mkdir(exist_ok=True)
        pd.concat(df_list, ignore_index=True).to_csv(
            save_path / f"{run}_{year}_downscaling_cluster.csv", index=False
        )


# %%
if __name__ == "__main__":
    downscale_cluster_cpas(
        solar_lcoe_path="",
        onshorewind_lcoe_path="",
        offshorewind_lcoe_path="",
        offshore_distance=15,
    )

# # %%
# cluster.KMeans(n_clusters=num_clusters[region][tech], random_state=6).fit(
#     preprocessing.StandardScaler().fit_transform(grouped[cluster_cols])
# )

# cluster_cols = ["latitude", "longitude"]
# n_clusters = int(len(sel_solar_) // 1.1)
# clusters = cluster.KMeans(n_clusters=n_clusters).fit(
#     preprocessing.StandardScaler().fit_transform(sel_solar_[cluster_cols])
# )
# sel_solar_["cluster"] = clusters.labels_

# df_list = []
# popden_low = 0
# for frac, popden_high in [(1, 5), (0.5, 10), (0.25, 30), (0.125, 40), (0.1, 60)]:
#     print(popden_low, popden_high)
#     selected_solar = sel_solar_.query("@popden_low < m_popden <= @popden_high")
#     selected_clust = (
#         pd.DataFrame(selected_solar["cluster"].unique())
#         .sample(frac=frac, replace=False)
#         .rename(columns={0: "cluster"})
#     )
#     final_solar = selected_solar.loc[
#         selected_solar["cluster"].isin(selected_clust["cluster"])
#     ]
#     start_cap = selected_solar["incap"].sum()
#     final_cap = final_solar["incap"].sum()
#     selected_frac = (final_cap / start_cap).round(3)
#     print(f"{selected_frac} selected, {frac} fraction expected")
#     df_list.append(final_solar)
#     popden_low = popden_high

# solar_pop_frac = [(1, 5), (0.5, 10), (0.25, 30), (0.125, 40), (0.1, 60)]
# wind_pop_frac = [(0.6, 1), (0.2, 3), (0.1, 5), (0.05, 10)]
# solar_cluster_factor = 1.1
# wind_cluster_factor = 6


# sel_S_0 = sel_solar_.query("m_popden <= 5")
# sel_S_0_c = (
#     pd.DataFrame(sel_S_0["cluster"].unique())
#     .sample(frac=0.2, replace=False)
#     .rename(columns={0: "cluster"})
# ).reset_index(drop=True)
# # sel_S_0_c["unique_id"] = sel_S_0_c.index + 1
# # sel_S_0_f_ = pd.merge(
# #     sel_S_0,
# #     sel_S_0_c,
# #     how="right"
# # )
# sel_S_0_f = sel_S_0.loc[sel_S_0["cluster"].isin(sel_S_0_c["cluster"]), :]

# # %%
