# Code description

This is code originally written to work with inputs from EER's RIO model. It should be modified to work with files created by PowerGenome and GenX.

## downscale_results.py

This script has functions for loading the model results (`load_capacity_results()`), loading the LCOE file for each technology (`load_lcoe_file()`), selecting the least cost CPAs that would satisfy the capacity built by a model using their derated capacity (`select_cpas()`), and then randomly selecting clusters of CPAs based on their derated capacity (`random_select_cluster_cpas()`).

Results for each run/year are saved in CSV files. The full pipeline is accessed using the function `downscale_cluster_cpas()`.

## plot_results.py

This script has functions that take downscaled results and plot them using either matplotlib (`mlp_maps()`) or Altair (`alt_maps()`). It was originally run in interactive mode -- downscaled results were generated and then the maps were generated.
