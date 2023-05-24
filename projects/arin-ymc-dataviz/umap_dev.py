#!/usr/bin/env python3

import matplotlib.pyplot as plt
import umapper
from exptio import (
    convert_df_to_aliby,
    data_directory_arinsbs,
    experiment_prefix,
    import_flavin,
    import_labels,
)
from postprocessor.core.processes.catch22 import catch22, catch22Parameters
from slicer import slice_interval, slice_strain

catch22processor = catch22(catch22Parameters.default())

# Load data, featurise.
data_directory = data_directory_arinsbs
filepath_prefix = data_directory + experiment_prefix
filename_targets = "categories_20016_detrend.csv"

signal_flavin = convert_df_to_aliby(import_flavin(filepath_prefix), filepath_prefix)
signal_flavin_processed = slice_strain(signal_flavin, ["by4741", "zwf1_Del"])
signal_flavin_processed = slice_interval(signal_flavin_processed, 25, 168)
features = catch22processor.run(signal_flavin_processed)

# Load targets (i.e. oscillatory/non-oscillatory ground truths)
targets = import_labels(
    "/home/arin/git/time-series-pipeline/" + filename_targets,
    signal_flavin_processed.index.get_level_values("cellID"),
)

# Drop non-oscillatory
features = features.iloc[targets == 1]

# Compute and draw UMAP
fig_umap, ax_umap = plt.subplots()
umapper.umap_plot(features, n_neighbors=20, min_dist=0.5, n_components=2, ax=ax_umap)
