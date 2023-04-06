#!/usr/bin/env python3

import igraph as ig
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
from postprocessor.core.processes.standardscaler import standardscaler
from sklearn.metrics.pairwise import euclidean_distances
from slicer import slice_interval, slice_strain
from utils import generate_palette_map, graph_prune

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

# Scale catch22 data matrix
features_scaled = standardscaler.as_function(features.T).T

# Compute distance matrix
distance_matrix = euclidean_distances(features_scaled)

# Prune distance matrix
distance_matrix_pruned = graph_prune(distance_matrix, 10)

# Creates graph
graph = ig.Graph.Weighted_Adjacency(distance_matrix_pruned.tolist(), mode="undirected")
graph.vs["strain"] = features.index.get_level_values("strain")

# Plot
palette_map = generate_palette_map(features)
# NOTE: variable edge widths do not work with the matplotlib backend;
# must use `cairo`.
visual_style = {
    "vertex_color": [palette_map[strain] for strain in graph.vs["strain"]],
    "layout": graph.layout("drl"),
    # "edge_width": [round(weight) for weight in graph.es["weight"]],
    "edge_width": 0.2,
}
fig, ax = plt.subplots()
ig.plot(graph, target=ax, **visual_style)
