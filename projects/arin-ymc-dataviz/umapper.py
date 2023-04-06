#!/usr/bin/env python3

"""UMAP-related wrappers"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import umap
from postprocessor.core.processes.standardscaler import standardscaler
from postprocessor.core.processes.umapembedding import (
    umapembedding,
    umapembeddingParameters,
)
from utils import generate_palette_map


class _UMAPEmbeddingPlotter:
    """Draw scatterplot of UMAP embedding"""

    def __init__(self, embedding, labels, palette_map):
        self.embedding = embedding
        self.labels = labels
        self.palette_map = palette_map

    def plot(self, ax):
        """Draw scatterplot on the provided Axes."""
        sns.scatterplot(
            x=self.embedding[:, 0],
            y=self.embedding[:, 1],
            hue=self.labels,
            palette=self.palette_map,
            ax=ax,
        )


def umap_embedding_plot(embedding, labels, palette_map, ax=None):
    """Plot UMAP embedding"""
    plotter = _UMAPEmbeddingPlotter(embedding, labels, palette_map)
    if ax is None:
        ax = plt.gca()
    plotter.plot(ax)
    return ax


# Wrapper for them all, with all the parameters, default for my most common
# use case
def umap_plot(
    data, n_neighbors, min_dist, n_components, scale=True, label_index="strain", ax=None
):
    """Compute embedding and draw UMAP"""
    params = umapembeddingParameters(
        n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components
    )

    runner = umapembedding(params)
    if scale:
        data_scaled = standardscaler.as_function(data.T).T
    else:
        data_scaled = data

    plotter = _UMAPEmbeddingPlotter(
        embedding=runner.run(data_scaled),
        labels=data.index.get_level_values(label_index),
        palette_map=generate_palette_map(data),
    )
    if ax is None:
        ax = plt.gca()
    plotter.plot(ax)
    return ax
