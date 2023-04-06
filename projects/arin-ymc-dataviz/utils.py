#!/usr/bin/env python3

import numpy as np
import pandas as pd
from matplotlib import cm, colors
from postprocessor.routines.median_plot import median_plot


def rearrange_by_sorted_list(df, my_list):
    """Rearrange rows of a DataFrame by a list"""
    return df.reindex(df.index[np.argsort(my_list)].to_list())


def pretty_format(string):
    """Replaces '_Del' with Delta character and formats CEN.PK"""
    # Strain names have '_Del' due to charset limitations of microscopy software
    # Similar problem with CEN.PK
    tmp = string.replace("CEN_PK_Mat_A_Koetter", "CEN.PK Mat A")
    tmp = tmp.replace("CEN_PK_Mat_A", "CEN.PK Mat A")
    tmp = tmp.replace("_Del_", "$\Delta$ ")
    tmp = tmp.replace("_Del", "$\Delta$")
    return tmp


def simple_median_plot(
    trace_df,
    ylabel,
    median_color="b",
    error_color="lightblue",
    xlabel="Time point",
    ax=None,
):
    """Wrapper for median plot to strip away hard-coding/irrelevant stuff"""
    return median_plot(
        trace_df=trace_df,
        trace_name="signal",
        label="signal",
        median_color=median_color,
        error_color=error_color,
        xlabel=xlabel,
        ylabel=ylabel,
        ax=ax,
    )


def graph_prune(distanceMatrix, neighbours):
    """
    Prunes a complete graph (input distance matrix), keeping at least a
    specified number of neighbours for each node.

    Parameters:
    -----------
    distanceMatrix = 2D numpy array
    neighbours = integer

    Return: Dij_pruned, a 2D numpy array, represents distance matrix of pruned
            graph
    """
    Dij_temp = distanceMatrix
    Adj = np.zeros(distanceMatrix.shape)
    for ii in range(distanceMatrix.shape[0]):
        idx = np.argsort(Dij_temp[ii, :])
        Adj[ii, idx[1]] = 1
        Adj[idx[1], ii] = 1
        for jj in range(neighbours):
            Adj[ii, idx[jj]] = 1
            Adj[idx[jj], ii] = 1
    Dij_pruned = Dij_temp * Adj
    return Dij_pruned


def generate_palette_map(df):
    """Create a palette map based on the strains in a dataframe"""
    strain_list = np.unique(df.index.get_level_values("strain"))
    palette_cm = cm.get_cmap("Set1", len(strain_list) + 1)
    palette_rgb = [
        colors.rgb2hex(palette_cm(index / len(strain_list))[:3])
        for index, _ in enumerate(strain_list)
    ]
    palette_map = dict(zip(strain_list, palette_rgb))
    return palette_map


# https://stackoverflow.com/a/2566508
def find_nearest(array, value):
    """find index of nearest value in numpy array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
