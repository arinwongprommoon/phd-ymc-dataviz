#!/usr/bin/env python3

import pandas as pd


def slice_signal(signal, cell_index_list, interval_start, interval_end):
    """Slice signal by list of cell indices and interval start & end."""
    out_signal = signal.loc[cell_index_list]
    out_signal = out_signal.iloc[:, interval_start:interval_end]
    return out_signal


def slice_strain(signal_df, strain_list):
    """Subset a signal DataFrame based on a list of strains to include.

    Parameters
    ----------
    signal_df : pandas.DataFrame
        signal DataFrame
    strain_list : list
        list of strains (strings)

    Examples
    --------
    FIXME: Add docs.

    """
    return signal_df.loc[strain_list]


def slice_interval(signal_df, interval_start, interval_end):
    """Slice signal DataFrame to a starting and ending timepoint (NaNs removed).

    Parameters
    ----------
    signal_df : pandas.DataFrame
        signal DataFrame
    interval_start : int
        starting timepoint
    interval_end : int
        ending timepoint

    Examples
    --------
    FIXME: Add docs.

    """
    return signal_df.iloc[:, interval_start:interval_end].dropna()
