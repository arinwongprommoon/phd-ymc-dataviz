#!/usr/bin/env python3

import numpy as np
import pandas as pd

data_directory_arinsbs = "/home/arin/git/time-series-pipeline/data/arin/"
experiment_prefix = "Omero20016_"

STRAINLOOKUP_SUFFIX = "strains.csv"
FLAVIN_SUFFIX = "flavin.csv"
MCHERRY_SUFFIX = "mCherry.csv"
BIRTHS_SUFFIX = "births.csv"
# LABELS_SUFFIX = "labels.csv" # TODO: potentially make handling labels similar to the others?

# mainly handles births
class SignalCSVHandler:
    def __init__(self, filepath_prefix):
        self.filepath_prefix = filepath_prefix

    def import_signal(self, filename_suffix):
        signal = pd.read_csv(self.filepath_prefix + filename_suffix)
        return signal


# to handle traces e.g. flavin, mCherry -- need to replace 0 with NaN
class TraceCSVHandler(SignalCSVHandler):
    def __init__(self, filepath_prefix):
        super().__init__(filepath_prefix)

    def import_trace(self, filename_suffix):
        signal = super().import_signal(filename_suffix)
        return signal.replace(0, np.nan)


# wrappers to simplify code
def import_births(filepath_prefix):
    handler = SignalCSVHandler(filepath_prefix)
    return handler.import_signal(BIRTHS_SUFFIX)


def import_flavin(filepath_prefix):
    handler = TraceCSVHandler(filepath_prefix)
    return handler.import_trace(FLAVIN_SUFFIX)


def import_mCherry(filepath_prefix):
    handler = TraceCSVHandler(filepath_prefix)
    return handler.import_trace(MCHERRY_SUFFIX)


# labels
def import_labels(filepath, list_cellIDs):
    targets_df = pd.read_csv(filepath, header=None, index_col=0)
    targets_df.index.names = ["cellID"]
    targets = targets_df.loc[list_cellIDs].to_numpy().flatten()
    return targets


def convert_df_to_aliby(signal, strainlookup_filepath_prefix):
    # Import look-up table for strains
    strainlookup_df = pd.read_csv(strainlookup_filepath_prefix + STRAINLOOKUP_SUFFIX)
    strainlookup_dict = dict(zip(strainlookup_df.position, strainlookup_df.strain))
    # Positions -> Strain (more informative)
    signal = signal.replace({"position": strainlookup_dict})
    signal.rename(columns={"position": "strain"}, inplace=True)
    signal = signal.drop(["distfromcentre"], axis=1)
    # Convert to multi-index dataframe
    signal_temp = signal.iloc[:, 2:]
    multiindex = pd.MultiIndex.from_frame(signal[["strain", "cellID"]])
    signal = pd.DataFrame(signal_temp.to_numpy(), index=multiindex)
    return signal


# filepath_prefix = data_directory + experiment_prefix
# signal_flavin = convert_df_to_aliby(import_flavin(filepath_prefix), filepath_prefix)
# signal_mCherry = convert_df_to_aliby(import_mCherry(filepath_prefix), filepath_prefix)
# signal_births = convert_df_to_aliby(import_births(filepath_prefix), filepath_prefix)


def df_to_biodare(df, filename):
    """Export any dataframe to CSV format preferred by BioDARE

    https://biodare2.ed.ac.uk/
    """
    df.T.to_csv(filename)
