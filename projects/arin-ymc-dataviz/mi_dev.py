#!/usr/bin/env python3

import numpy as np
import pandas as pd
from postprocessor.core.multisignal.mi import mi, miParameters
from postprocessor.core.multisignal.micatch22 import micatch22, micatch22Parameters
from postprocessor.core.processes.catch22 import catch22
from postprocessor.core.processes.standardscaler import standardscaler

# 1. Import Alejandro's data
filepath_prefix = "/home/arin/git/time-series-pipeline/data/agranados/infdata_rep"
MEDIUM_LIST = ["rich", "stress"]
replicate = 1

# convenience functions for a specific data structure
def convert_csv_agranados_to_aliby(replicate, medium):
    signal = pd.read_csv(filepath_prefix + str(replicate) + "_" + medium + ".csv")
    multiindex_array = [[medium] * len(signal), list(range(len(signal)))]
    multiindex = pd.MultiIndex.from_arrays(multiindex_array, names=("strain", "cellID"))
    signal = pd.DataFrame(signal.to_numpy(), multiindex)
    return signal


signal_rich = convert_csv_agranados_to_aliby(replicate, "rich")
signal_stress = convert_csv_agranados_to_aliby(replicate, "stress")

# 2. mi-catch22
results_mi = mi.as_function([signal_rich, signal_stress], overtime=False)
results_micatch22 = micatch22.as_function([signal_rich, signal_stress], overtime=False)
breakpoint()
print("hello")

# variance of catch22
def get_variances(df):
    features = catch22.as_function(df)
    scaled = standardscaler.as_function(features)
    return scaled.var().sort_values(ascending=False)
