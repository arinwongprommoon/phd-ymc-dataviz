#!/usr/bin/env python3

# TODO if have time: separate import and plotting

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from postprocessor.core.multisignal.align import alignParameters
from postprocessor.core.multisignal.crosscorr import crosscorr, crosscorrParameters
from postprocessor.core.processes.butter import butter, butterParameters
from postprocessor.core.processes.fft import fft, fftParameters
from postprocessor.core.processes.findpeaks import findpeaks, findpeaksParameters
from postprocessor.core.processes.standardscaler import (
    standardscaler,
    standardscalerParameters,
)
from postprocessor.grouper import NameGrouper
from postprocessor.routines.heatmap import heatmap
from postprocessor.routines.histogram import histogram
from postprocessor.routines.mean_plot import mean_plot
from postprocessor.routines.median_plot import median_plot
from postprocessor.routines.single_birth_plot import single_birth_plot

from exptio import import_labels
from process_pipeline import apply_all_processes
from signalcollection import RawSignalCollection, SlicedSignalCollection
from utils import pretty_format

### IMPORT AND SELECT DATA

# TODO:
# - Write dictionary with keys = OmeroIDs and values = strains so that script
#   can loop through.  Currently just looks at one strain in one experiment.

# HDF5 -- Import files, group by strain, choose signals

if True:
    data_directory = Path(
        "/home/jupyter-arin/data/1253_2023_03_23_flavin_by4742swain_by4742morgan_tsa1tsa2morgan_lysmedia_01_00/"
    )
    h5path_to_abbrev_dict = {
        "extraction/Flavin_bgsub/max/mean": ("flavin", "continuous"),
        #"extraction/mCherry2_bgsub/max/mean": ("mCherry", "continuous"),
        "extraction/general/None/area": ("area", "continuous"),
        "extraction/general/None/volume": ("volume", "continuous"),
        "postprocessing/buddings/extraction_general_None_volume": ("births", "binary"),
    }
    # Define oscillation classification, if any
    #classification_filename = "20016_by4741_labels.csv"

    grouper = NameGrouper(data_directory)
    # bodge: assumes that directory name begins with 5-digit OmeroID,
    # ideally should pull from metadata
    experimentID = grouper.name[0:4]
    raw = RawSignalCollection(grouper, h5path_to_abbrev_dict)


# Choose strain & interval
if True:
    # Choose strain and interval
    strain_name = "by4742swain"
    wildtype_name = "by4742swain"
    interval_start = 0+5
    interval_end = 240
    critical_freqs = 1 / 350

    remove_nosc = False
    if remove_nosc:
        classification_labels = pd.read_csv(classification_filename)
    else:
        classification_labels = None

    # Define the SignalCollection object that stores the signals from the strain
    # of interest, sliced
    strain = SlicedSignalCollection(
        raw, strain_name, interval_start, interval_end, classification_labels
    )
    wildtype = SlicedSignalCollection(
        raw, wildtype_name, interval_start, interval_end, classification_labels
    )

# For core post-processes
# Structure: appending string, process object, parameters object, operates on, gives out
# (Maybe better as an objects with inheritance as I keep repeating stuff?)
process_dict = {
    "/butter": {
        "runner": butter(
            butterParameters.from_dict(
                {
                    "order": 2,
                    "critical_freqs": critical_freqs,
                    "filter_type": "highpass",
                    "sampling_freq": 1 / 5,
                }
            )
        ),
        "signame_endswith": "",
        "input_sigtype": "continuous",
        "output_sigtype": "continuous",
    },
    # for htb2 levels 'stuck' in starvation
    # "/scale": {
    #    "runner": standardscaler(standardscalerParameters.default()),
    #    "signame_endswith": "",
    #    "input_sigtype": "continuous",
    #    "output_sigtype": "continuous",
    # },
    "_scale": {
        "runner": standardscaler(standardscalerParameters.default()),
        "signame_endswith": "/butter",
        "input_sigtype": "continuous",
        "output_sigtype": "continuous",
    },
    "_peaks": {
        "runner": findpeaks(findpeaksParameters.default()),
        "signame_endswith": "/butter",
        "input_sigtype": "continuous",
        "output_sigtype": "binary",
    },
    "_acf": {
        "runner": crosscorr(
            crosscorrParameters.from_dict(
                {
                    "normalised": True,
                    "only_pos": True,
                }
            )
        ),
        "signame_endswith": "/butter",
        "input_sigtype": "continuous",
        "output_sigtype": "continuous",
    },
    "_fft": {
        "runner": fft(fftParameters.default()),
        "signame_endswith": "/butter",
        "input_sigtype": "continuous",
        "output_sigtype": "continuous",
    },
}

# For multisignal post-processes
# Structure: appending string, process object, parameters object, operates on, gives out
multisignal_process_dict = {
    "/xcf": {
        "runner": crosscorr(crosscorrParameters.default()),
        "input0": "flavin/butter",
        "input1": "mCherry/butter",
        "output_sigtype": "continuous",
    }
}

# For align post-processes
# align is a special case because:
# (a) it gives out TWO outputs, either of which could be useful (but not always)
# (b) I re-use it with different parameters
# Structure: ID, process object, parameters object, operates on, gives out
align_process_dict = {
    "0": {
        "parameters": alignParameters.default(),
        "input0": "flavin/butter",
        "input1": "flavin/butter_peaks",
        "output0": True,
        "output0_sigtype": "continuous",
        "output1": False,
        "output1_sigtype": "binary",
        "infix": "_align_BY_",
    },
    "1": {
        "parameters": alignParameters.from_dict(
            {
                "slice_before_first_event": False,
                "events_at_least": 2,
            }
        ),
        "input0": "flavin/butter",
        "input1": "births",
        "output0": True,
        "output0_sigtype": "continuous",
        "output1": True,
        "output1_sigtype": "binary",
        "infix": "_align_BY_",
    },
    "2": {
        "parameters": alignParameters.from_dict(
            {
                "slice_before_first_event": True,
                "events_at_least": 2,
            }
        ),
        "input0": "flavin/butter",
        "input1": "births",
        "output0": True,
        "output0_sigtype": "continuous",
        "output1": True,
        "output1_sigtype": "binary",
        "infix": "_alignTrunc_BY_",
    },
}

if True:
    # Apply all processes to all signal collections
    signalcollection_list = [strain, wildtype]
    snr_cutoff = 0.01766784

    for signalcollection in signalcollection_list:
        apply_all_processes(
            signalcollection, process_dict, multisignal_process_dict, align_process_dict
        )
        # TODO: don't hard-code this, make it possible to do mCherry too
        # Peak intervals
        signalcollection.lists[
            "flavin/butter_peaks_intervals"
        ] = signalcollection.mask_to_intervals(
            signalcollection.signals["flavin/butter_peaks"]
        )
        # Birth intervals
        signalcollection.lists["births_intervals"] = signalcollection.mask_to_intervals(
            signalcollection.signals["births"]
        )
        # Principal frequencies from FFT
        signalcollection.lists[
            "flavin/butter_freqs"
        ] = signalcollection.get_principal_freqs(
            signalcollection.signals["flavin/butter_fft"][0],
            signalcollection.signals["flavin/butter_fft"][1],
        )
        if "mCherry" in signalcollection.signals.keys():
            signalcollection.lists[
                "mCherry/butter_freqs"
            ] = signalcollection.get_principal_freqs(
                signalcollection.signals["mCherry/butter_fft"][0],
                signalcollection.signals["mCherry/butter_fft"][1],
            )
        # Signal to noise ratio
        signalcollection.lists["flavin/butter_snr"] = signalcollection.get_snr(
            signalcollection.signals["flavin/butter_fft"][0],
            signalcollection.signals["flavin/butter_fft"][1],
            cutoff_freq=snr_cutoff,
        )


### PLOTTING

# Display
wild_type_on = True

# Convenience functions
# TODO: move them somewhere
def set_origin_axes(ax):
    # draw grid
    # ax.grid(True, which="both")

    # set the x-spine (see below for more info on `set_position`)
    # ax.spines["left"].set_position("zero")

    # turn off the right spine/ticks
    ax.spines["right"].set_color("none")
    ax.yaxis.tick_left()

    # set the y-spine
    # ax.spines["bottom"].set_position("zero")

    # turn off the top spine/ticks
    ax.spines["top"].set_color("none")
    ax.xaxis.tick_bottom()

    # draw axes
    ax.axhline(y=0, color="k", linewidth=0.5)
    ax.axvline(x=0, color="k", linewidth=0.5)


# TIME SERIES -- INDIVIDUAL

# Plot one time series
if False:
    cell_index = ("htb2mCherry_013", 39, 1)
    tick_spacing = 60

    fig_ts, ax_ts = plt.subplots(figsize=(20, 5))
    timepoints = strain.signals["births"].columns
    flavin_ts = strain.signals["flavin/butter_scale"].loc[cell_index, :]
    birth_mask = strain.signals["births"].loc[cell_index, :]
    plot_title = f"Sample time series"
    ylabel = "Normalised fluorescence (AU)"
    single_birth_plot(
        trace_timepoints=timepoints,
        trace_values=flavin_ts,
        trace_name="flavin",
        birth_mask=birth_mask,
        trace_color="#3714b0",
        plot_title=plot_title,
        ylabel=ylabel,
        ax=ax_ts,
    )
    if "mCherry" in strain.signals.keys():
        mCherry_ts = strain.signals["mCherry/butter_scale"].loc[cell_index, :]
        single_birth_plot(
            trace_timepoints=timepoints,
            trace_values=mCherry_ts,
            trace_name="HTB2::mCherry",
            birth_mask=birth_mask,
            trace_color="#cb0077",
            trace_linestyle="dotted",
            plot_title=plot_title,
            ylabel=ylabel,
            ax=ax_ts,
        )
    ax_ts.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

# Multicolour
if False:
    cellindex = 28
    tick_spacing = 60

    strain.signals["flavin/butter_scaled"] = standardscaler.as_function(
        strain.signals["flavin/butter"]
    )
    trace_time = strain.signals["flavin/butter_scaled"].columns
    y_plot = strain.signals["flavin/butter_scaled"].loc[cellindex]
    birth_mask_bool = strain.signals["births"].loc[cellindex].astype(bool)

    fig_ts_multicolour, ax_ts_multicolour = plt.subplots()
    fig_ts_multicolour.set_size_inches(10, 4)
    # Create line segments
    points = np.array([np.array(list(trace_time)), y_plot.to_numpy()]).T.reshape(
        -1, 1, 2
    )
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Map data points to colours
    norm = plt.Normalize(y_plot.min(), y_plot.max())
    lc = LineCollection(segments, cmap="RdBu", norm=norm)
    # Set values
    lc.set_array(y_plot)
    lc.set_linewidth(2)
    ax_ts_multicolour.add_collection(lc)
    ax_ts_multicolour.autoscale()
    # Add births
    for occurence, birth_time in enumerate(trace_time[birth_mask_bool]):
        if occurence == 0:
            label = "birth event"
        else:
            label = None
        ax_ts_multicolour.axvline(
            birth_time,
            color="k",
            linestyle="--",
            label=label,
        )
    # fig_ts.colorbar(line, ax=ax_ts)

# TIME SERIES -- POPULATION

# Plot mean time series (from raw, without filtering)
if True:
    fig_mean_ts_raw, ax_mean_ts_raw = plt.subplots()
    mean_plot(
        raw.signals["flavin"].loc[strain_name],
        unit_scaling=1 / 60,
        label=strain_name,
        mean_color="#3714b0",
        error_color="#725ac3",
        xlabel="Time (hours)",
        ylabel="Flavin fluorescence (AU)",
        plot_title="Mean time series (raw data)",
        ax=ax_mean_ts_raw,
    )

# Plot mean time series (after processing and filtering)
if True:
    fig_mean_ts_detrend, ax_mean_ts_detrend = plt.subplots()
    mean_plot(
        strain.signals["flavin/butter"],
        unit_scaling=1 / 60,
        label=strain_name,
        mean_color="#3714b0",
        error_color="#725ac3",
        xlabel="Time (hours)",
        plot_title="Mean time series (after processing and filtering)",
        ax=ax_mean_ts_detrend,
    )

# Plot median time series (after processing, filtering, and scaling)
# -- no alignment
if True:
    fig_median_ts_detrend, ax_median_ts_detrend = plt.subplots()
    median_plot(
        standardscaler.as_function(strain.signals["flavin/butter"]),
        unit_scaling=1 / 60,
        label=strain_name,
        median_color="#3714b0",
        error_color="#725ac3",
        xlabel="Time (hours)",
        plot_title="Median time series (after processing, filtering, scaling)",
        ax=ax_median_ts_detrend,
    )

# Plot heatmap
if True:
    fig_heatmap, ax_heatmap = plt.subplots()
    heatmap(
        strain.signals["flavin/butter"].dropna(how="all"),
        trace_name="flavin",
        buddings_df=strain.signals["births"].loc[
            strain.signals["flavin/butter"]
            .dropna()
            .index.intersection(strain.signals["births"].index)
        ],
        plot_title="Flavin fluorescence",
        unit_scaling=1 / 12,
        xtick_step=2,
        xlabel="Time (hours)",
        ax=ax_heatmap,
    )

# Plot MEDIAN flavin time series
if True:
    fig_median_flavin, ax_median_flavin = plt.subplots()
    num_cells_med = len(strain.signals["flavin/butter_align_BY_flavin_butter_peaks"])
    plot_title = f"Single-cell flavin fluorescence, aligned to first flavin peak\n(n = {num_cells_med})"
    if wild_type_on:
        median_plot(
            standardscaler.as_function(
                wildtype.signals["flavin/butter_align_BY_flavin_butter_peaks"]
            ),
            unit_scaling=1 / 60,
            label="wild type",
            median_color="grey",
            error_color="lightgrey",
            xlabel="Time (hours)",
            plot_title=plot_title,
            ax=ax_median_flavin,
        )
    median_plot(
        standardscaler.as_function(
            strain.signals["flavin/butter_align_BY_flavin_butter_peaks"]
        ),
        unit_scaling=1 / 60,
        label=pretty_format(strain_name),
        median_color="#3714b0",
        error_color="#725ac3",
        xlabel="Time (hours)",
        plot_title=plot_title,
        ax=ax_median_flavin,
    )

# Align to first birth
if True:
    wild_type_on = False
    if strain_name == wildtype_name:
        wild_type_on = False
    fig_median_sync, ax_median_sync = plt.subplots()
    max_time = 120
    num_cells_align = len(strain.signals["flavin/butter_alignTrunc_BY_births"])
    plot_title = f"Single-cell flavin fluorescence, aligned to first budding event\n(n = {num_cells_align})"
    if wild_type_on:
        median_plot(
            standardscaler.as_function(
                wildtype.signals["flavin/butter_alignTrunc_BY_births"]
            ).iloc[:, 0:max_time],
            unit_scaling=1 / 60,
            label="wild type",
            median_color="grey",
            error_color="lightgrey",
            xlabel="Time (hours)",
            plot_title=plot_title,
            ax=ax_median_sync,
        )
    median_plot(
        standardscaler.as_function(
            strain.signals["flavin/butter_alignTrunc_BY_births"]
        ).iloc[:, 0:max_time],
        unit_scaling=1 / 60,
        label=pretty_format(strain_name),
        median_color="#3714b0",
        error_color="#725ac3",
        xlabel="Time (hours)",
        plot_title=plot_title,
        ax=ax_median_sync,
    )

# Plot heatmap of flavin signals aligned by first birth
if True:
    fig_align, ax_align = plt.subplots()
    max_time = 120
    trace_df = strain.signals["flavin/butter_alignTrunc_BY_births"].dropna(how="all")
    buddings_df = strain.signals["births_alignTrunc_BY_births"].loc[
        trace_df.index.intersection(strain.signals["births_alignTrunc_BY_births"].index)
    ]
    trace_df, buddings_df = strain.sort_by_second_interval(trace_df, buddings_df)
    heatmap(
        trace_df=trace_df.iloc[:, 0:max_time],
        trace_name="flavin",
        buddings_df=buddings_df.iloc[:, 0:max_time],
        plot_title=f"Flavin fluorescence, \naligned by first budding",
        unit_scaling=1 / 12,
        xtick_step=2,
        xlabel="Time (hours)",
        ax=ax_align,
    )


# HISTOGRAMS

# Plot histogram of cycle lengths
if True:
    fig_hist, ax_hist = plt.subplots()
    plot_title = "Distribution of mean cycle lengths\n found by peak-to-peak distances"  # of each cell"
    num_cdcs = len(strain.lists["births_intervals"])
    histogram(
        5 * strain.lists["births_intervals"],
        label=f"cell division cycle (n = {num_cdcs})",
        color="#cb0077",
        lognormal=False,
        ax=ax_hist,
        plot_title=plot_title,
    )
    num_ymcs = len(strain.lists["flavin/butter_peaks_intervals"])
    histogram(
        5 * strain.lists["flavin/butter_peaks_intervals"],
        label=f"flavin cycle (n = {num_ymcs})",
        color="#3714b0",
        lognormal=False,
        ax=ax_hist,
        plot_title=plot_title,
    )
    tick_spacing = 60
    ax_hist.set_xlim((0, 420))
    ax_hist.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

# Histogram of peak periods from FFT
if True:
    fig_hist_pds, ax_hist_pds = plt.subplots()
    plot_title = "Distribution of mean cycle lengths\n derived from Fourier spectra"  # of each cell"
    num_ymcs = sum(~np.isnan(strain.lists["flavin/butter_freqs"]))

    def draw_hist_pds(ax):
        histogram(
            values=1 / (60 * strain.lists["flavin/butter_freqs"]),
            label=f"flavin cycle (n = {num_ymcs})",
            color="#3714b0",
            binsize=1 / 3,
            lognormal=False,
            xlabel="Period (hr)",
            ylabel="Cells",
            ax=ax,
            plot_title=plot_title,
        )
        if "mCherry" in strain.signals.keys():
            num_cdcs = sum(~np.isnan(strain.lists["mCherry/butter_freqs"]))
            histogram(
                values=1 / (60 * strain.lists["mCherry/butter_freqs"]),
                label=f"cell division cycle (n = {num_cdcs})",
                color="#cb0077",
                binsize=1 / 3,
                lognormal=False,
                xlabel="Period (hr)",
                ylabel="Cells",
                ax=ax,
                plot_title=plot_title,
            )

    draw_hist_pds(ax_hist_pds)

# Histogram of signal-to-noise ratio
if True:
    fig_snr, ax_snr = plt.subplots()
    num_snr = sum(~np.isnan(strain.lists["flavin/butter_snr"]))
    plot_title = f"Signal-to-noise ratios across cells \n(n = {num_snr}, cut-off frequency = {snr_cutoff:.3f})"
    histogram(
        values=strain.lists["flavin/butter_snr"],
        label=f"flavin time series",
        binsize=0.5,
        xlabel="Signal-to-noise ratio",
        ylabel="Cells",
        color="#3714b0",
        ax=ax_snr,
        plot_title=plot_title,
    )
    ax_snr.set_xlim((0, 10))


# SPECTRAL METHODS

# Drop non-oscillatory if haven't done
if False:
    targets = import_labels(
        filepath_targets,
        strain.signals["flavin/butter"].index.get_level_values("cellID"),
    )
    strain.signals["flavin/butter_osc"] = strain.signals["flavin/butter"].iloc[
        targets == 1
    ]
    strain.signals["births_osc"] = strain.signals["births"].loc[
        strain.signals["flavin/butter_osc"].index
    ]

# CROSS-CORRELATION
# Autocorrelation of flavin time series, compare strains if applicable
if True:
    fig_autocorr, ax_autocorr = plt.subplots()
    max_lag = 120
    num_cells_acf = len(strain.signals["flavin/butter_acf"])
    plot_title = f"Autocorrelation across population of flavin time series\n(n = {num_cells_acf})"
    if wild_type_on:
        median_plot(
            wildtype.signals["flavin/butter_acf"].iloc[:, 0:max_lag],
            unit_scaling=1 / 12,
            label="wild type",
            median_color="grey",
            error_color="lightgrey",
            xlabel="Lag (hours)",
            ylabel="Correlation",
            plot_title=plot_title,
            ax=ax_autocorr,
        )
    median_plot(
        strain.signals["flavin/butter_acf"].iloc[:, 0:max_lag],
        unit_scaling=1 / 12,
        label=pretty_format(strain_name),
        xlabel="Lag (hours)",
        ylabel="Correlation",
        median_color="#3714b0",
        error_color="#725ac3",
        plot_title=plot_title,
        ax=ax_autocorr,
    )
    set_origin_axes(ax_autocorr)

# Autocorrelation of flavin and mCherry time series, same strain
if True:
    fig_autocorr_comp, ax_autocorr_comp = plt.subplots(figsize=(6, 6))

    inset_histogram = True
    max_lag = 120

    num_cells_flavin_acf = len(strain.signals["flavin/butter_acf"])
    # plot_title = f"Autocorrelation across populations of time series (n = {num_cells_flavin_acf})"
    plot_title = ""

    if "mCherry" in strain.signals.keys():
        num_cells_mCherry_acf = len(strain.signals["mCherry/butter_acf"])
        median_plot(
            strain.signals["mCherry/butter_acf"].iloc[:, 0:max_lag],
            unit_scaling=1 / 12,
            label="mCherry",
            median_color="#cb0077",
            error_color="#d854a1",
            xlabel="Lag (hours)",
            ylabel="Population autocorrelation function",
            plot_title=plot_title,
            ax=ax_autocorr_comp,
        )
    else:
        print("Warning: no mCherry for autocorrelation plot")
    median_plot(
        strain.signals["flavin/butter_acf"].iloc[:, 0:max_lag],
        unit_scaling=1 / 12,
        label="flavin",
        median_color="#3714b0",
        error_color="#725ac3",
        xlabel="Lag (hours)",
        ylabel="Population autocorrelation function",
        plot_title=plot_title,
        ax=ax_autocorr_comp,
    )

    # draw x/y axes from origin
    set_origin_axes(ax_autocorr_comp)
    # shrink main axes height by 10% on bottom
    box = ax_autocorr_comp.get_position()
    ax_autocorr_comp.set_position(
        [box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.80]
    )
    # draw legend below main axes
    ax_autocorr_comp.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)

    if inset_histogram:
        # create inset
        ax_inset = ax_autocorr_comp.inset_axes([0.5, 0.6, 0.40, 0.40])
        # draw histogram
        draw_hist_pds(ax_inset)
        ax_inset.get_legend().remove()

# Crosscorrelation of flavin and mCherry time series, same strain
if False:
    fig_crosscorr, ax_crosscorr = plt.subplots()
    max_lag = 120
    num_cells_xcf = len(strain.signals["flavin_butterXmCherry_butter/xcf"])
    plot_title = f"Cross-correlation between flavin and mCherry \nacross population (n = {num_cells_xcf})"
    median_plot(
        strain.signals["flavin_butterXmCherry_butter/xcf"].loc[:, -max_lag:max_lag],
        unit_scaling=1 / 12,
        label="flavin vs mCherry cross-correlation",
        xlabel="Lag (hours)",
        ylabel="Correlation",
        median_color="grey",
        error_color="lightgrey",
        median_marker="o",
        plot_title=plot_title,
        ax=ax_crosscorr,
    )
    set_origin_axes(ax_crosscorr)


# FFT
# Mean plot
if True:
    fft_freqs = strain.signals["flavin/butter_fft"][0].dropna()
    fft_power = strain.signals["flavin/butter_fft"][1].dropna()
    fft_peak_freq = fft_freqs.iloc[0, np.argmax(np.mean(fft_power).to_numpy())]
    fig_fft, ax_fft = plt.subplots()
    mean_plot(
        fft_power,
        unit_scaling=fft_freqs.iloc[0, 1],
        label="",
        xlabel="Frequency ($\mathrm{min}^{-1}$)",
        ylabel="Power",
        mean_color="#3714b0",
        error_color="#725ac3",
        plot_title=f"Mean Fourier spectrum across all time series\n(n = {len(fft_power)})",
        ax=ax_fft,
    )
    ax_fft.axvline(
        fft_peak_freq,
        color="r",
        label=f"Peak power,\n corresponds to {1/fft_peak_freq:.2f} min period",
    )
    ax_fft.axvline(
        critical_freqs,
        color="g",
        label=f"Critical frequency for high-pass filter,\n corresponds to {1/critical_freqs:.2f} min period",
    )
    ax_fft.legend(loc="upper right")

# of mean aligned signal
if False:
    strain.signals[
        "flavin/butter_align_BY_flavin_butter_peaks_scaled"
    ] = standardscaler.as_function(
        strain.signals["flavin/butter_align_BY_flavin_butter_peaks"]
    )

    fft_freqs, fft_power = fft.as_function(
        strain.signals["flavin/butter_align_BY_flavin_butter_peaks_scaled"]
        .mean(axis=0)
        .to_frame()
        .T
    )
    fft_peak_freq = fft_freqs.iloc[0, np.argmax(np.mean(fft_power).to_numpy())]
    fig_fft, ax_fft = plt.subplots()
    num_cells_align = len(
        strain.signals["flavin/butter_align_BY_flavin_butter_peaks_scaled"]
    )
    ax_fft.plot(fft_freqs.iloc[0].to_numpy(), fft_power.iloc[0].to_numpy(), "b")
    ax_fft.set_xlabel("Frequency ($\mathrm{min}^{-1}$)")
    ax_fft.set_ylabel("Power")
    ax_fft.set_title(
        f"Fourier spectrum of mean flavin time series, \naligned by first flavin peak (n = {num_cells_align})"
    )
    ax_fft.axvline(
        fft_peak_freq,
        color="r",
        label=f"Peak power, corresponds to {1/fft_peak_freq:.2f} min period",
    )
    ax_fft.legend(loc="upper right")


# MULTIPLE STRAINS / CLUSTERING-RELATED
# FIXME: This best fits in a separate script.

if False:
    strain_list = ["FY4", "rim11_Del", "swe1_Del", "tsa1_Del_tsa2_Del"]
    all_flavin_processed, all["births"] = preprocessing(
        strain_name=strain_list,
        interval_start=25,
        interval_end=168,
    )

# PCA
if False:
    from postprocessor.core.processes.catch22 import catch22, catch22Parameters
    from sklearn.decomposition import PCA

    # 'borrowing' from umapper.py -- I can totally see me re-writing this into
    # a general embedding plotter thing later, definitely putting inheritance there
    from umapper import umap_embedding_plot
    from utils import generate_palette_map

    catch22processor = catch22(catch22Parameters.default())

    targets = import_labels(
        filepath_targets,
        all_flavin_processed.index.get_level_values("cellID"),
    )
    features = catch22processor.run(all_flavin_processed)

    # Drop non-oscillatory
    features = features.iloc[targets == 1]

    # Scale catch22 data matrix
    features_scaled = standardscaler.as_function(features.T).T

    fig_pca, ax_pca = plt.subplots()
    umap_embedding_plot(
        embedding=PCA().fit_transform(features_scaled),
        labels=features_scaled.index.get_level_values("strain"),
        palette_map=generate_palette_map(features_scaled),
        ax=ax_pca,
    )
    ax_pca.set_xlabel("PC1")
    ax_pca.set_ylabel("PC2")
    ax_pca.set_title("Pricipal component analysis")

# PDF
# shamelessly copied from https://stackoverflow.com/a/17788764

with PdfPages(f"{strain_name}_{experimentID}_plots.pdf") as pdf:
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
