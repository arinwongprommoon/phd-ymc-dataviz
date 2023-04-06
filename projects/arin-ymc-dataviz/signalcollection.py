#!/usr/bin/env python3

import numpy as np
from slicer import slice_signal
from utils import rearrange_by_sorted_list, find_nearest

# TODO: make this an abc?
class SignalCollection:
    """
    Collects processed dataframes (signals), types of data (sigtypes:
    'continuous' or 'binary'), along with relevant lists (lists: e.g.
    intervals between budding events).  Types of data informs the processes
    that are applied to them.
    """

    def __init__(self) -> None:
        self.signals = {}
        # 'continous' or 'binary'
        self.sigtypes = {}

        self.lists = {}

    def mask_to_intervals(self, mask_df):
        """Convert binary DataFrame, indicating events, to list of intervals

        Parameters
        ----------
        mask_df : pandas.DataFrame
            DataFrame containing elements of value 0 or 1, indicating whether an
            event occurs at each time point (columns) in multiple time series
            (rows).

        """
        temp_df = mask_df.apply(lambda x: np.diff(np.where(np.array(x) > 0)), axis=1)
        intervals = np.hstack(temp_df.to_list()).ravel()
        return np.array(intervals)

    def mask_to_mean_intervals(self, mask_df):
        """Convert binary DataFrame, indicating events, to list of mean intervals for each cell.

        Parameters
        ----------
        mask_df : pandas.DataFrame
            DataFrame containing elements of value 0 or 1, indicating whether an
            event occurs at each time point (columns) in multiple time series
            (rows).

        """
        temp_df = mask_df.apply(lambda x: np.diff(np.where(np.array(x) > 0)), axis=1)
        mean_intervals = []
        for el in temp_df.to_list():
            if not el.size:
                pass
            else:
                mean_intervals.append(np.mean(el))
        return np.array(mean_intervals)

    def get_second_intervals(self, mask_df):
        """Find intervals between first and second events from a mask DataFrame

        Find intervals between first and second events from a mask DataFrame
        Assumes: already aligned

        Parameters
        ----------
        mask_df : pandas.DataFrame
            binary DataFrame, functions as mask

        Return
        ------
        second_event_intervals : list
            list of intervals

        """
        second_event_intervals = []
        for index in mask_df.index:
            event_locs = np.where(mask_df.loc[index].to_numpy() == 1)[0]
            event_locs = np.delete(
                event_locs, 0
            )  # the first element of this list is always zero
            if event_locs.any():
                second_event_intervals.append(
                    event_locs[0]
                )  # this is when the 2nd event is in relation to the 1st
            else:
                second_event_intervals.append(None)
        # Absence of second event represented by nan
        second_event_intervals = np.array(second_event_intervals, dtype=float)

        return second_event_intervals

    def sort_by_second_interval(self, signal_df_aligned, mask_df_aligned):
        """Sort signal DataFrame and mask DataFrame by length of second interval of mask DataFrame

        Parameters
        ----------
        signal_df_aligned : pandas.DataFrame
            DataFrame with signals
        mask_df_aligned : pandas.DataFrame
            DataFrame with elements either 0 or 1, functions as mask

        """
        second_events_intervals = self.get_second_intervals(mask_df_aligned)
        signal_df_aligned_sorted = rearrange_by_sorted_list(
            signal_df_aligned, second_events_intervals
        )
        mask_df_aligned_sorted = rearrange_by_sorted_list(
            mask_df_aligned, second_events_intervals
        )
        return signal_df_aligned_sorted, mask_df_aligned_sorted

    # TODO: make this pythonic?
    def get_principal_freqs(self, fft_freqs_df, fft_power_df):
        """Get principal Fourier frequency from highest peak in periodogram

        Parameters
        ----------
        fft_freqs_df : pandas.DataFrame
            DataFrame showing in each row the frequency dimension of each
            Fourier spectrum
        fft_power_df : pandas.DataFrame
            DataFrame showing in each row the periodogram (Fourier spectrum)
        """
        fft_freqs_array = fft_freqs_df.to_numpy()
        fft_power_array = fft_power_df.to_numpy()
        peak_freqs = []
        for index, power_spectrum in enumerate(fft_power_array):
            argmax = np.argmax(power_spectrum)
            peak_freqs.append(fft_freqs_array[index, argmax])
        return np.array(peak_freqs)

    def get_snr(self, fft_freqs_df, fft_power_df, cutoff_freq):
        """Get signal-to-noise ratio from a Fourier spectrum

        Get signal-to-noise ratio from a Fourier spectrum. Defines a cut-off
        frequency; frequencies lower than this is considered signal, while
        frequencies higher than this is considered noise. The signal-to-noise
        ratio is defined as the area under the Fourier spectrum to the left of
        the cut-off divided by the area under the Fourier spectrum to the right
        of the cut-off. Follows:

        Parameters
        ----------
        fft_freqs_df : pandas.DataFrame
            DataFrame showing in each row the frequency dimension of each
            Fourier spectrum
        fft_power_df : pandas.DataFrame
            DataFrame showing in each row the periodogram (Fourier spectrum)
        cutoff_freq : float
            cut-off frequency to divide signal and noise
        """
        fft_freqs_array = fft_freqs_df.to_numpy()
        fft_power_array = fft_power_df.to_numpy()
        snr = []
        for rowindex, _ in enumerate(fft_power_array):
            cutoff_freq_nearest = find_nearest(
                fft_freqs_array[rowindex, :], cutoff_freq
            )
            # nans can occur if the origin time series has nans -- skip over these
            if np.isnan(cutoff_freq_nearest):
                snr.append(np.nan)
            else:
                cutoff_colindex = np.where(
                    fft_freqs_array[rowindex, :] == cutoff_freq_nearest
                )[0].item()
                area_all = np.trapz(
                    y=fft_power_array[rowindex, :], x=fft_freqs_array[rowindex, :]
                )
                area_signal = np.trapz(
                    y=fft_power_array[rowindex, 0:cutoff_colindex],
                    x=fft_freqs_array[rowindex, 0:cutoff_colindex],
                )
                area_noise = area_all - area_signal
                snr.append(area_signal / area_noise)
        return np.array(snr)


class RawSignalCollection(SignalCollection):
    """
    Collects raw signals from a grouper object, then defines names of
    signals based on abbrev dict.
    """

    def __init__(self, grouper, h5path_to_abbrev_dict) -> None:
        """define object

        Parameters
        ----------
        grouper : postprocessor.grouper.Grouper or NameGrouper object
            Grouper object that collects data from post-processed experiment
        h5path_to_abbrev_dict : dict
            Dictionary that specifies a lookup table -- keys: path in HDF5 file
            to get to signal of interest, values: abbreviation which is the
            signal name to write to in the RawSignalCollection object

        Examples
        --------
        FIXME: Add docs.

        """
        super().__init__()
        for h5path, (abbrev, sigtype) in h5path_to_abbrev_dict.items():
            self.signals[abbrev] = grouper.concat_signal(h5path)
            self.sigtypes[abbrev] = sigtype


class SlicedSignalCollection(SignalCollection):
    """
    Slices a raw signal collection, stores signals as appropriate
    """

    def __init__(
        self, raw, strain_name, interval_start, interval_end, classification_labels=None
    ) -> None:
        """defines how to slice

        Parameters
        ----------
        raw : signalcollection.RawSignalCollection object
            Source raw signal collection
        strain_name : string
            Name of strain of interest
        interval_start : int
            Starting time point
        interval_end : int
            Ending time point
        classification_labels : pandas.DataFrame
            DataFrame that shows classification labels, having columns:
            ['position', 'trap', 'cell_label', 'score'], and 'score'
            taking values of 0 or 1.  Pass a DataFrame here if you want
            0-scored time series to be removed.

        Examples
        --------
        FIXME: Add docs.

        """
        super().__init__()
        for signame, signal in raw.signals.items():
            # Slice by strain name and interval
            try:
                self.signals[signame] = slice_signal(
                    signal, strain_name, interval_start, interval_end
                )
                self.sigtypes[signame] = raw.sigtypes[signame]
            except KeyError:
                print(f"{signame} signal is not present in '{strain_name}' strain")
            if classification_labels is not None:
                self.signals[signame] = self.signals[signame].loc[
                    classification_labels[classification_labels.score == 1]
                    .set_index(["position", "trap", "cell_label"])
                    .index
                ]
