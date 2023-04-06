#!/usr/bin/env python3

from postprocessor.core.multisignal.align import align


def apply_postprocesses(signalcollection, process_dict):
    """Apply postprocesses

    Parameters
    ----------
    signalcollection : signalcollection.SignalCollection object
        Collect of signals to apply processes to
    process_dict : dict
        Dictionary of postprocesses to apply. Keys: string to append to name of
        source signal. Values: dictionary that specifies process, parameters,
        type of signal to operate on, type of signal of output.

    Examples
    --------
    FIXME: Add docs.

    """
    # Goes through each appending string, which represents a 'key' for each post-
    # process
    for process_appending_string, options_dict in process_dict.items():
        # Defines a list of signal (names) that I want my post-process to run on,
        # based on the ending of the signal name.
        applied_signames = [
            signame
            for signame in list(signalcollection.signals.keys())
            if signame.endswith(options_dict["signame_endswith"])
        ]
        for signame in applied_signames:
            # New signal name: append string that indicates which post-process
            # created it
            newname = signame + process_appending_string
            # Run if the signal type (continuous or binary) is correct
            if signalcollection.sigtypes[signame] == options_dict["input_sigtype"]:
                signalcollection.signals[newname] = options_dict["runner"].run(
                    signalcollection.signals[signame]
                )
                # Specify signal type of new signal
                signalcollection.sigtypes[newname] = options_dict["output_sigtype"]
            else:
                print(
                    f"{signame} is not {options_dict['input_sigtype']}, skipping {process_appending_string}"
                )


def apply_multisignal_postprocesses(signalcollection, multisignal_process_dict):
    """Apply multisignal postprocesses

    Parameters
    ----------
    signalcollection : signalcollection.SignalCollection object
        Collect of signals to apply processes to
    multisignal_process_dict : dict
        Dictionary of postprocesses to apply. Keys: string to append to name of
        source signal. Values: dictionary that specifies process, parameters,
        the names of the two signals to operate on, type of signal of output.

    Examples
    --------
    FIXME: Add docs.

    """

    # Goes through each appending string, which represents a 'key' for each post-
    # process
    for process_appending_string, options_dict in multisignal_process_dict.items():
        newname = (
            options_dict["input0"].replace("/", "_")
            + "X"
            + options_dict["input1"].replace("/", "_")
            + process_appending_string
        )
        # Run process
        signalcollection.signals[newname] = options_dict["runner"].run(
            signalcollection.signals[options_dict["input0"]],
            signalcollection.signals[options_dict["input1"]],
        )
        signalcollection.sigtypes[newname] = options_dict["output_sigtype"]


def apply_align_postprocesses(signalcollection, align_process_dict):
    """Apply align postprocesses

    Parameters
    ----------
    signalcollection : signalcollection.SignalCollection object
        Collect of signals to apply processes to
    process_dict : dict
        Dictionary of postprocesses to apply. Keys: ints.
        Values: dictionary that specifies
        - parameters
        - first input signal name
        - second input signal name
        - whether to store first output signal
        - type of first output signal
        - whether to store second output signal
        - type of second output signal
        - infix for resulting signal name

    Examples
    --------
    FIXME: Add docs.

    """
    for options_dict in align_process_dict.values():
        runner = align(options_dict["parameters"])
        output0, output1 = runner.run(
            signalcollection.signals[options_dict["input0"]],
            signalcollection.signals[options_dict["input1"]],
        )
        if options_dict["output0"]:
            newname0 = (
                options_dict["input0"]
                + options_dict["infix"]
                + options_dict["input1"].replace("/", "_")
            )
            signalcollection.signals[newname0] = output0
            signalcollection.sigtypes[newname0] = options_dict["output0_sigtype"]

        if options_dict["output1"]:
            newname1 = (
                options_dict["input1"]
                + options_dict["infix"]
                + options_dict["input1"].replace("/", "_")
            )
            signalcollection.signals[newname1] = output1
            signalcollection.sigtypes[newname1] = options_dict["output1_sigtype"]


def apply_all_processes(
    signalcollection, process_dict, multisignal_process_dict, align_process_dict
):
    apply_postprocesses(signalcollection, process_dict)
    try:
        apply_multisignal_postprocesses(signalcollection, multisignal_process_dict)
    except:
        pass
    apply_align_postprocesses(signalcollection, align_process_dict)
