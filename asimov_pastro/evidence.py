"""
Read Bayesian evidence values from gravitational wave PE result files.

Supports two formats:

* **Bilby native HDF5** (``{label}_result.hdf5``) — evidence values are
  stored as scalar datasets at the root level: ``/log_evidence``,
  ``/log_noise_evidence``, ``/log_bayes_factor``.

* **PESummary HDF5** (``IGWN-GWTC*_PEDataRelease_*.h5``) — evidence values
  are stored under ``{approximant_label}/meta_data/sampler/`` with ``ln_``
  prefixed names.

The public entry point is :func:`read_evidence`, which auto-detects the
format.
"""

import h5py


def read_bilby_evidence(filepath):
    """
    Read evidence values from a bilby native HDF5 result file.

    Parameters
    ----------
    filepath : str or Path
        Path to a bilby result HDF5 file (e.g. ``{label}_result.hdf5``).

    Returns
    -------
    dict with keys:
        ``log_evidence``, ``log_noise_evidence``, ``log_bayes_factor``,
        ``log_evidence_err`` (may be ``None`` if not present).
    """
    with h5py.File(filepath, "r") as f:
        log_evidence = float(f["log_evidence"][()])
        log_noise_evidence = float(f["log_noise_evidence"][()])
        log_bayes_factor = float(f["log_bayes_factor"][()])
        log_evidence_err = (
            float(f["log_evidence_err"][()]) if "log_evidence_err" in f else None
        )
    return {
        "log_evidence": log_evidence,
        "log_noise_evidence": log_noise_evidence,
        "log_bayes_factor": log_bayes_factor,
        "log_evidence_err": log_evidence_err,
    }


def read_pesummary_evidence(filepath, label):
    """
    Read evidence values from a PESummary HDF5 file.

    Parameters
    ----------
    filepath : str or Path
        Path to a PESummary HDF5 file.
    label : str
        The approximant label to read from (e.g. ``"C01:IMRPhenomXPHM"``).
        Must be a key at the top level of the file that has a
        ``meta_data/sampler`` sub-group containing evidence values.

    Returns
    -------
    dict with keys:
        ``log_evidence``, ``log_noise_evidence``, ``log_bayes_factor``,
        ``log_evidence_err`` (may be ``None`` if not present).

    Raises
    ------
    KeyError
        If the label or required datasets are not found in the file.
    ValueError
        If the specified label does not contain sampler evidence values
        (e.g. it is a ``Mixed`` combined result).
    """
    with h5py.File(filepath, "r") as f:
        if label not in f:
            available = [k for k in f.keys() if k not in ("history", "version")]
            raise KeyError(
                f"Label '{label}' not found in {filepath}. "
                f"Available: {available}"
            )
        try:
            sampler = f[f"{label}/meta_data/sampler"]
        except KeyError:
            raise KeyError(
                f"No meta_data/sampler group found under '{label}' in {filepath}"
            )

        required = {"ln_evidence", "ln_noise_evidence", "ln_bayes_factor"}
        missing = required - set(sampler.keys())
        if missing:
            raise ValueError(
                f"Label '{label}' in {filepath} is missing evidence keys: "
                f"{missing}. This may be a combined 'Mixed' result that does "
                f"not carry individual evidence values."
            )

        log_evidence = float(sampler["ln_evidence"][()])
        log_noise_evidence = float(sampler["ln_noise_evidence"][()])
        log_bayes_factor = float(sampler["ln_bayes_factor"][()])
        log_evidence_err = (
            float(sampler["ln_evidence_error"][()])
            if "ln_evidence_error" in sampler
            else None
        )

    return {
        "log_evidence": log_evidence,
        "log_noise_evidence": log_noise_evidence,
        "log_bayes_factor": log_bayes_factor,
        "log_evidence_err": log_evidence_err,
    }


def detect_ifos_in_bilby_result(filepath):
    """
    Detect which IFOs are present in a bilby result file.

    Looks for IFO-prefixed SNR parameters in the posterior samples
    (e.g. ``H1_optimal_snr``, ``L1_optimal_snr``).

    Parameters
    ----------
    filepath : str or Path
        Path to a bilby native HDF5 result file.

    Returns
    -------
    list of str
        Sorted list of IFO codes found (e.g. ``["H1", "L1"]``).
    """
    known_ifos = ["H1", "L1", "V1", "K1", "ET", "CE"]
    found = []
    with h5py.File(filepath, "r") as f:
        if "posterior" in f:
            keys = set(f["posterior"].keys())
            for ifo in known_ifos:
                # Any IFO-prefixed parameter is sufficient
                if any(k.startswith(f"{ifo}_") for k in keys):
                    found.append(ifo)
    return sorted(found)


def read_evidence(filepath, label=None):
    """
    Read evidence values from a PE result file, auto-detecting the format.

    Parameters
    ----------
    filepath : str or Path
        Path to either a bilby native HDF5 result file or a PESummary HDF5
        file.
    label : str, optional
        Approximant label to use when reading a PESummary file.  If the file
        is in bilby native format this is ignored.  If the file is PESummary
        format and ``label`` is ``None``, the first non-metadata top-level
        group that carries evidence values is used.

    Returns
    -------
    dict with keys:
        ``log_evidence``, ``log_noise_evidence``, ``log_bayes_factor``,
        ``log_evidence_err``.
    """
    with h5py.File(filepath, "r") as f:
        is_bilby_native = "log_evidence" in f

    if is_bilby_native:
        return read_bilby_evidence(filepath)

    # PESummary format
    if label is not None:
        return read_pesummary_evidence(filepath, label)

    # No label given: try each top-level group until one has evidence values
    with h5py.File(filepath, "r") as f:
        candidates = [
            k for k in f.keys()
            if k not in ("history", "version") and isinstance(f[k], h5py.Group)
        ]

    for candidate in candidates:
        try:
            return read_pesummary_evidence(filepath, candidate)
        except (KeyError, ValueError):
            continue

    raise ValueError(
        f"Could not find evidence values in {filepath}. "
        f"Tried labels: {candidates}. "
        f"For PESummary files, specify the 'approximant' in the pastro config."
    )
