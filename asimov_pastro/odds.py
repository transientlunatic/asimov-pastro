"""
Bayesian odds and p_astro calculation.

Implements the framework from:
  Ashton, Thrane & Smith (2019) arXiv:1909.11872
  Ashton & Thrane (2020)        arXiv:2006.05039

The noise hypothesis decomposes each detector into either containing a
glitch (modelled as an incoherent CBC) or Gaussian noise.  For N detectors
this gives 2^N terms (Eq. 9 of arXiv:1909.11872).  The odds ratio is then:

    O = xi * Z_S / Z_N

and p_astro = O / (1 + O).

Future work
-----------
* Replace the point-estimate xi with a proper marginalisation over the
  merger rate posterior using a log-normal prior (arXiv:2006.05039 §2).
* Infer xi_g from a full glitch population model (running PE on the top-N
  Omicron triggers and fitting a spin/mass hypermodel) rather than using a
  raw trigger count.
"""

import numpy as np
from itertools import product as cartesian_product
from scipy.special import logsumexp

# seconds per Julian year
_SECONDS_PER_YEAR = 365.25 * 24 * 3600


# ---------------------------------------------------------------------------
# Prior signal probability xi
# ---------------------------------------------------------------------------

def compute_xi(merger_rate, sensitive_volume, segment_duration):
    """
    Compute the prior probability xi that a segment contains a merger.

    Uses the formula (arXiv:2006.05039 Eq. 3 generalised):

        xi = R * V_T * delta_T

    where ``delta_T`` is the width of the geocentric coalescence-time prior
    used in the PE run (almost always 0.2 s for triggered analyses).

    .. note::
       Future improvement: marginalise over the rate posterior rather than
       using a point estimate for ``merger_rate``.  See arXiv:2006.05039 §2.

    Parameters
    ----------
    merger_rate : float
        Astrophysical merger rate in Gpc^-3 yr^-1.
    sensitive_volume : float
        Detector network sensitive volume in Gpc^3.
    segment_duration : float
        Width of the coalescence-time prior (seconds).  Default in blueprints
        is 0.2 s; can be read from the bilby ini file in future.

    Returns
    -------
    float
        Prior signal probability xi (dimensionless, << 1 for typical rates).
    """
    delta_t_years = segment_duration / _SECONDS_PER_YEAR
    return merger_rate * sensitive_volume * delta_t_years


# ---------------------------------------------------------------------------
# Per-detector glitch rate xi_g
# ---------------------------------------------------------------------------

def compute_xi_g_from_triggers(trigger_list, context_duration, segment_duration):
    """
    Estimate per-detector glitch probability per segment from Omicron triggers.

    Uses a Poisson point estimate:  xi_g = (N_triggers / T_context) * delta_T,
    i.e. the expected number of triggers in a single analysis segment.

    .. note::
       Future improvement: replace this with a full glitch population model
       (nested sampling on the top-N triggers + hyperparameter inference) as
       described in arXiv:1909.11872 Appendix B and arXiv:2006.05039 §3.

    Parameters
    ----------
    trigger_list : dict
        Mapping of IFO code to list of trigger file paths, as advertised by
        asimov-pyomicron under the ``"trigger list"`` asset key.
    context_duration : float
        Duration of the contextual data window in seconds (e.g. 86400 for
        24 hours as used in arXiv:2006.05039).
    segment_duration : float
        Duration of a single analysis segment in seconds (same ``delta_T``
        used for xi).

    Returns
    -------
    dict
        Mapping of IFO code to xi_g (float).
    """
    xi_g = {}
    for ifo, files in trigger_list.items():
        n_triggers = _count_triggers(files)
        rate = n_triggers / context_duration  # triggers per second
        xi_g[ifo] = rate * segment_duration
    return xi_g


def _count_triggers(files):
    """
    Return the total number of Omicron triggers across a list of files.

    Currently uses file count as a proxy.  A proper implementation would
    open each file (ROOT/HDF5/XML) and count rows above the SNR threshold.
    """
    # TODO: open trigger files and count rows, e.g. via gwpy.table.EventTable
    return len(files)


# ---------------------------------------------------------------------------
# Noise evidence (Eq. 9, arXiv:1909.11872)
# ---------------------------------------------------------------------------

def compute_log_noise_evidence(ifos, log_Z_glitch, log_Z_gaussian, xi_g):
    """
    Compute the log noise evidence by summing over all 2^N glitch subsets.

    For each possible assignment of detectors to {glitch, Gaussian noise},
    the contribution to the noise evidence is:

        term_S = prod_{x in S} xi_g[x] * Z_glitch[x]
               * prod_{x not in S} (1 - xi_g[x]) * Z_gaussian[x]

    The total noise evidence is the sum over all subsets S.  Uses
    ``logsumexp`` for numerical stability.

    Parameters
    ----------
    ifos : list of str
        IFO codes (e.g. ``["H1", "L1", "V1"]``).
    log_Z_glitch : dict
        Mapping IFO -> log single-detector evidence (glitch hypothesis, from
        single-IFO PE run).
    log_Z_gaussian : dict
        Mapping IFO -> log noise evidence (Gaussian hypothesis, i.e.
        ``log_noise_evidence`` from the same single-IFO PE run).
    xi_g : dict
        Mapping IFO -> per-segment glitch probability.

    Returns
    -------
    float
        Log noise evidence (log Z_N).
    """
    log_terms = []

    # Iterate over all 2^N assignments: True = glitch, False = Gaussian noise
    for glitch_flags in cartesian_product([False, True], repeat=len(ifos)):
        log_term = 0.0
        for ifo, is_glitch in zip(ifos, glitch_flags):
            p = xi_g[ifo]
            if is_glitch:
                log_term += np.log(p) + log_Z_glitch[ifo]
            else:
                log_term += np.log1p(-p) + log_Z_gaussian[ifo]
        log_terms.append(log_term)

    return float(logsumexp(log_terms))


# ---------------------------------------------------------------------------
# Final p_astro calculation
# ---------------------------------------------------------------------------

def compute_pastro(log_Z_S, log_Z_N, xi):
    """
    Compute the astrophysical probability p_astro.

    Implements:

        O = xi * Z_S / Z_N  =>  log O = log(xi) + log Z_S - log Z_N
        p_astro = O / (1 + O) = 1 / (1 + exp(-log O))

    Parameters
    ----------
    log_Z_S : float
        Log signal evidence from the coherent (multi-IFO) PE run.
    log_Z_N : float
        Log noise evidence from :func:`compute_log_noise_evidence`.
    xi : float
        Prior signal probability from :func:`compute_xi`.

    Returns
    -------
    p_astro : float
        Astrophysical probability in [0, 1].
    log_odds : float
        Natural log of the odds ratio O = xi * Z_S / Z_N.
    """
    log_odds = np.log(xi) + log_Z_S - log_Z_N
    p_astro = float(1.0 / (1.0 + np.exp(-log_odds)))
    return p_astro, float(log_odds)
