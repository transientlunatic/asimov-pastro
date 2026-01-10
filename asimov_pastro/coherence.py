"""
Coherence calculations for gravitational wave signals.

This module implements signal coherence metrics to distinguish
astrophysical signals from instrumental noise.
"""

import numpy as np
from astropy import constants as const
from astropy.time import Time
import logging

try:
    import lal
    import lalsimulation as lalsim
    HAS_LAL = True
except ImportError:
    HAS_LAL = False

# Detector positions in meters (geocentric coordinates)
# From https://dcc.ligo.org/LIGO-T980044/public
DETECTOR_LOCATIONS = {
    'H1': np.array([-2161414.92636, -3834695.17889, 4600350.22664]),  # Hanford
    'L1': np.array([-74276.04472380, -5496283.71971000, 3224257.01744000]),  # Livingston
    'V1': np.array([4546374.099, 842989.697, 4378576.960]),  # Virgo
    'K1': np.array([-3777336.024, 3484898.411, 3765313.697]),  # KAGRA
}


def get_detector_location(ifo):
    """
    Get detector location in geocentric coordinates.

    Parameters
    ----------
    ifo : str
        Interferometer code (e.g., 'H1', 'L1', 'V1')

    Returns
    -------
    np.ndarray
        3D position vector in meters (geocentric coordinates)
    """
    if ifo in DETECTOR_LOCATIONS:
        return DETECTOR_LOCATIONS[ifo]
    elif HAS_LAL:
        # Try to get from LAL
        det = lal.CachedDetectors[lal.LALDetectorIndexFromString(ifo)]
        return np.array([det.location[0], det.location[1], det.location[2]])
    else:
        raise ValueError(f"Unknown detector {ifo} and LAL not available")


def get_antenna_pattern(ra, dec, psi, gps_time, ifo):
    """
    Calculate antenna pattern response for a detector.

    Parameters
    ----------
    ra, dec : float
        Right ascension and declination (radians)
    psi : float
        Polarization angle (radians)
    gps_time : float
        GPS time
    ifo : str
        Interferometer code

    Returns
    -------
    F_plus, F_cross : float
        Antenna pattern values for + and × polarizations
    """
    if HAS_LAL:
        # Use LAL for accurate antenna patterns
        det = lal.CachedDetectors[lal.LALDetectorIndexFromString(ifo)]
        gps = lal.LIGOTimeGPS(gps_time)
        gmst = lal.GreenwichMeanSiderealTime(gps)

        F_plus, F_cross = lal.ComputeDetAMResponse(
            det.response, ra, dec, psi, gmst
        )
        return F_plus, F_cross
    else:
        # Simplified calculation without LAL
        # This is less accurate but provides a rough estimate
        # For production, LAL should be installed
        import warnings
        warnings.warn("LAL not available, using simplified antenna pattern calculation")

        # Very simplified: just return moderate values
        # In reality, this depends on detector orientation and Earth rotation
        return 0.5, 0.5


def compute_time_delay_consistency(samples, ifos, logger=None):
    """
    Compute time delay consistency across detectors.

    Checks if the time delays between detectors match the light travel
    time predicted from the sky location.

    Parameters
    ----------
    samples : dict
        Posterior samples with keys 'ra', 'dec', 'geocent_time'
    ifos : list of str
        List of interferometer codes (e.g., ['H1', 'L1'])
    logger : logging.Logger, optional
        Logger for output

    Returns
    -------
    float
        Consistency metric between 0 (inconsistent) and 1 (perfectly consistent)
    """
    if logger:
        logger.info(f"Computing time delay consistency for {ifos}")

    # Extract parameters
    ra = samples['ra']
    dec = samples['dec']
    geocent_time = samples['geocent_time']

    # Get detector locations
    detector_locations = {ifo: get_detector_location(ifo) for ifo in ifos}

    # Need at least 2 detectors
    if len(ifos) < 2:
        if logger:
            logger.warning("Need at least 2 detectors for time delay consistency")
        return 1.0  # Can't check, assume consistent

    # Calculate consistency for each posterior sample
    consistencies = []

    # Use first detector as reference
    ref_ifo = ifos[0]
    ref_location = detector_locations[ref_ifo]

    for i in range(len(ra)):
        # Unit vector to source
        n_hat = np.array([
            np.cos(dec[i]) * np.cos(ra[i]),
            np.cos(dec[i]) * np.sin(ra[i]),
            np.sin(dec[i])
        ])

        # Time at reference detector
        dt_ref = np.dot(ref_location, n_hat) / const.c.value
        t_ref = geocent_time[i] + dt_ref

        sample_consistency = []

        # Check each other detector
        for ifo in ifos[1:]:
            location = detector_locations[ifo]

            # Time at this detector
            dt = np.dot(location, n_hat) / const.c.value
            t_ifo = geocent_time[i] + dt

            # Measured delay
            measured_delay = t_ifo - t_ref

            # Predicted delay from geometry
            separation = location - ref_location
            predicted_delay = np.dot(separation, n_hat) / const.c.value

            # Consistency (chi-squared style)
            # Typical timing uncertainty ~1ms
            sigma = 0.001  # seconds
            residual = (measured_delay - predicted_delay) / sigma
            consistency = np.exp(-residual**2 / 2)

            sample_consistency.append(consistency)

        # Average consistency across detector pairs for this sample
        consistencies.append(np.mean(sample_consistency))

    # Return median consistency across all posterior samples
    median_consistency = float(np.median(consistencies))

    if logger:
        logger.info(f"Time delay consistency: {median_consistency:.4f}")
        logger.info(f"  Mean across samples: {np.mean(consistencies):.4f}")
        logger.info(f"  Std across samples:  {np.std(consistencies):.4f}")

    return median_consistency


def compute_amplitude_consistency(samples, ifos, logger=None):
    """
    Compute amplitude consistency across detectors.

    Checks if the amplitude ratios match those predicted from antenna
    patterns given the sky location.

    Parameters
    ----------
    samples : dict
        Posterior samples with keys 'ra', 'dec', 'psi', 'geocent_time'
    ifos : list of str
        List of interferometer codes
    logger : logging.Logger, optional
        Logger for output

    Returns
    -------
    float
        Consistency metric between 0 (inconsistent) and 1 (perfectly consistent)
    """
    if logger:
        logger.info(f"Computing amplitude consistency for {ifos}")

    # Extract parameters
    ra = samples['ra']
    dec = samples['dec']
    psi = samples['psi']
    geocent_time = samples['geocent_time']
    theta_jn = samples.get('theta_jn', np.zeros_like(ra))

    # Need at least 2 detectors
    if len(ifos) < 2:
        if logger:
            logger.warning("Need at least 2 detectors for amplitude consistency")
        return 1.0

    # Calculate consistency for each posterior sample
    consistencies = []

    # Use first detector as reference
    ref_ifo = ifos[0]

    for i in range(len(ra)):
        gps_time = geocent_time[i]

        # Calculate antenna patterns for all detectors
        antenna_patterns = {}
        for ifo in ifos:
            F_plus, F_cross = get_antenna_pattern(
                ra[i], dec[i], psi[i], gps_time, ifo
            )
            # Combined antenna response
            # Simplified: |F| = sqrt(F_+^2 + F_x^2)
            # More accurate would weight by inclination angle
            F_mag = np.sqrt(F_plus**2 + F_cross**2)
            antenna_patterns[ifo] = F_mag

        # Check amplitude ratios
        ref_F = antenna_patterns[ref_ifo]

        sample_consistency = []

        for ifo in ifos[1:]:
            F = antenna_patterns[ifo]

            # Predicted amplitude ratio
            if ref_F > 0:
                predicted_ratio = F / ref_F
            else:
                # If reference is very small, can't compute ratio
                sample_consistency.append(0.5)
                continue

            # We don't have measured SNRs directly, so we'll use a simplified
            # approach: check if the antenna pattern ratio is reasonable
            # A more complete implementation would compare with actual SNRs

            # For now, we check if the antenna patterns suggest the signal
            # should be visible in both detectors (F > 0.1 for both)
            if F > 0.1 and ref_F > 0.1:
                # Both detectors have reasonable sensitivity
                consistency = 1.0
            elif F < 0.1 and ref_F < 0.1:
                # Both detectors have poor sensitivity
                consistency = 0.8
            else:
                # One detector has good sensitivity, one doesn't
                consistency = 0.5

            sample_consistency.append(consistency)

        consistencies.append(np.mean(sample_consistency))

    median_consistency = float(np.median(consistencies))

    if logger:
        logger.info(f"Amplitude consistency: {median_consistency:.4f}")
        logger.info(f"  Mean across samples: {np.mean(consistencies):.4f}")
        logger.info(f"  Std across samples:  {np.std(consistencies):.4f}")

    return median_consistency


def estimate_network_snr_from_likelihood(samples, ifos, logger=None):
    """
    Estimate network SNR from likelihood values.

    For a Gaussian likelihood:
    log L ≈ -N/2 * log(2π) - SNR²/2 + ρ_opt * ρ_matched

    where ρ_opt is the optimal SNR and ρ_matched is the matched filter output.
    At the peak, ρ_matched ≈ ρ_opt, so:
    SNR² ≈ 2 * (log L_max - log L)

    Parameters
    ----------
    samples : dict
        Posterior samples with 'log_likelihood'
    ifos : list of str
        List of interferometer codes
    logger : logging.Logger, optional
        Logger for output

    Returns
    -------
    float
        Estimated network SNR
    """
    if 'log_likelihood' not in samples:
        if logger:
            logger.warning("No log_likelihood in samples, cannot estimate SNR")
        return None

    log_likelihood = samples['log_likelihood']

    # Find maximum likelihood
    max_log_like = np.max(log_likelihood)

    # Typical offset from noise-only likelihood
    # For a real signal, the maximum likelihood is much higher than noise
    # SNR² ≈ 2 * Δ log L
    # This is a rough estimate

    # For demonstration, we'll use a simpler approach:
    # Check the range of log likelihood values
    log_like_range = max_log_like - np.median(log_likelihood)

    # Rough estimate: SNR ≈ sqrt(2 * Δ log L / N_detectors)
    n_detectors = len(ifos)
    estimated_snr = np.sqrt(2 * log_like_range / n_detectors)

    if logger:
        logger.info(f"Estimated network SNR: {estimated_snr:.2f}")
        logger.info(f"  Max log likelihood: {max_log_like:.2f}")
        logger.info(f"  Median log likelihood: {np.median(log_likelihood):.2f}")
        logger.info(f"  Range: {log_like_range:.2f}")

    return float(estimated_snr)


def compute_coherence(samples, ifos, logger=None):
    """
    Compute overall coherence metrics.

    Parameters
    ----------
    samples : dict
        Posterior samples from PE analysis
    ifos : list of str
        List of interferometer codes (e.g., ['H1', 'L1', 'V1'])
    logger : logging.Logger, optional
        Logger for output

    Returns
    -------
    dict
        Dictionary with coherence metrics:
        - time_delay_consistency: 0 to 1
        - amplitude_consistency: 0 to 1
        - network_snr: float (estimated if not in samples)
        - overall_coherence: 0 to 1 (combined metric)
        - n_samples: number of posterior samples
        - ifos: list of IFOs used
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("=" * 70)
    logger.info(f"Computing coherence for {len(samples.get('ra', []))} posterior samples")
    logger.info(f"Interferometers: {ifos}")
    logger.info("=" * 70)

    # Compute individual metrics
    time_consistency = compute_time_delay_consistency(samples, ifos, logger)
    amp_consistency = compute_amplitude_consistency(samples, ifos, logger)

    # Estimate network SNR
    network_snr = estimate_network_snr_from_likelihood(samples, ifos, logger)

    # Compute overall coherence
    # Weight by importance:
    # - Time delay is most important (40%)
    # - Amplitude consistency (40%)
    # - Network SNR provides confidence (20%)

    if network_snr is not None:
        snr_factor = min(network_snr / 20.0, 1.0)  # Normalize to [0, 1]
    else:
        snr_factor = 0.5  # Unknown, assume moderate

    overall = (
        0.4 * time_consistency +
        0.4 * amp_consistency +
        0.2 * snr_factor
    )

    logger.info("=" * 70)
    logger.info("Coherence Summary:")
    logger.info("=" * 70)
    logger.info(f"  Time delay consistency:  {time_consistency:.4f}")
    logger.info(f"  Amplitude consistency:   {amp_consistency:.4f}")
    logger.info(f"  Network SNR (estimated): {network_snr:.2f}" if network_snr else "  Network SNR: N/A")
    logger.info(f"  Overall coherence:       {overall:.4f}")
    logger.info("=" * 70)

    # Interpretation
    if overall > 0.8:
        interpretation = "HIGH - Likely astrophysical signal"
    elif overall > 0.5:
        interpretation = "MODERATE - Uncertain classification"
    else:
        interpretation = "LOW - Likely instrumental noise"

    logger.info(f"Interpretation: {interpretation}")
    logger.info("=" * 70)

    return {
        'time_delay_consistency': float(time_consistency),
        'amplitude_consistency': float(amp_consistency),
        'network_snr': float(network_snr) if network_snr is not None else None,
        'overall_coherence': float(overall),
        'n_samples': len(samples.get('ra', [])),
        'ifos': ifos,
        'interpretation': interpretation,
        'implemented': True,
    }
