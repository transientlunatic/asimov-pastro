Coherence Module
================

.. automodule:: asimov_pastro.coherence
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

Detector Information
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: get_detector_location

.. autofunction:: get_antenna_pattern

Coherence Checks
^^^^^^^^^^^^^^^^

.. autofunction:: compute_time_delay_consistency

.. autofunction:: compute_amplitude_consistency

.. autofunction:: estimate_network_snr_from_likelihood

.. autofunction:: compute_coherence_checks

Module Constants
----------------

.. data:: DETECTOR_LOCATIONS

   Dictionary mapping interferometer codes to geocentric positions (meters).

   Available detectors:

   * ``H1``: LIGO Hanford
   * ``L1``: LIGO Livingston
   * ``V1``: Virgo
   * ``K1``: KAGRA

   Positions from LIGO-T980044.

.. data:: HAS_LAL

   Boolean flag indicating if LAL is available for antenna pattern calculations.

   If ``False``, simplified antenna patterns are used (less accurate).
