Quick Start
===========

This guide will get you started with asimov-pastro in 5 minutes.

Prerequisites
-------------

* Python 3.9 or later
* Posterior samples from a gravitational wave parameter estimation analysis
* Basic familiarity with GW data analysis

Installation
------------

.. code-block:: bash

   cd asimov-pastro
   pip install -e .

First Analysis
--------------

Using the Command-Line Tool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The fastest way to analyze posterior samples:

.. code-block:: bash

   ./pastro-coherence GW191215_posterior.h5 -i H1 -i L1 -i V1

You'll see output like:

.. code-block:: text

   ======================================================================
   Coherence Calculation
   ======================================================================
   Input:  GW191215_posterior.h5
   Output: GW191215_posterior_coherence.json

   Loading samples...
   ✓ Loaded 8884 samples
   ✓ Parameters: 60

   Using IFOs: H1, L1, V1

   ======================================================================
   RESULTS
   ======================================================================

   Time delay consistency: 1.0000  [GREEN]
   Amplitude consistency:  1.0000  [GREEN]
   Network SNR:            2.61
   Overall coherence:      0.8261  [GREEN]

   Interpretation: HIGH - Likely astrophysical signal

   ✓ Saved to GW191215_posterior_coherence.json

Using the Python API
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from asimov_pastro.calculator import PastroCalculator

   # Initialize calculator
   calc = PastroCalculator()

   # Load posterior samples
   samples = calc.load_samples(['GW191215_posterior.h5'])

   # Compute coherence checks
   result = calc.compute_coherence_checks(samples, ifos=['H1', 'L1', 'V1'])

   # Print results
   print(f"Time delay consistency: {result['time_delay_consistency']:.4f}")
   print(f"Amplitude consistency:  {result['amplitude_consistency']:.4f}")
   print(f"Overall coherence:      {result['overall_coherence']:.4f}")
   print(f"Interpretation:         {result['interpretation']}")

Understanding the Output
------------------------

The analysis produces several metrics:

Time Delay Consistency
^^^^^^^^^^^^^^^^^^^^^^

* **What it checks**: Whether timing differences between detectors match light travel time
* **High (>0.9)**: Consistent with gravitational wave propagation
* **Low (<0.5)**: Random timing, likely glitch

Amplitude Consistency
^^^^^^^^^^^^^^^^^^^^^

* **What it checks**: Whether amplitudes match antenna pattern predictions
* **High (>0.9)**: Consistent amplitude patterns across detectors
* **Low (<0.5)**: Inconsistent amplitudes, likely single-detector noise

Network SNR
^^^^^^^^^^^

* **What it is**: Combined signal-to-noise ratio from all detectors
* **>10**: Strong detection
* **5-10**: Moderate detection
* **<5**: Weak detection

Overall Coherence
^^^^^^^^^^^^^^^^^

Combined diagnostic score (weighted average of above metrics):

* **>0.8 (HIGH)**: Likely astrophysical signal ✓
* **0.5-0.8 (MODERATE)**: Uncertain classification
* **<0.5 (LOW)**: Likely instrumental noise ✗

Important Notes
---------------

.. warning::
   **These are diagnostic checks**, not the full Bayesian p_astro calculation
   described in Ashton & Thrane (2019, 2020). They provide preliminary quality
   assessment. The full Bayesian implementation is planned for future development.

Working with Different File Formats
------------------------------------

PESummary Files (GWTC)
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   ./pastro-coherence GWTC_file.h5 -l "C01:Mixed" -i H1 -i L1 -i V1

Bilby Files
^^^^^^^^^^^

.. code-block:: bash

   ./pastro-coherence bilby_posterior.h5 -i H1 -i L1

PyCBC Files
^^^^^^^^^^^

.. code-block:: bash

   ./pastro-coherence pycbc_posterior.h5 -i H1 -i L1

Next Steps
----------

* :doc:`cli_usage` - Full CLI documentation
* :doc:`examples` - More detailed examples
* :doc:`methodology` - Understanding the science
* :doc:`api/calculator` - Python API reference
