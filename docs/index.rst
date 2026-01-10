asimov-pastro: Astrophysical Probability for Gravitational Waves
===================================================================

**asimov-pastro** is an Asimov pipeline plugin for computing the astrophysical probability (p_astro) of gravitational wave candidates. It implements the Bayesian framework from:

* `Ashton & Thrane (2019) - arXiv:1909.11872 <https://arxiv.org/abs/1909.11872>`_: *Gravitational wave detection without boot straps: a Bayesian approach*
* `Ashton & Thrane (2020) - arXiv:2006.05039 <https://arxiv.org/abs/2006.05039>`_: *The astrophysical odds of GW151216*

The pipeline computes the probability that a gravitational wave candidate is of astrophysical origin versus terrestrial noise, enabling systematic classification of detections at scale.

.. note::
   **Current Status**: Diagnostic coherence checks are fully implemented and tested. The full Bayesian p_astro calculation is planned for future development.

Key Features
------------

* **PE-Tool Agnostic**: Works with posterior samples from any parameter estimation pipeline (Bilby, PyCBC, RIFT, etc.)
* **Format Support**: Reads PESummary (GWTC), Bilby, and PyCBC HDF5 formats
* **Diagnostic Checks**: Time delay and amplitude consistency checks across detectors
* **Standalone CLI**: Command-line tool for quick analysis
* **Asimov Integration**: Full pipeline plugin for workflow management
* **Scalable**: Designed for processing entire observing run catalogs

Quick Start
-----------

Installation
^^^^^^^^^^^^

.. code-block:: bash

   pip install -e .

Basic Usage
^^^^^^^^^^^

Using the standalone CLI:

.. code-block:: bash

   ./pastro-coherence GW191215_posterior.h5 -i H1 -i L1 -i V1

Using the Python API:

.. code-block:: python

   from asimov_pastro.calculator import PastroCalculator

   calc = PastroCalculator()
   samples = calc.load_samples(['posterior.h5'])
   checks = calc.compute_coherence_checks(samples, ifos=['H1', 'L1', 'V1'])

   print(f"Overall coherence: {checks['overall_coherence']:.4f}")
   print(f"Interpretation: {checks['interpretation']}")

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   cli_usage
   methodology
   examples

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/calculator
   api/coherence
   api/pipeline

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   roadmap
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
