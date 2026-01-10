Pipeline Module
===============

.. automodule:: asimov_pastro
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: cli

Pastro Pipeline Class
---------------------

.. autoclass:: Pastro
   :members:
   :undoc-members:
   :show-inheritance:

   Asimov pipeline plugin for p_astro calculation.

   This class integrates with the Asimov workflow management system to:

   * Submit HTCondor jobs for p_astro computation
   * Monitor job completion
   * Collect results and metadata
   * Interface with other pipelines (PE, Omicron, etc.)

   Properties
   ^^^^^^^^^^

   .. autoattribute:: interferometers

   Pipeline Methods
   ^^^^^^^^^^^^^^^^

   .. automethod:: __init__

   .. automethod:: get_pe_samples

   .. automethod:: get_omicron_triggers

   .. automethod:: compute_coherence_checks

   .. automethod:: compute_pastro

   HTCondor Integration
   ^^^^^^^^^^^^^^^^^^^^

   .. automethod:: build_dag

   .. automethod:: submit_dag

   .. automethod:: detect_completion

   Asset Management
   ^^^^^^^^^^^^^^^^

   .. automethod:: collect_assets

   .. automethod:: collect_logs

   .. automethod:: after_completion
