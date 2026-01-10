Methodology
===========

Overview
--------

asimov-pastro implements a Bayesian framework for computing the astrophysical probability (p_astro) of gravitational wave candidates. The method distinguishes between:

* **Astrophysical signals**: Coherent gravitational waves from compact binary mergers
* **Terrestrial noise**: Instrumental glitches and detector artifacts

Theoretical Framework
---------------------

The astrophysical probability is defined as:

.. math::

   p_{\text{astro}} = \frac{1}{1 + \mathcal{O}_N^A}

where :math:`\mathcal{O}_N^A` is the odds ratio of the terrestrial noise hypothesis to the astrophysical signal hypothesis:

.. math::

   \mathcal{O}_N^A = \frac{P(\text{data} | \text{noise})}{P(\text{data} | \text{signal})}

This framework is described in:

* **Ashton & Thrane (2019)** - `arXiv:1909.11872 <https://arxiv.org/abs/1909.11872>`_
* **Ashton & Thrane (2020)** - `arXiv:2006.05039 <https://arxiv.org/abs/2006.05039>`_

Current Implementation
----------------------

.. note::
   **Implementation Status**: The current version provides **diagnostic coherence checks** that assess signal consistency across detectors. The full Bayesian p_astro calculation is planned for future development.

Diagnostic Coherence Checks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pipeline currently implements three diagnostic checks:

1. Time Delay Consistency
~~~~~~~~~~~~~~~~~~~~~~~~~~

Checks if timing differences between detectors match the light travel time predicted from the sky location.

**Method**:

For each posterior sample with sky location :math:`(\alpha, \delta)`:

1. Calculate unit vector to source:

   .. math::

      \hat{n} = \begin{pmatrix}
      \cos\delta \cos\alpha \\
      \cos\delta \sin\alpha \\
      \sin\delta
      \end{pmatrix}

2. Calculate predicted time delay between detectors:

   .. math::

      \Delta t_{\text{pred}} = \frac{(\vec{r}_2 - \vec{r}_1) \cdot \hat{n}}{c}

   where :math:`\vec{r}_1, \vec{r}_2` are detector positions.

3. Compare with measured delay from posterior samples.

4. Consistency metric:

   .. math::

      \chi_{\text{time}} = \exp\left(-\frac{(\Delta t_{\text{meas}} - \Delta t_{\text{pred}})^2}{2\sigma^2}\right)

   with :math:`\sigma \approx 1` ms (typical timing uncertainty).

**Interpretation**:

* **High** (:math:`> 0.9`): Time delays consistent with GW propagation at speed of light
* **Low** (:math:`< 0.5`): Random timing, likely glitch

**Reference**: Standard technique used in coherent GW analysis pipelines (cWB, LIB, etc.)

2. Amplitude Consistency
~~~~~~~~~~~~~~~~~~~~~~~~~

Checks if amplitude ratios across detectors match antenna pattern predictions.

**Method**:

1. Calculate antenna patterns for each detector:

   .. math::

      F_+ = \frac{1}{2}(1 + \cos^2\iota)\cos 2\psi D_+(\alpha, \delta, t)

   .. math::

      F_\times = \cos\iota \sin 2\psi D_\times(\alpha, \delta, t)

   where :math:`D_+, D_\times` depend on detector orientation and Earth rotation.

2. Combined antenna response:

   .. math::

      |F| = \sqrt{F_+^2 + F_\times^2}

3. Check if all detectors have reasonable sensitivity for the given sky location.

**Note**: Full amplitude consistency would compare predicted SNR ratios:

.. math::

   \frac{\rho_2}{\rho_1} \approx \frac{|F_2|}{|F_1|}

This requires per-detector SNR values from the PE analysis.

**Interpretation**:

* **High**: Amplitude patterns consistent across detectors
* **Low**: Inconsistent amplitudes, likely single-detector glitch

**Reference**: Antenna pattern theory from LIGO-T010110

3. Network SNR Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Estimates the combined network signal-to-noise ratio.

**Method**:

From likelihood values in posterior samples:

.. math::

   \text{SNR}_{\text{net}} \approx \sqrt{\frac{2\Delta\log\mathcal{L}}{N_{\text{det}}}}

where :math:`\Delta\log\mathcal{L}` is the likelihood range.

Alternatively, if per-detector SNRs are available:

.. math::

   \text{SNR}_{\text{net}} = \sqrt{\sum_i \rho_i^2}

**Interpretation**:

* **> 10**: Strong detection
* **5-10**: Moderate detection
* **< 5**: Weak detection

Overall Coherence Score
~~~~~~~~~~~~~~~~~~~~~~~~

The diagnostic checks are combined into an overall score:

.. math::

   C_{\text{overall}} = 0.4\chi_{\text{time}} + 0.4\chi_{\text{amp}} + 0.2\min\left(\frac{\text{SNR}}{20}, 1\right)

**Classification**:

* **HIGH** (:math:`> 0.8`): Likely astrophysical signal
* **MODERATE** (:math:`0.5 - 0.8`): Uncertain classification
* **LOW** (:math:`< 0.5`): Likely instrumental noise

Important Distinctions
----------------------

Diagnostic Checks vs Full p_astro
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
   The current coherence checks are **diagnostic tools**, not the full Bayesian p_astro calculation described in the papers.

**What we have**: Heuristic checks that qualitatively assess coherence

**What the papers describe**: Full Bayesian model comparison with proper marginalization

The full implementation will compute:

.. math::

   \mathcal{O}_N^A = \frac{\int P(\text{data}|\vec{\theta}_N) P(\vec{\theta}_N) d\vec{\theta}_N}{\int P(\text{data}|\vec{\theta}_A) P(\vec{\theta}_A) d\vec{\theta}_A}

where:

* :math:`\vec{\theta}_A`: Astrophysical signal parameters
* :math:`\vec{\theta}_N`: Noise/glitch parameters
* Full likelihood evaluation in both hypotheses

Single vs Multi-Detector PE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. important::
   The pipeline uses **one multi-detector PE run**, not separate single-detector runs.

The multi-detector posterior samples already contain all necessary information:

* Individual detector arrival times (from geometry + geocent_time)
* Expected amplitudes (from antenna patterns)
* Coherence is implicit in how well parameters explain all detectors

See :doc:`FAQ <faq>` for detailed explanation.

Future Development
------------------

Planned Features
^^^^^^^^^^^^^^^^

1. **Full Bayesian p_astro**:

   * Implement complete likelihood models
   * Astrophysical signal hypothesis (coherent)
   * Terrestrial noise hypothesis (incoherent)
   * Proper marginalization over parameters

2. **Glitch Characterization**:

   * Integrate Omicron trigger information
   * Time-frequency glitch modeling
   * Detector-specific noise profiles

3. **Population Analysis**:

   * Multi-event p_astro computation
   * Selection effects
   * Astrophysical rate inference

4. **Null Stream Analysis**:

   * Advanced coherence test
   * Combine data in null configuration
   * Check for excess power

References
----------

**Primary Papers**:

* Ashton & Thrane (2019), *Gravitational wave detection without boot straps: a Bayesian approach*, `arXiv:1909.11872 <https://arxiv.org/abs/1909.11872>`_
* Ashton & Thrane (2020), *The astrophysical odds of GW151216*, `arXiv:2006.05039 <https://arxiv.org/abs/2006.05039>`_

**Related Methods**:

* Biscoveanu et al. (2023), *Unified p_astro for gravitational waves*, `arXiv:2305.00071 <https://arxiv.org/abs/2305.00071>`_
* Anderson et al. (2001), *Excess power statistic for detection of burst sources*, PRD 63, 042003
* Klimenko et al. (2016), *Method for detection and reconstruction of gravitational wave transients*, PRD 93, 042004

**Technical References**:

* LIGO-T010110: *Antenna patterns and detector responses*
* LIGO-T980044: *Detector positions and orientations*
