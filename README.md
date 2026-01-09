# asimov-pastro

An Asimov pipeline plugin for computing the astrophysical probability (p_astro) of gravitational wave candidates.

## Overview

This pipeline implements the Bayesian framework described in:
- [arXiv:1909.11872](https://arxiv.org/abs/1909.11872) - Gravitational wave detection without boot straps: a Bayesian approach
- [arXiv:2006.05039](https://arxiv.org/abs/2006.05039) - The astrophysical odds of GW151216

The pipeline computes p_astro, the probability that a gravitational wave candidate is of astrophysical origin (vs. terrestrial noise), using:
1. **Posterior samples** from parameter estimation (PE) analyses
2. **Signal coherence** across multiple detectors
3. **Glitch characterization** from Omicron triggers (optional)
4. **Bayesian odds calculation** comparing astrophysical and terrestrial hypotheses

## Key Features

- **PE-tool agnostic**: Works with posterior samples from any PE pipeline
  - bilby / bilby_pipe
  - pycbc_inference
  - RIFT
  - LALInference
  - Any tool that outputs standard formats (HDF5, JSON)

- **Omicron integration**: Can incorporate instrumental glitch information from asimov-pyomicron

- **Scalable**: Integrates with Asimov's workflow management for processing many events

- **Multi-event support**: Can be used with Asimov's ProjectAnalysis for population-level studies

## Installation

```bash
# Clone the repository
cd /path/to/asimov-pipelines/asimov-pastro

# Install in development mode
pip install -e .
```

The pipeline will automatically register with Asimov via entry points.

## Quick Start

### 1. Set up your Asimov project

```bash
mkdir my_pastro_project
cd my_pastro_project
asimov init "P_astro Analysis"
```

### 2. Add an event and PE analysis

First, run your parameter estimation analysis (e.g., with bilby_pipe or pycbc_inference):

```yaml
# In your ledger file
productions:
  - name: "GW150914_bilby_pe"
    pipeline: "bilby"
    status: "complete"
    # ... PE configuration ...
```

### 3. Add p_astro analysis

Add a new production that depends on the PE analysis:

```yaml
productions:
  - name: "GW150914_pastro"
    pipeline: "pastro"
    needs:
      - "GW150914_bilby_pe"  # PE analysis dependency
    interferometers:
      - H1
      - L1
    scheduler:
      request_memory: "2GB"
      request_disk: "1GB"
      accounting_group: "ligo.dev.o4.cbc.pe.bayesian"
```

### 4. Run the analysis

```bash
asimov apply
```

The pipeline will:
1. Find posterior samples from the dependent PE analysis
2. Compute signal coherence across detectors
3. Calculate p_astro
4. Store results in `pastro_results.json`

## Architecture

### Pipeline Flow

```
Input: PE Posterior Samples + (Optional) Omicron Triggers
        ↓
Step 1: Load samples (format-agnostic reader)
        ↓
Step 2: Compute signal coherence
        ↓
Step 3: Characterize noise environment
        ↓
Step 4: Compute Bayesian odds ratio
        ↓
Output: p_astro value + diagnostics
```

### Integration with Other Pipelines

```
Event Detection
    ↓
gwdata (Frame files, PSDs, Calibration)
    ↓
[Optional] asimov-pyomicron (Glitch triggers) ──┐
    ↓                                             │
PE Pipeline (bilby/pycbc/RIFT/etc.)             │
    ↓                                             │
Posterior Samples                                │
    ↓                                             │
asimov-pastro ←──────────────────────────────────┘
    ↓
p_astro Classification
```

## Output Format

The pipeline produces `pastro_results.json`:

```json
{
  "pastro": 0.95,
  "coherence": {
    "time_delay_consistency": 0.98,
    "amplitude_consistency": 0.92,
    "phase_consistency": 0.89,
    "network_snr": 23.4
  },
  "samples_info": {
    "n_samples": 10000,
    "parameters": ["mass_1", "mass_2", "luminosity_distance", ...]
  }
}
```

The p_astro value is also stored in the production metadata for easy access:

```python
production.meta['pastro']  # Access in Python
```

## Configuration Options

### Required Metadata

- `interferometers`: List of detectors (e.g., `["H1", "L1", "V1"]`)
- `needs`: Dependency on PE analysis that provides posterior samples

### Optional Metadata

- `scheduler.request_memory`: Memory allocation (default: "2GB")
- `scheduler.request_disk`: Disk space (default: "1GB")
- `scheduler.accounting_group`: HTCondor accounting group

### Advanced Configuration

For custom coherence calculations or noise models, you can extend the `PastroCalculator` class:

```python
from asimov_pastro.calculator import PastroCalculator

class CustomPastroCalculator(PastroCalculator):
    def compute_coherence(self, samples):
        # Custom coherence implementation
        pass

    def compute_pastro(self, samples, coherence, glitch_info):
        # Custom p_astro calculation
        pass
```

## Development Status

### Implemented ✓

- [x] Asimov pipeline interface
- [x] Format-agnostic sample reader (HDF5, JSON)
- [x] Omicron trigger ingestion
- [x] HTCondor job submission
- [x] Results collection and metadata storage

### In Progress / TODO

- [ ] **Coherence calculations**
  - [ ] Time delay consistency across IFOs
  - [ ] Amplitude consistency given antenna patterns
  - [ ] Phase relationship verification
  - [ ] Network SNR computation

- [ ] **Bayesian odds calculation**
  - [ ] Astrophysical signal likelihood model
  - [ ] Terrestrial noise likelihood model
  - [ ] Prior on astrophysical rate
  - [ ] Integration over parameters

- [ ] **Noise characterization**
  - [ ] Glitch parameterization from Omicron triggers
  - [ ] Time-frequency glitch modeling
  - [ ] Detector-specific noise profiles

- [ ] **Validation and testing**
  - [ ] Unit tests for all components
  - [ ] Integration tests with real PE outputs
  - [ ] Comparison with published p_astro values

- [ ] **Documentation**
  - [ ] Mathematical derivation documentation
  - [ ] API reference
  - [ ] Tutorial notebooks

## Contributing

This is an active development project. Contributions are welcome!

### Priority Areas

1. **Coherence calculations**: Implement the methods from arXiv:1909.11872
2. **Bayesian framework**: Full implementation of the odds ratio calculation
3. **Validation**: Test against known results (GW150914, GW151216, etc.)
4. **Documentation**: Mathematical background and usage examples

### Development Setup

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run tests with coverage
pytest --cov=asimov_pastro --cov-report=html
```

## References

- Ashton, G., & Thrane, E. (2019). *Gravitational wave detection without boot straps: a Bayesian approach*. [arXiv:1909.11872](https://arxiv.org/abs/1909.11872)
- Ashton, G., & Thrane, E. (2020). *The astrophysical odds of GW151216*. [arXiv:2006.05039](https://arxiv.org/abs/2006.05039)

## License

This project follows the licensing of the parent Asimov project.

## Contact

For questions or issues, please open an issue on the repository or contact the developers.
