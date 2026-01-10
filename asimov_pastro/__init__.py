"""
Asimov pipeline for computing astrophysical probability (p_astro) of
gravitational wave candidates.

This pipeline implements the Bayesian framework described in:
- arXiv:1909.11872 - Gravitational wave detection without boot straps
- arXiv:2006.05039 - The astrophysical odds of GW151216

The pipeline can consume posterior samples from any parameter estimation
tool (bilby, pycbc_inference, RIFT, etc.) and compute the probability
that a candidate is of astrophysical origin vs. terrestrial noise.
"""

import os
import glob
import json
from pathlib import Path

import asimov.pipeline
from asimov.utils import set_directory

try:
    import htcondor2 as htcondor
    import classad2 as classad
except ImportError:
    import htcondor
    import classad

from asimov import config


class Pastro(asimov.pipeline.Pipeline):
    """
    Asimov pipeline for computing p_astro (astrophysical probability).

    This pipeline:
    1. Ingests posterior samples from any PE pipeline (bilby, pycbc, etc.)
    2. Optionally ingests Omicron glitch triggers for noise characterization
    3. Computes signal coherence across detectors
    4. Calculates Bayesian odds: P(astrophysical) / P(terrestrial)
    5. Outputs p_astro value and diagnostic information

    The pipeline is agnostic to the source of posterior samples and can
    work with outputs from bilby_pipe, pycbc_inference, RIFT, or any
    other parameter estimation tool that produces standard formats.
    """

    name = "pastro"

    def __init__(self, production, category=None):
        """
        Initialize the p_astro pipeline.

        Parameters
        ----------
        production : asimov.Production
            The production object containing analysis metadata
        category : str, optional
            The category of this analysis
        """
        super().__init__(production, category)
        self.logger.info(f"Initializing p_astro pipeline for {production.name}")

    @property
    def interferometers(self):
        """Get list of interferometers for this analysis."""
        ifos = self.production.meta.get("interferometers", [])
        if isinstance(ifos, dict):
            return list(ifos.keys())
        return ifos

    def get_pe_samples(self):
        """
        Retrieve posterior samples from dependent PE analyses.

        This method is agnostic to the PE tool used - it can read
        samples from bilby, pycbc, RIFT, or any other pipeline.

        Returns
        -------
        list of str
            Paths to posterior sample files
        """
        dependencies = self.production._previous_assets()
        sample_files = dependencies.get('samples', [])

        if not sample_files:
            self.logger.warning("No posterior samples found in dependencies")
            return []

        self.logger.info(f"Found {len(sample_files)} posterior sample file(s)")
        for sample_file in sample_files:
            self.logger.info(f"  - {sample_file}")

        return sample_files

    def get_omicron_triggers(self):
        """
        Retrieve Omicron glitch triggers if available.

        Returns
        -------
        dict or None
            Dictionary mapping IFOs to trigger files, or None if not available
        """
        dependencies = self.production._previous_assets()

        # Look for Omicron outputs in various formats
        omicron_triggers = {}

        # Check for direct omicron trigger files
        if 'omicron_triggers' in dependencies:
            omicron_triggers = dependencies['omicron_triggers']

        # Check for generic trigger files
        elif 'triggers' in dependencies:
            omicron_triggers = dependencies['triggers']

        if omicron_triggers:
            self.logger.info(f"Found Omicron triggers for IFOs: {list(omicron_triggers.keys())}")
        else:
            self.logger.info("No Omicron triggers found (optional)")

        return omicron_triggers if omicron_triggers else None

    def compute_coherence(self, samples):
        """
        Compute signal coherence across detectors.

        Compares signal parameters across IFOs to determine if they're
        consistent with a single astrophysical source.

        Parameters
        ----------
        samples : dict
            Posterior samples organized by parameter

        Returns
        -------
        dict
            Coherence metrics (time delays, amplitude ratios, network SNR, etc.)
        """
        from .coherence import compute_coherence as compute_coherence_impl

        self.logger.info("Computing signal coherence across detectors")

        # Get IFOs from production metadata
        ifos = self.interferometers

        # Call the implementation
        coherence_metrics = compute_coherence_impl(samples, ifos, self.logger)

        return coherence_metrics

    def compute_pastro(self, samples, coherence, glitch_triggers=None):
        """
        Compute the astrophysical probability using Bayesian odds.

        Implements the method from arXiv:1909.11872 and arXiv:2006.05039:

        p_astro = 1 / (1 + O_N^A)

        where O_N^A is the odds ratio of terrestrial noise to astrophysical signal.

        Parameters
        ----------
        samples : dict
            Posterior samples from PE analysis
        coherence : dict
            Coherence metrics from compute_coherence()
        glitch_triggers : dict, optional
            Omicron trigger information for noise characterization

        Returns
        -------
        float
            The astrophysical probability (0 to 1)
        """
        self.logger.info("Computing p_astro using Bayesian odds")

        # TODO: Implement full Bayesian calculation
        # 1. Model astrophysical signal hypothesis (coherent signal)
        # 2. Model terrestrial noise hypothesis (glitches, incoherent)
        # 3. Compute likelihood ratio
        # 4. Convert to probability

        # Placeholder: return 0.5 (maximum uncertainty)
        pastro = 0.5

        self.logger.warning("p_astro calculation not yet implemented, returning placeholder")
        return pastro

    def build_dag(self, dryrun=False):
        """
        Build the HTCondor job for p_astro calculation.

        Parameters
        ----------
        dryrun : bool
            If True, don't actually submit the job
        """
        self.logger.info("Building p_astro calculation job")

        name = self.production.name
        rundir = self.production.rundir

        # Get input data
        sample_files = self.get_pe_samples()
        omicron_triggers = self.get_omicron_triggers()

        if not sample_files:
            raise ValueError("No posterior samples available for p_astro calculation")

        # Create a Python script that will run the calculation
        script_path = os.path.join(rundir, f"{name}_calculate_pastro.py")

        script_content = f"""#!/usr/bin/env python
\"\"\"
P_astro calculation script for {name}
Generated by asimov-pastro pipeline
\"\"\"

import sys
import json
from asimov_pastro.calculator import PastroCalculator

# Input files
sample_files = {sample_files}
omicron_triggers = {omicron_triggers}

# Run calculation
calculator = PastroCalculator()
results = calculator.calculate(
    sample_files=sample_files,
    omicron_triggers=omicron_triggers,
)

# Save results
with open('pastro_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"p_astro = {{results['pastro']:.4f}}")
sys.exit(0)
"""

        os.makedirs(rundir, exist_ok=True)
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)

        # Build HTCondor submit description
        executable = script_path

        description = {
            "executable": executable,
            "output": f"{name}.out",
            "error": f"{name}.err",
            "log": f"{name}.log",
            "request_disk": self.production.meta.get("scheduler", {}).get("request disk", "1GB"),
            "request_memory": self.production.meta.get("scheduler", {}).get("request memory", "2GB"),
            "request_cpus": self.production.meta.get("scheduler", {}).get("cpus", "1"),
            "batch_name": f"pastro/{name}",
            "+flock_local": "True",
        }

        # Add accounting group if provided
        accounting_group = self.production.meta.get("scheduler", {}).get("accounting group", None)
        if accounting_group:
            description["accounting_group_user"] = config.get("condor", "user")
            description["accounting_group"] = accounting_group

        self.job = htcondor.Submit(description)

        with set_directory(rundir):
            with open(f"{name}.sub", "w") as subfile:
                subfile.write(self.job.__str__())

        self.logger.info(f"HTCondor job description created at {rundir}/{name}.sub")

    def submit_dag(self, dryrun=False):
        """
        Submit the p_astro calculation job to HTCondor.

        Parameters
        ----------
        dryrun : bool
            If True, don't actually submit the job

        Returns
        -------
        int
            HTCondor cluster ID
        """
        if dryrun:
            self.logger.info("Dry run mode - not submitting job")
            return -1

        self.logger.info("Submitting p_astro job to HTCondor")

        with set_directory(self.production.rundir):
            try:
                schedd = htcondor.Schedd()
            except:
                schedulers = htcondor.Collector().locate(
                    htcondor.DaemonTypes.Schedd,
                    config.get("condor", "scheduler")
                )
                schedd = htcondor.Schedd(schedulers)

            result = schedd.submit(self.job)
            cluster_id = result.cluster()

            self.logger.info(f"Submitted job {cluster_id} to HTCondor")

        self.production.job_id = int(cluster_id)
        self.production.status = "running"

        return cluster_id

    def detect_completion(self):
        """
        Check if the p_astro calculation has completed.

        Returns
        -------
        bool
            True if the job has completed successfully
        """
        rundir = self.production.rundir
        results_file = os.path.join(rundir, "pastro_results.json")

        if os.path.exists(results_file):
            self.logger.info("P_astro calculation completed - results file found")
            return True

        return False

    def collect_assets(self):
        """
        Collect output files from the p_astro calculation.

        Returns
        -------
        dict
            Dictionary of output files and results
        """
        rundir = self.production.rundir
        assets = {}

        # Main results file
        results_file = os.path.join(rundir, "pastro_results.json")
        if os.path.exists(results_file):
            assets['pastro_results'] = results_file

            # Load and log the p_astro value
            try:
                with open(results_file, 'r') as f:
                    results = json.load(f)
                    pastro = results.get('pastro', None)
                    if pastro is not None:
                        self.logger.info(f"p_astro = {pastro:.4f}")
                        self.production.meta['pastro'] = pastro
            except Exception as e:
                self.logger.error(f"Failed to read results: {e}")

        # Look for diagnostic plots
        plot_files = glob.glob(os.path.join(rundir, "*.png"))
        if plot_files:
            assets['plots'] = plot_files

        return assets

    def collect_logs(self):
        """
        Collect log files from the job.

        Returns
        -------
        dict
            Dictionary of log file contents
        """
        rundir = self.production.rundir
        name = self.production.name
        logs = {}

        for ext in ['out', 'err', 'log']:
            log_file = os.path.join(rundir, f"{name}.{ext}")
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    logs[ext] = f.read()

        return logs

    def after_completion(self):
        """
        Post-processing after successful completion.
        """
        self.logger.info("P_astro calculation completed successfully")
        self.production.status = "finished"

        # Store results in production metadata
        assets = self.collect_assets()
        if 'pastro_results' in assets:
            with open(assets['pastro_results'], 'r') as f:
                results = json.load(f)
                self.production.meta['pastro_results'] = results
