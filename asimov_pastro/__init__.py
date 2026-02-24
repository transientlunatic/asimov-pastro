"""
Asimov pipeline for computing astrophysical probability (p_astro) of
gravitational wave candidates.

Implements the Bayesian framework from:
  Ashton, Thrane & Smith (2019)  arXiv:1909.11872
  Ashton & Thrane (2020)         arXiv:2006.05039

Expected upstream dependencies (declared in the blueprint via ``needs``):

* **Per-IFO PE runs** — one bilby analysis per detector, run incoherently.
  Each must advertise a ``"samples"`` asset (bilby result HDF5 file).
  The IFOs present are detected automatically from the posterior parameters.
* **Coherent PE run** — a multi-IFO bilby analysis.
  Also advertised as ``"samples"``; distinguished by having >1 IFO.
* **pyomicron run** (optional) — advertises a ``"trigger list"`` asset
  mapping IFO codes to lists of Omicron trigger files.

Blueprint configuration (under the ``pastro`` key):

.. code-block:: yaml

    pastro:
      merger_rate: 28.3        # R, Gpc^-3 yr^-1
      sensitive_volume: 1.5    # V_T, Gpc^3
      segment_duration: 0.2    # delta-T, seconds (default 0.2)
      approximant: C01:IMRPhenomXPHM  # only needed for PESummary input files
      context_duration: 86400  # Omicron trigger window, seconds (default 24 h)
      default_xi_g: 0.01       # fallback glitch probability if no triggers available
"""

import json
import os

try:
    import asimov.pipeline
    from asimov.utils import set_directory
    from asimov import config as asimov_config
except (ImportError, ModuleNotFoundError):
    asimov = None
    asimov_config = None

from .evidence import detect_ifos_in_bilby_result, read_evidence
from .odds import compute_log_noise_evidence, compute_pastro, compute_xi, compute_xi_g_from_triggers

try:
    import htcondor2 as htcondor
except ImportError:
    try:
        import htcondor
    except ImportError:
        htcondor = None


_DEFAULT_SEGMENT_DURATION = 0.2    # seconds
_DEFAULT_CONTEXT_DURATION = 86400  # 24 hours, as used in arXiv:2006.05039
_DEFAULT_XI_G = 0.01               # conservative fallback glitch probability


_Pipeline = asimov.pipeline.Pipeline if asimov is not None else object


class Pastro(_Pipeline):
    """
    Asimov pipeline for computing p_astro (astrophysical probability).

    See module docstring for usage details.
    """

    name = "pastro"

    def __init__(self, production, category=None):
        super().__init__(production, category)
        self.logger.info(
            f"Initialising p_astro pipeline for {production.name}"
        )

    # ------------------------------------------------------------------
    # Configuration helpers
    # ------------------------------------------------------------------

    @property
    def _pastro_config(self):
        return self.production.meta.get("pastro", {})

    @property
    def merger_rate(self):
        return float(self._pastro_config["merger_rate"])

    @property
    def sensitive_volume(self):
        return float(self._pastro_config["sensitive_volume"])

    @property
    def segment_duration(self):
        return float(
            self._pastro_config.get("segment_duration", _DEFAULT_SEGMENT_DURATION)
        )

    @property
    def context_duration(self):
        return float(
            self._pastro_config.get("context_duration", _DEFAULT_CONTEXT_DURATION)
        )

    @property
    def approximant(self):
        return self._pastro_config.get("approximant", None)

    @property
    def default_xi_g(self):
        return float(
            self._pastro_config.get("default_xi_g", _DEFAULT_XI_G)
        )

    # ------------------------------------------------------------------
    # Asset retrieval from upstream dependencies
    # ------------------------------------------------------------------

    def _get_all_pe_sample_files(self):
        """Return all HDF5/JSON sample files from upstream PE dependencies."""
        assets = self.production._previous_assets()
        sample_files = assets.get("samples", [])
        if isinstance(sample_files, str):
            sample_files = [sample_files]
        return [f for f in sample_files if f.endswith((".hdf5", ".h5", ".json"))]

    def _partition_pe_files(self, sample_files):
        """
        Split PE result files into one coherent file and per-IFO files.

        Uses :func:`~asimov_pastro.evidence.detect_ifos_in_bilby_result` to
        inspect each file.  The file with the largest IFO set is taken as the
        coherent result; files with exactly one IFO are taken as single-detector
        results.

        Returns
        -------
        coherent_file : str or None
        per_ifo_files : dict mapping IFO code to file path
        """
        ifo_map = {}  # file -> list of IFOs
        for f in sample_files:
            try:
                ifos = detect_ifos_in_bilby_result(f)
            except Exception as e:
                self.logger.warning(
                    f"Could not detect IFOs in {f}: {e} — skipping"
                )
                continue
            if ifos:
                ifo_map[f] = ifos

        if not ifo_map:
            return None, {}

        # Coherent file: most IFOs
        coherent_file = max(ifo_map, key=lambda f: len(ifo_map[f]))
        coherent_ifos = ifo_map[coherent_file]

        per_ifo_files = {}
        for f, ifos in ifo_map.items():
            if f == coherent_file:
                continue
            if len(ifos) == 1:
                per_ifo_files[ifos[0]] = f
            else:
                self.logger.warning(
                    f"{f} has {ifos} IFOs but is not the coherent file — ignoring"
                )

        self.logger.info(
            f"Coherent PE file ({coherent_ifos}): {coherent_file}"
        )
        for ifo, f in per_ifo_files.items():
            self.logger.info(f"Single-IFO PE file ({ifo}): {f}")

        return coherent_file, per_ifo_files

    def _get_trigger_list(self):
        """Return the Omicron trigger list from upstream assets, or None."""
        assets = self.production._previous_assets()
        trigger_list = assets.get("trigger list", None)
        if trigger_list:
            self.logger.info(
                f"Found Omicron triggers for IFOs: {list(trigger_list.keys())}"
            )
        else:
            self.logger.info(
                "No Omicron trigger list found — will use default xi_g"
            )
        return trigger_list

    # ------------------------------------------------------------------
    # Core calculation
    # ------------------------------------------------------------------

    def run_calculation(self):
        """
        Perform the p_astro calculation and return results as a dict.

        This method is called from the HTCondor job script.  It is also
        usable directly (e.g. in tests) without an HTCondor environment.

        Returns
        -------
        dict with keys:
            ``p_astro``, ``log_odds``, ``log_Z_S``, ``log_Z_N``, ``xi``,
            ``xi_g``, ``ifos``, ``log_evidence_coherent``,
            ``log_noise_evidence_per_ifo``.
        """
        # --- locate and classify PE result files -------------------------
        sample_files = self._get_all_pe_sample_files()
        if not sample_files:
            raise RuntimeError(
                "No PE sample files found in upstream dependencies. "
                "Ensure coherent and per-IFO bilby analyses are listed "
                "in the blueprint 'needs' field and have completed."
            )

        coherent_file, per_ifo_files = self._partition_pe_files(sample_files)

        if coherent_file is None:
            raise RuntimeError(
                "Could not identify a coherent (multi-IFO) PE result file. "
                "Check that the coherent bilby analysis has completed and "
                "that its result file is being advertised."
            )
        if not per_ifo_files:
            raise RuntimeError(
                "No single-IFO PE result files found. "
                "Ensure per-IFO bilby analyses are listed in 'needs' and "
                "have completed."
            )

        ifos = sorted(per_ifo_files.keys())
        self.logger.info(f"IFOs for noise model: {ifos}")

        # --- read evidences ----------------------------------------------
        self.logger.info("Reading coherent evidence from PE result file")
        coherent_ev = read_evidence(coherent_file, label=self.approximant)
        log_Z_S = coherent_ev["log_evidence"]
        self.logger.info(f"  log Z_S = {log_Z_S:.3f}")

        log_Z_glitch = {}
        log_Z_gaussian = {}
        for ifo, f in per_ifo_files.items():
            self.logger.info(f"Reading single-IFO evidence for {ifo}")
            ev = read_evidence(f, label=self.approximant)
            log_Z_glitch[ifo] = ev["log_evidence"]
            log_Z_gaussian[ifo] = ev["log_noise_evidence"]
            self.logger.info(
                f"  {ifo}: log Z_glitch = {log_Z_glitch[ifo]:.3f}, "
                f"log Z_gaussian = {log_Z_gaussian[ifo]:.3f}"
            )

        # --- xi: prior signal probability --------------------------------
        xi = compute_xi(self.merger_rate, self.sensitive_volume, self.segment_duration)
        self.logger.info(
            f"xi = {xi:.3e}  "
            f"(R={self.merger_rate} Gpc^-3 yr^-1, "
            f"V={self.sensitive_volume} Gpc^3, "
            f"dT={self.segment_duration} s)"
        )

        # --- xi_g: per-detector glitch probability -----------------------
        trigger_list = self._get_trigger_list()
        if trigger_list:
            xi_g = compute_xi_g_from_triggers(
                trigger_list, self.context_duration, self.segment_duration
            )
            # Ensure we have xi_g for every IFO in the noise model
            for ifo in ifos:
                if ifo not in xi_g:
                    self.logger.warning(
                        f"No triggers found for {ifo} — using default xi_g"
                    )
                    xi_g[ifo] = self.default_xi_g
        else:
            xi_g = {ifo: self.default_xi_g for ifo in ifos}

        for ifo in ifos:
            self.logger.info(f"  xi_g[{ifo}] = {xi_g[ifo]:.3e}")

        # --- noise evidence (Eq. 9) --------------------------------------
        log_Z_N = compute_log_noise_evidence(
            ifos, log_Z_glitch, log_Z_gaussian, xi_g
        )
        self.logger.info(f"log Z_N = {log_Z_N:.3f}")

        # --- p_astro -----------------------------------------------------
        p_astro, log_odds = compute_pastro(log_Z_S, log_Z_N, xi)
        self.logger.info(f"log odds = {log_odds:.3f}")
        self.logger.info(f"p_astro  = {p_astro:.6f}")

        return {
            "p_astro": p_astro,
            "log_odds": log_odds,
            "log_Z_S": log_Z_S,
            "log_Z_N": log_Z_N,
            "xi": xi,
            "xi_g": xi_g,
            "ifos": ifos,
            "log_evidence_per_ifo": log_Z_glitch,
            "log_noise_evidence_per_ifo": log_Z_gaussian,
            "inputs": {
                "merger_rate": self.merger_rate,
                "sensitive_volume": self.sensitive_volume,
                "segment_duration": self.segment_duration,
                "context_duration": self.context_duration,
                "coherent_file": coherent_file,
                "per_ifo_files": per_ifo_files,
            },
        }

    # ------------------------------------------------------------------
    # HTCondor job management
    # ------------------------------------------------------------------

    def build_dag(self, dryrun=False):
        """
        Write the calculation script and HTCondor submit file.
        """
        self.logger.info("Building p_astro calculation job")
        name = self.production.name
        rundir = self.production.rundir
        os.makedirs(rundir, exist_ok=True)

        # Write a small driver script that imports and calls run_calculation
        script_path = os.path.join(rundir, f"{name}_pastro.py")
        script = f"""\
#!/usr/bin/env python
# Auto-generated by asimov-pastro — do not edit by hand.
import json, logging, sys
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

from asimov.ledger import Ledger
from asimov.event import Event

ledger = Ledger("{self.production.event.repository.directory}")
event = Event.from_ledger(ledger, "{self.production.event.name}")
production = event.get_production("{self.production.name}")
pipeline = production.pipeline

results = pipeline.run_calculation()

with open("pastro_results.json", "w") as fh:
    json.dump(results, fh, indent=2)

print(f"p_astro = {{results['p_astro']:.6f}}")
sys.exit(0)
"""
        with open(script_path, "w") as fh:
            fh.write(script)
        os.chmod(script_path, 0o755)

        scheduler_meta = self.production.meta.get("scheduler", {})
        description = {
            "executable": script_path,
            "output": os.path.join(rundir, f"{name}.out"),
            "error": os.path.join(rundir, f"{name}.err"),
            "log": os.path.join(rundir, f"{name}.log"),
            "request_disk": scheduler_meta.get("request disk", "1GB"),
            "request_memory": scheduler_meta.get("request memory", "2GB"),
            "request_cpus": str(scheduler_meta.get("cpus", 1)),
            "batch_name": f"pastro/{name}",
            "+flock_local": "True",
        }

        accounting_group = scheduler_meta.get("accounting group")
        if accounting_group:
            description["accounting_group_user"] = asimov_config.get("condor", "user")
            description["accounting_group"] = accounting_group

        self.job = htcondor.Submit(description)

        with set_directory(rundir):
            with open(f"{name}.sub", "w") as fh:
                fh.write(str(self.job))

        self.logger.info(f"HTCondor submit file written to {rundir}/{name}.sub")

    def submit_dag(self, dryrun=False):
        """Submit the p_astro job to HTCondor."""
        if dryrun:
            self.logger.info("Dry run — not submitting")
            return -1

        with set_directory(self.production.rundir):
            try:
                schedd = htcondor.Schedd()
            except Exception:
                schedulers = htcondor.Collector().locate(
                    htcondor.DaemonTypes.Schedd,
                    asimov_config.get("condor", "scheduler"),
                )
                schedd = htcondor.Schedd(schedulers)

            result = schedd.submit(self.job)
            cluster_id = result.cluster()

        self.production.job_id = int(cluster_id)
        self.production.status = "running"
        self.logger.info(f"Submitted HTCondor job {cluster_id}")
        return cluster_id

    def detect_completion(self):
        """Return True when ``pastro_results.json`` exists in the run directory."""
        results_file = os.path.join(self.production.rundir, "pastro_results.json")
        if os.path.exists(results_file):
            self.logger.info("p_astro results file found — job complete")
            return True
        return False

    def collect_assets(self):
        """Return a dict of output assets, including the results JSON."""
        rundir = self.production.rundir
        assets = {}

        results_file = os.path.join(rundir, "pastro_results.json")
        if os.path.exists(results_file):
            assets["pastro results"] = results_file
            try:
                with open(results_file) as fh:
                    results = json.load(fh)
                p_astro = results.get("p_astro")
                if p_astro is not None:
                    self.logger.info(f"p_astro = {p_astro:.6f}")
                    self.production.meta["p_astro"] = p_astro
            except Exception as e:
                self.logger.error(f"Failed to read results file: {e}")

        return assets

    def collect_logs(self):
        """Return contents of HTCondor log files."""
        rundir = self.production.rundir
        name = self.production.name
        logs = {}
        for ext in ("out", "err", "log"):
            path = os.path.join(rundir, f"{name}.{ext}")
            if os.path.exists(path):
                with open(path) as fh:
                    logs[ext] = fh.read()
        return logs

    def after_completion(self):
        """Post-processing hook: store results in production metadata."""
        self.production.status = "finished"
        assets = self.collect_assets()
        if "pastro results" in assets:
            with open(assets["pastro results"]) as fh:
                self.production.meta["pastro_results"] = json.load(fh)
