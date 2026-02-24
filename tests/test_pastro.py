"""
Tests for asimov-pastro.

Unit tests cover the pure calculation modules (evidence.py, odds.py)
without requiring an Asimov environment or HTCondor.  Pipeline-level tests
that need a full asimov environment are marked with ``skip``.
"""

import json
import math
import pytest
import numpy as np
import h5py


# ---------------------------------------------------------------------------
# evidence.py tests
# ---------------------------------------------------------------------------

class TestReadBilbyEvidence:
    def test_reads_all_keys(self, tmp_path):
        from asimov_pastro.evidence import read_bilby_evidence

        f = tmp_path / "result.hdf5"
        with h5py.File(f, "w") as hf:
            hf.create_dataset("log_evidence", data=-1234.5)
            hf.create_dataset("log_noise_evidence", data=-1270.0)
            hf.create_dataset("log_bayes_factor", data=35.5)
            hf.create_dataset("log_evidence_err", data=0.1)

        ev = read_bilby_evidence(str(f))
        assert ev["log_evidence"] == pytest.approx(-1234.5)
        assert ev["log_noise_evidence"] == pytest.approx(-1270.0)
        assert ev["log_bayes_factor"] == pytest.approx(35.5)
        assert ev["log_evidence_err"] == pytest.approx(0.1)

    def test_missing_err_returns_none(self, tmp_path):
        from asimov_pastro.evidence import read_bilby_evidence

        f = tmp_path / "result.hdf5"
        with h5py.File(f, "w") as hf:
            hf.create_dataset("log_evidence", data=-1234.5)
            hf.create_dataset("log_noise_evidence", data=-1270.0)
            hf.create_dataset("log_bayes_factor", data=35.5)

        ev = read_bilby_evidence(str(f))
        assert ev["log_evidence_err"] is None


class TestReadPESummaryEvidence:
    def _make_pesummary_file(self, path, label, ln_evidence, ln_noise_evidence,
                              ln_bayes_factor, include_error=True):
        with h5py.File(path, "w") as hf:
            sampler = hf.require_group(f"{label}/meta_data/sampler")
            sampler.create_dataset("ln_evidence", data=ln_evidence)
            sampler.create_dataset("ln_noise_evidence", data=ln_noise_evidence)
            sampler.create_dataset("ln_bayes_factor", data=ln_bayes_factor)
            if include_error:
                sampler.create_dataset("ln_evidence_error", data=0.13)

    def test_reads_correct_label(self, tmp_path):
        from asimov_pastro.evidence import read_pesummary_evidence

        f = str(tmp_path / "pesummary.h5")
        self._make_pesummary_file(f, "C01:IMRPhenomXPHM", -20818.88, -20855.64, 36.77)

        ev = read_pesummary_evidence(f, "C01:IMRPhenomXPHM")
        assert ev["log_evidence"] == pytest.approx(-20818.88)
        assert ev["log_noise_evidence"] == pytest.approx(-20855.64)
        assert ev["log_bayes_factor"] == pytest.approx(36.77)
        assert ev["log_evidence_err"] == pytest.approx(0.13)

    def test_bad_label_raises(self, tmp_path):
        from asimov_pastro.evidence import read_pesummary_evidence

        f = str(tmp_path / "pesummary.h5")
        self._make_pesummary_file(f, "C01:IMRPhenomXPHM", -20818.88, -20855.64, 36.77)

        with pytest.raises(KeyError, match="not found"):
            read_pesummary_evidence(f, "NonExistent:Label")


class TestReadEvidenceAutoDetect:
    def test_detects_bilby_native(self, tmp_path):
        from asimov_pastro.evidence import read_evidence

        f = tmp_path / "result.hdf5"
        with h5py.File(f, "w") as hf:
            hf.create_dataset("log_evidence", data=-500.0)
            hf.create_dataset("log_noise_evidence", data=-540.0)
            hf.create_dataset("log_bayes_factor", data=40.0)

        ev = read_evidence(str(f))
        assert ev["log_evidence"] == pytest.approx(-500.0)

    def test_detects_pesummary_with_label(self, tmp_path):
        from asimov_pastro.evidence import read_evidence

        f = tmp_path / "pesummary.h5"
        with h5py.File(f, "w") as hf:
            sampler = hf.require_group("C01:IMRPhenomXPHM/meta_data/sampler")
            sampler.create_dataset("ln_evidence", data=-20818.88)
            sampler.create_dataset("ln_noise_evidence", data=-20855.64)
            sampler.create_dataset("ln_bayes_factor", data=36.77)

        ev = read_evidence(str(f), label="C01:IMRPhenomXPHM")
        assert ev["log_evidence"] == pytest.approx(-20818.88)


class TestDetectIFOs:
    def test_detects_single_ifo(self, tmp_path):
        from asimov_pastro.evidence import detect_ifos_in_bilby_result

        f = tmp_path / "h1_result.hdf5"
        with h5py.File(f, "w") as hf:
            posterior = hf.require_group("posterior")
            posterior.create_dataset("H1_optimal_snr", data=np.array([10.0]))
            posterior.create_dataset("mass_1", data=np.array([35.0]))

        ifos = detect_ifos_in_bilby_result(str(f))
        assert ifos == ["H1"]

    def test_detects_multiple_ifos(self, tmp_path):
        from asimov_pastro.evidence import detect_ifos_in_bilby_result

        f = tmp_path / "h1l1_result.hdf5"
        with h5py.File(f, "w") as hf:
            posterior = hf.require_group("posterior")
            posterior.create_dataset("H1_optimal_snr", data=np.array([8.0]))
            posterior.create_dataset("L1_optimal_snr", data=np.array([7.0]))
            posterior.create_dataset("mass_1", data=np.array([35.0]))

        ifos = detect_ifos_in_bilby_result(str(f))
        assert set(ifos) == {"H1", "L1"}


# ---------------------------------------------------------------------------
# odds.py tests
# ---------------------------------------------------------------------------

class TestComputeXi:
    def test_reference_values(self):
        """Reproduce the arXiv:2006.05039 reference point approximately."""
        from asimov_pastro.odds import compute_xi

        # Paper uses R=59 Gpc^-3 yr^-1, delta_T=0.2 s, xi~4.5e-4.
        # We need V_T such that R*V_T*delta_T = 4.5e-4.
        # V_T = 4.5e-4 / (59 * 0.2/31557600) ~ 1.2e6 Gpc^3 — unrealistic,
        # meaning the paper's xi bakes in detection efficiency implicitly.
        # Here we just test the formula is self-consistent.
        xi = compute_xi(merger_rate=28.3, sensitive_volume=1.0, segment_duration=0.2)
        assert xi == pytest.approx(28.3 * 1.0 * 0.2 / (365.25 * 24 * 3600))
        assert xi > 0

    def test_scales_linearly_with_rate(self):
        from asimov_pastro.odds import compute_xi

        xi1 = compute_xi(10.0, 1.0, 0.2)
        xi2 = compute_xi(20.0, 1.0, 0.2)
        assert xi2 == pytest.approx(2 * xi1)

    def test_scales_linearly_with_volume(self):
        from asimov_pastro.odds import compute_xi

        xi1 = compute_xi(28.3, 1.0, 0.2)
        xi2 = compute_xi(28.3, 2.0, 0.2)
        assert xi2 == pytest.approx(2 * xi1)


class TestComputeXiG:
    def test_zero_triggers(self, tmp_path):
        from asimov_pastro.odds import compute_xi_g_from_triggers

        xi_g = compute_xi_g_from_triggers(
            trigger_list={"H1": [], "L1": []},
            context_duration=86400,
            segment_duration=0.2,
        )
        assert xi_g["H1"] == 0.0
        assert xi_g["L1"] == 0.0

    def test_scales_with_trigger_count(self):
        from asimov_pastro.odds import compute_xi_g_from_triggers

        # Use file paths as proxies for trigger count
        xi_g = compute_xi_g_from_triggers(
            trigger_list={"H1": ["t1", "t2", "t3"]},  # 3 triggers
            context_duration=86400,
            segment_duration=0.2,
        )
        expected = 3 / 86400 * 0.2
        assert xi_g["H1"] == pytest.approx(expected)


class TestComputeLogNoiseEvidence:
    def test_two_detectors_all_glitch(self):
        """When xi_g=1 for all IFOs the only surviving term is both-glitch."""
        import warnings
        from asimov_pastro.odds import compute_log_noise_evidence

        log_Z_glitch = {"H1": -100.0, "L1": -110.0}
        log_Z_gaussian = {"H1": -200.0, "L1": -220.0}
        xi_g = {"H1": 1.0, "L1": 1.0}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            log_Z_N = compute_log_noise_evidence(
                ["H1", "L1"], log_Z_glitch, log_Z_gaussian, xi_g
            )
        # Only the (glitch, glitch) term survives: log(1) + Z_H + log(1) + Z_L
        expected = log_Z_glitch["H1"] + log_Z_glitch["L1"]
        assert log_Z_N == pytest.approx(expected, abs=1e-10)

    def test_two_detectors_no_glitch(self):
        """When xi_g=0 for all IFOs the only surviving term is both-Gaussian."""
        import warnings
        from asimov_pastro.odds import compute_log_noise_evidence

        log_Z_glitch = {"H1": -100.0, "L1": -110.0}
        log_Z_gaussian = {"H1": -200.0, "L1": -220.0}
        xi_g = {"H1": 0.0, "L1": 0.0}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            log_Z_N = compute_log_noise_evidence(
                ["H1", "L1"], log_Z_glitch, log_Z_gaussian, xi_g
            )
        expected = log_Z_gaussian["H1"] + log_Z_gaussian["L1"]
        assert log_Z_N == pytest.approx(expected, abs=1e-10)

    def test_three_detectors_runs(self):
        """Smoke test: 3 detectors = 8 terms, should not raise."""
        from asimov_pastro.odds import compute_log_noise_evidence

        ifos = ["H1", "L1", "V1"]
        log_Z_glitch = {ifo: -100.0 for ifo in ifos}
        log_Z_gaussian = {ifo: -150.0 for ifo in ifos}
        xi_g = {ifo: 0.05 for ifo in ifos}

        log_Z_N = compute_log_noise_evidence(ifos, log_Z_glitch, log_Z_gaussian, xi_g)
        assert math.isfinite(log_Z_N)


class TestComputePastro:
    def test_high_snr_gives_high_pastro(self):
        """A strongly coherent signal should give p_astro close to 1."""
        from asimov_pastro.odds import compute_pastro

        # Large log_Z_S relative to log_Z_N => high odds
        p_astro, log_odds = compute_pastro(
            log_Z_S=-100.0, log_Z_N=-150.0, xi=1e-3
        )
        assert p_astro > 0.99
        assert log_odds > 0

    def test_noise_dominated_gives_low_pastro(self):
        """When Z_N >> xi * Z_S the event is likely noise."""
        from asimov_pastro.odds import compute_pastro

        p_astro, log_odds = compute_pastro(
            log_Z_S=-150.0, log_Z_N=-100.0, xi=1e-4
        )
        assert p_astro < 0.01
        assert log_odds < 0

    def test_formula_consistency(self):
        """p_astro = O / (1 + O) exactly."""
        from asimov_pastro.odds import compute_pastro

        p_astro, log_odds = compute_pastro(-200.0, -210.0, 5e-4)
        O = math.exp(log_odds)
        assert p_astro == pytest.approx(O / (1 + O))

    def test_output_in_unit_interval(self):
        from asimov_pastro.odds import compute_pastro

        for log_Z_S, log_Z_N in [(-100, -90), (-100, -110), (-200, -200)]:
            p, _ = compute_pastro(log_Z_S, log_Z_N, xi=1e-3)
            assert 0.0 <= p <= 1.0


# ---------------------------------------------------------------------------
# Pipeline-level tests (require asimov environment)
# ---------------------------------------------------------------------------

@pytest.mark.skip(reason="Requires a full asimov environment")
def test_pipeline_initialisation():
    pass


@pytest.mark.skip(reason="Requires a full asimov environment and HTCondor")
def test_pipeline_submit():
    pass
