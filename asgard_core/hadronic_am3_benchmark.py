from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AM3BenchmarkChannel:
    am3_process: str
    ka2008_channel: str
    description: str
    benchmark_quantity: str


def am3_photopion_benchmark_catalog() -> tuple[AM3BenchmarkChannel, ...]:
    """Reference mapping between AM3 process groups and KA2008 reference channels.

    This is metadata only: it freezes the comparison targets used in regression tests.
    """
    return (
        AM3BenchmarkChannel(
            am3_process="photopion",
            ka2008_channel="gamma",
            description="Neutral-pion decay photons / photopion gamma production",
            benchmark_quantity="Q_gamma_pi0",
        ),
        AM3BenchmarkChannel(
            am3_process="muon_decay",
            ka2008_channel="e_plus",
            description="Positron channel from positive-muon decay in the charged-pion chain",
            benchmark_quantity="Q_eplus",
        ),
        AM3BenchmarkChannel(
            am3_process="muon_decay",
            ka2008_channel="nu_mu_bar",
            description="Antimuon-neutrino channel from positive-muon decay",
            benchmark_quantity="Q_numu_bar",
        ),
        AM3BenchmarkChannel(
            am3_process="pion_decay",
            ka2008_channel="nu_mu",
            description="Prompt muon-neutrino channel; AM3 muon-decay contribution must be summed separately",
            benchmark_quantity="Q_numu",
        ),
        AM3BenchmarkChannel(
            am3_process="muon_decay",
            ka2008_channel="nu_e",
            description="Electron-neutrino channel from positive-muon decay",
            benchmark_quantity="Q_nue",
        ),
        AM3BenchmarkChannel(
            am3_process="muon_decay",
            ka2008_channel="e_minus",
            description="Electron channel from negative-muon decay",
            benchmark_quantity="Q_eminus",
        ),
        AM3BenchmarkChannel(
            am3_process="muon_decay",
            ka2008_channel="nu_e_bar",
            description="Electron-antineutrino channel from negative-muon decay",
            benchmark_quantity="Q_nue_bar",
        ),
    )
