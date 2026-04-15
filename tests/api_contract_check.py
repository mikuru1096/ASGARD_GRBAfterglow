from __future__ import annotations

from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_models import ReverseShockConfig, SpectrumOutputConfig
from asgard_presets import build_baseline_config
from mergered import (
    plot_characteristic_frequencies,
    plot_spectrum,
    run_fit,
)


def build_contract_config():
    return build_baseline_config(
        num_gam_e=161,
        num_nu=161,
        num_r=240,
        num_theta=240,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=1.0e-1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
        ),
        spectrum_output=SpectrumOutputConfig(enabled=True, num_nu_obs=96),
    )


def main() -> None:
    result = run_fit(build_contract_config())

    assert result.characteristic_time_s.shape == result.nu_m.shape
    assert result.characteristic_time_s.shape == result.nu_c.shape
    assert result.characteristic_time_s.shape == result.nu_a.shape
    assert result.rs_nu_m is not None
    assert result.rs_nu_c is not None
    assert result.rs_nu_a is not None
    assert result.characteristic_time_s.shape == result.rs_nu_m.shape
    assert result.characteristic_time_s.shape == result.rs_nu_c.shape
    assert result.characteristic_time_s.shape == result.rs_nu_a.shape
    assert result.spectrum_freq_hz is not None
    assert result.spectrum_fnu is not None
    assert result.spectrum_fnu.shape == (result.spectrum_freq_hz.shape[0], result.t_obs_s.shape[0])

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        freq_plot = tmpdir_path / "characteristic.pdf"
        spec_plot = tmpdir_path / "spectrum.pdf"

        plot_characteristic_frequencies(result, include_reverse=True, outfile=str(freq_plot), show=False)
        plot_spectrum(result, times_s=[1.0e3, 1.0e4], quantity="nufnu", outfile=str(spec_plot), show=False)

        assert freq_plot.exists() and freq_plot.stat().st_size > 0
        assert spec_plot.exists() and spec_plot.stat().st_size > 0

    print("PASS: API contract check succeeded.")


if __name__ == "__main__":
    main()
