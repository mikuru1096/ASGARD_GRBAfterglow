from asgard_presets import build_spectrum_demo_config
from mergered import plot_light_curve, plot_spectrum, run_fit


def build_demo_config():
    return build_spectrum_demo_config()


def main() -> None:
    result = run_fit(build_demo_config())
    plot_light_curve(result, outfile="Radiation_Lightcurves.pdf", show=False)
    plot_spectrum(result, times_s=[1e3, 1e4, 1e5], quantity="nufnu", outfile="Radiation_Spectra.pdf", show=False)


if __name__ == "__main__":
    main()
