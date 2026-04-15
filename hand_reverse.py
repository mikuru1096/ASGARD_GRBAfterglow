from asgard_presets import build_reverse_demo_config
from mergered import plot_characteristic_frequencies, plot_light_curve, run_fit


def build_demo_config():
    return build_reverse_demo_config()


def main() -> None:
    result = run_fit(build_demo_config())
    plot_light_curve(result, outfile="Radiation_Lightcurves_RS.pdf", show=False)
    plot_characteristic_frequencies(result, include_reverse=True, outfile="Characteristic_Frequencies_RS.pdf", show=False)


if __name__ == "__main__":
    main()
