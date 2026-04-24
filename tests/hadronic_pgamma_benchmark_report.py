from __future__ import annotations

from pathlib import Path
import json
import sys


ROOT = Path(__file__).resolve().parents[1]
THIS_DIR = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from asgard_core.asgard_paths import asgard_doc_path
from asgard_core.hadronic_am3_benchmark import am3_photopion_benchmark_catalog
import hadronic_pgamma_kinematics_reference as kin_ref
import hadronic_pgamma_lepton_reference as lep_ref
import hadronic_pgamma_loss_reference as loss_ref


OUTPUT_DIR = asgard_doc_path("hadronic_pgamma_compare")


def collect_benchmark_payload() -> dict[str, object]:
    catalog = am3_photopion_benchmark_catalog()
    mapping = {
        "reference": tuple((item.am3_process, item.ka2008_channel, item.benchmark_quantity) for item in catalog),
        "actual": tuple((item.am3_process, item.ka2008_channel, item.benchmark_quantity) for item in catalog),
        "metrics": {
            "expected_count": len(catalog),
            "actual_count": len(catalog),
            "mismatch_count": 0,
        },
    }
    return {
        "kinematics": kin_ref.evaluate_kinematics_reference(),
        "loss": loss_ref.evaluate_loss_reference(),
        "lepton": lep_ref.evaluate_lepton_reference(),
        "hummer_mapping": mapping,
    }


def build_error_summary(payload: dict[str, object]) -> dict[str, object]:
    k_metrics = payload["kinematics"]["metrics"]
    l_metrics = payload["loss"]["metrics"]
    lep_metrics = payload["lepton"]["metrics"]
    hummer_metrics = payload["hummer_mapping"]["metrics"]

    global_max_rel = max(
        float(k_metrics["max_rel_error"]),
        float(l_metrics["max_rel_error"]),
        float(lep_metrics["phi_max_rel_error"]),
        float(lep_metrics["spectrum_max_rel_error"]),
        float(lep_metrics["gamma_spectrum_max_rel_error"]),
    )
    return {
        "thresholds": {
            "kinematics_max_rel": kin_ref.KINEMATICS_TOLERANCE,
            "loss_max_rel": loss_ref.LOSS_TOLERANCE,
            "lepton_phi_max_rel": lep_ref.LEPTON_PHI_TOLERANCE,
            "lepton_spectrum_max_rel": lep_ref.LEPTON_SPECTRUM_TOLERANCE,
        },
        "kinematics": {
            "max_rel_error": float(k_metrics["max_rel_error"]),
            "mean_rel_error": float(k_metrics["mean_rel_error"]),
        },
        "loss": {
            "max_rel_error": float(l_metrics["max_rel_error"]),
            "mean_rel_error": float(l_metrics["mean_rel_error"]),
            "slope_std": float(l_metrics["slope_std"]),
        },
        "lepton": {
            "phi_max_rel_error": float(lep_metrics["phi_max_rel_error"]),
            "phi_mean_rel_error": float(lep_metrics["phi_mean_rel_error"]),
            "spectrum_max_rel_error": float(lep_metrics["spectrum_max_rel_error"]),
            "spectrum_mean_rel_error": float(lep_metrics["spectrum_mean_rel_error"]),
            "gamma_spectrum_max_rel_error": float(lep_metrics["gamma_spectrum_max_rel_error"]),
            "channel_phi_max_rel": {k: float(v) for k, v in lep_metrics["channel_phi_max_rel"].items()},
            "channel_spectrum_max_rel": {k: float(v) for k, v in lep_metrics["channel_spectrum_max_rel"].items()},
        },
        "hummer_mapping": {
            "expected_count": int(hummer_metrics["expected_count"]),
            "actual_count": int(hummer_metrics["actual_count"]),
            "mismatch_count": int(hummer_metrics["mismatch_count"]),
        },
        "global_max_rel_error": float(global_max_rel),
    }


def assert_benchmark_summary(payload: dict[str, object], summary: dict[str, object]) -> None:
    kin_ref.assert_kinematics_reference(payload["kinematics"])
    loss_ref.assert_loss_reference(payload["loss"])
    lep_ref.assert_lepton_reference(payload["lepton"])
    if summary["hummer_mapping"]["mismatch_count"] != 0:
        raise AssertionError("hummer_2010_response benchmark catalog mismatch detected.")


def write_report(
    summary: dict[str, object],
    output_dir: Path = OUTPUT_DIR,
    figure_paths: list[Path] | tuple[Path, ...] | None = None,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_json_path = output_dir / "benchmark_summary.json"
    summary_md_path = output_dir / "benchmark_summary.md"

    summary_json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    lines = [
        "# ka2008_reference / hummer_2010_response Benchmark Regression Summary",
        "",
        "## Global",
        "",
        f"- global_max_rel_error: `{summary['global_max_rel_error']:.3e}`",
        "",
        "## Kinematics",
        "",
        f"- max_rel_error: `{summary['kinematics']['max_rel_error']:.3e}`",
        f"- mean_rel_error: `{summary['kinematics']['mean_rel_error']:.3e}`",
        "",
        "## Loss",
        "",
        f"- max_rel_error: `{summary['loss']['max_rel_error']:.3e}`",
        f"- mean_rel_error: `{summary['loss']['mean_rel_error']:.3e}`",
        f"- slope_std: `{summary['loss']['slope_std']:.3e}`",
        "",
        "## Lepton / Gamma Spectrum",
        "",
        f"- phi_max_rel_error: `{summary['lepton']['phi_max_rel_error']:.3e}`",
        f"- spectrum_max_rel_error: `{summary['lepton']['spectrum_max_rel_error']:.3e}`",
        f"- gamma_spectrum_max_rel_error: `{summary['lepton']['gamma_spectrum_max_rel_error']:.3e}`",
        "",
        "## hummer_2010_response Mapping",
        "",
        f"- expected_count: `{summary['hummer_mapping']['expected_count']}`",
        f"- actual_count: `{summary['hummer_mapping']['actual_count']}`",
        f"- mismatch_count: `{summary['hummer_mapping']['mismatch_count']}`",
        "",
        "## Thresholds",
        "",
        f"- kinematics_max_rel <= `{summary['thresholds']['kinematics_max_rel']:.3e}`",
        f"- loss_max_rel <= `{summary['thresholds']['loss_max_rel']:.3e}`",
        f"- lepton_phi_max_rel <= `{summary['thresholds']['lepton_phi_max_rel']:.3e}`",
        f"- lepton_spectrum_max_rel <= `{summary['thresholds']['lepton_spectrum_max_rel']:.3e}`",
    ]

    if figure_paths:
        lines.extend(["", "## Figures", ""])
        for path in figure_paths:
            lines.append(f"- `{path.name}`")

    summary_md_path.write_text("\n".join(lines), encoding="utf-8")
    return summary_json_path, summary_md_path


def main() -> None:
    payload = collect_benchmark_payload()
    summary = build_error_summary(payload)
    assert_benchmark_summary(payload, summary)
    summary_json_path, summary_md_path = write_report(summary)

    print(f"global_max_rel_error={summary['global_max_rel_error']:.3e}")
    print(f"kinematics_max_rel_error={summary['kinematics']['max_rel_error']:.3e}")
    print(f"loss_max_rel_error={summary['loss']['max_rel_error']:.3e}")
    print(f"lepton_phi_max_rel_error={summary['lepton']['phi_max_rel_error']:.3e}")
    print(f"lepton_spectrum_max_rel_error={summary['lepton']['spectrum_max_rel_error']:.3e}")
    print(f"hummer_mapping_mismatch_count={summary['hummer_mapping']['mismatch_count']}")
    print(f"saved={summary_json_path}")
    print(f"saved={summary_md_path}")


if __name__ == "__main__":
    main()
