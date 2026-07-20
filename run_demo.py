from pathlib import Path
import json

from thz_field_analysis import (
    calculate_scan_metrics,
    make_synthetic_vector_pulse,
    analyze_field,
    create_all_visuals,
    write_markdown_report,
)
from thz_field_analysis.core import save_dataset_npz


def main() -> None:
    # Default example:
    # 201 images, 10 µm stage step, retroreflector multiplier = 2.
    # This gives approximately 66.7 fs sampling and 0.0746 THz FFT spacing.
    stage_step_um = 10.0
    n_time_samples = 201
    delay_multiplier = 2.0

    metrics = calculate_scan_metrics(
        stage_step_um=stage_step_um,
        n_time_samples=n_time_samples,
        delay_multiplier=delay_multiplier,
    )

    t_s, ex, ey = make_synthetic_vector_pulse(
        n_time_samples=n_time_samples,
        ny=96,
        nx=96,
        stage_step_um=stage_step_um,
        delay_multiplier=delay_multiplier,
    )

    analysis, arrays = analyze_field(
        ex=ex,
        ey=ey,
        dt_s=metrics.time_step_s,
        threshold_fraction=0.10,
    )

    output_dir = Path("outputs/demo")
    create_all_visuals(
        ex=ex,
        ey=ey,
        metrics=metrics,
        analysis=analysis,
        arrays=arrays,
        output_dir=output_dir,
        max_spectrum_thz=3.0,
    )

    save_dataset_npz("data/example/synthetic_radial_thz.npz", t_s, ex, ey)

    write_markdown_report(
        metrics=metrics,
        analysis=analysis,
        output_dir=output_dir,
        report_path="DEMO_REPORT.md",
        dataset_name="Synthetic radially polarized single-cycle THz pulse",
    )

    print("\nDemo completed.")
    print(f"Time step: {metrics.time_step_s * 1e15:.3f} fs")
    print(f"Frequency resolution: {metrics.frequency_resolution_hz / 1e12:.5f} THz")
    print(f"Nyquist frequency: {metrics.nyquist_frequency_hz / 1e12:.3f} THz")
    print(f"Positive-frequency frames: {metrics.n_positive_frequency_frames}")
    print(f"Pulse FWHM: {analysis.pulse_fwhm_s * 1e12:.3f} ps")
    print(f"Spectral FWHM: {analysis.spectral_fwhm_hz / 1e12:.3f} THz")
    print("Open DEMO_REPORT.md and outputs/demo/ to see the results.")


if __name__ == "__main__":
    main()
