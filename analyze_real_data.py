from __future__ import annotations

import argparse
from pathlib import Path

from thz_field_analysis import (
    calculate_scan_metrics,
    load_field_cube,
    analyze_field,
    create_all_visuals,
    write_markdown_report,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze real spatiotemporal THz electric-field image cubes."
    )
    parser.add_argument(
        "--ex",
        required=True,
        help="Path to Ex .npy, or a .npz file containing Ex and Ey.",
    )
    parser.add_argument(
        "--ey",
        default=None,
        help="Path to Ey .npy. Omit when --ex is a combined .npz file.",
    )
    parser.add_argument(
        "--time-axis",
        type=int,
        default=0,
        choices=[0, 1, 2],
        help="Axis containing time in the stored cube. Default: 0.",
    )
    parser.add_argument(
        "--stage-step-um",
        type=float,
        required=True,
        help="Mechanical delay-stage step in micrometres.",
    )
    parser.add_argument(
        "--delay-multiplier",
        type=float,
        default=2.0,
        help="Use 2 for a retroreflector delay line; otherwise use the correct optical-path multiplier.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.10,
        help="Fraction of peak spectrum used to define usable bandwidth. Default: 0.10.",
    )
    parser.add_argument(
        "--max-spectrum-thz",
        type=float,
        default=3.0,
        help="Maximum frequency shown in the spectrum plot.",
    )
    parser.add_argument(
        "--output",
        default="outputs/real_data",
        help="Output directory.",
    )
    parser.add_argument(
        "--skip-animation",
        action="store_true",
        help="Skip GIF creation for a faster analysis run.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    ex, ey = load_field_cube(
        ex_path=args.ex,
        ey_path=args.ey,
        time_axis=args.time_axis,
    )

    metrics = calculate_scan_metrics(
        stage_step_um=args.stage_step_um,
        n_time_samples=ex.shape[0],
        delay_multiplier=args.delay_multiplier,
    )

    analysis, arrays = analyze_field(
        ex=ex,
        ey=ey,
        dt_s=metrics.time_step_s,
        threshold_fraction=args.threshold,
    )

    output_dir = Path(args.output)
    create_all_visuals(
        ex=ex,
        ey=ey,
        metrics=metrics,
        analysis=analysis,
        arrays=arrays,
        output_dir=output_dir,
        max_spectrum_thz=args.max_spectrum_thz,
        create_animation=not args.skip_animation,
    )

    report_path = output_dir.parent / f"{output_dir.name}_REPORT.md"
    write_markdown_report(
        metrics=metrics,
        analysis=analysis,
        output_dir=output_dir,
        report_path=report_path,
        dataset_name=Path(args.ex).name,
    )

    print(f"Analysis complete. Report: {report_path}")


if __name__ == "__main__":
    main()
