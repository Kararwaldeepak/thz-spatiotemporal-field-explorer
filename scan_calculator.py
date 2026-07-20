from __future__ import annotations

import argparse
import json
import math

from thz_field_analysis import calculate_scan_metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate THz delay-scan sampling and FFT quantities."
    )
    parser.add_argument("--stage-step-um", type=float, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--n-images", type=int)
    group.add_argument(
        "--scan-span-um",
        type=float,
        help="Mechanical first-to-last stage span in micrometres.",
    )
    parser.add_argument(
        "--delay-multiplier",
        type=float,
        default=2.0,
        help="Use 2 for a retroreflector.",
    )
    parser.add_argument(
        "--max-thz",
        type=float,
        default=None,
        help="Optionally count non-negative FFT frames from 0 to this frequency.",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        default=None,
        help="Optional path for saving the calculated quantities as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.n_images is not None:
        n_images = args.n_images
    else:
        ratio = args.scan_span_um / args.stage_step_um
        nearest = round(ratio)
        if not math.isclose(ratio, nearest, rel_tol=0, abs_tol=1e-8):
            raise ValueError(
                "scan-span-um must be an integer multiple of stage-step-um."
            )
        n_images = nearest + 1

    metrics = calculate_scan_metrics(
        stage_step_um=args.stage_step_um,
        n_time_samples=n_images,
        delay_multiplier=args.delay_multiplier,
    )

    print("\nTHz delay-scan calculator")
    print("-------------------------")
    print(f"Stage step:                    {metrics.stage_step_um:.6f} µm")
    print(f"Number of time images:         {metrics.n_time_samples}")
    print(f"Mechanical scan span:          {metrics.stage_scan_span_um:.6f} µm")
    print(f"Time sampling interval:        {metrics.time_step_s * 1e15:.6f} fs")
    print(f"First-to-last time span:       {metrics.time_span_s * 1e12:.6f} ps")
    print(f"FFT record duration:           {metrics.fft_record_duration_s * 1e12:.6f} ps")
    print(f"Frequency-bin spacing:         {metrics.frequency_resolution_hz / 1e12:.8f} THz")
    print(f"Nyquist frequency:             {metrics.nyquist_frequency_hz / 1e12:.6f} THz")
    print(f"Non-negative frequency frames: {metrics.n_positive_frequency_frames}")

    output = metrics.to_dict()

    if args.max_thz is not None:
        if args.max_thz < 0:
            raise ValueError("--max-thz must be non-negative.")
        max_hz = min(args.max_thz * 1e12, metrics.nyquist_frequency_hz)
        frames = int(math.floor(max_hz / metrics.frequency_resolution_hz + 1e-12)) + 1
        frames = min(frames, metrics.n_positive_frequency_frames)
        print(f"Frames from 0 to {args.max_thz:g} THz:      {frames}")
        output["requested_max_frequency_hz"] = args.max_thz * 1e12
        output["frames_from_zero_to_requested_max"] = frames

    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as handle:
            json.dump(output, handle, indent=2)
        print(f"Saved JSON: {args.json_path}")


if __name__ == "__main__":
    main()
