from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional
import json
import math

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c
from scipy.signal import hilbert


@dataclass
class ScanMetrics:
    """Sampling quantities calculated from a mechanical delay scan."""

    stage_step_um: float
    n_time_samples: int
    delay_multiplier: float
    stage_scan_span_um: float
    time_step_s: float
    time_span_s: float
    fft_record_duration_s: float
    sampling_rate_hz: float
    frequency_resolution_hz: float
    nyquist_frequency_hz: float
    n_positive_frequency_frames: int

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class FieldAnalysis:
    """Measured quantities extracted from E_x(x,y,t) and E_y(x,y,t)."""

    reference_y: int
    reference_x: int
    dominant_component: str
    pulse_peak_time_s: float
    pulse_fwhm_s: float
    spectrum_peak_frequency_hz: float
    spectral_fwhm_low_hz: float
    spectral_fwhm_high_hz: float
    spectral_fwhm_hz: float
    threshold_fraction: float
    threshold_band_low_hz: float
    threshold_band_high_hz: float
    threshold_bandwidth_hz: float
    usable_frequency_frames: int

    def to_dict(self) -> dict:
        return asdict(self)


def calculate_scan_metrics(
    stage_step_um: float,
    n_time_samples: int,
    delay_multiplier: float = 2.0,
) -> ScanMetrics:
    """
    Calculate scan and FFT sampling quantities.

    Parameters
    ----------
    stage_step_um:
        Mechanical translation-stage step in micrometres.
    n_time_samples:
        Number of recorded time-domain images.
    delay_multiplier:
        Optical path multiplier. Use 2 for a retroreflector delay line because
        moving the stage by dL changes the optical path by 2*dL. Use 1 when the
        entered displacement already represents optical path change.

    Notes
    -----
    dt = delay_multiplier * stage_step / c
    measured time span = (N - 1) * dt
    NumPy FFT record duration = N * dt
    df = 1 / (N * dt)
    f_Nyquist = 1 / (2 * dt)
    """
    if stage_step_um <= 0:
        raise ValueError("stage_step_um must be positive.")
    if n_time_samples < 2:
        raise ValueError("n_time_samples must be at least 2.")
    if delay_multiplier <= 0:
        raise ValueError("delay_multiplier must be positive.")

    step_m = stage_step_um * 1e-6
    dt = delay_multiplier * step_m / c
    span = (n_time_samples - 1) * dt
    fft_duration = n_time_samples * dt
    fs = 1.0 / dt
    df = 1.0 / fft_duration
    f_nyquist = fs / 2.0
    stage_span_um = (n_time_samples - 1) * stage_step_um
    n_positive = n_time_samples // 2 + 1

    return ScanMetrics(
        stage_step_um=float(stage_step_um),
        n_time_samples=int(n_time_samples),
        delay_multiplier=float(delay_multiplier),
        stage_scan_span_um=float(stage_span_um),
        time_step_s=float(dt),
        time_span_s=float(span),
        fft_record_duration_s=float(fft_duration),
        sampling_rate_hz=float(fs),
        frequency_resolution_hz=float(df),
        nyquist_frequency_hz=float(f_nyquist),
        n_positive_frequency_frames=int(n_positive),
    )


def make_synthetic_vector_pulse(
    n_time_samples: int = 201,
    ny: int = 96,
    nx: int = 96,
    stage_step_um: float = 10.0,
    delay_multiplier: float = 2.0,
    pulse_center_ps: float = 6.0,
    field_sigma_ps: float = 0.32,
    beam_sigma_pixels: float = 22.0,
    wavefront_delay_ps: float = 0.25,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create a simple radially polarized single-cycle THz field movie.

    Returns
    -------
    t_s:
        Time axis, shape (Nt,).
    ex:
        Horizontal electric-field cube, shape (Nt, Ny, Nx).
    ey:
        Vertical electric-field cube, shape (Nt, Ny, Nx).

    The model contains:
    - a Gaussian transverse beam,
    - a derivative-of-Gaussian single-cycle pulse,
    - radial polarization,
    - a small radius-dependent delay to make the movie spatiotemporal.
    """
    metrics = calculate_scan_metrics(
        stage_step_um=stage_step_um,
        n_time_samples=n_time_samples,
        delay_multiplier=delay_multiplier,
    )
    t_s = np.arange(n_time_samples) * metrics.time_step_s

    y = np.arange(ny) - (ny - 1) / 2
    x = np.arange(nx) - (nx - 1) / 2
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    spatial = np.exp(-(rr**2) / (2.0 * beam_sigma_pixels**2))

    radial_x = np.divide(xx, rr, out=np.zeros_like(xx, dtype=float), where=rr > 0)
    radial_y = np.divide(yy, rr, out=np.zeros_like(yy, dtype=float), where=rr > 0)

    radius_norm = rr / max(rr.max(), 1.0)
    local_delay_s = wavefront_delay_ps * 1e-12 * radius_norm**2

    tau = t_s[:, None, None] - pulse_center_ps * 1e-12 - local_delay_s[None, :, :]
    sigma_s = field_sigma_ps * 1e-12

    # Derivative of a Gaussian: a basic single-cycle THz pulse.
    pulse = -(tau / sigma_s) * np.exp(-(tau**2) / (2.0 * sigma_s**2))
    pulse /= np.max(np.abs(pulse))

    ex = pulse * spatial[None, :, :] * radial_x[None, :, :]
    ey = pulse * spatial[None, :, :] * radial_y[None, :, :]

    # Add a weak deterministic spatial-temporal ripple for realism.
    ripple = 0.015 * np.sin(2 * np.pi * 0.42e12 * t_s)[:, None, None]
    ripple = ripple * np.exp(-(rr[None, :, :] ** 2) / (2.0 * (1.3 * beam_sigma_pixels) ** 2))
    ex = ex + ripple * np.cos(np.arctan2(yy, xx))[None, :, :]
    ey = ey + ripple * np.sin(np.arctan2(yy, xx))[None, :, :]

    return t_s, ex.astype(np.float64), ey.astype(np.float64)


def _standardize_cube(cube: np.ndarray, time_axis: int) -> np.ndarray:
    cube = np.asarray(cube, dtype=float)
    if cube.ndim != 3:
        raise ValueError(f"Expected a 3D cube, received shape {cube.shape}.")
    if time_axis not in (0, 1, 2):
        raise ValueError("time_axis must be 0, 1, or 2.")
    return np.moveaxis(cube, time_axis, 0)


def load_field_cube(
    ex_path: str | Path,
    ey_path: Optional[str | Path] = None,
    time_axis: int = 0,
    ex_key: str = "Ex",
    ey_key: str = "Ey",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load electric-field cubes from .npy or .npz files.

    Accepted forms
    --------------
    1. ex.npy and ey.npy
    2. one .npz file containing arrays named Ex and Ey

    The returned shape is always (Nt, Ny, Nx).
    """
    ex_path = Path(ex_path)
    if not ex_path.exists():
        raise FileNotFoundError(ex_path)

    if ex_path.suffix.lower() == ".npz":
        with np.load(ex_path) as data:
            if ex_key not in data:
                raise KeyError(f"{ex_key!r} not found in {ex_path.name}.")
            ex = data[ex_key]
            if ey_key in data:
                ey = data[ey_key]
            else:
                ey = np.zeros_like(ex)
    elif ex_path.suffix.lower() == ".npy":
        ex = np.load(ex_path)
        if ey_path is None:
            ey = np.zeros_like(ex)
        else:
            ey = np.load(Path(ey_path))
    else:
        raise ValueError("Use .npy or .npz input files.")

    ex = _standardize_cube(ex, time_axis)
    ey = _standardize_cube(ey, time_axis)

    if ex.shape != ey.shape:
        raise ValueError(f"Ex shape {ex.shape} and Ey shape {ey.shape} do not match.")

    return ex, ey


def _fwhm_interval(axis: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or values.ndim != 1 or len(axis) != len(values):
        raise ValueError("axis and values must be one-dimensional arrays of equal length.")
    if not np.any(np.isfinite(values)):
        return math.nan, math.nan, math.nan

    peak = np.nanmax(values)
    if not np.isfinite(peak) or peak <= 0:
        return math.nan, math.nan, math.nan

    mask = values >= 0.5 * peak
    indices = np.flatnonzero(mask)
    if len(indices) == 0:
        return math.nan, math.nan, math.nan

    low = float(axis[indices[0]])
    high = float(axis[indices[-1]])
    return low, high, high - low


def _threshold_interval(
    axis: np.ndarray,
    values: np.ndarray,
    threshold_fraction: float,
) -> tuple[float, float, float, int]:
    if not 0 < threshold_fraction < 1:
        raise ValueError("threshold_fraction must lie between 0 and 1.")
    values = np.asarray(values, dtype=float)
    peak = np.nanmax(values)
    if not np.isfinite(peak) or peak <= 0:
        return math.nan, math.nan, math.nan, 0
    mask = values >= threshold_fraction * peak
    indices = np.flatnonzero(mask)
    if len(indices) == 0:
        return math.nan, math.nan, math.nan, 0
    low = float(axis[indices[0]])
    high = float(axis[indices[-1]])
    return low, high, high - low, int(mask.sum())


def analyze_field(
    ex: np.ndarray,
    ey: np.ndarray,
    dt_s: float,
    threshold_fraction: float = 0.10,
) -> tuple[FieldAnalysis, dict]:
    """
    Analyze a vector THz electric-field movie.

    The brightest pixel is selected from total time-integrated energy. At that
    pixel, the stronger of Ex or Ey is used as the signed reference waveform.
    Spectral bandwidth is calculated from the spatially averaged spectral
    amplitude, avoiding cancellation in vector beams.
    """
    ex = np.asarray(ex, dtype=float)
    ey = np.asarray(ey, dtype=float)

    if ex.shape != ey.shape or ex.ndim != 3:
        raise ValueError("Ex and Ey must be matching 3D arrays shaped (Nt, Ny, Nx).")
    if dt_s <= 0:
        raise ValueError("dt_s must be positive.")

    nt, ny, nx = ex.shape
    t_s = np.arange(nt) * dt_s

    energy_map = np.sum(ex**2 + ey**2, axis=0)
    ref_y, ref_x = np.unravel_index(np.argmax(energy_map), energy_map.shape)
    ex_trace = ex[:, ref_y, ref_x]
    ey_trace = ey[:, ref_y, ref_x]

    if np.max(np.abs(ex_trace)) >= np.max(np.abs(ey_trace)):
        trace = ex_trace
        component = "Ex"
    else:
        trace = ey_trace
        component = "Ey"

    envelope = np.abs(hilbert(trace - np.mean(trace)))
    _, _, pulse_fwhm_s = _fwhm_interval(t_s, envelope)
    pulse_peak_time_s = float(t_s[int(np.argmax(envelope))])

    ex_f = np.fft.rfft(ex - np.mean(ex, axis=0, keepdims=True), axis=0)
    ey_f = np.fft.rfft(ey - np.mean(ey, axis=0, keepdims=True), axis=0)
    freq_hz = np.fft.rfftfreq(nt, d=dt_s)

    spectral_amplitude = np.sqrt(np.mean(np.abs(ex_f) ** 2 + np.abs(ey_f) ** 2, axis=(1, 2)))
    if spectral_amplitude.max() > 0:
        spectral_amplitude_norm = spectral_amplitude / spectral_amplitude.max()
    else:
        spectral_amplitude_norm = spectral_amplitude

    # Ignore the DC bin when locating the THz spectral peak.
    if len(freq_hz) > 1:
        peak_index = 1 + int(np.argmax(spectral_amplitude_norm[1:]))
    else:
        peak_index = 0
    peak_frequency_hz = float(freq_hz[peak_index])

    fwhm_low, fwhm_high, fwhm_bw = _fwhm_interval(freq_hz, spectral_amplitude_norm)
    th_low, th_high, th_bw, usable_frames = _threshold_interval(
        freq_hz, spectral_amplitude_norm, threshold_fraction
    )

    analysis = FieldAnalysis(
        reference_y=int(ref_y),
        reference_x=int(ref_x),
        dominant_component=component,
        pulse_peak_time_s=pulse_peak_time_s,
        pulse_fwhm_s=float(pulse_fwhm_s),
        spectrum_peak_frequency_hz=peak_frequency_hz,
        spectral_fwhm_low_hz=float(fwhm_low),
        spectral_fwhm_high_hz=float(fwhm_high),
        spectral_fwhm_hz=float(fwhm_bw),
        threshold_fraction=float(threshold_fraction),
        threshold_band_low_hz=float(th_low),
        threshold_band_high_hz=float(th_high),
        threshold_bandwidth_hz=float(th_bw),
        usable_frequency_frames=int(usable_frames),
    )

    arrays = {
        "t_s": t_s,
        "trace": trace,
        "envelope": envelope,
        "energy_map": energy_map,
        "freq_hz": freq_hz,
        "ex_f": ex_f,
        "ey_f": ey_f,
        "spectral_amplitude_norm": spectral_amplitude_norm,
        "field_magnitude": np.sqrt(ex**2 + ey**2),
    }
    return analysis, arrays


def _save_and_close(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_scan_metrics(metrics: ScanMetrics, output_path: Path) -> None:
    fig = plt.figure(figsize=(9, 6))
    fig.suptitle("THz Delay-Scan and FFT Sampling Calculator", fontsize=16)
    text = (
        "Mechanical delay scan\n"
        f"Stage step: {metrics.stage_step_um:.3f} µm\n"
        f"Number of time images: {metrics.n_time_samples}\n"
        f"Stage scan span: {metrics.stage_scan_span_um:.3f} µm\n"
        f"Optical-path multiplier: {metrics.delay_multiplier:g}\n\n"
        "Time-domain sampling\n"
        "Δt = multiplier × ΔL / c\n"
        f"Time step, Δt: {metrics.time_step_s * 1e15:.3f} fs\n"
        f"First-to-last time span: {metrics.time_span_s * 1e12:.3f} ps\n"
        f"FFT record duration, NΔt: {metrics.fft_record_duration_s * 1e12:.3f} ps\n"
        f"Sampling rate: {metrics.sampling_rate_hz / 1e12:.3f} THz\n\n"
        "Frequency-domain sampling\n"
        "Δf = 1 / (NΔt)\n"
        f"Frequency resolution, Δf: {metrics.frequency_resolution_hz / 1e12:.5f} THz\n"
        "fNyquist = 1 / (2Δt)\n"
        f"Nyquist frequency: {metrics.nyquist_frequency_hz / 1e12:.3f} THz\n"
        f"Positive-frequency frames: {metrics.n_positive_frequency_frames}"
    )
    fig.text(0.08, 0.88, text, va="top", family="monospace", fontsize=12)
    _save_and_close(fig, output_path)


def plot_temporal_trace(arrays: dict, analysis: FieldAnalysis, output_path: Path) -> None:
    t_ps = arrays["t_s"] * 1e12
    trace = arrays["trace"]
    envelope = arrays["envelope"]
    norm = max(np.max(np.abs(trace)), 1e-30)

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.plot(t_ps, trace / norm, label=f"{analysis.dominant_component} field")
    ax.plot(t_ps, envelope / max(envelope.max(), 1e-30), label="Envelope")
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Normalized electric field")
    ax.set_title(
        f"Time-Domain THz Pulse at Pixel (x={analysis.reference_x}, y={analysis.reference_y})"
    )
    ax.legend()
    _save_and_close(fig, output_path)


def plot_spectrum(arrays: dict, analysis: FieldAnalysis, output_path: Path, max_thz: Optional[float] = None) -> None:
    freq_thz = arrays["freq_hz"] / 1e12
    spectrum = arrays["spectral_amplitude_norm"]

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.plot(freq_thz, spectrum)
    ax.axhline(0.5, linestyle="--", linewidth=0.9, label="FWHM level")
    ax.axhline(analysis.threshold_fraction, linestyle=":", linewidth=0.9, label="Usable-band threshold")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Normalized spectral amplitude")
    ax.set_title("Spatially Averaged THz Spectrum")
    if max_thz is not None:
        ax.set_xlim(0, max_thz)
    ax.legend()
    _save_and_close(fig, output_path)


def plot_energy_map(arrays: dict, output_path: Path) -> None:
    energy = arrays["energy_map"]
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    image = ax.imshow(energy, origin="lower")
    fig.colorbar(image, ax=ax, label="Integrated field energy (a.u.)")
    ax.set_xlabel("x pixel")
    ax.set_ylabel("y pixel")
    ax.set_title("Time-Integrated THz Field Energy")
    _save_and_close(fig, output_path)


def plot_time_snapshot(arrays: dict, time_index: int, output_path: Path) -> None:
    field = arrays["field_magnitude"]
    t_ps = arrays["t_s"][time_index] * 1e12
    frame = field[time_index]
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    image = ax.imshow(frame, origin="lower")
    fig.colorbar(image, ax=ax, label="|E| (a.u.)")
    ax.set_xlabel("x pixel")
    ax.set_ylabel("y pixel")
    ax.set_title(f"THz Spatial Frame at t = {t_ps:.3f} ps")
    _save_and_close(fig, output_path)


def plot_signed_component_snapshot(
    ex: np.ndarray,
    ey: np.ndarray,
    arrays: dict,
    time_index: int,
    component: str,
    output_path: Path,
) -> None:
    cube = ex if component == "Ex" else ey
    frame = cube[time_index]
    vmax = max(np.max(np.abs(frame)), 1e-30)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    image = ax.imshow(frame, origin="lower", vmin=-vmax, vmax=vmax)
    fig.colorbar(image, ax=ax, label=f"{component} (a.u.)")
    ax.set_xlabel("x pixel")
    ax.set_ylabel("y pixel")
    ax.set_title(f"{component} at t = {arrays['t_s'][time_index] * 1e12:.3f} ps")
    _save_and_close(fig, output_path)


def plot_frequency_map(arrays: dict, frequency_index: int, output_path: Path) -> None:
    ex_f = arrays["ex_f"]
    ey_f = arrays["ey_f"]
    freq_thz = arrays["freq_hz"][frequency_index] / 1e12
    amplitude = np.sqrt(np.abs(ex_f[frequency_index]) ** 2 + np.abs(ey_f[frequency_index]) ** 2)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    image = ax.imshow(amplitude, origin="lower")
    fig.colorbar(image, ax=ax, label="Spectral field amplitude (a.u.)")
    ax.set_xlabel("x pixel")
    ax.set_ylabel("y pixel")
    ax.set_title(f"THz Frequency Frame at f = {freq_thz:.4f} THz")
    _save_and_close(fig, output_path)


def plot_xt_map(arrays: dict, y_index: int, output_path: Path) -> None:
    field = arrays["field_magnitude"][:, y_index, :]
    extent = [
        0,
        field.shape[1] - 1,
        arrays["t_s"][0] * 1e12,
        arrays["t_s"][-1] * 1e12,
    ]
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    image = ax.imshow(field, origin="lower", aspect="auto", extent=extent)
    fig.colorbar(image, ax=ax, label="|E| (a.u.)")
    ax.set_xlabel("x pixel")
    ax.set_ylabel("Time (ps)")
    ax.set_title(f"Spatiotemporal x–t Map at y = {y_index}")
    _save_and_close(fig, output_path)


def create_time_animation(arrays: dict, output_path: Path, max_frames: int = 80) -> None:
    field = arrays["field_magnitude"]
    nt = field.shape[0]
    indices = np.unique(np.linspace(0, nt - 1, min(nt, max_frames)).astype(int))
    global_max = max(np.max(field), 1e-30)

    images = []
    for index in indices:
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        image = ax.imshow(field[index], origin="lower", vmin=0, vmax=global_max)
        fig.colorbar(image, ax=ax, label="|E| (a.u.)")
        ax.set_xlabel("x pixel")
        ax.set_ylabel("y pixel")
        ax.set_title(f"t = {arrays['t_s'][index] * 1e12:.3f} ps")
        fig.tight_layout()
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba())
        images.append(rgba[:, :, :3].copy())
        plt.close(fig)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    imageio.mimsave(output_path, images, duration=0.10, loop=0)


def create_all_visuals(
    ex: np.ndarray,
    ey: np.ndarray,
    metrics: ScanMetrics,
    analysis: FieldAnalysis,
    arrays: dict,
    output_dir: str | Path,
    max_spectrum_thz: Optional[float] = 3.0,
    create_animation: bool = True,
) -> list[Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    path = output_dir / "01_scan_calculator.png"
    plot_scan_metrics(metrics, path)
    written.append(path)

    path = output_dir / "02_time_domain_pulse.png"
    plot_temporal_trace(arrays, analysis, path)
    written.append(path)

    path = output_dir / "03_frequency_spectrum.png"
    plot_spectrum(arrays, analysis, path, max_thz=max_spectrum_thz)
    written.append(path)

    path = output_dir / "04_integrated_energy_map.png"
    plot_energy_map(arrays, path)
    written.append(path)

    peak_index = int(np.argmax(arrays["envelope"]))
    candidate_indices = sorted(
        set(
            int(np.clip(i, 0, ex.shape[0] - 1))
            for i in (peak_index - 8, peak_index, peak_index + 8)
        )
    )
    for number, index in enumerate(candidate_indices, start=1):
        path = output_dir / f"05_time_frame_{number:02d}.png"
        plot_time_snapshot(arrays, index, path)
        written.append(path)

    path = output_dir / "06_signed_component_at_peak.png"
    plot_signed_component_snapshot(
        ex, ey, arrays, peak_index, analysis.dominant_component, path
    )
    written.append(path)

    positive_freq = arrays["freq_hz"]
    requested_thz = [0.3, 0.6, 1.0]
    for number, target_thz in enumerate(requested_thz, start=1):
        index = int(np.argmin(np.abs(positive_freq - target_thz * 1e12)))
        path = output_dir / f"07_frequency_frame_{number:02d}_{positive_freq[index]/1e12:.4f}THz.png"
        plot_frequency_map(arrays, index, path)
        written.append(path)

    path = output_dir / "08_spatiotemporal_xt_map.png"
    plot_xt_map(arrays, ex.shape[1] // 2, path)
    written.append(path)

    if create_animation:
        path = output_dir / "09_time_evolution.gif"
        create_time_animation(arrays, path)
        written.append(path)

    with open(output_dir / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(
            {"scan_metrics": metrics.to_dict(), "field_analysis": analysis.to_dict()},
            handle,
            indent=2,
        )
    written.append(output_dir / "summary.json")

    return written


def _fmt(value: float, scale: float, digits: int = 4) -> str:
    if not np.isfinite(value):
        return "not available"
    return f"{value / scale:.{digits}f}"


def write_markdown_report(
    metrics: ScanMetrics,
    analysis: FieldAnalysis,
    output_dir: str | Path,
    report_path: str | Path,
    dataset_name: str = "THz electric-field dataset",
) -> None:
    output_dir = Path(output_dir)
    report_path = Path(report_path)
    import os
    rel = Path(os.path.relpath(output_dir, report_path.parent)).as_posix()

    lines = [
        "# Spatiotemporal THz Electric-Field Report",
        "",
        f"Dataset: **{dataset_name}**",
        "",
        "## 1. What was recorded?",
        "",
        "A camera records one two-dimensional electric-field image at each delay position.",
        "Stacking the images creates a three-dimensional movie:",
        "",
        r"\[",
        r"E_x(x,y,t), \qquad E_y(x,y,t)",
        r"\]",
        "",
        r"- \(x,y\): position in the camera image",
        r"- \(t\): pump–probe delay",
        "- each time image is one **time-domain frame**",
        r"- the FFT along \(t\) creates frequency images \(E_x(x,y,f)\) and \(E_y(x,y,f)\)",
        "",
        "## 2. Delay-scan calculations",
        "",
        "| Quantity | Result |",
        "|---|---:|",
        f"| Stage step | {metrics.stage_step_um:.3f} µm |",
        f"| Number of time-domain images | {metrics.n_time_samples} |",
        f"| Stage scan span | {metrics.stage_scan_span_um:.3f} µm |",
        f"| Optical-path multiplier | {metrics.delay_multiplier:g} |",
        f"| Time step, Δt | {_fmt(metrics.time_step_s, 1e-15, 3)} fs |",
        f"| First-to-last time span | {_fmt(metrics.time_span_s, 1e-12, 4)} ps |",
        f"| FFT record duration, NΔt | {_fmt(metrics.fft_record_duration_s, 1e-12, 4)} ps |",
        f"| Frequency resolution, Δf | {_fmt(metrics.frequency_resolution_hz, 1e12, 5)} THz |",
        f"| Nyquist frequency | {_fmt(metrics.nyquist_frequency_hz, 1e12, 4)} THz |",
        f"| Positive-frequency frames | {metrics.n_positive_frequency_frames} |",
        "",
        "The equations are:",
        "",
        r"\[",
        r"\Delta t = \frac{m\,\Delta L}{c}",
        r"\]",
        "",
        r"\[",
        r"\Delta f = \frac{1}{N\Delta t}",
        r"\]",
        "",
        r"\[",
        r"f_{\mathrm{Nyquist}} = \frac{1}{2\Delta t}",
        r"\]",
        "",
        r"where \(m=2\) for a retroreflector delay line.",
        "",
        "## 3. Pulse and bandwidth results",
        "",
        "| Quantity | Result |",
        "|---|---:|",
        f"| Reference pixel | x={analysis.reference_x}, y={analysis.reference_y} |",
        f"| Dominant field component | {analysis.dominant_component} |",
        f"| Pulse peak time | {_fmt(analysis.pulse_peak_time_s, 1e-12, 4)} ps |",
        f"| Time-domain pulse width (envelope FWHM) | {_fmt(analysis.pulse_fwhm_s, 1e-12, 4)} ps |",
        f"| Spectral peak | {_fmt(analysis.spectrum_peak_frequency_hz, 1e12, 4)} THz |",
        f"| Spectral FWHM range | {_fmt(analysis.spectral_fwhm_low_hz, 1e12, 4)}–{_fmt(analysis.spectral_fwhm_high_hz, 1e12, 4)} THz |",
        f"| Spectral FWHM bandwidth | {_fmt(analysis.spectral_fwhm_hz, 1e12, 4)} THz |",
        f"| {analysis.threshold_fraction*100:.0f}% usable range | {_fmt(analysis.threshold_band_low_hz, 1e12, 4)}–{_fmt(analysis.threshold_band_high_hz, 1e12, 4)} THz |",
        f"| Frequency frames above threshold | {analysis.usable_frequency_frames} |",
        "",
        "## 4. Visual explanation",
        "",
        "### Scan calculator",
        f"![Scan calculator]({rel}/01_scan_calculator.png)",
        "",
        "### Time-domain pulse",
        f"![Time-domain pulse]({rel}/02_time_domain_pulse.png)",
        "",
        "### Frequency spectrum",
        f"![Frequency spectrum]({rel}/03_frequency_spectrum.png)",
        "",
        "### Time-integrated energy",
        f"![Energy map]({rel}/04_integrated_energy_map.png)",
        "",
        "### Three time-domain image frames",
        f"![Time frame 1]({rel}/05_time_frame_01.png)",
        "",
        f"![Time frame 2]({rel}/05_time_frame_02.png)",
        "",
        f"![Time frame 3]({rel}/05_time_frame_03.png)",
        "",
        "### Spatiotemporal x–t map",
        f"![x-t map]({rel}/08_spatiotemporal_xt_map.png)",
        "",
        "### Animated field movie",
        f"Open `{rel}/09_time_evolution.gif`.",
        "",
        "## 5. Important interpretation",
        "",
        "**Time resolution is not the same as pulse width.**",
        "",
        r"- Time step \(\Delta t\) tells how closely two measurements are separated.",
        "- Pulse width tells how long the measured THz transient lasts.",
        r"- Frequency resolution \(\Delta f\) improves when the time window becomes longer.",
        "- Nyquist frequency increases when the time step becomes smaller.",
        "- Nyquist frequency is only a mathematical sampling limit; the experimentally",
        "  useful bandwidth may be much smaller because of source, detector, optics,",
        "  noise, and signal-to-noise ratio.",
        "",
        "**Number of frequency frames**",
        "",
        r"For \(N\) real time-domain images, a real FFT produces:",
        "",
        r"\[",
        r"N_f = \left\lfloor \frac{N}{2} \right\rfloor + 1",
        r"\]",
        "",
        "positive-frequency frames, including DC and the Nyquist bin when applicable.",
        "Only a subset may contain useful THz signal.",
        "",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")


def save_dataset_npz(
    output_path: str | Path,
    t_s: np.ndarray,
    ex: np.ndarray,
    ey: np.ndarray,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, t_s=t_s, Ex=ex, Ey=ey)
