# Spatiotemporal THz Electric-Field Explorer

A beginner-friendly Python repository for understanding and calculating a
spatiotemporal terahertz electric field:

\[
E_x(x,y,t), \qquad E_y(x,y,t)
\]

The project explains, calculates, and visualizes:

- mechanical delay-stage scan length,
- time step and time resolution,
- total time window,
- sampling rate,
- frequency resolution,
- Nyquist frequency,
- number of time-domain image frames,
- number of positive-frequency image frames,
- time-domain THz pulse width,
- spectral peak and bandwidth,
- spatial THz field maps,
- frequency-resolved THz field maps,
- spatiotemporal \(x-t\) maps,
- an animated THz electric-field movie.

The default example uses **201 time-domain images**, a **10 µm stage step**, and
a **retroreflector**, closely matching a common THz imaging scan.

---

## 1. The simplest possible explanation

Imagine a camera taking a picture of a moving wave.

1. The first picture is taken at delay \(t_0\).
2. The delay stage moves a small distance.
3. The next picture is taken at \(t_1\).
4. This continues until many pictures are collected.

Each picture is a two-dimensional THz electric-field map. Putting all pictures
in order makes a movie:

```text
image 1, image 2, image 3, ... image N
                 ↓
             E(x, y, t)
```

A Fourier transform is then applied along the time direction at every camera
pixel:

```text
E(x, y, t)  --FFT along time-->  E(x, y, f)
```

Therefore, a time-domain image movie becomes a set of frequency-domain images.

---

## 2. Core equations

For a delay stage step \(\Delta L\):

\[
\Delta t = \frac{m\Delta L}{c}
\]

where:

- \(c\) is the speed of light,
- \(m=2\) for a retroreflector because the optical path changes twice the
  mechanical stage movement,
- \(m=1\) when the entered displacement already represents optical-path change.

For \(N\) recorded time images:

\[
T_{\mathrm{span}}=(N-1)\Delta t
\]

NumPy's FFT bin spacing uses the record duration \(N\Delta t\):

\[
\Delta f = \frac{1}{N\Delta t}
\]

The Nyquist frequency is:

\[
f_{\mathrm{Nyquist}}=\frac{1}{2\Delta t}
\]

For real-valued time-domain data, the number of non-negative frequency frames is:

\[
N_f=\left\lfloor\frac{N}{2}\right\rfloor+1
\]

---

## 3. Numerical example: 201 images and 10 µm stage step

For a retroreflector:

```text
stage step             = 10 µm
number of images       = 201
stage scan span        = 2000 µm = 2 mm
time step              ≈ 66.713 fs
first-to-last span     ≈ 13.343 ps
FFT frequency spacing  ≈ 0.07458 THz
Nyquist frequency      ≈ 7.495 THz
positive-frequency bins = 101
```

The Nyquist frequency is only the maximum frequency that can be represented
without aliasing. It is not automatically the useful experimental bandwidth.
The useful bandwidth is usually smaller because of the THz source, detector,
optics, absorption, and noise.

---

## 4. Time resolution, pulse width, and frequency resolution are different

### Time step or sampling interval

\(\Delta t\) is the time separation between neighboring images.

### Pulse width

The pulse width is measured from the THz waveform itself. This repository uses
the full width at half maximum of the analytic-signal envelope.

### Frequency resolution

\(\Delta f\) is the spacing between FFT frequency bins. A longer time record
gives finer frequency spacing.

### Nyquist frequency

The Nyquist frequency is controlled by the time step. A smaller time step gives
a larger Nyquist frequency.

A useful rule is:

```text
longer scan  → better frequency resolution
smaller step → higher Nyquist frequency
```

---

## 5. Repository layout

```text
spatiotemporal-thz-electric-field/
├── README.md
├── DEMO_REPORT.md
├── DATA_FORMAT.md
├── GITHUB_UPLOAD.md
├── requirements.txt
├── run_demo.py
├── analyze_real_data.py
├── thz_field_analysis/
│   ├── __init__.py
│   └── core.py
├── data/
│   └── example/
├── outputs/
│   └── demo/
└── notebooks/
    └── beginner_tutorial.ipynb
```

---

## 6. Installation

Create and activate a Python environment, then run:

```bash
pip install -r requirements.txt
```

---

## 7. Run the complete demonstration

```bash
python run_demo.py
```

The program will:

1. calculate all delay-scan quantities,
2. generate a synthetic radially polarized THz electric-field movie,
3. calculate the time-domain pulse width,
4. calculate the frequency spectrum and bandwidth,
5. generate time-domain and frequency-domain image frames,
6. create an \(x-t\) map,
7. create an animated GIF,
8. write `DEMO_REPORT.md`,
9. save the example dataset as an `.npz` file.

---

## 8. Analyze real data

For one `.npz` file containing arrays named `Ex` and `Ey`:

```bash
python analyze_real_data.py \
  --ex my_thz_field.npz \
  --stage-step-um 10 \
  --delay-multiplier 2
```

For a faster run without creating the animated GIF:

```bash
python analyze_real_data.py \
  --ex my_thz_field.npz \
  --stage-step-um 10 \
  --delay-multiplier 2 \
  --skip-animation
```

For two `.npy` files:

```bash
python analyze_real_data.py \
  --ex Ex.npy \
  --ey Ey.npy \
  --stage-step-um 10 \
  --delay-multiplier 2
```

See `DATA_FORMAT.md` for details.

---

## 9. How the spatiotemporal field is calculated

At each pixel \((x,y)\), the camera provides a time waveform:

\[
E_x(t),\quad E_y(t)
\]

The vector field magnitude is:

\[
|E(x,y,t)|=\sqrt{E_x^2(x,y,t)+E_y^2(x,y,t)}
\]

The complex frequency-domain fields are:

\[
\tilde{E}_x(x,y,f)=\mathcal{F}_t\{E_x(x,y,t)\}
\]

\[
\tilde{E}_y(x,y,f)=\mathcal{F}_t\{E_y(x,y,t)\}
\]

The spectral field magnitude is:

\[
|\tilde{E}(x,y,f)|=
\sqrt{|\tilde{E}_x|^2+|\tilde{E}_y|^2}
\]

Because the FFT is calculated along time independently at every pixel, each
frequency bin becomes a complete two-dimensional THz image.

---

## 10. How bandwidth is reported

This repository reports two bandwidths:

1. **Spectral FWHM bandwidth**: the frequency range above half of the maximum
   spectral amplitude.
2. **Threshold bandwidth**: the frequency range above a user-defined fraction
   of the maximum, default 10%.

For rigorous experiments, a signal-to-noise threshold is often better than a
simple percentage threshold. The current code is intentionally transparent and
easy to modify.

---

## 11. What a school student should remember

- A time-domain THz camera dataset is a movie.
- The stage step controls how often the movie is sampled.
- The total scan length controls how finely frequency can be separated.
- The FFT changes a time movie into frequency pictures.
- More time images produce more frequency bins.
- Not every frequency bin contains useful signal.
- The source and detector determine the real experimental bandwidth.

---

## 12. Citation suggestion

When using this educational code in research, cite the repository and describe
the exact sampling formulas, threshold definition, and pulse-width definition
used in your analysis.

---

## 13. Verify the calculations

Run the included tests:

```bash
python -m unittest
```

The test checks the 201-image, 10 µm, retroreflector example against the expected
time step, frequency resolution, Nyquist frequency, scan span, and FFT frame count.

## 14. Scientific caution about “time resolution”

This repository calculates the **time sampling interval** from the delay-stage
step. The true experimental temporal resolution can also depend on the optical
gate-pulse duration, detector response, phase matching, electronics, and other
instrument effects. Therefore, do not claim that the stage-derived time step
alone is the complete instrument response.

---

## 15. Use only the scan calculator

No THz dataset is required for this command:

```bash
python scan_calculator.py \
  --stage-step-um 10 \
  --n-images 201 \
  --delay-multiplier 2 \
  --max-thz 3
```

The same scan can be entered through its mechanical span:

```bash
python scan_calculator.py \
  --stage-step-um 10 \
  --scan-span-um 2000 \
  --delay-multiplier 2 \
  --max-thz 3
```

For this example, there are 41 non-negative FFT bins from 0 to 3 THz,
including the DC frame.

---

## Live browser simulation

This repository now includes a self-contained `index.html` web application.
It uses only HTML, CSS, and JavaScript, so no Python installation or external
JavaScript library is required for the live GitHub Pages version.

To enable GitHub Pages:

1. Push `index.html` and `.nojekyll` to the root of the `main` branch.
2. Open **Repository → Settings → Pages**.
3. Under **Build and deployment**, select **Deploy from a branch**.
4. Select **main** and **/(root)**.
5. Save and wait for GitHub Pages to redeploy.

The expected live URL is:

```text
https://kararwaldeepak.github.io/spatiotemporal-thz-electric-field/
```

The page includes interactive controls for scan step, number of frames,
retroreflector multiplier, pulse duration, bandwidth threshold, polarization,
time-domain frames, and frequency-domain frames.
