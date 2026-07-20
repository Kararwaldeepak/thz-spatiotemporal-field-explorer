# Real-data format

The analysis expects two electric-field movies:

- `Ex`: horizontal THz electric-field component
- `Ey`: vertical THz electric-field component

The recommended array shape is:

```text
(Nt, Ny, Nx)
```

where:

- `Nt` is the number of delay positions or time-domain images,
- `Ny` is the image height,
- `Nx` is the image width.

## Option 1: One compressed NPZ file

```python
import numpy as np

np.savez_compressed(
    "my_thz_field.npz",
    Ex=Ex,
    Ey=Ey,
)
```

Run:

```bash
python analyze_real_data.py \
  --ex my_thz_field.npz \
  --stage-step-um 10 \
  --delay-multiplier 2
```

## Option 2: Two NPY files

```python
np.save("Ex.npy", Ex)
np.save("Ey.npy", Ey)
```

Run:

```bash
python analyze_real_data.py \
  --ex Ex.npy \
  --ey Ey.npy \
  --stage-step-um 10 \
  --delay-multiplier 2
```

## Different axis order

For an array stored as `(Ny, Nx, Nt)`, time is axis 2:

```bash
python analyze_real_data.py \
  --ex my_thz_field.npz \
  --time-axis 2 \
  --stage-step-um 10
```

## Important physical requirement

The values must be signed electric-field samples, not ordinary intensity images.
A signed THz field can be positive or negative. If the camera data are only
intensity images, phase-sensitive electric-field reconstruction is not possible
without additional measurements or calibration.
