# Sweeping Turbulence Generator

Creates evolving, multi-scale turbulent fields by continuously injecting energy at different spatial scales while allowing old structures to decay.

## The Algorithm

1. **Initialize** an empty real-space field and start at wavenumber `k = k_min`

2. **Each timestep**:
   - **Decay**: Multiply existing field by `exp(-dt/t_decay)` — old structures fade exponentially
   - **Inject**: Generate random-phase Fourier modes in a narrow k-band `[k, k+dk]`
   - **Transform**: Inverse FFT to get the new real-space contribution
   - **Accumulate**: Add `|F|²` (power) to the decayed field
   - **Sweep**: Update k according to a pattern (linear bounce or sinusoidal)

3. **Output**: Save each frame as binary with the current driven wavenumber

## Pseudocode

```
INITIALIZE:
    F_real[nx, ny] = 0          # Accumulated real-space field
    k_current = k_min           # Starting wavenumber
    decay_factor = exp(-dt / t_decay)

FOR step = 1 TO n_steps:
    
    # 1. DECAY existing field
    F_real *= decay_factor
    
    # 2. INJECT new k-space modes
    F_kspace[nx, ny] = 0
    FOR each (kx, ky) where k_current <= |k| < k_current + dk:
        phase = random(0, 2π)
        F_kspace[kx, ky] = amplitude * exp(i * phase)
    
    # 3. TRANSFORM to real space
    F_new = IFFT(F_kspace)
    
    # 4. ACCUMULATE power
    F_real += |F_new|²
    
    # 5. SWEEP wavenumber
    k_current = update_k(step)  # bounce, sine, etc.
    
    # 6. OUTPUT
    SAVE(F_real, k_current, step)
```

## Flowchart

```
┌─────────────────────────────────────┐
│           INITIALIZE                │
│  F_real = 0, k = k_min              │
│  decay = exp(-dt/t_decay)           │
└─────────────────┬───────────────────┘
                  │
                  ▼
        ┌─────────────────┐
        │  step = 1..N    │◄──────────────────┐
        └────────┬────────┘                   │
                 │                            │
                 ▼                            │
┌─────────────────────────────────────┐       │
│  DECAY: F_real *= decay_factor      │       │
└─────────────────┬───────────────────┘       │
                  │                           │
                  ▼                           │
┌─────────────────────────────────────┐       │
│  INJECT: Fill F_kspace              │       │
│  for |k| ∈ [k_current, k_current+dk]│       │
│  with random phases                 │       │
└─────────────────┬───────────────────┘       │
                  │                           │
                  ▼                           │
┌─────────────────────────────────────┐       │
│  TRANSFORM: F_new = IFFT(F_kspace)  │       │
└─────────────────┬───────────────────┘       │
                  │                           │
                  ▼                           │
┌─────────────────────────────────────┐       │
│  ACCUMULATE: F_real += |F_new|²     │       │
└─────────────────┬───────────────────┘       │
                  │                           │
                  ▼                           │
┌─────────────────────────────────────┐       │
│  SWEEP: k_current = f(step)         │       │
│  (sine, bounce, etc.)               │       │
└─────────────────┬───────────────────┘       │
                  │                           │
                  ▼                           │
┌─────────────────────────────────────┐       │
│  SAVE frame                         │       │
└─────────────────┬───────────────────┘       │
                  │                           │
                  ▼                           │
              ┌───────┐    NO                 │
              │ done? ├──────────────────────►┘
              └───┬───┘
                  │ YES
                  ▼
               [ END ]
```

## Physical Intuition

- **Low k** = large-scale structures (big blobs)
- **High k** = small-scale structures (fine grain)
- **Decay** = memory — older injections persist but fade
- **Sweeping** = the "driven" scale changes over time

The result is a field with overlapping scales: fresh structures at the current k, plus ghostly remnants of previously-driven scales. Like a driven turbulent system with finite correlation time.

## Parameters

| Parameter | Meaning |
|-----------|------------------------------------------|
| `k_min`   | Minimum wavenumber                       |
| `k_max`   | Maximum wavenumber                       |
| `dk`      | Width of the active k-band               |
| `t_decay` | Exponential decay timescale              |
| `dt`      | Timestep                                 |
| `n_steps` | Number of frames                         |
| `amplitude` | Injection strength                     |

## Building

### With CMake

```bash
cmake -B build
cmake --build build
```

### Manual

```bash
mpif90 -O2 -I/usr/include sweep.f90 -o sweep -lfftw3_mpi -lfftw3 -lm
```

## Running

```bash
# Generate frames (k_min k_max dk t_decay dt n_steps amplitude)
mpirun -np 1 ./build/sweep 4 64 8 5.0 1.0 250 1.0

# Visualize
pip install -r requirements.txt
python visualize_sweep.py
```

## Development

Use the included devcontainer for a complete development environment with:
- GFortran
- OpenMPI
- Parallel FFTW3
- Parallel HDF5
- Python with visualization dependencies

## Output Format

Each frame is saved as `frame_NNNN.bin` with:
- Header: `nx, ny` (int32)
- Metadata: `k_current` (float64)
- Data: `|F|²` field (float64, ny × nx)

Visualization produces:
- `sweep_animation.gif` - animated GIF of all frames
- `sweep_first.png`, `sweep_middle.png`, `sweep_last.png` - static frames
