# fir-filter

Finite impulse response (FIR) filter design in pure Rust.

## Features

- **Windowed sinc design** for lowpass, highpass, bandpass, and bandstop filters
- **Window functions**: Hamming, Hann, Blackman, Rectangular
- **Frequency response** computation and magnitude analysis
- **Pure `std`** — no external dependencies

## Modules

| Module | Description |
|--------|-------------|
| `sinc` | Sinc function, ideal impulse responses |
| `window` | Window function generation and application |
| `lowpass` | Lowpass FIR filter design |
| `highpass` | Highpass FIR filter design |
| `bandpass` | Bandpass and bandstop FIR filter design |

## Quick Start

```rust
use fir_filter::{design_lowpass, filter, WindowType};

// Design a 50th-order lowpass filter with 0.25 Nyquist cutoff
let coeffs = design_lowpass(50, 0.25, WindowType::Hamming);

// Apply to a signal
let signal = vec![0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0];
let filtered = filter(&coeffs, &signal);
```

## License

MIT OR Apache-2.0
