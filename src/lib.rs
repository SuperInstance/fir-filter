//! # FIR Filter Design
//!
//! Finite impulse response filter design using windowed sinc method.
//! Supports Hamming, Hann, and Blackman windows with lowpass, highpass,
//! and bandpass filter configurations.

pub mod sinc;
pub mod window;
pub mod lowpass;
pub mod highpass;
pub mod bandpass;

/// Re-export of the main filter design function
pub use sinc::{sinc, normalized_sinc};
pub use window::{WindowType, apply_window};
pub use lowpass::design_lowpass;
pub use highpass::design_highpass;
pub use bandpass::design_bandpass;

/// Apply a FIR filter to a signal using direct convolution
pub fn filter(coeffs: &[f64], signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let m = coeffs.len();
    let mut output = vec![0.0; n];
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..m {
            if i >= j {
                sum += coeffs[j] * signal[i - j];
            }
        }
        output[i] = sum;
    }
    output
}

/// Compute the frequency response of a filter at a given frequency
pub fn frequency_response(coeffs: &[f64], omega: f64) -> (f64, f64) {
    let mut re = 0.0;
    let mut im = 0.0;
    for (n, &c) in coeffs.iter().enumerate() {
        re += c * (omega * n as f64).cos();
        im -= c * (omega * n as f64).sin();
    }
    (re, im)
}

/// Compute the magnitude response at a given normalized frequency (0 to π)
pub fn magnitude_response(coeffs: &[f64], omega: f64) -> f64 {
    let (re, im) = frequency_response(coeffs, omega);
    re.hypot(im)
}
