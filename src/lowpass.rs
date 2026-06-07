//! Lowpass FIR filter design using windowed sinc method.

use crate::sinc::ideal_lowpass_impulse;
use crate::window::{WindowType, generate_window, apply_window};

/// Design a lowpass FIR filter.
///
/// # Arguments
/// * `order` - Filter order (must be even), produces `order + 1` coefficients
/// * `cutoff` - Normalized cutoff frequency (0.0 to 1.0, where 1.0 = Nyquist/fs/2)
/// * `window_type` - Window function to apply
///
/// # Returns
/// Vector of FIR filter coefficients
pub fn design_lowpass(order: usize, cutoff: f64, window_type: WindowType) -> Vec<f64> {
    assert!(order.is_multiple_of(2), "Filter order must be even");
    assert!(cutoff > 0.0 && cutoff < 1.0, "Cutoff must be in (0, 1)");
    let ideal = ideal_lowpass_impulse(order, cutoff);
    let window = generate_window(window_type, ideal.len());
    apply_window(&ideal, &window)
}

/// Compute the -3dB cutoff frequency of a designed lowpass filter.
/// Searches for the frequency where magnitude drops to 1/sqrt(2) ≈ 0.7071.
pub fn compute_3db_frequency(coeffs: &[f64]) -> f64 {
    let steps = 1000;
    let mut omega = 0.0;
    let target = 1.0 / std::f64::consts::SQRT_2;
    for i in 1..steps {
        let w = std::f64::consts::PI * i as f64 / steps as f64;
        let mag = crate::magnitude_response(coeffs, w);
        if mag < target && omega > 0.0 {
            // Linear interpolation
            return w;
        }
        if mag >= target {
            omega = w;
        }
    }
    omega
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::magnitude_response;
    use std::f64::consts::PI;

    #[test]
    fn test_lowpass_dc_gain() {
        // DC gain should be approximately 1.0 (unified by sinc integral)
        let coeffs = design_lowpass(50, 0.25, WindowType::Hamming);
        let dc = coeffs.iter().sum::<f64>();
        assert!((dc - 1.0).abs() < 0.1, "DC gain {} expected near 1.0", dc);
    }

    #[test]
    fn test_lowpass_stopband_attenuation() {
        let coeffs = design_lowpass(50, 0.25, WindowType::Hamming);
        // At high frequency (near π), magnitude should be very small
        let mag = magnitude_response(&coeffs, PI * 0.9);
        assert!(mag < 0.05, "Stopband magnitude too high: {}", mag);
    }

    #[test]
    fn test_lowpass_passband_unity() {
        let coeffs = design_lowpass(50, 0.25, WindowType::Hamming);
        // Normalize to unit DC gain
        let dc: f64 = coeffs.iter().sum();
        let norm_coeffs: Vec<f64> = coeffs.iter().map(|&c| c / dc).collect();
        let passband = magnitude_response(&norm_coeffs, PI * 0.1);
        assert!((passband - 1.0).abs() < 0.1, "Passband not near unity: {}", passband);
    }

    #[test]
    fn test_lowpass_symmetry() {
        let coeffs = design_lowpass(50, 0.25, WindowType::Blackman);
        for i in 0..coeffs.len() / 2 {
            assert!((coeffs[i] - coeffs[coeffs.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_lowpass_length() {
        let coeffs = design_lowpass(30, 0.5, WindowType::Hann);
        assert_eq!(coeffs.len(), 31);
    }

    #[test]
    #[should_panic(expected = "Filter order must be even")]
    fn test_lowpass_odd_order_panics() {
        design_lowpass(31, 0.5, WindowType::Hamming);
    }

    #[test]
    #[should_panic(expected = "Cutoff must be in")]
    fn test_lowpass_invalid_cutoff_panics() {
        design_lowpass(20, 1.5, WindowType::Hamming);
    }

    #[test]
    fn test_lowpass_different_windows() {
        let cutoff = 0.3;
        for wt in [WindowType::Hamming, WindowType::Hann, WindowType::Blackman] {
            let coeffs = design_lowpass(40, cutoff, wt);
            assert_eq!(coeffs.len(), 41);
        }
    }
}
