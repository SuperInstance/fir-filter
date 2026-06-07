//! Highpass FIR filter design using spectral inversion.

use crate::sinc::ideal_highpass_impulse;
use crate::window::{WindowType, generate_window, apply_window};

/// Design a highpass FIR filter.
///
/// # Arguments
/// * `order` - Filter order (must be even), produces `order + 1` coefficients
/// * `cutoff` - Normalized cutoff frequency (0.0 to 1.0, where 1.0 = Nyquist)
/// * `window_type` - Window function to apply
pub fn design_highpass(order: usize, cutoff: f64, window_type: WindowType) -> Vec<f64> {
    assert!(order.is_multiple_of(2), "Filter order must be even");
    assert!(cutoff > 0.0 && cutoff < 1.0, "Cutoff must be in (0, 1)");
    let ideal = ideal_highpass_impulse(order, cutoff);
    let window = generate_window(window_type, ideal.len());
    apply_window(&ideal, &window)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::magnitude_response;
    use std::f64::consts::PI;

    #[test]
    fn test_highpass_nyquist_gain() {
        let coeffs = design_highpass(50, 0.25, WindowType::Hamming);
        let mag = magnitude_response(&coeffs, PI * 0.9);
        assert!(mag > 0.3, "Highpass should pass high frequencies, got {}", mag);
    }

    #[test]
    fn test_highpass_symmetry() {
        let coeffs = design_highpass(50, 0.25, WindowType::Hamming);
        for i in 0..coeffs.len() / 2 {
            assert!((coeffs[i] - coeffs[coeffs.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_highpass_dc_zero() {
        // DC response: sum of all coefficients (should be near zero for highpass)
        let coeffs = design_highpass(50, 0.25, WindowType::Blackman);
        let dc: f64 = coeffs.iter().sum();
        assert!(dc.abs() < 0.05, "Highpass DC should be near zero, got {}", dc);
    }

    #[test]
    fn test_highpass_length() {
        let coeffs = design_highpass(30, 0.5, WindowType::Hann);
        assert_eq!(coeffs.len(), 31);
    }

    #[test]
    #[should_panic(expected = "Filter order must be even")]
    fn test_highpass_odd_order_panics() {
        design_highpass(31, 0.5, WindowType::Hamming);
    }
}
