//! Bandpass FIR filter design using windowed sinc method.
//! A bandpass filter is constructed by convolving a highpass and lowpass filter,
//! or equivalently, subtracting a lowpass from a higher-cutoff lowpass.

use crate::sinc;
use crate::window::{WindowType, generate_window, apply_window};

/// Design a bandpass FIR filter.
///
/// # Arguments
/// * `order` - Filter order (must be even), produces `order + 1` coefficients
/// * `low_cutoff` - Normalized low cutoff frequency (0.0 to 1.0)
/// * `high_cutoff` - Normalized high cutoff frequency (0.0 to 1.0), must be > low_cutoff
/// * `window_type` - Window function to apply
pub fn design_bandpass(
    order: usize,
    low_cutoff: f64,
    high_cutoff: f64,
    window_type: WindowType,
) -> Vec<f64> {
    assert!(order.is_multiple_of(2), "Filter order must be even");
    assert!(low_cutoff > 0.0 && high_cutoff < 1.0, "Cutoffs must be in (0, 1)");
    assert!(high_cutoff > low_cutoff, "High cutoff must exceed low cutoff");

    // Bandpass = ideal_lp(high) - ideal_lp(low)
    let ideal_low = sinc::ideal_lowpass_impulse(order, low_cutoff);
    let ideal_high = sinc::ideal_lowpass_impulse(order, high_cutoff);
    let ideal_bp: Vec<f64> = ideal_high.iter().zip(ideal_low.iter()).map(|(&h, &l)| h - l).collect();
    let window = generate_window(window_type, ideal_bp.len());
    apply_window(&ideal_bp, &window)
}

/// Design a bandstop (notch) FIR filter.
///
/// # Arguments
/// * `order` - Filter order (must be even)
/// * `low_cutoff` - Normalized low cutoff frequency
/// * `high_cutoff` - Normalized high cutoff frequency
/// * `window_type` - Window function to apply
pub fn design_bandstop(
    order: usize,
    low_cutoff: f64,
    high_cutoff: f64,
    window_type: WindowType,
) -> Vec<f64> {
    assert!(order.is_multiple_of(2), "Filter order must be even");
    assert!(low_cutoff > 0.0 && high_cutoff < 1.0, "Cutoffs must be in (0, 1)");
    assert!(high_cutoff > low_cutoff, "High cutoff must exceed low cutoff");

    // Bandstop = allpass - bandpass
    let bp = design_bandpass(order, low_cutoff, high_cutoff, window_type);
    let mid = order / 2;
    bp.iter().enumerate().map(|(i, &v)| if i == mid { 1.0 - v } else { -v }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::magnitude_response;
    use std::f64::consts::PI;

    #[test]
    fn test_bandpass_symmetry() {
        let coeffs = design_bandpass(50, 0.2, 0.4, WindowType::Hamming);
        for i in 0..coeffs.len() / 2 {
            assert!((coeffs[i] - coeffs[coeffs.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_bandpass_length() {
        let coeffs = design_bandpass(30, 0.2, 0.4, WindowType::Hann);
        assert_eq!(coeffs.len(), 31);
    }

    #[test]
    fn test_bandpass_passes_midband() {
        let coeffs = design_bandpass(50, 0.2, 0.4, WindowType::Hamming);
        // The bandpass should have higher response in passband than stopband
        let passband = magnitude_response(&coeffs, PI * 0.3);
        let stopband_low = magnitude_response(&coeffs, PI * 0.05);
        assert!(passband > stopband_low, "Passband {} should exceed stopband {}", passband, stopband_low);
    }

    #[test]
    fn test_bandpass_attenuates_low() {
        let coeffs = design_bandpass(50, 0.2, 0.4, WindowType::Blackman);
        let passband = magnitude_response(&coeffs, PI * 0.3);
        let stopband = magnitude_response(&coeffs, PI * 0.05);
        let ratio = stopband / passband;
        assert!(ratio < 0.5, "Stopband should be attenuated relative to passband, ratio={}", ratio);
    }

    #[test]
    fn test_bandpass_attenuates_high() {
        let coeffs = design_bandpass(100, 0.1, 0.4, WindowType::Blackman);
        let passband = magnitude_response(&coeffs, PI * 0.25);
        let stopband = magnitude_response(&coeffs, PI * 0.8);
        let ratio = stopband / passband;
        assert!(ratio < 0.55, "High stopband should be attenuated relative to passband, ratio={}", ratio);
    }

    #[test]
    #[should_panic(expected = "High cutoff must exceed")]
    fn test_bandpass_reversed_cutoffs_panics() {
        design_bandpass(20, 0.4, 0.2, WindowType::Hamming);
    }

    #[test]
    #[should_panic(expected = "Filter order must be even")]
    fn test_bandpass_odd_order_panics() {
        design_bandpass(21, 0.2, 0.4, WindowType::Hamming);
    }

    #[test]
    fn test_bandstop_dc_near_unity() {
        let coeffs = design_bandstop(50, 0.2, 0.4, WindowType::Hamming);
        let dc: f64 = coeffs.iter().sum();
        // Bandstop should pass DC ≈ 1.0
        assert!((dc - 1.0).abs() < 0.2, "Bandstop DC near unity, got {}", dc);
    }

    #[test]
    fn test_bandpass_dc_near_zero() {
        let coeffs = design_bandpass(50, 0.2, 0.4, WindowType::Blackman);
        let dc: f64 = coeffs.iter().sum();
        assert!(dc.abs() < 0.1, "Bandpass DC should be near zero, got {}", dc);
    }
}
