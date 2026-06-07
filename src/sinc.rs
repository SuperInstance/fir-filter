//! Sinc function and normalized sinc function for FIR filter design.

/// Compute the sinc function: sin(π·x) / (π·x)
/// Returns 1.0 when x is zero (handling the removable singularity).
pub fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-15 {
        1.0
    } else {
        let px = std::f64::consts::PI * x;
        px.sin() / px
    }
}

/// Compute the normalized sinc function: sin(x) / x
/// Returns 1.0 when x is zero.
pub fn normalized_sinc(x: f64) -> f64 {
    if x.abs() < 1e-15 {
        1.0
    } else {
        x.sin() / x
    }
}

/// Generate a symmetric time vector centered at zero for filter design.
/// `order` must be even; produces `order + 1` coefficients.
pub fn time_vector(order: usize) -> Vec<f64> {
    let half = order as f64 / 2.0;
    (0..=order).map(|n| n as f64 - half).collect()
}

/// Compute the ideal lowpass impulse response (sinc in time domain).
/// `cutoff` is the normalized cutoff frequency (0.0 to 1.0, where 1.0 = Nyquist).
pub fn ideal_lowpass_impulse(order: usize, cutoff: f64) -> Vec<f64> {
    let t = time_vector(order);
    t.iter().map(|&n| 2.0 * cutoff * sinc(2.0 * cutoff * n)).collect()
}

/// Compute the ideal highpass impulse response.
/// Derived from delta minus lowpass: h[n] = δ(n) - h_lp[n]
pub fn ideal_highpass_impulse(order: usize, cutoff: f64) -> Vec<f64> {
    let lp = ideal_lowpass_impulse(order, cutoff);
    let mid = order / 2;
    lp.iter().enumerate().map(|(i, &v)| {
        if i == mid { 1.0 - v } else { -v }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_sinc_at_zero() {
        assert!((sinc(0.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_sinc_at_integer() {
        // sinc(n) = 0 for non-zero integer n
        assert!(sinc(1.0).abs() < 1e-15);
        assert!(sinc(2.0).abs() < 1e-15);
        assert!(sinc(-1.0).abs() < 1e-15);
    }

    #[test]
    fn test_sinc_symmetry() {
        // sinc(x) = sinc(-x)
        for x in [0.1, 0.5, 1.5, 3.7] {
            assert!((sinc(x) - sinc(-x)).abs() < 1e-14);
        }
    }

    #[test]
    fn test_normalized_sinc_at_zero() {
        assert!((normalized_sinc(0.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_normalized_sinc_pi() {
        assert!(normalized_sinc(PI).abs() < 1e-15);
    }

    #[test]
    fn test_time_vector_centering() {
        let t = time_vector(10);
        assert!((t[5]).abs() < 1e-15); // center element should be 0
        assert_eq!(t.len(), 11);
    }

    #[test]
    fn test_time_vector_symmetry() {
        let t = time_vector(8);
        for i in 0..t.len() / 2 {
            assert!((t[i] + t[t.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_ideal_lowpass_impulse_symmetry() {
        let h = ideal_lowpass_impulse(20, 0.25);
        for i in 0..h.len() / 2 {
            assert!((h[i] - h[h.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_ideal_highpass_impulse_symmetry() {
        let h = ideal_highpass_impulse(20, 0.25);
        for i in 0..h.len() / 2 {
            assert!((h[i] - h[h.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_ideal_lowpass_peak_at_center() {
        let h = ideal_lowpass_impulse(20, 0.5);
        let mid = h.len() / 2;
        // Center value should be 2*cutoff = 1.0
        assert!((h[mid] - 1.0).abs() < 1e-14);
    }
}
