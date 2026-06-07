//! Window functions for FIR filter design.
//! Supports Hamming, Hann, Blackman, and rectangular windows.

/// Supported window function types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowType {
    /// Rectangular window (no windowing)
    Rectangular,
    /// Hamming window
    Hamming,
    /// Hann (Hanning) window
    Hann,
    /// Blackman window
    Blackman,
}

/// Generate a window function of the given type and length.
pub fn generate_window(window_type: WindowType, length: usize) -> Vec<f64> {
    match window_type {
        WindowType::Rectangular => vec![1.0; length],
        WindowType::Hamming => (0..length)
            .map(|n| 0.54 - 0.46 * (2.0 * std::f64::consts::PI * n as f64 / (length - 1) as f64).cos())
            .collect(),
        WindowType::Hann => (0..length)
            .map(|n| 0.5 * (1.0 - (2.0 * std::f64::consts::PI * n as f64 / (length - 1) as f64).cos()))
            .collect(),
        WindowType::Blackman => (0..length)
            .map(|n| {
                let x = 2.0 * std::f64::consts::PI * n as f64 / (length - 1) as f64;
                0.42 - 0.5 * x.cos() + 0.08 * (2.0 * x).cos()
            })
            .collect(),
    }
}

/// Apply a window function to a set of coefficients (element-wise multiplication).
pub fn apply_window(coeffs: &[f64], window: &[f64]) -> Vec<f64> {
    coeffs.iter().zip(window.iter()).map(|(&c, &w)| c * w).collect()
}

/// Design a windowed filter: applies the specified window to ideal impulse response.
pub fn windowed_filter(ideal: &[f64], window_type: WindowType) -> Vec<f64> {
    let window = generate_window(window_type, ideal.len());
    apply_window(ideal, &window)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rectangular_window() {
        let w = generate_window(WindowType::Rectangular, 10);
        assert_eq!(w, vec![1.0; 10]);
    }

    #[test]
    fn test_hann_endpoint_zero() {
        let w = generate_window(WindowType::Hann, 101);
        assert!(w[0].abs() < 1e-14);
        assert!(w[100].abs() < 1e-14);
    }

    #[test]
    fn test_hann_peak() {
        let w = generate_window(WindowType::Hann, 101);
        let mid = 50;
        assert!((w[mid] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_hamming_symmetry() {
        let w = generate_window(WindowType::Hamming, 51);
        for i in 0..w.len() / 2 {
            assert!((w[i] - w[w.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_blackman_symmetry() {
        let w = generate_window(WindowType::Blackman, 51);
        for i in 0..w.len() / 2 {
            assert!((w[i] - w[w.len() - 1 - i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_blackman_low_sidelobe() {
        // Blackman window values should be mostly non-negative
        // (edges can be very slightly negative due to floating point)
        let w = generate_window(WindowType::Blackman, 51);
        let neg_count = w.iter().filter(|&&v| v < -1e-10).count();
        assert!(neg_count == 0, "Blackman should not have significant negative values");
    }

    #[test]
    fn test_hann_values_range() {
        let w = generate_window(WindowType::Hann, 51);
        for &v in &w {
            assert!(v >= 0.0 && v <= 1.0);
        }
    }

    #[test]
    fn test_apply_window_identity() {
        let coeffs = vec![1.0, 2.0, 3.0, 4.0];
        let window = vec![1.0, 1.0, 1.0, 1.0];
        let result = apply_window(&coeffs, &window);
        assert_eq!(result, coeffs);
    }

    #[test]
    fn test_apply_window_zero() {
        let coeffs = vec![1.0, 2.0, 3.0];
        let window = vec![0.0, 0.0, 0.0];
        let result = apply_window(&coeffs, &window);
        assert!(result.iter().all(|&v| v.abs() < 1e-15));
    }

    #[test]
    fn test_windowed_filter_length() {
        let ideal = vec![1.0; 21];
        let result = windowed_filter(&ideal, WindowType::Hamming);
        assert_eq!(result.len(), 21);
    }
}
