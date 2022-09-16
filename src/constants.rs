/// This module contains constants.
///
/// Fundamental constants use the [NIST] [CODATA].
///
/// [NIST]: <https://physics.nist.gov/cgi-bin/cuu/Value?u|search_for=abbr_in!>
/// [CODATA 2018]: <https://www.nist.gov/publications/codata-recommended-values-fundamental-physical-constants-2018>
use once_cell::sync::Lazy;

use crate::{
    eV,
    units::{
        JoulePerKelvin, JoulePerMeter3Kelvin4, JouleSeconds, Meters3PerKgPerSecond2,
        MetersPerSecond, PositiveFloat, WattsPerMeters2Kelvin4,
    },
    DimensionlessPositiveFloat, FloatingPointUnit,
};

/// PI
pub const PI: f64 = std::f64::consts::PI;

/// Speed of light [CODATA 2018]
pub static C_M_PER_S: MetersPerSecond = 299792458.;

/// Gravitational constant [CODATA 2018]
pub static G: Meters3PerKgPerSecond2 = 6.67430e-11;

/// Boltzmann constant [CODATA 2018]
pub static BOLTZMANN: JoulePerKelvin = 1.380649e-23;

/// Stefan-Boltzmann constant [CODATA 2018]
pub static STEFAN_BOLTZMANN: WattsPerMeters2Kelvin4 = 5.6703744194e-8;

/// Reduced Planck constant [CODATA 2018]
pub static H_BAR: JouleSeconds = 1.054571817e-34;

/// Vector of neutrino masses (defaults to 3 massless neutrinos)
pub static DEFAULT_NEUTRINO_MASSES: Lazy<[eV; 3]> =
    Lazy::new(|| [eV::zero(), eV::zero(), eV::zero()]);

/// Effective number of neutrinos
pub static DEFAULT_N_EFF: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| DimensionlessPositiveFloat::new(3.04).unwrap()); // WMAP (ApJ Spergel et al 2007)

/// Ratio of neutrino to photon temperature
pub static T_NU_TO_T_GAMMA_RATIO: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| PositiveFloat((4. / 11_f64).powf(1. / 3.)));

/// $\alpha = \frac{\pi^2 k^4}{15 \hbar^3 c^3}$
///
/// Eqn 2.29 from Ryden
pub static ALPHA: Lazy<JoulePerMeter3Kelvin4> =
    Lazy::new(|| PI.powi(2) * BOLTZMANN.powi(4) / (15. * H_BAR.powi(3) * C_M_PER_S.powi(3)));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alpha() {
        assert!(*ALPHA > 7.0e-16);
        assert!(*ALPHA < 8.0e-16);
    }
}
