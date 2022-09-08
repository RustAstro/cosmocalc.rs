use once_cell::sync::Lazy;

use crate::{
    eV,
    units::{
        JoulePerKelvin, JoulePerMeter3Kelvin4, JouleSeconds, Meters3PerKgPerSecond2,
        MetersPerSecond, PositiveFloat, WattsPerMeters2Kelvin4,
    },
    DimensionlessPositiveFloat,
};

pub const PI: f64 = std::f64::consts::PI;

// Neutrinos
pub static DEFAULT_NEUTRINO_MASSES: Lazy<[eV; 3]> = Lazy::new(|| {
    [
        DimensionlessPositiveFloat::zero(),
        DimensionlessPositiveFloat::zero(),
        DimensionlessPositiveFloat::zero(),
    ]
});
pub static DEFAULT_N_EFF: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| DimensionlessPositiveFloat::new(3.04).unwrap()); // WMAP (ApJ Spergel et al 2007)

/// Ratio of neutrino to photon temperature
pub static T_NU_TO_T_GAMMA_RATIO: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| PositiveFloat((4. / 11_f64).powf(1. / 3.)));

// SI units: https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units

/// Speed of light
pub static C_M_PER_S: MetersPerSecond = 299792458.;

/// Gravitational constant
pub static G: Meters3PerKgPerSecond2 = 6.6743e-11;

/// Boltzmann constant
pub static BOLTZMANN: JoulePerKelvin = 1.380649e-23;

/// Stefan-Boltzmann constant
pub static STEFAN_BOLTZMANN: WattsPerMeters2Kelvin4 = 5.670374e-8;

/// Reduced Planck constant
pub static H_BAR: JouleSeconds = 1.05457e-34;

/// $\alpha = \frac{\pi^2 k^4}{15 \hbar^3 c^3}$
///
/// Eqn 2.29 from Ryden
pub static ALPHA: Lazy<JoulePerMeter3Kelvin4> =
    Lazy::new(|| PI.powf(2.) * BOLTZMANN.powf(4.) / (15. * H_BAR.powf(3.) * C_M_PER_S.powf(3.)));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alpha() {
        assert!(*ALPHA > 7.0e-16);
        assert!(*ALPHA < 8.0e-16);
    }
}
