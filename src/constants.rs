use once_cell::sync::Lazy;

use crate::{
    eV,
    units::{Meters3PerKgPerSecond2, MetersPerSecond, PositiveFloat},
    DimensionlessPositiveFloat,
};

pub const PI: f64 = std::f64::consts::PI;

pub static ZERO: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| DimensionlessPositiveFloat::new(0.0).unwrap());
pub static ONE: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| DimensionlessPositiveFloat::new(1.0).unwrap());

// Neutrinos
pub static DEFAULT_NEUTRINO_MASSES: Lazy<[eV; 3]> = Lazy::new(|| [*ZERO, *ZERO, *ZERO]);
pub static DEFAULT_N_EFF: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| DimensionlessPositiveFloat::new(3.04).unwrap()); // WMAP (ApJ Spergel et al 2007)

/// Ratio of neutrino to photon temperature
pub static T_NU_TO_T_GAMMA_RATIO: Lazy<DimensionlessPositiveFloat> =
    Lazy::new(|| PositiveFloat((4. / 11_f64).powf(1. / 3.)));

// SI units: https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units
pub static C_M_PER_S: MetersPerSecond = 299792458.;

/// Gravitational constant
pub static G: Meters3PerKgPerSecond2 = 6.6743e-11;
