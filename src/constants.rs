use once_cell::sync::Lazy;

use crate::PositiveFloat;

pub const ZERO: Lazy<PositiveFloat> = Lazy::new(|| PositiveFloat::new(0.0).unwrap());

// Neutrinos
pub const DEFAULT_NEUTRINO_MASSES: Lazy<[PositiveFloat; 3]> = Lazy::new(|| [*ZERO, *ZERO, *ZERO]);
pub const DEFAULT_N_EFF: Lazy<PositiveFloat> = Lazy::new(|| PositiveFloat::new(3.04).unwrap()); // WMAP (ApJ Spergel et al 2007)
