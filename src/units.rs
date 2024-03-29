use anyhow::anyhow;

pub mod dimensionless;
pub mod energy;
pub mod length;
pub(crate) mod macros;
pub mod mass;
pub mod temperature;
pub mod time;
pub mod traits;

pub use dimensionless::DimensionlessFloat;
pub use traits::FloatingPointUnit;

// Continuous positive quantities that are dimensionless (e.g. ratios like the omegas)
pub type DimensionlessPositiveFloat = PositiveFloat;

// Hubble parameter units
pub type KmPerSecPerMpc = f64;
pub type HInvKmPerSecPerMpc = f64;

// Densities
pub type KilogramsPerMeter3 = PositiveFloat;

// Constant units
// TODO: Work out a better way to handle types for composite unit information
pub type MetersPerSecond = f64;
pub type Meters3PerKgPerSecond2 = f64;
pub type Meters2KgPerSecond2Kelvin = f64;
pub type JouleSeconds = f64;
pub type JoulePerMeter3Kelvin4 = f64;
pub type WattsPerMeters2Kelvin4 = f64;
pub type JoulePerKelvin = f64;
pub type Mpc3 = f64;
pub type HInvMpc = f64;

/// Represents continuous physical quantities that _cannot_ be negative.
#[derive(Clone, Copy, Debug, PartialOrd, PartialEq)]
pub struct PositiveFloat(pub f64);

impl PositiveFloat {
    pub fn new(x: f64) -> Result<Self, anyhow::Error> {
        if x < 0. {
            return Err(anyhow!("expected positive number"));
        }
        Ok(Self(x))
    }

    pub fn zero() -> Self {
        Self(0.)
    }

    pub fn one() -> Self {
        Self(1.)
    }

    // Passthrough methods for convenience
    pub fn floor(&self) -> f64 {
        self.0.floor()
    }

    pub fn powf(&self, exp: f64) -> f64 {
        self.0.powf(exp)
    }
}

impl std::ops::Sub for PositiveFloat {
    type Output = PositiveFloat;

    fn sub(self, rhs: Self) -> Self::Output {
        PositiveFloat(self.0 - rhs.0)
    }
}

impl std::ops::Add for PositiveFloat {
    type Output = PositiveFloat;

    fn add(self, rhs: Self) -> Self::Output {
        PositiveFloat(self.0 + rhs.0)
    }
}
