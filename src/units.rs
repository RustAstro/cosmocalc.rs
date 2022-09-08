use anyhow::anyhow;

// Continuous positive quantities that are dimensionless (e.g. ratios like the omegas)
pub type DimensionlessPositiveFloat = PositiveFloat;

// Time
pub type Gyr = PositiveFloat;
pub type Seconds = PositiveFloat;

// Energy
#[allow(non_camel_case_types)]
pub type eV = PositiveFloat;
pub type Joule = PositiveFloat;

// Mass
pub type Kilogram = PositiveFloat;
pub type Gram = PositiveFloat;

// Temperatures
pub type Kelvin = PositiveFloat;

// Hubble parameter units
pub type KmPerSecPerMpc = f64;
pub type HInvKmPerSecPerMpc = f64;

// Distances
pub type Meter = PositiveFloat;
pub type Km = PositiveFloat;
pub type Mpc = PositiveFloat;

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

// Conversions
pub const KILOMETER_TO_METER: f64 = 1000.;
pub const MPC_TO_METERS: f64 = 3.086e+22;
pub const MPC_TO_KILOMETERS: f64 = 3.086e+19;

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
