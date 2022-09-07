use anyhow::anyhow;

// Continuous positive quantities that are dimensionless (e.g. ratios like the omegas)
pub type DimensionlessPositiveFloat = PositiveFloat;

// Time
pub type Gyr = PositiveFloat;
pub type Seconds = PositiveFloat;

// Energy
#[allow(non_camel_case_types)]
pub type eV = PositiveFloat;

// Temperatures
pub type Kelvin = PositiveFloat;

// Hubble parameter units
pub type KmPerSecPerMpc = f32;
pub type HInvKmPerSecPerMpc = f32;

// Distances
pub type Meter = PositiveFloat;
pub type Km = PositiveFloat;
pub type Mpc = PositiveFloat;

// Conversions
pub const KILOMETER_TO_METER: f32 = 1000.;
pub const MPC_TO_METERS: f32 = 3.086e+22;
pub const MPC_TO_KILOMETERS: f32 = 3.086e+19;

/// Represents continuous physical quantities that _cannot_ be negative.
#[derive(Clone, Copy, Debug, PartialOrd, PartialEq)]
pub struct PositiveFloat(pub f32);

impl PositiveFloat {
    pub fn new(x: f32) -> Result<Self, anyhow::Error> {
        if x < 0. {
            return Err(anyhow!("expected positive number"));
        }
        Ok(Self(x))
    }

    // Passthrough methods for convenience
    pub fn floor(&self) -> f32 {
        self.0.floor()
    }

    pub fn powf(&self, exp: f32) -> f32 {
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
