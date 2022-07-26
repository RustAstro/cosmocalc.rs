use anyhow::anyhow;

pub type DimensionlessPositiveFloat = PositiveFloat;

// Energy
#[allow(non_camel_case_types)]
pub type eV = PositiveFloat;

// Temperatures
pub type TemperatureKelvin = PositiveFloat;

// Hubble parameter units
pub type KmPerSecPerMpc = f32;
pub type HInvKmPerSecPerMpc = f32;

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
