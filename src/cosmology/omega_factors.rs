use crate::{units::PositiveFloat, DimensionlessPositiveFloat};

/// Represents a collection of dimensionless density parameters.
pub struct OmegaFactors {
    /// Ratio of non-relativistic matter to critical density at `z=0`.
    pub Omega_M0: DimensionlessPositiveFloat,
    /// Ratio of dark energy density to critical density at `z=0`.
    pub Omega_DE0: DimensionlessPositiveFloat,
    /// Ratio of baryon density to critical density at `z=0`.
    pub Omega_b0: DimensionlessPositiveFloat,
}

impl OmegaFactors {
    pub fn new(Omega_M0: f64, Omega_DE0: f64, Omega_b0: f64) -> Result<Self, anyhow::Error> {
        if Omega_b0 > Omega_M0 {
            return Err(anyhow::anyhow!("cannot have more baryons than matter"));
        }

        Ok(OmegaFactors {
            Omega_M0: DimensionlessPositiveFloat::new(Omega_M0)?,
            Omega_DE0: DimensionlessPositiveFloat::new(Omega_DE0)?,
            Omega_b0: DimensionlessPositiveFloat::new(Omega_b0)?,
        })
    }

    /// Dark matter density at `z=0` is matter at `z=0` minus baryons at `z=0`.
    pub fn omega_dark_matter_density_0(&self) -> DimensionlessPositiveFloat {
        self.Omega_M0 - self.Omega_b0
    }

    /// Curvature density at `z=0` given densities of relativistic particles at `z=0`.
    pub fn curvature_density_0(
        &self,
        omega_nu0: DimensionlessPositiveFloat,
        omega_gamma0: DimensionlessPositiveFloat,
    ) -> DimensionlessPositiveFloat {
        PositiveFloat(1.0) - self.Omega_M0 - self.Omega_DE0 - omega_nu0 - omega_gamma0
    }
}
