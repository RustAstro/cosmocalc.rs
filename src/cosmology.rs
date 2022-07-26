use anyhow::anyhow;

use crate::{
    constants::{C_M_PER_S, DEFAULT_NEUTRINO_MASSES, DEFAULT_N_EFF, ONE, ZERO},
    eV, DimensionlessPositiveFloat, HInvKmPerSecPerMpc, KmPerSecPerMpc, TemperatureKelvin,
};

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
    pub fn new(Omega_M0: f32, Omega_DE0: f32, Omega_b0: f32) -> Result<Self, anyhow::Error> {
        if Omega_b0 < Omega_M0 {
            return Err(anyhow!("cannot have more baryons than matter"));
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

    /// Curvature density at `z=0` given neutrino density at `z=0`.
    pub fn curvature_density_0(
        &self,
        omega_nu0: DimensionlessPositiveFloat,
    ) -> DimensionlessPositiveFloat {
        todo!()
    }
}

/// Represents an FLRW cosmology.
///
/// This represents an homogenous and isotropic cosmology based
/// on the FLRW (Friedmann-Lemaitre-Robertson-Walker) metric.
pub struct FLRWCosmology {
    /// A descriptive name.
    pub name: Option<String>,
    /// Literature reference.
    pub reference: Option<String>,

    /// Hubble constant at `z=0` (km/(s/Mpc)).
    pub H_0: KmPerSecPerMpc,

    /// Omega factors for this cosmology.
    pub omega: OmegaFactors,

    /// Temperature of the CMB at `z=0`.
    pub T_CMB0: TemperatureKelvin,
    /// Number of effective neutrino species.
    pub N_eff: DimensionlessPositiveFloat,
    /// Mass of neutrino species in eV.
    pub m_nu: Vec<eV>,
}

impl FLRWCosmology {
    /// Instantiate a new FLRW cosmology.
    pub fn new(
        name: Option<String>,
        reference: Option<String>,
        H_0: f32,
        omega: OmegaFactors,
        T_CMB0: Option<TemperatureKelvin>,
        N_eff: Option<DimensionlessPositiveFloat>,
        m_nu: Option<Vec<eV>>,
    ) -> Result<Self, anyhow::Error> {
        if N_eff.unwrap_or(*DEFAULT_N_EFF).floor()
            != m_nu
                .clone()
                .unwrap_or(DEFAULT_NEUTRINO_MASSES.to_vec())
                .len() as f32
        {
            return Err(anyhow!(
                "number of neutrino masses must match the number of effective neutrino species"
            ));
        }

        Ok(Self {
            name,
            reference,
            H_0,
            omega,
            T_CMB0: T_CMB0.unwrap_or(*ZERO),
            N_eff: N_eff.unwrap_or(*DEFAULT_N_EFF),
            m_nu: m_nu.unwrap_or(DEFAULT_NEUTRINO_MASSES.to_vec()),
        })
    }

    /// Dimensionless hubble parameter h where 100 km/s/Mpc * h = H0
    pub fn little_h(&self) -> DimensionlessPositiveFloat {
        DimensionlessPositiveFloat::new(self.H_0 / 100.0).unwrap()
    }

    /// Hubble time: Inverse of the Hubble constant H_0
    pub fn hubble_time(&self) {
        todo!()
    }

    /// Hubble distance in Mpc: $D_H = c / H_0$.
    pub fn hubble_distance(&self) -> KmPerSecPerMpc {
        // Factor of 1000 to convert c in m/s to c in km/s so that
        // the units cancel.
        //*C_M_PER_S / (1000. * self.H_0)
        *C_M_PER_S / (self.H_0)
    }

    /// Hubble distance in h^{-1} Mpc.
    pub fn hubble_distance_little_h(&self) -> HInvKmPerSecPerMpc {
        *C_M_PER_S / (10. * self.little_h().0)
    }

    /// Critical density at `z=0`.
    pub fn critical_density0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless photon density (density/critical density) at `z=0`.
    pub fn omega_gamma0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless neutrino density (density/critical density) at `z=0`
    pub fn omega_nu0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless dark matter density (density/critical density) at `z=0`
    pub fn omega_dm0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless effective curvature density (density/critical density) at `z=0`
    pub fn omega_k0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless matter density (density/critical density) at `z=0`
    pub fn omega_m0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless dark energy density (density/critical density) at `z=0`
    pub fn omega_de0(&self) -> DimensionlessPositiveFloat {
        todo!()
    }

    /// Dimensionless total density (density/critical density) at `z=0`.
    pub fn omega_tot0(&self) -> DimensionlessPositiveFloat {
        self.omega_m0()
            + self.omega_gamma0()
            + self.omega_nu0()
            + self.omega_de0()
            + self.omega_k0()
    }

    /// Whether this cosmology is spatially flat
    pub fn is_flat(&self) -> bool {
        self.omega_k0() == *ZERO && self.omega_tot0() == *ONE
    }

    pub fn neutrino_temperature(&self) {
        todo!()
    }
}
