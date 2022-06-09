use anyhow::anyhow;

use crate::{
    constants::{DEFAULT_NEUTRINO_MASSES, DEFAULT_N_EFF, ZERO},
    PositiveFloat,
};

/// Represents an FLRW cosmology.
///
/// This represents an homogenous and isotropic cosmology based
/// on the FLRW (Friedmann-Lemaitre-Robertson-Walker) metric.
struct FLRWCosmology {
    /// A descriptive name.
    name: Option<String>,
    /// Literature reference.
    reference: Option<String>,

    /// Hubble constant at `z=0`.
    H_0: f32,

    /// Ratio of non-relativistic matter to critical density at `z=0`.
    Omega_M0: PositiveFloat,
    /// Ratio of dark matter density to critical density at `z=0`.
    Omega_DE0: PositiveFloat,
    /// Ratio of baryon density to critical density at `z=0`.
    Omega_b0: PositiveFloat,

    /// Temperature of the CMB at `z=0`.
    T_CMB0: PositiveFloat,
    /// Number of effective neutrino species.
    N_eff: PositiveFloat,
    /// Mass of neutrino species in eV.
    m_nu: Vec<PositiveFloat>,
}

impl FLRWCosmology {
    /// Instantiate a new FLRW cosmology.
    pub fn new(
        name: Option<String>,
        reference: Option<String>,
        H_0: f32,
        Omega_M0: PositiveFloat,
        Omega_DE0: PositiveFloat,
        Omega_b0: PositiveFloat,
        T_CMB0: Option<PositiveFloat>,
        N_eff: Option<PositiveFloat>,
        m_nu: Option<Vec<PositiveFloat>>,
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

        if Omega_b0 < Omega_M0 {
            return Err(anyhow!("cannot have more baryons than matter"));
        }

        Ok(Self {
            name,
            reference,
            H_0,
            Omega_M0,
            Omega_DE0,
            Omega_b0,
            T_CMB0: T_CMB0.unwrap_or(*ZERO),
            N_eff: N_eff.unwrap_or(*DEFAULT_N_EFF),
            m_nu: m_nu.unwrap_or(DEFAULT_NEUTRINO_MASSES.to_vec()),
        })
    }
}
