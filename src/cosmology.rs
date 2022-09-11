use anyhow::anyhow;

mod omega_factors;

pub use omega_factors::OmegaFactors;

use crate::{
    constants::{self, C_M_PER_S, DEFAULT_NEUTRINO_MASSES, DEFAULT_N_EFF},
    eV, units,
    units::{FloatingPointUnit, PositiveFloat, Seconds},
    DimensionlessFloat, DimensionlessPositiveFloat, HInvKmPerSecPerMpc, Kelvin, KmPerSecPerMpc,
    Mpc, Redshift,
};

/// Represents an FLRW cosmology.
///
/// This represents an homogenous and isotropic cosmology based
/// on the FLRW (Friedmann-Lemaitre-Robertson-Walker) metric.
///
/// # Examples
///
/// ```
/// use cosmocalc::{Distances, FLRWCosmology};
///
/// let cosmology = FLRWCosmology::two_component(0.286, 0.714, 69.6);
/// assert!(cosmology.radial_comoving_distance(2.0).0 > 5273.);
/// assert!(cosmology.radial_comoving_distance(2.0).0 < 5274.);
/// ```
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
    pub T_CMB0: Option<Kelvin>,
    /// Number of effective neutrino species.
    pub N_eff: DimensionlessPositiveFloat,
    /// Mass of neutrino species in eV.
    pub m_nu: Vec<eV>,
}

impl FLRWCosmology {
    /// Instantiate a simple two component cosmology.
    pub fn two_component(Omega_M0: f64, Omega_DE0: f64, H_0: f64) -> Self {
        let omega = OmegaFactors::new(Omega_M0, Omega_DE0, Omega_M0).unwrap();
        Self {
            name: None,
            reference: None,
            H_0,
            omega,
            T_CMB0: Some(PositiveFloat::zero()),
            N_eff: PositiveFloat::zero(),
            m_nu: vec![],
        }
    }

    /// Instantiate a new FLRW cosmology.
    pub fn new(
        name: Option<String>,
        reference: Option<String>,
        H_0: f64,
        omega: OmegaFactors,
        T_CMB0: Option<f64>,
        N_eff: Option<DimensionlessPositiveFloat>,
        m_nu: Option<Vec<eV>>,
    ) -> Result<Self, anyhow::Error> {
        if N_eff.unwrap_or(*DEFAULT_N_EFF).floor()
            != m_nu
                .clone()
                .unwrap_or_else(|| DEFAULT_NEUTRINO_MASSES.to_vec())
                .len() as f64
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
            T_CMB0: T_CMB0.map(PositiveFloat),
            N_eff: N_eff.unwrap_or(*DEFAULT_N_EFF),
            m_nu: m_nu.unwrap_or_else(|| DEFAULT_NEUTRINO_MASSES.to_vec()),
        })
    }

    pub fn E(&self, z: Redshift) -> Mpc {
        PositiveFloat(
            (self.omega_m0().0 * (1. + z).powf(3.)
                + self.omega_k0().0 * (1. + z).powf(2.)
                + self.omega_de0().0
                + (self.omega_gamma0().0 + self.omega_nu0().0) * (1. + z).powf(4.))
            .sqrt(),
        )
    }

    /// Hubble expansion rate (km/s/Mpc) at redshift z.
    pub fn H(&self, z: Redshift) -> KmPerSecPerMpc {
        self.H_0 * self.E(z).0
    }

    /// Scale factor at redshift z.
    pub fn scale_factor(&self, z: Redshift) -> DimensionlessPositiveFloat {
        PositiveFloat(1.0 / (z + 1.0))
    }

    /// Dimensionless hubble parameter h where 100 km/s/Mpc * h = H0
    pub fn little_h(&self) -> DimensionlessPositiveFloat {
        DimensionlessPositiveFloat::new(self.H_0 / 100.0).unwrap()
    }

    /// Hubble time: Inverse of the Hubble constant H_0
    pub fn hubble_time(&self) -> Seconds {
        // H_0 units are km/s/Mpc so we need to convert Mpc to km
        // such that the distance units cancel
        PositiveFloat(1. / self.H_0 * units::MPC_TO_KILOMETERS)
    }

    /// Hubble distance in Mpc: $D_H = c / H_0$.
    pub fn hubble_distance(&self) -> KmPerSecPerMpc {
        // Factor of 1000 to convert c in m/s to c in km/s so that
        // the units cancel.
        C_M_PER_S / (self.H_0 * units::KILOMETER_TO_METER)
    }

    /// Hubble distance in h^{-1} Mpc.
    pub fn hubble_distance_little_h(&self) -> HInvKmPerSecPerMpc {
        C_M_PER_S / (1.0e5)
    }

    /// Critical mass density at redshift z.
    pub fn critical_density(&self, z: Redshift) -> units::KilogramsPerMeter3 {
        if z == 0.0 {
            PositiveFloat(
                3. * self.H_0.powf(2.)
                    / (8. * constants::PI * constants::G * units::MPC_TO_KILOMETERS.powf(2.)),
            )
        } else {
            PositiveFloat(
                3. * self.H(z).powf(2.)
                    / (8. * constants::PI * constants::G * units::MPC_TO_KILOMETERS.powf(2.)),
            )
        }
    }

    /// Dimensionless photon density (density/critical density) at `z=0`.
    ///
    /// Eqn. 2.28 from Ryden divided by the critical density at `z=0`
    pub fn omega_gamma0(&self) -> DimensionlessFloat {
        match self.T_CMB0 {
            Some(T_CMB0) => DimensionlessFloat(
                *constants::ALPHA * T_CMB0.powf(4.)
                    / (self.critical_density(0.).0 * C_M_PER_S.powf(2.)),
            ),
            None => DimensionlessFloat::zero(),
        }
    }

    /// Dimensionless photon density (density/critical density) at `z>0`
    pub fn omega_gamma(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_gamma0().0 * (1.0 + z).powf(4.) * 1. / self.E(z).0.powf(2.))
    }

    /// Dimensionless neutrino density (density/critical density) at `z=0`
    pub fn omega_nu0(&self) -> DimensionlessFloat {
        match self.T_CMB0 {
            Some(_) => DimensionlessFloat(
                7. / 8. * (4.0f64 / 11.).powf(4. / 3.) * self.N_eff.0 * self.omega_gamma0().0,
            ),
            None => DimensionlessFloat::zero(),
        }
    }

    /// Dimensionless neutrino density (density/critical density) at `z>0`
    pub fn omega_nu(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_nu0().0 * (1.0 + z).powf(4.) * 1. / self.E(z).0.powf(2.))
    }

    /// Dimensionless dark matter density (density/critical density) at `z=0`
    pub fn omega_dm0(&self) -> DimensionlessFloat {
        self.omega.omega_dark_matter_density_0()
    }

    /// Dimensionless dark matter density (density/critical density) at `z>0`
    pub fn omega_dm(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_dm0().0 * (1.0 + z).powf(3.) * 1. / self.E(z).0.powf(2.))
    }

    /// Dimensionless effective curvature density (density/critical density) at `z=0`
    pub fn omega_k0(&self) -> DimensionlessFloat {
        self.omega
            .curvature_density_0(self.omega_nu0(), self.omega_gamma0())
    }

    /// Dimensionless effective curvature density (density/critical density) at `z>0`
    pub fn omega_k(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_k0().0 * (1.0 + z).powf(2.) * 1. / self.E(z).0.powf(2.))
    }

    /// Dimensionless matter density (density/critical density) at `z=0`
    pub fn omega_m0(&self) -> DimensionlessFloat {
        self.omega.Omega_M0
    }

    /// Dimensionless matter density (density/critical density) at `z>0`
    pub fn omega_m(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_m0().0 * (1.0 + z).powf(3.) * 1. / self.E(z).0.powf(2.))
    }

    /// Dimensionless baryon density (density/critical density) at `z=0`
    pub fn omega_b0(&self) -> DimensionlessFloat {
        self.omega.Omega_b0
    }

    /// Dimensionless baryon density (density/critical density) at `z>0`
    pub fn omega_b(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_b0().0 * (1.0 + z).powf(3.) * 1. / self.E(z).0.powf(2.))
    }

    /// Dimensionless dark energy density (density/critical density) at `z=0`
    pub fn omega_de0(&self) -> DimensionlessFloat {
        self.omega.Omega_DE0
    }

    /// Dimensionless dark energy density (density/critical density) at `z>0`.
    pub fn omega_de(&self, z: Redshift) -> DimensionlessFloat {
        DimensionlessFloat(self.omega_de0().0 / self.E(z).0.powf(2.))
    }

    /// Dimensionless total density (density/critical density) at `z=0`.
    pub fn omega_tot0(&self) -> DimensionlessFloat {
        self.omega_m0()
            + self.omega_gamma0()
            + self.omega_nu0()
            + self.omega_de0()
            + self.omega_k0()
    }

    /// Dimensionless total density (density/critical density) at `z>0`.
    pub fn omega_tot(&self, z: Redshift) -> DimensionlessFloat {
        self.omega_m(z)
            + self.omega_gamma(z)
            + self.omega_nu(z)
            + self.omega_de(z)
            + self.omega_k(z)
    }

    /// Whether this cosmology is spatially flat
    pub fn is_flat(&self) -> bool {
        self.omega_k0() == DimensionlessFloat::zero()
            && self.omega_tot0() == DimensionlessFloat::one()
    }

    /// Neutrino temperature at `z=0`.
    pub fn neutrino_temperature0(&self) -> Kelvin {
        match self.T_CMB0 {
            Some(T_cmb) => PositiveFloat(T_cmb.0 * (*constants::T_NU_TO_T_GAMMA_RATIO).0),
            None => DimensionlessPositiveFloat::zero(),
        }
    }
}
