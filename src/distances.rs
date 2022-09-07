use crate::{units::PositiveFloat, FLRWCosmology, Mpc, Redshift};

/// Bin width in redshift integrals.
const DZ: f32 = 0.0001;

/// Cosmological distances following Hogg 2000
/// https://arxiv.org/pdf/astro-ph/9905116.pdf
pub trait Distances {
    /// Line of sight (radial) comoving distance in Megaparsecs.
    fn radial_comoving_distance(&self, z: Redshift) -> Mpc;
    /// Transverse comoving distance in Megaparsecs.
    fn transverse_comoving_distance(&self, z: Redshift) -> Mpc;
    /// Angular diameter distance in Megaparsecs.
    fn angular_diameter_distance(&self, z: Redshift) -> Mpc;
    /// Luminosity distance in Megaparsecs.
    fn luminosity_distance(&self, z: Redshift) -> Mpc;
}

impl Distances for FLRWCosmology {
    fn radial_comoving_distance(&self, z: Redshift) -> Mpc {
        let mut integrand: f32 = 0.0;
        let mut z_prime = 0.0;
        while z_prime < z {
            z_prime += DZ / 2.;
            integrand += (DZ / 2.) / self.E(z_prime).0;
        }
        PositiveFloat(self.hubble_distance() * integrand)
    }

    fn transverse_comoving_distance(&self, _z: Redshift) -> Mpc {
        todo!()
    }

    fn angular_diameter_distance(&self, _z: Redshift) -> Mpc {
        todo!()
    }

    fn luminosity_distance(&self, _z: Redshift) -> Mpc {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::cosmology::OmegaFactors;

    use super::*;

    #[test]
    fn flat_universe_distances() {
        let omegas = OmegaFactors::new(0.286, 0.714, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        // Megaparsecs
        // Ned Wright calculator yields 6481.1
        assert_eq!(cosmology.radial_comoving_distance(3.0).0, 6482.5117);
    }
}
