use crate::{units::PositiveFloat, DimensionlessFloat, FLRWCosmology, Mpc, Redshift};

/// Bin width in redshift integrals.
const DZ: f64 = 0.0001;

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
        let mut integrand: f64 = 0.0;
        let mut z_prime = 0.0;
        while z_prime < z {
            z_prime += DZ / 2.;
            integrand += (DZ / 2.) / self.E(z_prime).0;
        }
        PositiveFloat(self.hubble_distance() * integrand)
    }

    fn transverse_comoving_distance(&self, z: Redshift) -> Mpc {
        let radial_comoving = self.radial_comoving_distance(z);
        let omega_k = self.omega_k(z);
        if omega_k > DimensionlessFloat::zero() {
            // Negative curvature (open)
            let sqrt_omega_k = (omega_k.0).sqrt();
            PositiveFloat(
                self.hubble_distance() * 1. / sqrt_omega_k
                    * f64::sinh(sqrt_omega_k * radial_comoving.0 / self.hubble_distance()),
            )
        } else if omega_k == DimensionlessFloat::zero() {
            // Flat
            radial_comoving
        } else {
            // Positive curvature (closed)
            let abs_sqrt_omega_k = (-1. * omega_k.0).sqrt();
            PositiveFloat(
                self.hubble_distance() * 1. / abs_sqrt_omega_k
                    * f64::sin(abs_sqrt_omega_k * radial_comoving.0 / self.hubble_distance()),
            )
        }
    }

    fn angular_diameter_distance(&self, z: Redshift) -> Mpc {
        PositiveFloat(self.transverse_comoving_distance(z).0 / (1. + z))
    }

    fn luminosity_distance(&self, z: Redshift) -> Mpc {
        // TODO: K-CORRECTIONS
        PositiveFloat(self.transverse_comoving_distance(z).0 * (1. + z))
    }
}

#[cfg(test)]
mod tests {
    use crate::{cosmology::OmegaFactors, DimensionlessPositiveFloat};

    use super::*;

    #[test]
    fn flat_universe_distances_no_relativistic_contribution() {
        // TESTED vs: astro.py 5.1 FlatLambdaCDM
        let omegas = OmegaFactors::new(0.286, 0.714, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(3.0).0 > 6482.5);
        assert!(cosmology.radial_comoving_distance(3.0).0 < 6482.6);
        assert!(cosmology.angular_diameter_distance(3.0).0 > 1620.6);
        assert!(cosmology.angular_diameter_distance(3.0).0 < 1620.7);
        assert!(cosmology.luminosity_distance(3.0).0 > 25930.0);
        assert!(cosmology.luminosity_distance(3.0).0 < 25930.2);
    }

    #[test]
    fn flat_universe_distances_with_radiation_but_no_neutrinos() {
        // TESTED vs: astro.py 5.1 FlatLambdaCDM
        let omegas = OmegaFactors::new(0.299, 0.7, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(
            None,
            None,
            69.6,
            omegas,
            Some(2.7255),
            Some(PositiveFloat(0.)),
            Some(vec![]),
        )
        .unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(3.0).0 > 6395.0);
        assert!(cosmology.radial_comoving_distance(3.0).0 < 6399.0);
        assert!(cosmology.angular_diameter_distance(3.0).0 > 1599.0);
        assert!(cosmology.angular_diameter_distance(3.0).0 < 1600.0);
        assert!(cosmology.luminosity_distance(3.0).0 > 25588.);
        assert!(cosmology.luminosity_distance(3.0).0 < 25589.);
    }

    #[test]
    fn flat_universe_distances_with_radiation_and_neutrinos() {
        let omegas = OmegaFactors::new(0.25, 0.7, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(
            None,
            None,
            69.6,
            omegas,
            Some(2.7255),
            Some(PositiveFloat(3.04)),
            Some(vec![
                DimensionlessPositiveFloat::zero(),
                DimensionlessPositiveFloat::zero(),
                DimensionlessPositiveFloat::zero(),
            ]),
        )
        .unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(3.0).0 > 6598.);
        assert!(cosmology.radial_comoving_distance(3.0).0 < 6598.5);
        assert!(cosmology.angular_diameter_distance(3.0).0 > 1600.5);
        assert!(cosmology.angular_diameter_distance(3.0).0 < 1700.0);
        assert!(cosmology.luminosity_distance(3.0).0 > 25000.);
        assert!(cosmology.luminosity_distance(3.0).0 < 27000.);
    }

    #[test]
    fn open_universe_distances_no_relativistic_contribution() {
        let omegas = OmegaFactors::new(0.286, 0.0, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(3.0).0 > 5200.);
        assert!(cosmology.radial_comoving_distance(3.0).0 < 5300.);
        assert!(cosmology.angular_diameter_distance(3.0).0 > 1250.);
        assert!(cosmology.angular_diameter_distance(3.0).0 < 1600.);
        // No k-corrections here
        assert!(cosmology.luminosity_distance(3.0).0 > 22000.);
        assert!(cosmology.luminosity_distance(3.0).0 < 24000.);
    }

    #[test]
    fn closed_universe_distances_no_relativistic_contribution() {
        let omegas = OmegaFactors::new(0.286, 0.8, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(2.0).0 > 5000.);
        assert!(cosmology.radial_comoving_distance(2.0).0 < 6000.);
        assert!(cosmology.angular_diameter_distance(2.0).0 > 1500.);
        assert!(cosmology.angular_diameter_distance(2.0).0 < 2000.);
        // No k-corrections here
        assert!(cosmology.luminosity_distance(2.0).0 > 14000.);
        assert!(cosmology.luminosity_distance(2.0).0 < 16000.);
    }

    #[test]
    fn simple_two_component() {
        let cosmology = FLRWCosmology::two_component(0.286, 0.714, 69.6);
        assert!(cosmology.radial_comoving_distance(2.0).0 > 5273.);
        assert!(cosmology.radial_comoving_distance(2.0).0 < 5274.);
    }
}
