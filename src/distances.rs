use crate::{
    constants,
    units::{FloatingPointUnit, Mpc3},
    DimensionlessFloat, FLRWCosmology, Mpc, Redshift,
};

/// Bin width in redshift integrals.
const DZ: f64 = 0.0001;

/// Cosmological distances following [Hogg 2000]
/// [Hogg 2000]: <https://arxiv.org/pdf/astro-ph/9905116.pdf>
pub trait Distances {
    /// Line of sight (radial) comoving distance in Megaparsecs.
    fn radial_comoving_distance(&self, z: Redshift) -> Mpc;
    /// Transverse comoving distance in Megaparsecs.
    fn transverse_comoving_distance(&self, z: Redshift) -> Mpc;
    /// Angular diameter distance in Megaparsecs.
    fn angular_diameter_distance(&self, z: Redshift) -> Mpc;
    /// Luminosity distance in Megaparsecs.
    ///
    /// This should be used with bolometric quantities, i.e.
    /// it does not include K-corrections.
    fn luminosity_distance(&self, z: Redshift) -> Mpc;
    /// Comoving volume.
    fn comoving_volume(&self, z: Redshift) -> Mpc3;
}

impl Distances for FLRWCosmology {
    fn radial_comoving_distance(&self, z: Redshift) -> Mpc {
        // We operate over 1e4
        let max_range = (z.0 as i64) * 10000;
        let step = DZ;
        let integrand: f64 = (0..max_range)
            .map(|z_prime_e4| step / self.E(Redshift::new(z_prime_e4 as f64 / 10000.)).0)
            .sum();
        /*
        Function body of E():
        (self.omega.Omega_M0.0 * (1. + z.0).powi(3)
                + self.omega_k0.0 * (1. + z.0).powi(2)
                + self.omega.Omega_DE0.0
                + (self.omega_gamma0.0 + self.omega_nu0.0) * (1. + z.0).powi(4))
            .sqrt()
        */
        Mpc::new(self.hubble_distance().0 * integrand)
    }

    fn transverse_comoving_distance(&self, z: Redshift) -> Mpc {
        let radial_comoving = self.radial_comoving_distance(z);
        let omega_k = self.omega_k(z);
        if omega_k > DimensionlessFloat::zero() {
            // Negative curvature (open)
            let sqrt_omega_k = (omega_k.0).sqrt();
            Mpc::new(
                self.hubble_distance().0 * 1. / sqrt_omega_k
                    * f64::sinh(sqrt_omega_k * radial_comoving.0 / self.hubble_distance().0),
            )
        } else if omega_k == DimensionlessFloat::zero() {
            // Flat
            radial_comoving
        } else {
            // Positive curvature (closed)
            let abs_sqrt_omega_k = (-1. * omega_k.0).sqrt();
            Mpc::new(
                self.hubble_distance().0 * 1. / abs_sqrt_omega_k
                    * f64::sin(abs_sqrt_omega_k * radial_comoving.0 / self.hubble_distance().0),
            )
        }
    }

    fn angular_diameter_distance(&self, z: Redshift) -> Mpc {
        Mpc::new(self.transverse_comoving_distance(z).0 / (1. + z.0))
    }

    fn luminosity_distance(&self, z: Redshift) -> Mpc {
        // TODO: K-CORRECTIONS
        Mpc::new(self.transverse_comoving_distance(z).0 * (1. + z.0))
    }

    /// Comoving volume
    fn comoving_volume(&self, z: Redshift) -> Mpc3 {
        // KmPerSecPerMpc
        let d_H = self.hubble_distance().0;
        let d_H_cubed = d_H.powi(3);

        let transverse_comoving = self.transverse_comoving_distance(z);
        let omega_k = self.omega_k(z);
        if omega_k > DimensionlessFloat::zero() {
            // Negative curvature (open)
            let sqrt_omega_k = (omega_k.0).sqrt();
            let coefficient = 4. * constants::PI * d_H_cubed / (2. * omega_k.0);
            let term_1_in_parens = transverse_comoving.0 / d_H
                * (1. + omega_k.0 * transverse_comoving.powi(2) / d_H.powi(2)).sqrt();
            let term_2_in_parens =
                1. / sqrt_omega_k * f64::asinh(sqrt_omega_k * transverse_comoving.0 / d_H);

            coefficient * (term_1_in_parens - term_2_in_parens)
        } else if omega_k == DimensionlessFloat::zero() {
            // Flat
            4. * constants::PI * transverse_comoving.powi(3) / 3.
        } else {
            // Positive curvature (closed)
            let sqrt_omega_k = (-1. * omega_k.0).sqrt();
            let coefficient = 4. * constants::PI * d_H_cubed / (2. * omega_k.0);
            let term_1_in_parens = transverse_comoving.0 / d_H
                * (1. + omega_k.0 * transverse_comoving.powi(2) / d_H.powi(2)).sqrt();
            let term_2_in_parens =
                1. / sqrt_omega_k * f64::asin(sqrt_omega_k * transverse_comoving.0 / d_H);

            coefficient * (term_1_in_parens - term_2_in_parens)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{cosmology::OmegaFactors, eV, units::PositiveFloat};

    use super::*;

    #[test]
    fn flat_universe_distances_no_relativistic_contribution() {
        // TESTED vs: astro.py 5.1 FlatLambdaCDM
        let omegas = OmegaFactors::new(0.286, 0.714, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) > Mpc::new(6482.5));
        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) < Mpc::new(6482.8));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) > Mpc::new(1620.6));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) < Mpc::new(1620.7));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) > Mpc::new(25930.0));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) < Mpc::new(25931.0));
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

        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) > Mpc::new(6395.0));
        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) < Mpc::new(6399.0));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) > Mpc::new(1599.0));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) < Mpc::new(1600.0));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) > Mpc::new(25589.));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) < Mpc::new(25594.));
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
            Some(vec![eV::zero(), eV::zero(), eV::zero()]),
        )
        .unwrap();

        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) > Mpc::new(6598.));
        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) < Mpc::new(6598.5));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) > Mpc::new(1600.5));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) < Mpc::new(1700.0));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) > Mpc::new(25000.));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) < Mpc::new(27000.));
    }

    #[test]
    fn open_universe_distances_no_relativistic_contribution() {
        let omegas = OmegaFactors::new(0.286, 0.0, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) > Mpc::new(5200.));
        assert!(cosmology.radial_comoving_distance(Redshift::new(3.0)) < Mpc::new(5300.));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) > Mpc::new(1250.));
        assert!(cosmology.angular_diameter_distance(Redshift::new(3.0)) < Mpc::new(1600.));
        // No k-corrections here
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) > Mpc::new(22000.));
        assert!(cosmology.luminosity_distance(Redshift::new(3.0)) < Mpc::new(24000.));
    }

    #[test]
    fn closed_universe_distances_no_relativistic_contribution() {
        let omegas = OmegaFactors::new(0.286, 0.8, 0.05).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 69.6, omegas, None, None, None).unwrap();

        // Megaparsecs
        assert!(cosmology.radial_comoving_distance(Redshift::new(2.0)) > Mpc::new(5000.));
        assert!(cosmology.radial_comoving_distance(Redshift::new(2.0)) < Mpc::new(6000.));
        assert!(cosmology.angular_diameter_distance(Redshift::new(2.0)) > Mpc::new(1500.));
        assert!(cosmology.angular_diameter_distance(Redshift::new(2.0)) < Mpc::new(2000.));
        // No k-corrections here
        assert!(cosmology.luminosity_distance(Redshift::new(2.0)) > Mpc::new(14000.));
        assert!(cosmology.luminosity_distance(Redshift::new(2.0)) < Mpc::new(16000.));
    }

    #[test]
    fn simple_two_component() {
        let cosmology = FLRWCosmology::two_component(0.286, 0.714, 69.6);
        assert!(cosmology.radial_comoving_distance(Redshift::new(2.0)) > Mpc::new(5273.));
        assert!(cosmology.radial_comoving_distance(Redshift::new(2.0)) < Mpc::new(5274.));
    }

    /// This test is here for profilng purposes. It runs the luminosity distance computation
    /// for many iterations such that profiling data can be collected.
    #[ignore]
    #[test]
    fn simple_luminosity() {
        let omegas = OmegaFactors::new(0.27, 0.73, 0.044).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();
        //for _ in 0..10000000 {
        for _ in 0..1 {
            cosmology.luminosity_distance(Redshift::new(2.0));
        }
    }

    #[test]
    fn comoving_volume() {
        // TESTED vs: astro.py 5.1 FlatLambdaCDM. Within 10e8 Mpc3.
        let omegas = OmegaFactors::new(0.27, 0.73, 0.044).unwrap();
        let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();
        assert!(cosmology.comoving_volume(Redshift::new(3.0)) > 1179361698730.);
        assert!(cosmology.comoving_volume(Redshift::new(3.0)) < 1179470000000.);
    }
}
