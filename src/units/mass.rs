use std::ops::{Add, Sub};

use crate::{
    constants,
    units::{macros::floating_point_unit_impl, traits::FloatingPointUnit},
    Joule,
};

floating_point_unit_impl! { Kilogram }
floating_point_unit_impl! { Gram }

impl From<Joule> for Kilogram {
    fn from(energy: Joule) -> Self {
        // Convert between energy and mass using $E=mc^2$
        Kilogram::new(energy.0 / (constants::C_M_PER_S).powi(2))
    }
}

impl From<Kilogram> for Joule {
    fn from(mass: Kilogram) -> Self {
        // Convert between energy and mass using $E=mc^2$
        Joule::new(mass.0 * (constants::C_M_PER_S).powi(2))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_energy_equivalance() {
        // ~1kg
        let equivalent_mass: Kilogram = Joule::new(8.9e16).into();
        assert!(equivalent_mass > Kilogram::new(0.99));
        assert!(equivalent_mass < Kilogram::new(1.01));

        let equivalent_energy: Joule = Kilogram::new(1.).into();
        assert!(equivalent_energy > Joule::new(8.9e16));
        assert!(equivalent_energy < Joule::new(9.05e16));
    }
}
