use crate::{constants, units::mass::Kilogram, FloatingPointUnit, Joule};

/// Convert energy to mass using $E=mc^2$
pub fn energy_to_mass(energy: Joule) -> Kilogram {
    Kilogram::new(energy.0 / (constants::C_M_PER_S).powf(2.))
}

/// Convert mass to energy using $E=mc^2$
pub fn mass_to_energy(mass: Kilogram) -> Joule {
    Joule::new(mass.0 * (constants::C_M_PER_S).powf(2.))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_energy_equivalance() {
        // ~1kg
        assert!(energy_to_mass(Joule::new(8.9e16)).0 > 0.99);
        assert!(energy_to_mass(Joule::new(8.9e16)).0 < 1.01);

        assert!(mass_to_energy(Kilogram::new(1.)).0 > 8.9e16);
        assert!(mass_to_energy(Kilogram::new(1.)).0 < 9.05e16);
    }
}
