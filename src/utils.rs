use crate::{constants, units::PositiveFloat, Joule, Kilogram};

/// Convert energy to mass using $E=mc^2$
pub fn energy_to_mass(energy: Joule) -> Kilogram {
    PositiveFloat(energy.0 / (constants::C_M_PER_S).powf(2.))
}

/// Convert mass to energy using $E=mc^2$
pub fn mass_to_energy(mass: Kilogram) -> Joule {
    PositiveFloat(mass.0 * (constants::C_M_PER_S).powf(2.))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_energy_equivalance() {
        // ~1kg
        assert!(energy_to_mass(PositiveFloat(8.9e16)).0 > 0.99);
        assert!(energy_to_mass(PositiveFloat(8.9e16)).0 < 1.01);

        assert!(mass_to_energy(PositiveFloat(1.)).0 > 8.9e16);
        assert!(mass_to_energy(PositiveFloat(1.)).0 < 9.05e16);
    }
}
