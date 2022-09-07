use crate::Redshift;

pub trait DarkEnergyEquationOfState {
    fn w(&self, z: Vec<Redshift>) -> Vec<f64>;
}

// TODO: impl DarkEnergyEquationOfState for Cosmology
