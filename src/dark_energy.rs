use crate::Redshift;

pub trait DarkEnergyEquationOfState {
    fn w(&self, z: Vec<Redshift>) -> Vec<f32>;
}

// TODO: impl DarkEnergyEquationOfState for Cosmology
