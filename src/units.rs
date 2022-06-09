use anyhow::anyhow;

#[derive(Clone, Copy, Debug, PartialOrd, PartialEq)]
pub struct PositiveFloat(f32);

impl PositiveFloat {
    pub fn new(x: f32) -> Result<Self, anyhow::Error> {
        if x < 0. {
            return Err(anyhow!("expected positive number"));
        }
        Ok(Self(x))
    }

    // Passthrough methods for convenience
    pub fn floor(&self) -> f32 {
        self.0.floor()
    }
}
