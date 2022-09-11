pub trait FloatingPointUnit {
    /// Create the value.
    fn new(inner: f64) -> Self;

    /// Get the inner unit.
    fn inner(&self) -> f64;

    /// Default implementations

    /// Get the zero value for this unit.
    fn zero() -> Self
    where
        Self: std::marker::Sized,
    {
        Self::new(0.)
    }

    /// Get the one value for this unit.
    fn one() -> Self
    where
        Self: std::marker::Sized,
    {
        Self::new(1.)
    }

    /// Round the value down.
    fn floor(&self) -> f64 {
        self.inner().floor()
    }

    /// The inner part of this value raised to a power.
    fn powf(&self, exp: f64) -> f64 {
        self.inner().powf(exp)
    }
}
