use std::ops::{Add, Sub};

use crate::units::{macros::floating_point_unit_impl, traits::FloatingPointUnit};

floating_point_unit_impl! { Meter }
floating_point_unit_impl! { Kilometer }
floating_point_unit_impl! { Mpc }

// Conversions
pub const KILOMETER_TO_METER: f64 = 1000.;
pub const MPC_TO_METERS: f64 = 3.086e+22;
pub const MPC_TO_KILOMETERS: f64 = 3.086e+19;

impl From<Kilometer> for Meter {
    fn from(km: Kilometer) -> Meter {
        Meter(1000. * km.0)
    }
}

impl From<Mpc> for Meter {
    fn from(mpc: Mpc) -> Self {
        Meter(mpc.0 * MPC_TO_METERS)
    }
}

impl From<Mpc> for Kilometer {
    fn from(mpc: Mpc) -> Self {
        Kilometer(mpc.0 * MPC_TO_KILOMETERS)
    }
}
