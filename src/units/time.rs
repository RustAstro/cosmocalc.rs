use std::ops::{Add, Sub};

use crate::units::{macros::floating_point_unit_impl, traits::FloatingPointUnit};

floating_point_unit_impl! { Gyr }
floating_point_unit_impl! { Seconds }

pub const SECONDS_PER_YR: f64 = 3.154e+7;
pub const SECONDS_PER_GYR: f64 = 3.154e+16;

impl From<Seconds> for Gyr {
    fn from(seconds: Seconds) -> Self {
        Gyr::new(seconds.0 / SECONDS_PER_GYR)
    }
}
