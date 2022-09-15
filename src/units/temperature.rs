use std::{
    default::Default,
    ops::{Add, Sub},
};

use crate::units::{macros::floating_point_unit_impl, traits::FloatingPointUnit};

floating_point_unit_impl! { Kelvin }

impl Default for Kelvin {
    fn default() -> Self {
        Self(0.)
    }
}
