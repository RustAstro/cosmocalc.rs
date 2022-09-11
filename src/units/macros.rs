macro_rules! floating_point_unit_impl {
    ($outer : ident) => {
        #[derive(Clone, Copy, Debug, PartialOrd, PartialEq)]
        pub struct $outer(pub f64);

        impl FloatingPointUnit for $outer {
            fn new(x: f64) -> Self {
                Self(x)
            }

            fn inner(&self) -> f64 {
                self.0
            }
        }

        impl Add<$outer> for $outer {
            type Output = $outer;

            fn add(self, b: $outer) -> $outer {
                $outer((self.0.add(&b.0)))
            }
        }

        impl Sub<$outer> for $outer {
            type Output = $outer;

            fn sub(self, b: $outer) -> $outer {
                $outer((self.0.sub(&b.0)))
            }
        }
    };
}

pub(crate) use floating_point_unit_impl;
