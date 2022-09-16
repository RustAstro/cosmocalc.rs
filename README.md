# cosmocalc.rs

[![Crates.io][crates-badge]][crates-url]
[![Documentation][docs-badge]][docs-url]

[crates-badge]: https://img.shields.io/crates/v/cosmocalc.svg
[crates-url]: https://crates.io/crates/cosmocalc
[docs-badge]: https://docs.rs/cosmocalc/badge.svg
[docs-url]: https://docs.rs/cosmocalc

A library for computing quantities in cosmology in the Rust programming language

# Features

## Cosmological distance calculations

```rust
let cosmology = FLRWCosmology::two_component(0.286, 0.714, 69.6);

z = Redshift::new(2.0);

let d_c = cosmology.radial_comoving_distance(z);
let d_m = cosmology.transverse_comoving_distance(z);
let d_a = cosmology.angular_diameter_distance(z);
let d_l = cosmology.luminosity_distance(z);
let v = cosmology.comoving_volume(z);
```

## Set contributions from matter, dark energy and relativistic particles for flat or non-flat cosmology

```rust
let omega_m = 0.299;
let omega_de = 0.7;
let omega_baryon = 0.05;
let H_0 = 69.6;
let T_CMB0 = 2.7255;
let omegas = OmegaFactors::new(omega_m, omega_de, omega_baryon).unwrap();
let cosmology = FLRWCosmology::new(
    None,
    None,
    H_0,
    omegas,
    Some(T_CMB0),
    Some(PositiveFloat(0.)),
    Some(vec![]),
)
.unwrap();

z = Redshift::new(2.0);

let t = cosmology.lookback_time(z);
let omega_at_z = cosmology.omega_tot(z);
let omega_de_at_z = cosmology.omega_de(z);
let critical_density_at_z = cosmology.critical_density(z);
let T_CMB_at_z = cosmology.T_CMB(z);
let T_nu_at_z = cosmology.T_nu(z);
let d_H = cosmology.hubble_distance();
let t_H = cosmology.hubble_time();
let expansion_rate_at_z = cosmology.H(z);
let a_z = cosmology.scale_factor(z);
```

# Developers

## Dev setup

This project requires [`Rust`](https://www.rust-lang.org/tools/install). Once installed:

```
cargo build
cargo test
```

## Benchmarks

Run `criterion` benchmarks using:

```
cargo bench
```

This will generate a report at `target/criterion/report/index.html`.
