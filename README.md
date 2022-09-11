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
assert!(cosmology.radial_comoving_distance(Redshift::new(2.0)) > Mpc::new(5273.));
```

# Developers

## Dev setup

This project requires [`Rust`](https://www.rust-lang.org/tools/install). Once installed:

```
cargo build
cargo test
```
