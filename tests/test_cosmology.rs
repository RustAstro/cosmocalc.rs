use cosmocalc::{
    cosmology::{FLRWCosmology, OmegaFactors},
    DimensionlessPositiveFloat, FloatingPointUnit, Gyr, Redshift,
};

#[test]
fn hubble_distance_and_time() {
    let omegas = OmegaFactors::new(0.27, 0.73, 0.044).unwrap();
    let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();

    assert!(cosmology.is_flat());

    // Little h for 70 km s^{-1} Mpc^{-1} should be 0.70
    assert_eq!(
        cosmology.little_h(),
        DimensionlessPositiveFloat::new(0.70).unwrap()
    );

    // Should be around 3000 h^-1 Mpc
    assert!(cosmology.hubble_distance_little_h() > 2950.0);
    assert!(cosmology.hubble_distance_little_h() < 3000.0);

    // For H_0 = 70, should be 4285.7 Mpc
    assert!(cosmology.hubble_distance() > 4000.0);
    assert!(cosmology.hubble_distance() < 4500.0);

    // D_H in units of h^{-1} Mpc should be equal to D_H in units of Mpc
    assert!(
        cosmology.hubble_distance_little_h()
            - (cosmology.hubble_distance() * cosmology.little_h().0)
            < 0.01
    );

    // t_H should be ~3e17 h^-1 seconds so t_H = 4e17 if h=0.70
    assert!(cosmology.hubble_time().0 > 4.4e17);
    assert!(cosmology.hubble_time().0 < 4.5e17);
}

#[test]
fn densities() {
    let omegas = OmegaFactors::new(0.27, 0.73, 0.044).unwrap();
    let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();

    assert!(cosmology.critical_density(Redshift::zero()).0 > 8.7e-27);
    assert!(cosmology.critical_density(Redshift::zero()).0 < 9.5e-27);
}

#[test]
fn lookback_time() {
    // TESTED vs: astro.py 5.1 FlatLambdaCDM
    // Accurate to 0.1 Gyr
    let omegas = OmegaFactors::new(0.27, 0.73, 0.044).unwrap();
    let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();
    assert!(cosmology.lookback_time(Redshift::zero()) == Gyr::zero());
    assert!(cosmology.lookback_time(Redshift::new(3.0)) > Gyr::new(11.64));
    assert!(cosmology.lookback_time(Redshift::new(3.0)) < Gyr::new(11.65));
}
