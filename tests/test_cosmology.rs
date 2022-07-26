use cosmocalc::{
    cosmology::{FLRWCosmology, OmegaFactors},
    DimensionlessPositiveFloat,
};

#[test]
fn hubble_distance() {
    let omegas = OmegaFactors::new(0.27, 0.044, 0.73).unwrap();
    let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();

    assert_eq!(
        cosmology.little_h(),
        DimensionlessPositiveFloat::new(0.70).unwrap()
    );
    // Should be around 3000 h^-1 Mpc
    assert!(cosmology.hubble_distance_little_h() > 2999.0);
    assert!(cosmology.hubble_distance_little_h() > 3001.0);

    // D_H in units of h^{-1} Mpc should be equal to D_H in units of Mpc
    // TODO: Fix
    // assert_eq!(
    //     cosmology.hubble_distance_little_h(),
    //     cosmology.hubble_distance() / cosmology.little_h()
    // );
}
