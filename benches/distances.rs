use cosmocalc::{cosmology::OmegaFactors, Distances, FLRWCosmology, FloatingPointUnit, Redshift};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

pub fn bench_luminosity_distance(c: &mut Criterion) {
    let mut group = c.benchmark_group("d_L");
    let omegas = OmegaFactors::new(0.27, 0.73, 0.044).unwrap();
    let cosmology = FLRWCosmology::new(None, None, 70.0, omegas, None, None, None).unwrap();

    let z_1 = Redshift::new(1.);
    group.bench_with_input(
        BenchmarkId::new("d_L(z=1)", format!("{:?}", z_1)),
        &z_1,
        |b, z| b.iter(|| cosmology.luminosity_distance(*z)),
    );
    let z_2 = Redshift::new(2.);
    group.bench_with_input(
        BenchmarkId::new("d_L(z=2)", format!("{:?}", z_2)),
        &z_2,
        |b, z| b.iter(|| cosmology.luminosity_distance(*z)),
    );
    group.finish();
}

criterion_group!(benches, bench_luminosity_distance);
criterion_main!(benches);
