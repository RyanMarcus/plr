use criterion::{black_box, criterion_group, criterion_main, Criterion};
use plr::OptimalPLR;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut data = Vec::new();

    for idx in 0..1000 {
        let x = idx as f64 / 1000.0;
        let y = f64::sin(x);

        data.push((x, y));
    }
    
    c.bench_function("optimal 1000", |b| b.iter(|| {
        let mut plr = OptimalPLR::new(0.005);

        for &(x, y) in data.iter() {
            plr.process(x, y);
        }

        plr.finish();
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
