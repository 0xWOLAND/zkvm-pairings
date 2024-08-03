use field::{
    common::{Bls12381Curve, Curve as FieldCurve},
    fp12::Fp12,
};

pub trait Curve: FieldCurve {
    const LOG_ATE_LOOP_COUNT: usize;
    const LOOP_COUNTER: [u64; 64];
}

impl Curve for Bls12381Curve {
    const LOG_ATE_LOOP_COUNT: usize = 62;
    const LOOP_COUNTER: [u64; 64] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        1, 0, 1, 1,
    ];
}
