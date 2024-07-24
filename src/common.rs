use crate::fp::Fp;

pub trait Curve: Clone + Copy {
    const A: Fp;
    const MAX_BITS: u64;
}

#[derive(Clone, Copy)]
pub struct Bls12381Curve {}

impl Curve for Bls12381Curve {
    const A: Fp = Fp::zero();
    const MAX_BITS: u64 = 384;
}
