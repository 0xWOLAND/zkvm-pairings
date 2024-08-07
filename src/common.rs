use crate::fp::{Bls12381, Bn254, FpElement};

pub trait AffinePoint: FpElement {
    const X: u64;
    // const BETA: Self;
    // const G1_X: Self;
    // const G1_Y: Self;
    // const G2_X0: Self;
    // const G2_X1: Self;
    // const G2_Y0: Self;
    // const G2_Y1: Self;
}

impl AffinePoint for Bls12381 {
    const X: u64 = 0xd201000000010000; // -0xd201000000010000
}

impl AffinePoint for Bn254 {
    const X: u64 = 0x7;
}
