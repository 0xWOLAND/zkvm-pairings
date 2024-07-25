use crate::fp::Fp;

pub trait Curve: Clone + Copy {
    const A: [u64; 6];
    const MAX_BITS: u64;
    const MODULUS: [u64; 6];
    const INV: u64;
    const R: [u64; 6];
}

#[derive(Clone, Copy)]
pub struct Bls12381Curve {}

impl Curve for Bls12381Curve {
    const A: [u64; 6] = [0; 6];
    const MAX_BITS: u64 = 384;
    const MODULUS: [u64; 6] = [
        0xb9fe_ffff_ffff_aaab,
        0x1eab_fffe_b153_ffff,
        0x6730_d2a0_f6b0_f624,
        0x6477_4b84_f385_12bf,
        0x4b1b_a7b6_434b_acd7,
        0x1a01_11ea_397f_e69a,
    ];

    /// INV = -(p^{-1} mod 2^64) mod 2^64
    const INV: u64 = 0x89f3_fffc_fffc_fffd;

    /// R = 2^384 mod p
    const R: [u64; 6] = [
        0x7609_0000_0002_fffd,
        0xebf4_000b_c40c_0002,
        0x5f48_9857_53c7_58ba,
        0x77ce_5853_7052_5745,
        0x5c07_1a97_a256_ec6d,
        0x15f6_5ec3_fa80_e493,
    ];
}
