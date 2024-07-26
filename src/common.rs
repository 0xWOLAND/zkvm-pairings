use core::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::fr::Fr;

pub(crate) trait FieldElement:
    Add + AddAssign + Sub + SubAssign + Mul + MulAssign + Div + DivAssign + Debug + Sync + Sized
{
}

pub(crate) trait AffinePoint<C: Curve>:
    Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Fr<C>, Output = Self>
    + Copy
    + Debug
{
    type Dtype: FieldElement;
    fn new(x: Self::Dtype, y: Self::Dtype, is_infinity: bool) -> Self;
    fn identity() -> Self;
    fn is_zero(&self) -> bool;
    fn generator() -> Self;
    fn is_valid(&self) -> Result<(), String>;
    fn random(rng: impl rand::Rng) -> Self;
    fn double(&self) -> Self;
}

pub trait Curve: Clone + Copy + Sync + Send + Debug {
    const B: [u64; 6];
    const B2_X: [u64; 6];
    const B2_Y: [u64; 6];
    const X: u64;
    const MAX_BITS: u64;
    const MODULUS: [u64; 6];
    const BETA: [u64; 6];
    const G1_X: [u64; 6];
    const G1_Y: [u64; 6];
    const G2_X0: [u64; 6];
    const G2_X1: [u64; 6];
    const G2_Y0: [u64; 6];
    const G2_Y1: [u64; 6];
    const INV: u64;
    const R: [u64; 6];

    const FR_MODULUS: [u64; 4];
    const FR_BITS: u64;
    const FR_GENERATOR: [u64; 4] = [7, 0, 0, 0];
    const FR_R: [u64; 4];
    const FR_R2: [u64; 4];
    const FR_R3: [u64; 4];
    const FR_INV: u64;
    const FR_TWO_INV: [u64; 4];
    const FR_S: u64;
    const FR_ROOT_OF_UNITY: [u64; 4];
    const FR_ROOT_OF_UNITY_INV: [u64; 4];
    const FR_DELTA: [u64; 4];
}

// BN254 curve R-value
// 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001

#[derive(Clone, Copy, Debug)]
pub struct Bls12381Curve {}

impl Curve for Bls12381Curve {
    const B: [u64; 6] = [4, 0, 0, 0, 0, 0];
    const B2_X: [u64; 6] = [4, 0, 0, 0, 0, 0];
    const B2_Y: [u64; 6] = [4, 0, 0, 0, 0, 0];
    const X: u64 = 0xd201000000010000; // -0xd201000000010000
    const MAX_BITS: u64 = 384;
    const MODULUS: [u64; 6] = [
        0xb9fe_ffff_ffff_aaab,
        0x1eab_fffe_b153_ffff,
        0x6730_d2a0_f6b0_f624,
        0x6477_4b84_f385_12bf,
        0x4b1b_a7b6_434b_acd7,
        0x1a01_11ea_397f_e69a,
    ];
    /// A nontrivial third root of unity in Fp
    const BETA: [u64; 6] = [
        0x2e01fffffffefffe,
        0xde17d813620a0002,
        0xddb3a93be6f89688,
        0xba69c6076a0f77ea,
        0x5f19672fdf76ce51,
        0x0,
    ];

    const G1_X: [u64; 6] = [
        0xfb3af00adb22c6bb,
        0x6c55e83ff97a1aef,
        0xa14e3a3f171bac58,
        0xc3688c4f9774b905,
        0x2695638c4fa9ac0f,
        0x17f1d3a73197d794,
    ];

    const G1_Y: [u64; 6] = [
        0x0caa232946c5e7e1,
        0xd03cc744a2888ae4,
        0x00db18cb2c04b3ed,
        0xfcf5e095d5d00af6,
        0xa09e30ed741d8ae4,
        0x08b3f481e3aaa0f1,
    ];

    const G2_X0: [u64; 6] = [
        0xd48056c8c121bdb8,
        0x0bac0326a805bbef,
        0xb4510b647ae3d177,
        0xc6e47ad4fa403b02,
        0x260805272dc51051,
        0x024aa2b2f08f0a91,
    ];

    const G2_X1: [u64; 6] = [
        0xe5ac7d055d042b7e,
        0x334cf11213945d57,
        0xb5da61bbdc7f5049,
        0x596bd0d09920b61a,
        0x7dacd3a088274f65,
        0x13e02b6052719f60,
    ];

    const G2_Y0: [u64; 6] = [
        0xe193548608b82801,
        0x923ac9cc3baca289,
        0x6d429a695160d12c,
        0xadfd9baa8cbdd3a7,
        0x8cc9cdc6da2e351a,
        0x0ce5d527727d6e11,
    ];

    const G2_Y1: [u64; 6] = [
        0xaaa9075ff05f79be,
        0x3f370d275cec1da1,
        0x267492ab572e99ab,
        0xcb3e287e85a763af,
        0x32acd2b02bc28b99,
        0x0606c4a02ea734cc,
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
    /// Constant representing the modulus
    /// q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
    const FR_MODULUS: [u64; 4] = [
        0xffff_ffff_0000_0001,
        0x53bd_a402_fffe_5bfe,
        0x3339_d808_09a1_d805,
        0x73ed_a753_299d_7d48,
    ];

    // The number of bits needed to represent the modulus.
    const FR_BITS: u64 = 255;

    // GENERATOR = 7 (multiplicative generator of r-1 order, that is also quadratic nonresidue)
    const FR_GENERATOR: [u64; 4] = [
        // 0x0000_000e_ffff_fff1, 0x17e3_63d3_0018_9c0f,
        // 0xff9c_5787_6f84_57b0,
        // 0x3513_3220_8fc5_a8c4,
        7, 0, 0, 0,
    ];

    /// INV = -(q^{-1} mod 2^64) mod 2^64
    const FR_INV: u64 = 0xffff_fffe_ffff_ffff;

    /// R = 2^256 mod q
    const FR_R: [u64; 4] = [
        0x0000_0001_ffff_fffe,
        0x5884_b7fa_0003_4802,
        0x998c_4fef_ecbc_4ff5,
        0x1824_b159_acc5_056f,
    ];

    /// R^2 = 2^512 mod q
    const FR_R2: [u64; 4] = [
        0xc999_e990_f3f2_9c6d,
        0x2b6c_edcb_8792_5c23,
        0x05d3_1496_7254_398f,
        0x0748_d9d9_9f59_ff11,
    ];

    /// R^3 = 2^768 mod q
    const FR_R3: [u64; 4] = [
        0xc62c_1807_439b_73af,
        0x1b3e_0d18_8cf0_6990,
        0x73d1_3c71_c7b5_f418,
        0x6e2a_5bb9_c8db_33e9,
    ];

    /// 2^-1
    const FR_TWO_INV: [u64; 4] = [
        0x7fffffff80000001,
        0xa9ded2017fff2dff,
        0x199cec0404d0ec02,
        0x39f6d3a994cebea4,
    ];

    // 2^S * t = MODULUS - 1 with t odd
    const FR_S: u64 = 32;

    /// GENERATOR^t where t * 2^s + 1 = q
    /// with t odd. In other words, this
    /// is a 2^s root of unity.
    ///
    /// `GENERATOR = 7 mod q` is a generator
    /// of the q - 1 order multiplicative
    /// subgroup.
    const FR_ROOT_OF_UNITY: [u64; 4] = [
        0x3829971f439f0d2b,
        0xb63683508c2280b9,
        0xd09b681922c813b4,
        0x16a2a19edfe81f20,
    ];

    /// ROOT_OF_UNITY^-1
    const FR_ROOT_OF_UNITY_INV: [u64; 4] = [
        0xFB4D6E13CF19A78,
        0x6F67D4A2B566F833,
        0xED4F2F74A35D0168,
        0x538A6F66E19C653,
    ];

    /// GENERATOR^{2^s} where t * 2^s + 1 = q with t odd.
    /// In other words, this is a t root of unity.
    const FR_DELTA: [u64; 4] = [
        0x6c083479590189d7,
        0xf6502437c6a09c00,
        0x43cab354fabb0062,
        0x8634d0aa021aaf8,
    ];
}
