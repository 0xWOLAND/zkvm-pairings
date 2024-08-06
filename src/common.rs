use core::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::fp::{Bls12381, Bn254, FpElement};

pub trait AffinePoint: FpElement {
    const X: u64;
    // const BETA: Self;
    const G1_X: Self;
    const G1_Y: Self;
    const G2_X0: Self;
    const G2_X1: Self;
    const G2_Y0: Self;
    const G2_Y1: Self;
}

impl AffinePoint for Bls12381 {
    const X: u64 = 0xd201000000010000; // -0xd201000000010000
    const G1_X: Self = Self::from_raw_unchecked([
        0xfb3af00adb22c6bb,
        0x6c55e83ff97a1aef,
        0xa14e3a3f171bac58,
        0xc3688c4f9774b905,
        0x2695638c4fa9ac0f,
        0x17f1d3a73197d794,
    ]);

    const G1_Y: Self = Self::from_raw_unchecked([
        0x0caa232946c5e7e1,
        0xd03cc744a2888ae4,
        0x00db18cb2c04b3ed,
        0xfcf5e095d5d00af6,
        0xa09e30ed741d8ae4,
        0x08b3f481e3aaa0f1,
    ]);

    const G2_X0: Self = Self::from_raw_unchecked([
        0xd48056c8c121bdb8,
        0x0bac0326a805bbef,
        0xb4510b647ae3d177,
        0xc6e47ad4fa403b02,
        0x260805272dc51051,
        0x024aa2b2f08f0a91,
    ]);

    const G2_X1: Self = Self::from_raw_unchecked([
        0xe5ac7d055d042b7e,
        0x334cf11213945d57,
        0xb5da61bbdc7f5049,
        0x596bd0d09920b61a,
        0x7dacd3a088274f65,
        0x13e02b6052719f60,
    ]);

    const G2_Y0: Self = Self::from_raw_unchecked([
        0xe193548608b82801,
        0x923ac9cc3baca289,
        0x6d429a695160d12c,
        0xadfd9baa8cbdd3a7,
        0x8cc9cdc6da2e351a,
        0x0ce5d527727d6e11,
    ]);

    const G2_Y1: Self = Self::from_raw_unchecked([
        0xaaa9075ff05f79be,
        0x3f370d275cec1da1,
        0x267492ab572e99ab,
        0xcb3e287e85a763af,
        0x32acd2b02bc28b99,
        0x0606c4a02ea734cc,
    ]);
}

impl AffinePoint for Bn254 {
    const X: u64 = 0x7;

    const G1_X: Self = Self::from_raw_unchecked([0x1, 0x0, 0x0, 0x0, 0x0, 0x0]);

    const G1_Y: Self = Self::from_raw_unchecked([0x2, 0x0, 0x0, 0x0, 0x0, 0x0]);

    const G2_X0: Self = Self::from_raw_unchecked([
        0x46debd5cd992f6ed,
        0x674322d4f75edadd,
        0x426a00665e5c4479,
        0x1800deef121f1e76,
    ]);

    const G2_X1: Self = Self::from_raw_unchecked([
        0x97e485b7aef312c2,
        0xf1aa493335a9e712,
        0x7260bfb731fb5d25,
        0x198e9393920d483a,
    ]);

    const G2_Y0: Self = Self::from_raw_unchecked([
        0x4ce6cc0166fa7daa,
        0xe3d1e7690c43d37b,
        0x4aab71808dcb408f,
        0x12c85ea5db8c6deb,
    ]);

    const G2_Y1: Self = Self::from_raw_unchecked([
        0x55acdadcd122975b,
        0xbc4b313370b38ef3,
        0xec9e99ad690c3395,
        0x90689d0585ff075,
    ]);
}
