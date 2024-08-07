//! This module provides an implementation of the BLS12-381 base field `GF(p)`
//! where `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
use crate::utils::*;
use core::fmt;
use core::mem::transmute;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use num_bigint::BigUint;
use rand::RngCore;
use sp1_zkvm::io;
use sp1_zkvm::lib::unconstrained;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::str::FromStr;

cfg_if::cfg_if! {
    if #[cfg(target_os = "zkvm")] {
        use sp1_zkvm::syscalls::{syscall_bls12381_fp_mulmod, syscall_bls12381_fp_addmod, syscall_bls12381_fp_submod};
        use sp1_zkvm::lib::syscall_bn254_fp_addmod;
        use sp1_zkvm::lib::{io, unconstrained};
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Bls12381([u64; 6]);

#[derive(Clone, Copy, Debug)]
pub struct Bn254([u64; 4]);

pub trait FpElement:
    for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + for<'a> Neg<Output = Self>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + PartialEq
    + Copy
    + Clone
    + Sized
    + Debug
    + Sync
    + Send
    + From<u64>
{
    const LIMBS: usize;
    fn zero() -> Self;
    fn one() -> Self;
    fn modulus() -> String;
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
    fn is_one(&self) -> bool {
        self == &Self::one()
    }
    fn invert(&self) -> Option<Self>;
    fn sqrt(&self) -> Option<Self>;
    fn random(rng: impl RngCore) -> Self;
    fn pow_vartime(&self, by: &[u64]) -> Self;
    fn is_lexicographically_largest(&self) -> bool;

    // fn div(&self, rhs: &Self) -> Self {
    //     assert!(!rhs.is_zero(), "Division by zero");
    //     *self * rhs.invert().unwrap()
    // }
    fn square(&self) -> Self {
        *self * *self
    }
}

impl Bls12381 {
    pub const MODULUS: [u64; 6] = [
        // This should be the actual modulus values for BLS12-381
        0xb9fe_ffff_ffff_aaab,
        0x1eab_fffe_b153_ffff,
        0x6730_d2a0_f6b0_f624,
        0x6477_4b84_f385_12bf,
        0x4b1b_a7b6_434b_acd7,
        0x1a01_11ea_397f_e69a,
    ];

    pub(crate) const fn from_raw_unchecked(v: [u64; 6]) -> Self {
        Bls12381(v)
    }

    pub(crate) fn from_bytes(bytes: &[u8; 48]) -> Self {
        let bytes: [u8; 48] = (*bytes).try_into().unwrap();
        unsafe { transmute::<[u8; 48], Bls12381>(bytes) }
    }

    pub(crate) fn to_bytes(&self) -> [u8; 48] {
        unsafe { transmute::<Bls12381, [u8; 48]>(*self) }
    }

    pub(crate) fn _invert(&self) -> Self {
        // Exponentiate by p - 2
        self.pow_vartime(&[
            0xb9fe_ffff_ffff_aaa9,
            0x1eab_fffe_b153_ffff,
            0x6730_d2a0_f6b0_f624,
            0x6477_4b84_f385_12bf,
            0x4b1b_a7b6_434b_acd7,
            0x1a01_11ea_397f_e69a,
        ])
    }

    pub(crate) fn _sqrt(&self) -> Self {
        self.pow_vartime(&[
            0xee7f_bfff_ffff_eaab,
            0x07aa_ffff_ac54_ffff,
            0xd9cc_34a8_3dac_3d89,
            0xd91d_d2e1_3ce1_44af,
            0x92c6_e9ed_90d2_eb35,
            0x0680_447a_8e5f_f9a6,
        ])
    }
}

impl FpElement for Bls12381 {
    const LIMBS: usize = 6;

    fn zero() -> Self {
        Self([0; Self::LIMBS])
    }

    fn one() -> Self {
        Self([1, 0, 0, 0, 0, 0])
    }

    fn modulus() -> String {
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787".to_string()
    }

    #[cfg(not(target_os = "zkvm"))]
    fn invert(&self) -> Option<Self> {
        Some(self._invert()).filter(|_| !self.is_zero())
    }

    #[cfg(target_os = "zkvm")]
    fn invert(&self) -> Option<Self> {
        use sp1_zkvm::io::FD_HINT;

        // Compute the inverse using the zkvm syscall
        unconstrained! {
            let mut buf = [0u8; 48];
            buf.copy_from_slice(&self._invert().to_bytes());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 48] = byte_vec.try_into().unwrap();
        unsafe {
            let inv = Self::from_bytes(&bytes);
            Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
        }
    }

    #[cfg(not(target_os = "zkvm"))]
    fn sqrt(&self) -> Option<Self> {
        Some(self._sqrt()).filter(|s| s.square() == *self)
    }

    #[cfg(target_os = "zkvm")]
    fn sqrt(&self) -> Option<Self> {
        use sp1_zkvm::io::FD_HINT;

        // Compute the square root using the zkvm syscall
        unconstrained! {
            let mut buf = [0u8; 48];
            buf.copy_from_slice(&self._sqrt().to_bytes());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 48] = byte_vec.try_into().unwrap();
        let sqrt = Self::from_bytes(&bytes);
        Some(sqrt).filter(|s| s.square() == *self)
    }

    fn random(mut rng: impl RngCore) -> Self {
        let bytes = Self::MODULUS
            .iter()
            .map(|p| rng.next_u64() % p)
            .collect::<Vec<u64>>();

        Self::from_raw_unchecked(bytes.try_into().unwrap())
    }

    fn pow_vartime(&self, by: &[u64]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res *= *self;
                }
            }
        }
        res
    }

    fn is_lexicographically_largest(&self) -> bool {
        let lhs = self.to_bytes();
        let rhs = (-*self).to_bytes();

        for (l, r) in lhs.iter().zip(rhs.iter()).rev() {
            if l > r {
                return true;
            } else if l < r {
                return false;
            }
        }
        false
    }
}

impl<'a> Add<&'a Bls12381> for Bls12381 {
    type Output = Bls12381;

    #[cfg(not(target_os = "zkvm"))]
    fn add(self, rhs: &'a Self) -> Self {
        use num_bigint::BigUint;

        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let lhs = BigUint::from_bytes_le(&self.to_bytes());
            let rhs = BigUint::from_bytes_le(&rhs.to_bytes());

            let sum = (lhs + rhs) % BigUint::from_str(&<Bls12381 as FpElement>::modulus()).unwrap();

            let mut sum_slice = sum.to_u32_digits();
            sum_slice.resize(2 * LIMBS, 0);
            Self::from_raw_unchecked(transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
                sum_slice.try_into().unwrap(),
            ))
        }
    }

    #[cfg(target_os = "zkvm")]
    fn add(&self, rhs: &'a Self) -> Self {
        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(rhs.0);
            syscall_bls12381_fp_addmod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl Add for Bls12381 {
    type Output = Bls12381;

    fn add(self, rhs: Self) -> Self {
        (&rhs).add(&self)
    }
}

impl<'a> Sub<&'a Bls12381> for Bls12381 {
    type Output = Bls12381;

    #[cfg(not(target_os = "zkvm"))]
    fn sub(self, rhs: &'a Self) -> Self {
        (&rhs).neg().add(&self)
    }

    #[cfg(target_os = "zkvm")]
    fn sub(&self, rhs: &Self) -> Self {
        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(rhs.0);
            syscall_bls12381_fp_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl Sub for Bls12381 {
    type Output = Bls12381;

    fn sub(self, rhs: Self) -> Self {
        (&rhs).neg().add(&self)
    }
}

impl Mul for Bls12381 {
    type Output = Bls12381;

    #[cfg(not(target_os = "zkvm"))]
    fn mul(self, rhs: Self) -> Self {
        use num_bigint::BigUint;

        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let slice_lhs = transmute::<&[u64; LIMBS], &[u32; 2 * LIMBS]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs);
            let slice_rhs = transmute::<&[u64; LIMBS], &[u32; 2 * LIMBS]>(&rhs.0);
            let rhs = BigUint::from_slice(slice_rhs);

            let prod =
                (lhs * rhs) % BigUint::from_str(&<Bls12381 as FpElement>::modulus()).unwrap();

            let mut prod_slice = prod.to_u32_digits();
            prod_slice.resize(2 * LIMBS, 0);
            Self::from_raw_unchecked(transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
                prod_slice.try_into().unwrap(),
            ))
        }
    }

    #[cfg(target_os = "zkvm")]
    fn mul(&self, rhs: Self) -> Self {
        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(rhs.0);
            syscall_bls12381_fp_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl<'a> Mul<&'a Bls12381> for Bls12381 {
    type Output = Bls12381;

    fn mul(self, rhs: &'a Self) -> Self {
        self.mul(*rhs)
    }
}

impl Div for Bls12381 {
    type Output = Bls12381;

    fn div(self, rhs: Self) -> Self {
        rhs.invert().unwrap() * self
    }
}

impl<'a> Div<&'a Bls12381> for Bls12381 {
    type Output = Bls12381;

    fn div(self, rhs: &'a Self) -> Self {
        self.div(*rhs)
    }
}

impl Neg for Bls12381 {
    type Output = Bls12381;

    #[cfg(not(target_os = "zkvm"))]
    fn neg(self) -> Self {
        let (d0, borrow) = sbb(Self::MODULUS[0], self.0[0], 0);
        let (d1, borrow) = sbb(Self::MODULUS[1], self.0[1], borrow);
        let (d2, borrow) = sbb(Self::MODULUS[2], self.0[2], borrow);
        let (d3, borrow) = sbb(Self::MODULUS[3], self.0[3], borrow);
        let (d4, borrow) = sbb(Self::MODULUS[4], self.0[4], borrow);
        let (d5, _) = sbb(Self::MODULUS[5], self.0[5], borrow);

        // Use a mask if `self` was zero, which would mean
        // the result of the subtraction is p.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3] | self.0[4] | self.0[5]) == 0)
            as u64)
            .wrapping_sub(1);

        Self::from_raw_unchecked([
            d0 & mask,
            d1 & mask,
            d2 & mask,
            d3 & mask,
            d4 & mask,
            d5 & mask,
        ])
    }

    #[cfg(target_os = "zkvm")]
    fn neg(&self) -> Self {
        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(Self::MODULUS);
            syscall_bls12381_fp_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl PartialEq for Bls12381 {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl AddAssign for Bls12381 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Bls12381 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for Bls12381 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl DivAssign for Bls12381 {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Bn254 {
    const MODULUS: [u64; 4] = [
        0x3c208c16d87cfd47,
        0x97816a916871ca8d,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ];

    pub(crate) const fn from_raw_unchecked(v: [u64; 4]) -> Self {
        Bn254(v)
    }

    pub(crate) fn from_bytes(bytes: &[u8; 32]) -> Self {
        let bytes: [u8; 32] = (*bytes).try_into().unwrap();
        unsafe { transmute::<[u8; 32], Bn254>(bytes) }
    }

    pub(crate) fn to_bytes(&self) -> [u8; 32] {
        unsafe { transmute::<Bn254, [u8; 32]>(*self) }
    }

    pub(crate) const fn fr_modulus() -> &'static str {
        todo!()
    }

    pub(crate) fn _invert(&self) -> Self {
        self.pow_vartime(&[
            0x3c208c16d87cfd45,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ])
    }

    pub(crate) fn _sqrt(&self) -> Self {
        self.pow_vartime(&[
            0x4f082305b61f3f52,
            0x65e05aa45a1c72a3,
            0x6e14116da0605617,
            0xc19139cb84c680a,
        ])
    }
}

impl FpElement for Bn254 {
    const LIMBS: usize = 4;

    fn zero() -> Self {
        Self([0; Self::LIMBS])
    }

    fn one() -> Self {
        Self([1, 0, 0, 0])
    }

    fn modulus() -> String {
        "21888242871839275222246405745257275088696311157297823662689037894645226208583".to_string()
    }

    #[cfg(not(target_os = "zkvm"))]
    fn sqrt(&self) -> Option<Self> {
        Some(self._sqrt()).filter(|s| s.square() == *self)
    }

    #[cfg(target_os = "zkvm")]
    fn sqrt(&self) -> Option<Self> {
        use sp1_zkvm::io::FD_HINT;

        // Compute the square root using the zkvm syscall
        unconstrained! {
            let mut buf = [0u8; 32];
            buf.copy_from_slice(&self._sqrt().to_bytes());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 32] = byte_vec.try_into().unwrap();
        let sqrt = Self::from_bytes(&bytes);
        Some(sqrt).filter(|s| s.square() == *self)
    }

    #[cfg(not(target_os = "zkvm"))]
    fn invert(&self) -> Option<Self> {
        Some(self._invert()).filter(|_| !self.is_zero())
    }

    #[cfg(target_os = "zkvm")]
    fn invert(&self) -> Option<Self> {
        use sp1_zkvm::io::FD_HINT;

        // Compute the inverse using the zkvm syscall
        unconstrained! {
            let mut buf = [0u8; 32];
            buf.copy_from_slice(&self._invert().to_bytes());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 32] = byte_vec.try_into().unwrap();
        unsafe {
            let inv = Self::from_bytes(&bytes);
            Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
        }
    }

    fn random(mut rng: impl RngCore) -> Self {
        Self::from_raw_unchecked(
            Self::MODULUS
                .iter()
                .map(|p| rng.next_u64() % p)
                .collect::<Vec<u64>>()
                .try_into()
                .unwrap(),
        )
    }

    fn pow_vartime(&self, by: &[u64]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res *= *self;
                }
            }
        }
        res
    }

    fn is_lexicographically_largest(&self) -> bool {
        let lhs = self.to_bytes();
        let rhs = (-*self).to_bytes();

        for (l, r) in lhs.iter().zip(rhs.iter()).rev() {
            if l > r {
                return true;
            } else if l < r {
                return false;
            }
        }
        false
    }
}

impl<'a> Add<&'a Bn254> for Bn254 {
    type Output = Bn254;

    #[cfg(not(target_os = "zkvm"))]
    fn add(self, rhs: &'a Self) -> Self {
        use num_bigint::BigUint;

        const LIMBS: usize = <Bn254 as FpElement>::LIMBS;

        unsafe {
            let lhs = BigUint::from_bytes_le(&self.to_bytes());
            let rhs = BigUint::from_bytes_le(&rhs.to_bytes());

            let sum = (lhs + rhs) % BigUint::from_str(&<Bn254 as FpElement>::modulus()).unwrap();

            let mut sum_slice = sum.to_u32_digits();
            sum_slice.resize(2 * LIMBS, 0);
            Self::from_raw_unchecked(transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
                sum_slice.try_into().unwrap(),
            ))
        }
    }

    #[cfg(target_os = "zkvm")]
    fn add(&self, rhs: &Self) -> Self {
        const LIMBS: usize = <Bls12381 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(rhs.0);
            syscall_bn254_fp_addmod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl Add for Bn254 {
    type Output = Bn254;

    fn add(self, rhs: Self) -> Self {
        (&rhs).add(&self)
    }
}

impl<'a> Sub<&'a Bn254> for Bn254 {
    type Output = Bn254;

    #[cfg(not(target_os = "zkvm"))]
    fn sub(self, rhs: &'a Self) -> Self {
        (&rhs).neg().add(&self)
    }

    #[cfg(target_os = "zkvm")]
    fn sub(&self, rhs: &Self) -> Self {
        const LIMBS: usize = <Bn254 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(rhs.0);
            syscall_bn254_fp_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl Sub for Bn254 {
    type Output = Bn254;

    fn sub(self, rhs: Self) -> Self {
        (&rhs).neg().add(&self)
    }
}

impl Mul for Bn254 {
    type Output = Bn254;

    #[cfg(not(target_os = "zkvm"))]
    fn mul(self, rhs: Self) -> Self {
        use num_bigint::BigUint;

        const LIMBS: usize = <Bn254 as FpElement>::LIMBS;

        unsafe {
            let slice_lhs = transmute::<&[u64; LIMBS], &[u32; 2 * LIMBS]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs);
            let slice_rhs = transmute::<&[u64; LIMBS], &[u32; 2 * LIMBS]>(&rhs.0);
            let rhs = BigUint::from_slice(slice_rhs);

            let prod = (lhs * rhs) % BigUint::from_str(&<Bn254 as FpElement>::modulus()).unwrap();

            let mut prod_slice = prod.to_u32_digits();
            prod_slice.resize(2 * LIMBS, 0);
            Self::from_raw_unchecked(transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
                prod_slice.try_into().unwrap(),
            ))
        }
    }

    #[cfg(target_os = "zkvm")]
    fn mul(&self, rhs: Self) -> Self {
        const LIMBS: usize = <Bn254 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(rhs.0);
            syscall_bn254_fp_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl<'a> Mul<&'a Bn254> for Bn254 {
    type Output = Bn254;

    fn mul(self, rhs: &'a Self) -> Self {
        self.mul(*rhs)
    }
}

impl Div for Bn254 {
    type Output = Bn254;

    fn div(self, rhs: Self) -> Self {
        rhs.invert().unwrap() * self
    }
}

impl<'a> Div<&'a Bn254> for Bn254 {
    type Output = Bn254;

    fn div(self, rhs: &'a Self) -> Self {
        self.div(*rhs)
    }
}

impl Neg for Bn254 {
    type Output = Bn254;

    #[cfg(not(target_os = "zkvm"))]
    fn neg(self) -> Self {
        let (d0, borrow) = sbb(Self::MODULUS[0], self.0[0], 0);
        let (d1, borrow) = sbb(Self::MODULUS[1], self.0[1], borrow);
        let (d2, borrow) = sbb(Self::MODULUS[2], self.0[2], borrow);
        let (d3, _) = sbb(Self::MODULUS[3], self.0[3], borrow);

        // Use a mask if `self` was zero, which would mean
        // the result of the subtraction is p.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Self::from_raw_unchecked([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }

    #[cfg(target_os = "zkvm")]
    fn neg(&self) -> Self {
        const LIMBS: usize = <Bn254 as FpElement>::LIMBS;

        unsafe {
            let mut lhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(self.0);
            let rhs = transmute::<[u64; LIMBS], [u32; 2 * LIMBS]>(Self::MODULUS);
            syscall_bn254_fp_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            Self::from_raw_unchecked(transmute::<&mut [u32; 2 * LIMBS], &mut [u64; LIMBS]>(
                &mut lhs,
            ))
        }
    }
}

impl PartialEq for Bn254 {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl AddAssign for Bn254 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Bn254 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for Bn254 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl DivAssign for Bn254 {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl From<u64> for Bls12381 {
    fn from(n: u64) -> Self {
        Bls12381([n, 0, 0, 0, 0, 0])
    }
}

impl From<u64> for Bn254 {
    fn from(n: u64) -> Self {
        Bn254([n, 0, 0, 0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn bls_rand() -> Bls12381 {
        let mut rng = rand::thread_rng();
        Bls12381::random(&mut rng)
    }

    fn bn_rand() -> Bn254 {
        let mut rng = rand::thread_rng();
        Bn254::random(&mut rng)
    }

    #[test]
    fn test_bls12381_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bls12381::from_raw_unchecked(x.clone().try_into().unwrap());
            let b = Bls12381::from_raw_unchecked(x.try_into().unwrap());
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_bls12381_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bls12381::from_raw_unchecked(x.try_into().unwrap());
            let b = Bls12381::from_raw_unchecked(y.try_into().unwrap());
            assert_ne!(a, b);
        }
    }

    #[test]
    fn test_bn254_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..4).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bn254::from_raw_unchecked(x.clone().try_into().unwrap());
            let b = Bn254::from_raw_unchecked(x.try_into().unwrap());
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_bn254_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..4).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..4).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bn254::from_raw_unchecked(x.try_into().unwrap());
            let b = Bn254::from_raw_unchecked(y.try_into().unwrap());
            assert_ne!(a, b);
        }
    }

    #[test]
    fn test_bls12381_addition_subtraction() {
        for _ in 0..10 {
            let a = bls_rand();
            let b = bls_rand();
            let c = bls_rand();

            assert_eq!(a + b, b + a);
            assert_eq!(a + (b + c), (a + b) + c);
            assert_eq!(a + Bls12381::zero(), a);
            assert_eq!(a - Bls12381::zero(), a);
            assert_eq!(Bls12381::zero() - a, -a);
            assert_eq!(a - b, a + (-b));
            assert_eq!(a - b, a + (b * -Bls12381::one()));
            assert_eq!(-a, Bls12381::zero() - a);
            assert_eq!(-a, a * -Bls12381::one());
        }
    }

    #[test]
    fn test_bn254_addition_subtraction() {
        for _ in 0..10 {
            let a = bn_rand();
            let b = bn_rand();
            let c = bn_rand();

            assert_eq!(a + b, b + a);
            assert_eq!(a + (b + c), (a + b) + c);
            assert_eq!(a + Bn254::zero(), a);
            assert_eq!(a - Bn254::zero(), a);
            assert_eq!(Bn254::zero() - a, -a);
            assert_eq!(a - b, a + (-b));
            assert_eq!(a - b, a + (b * -Bn254::one()));
            assert_eq!(-a, Bn254::zero() - a);
            assert_eq!(-a, a * -Bn254::one());
        }
    }

    #[test]
    fn test_bls12381_multiplication() {
        for _ in 0..10 {
            let a = bls_rand();
            let b = bls_rand();
            let c = bls_rand();

            assert_eq!(a * b, b * a);
            assert_eq!(a * (b * c), (a * b) * c);
            assert_eq!(a * (b + c), a * b + a * c);
        }
    }

    #[test]
    fn test_bn254_multiplication() {
        for _ in 0..10 {
            let a = bn_rand();
            let b = bn_rand();
            let c = bn_rand();

            assert_eq!(a * b, b * a);
            assert_eq!(a * (b * c), (a * b) * c);
            assert_eq!(a * (b + c), a * b + a * c);
        }
    }

    #[test]
    fn test_bls12381_division() {
        for _ in 0..10 {
            let a = bls_rand();

            assert_eq!(a / Bls12381::one(), a);
            assert_eq!(a / a, Bls12381::one());
            assert_eq!(Bls12381::zero() / a, Bls12381::zero());

            let a = bls_rand();
            let b = bls_rand();
            let c = bls_rand();

            assert_eq!((a + b) / c, a / c + b / c);

            let a = bls_rand();
            let b = bls_rand();
            assert_eq!(a / b, a * b.invert().unwrap());
        }
    }

    #[test]
    fn test_bn254_division() {
        for _ in 0..10 {
            let a = bn_rand();

            assert_eq!(a / Bn254::one(), a);
            assert_eq!(a / a, Bn254::one());
            assert_eq!(Bn254::zero() / a, Bn254::zero());

            let a = bn_rand();
            let b = bn_rand();
            let c = bn_rand();

            assert_eq!((a + b) / c, a / c + b / c);

            let a = bn_rand();
            let b = bn_rand();
            assert_eq!(a / b, a * b.invert().unwrap());
        }
    }

    #[test]
    fn test_bls12381_inversion() {
        for _ in 0..10 {
            let a = bls_rand();

            assert_eq!(a * a.invert().unwrap(), Bls12381::one());
            assert_eq!(a.invert().unwrap() * a, Bls12381::one());
            assert_eq!(Bls12381::one().invert().unwrap(), Bls12381::one());
            assert!(Bls12381::zero().invert().is_none());
            let a = bls_rand();
            assert_eq!(a * a.invert().unwrap(), Bls12381::one());
            assert_eq!(a * a._invert(), Bls12381::one());
        }
    }

    #[test]
    fn test_bn254_inversion() {
        for _ in 0..10 {
            let a = bn_rand();

            assert_eq!(a * a.invert().unwrap(), Bn254::one());
            assert_eq!(a.invert().unwrap() * a, Bn254::one());
            assert_eq!(Bn254::one().invert().unwrap(), Bn254::one());
            assert!(Bn254::zero().invert().is_none());
            let a = bn_rand();
            assert_eq!(a * a.invert().unwrap(), Bn254::one());
            assert_eq!(a * a._invert(), Bn254::one());
        }
    }

    #[test]
    fn test_bls12381_sqrt() {
        for _ in 0..10 {
            let a = bls_rand();
            let sqrt = a.sqrt();
            if let Some(s) = sqrt {
                assert_eq!(s.square(), a);
            }
        }
    }

    #[test]
    fn test_bn254_sqrt() {
        for _ in 0..10 {
            let a = bn_rand();
            let sqrt = a.sqrt();
            if let Some(s) = sqrt {
                assert_eq!(s.square(), a);
            }
        }
    }

    #[test]
    fn test_bls12381_random() {
        for _ in 0..100 {
            let a = Bls12381::random(&mut rand::thread_rng());
            let b = Bls12381::random(&mut rand::thread_rng());
            assert_ne!(a, b);
        }
    }

    #[test]
    fn test_bn254_random() {
        for _ in 0..100 {
            let a = Bn254::random(&mut rand::thread_rng());
            let b = Bn254::random(&mut rand::thread_rng());
            assert_ne!(a, b);
        }
    }

    #[test]
    fn test_bytes() {
        for _ in 0..100 {
            let a = bls_rand();
            assert_eq!(a, Bls12381::from_bytes(&a.to_bytes()));
        }
    }
}
