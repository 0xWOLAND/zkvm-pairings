//! This module provides an implementation of the BLS12-381 scalar field $\mathbb{F}_q$
//! where `q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::RngCore;
use std::marker::PhantomData;
use std::mem::transmute;
use std::ops::Div;

use ff::{Field, PrimeField};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::fp::{Bls12381, Bn254, FpElement};
use crate::utils::{adc, sbb};

pub(crate) trait FrElement: FpElement {
    const FR_BITS: u32;
    const FR_MODULUS: [u64; 4];
    const FR_R: [u64; 4];
    const FR_S: u32;
    const FR_TWO_INV: [u64; 4];
    const FR_GENERATOR: [u64; 4];
    const FR_ROOT_OF_UNITY: [u64; 4];
    const FR_ROOT_OF_UNITY_INV: [u64; 4];
    const FR_DELTA: [u64; 4];
    const FR_MODULUS_STR: &'static str;
}

impl FrElement for Bls12381 {
    // The number of bits needed to represent the modulus.
    const FR_BITS: u32 = 255;

    const FR_MODULUS: [u64; 4] = [
        0xffff_ffff_0000_0001,
        0x53bd_a402_fffe_5bfe,
        0x3339_d808_09a1_d805,
        0x73ed_a753_299d_7d48,
    ];

    // GENERATOR = 7 (multiplicative generator of r-1 order, that is also quadratic nonresidue)
    const FR_GENERATOR: [u64; 4] = [
        // 0x0000_000e_ffff_fff1, 0x17e3_63d3_0018_9c0f,
        // 0xff9c_5787_6f84_57b0,
        // 0x3513_3220_8fc5_a8c4,
        7, 0, 0, 0,
    ];

    /// R = 2^256 mod q
    const FR_R: [u64; 4] = [
        0x0000_0001_ffff_fffe,
        0x5884_b7fa_0003_4802,
        0x998c_4fef_ecbc_4ff5,
        0x1824_b159_acc5_056f,
    ];

    /// 2^-1
    const FR_TWO_INV: [u64; 4] = [
        0x7fffffff80000001,
        0xa9ded2017fff2dff,
        0x199cec0404d0ec02,
        0x39f6d3a994cebea4,
    ];

    // 2^S * t = MODULUS - 1 with t odd
    const FR_S: u32 = 32;

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

    const FR_MODULUS_STR: &'static str =
        "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
}

impl FrElement for Bn254 {
    const FR_BITS: u32 = 256;

    const FR_MODULUS: [u64; 4] = [
        0x43e1f593f0000001,
        0x2833e84879b97091,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ];

    const FR_R: [u64; 4] = [
        0xac96341c4ffffffb,
        0x36fc76959f60cd29,
        0x666ea36f7879462e,
        0x0e0a77c19a07df2f,
    ];

    const FR_S: u32 = 28;

    const FR_TWO_INV: [u64; 4] = [
        0xa1f0fac9f8000001,
        0x9419f4243cdcb848,
        0xdc2822db40c0ac2e,
        0x183227397098d014,
    ];

    const FR_GENERATOR: [u64; 4] = [7, 0, 0, 0];

    const FR_ROOT_OF_UNITY: [u64; 4] = [
        0xd34f1ed960c37c9c,
        0x3215cf6dd39329c8,
        0x98865ea93dd31f74,
        0x03ddb9f5166d18b7,
    ];

    const FR_ROOT_OF_UNITY_INV: [u64; 4] = [
        0x0ed3e50a414e6dba,
        0xb22625f59115aba7,
        0x1bbe587180f34361,
        0x048127174daabc26,
    ];

    const FR_DELTA: [u64; 4] = [
        0x870e56bbe533e9a2,
        0x5b5f898e5e963f25,
        0x64ec26aad4c86e71,
        0x09226b6e22c6f0ca,
    ];

    const FR_MODULUS_STR: &'static str = "";
}

/// Represents an element of the scalar field $\mathbb{F}_q$ of the BLS12-381 elliptic
/// curve construction.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Scalar` values are always in
#[derive(Clone, Copy)]
pub struct Fr<F: FrElement>(pub [u64; 4], PhantomData<F>);

impl<F: FrElement> fmt::Debug for Fr<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl<F: FrElement> fmt::Display for Fr<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<F: FrElement> From<u64> for Fr<F> {
    fn from(val: u64) -> Fr<F> {
        Fr::from_raw([val, 0, 0, 0])
    }
}

impl<F: FrElement> ConstantTimeEq for Fr<F> {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
    }
}

impl<F: FrElement> PartialEq for Fr<F> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<F: FrElement> ConditionallySelectable for Fr<F> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fr::from_raw([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
        ])
    }
}

impl<'a, F: FrElement> Neg for &'a Fr<F> {
    type Output = Fr<F>;

    #[inline]
    fn neg(self) -> Fr<F> {
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sbb(F::FR_MODULUS[0], self.0[0], 0);
        let (d1, borrow) = sbb(F::FR_MODULUS[1], self.0[1], borrow);
        let (d2, borrow) = sbb(F::FR_MODULUS[2], self.0[2], borrow);
        let (d3, _) = sbb(F::FR_MODULUS[3], self.0[3], borrow);

        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Fr::from_raw([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }
}

impl<F: FrElement> Neg for Fr<F> {
    type Output = Fr<F>;

    #[inline]
    fn neg(self) -> Fr<F> {
        -&self
    }
}

impl<F: FrElement> Add<Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    #[inline]
    fn add(self, rhs: Fr<F>) -> Fr<F> {
        let (d0, carry) = adc(self.0[0], rhs.0[0], 0);
        let (d1, carry) = adc(self.0[1], rhs.0[1], carry);
        let (d2, carry) = adc(self.0[2], rhs.0[2], carry);
        let (d3, _) = adc(self.0[3], rhs.0[3], carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Fr::from_raw([d0, d1, d2, d3])).sub(Fr::from_raw(F::FR_MODULUS))
    }
}

impl<F: FrElement> Sub<Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    #[inline]
    fn sub(self, rhs: Fr<F>) -> Fr<F> {
        let (d0, borrow) = sbb(self.0[0], rhs.0[0], 0);
        let (d1, borrow) = sbb(self.0[1], rhs.0[1], borrow);
        let (d2, borrow) = sbb(self.0[2], rhs.0[2], borrow);
        let (d3, borrow) = sbb(self.0[3], rhs.0[3], borrow);

        // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
        // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
        let (d0, carry) = adc(d0, F::FR_MODULUS[0] & borrow, 0);
        let (d1, carry) = adc(d1, F::FR_MODULUS[1] & borrow, carry);
        let (d2, carry) = adc(d2, F::FR_MODULUS[2] & borrow, carry);
        let (d3, _) = adc(d3, F::FR_MODULUS[3] & borrow, carry);

        Fr::from_raw([d0, d1, d2, d3])
    }
}

impl<F: FrElement> Mul<Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    /// Multiplies `rhs` by `self`, returning the result.
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        use num_bigint::BigUint;

        unsafe {
            let modulus = BigUint::from_slice(&transmute::<[u64; 4], [u32; 8]>(F::FR_MODULUS));
            let slice_lhs = transmute::<&[u64; 4], &[u32; 8]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs) % &modulus;
            let rhs = BigUint::from_bytes_le(&rhs.to_bytes()) % &modulus;

            let prod = (lhs * rhs) % modulus;

            let mut prod_slice = prod.to_bytes_le();
            prod_slice.resize(32, 0);
            Fr::from_bytes(&prod_slice.try_into().unwrap()).unwrap()
        }
    }
}

impl<F: FrElement> Div<Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    fn div(self, rhs: Self) -> Self {
        self * rhs.invert().unwrap()
    }
}

impl<F: FrElement> AddAssign for Fr<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.add(&rhs);
    }
}

impl<F: FrElement> AddAssign<&Fr<F>> for Fr<F> {
    fn add_assign(&mut self, rhs: &Self) {
        *self = self.add(rhs);
    }
}

impl<F: FrElement> SubAssign for Fr<F> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.sub(rhs);
    }
}

impl<F: FrElement> SubAssign<&Fr<F>> for Fr<F> {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = self.sub(*rhs);
    }
}

impl<F: FrElement> MulAssign for Fr<F> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<F: FrElement> MulAssign<&Fr<F>> for Fr<F> {
    fn mul_assign(&mut self, rhs: &Self) {
        *self = *self * *rhs;
    }
}

impl<'a, F: FrElement> Add<&'a Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    fn add(self, rhs: &'a Fr<F>) -> Fr<F> {
        self.add(*rhs)
    }
}

impl<'a, F: FrElement> Sub<&'a Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    fn sub(self, rhs: &'a Fr<F>) -> Fr<F> {
        self.sub(*rhs)
    }
}

impl<'a, F: FrElement> Mul<&'a Fr<F>> for Fr<F> {
    type Output = Fr<F>;

    fn mul(self, rhs: &'a Fr<F>) -> Fr<F> {
        self * *rhs
    }
}

impl<F: FrElement> Default for Fr<F> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<F: FrElement> zeroize::DefaultIsZeroes for Fr<F> {}

impl<F: FrElement> Fr<F> {
    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Fr<F> {
        Fr::from_raw([0, 0, 0, 0])
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Fr<F> {
        Fr::from_raw([1, 0, 0, 0])
    }

    /// Doubles this field element.
    #[inline]
    pub fn double(&self) -> Fr<F> {
        // TODO: This can be achieved more efficiently with a bitshift.
        *self + *self
    }

    pub fn random(mut rng: impl RngCore) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<Fr<F>> {
        let mut tmp = Fr::from_raw([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], F::FR_MODULUS[0], 0);
        let (_, borrow) = sbb(tmp.0[1], F::FR_MODULUS[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], F::FR_MODULUS[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], F::FR_MODULUS[3], borrow);

        // If the element is smaller than MODULUS then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        CtOption::new(tmp, Choice::from(is_some))
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 32] {
        // Turn into canonical form by computing
        // (a.R) / R = a

        let mut res = [0; 32];
        res[0..8].copy_from_slice(&self.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&self.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&self.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&self.0[3].to_le_bytes());

        res
    }

    /// Converts a 512-bit little endian integer into
    /// a `Scalar` by reducing by the modulus.
    pub fn from_bytes_wide(bytes: &[u8; 64]) -> Fr<F> {
        Fr::from_u512([
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[48..56]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[56..64]).unwrap()),
        ])
    }

    fn from_u512(limbs: [u64; 8]) -> Fr<F> {
        // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
        // with the higher bits multiplied by 2^256. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^256 = R^3

        let d0 = Fr::<F>::from_raw([limbs[0], limbs[1], limbs[2], limbs[3]]);
        let d1 = Fr::<F>::from_raw([limbs[4], limbs[5], limbs[6], limbs[7]]);
        d0 + d1 * Fr::from_raw(F::FR_R)
    }

    /// Converts from an integer represented in little endian
    /// into its (congruent) `Scalar` representation.
    pub const fn from_raw(val: [u64; 4]) -> Self {
        Fr(val, PhantomData::<F>)
    }

    /// Squares this element.
    #[inline]
    pub fn square(&self) -> Fr<F> {
        *self * *self
    }

    /// Exponentiates `self` by `by`, where `by` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub fn pow_vartime(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res = res * *self;
                }
            }
        }
        res
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(&self) -> CtOption<Self> {
        #[inline(always)]
        fn square_assign_multi<F: FrElement>(n: &mut Fr<F>, num_times: usize) {
            for _ in 0..num_times {
                *n = n.square();
            }
        }
        // found using https://github.com/kwantam/addchain
        let mut t0 = self.square();
        let mut t1 = t0 * *self;
        let mut t16 = t0.square();
        let mut t6 = t16.square();
        let mut t5 = t6 * t0;
        t0 = t6 * t16;
        let mut t12 = t5 * t16;
        let mut t2 = t6.square();
        let mut t7 = t5 * t6;
        let mut t15 = t0 * t5;
        let mut t17 = t12.square();
        t1 = t1 * t17;
        let mut t3 = t7 * t2;
        let t8 = t1 * t17;
        let t4 = t8 * t2;
        let t9 = t8 * t7;
        t7 = t4 * t5;
        let t11 = t4 * t17;
        t5 = t9 * t17;
        let t14 = t7 * t15;
        let t13 = t11 * t12;
        t12 = t11 * t17;
        t15 = t15 * t12;
        t16 = t16 * t15;
        t3 = t3 * t16;
        t17 = t17 * t3;
        t0 = t0 * t17;
        t6 = t6 * t0;
        t2 = t2 * t6;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t17;
        square_assign_multi(&mut t0, 9);
        t0 = t0 * t16;
        square_assign_multi(&mut t0, 9);
        t0 = t0 * t15;
        square_assign_multi(&mut t0, 9);
        t0 = t0 * t15;
        square_assign_multi(&mut t0, 7);
        t0 = t0 * t14;
        square_assign_multi(&mut t0, 7);
        t0 = t0 * t13;
        square_assign_multi(&mut t0, 10);
        t0 = t0 * t12;
        square_assign_multi(&mut t0, 9);
        t0 = t0 * t11;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t8;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * *self;
        square_assign_multi(&mut t0, 14);
        t0 = t0 * t9;
        square_assign_multi(&mut t0, 10);
        t0 = t0 * t8;
        square_assign_multi(&mut t0, 15);
        t0 = t0 * t7;
        square_assign_multi(&mut t0, 10);
        t0 = t0 * t6;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t5;
        square_assign_multi(&mut t0, 16);
        t0 = t0 * t3;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 7);
        t0 = t0 * t4;
        square_assign_multi(&mut t0, 9);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t3;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t3;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 8);
        t0 = t0 * t2;
        square_assign_multi(&mut t0, 5);
        t0 = t0 * t1;
        square_assign_multi(&mut t0, 5);
        t0 = t0 * t1;

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }
}

impl<F: FrElement> From<Fr<F>> for [u8; 32] {
    fn from(value: Fr<F>) -> [u8; 32] {
        value.to_bytes()
    }
}

impl<'a, F: FrElement> From<&'a Fr<F>> for [u8; 32] {
    fn from(value: &'a Fr<F>) -> [u8; 32] {
        value.to_bytes()
    }
}

impl<F: FrElement> Eq for Fr<F> {}
impl<F: FrElement + 'static> Field for Fr<F> {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn random(mut rng: impl RngCore) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }

    #[must_use]
    fn square(&self) -> Self {
        self.square()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff::helpers::sqrt_ratio_generic(num, div)
    }

    fn sqrt(&self) -> CtOption<Self> {
        // (t - 1) // 2 = 6104339283789297388802252303364915521546564123189034618274734669823
        ff::helpers::sqrt_tonelli_shanks(
            self,
            &[
                0x7fff_2dff_7fff_ffff,
                0x04d0_ec02_a9de_d201,
                0x94ce_bea4_199c_ec04,
                0x0000_0000_39f6_d3a9,
            ],
        )
    }

    fn is_zero_vartime(&self) -> bool {
        self.0 == Self::zero().0
    }
}

impl<F: FrElement + 'static> PrimeField for Fr<F> {
    type Repr = [u8; 32];

    fn from_repr(r: Self::Repr) -> CtOption<Self> {
        Self::from_bytes(&r)
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes()
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_bytes()[0] & 1)
    }

    const MODULUS: &'static str = F::FR_MODULUS_STR;
    const NUM_BITS: u32 = F::FR_BITS as u32;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const TWO_INV: Self = Fr::from_raw(F::FR_TWO_INV);
    const MULTIPLICATIVE_GENERATOR: Self = Fr::from_raw(F::FR_GENERATOR);
    const S: u32 = F::FR_S;
    const ROOT_OF_UNITY: Self = Fr::from_raw(F::FR_ROOT_OF_UNITY);
    const ROOT_OF_UNITY_INV: Self = Fr::from_raw(F::FR_ROOT_OF_UNITY_INV);
    const DELTA: Self = Fr::from_raw(F::FR_DELTA);
}

impl<T, F: FrElement> core::iter::Sum<T> for Fr<F>
where
    T: core::borrow::Borrow<Fr<F>>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

impl<T, F: FrElement> core::iter::Product<T> for Fr<F>
where
    T: core::borrow::Borrow<Fr<F>>,
{
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::one(), |acc, item| acc * item.borrow())
    }
}

mod tests {
    use crate::fp::Bn254;

    use super::*;

    macro_rules! fr_tests {
        ($curve:ident, $rand_fn:ident, $curve_test: ident) => {
            mod $curve_test {
                use super::*;

                #[test]
                fn test_equality() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = a;
                        assert_eq!(a, b);
                    }
                }

                #[test]
                fn test_inequality() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = $rand_fn();
                        if a != b {
                            assert_ne!(a, b);
                        }
                    }
                }

                #[test]
                fn test_addition_subtraction() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = $rand_fn();
                        let c = $rand_fn();

                        // commutative
                        assert_eq!(a + b, b + a);
                        assert_eq!(a + (b + c), (a + b) + c);

                        // additive identity
                        assert_eq!(a + Fr::<$curve>::zero(), a);
                        assert_eq!(a - Fr::<$curve>::zero(), a);

                        assert_eq!(Fr::<$curve>::zero() - a, -a);
                        assert_eq!(a - b, a + (-b));
                        assert_eq!(a - b, a + (b * -Fr::<$curve>::one()));

                        assert_eq!(-a, Fr::<$curve>::zero() - a);
                        assert_eq!(-a, a * -Fr::<$curve>::one());
                    }
                }

                #[test]
                fn test_multiplication() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = $rand_fn();
                        let c = $rand_fn();

                        // commutative
                        assert_eq!(a * b, b * a);

                        // associative
                        assert_eq!(a * (b * c), (a * b) * c);

                        // distributive
                        assert_eq!(a * (b + c), a * b + a * c);
                    }
                }

                #[test]
                fn test_mul_equality() {
                    for _ in 0..10 {
                        let a = $rand_fn();

                        assert_eq!(a * Fr::<$curve>::zero(), Fr::<$curve>::zero());
                        assert_eq!(a * Fr::<$curve>::one(), a);
                        assert_eq!(a * Fr::<$curve>::from(2u64), a + a);
                        assert_eq!(a * Fr::<$curve>::from(3u64), a + a + a);
                        assert_eq!(a * Fr::<$curve>::from(4u64), a + a + a + a);
                    }
                }

                #[test]
                fn test_square_equality() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        assert_eq!(a.square(), a * a);
                    }
                }

                #[test]
                fn test_div() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = $rand_fn();
                        let c = $rand_fn();

                        // division by one
                        assert_eq!(a / Fr::<$curve>::one(), a);
                        assert_eq!(a / a, Fr::<$curve>::one());

                        // division by zero
                        assert_eq!(Fr::<$curve>::zero() / a, Fr::<$curve>::zero());

                        // division distributivity
                        assert_eq!((a + b) / c, a / c + b / c);

                        // division and multiplication equality
                        if b.is_zero().unwrap_u8() == 0 {
                            assert_eq!(a / b, a * b.invert().unwrap());
                        }
                    }
                }

                #[test]
                fn test_inversion() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        if !a.is_zero().unwrap_u8() == 0 {
                            assert_eq!(a * a.invert().unwrap(), Fr::<$curve>::one());
                            assert_eq!(a.invert().unwrap().invert().unwrap(), a);
                        }
                    }
                }
            }
        };
    }

    fn bls12381_fr_rand() -> Fr<Bls12381> {
        let mut rng = rand::thread_rng();
        Fr::<Bls12381>::random(&mut rng)
    }

    fn bn254_fr_rand() -> Fr<Bn254> {
        let mut rng = rand::thread_rng();
        Fr::<Bn254>::random(&mut rng)
    }

    fr_tests!(Bls12381, bls12381_fr_rand, bls12381_fr_test);
    fr_tests!(Bn254, bn254_fr_rand, bn254_fr_test);
}
