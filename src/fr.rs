//! This module provides an implementation of the BLS12-381 scalar field $\mathbb{F}_q$
//! where `q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::RngCore;
use std::marker::PhantomData;
use std::mem::transmute;

use ff::{Field, PrimeField};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::common::{Bls12381Curve, Curve};
use crate::utils::{adc, sbb};

/// Represents an element of the scalar field $\mathbb{F}_q$ of the BLS12-381 elliptic
/// curve construction.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Scalar` values are always in
#[derive(Clone, Copy)]
pub struct Fr<C: Curve>(pub [u64; 4], PhantomData<C>);

impl<C: Curve> fmt::Debug for Fr<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl<C: Curve> fmt::Display for Fr<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<C: Curve> From<u64> for Fr<C> {
    fn from(val: u64) -> Fr<C> {
        Fr::from_raw([val, 0, 0, 0])
    }
}

impl<C: Curve> ConstantTimeEq for Fr<C> {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
    }
}

impl<C: Curve> PartialEq for Fr<C> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<C: Curve> ConditionallySelectable for Fr<C> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fr::from_raw([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
        ])
    }
}

impl<'a, C: Curve> Neg for &'a Fr<C> {
    type Output = Fr<C>;

    #[inline]
    fn neg(self) -> Fr<C> {
        self.neg()
    }
}

impl<C: Curve> Neg for Fr<C> {
    type Output = Fr<C>;

    #[inline]
    fn neg(self) -> Fr<C> {
        -&self
    }
}

impl<'a, 'b, C: Curve> Sub<&'b Fr<C>> for &'a Fr<C> {
    type Output = Fr<C>;

    #[inline]
    fn sub(self, rhs: &'b Fr<C>) -> Fr<C> {
        self.sub(rhs)
    }
}

impl<'a, 'b, C: Curve> Add<&'b Fr<C>> for &'a Fr<C> {
    type Output = Fr<C>;

    #[inline]
    fn add(self, rhs: &'b Fr<C>) -> Fr<C> {
        self.add(rhs)
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fr<C>> for &'a Fr<C> {
    type Output = Fr<C>;

    #[inline]
    fn mul(self, rhs: &'b Fr<C>) -> Fr<C> {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fr<C>, Fr<C>);
impl_binops_multiplicative!(Fr<C>, Fr<C>);

impl<C: Curve> Default for Fr<C> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<C: Curve> zeroize::DefaultIsZeroes for Fr<C> {}

impl<C: Curve> Fr<C> {
    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Fr<C> {
        Fr::from_raw([0, 0, 0, 0])
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Fr<C> {
        Fr::from_raw([1, 0, 0, 0])
    }

    /// Doubles this field element.
    #[inline]
    pub const fn double(&self) -> Fr<C> {
        // TODO: This can be achieved more efficiently with a bitshift.
        self.add(self)
    }

    pub fn random(mut rng: impl RngCore) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<Fr<C>> {
        let mut tmp = Fr::from_raw([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], C::FR_MODULUS[0], 0);
        let (_, borrow) = sbb(tmp.0[1], C::FR_MODULUS[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], C::FR_MODULUS[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], C::FR_MODULUS[3], borrow);

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
    pub fn from_bytes_wide(bytes: &[u8; 64]) -> Fr<C> {
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

    fn from_u512(limbs: [u64; 8]) -> Fr<C> {
        // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
        // with the higher bits multiplied by 2^256. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^256 = R^3

        let d0 = Fr::from_raw([limbs[0], limbs[1], limbs[2], limbs[3]]);
        let d1 = Fr::from_raw([limbs[4], limbs[5], limbs[6], limbs[7]]);
        d0 + d1 * Fr::from_raw(C::FR_R)
    }

    /// Converts from an integer represented in little endian
    /// into its (congruent) `Scalar` representation.
    pub const fn from_raw(val: [u64; 4]) -> Self {
        Fr(val, PhantomData::<C>)
    }

    /// Squares this element.
    #[inline]
    pub fn square(&self) -> Fr<C> {
        self * self
    }

    /// Exponentiates `self` by `by`, where `by` is a
    /// little-endian order integer exponent.
    pub fn pow(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 0x1) as u8).into());
            }
        }
        res
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
                    res.mul_assign(self);
                }
            }
        }
        res
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(&self) -> CtOption<Self> {
        #[inline(always)]
        fn square_assign_multi<C: Curve>(n: &mut Fr<C>, num_times: usize) {
            for _ in 0..num_times {
                *n = n.square();
            }
        }
        // found using https://github.com/kwantam/addchain
        let mut t0 = self.square();
        let mut t1 = t0 * self;
        let mut t16 = t0.square();
        let mut t6 = t16.square();
        let mut t5 = t6 * t0;
        t0 = t6 * t16;
        let mut t12 = t5 * t16;
        let mut t2 = t6.square();
        let mut t7 = t5 * t6;
        let mut t15 = t0 * t5;
        let mut t17 = t12.square();
        t1 *= t17;
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
        t15 *= &t12;
        t16 *= &t15;
        t3 *= &t16;
        t17 *= &t3;
        t0 *= &t17;
        t6 *= &t0;
        t2 *= &t6;
        square_assign_multi(&mut t0, 8);
        t0 *= &t17;
        square_assign_multi(&mut t0, 9);
        t0 *= &t16;
        square_assign_multi(&mut t0, 9);
        t0 *= &t15;
        square_assign_multi(&mut t0, 9);
        t0 *= &t15;
        square_assign_multi(&mut t0, 7);
        t0 *= &t14;
        square_assign_multi(&mut t0, 7);
        t0 *= &t13;
        square_assign_multi(&mut t0, 10);
        t0 *= &t12;
        square_assign_multi(&mut t0, 9);
        t0 *= &t11;
        square_assign_multi(&mut t0, 8);
        t0 *= &t8;
        square_assign_multi(&mut t0, 8);
        t0 *= self;
        square_assign_multi(&mut t0, 14);
        t0 *= &t9;
        square_assign_multi(&mut t0, 10);
        t0 *= &t8;
        square_assign_multi(&mut t0, 15);
        t0 *= &t7;
        square_assign_multi(&mut t0, 10);
        t0 *= &t6;
        square_assign_multi(&mut t0, 8);
        t0 *= &t5;
        square_assign_multi(&mut t0, 16);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 7);
        t0 *= &t4;
        square_assign_multi(&mut t0, 9);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 5);
        t0 *= &t1;
        square_assign_multi(&mut t0, 5);
        t0 *= &t1;

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Multiplies `rhs` by `self`, returning the result.
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        use num_bigint::BigUint;

        unsafe {
            let modulus = BigUint::from_slice(&transmute::<[u64; 4], [u32; 8]>(C::FR_MODULUS));
            let slice_lhs = transmute::<&[u64; 4], &[u32; 8]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs) % &modulus;
            let rhs = BigUint::from_bytes_le(&rhs.to_bytes()) % &modulus;

            let prod = (lhs * rhs) % modulus;

            let mut prod_slice = prod.to_bytes_le();
            prod_slice.resize(32, 0);
            Fr::from_bytes(&prod_slice.try_into().unwrap()).unwrap()
        }
    }

    /// Subtracts `rhs` from `self`, returning the result.
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sbb(self.0[0], rhs.0[0], 0);
        let (d1, borrow) = sbb(self.0[1], rhs.0[1], borrow);
        let (d2, borrow) = sbb(self.0[2], rhs.0[2], borrow);
        let (d3, borrow) = sbb(self.0[3], rhs.0[3], borrow);

        // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
        // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
        let (d0, carry) = adc(d0, C::FR_MODULUS[0] & borrow, 0);
        let (d1, carry) = adc(d1, C::FR_MODULUS[1] & borrow, carry);
        let (d2, carry) = adc(d2, C::FR_MODULUS[2] & borrow, carry);
        let (d3, _) = adc(d3, C::FR_MODULUS[3] & borrow, carry);

        Fr::from_raw([d0, d1, d2, d3])
    }

    /// Adds `rhs` to `self`, returning the result.
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, carry) = adc(self.0[0], rhs.0[0], 0);
        let (d1, carry) = adc(self.0[1], rhs.0[1], carry);
        let (d2, carry) = adc(self.0[2], rhs.0[2], carry);
        let (d3, _) = adc(self.0[3], rhs.0[3], carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Fr::from_raw([d0, d1, d2, d3])).sub(&Fr::from_raw(C::FR_MODULUS))
    }

    /// Negates `self`.
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sbb(C::FR_MODULUS[0], self.0[0], 0);
        let (d1, borrow) = sbb(C::FR_MODULUS[1], self.0[1], borrow);
        let (d2, borrow) = sbb(C::FR_MODULUS[2], self.0[2], borrow);
        let (d3, _) = sbb(C::FR_MODULUS[3], self.0[3], borrow);

        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Fr::from_raw([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }
}

impl<C: Curve> From<Fr<C>> for [u8; 32] {
    fn from(value: Fr<C>) -> [u8; 32] {
        value.to_bytes()
    }
}

impl<'a, C: Curve> From<&'a Fr<C>> for [u8; 32] {
    fn from(value: &'a Fr<C>) -> [u8; 32] {
        value.to_bytes()
    }
}

impl<C: Curve> Eq for Fr<C> {}
impl<C: Curve + 'static> Field for Fr<C> {
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

impl<C: Curve + 'static> PrimeField for Fr<C> {
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

    const MODULUS: &'static str =
        "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
    const NUM_BITS: u32 = C::FR_BITS as u32;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const TWO_INV: Self = Fr::from_raw(C::FR_TWO_INV);
    const MULTIPLICATIVE_GENERATOR: Self = Fr::from_raw(C::FR_GENERATOR);
    const S: u32 = C::FR_S as u32;
    const ROOT_OF_UNITY: Self = Fr::from_raw(C::FR_ROOT_OF_UNITY);
    const ROOT_OF_UNITY_INV: Self = Fr::from_raw(C::FR_ROOT_OF_UNITY_INV);
    const DELTA: Self = Fr::from_raw(C::FR_DELTA);
}

impl<T, C: Curve> core::iter::Sum<T> for Fr<C>
where
    T: core::borrow::Borrow<Fr<C>>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

impl<T, C: Curve> core::iter::Product<T> for Fr<C>
where
    T: core::borrow::Borrow<Fr<C>>,
{
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::one(), |acc, item| acc * item.borrow())
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use num_bigint::BigUint;

    use crate::common::Bls12381Curve;

    use super::*;

    #[test]
    fn test_constants() {
        assert_eq!(
            Fr::<Bls12381Curve>::MODULUS,
            "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
        );

        assert_eq!(
            Fr::<Bls12381Curve>::from(2) * Fr::<Bls12381Curve>::TWO_INV,
            Fr::<Bls12381Curve>::ONE
        );

        assert_eq!(
            Fr::<Bls12381Curve>::ROOT_OF_UNITY * Fr::<Bls12381Curve>::ROOT_OF_UNITY_INV,
            Fr::<Bls12381Curve>::ONE,
        );

        // ROOT_OF_UNITY^{2^s} mod m == 1
        assert_eq!(
            Fr::<Bls12381Curve>::ROOT_OF_UNITY.pow(&[1u64 << Bls12381Curve::FR_S, 0, 0, 0]),
            Fr::<Bls12381Curve>::ONE,
        );

        // DELTA^{t} mod m == 1
        assert_eq!(
            Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_DELTA).pow(&[
                0xfffe_5bfe_ffff_ffff,
                0x09a1_d805_53bd_a402,
                0x299d_7d48_3339_d808,
                0x0000_0000_73ed_a753,
            ]),
            Fr::<Bls12381Curve>::one(),
        );
    }

    #[test]
    fn test_inv() {
        // Compute -(q^{-1} mod 2^64) mod 2^64 by exponentiating
        // by totient(2**64) - 1

        let mut inv = 1u64;
        for _ in 0..63 {
            inv = inv.wrapping_mul(inv);
            inv = inv.wrapping_mul(Bls12381Curve::FR_MODULUS[0]);
        }
        inv = inv.wrapping_neg();

        assert_eq!(inv, Bls12381Curve::FR_INV);
    }

    #[cfg(feature = "std")]
    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Fr::zero()),
            "0x0000000000000000000000000000000000000000000000000000000000000000"
        );
        assert_eq!(
            format!("{:?}", Fr::one()),
            "0x0000000000000000000000000000000000000000000000000000000000000001"
        );
        assert_eq!(
            format!("{:?}", R2),
            "0x1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffe"
        );
    }

    #[test]
    fn test_equality() {
        assert_eq!(Fr::<Bls12381Curve>::zero(), Fr::<Bls12381Curve>::zero());
        assert_eq!(Fr::<Bls12381Curve>::one(), Fr::<Bls12381Curve>::one());
        assert_eq!(
            Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R2),
            Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R2)
        );

        assert!(Fr::<Bls12381Curve>::zero() != Fr::<Bls12381Curve>::one());
        assert!(Fr::<Bls12381Curve>::one() != Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R2));
    }

    #[test]
    fn test_to_bytes() {
        assert_eq!(
            Fr::<Bls12381Curve>::zero().to_bytes(),
            [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            Fr::<Bls12381Curve>::one().to_bytes(),
            [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            (-&Fr::<Bls12381Curve>::one()).to_bytes(),
            [
                0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ]
        );
    }

    #[test]
    fn test_from_bytes() {
        assert_eq!(
            Fr::<Bls12381Curve>::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Fr::zero()
        );

        assert_eq!(
            Fr::<Bls12381Curve>::from_bytes(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Fr::one()
        );

        assert_eq!(
            Fr(Bls12381Curve::FR_R2, PhantomData::<Bls12381Curve>),
            Fr::from_bytes(&Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R2).to_bytes())
                .unwrap()
        );

        // -1 should work
        assert!(bool::from(
            Fr::<Bls12381Curve>::from_bytes(&[
                0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_some()
        ));

        // modulus is invalid
        assert!(bool::from(
            Fr::<Bls12381Curve>::from_bytes(&[
                1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_none()
        ));

        // Anything larger than the modulus is invalid
        assert!(bool::from(
            Fr::<Bls12381Curve>::from_bytes(&[
                2, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_none()
        ));
        assert!(bool::from(
            Fr::<Bls12381Curve>::from_bytes(&[
                1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 58, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_none()
        ));
        assert!(bool::from(
            Fr::<Bls12381Curve>::from_bytes(&[
                1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 116
            ])
            .is_none()
        ));
    }

    #[test]
    fn test_from_u512_zero() {
        assert_eq!(
            Fr::<Bls12381Curve>::zero(),
            Fr::<Bls12381Curve>::from_u512([
                Bls12381Curve::FR_MODULUS[0],
                Bls12381Curve::FR_MODULUS[1],
                Bls12381Curve::FR_MODULUS[2],
                Bls12381Curve::FR_MODULUS[3],
                0,
                0,
                0,
                0
            ])
        );
    }

    #[test]
    fn test_from_u512_r() {
        assert_eq!(
            Fr::<Bls12381Curve>::from_raw([1, 0, 0, 0]),
            Fr::<Bls12381Curve>::from_u512([1, 0, 0, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_bytes_wide_r() {
        assert_eq!(
            Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R),
            Fr::from_bytes_wide(&[
                254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88, 245, 79, 188, 236,
                239, 79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Fr::<Bls12381Curve>::one(),
            Fr::<Bls12381Curve>::from_bytes_wide(&[
                0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ])
        );
    }

    #[test]
    fn test_zero() {
        assert_eq!(Fr::<Bls12381Curve>::zero(), -&Fr::<Bls12381Curve>::zero());
        assert_eq!(
            Fr::<Bls12381Curve>::zero(),
            Fr::<Bls12381Curve>::zero() + Fr::<Bls12381Curve>::zero()
        );
        assert_eq!(
            Fr::<Bls12381Curve>::zero(),
            Fr::<Bls12381Curve>::zero() - Fr::<Bls12381Curve>::zero()
        );
        assert_eq!(
            Fr::<Bls12381Curve>::zero(),
            Fr::<Bls12381Curve>::zero() * Fr::<Bls12381Curve>::zero()
        );
    }

    #[cfg(test)]
    const LARGEST: Fr<Bls12381Curve> = Fr::<Bls12381Curve>::from_raw([
        0xffff_ffff_0000_0000,
        0x53bd_a402_fffe_5bfe,
        0x3339_d808_09a1_d805,
        0x73ed_a753_299d_7d48,
    ]);

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(
            tmp,
            Fr::<Bls12381Curve>::from_raw([
                0xffff_fffe_ffff_ffff,
                0x53bd_a402_fffe_5bfe,
                0x3339_d808_09a1_d805,
                0x73ed_a753_299d_7d48,
            ])
        );

        let mut tmp = LARGEST;
        tmp += &Fr::<Bls12381Curve>::from_raw([1, 0, 0, 0]);

        assert_eq!(tmp, Fr::zero());
    }

    #[test]
    fn test_negation() {
        let tmp = -&LARGEST;

        assert_eq!(tmp, Fr::<Bls12381Curve>::from_raw([1, 0, 0, 0]));

        let tmp = -&Fr::<Bls12381Curve>::zero();
        assert_eq!(tmp, Fr::zero());
        let tmp = -&Fr::<Bls12381Curve>::from_raw([1, 0, 0, 0]);
        assert_eq!(tmp, LARGEST);
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Fr::zero());

        let mut tmp = Fr::zero();
        tmp -= &LARGEST;

        let mut tmp2 = Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_MODULUS);
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_multiplication() {
        let mut cur = LARGEST;

        for _ in 0..10 {
            let mut tmp = cur;
            tmp *= &cur;

            let mut tmp2 = Fr::zero();
            for b in cur
                .to_bytes()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_squaring() {
        let mut cur = LARGEST;

        for _ in 0..10 {
            let mut tmp = cur;
            tmp = tmp.square();

            let mut tmp2 = Fr::zero();
            for b in cur
                .to_bytes()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_inversion() {
        assert!(bool::from(Fr::<Bls12381Curve>::zero().invert().is_none()));
        assert_eq!(Fr::<Bls12381Curve>::one().invert().unwrap(), Fr::one());
        assert_eq!(
            (-&Fr::<Bls12381Curve>::one()).invert().unwrap(),
            -&Fr::one()
        );

        let mut tmp = Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R2);

        for _ in 0..10 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Fr::one());

            tmp.add_assign(&Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R2));
        }
    }

    #[test]
    fn test_invert_is_pow() {
        let q_minus_2 = [
            0xffff_fffe_ffff_ffff,
            0x53bd_a402_fffe_5bfe,
            0x3339_d808_09a1_d805,
            0x73ed_a753_299d_7d48,
        ];

        let mut r1 = Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R);
        let mut r2 = Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R);
        let mut r3 = Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R);

        for _ in 0..10 {
            r1 = r1.invert().unwrap();
            r2 = r2.pow_vartime(&q_minus_2);
            r3 = r3.pow(&q_minus_2);

            assert_eq!(r1, r2);
            assert_eq!(r2, r3);
            // Add R so we check something different next time around
            r1.add_assign(&Fr::<Bls12381Curve>::from_raw(Bls12381Curve::FR_R));
            r2 = r1;
            r3 = r1;
        }
    }

    #[test]
    fn test_sqrt() {
        {
            assert_eq!(Fr::<Bls12381Curve>::zero().sqrt().unwrap(), Fr::zero());
        }

        let mut none_count = 0;
        for i in 0..100 {
            let square = Fr::<Bls12381Curve>::from_u128(i);
            let square_root = square.sqrt();
            if bool::from(square_root.is_none()) {
                none_count += 1;
            } else {
                assert_eq!(square_root.unwrap() * square_root.unwrap(), square);
            }
        }

        // There are 46 quadratic non-residues mod p = 3 mod 4
        assert_eq!(46, none_count);
    }

    #[test]
    fn test_double() {
        let a = Fr::<Bls12381Curve>::from_raw([
            0x1fff_3231_233f_fffd,
            0x4884_b7fa_0003_4802,
            0x998c_4fef_ecbc_4ff3,
            0x1824_b159_acc5_0562,
        ]);

        assert_eq!(a.double(), a + a);
    }

    #[cfg(feature = "zeroize")]
    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fr::from_raw([
            0x1fff_3231_233f_fffd,
            0x4884_b7fa_0003_4802,
            0x998c_4fef_ecbc_4ff3,
            0x1824_b159_acc5_0562,
        ]);
        a.zeroize();
        assert!(bool::from(a.is_zero()));
    }
}
