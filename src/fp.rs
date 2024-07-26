//! This module provides an implementation of the BLS12-381 base field `GF(p)`
//! where `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
use crate::common::{Curve, FieldElement};
use crate::utils::*;
use core::fmt;
use core::mem::transmute;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::RngCore;
use std::marker::PhantomData;

cfg_if::cfg_if! {
    if #[cfg(target_os = "zkvm")] {
        use sp1_zkvm::syscalls::{bls12381_sys_bigint, syscall_bls12381_fp_mulmod};
        use num_bigint::BigUint;
        use sp1_zkvm::lib::{io, unconstrained};
    }
}

// The internal representation of this type is six 64-bit unsigned
// integers in little-endian order. `Fp` values are always in

#[derive(Copy, Clone)]
/// Represents an element in the finite field Fp.
pub struct Fp<C: Curve>(pub(crate) [u64; 6], PhantomData<C>);

impl<C: Curve> fmt::Debug for Fp<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl<C: Curve> Default for Fp<C> {
    fn default() -> Self {
        Fp::<C>::zero()
    }
}

impl<C: Curve> From<u64> for Fp<C> {
    fn from(value: u64) -> Self {
        Fp::<C>::from_raw_unchecked([value, 0, 0, 0, 0, 0])
    }
}

#[cfg(feature = "zeroize")]
impl zeroize::DefaultIsZeroes for Fp {}

impl<C: Curve> Eq for Fp<C> {}
impl<C: Curve> PartialEq for Fp<C> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.iter().zip(other.0.iter()).all(|(a, b)| a == b)
    }
}

impl<'a, C: Curve> Neg for &'a Fp<C> {
    type Output = Fp<C>;

    #[inline]
    fn neg(self) -> Fp<C> {
        self.neg()
    }
}

impl<C: Curve> Neg for Fp<C> {
    type Output = Fp<C>;

    #[inline]
    fn neg(self) -> Fp<C> {
        -&self
    }
}

impl<'a, 'b, C: Curve> Sub<&'b Fp<C>> for &'a Fp<C> {
    type Output = Fp<C>;

    #[inline]
    fn sub(self, rhs: &'b Fp<C>) -> Fp<C> {
        self.sub(rhs)
    }
}

impl<'a, 'b, C: Curve> Add<&'b Fp<C>> for &'a Fp<C> {
    type Output = Fp<C>;

    #[inline]
    fn add(self, rhs: &'b Fp<C>) -> Fp<C> {
        self.add(rhs)
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fp<C>> for &'a Fp<C> {
    type Output = Fp<C>;

    #[inline]
    fn mul(self, rhs: &'b Fp<C>) -> Fp<C> {
        self.mul(rhs)
    }
}

impl<'a, 'b, C: Curve> Div<&'b Fp<C>> for &'a Fp<C> {
    type Output = Fp<C>;

    #[inline]
    fn div(self, rhs: &'b Fp<C>) -> Fp<C> {
        self.div(rhs)
    }
}

cfg_if::cfg_if! {
    if #[cfg(not(target_os = "zkvm"))] {
        impl_binops_multiplicative!(Fp<C>, Fp<C>);
    }
    else {
        impl_binops_multiplicative_mixed!(Fp<C>, Fp<C>, Fp<C>);
        impl<C: Curve> MulAssign<Fp<C>> for Fp<C> {
            #[inline]
            fn mul_assign(&mut self, rhs: Fp<C>) {
                unsafe {
                    let mut lhs = transmute::<[u64; 6], [u32; 12]>(self.0);
                    let rhs = transmute::<[u64; 6], [u32; 12]>(rhs.0);
                    syscall_bls12381_fp_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());

                    *self = Fp::<C>::from_raw_unchecked(transmute::<[u32; 12], [u64; 6]>(lhs));
                }
            }
        }

        impl<'b, C: Curve> MulAssign<&'b Fp<C>> for Fp<C>{
            #[inline]
            fn mul_assign(&mut self, rhs: &'b Fp<C>) {
                *self = &*self * rhs;
            }
        }
    }
}

impl_binops_additive!(Fp<C>, Fp<C>);
impl_binops_divisible!(Fp<C>, Fp<C>);

impl<C: Curve> Fp<C> {
    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Fp<C> {
        Fp::from_raw_unchecked([0, 0, 0, 0, 0, 0])
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Fp<C> {
        Fp::from_raw_unchecked([1, 0, 0, 0, 0, 0])
    }

    /// Checks if this element is zero.
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|&e| e == 0)
    }

    /// Attempts to convert a big-endian byte representation of
    /// a scalar into an `Fp`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 48]) -> Result<Fp<C>, ()> {
        let mut tmp = Fp::from_raw_unchecked([0, 0, 0, 0, 0, 0]);

        tmp.0[5] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[4] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[3] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[2] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());
        tmp.0[1] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap());
        tmp.0[0] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], C::MODULUS[0], 0);
        let (_, borrow) = sbb(tmp.0[1], C::MODULUS[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], C::MODULUS[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], C::MODULUS[3], borrow);
        let (_, borrow) = sbb(tmp.0[4], C::MODULUS[4], borrow);
        let (_, borrow) = sbb(tmp.0[5], C::MODULUS[5], borrow);

        // If the element is smaller than MODULUS then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        if (borrow as u8) & 1 > 0 {
            Err(())
        } else {
            Ok(tmp)
        }
    }

    /// Converts an element of `Fp` into a byte representation in
    /// big-endian byte order.
    pub fn to_bytes(self) -> [u8; 48] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let mut res = [0; 48];
        res[0..8].copy_from_slice(&self.0[5].to_be_bytes());
        res[8..16].copy_from_slice(&self.0[4].to_be_bytes());
        res[16..24].copy_from_slice(&self.0[3].to_be_bytes());
        res[24..32].copy_from_slice(&self.0[2].to_be_bytes());
        res[32..40].copy_from_slice(&self.0[1].to_be_bytes());
        res[40..48].copy_from_slice(&self.0[0].to_be_bytes());

        res
    }

    pub fn from_bytes_unsafe(bytes: &[u8; 48]) -> Fp<C> {
        unsafe { transmute::<[u8; 48], Fp<C>>(*bytes) }
    }

    pub fn to_bytes_unsafe(self) -> [u8; 48] {
        unsafe { transmute::<[u64; 6], [u8; 48]>(self.0) }
    }

    /// Reduces a big-endian 64-bit limb representation of a 768-bit number.
    pub fn from_u768(limbs: [u64; 12]) -> Fp<C> {
        // We reduce an arbitrary 768-bit number by decomposing it into two 384-bit digits
        // with the higher bits multiplied by 2^384. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^384 = R^3

        let d1 = Fp::<C>::from_raw_unchecked([
            limbs[11], limbs[10], limbs[9], limbs[8], limbs[7], limbs[6],
        ]);
        let d0 = Fp::<C>::from_raw_unchecked([
            limbs[5], limbs[4], limbs[3], limbs[2], limbs[1], limbs[0],
        ]);
        d0 + d1 * Fp::<C>::from_raw_unchecked(C::R)
    }

    pub(crate) fn random(mut rng: impl RngCore) -> Fp<C> {
        let mut bytes = [0u8; 96];
        rng.fill_bytes(&mut bytes);

        // Parse the random bytes as a big-endian number, to match Fp encoding order.
        Fp::from_u768([
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[48..56]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[56..64]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[64..72]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[72..80]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[80..88]).unwrap()),
            u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[88..96]).unwrap()),
        ])
    }

    /// Constructs an element of `Fp` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: [u64; 6]) -> Fp<C> {
        Fp(v, PhantomData::<C>)
    }

    /// Although this is labeled "vartime", it is only
    /// variable time with respect to the exponent. It
    /// is also not exposed in the public API.
    pub fn pow_vartime(&self, by: &[u64; 6]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res *= self;
                }
            }
        }
        res
    }

    #[inline]
    /// Computes the square root of this field element.
    pub fn sqrt(&self) -> Result<Self, ()> {
        // We use Shank's method, as p = 3 (mod 4). This means
        // we only need to exponentiate by (p+1)/4. This only
        // works for elements that are actually quadratic residue,
        // so we check that we got the correct result at the end.

        let sqrt = self.pow_vartime(&[
            0xee7f_bfff_ffff_eaab,
            0x07aa_ffff_ac54_ffff,
            0xd9cc_34a8_3dac_3d89,
            0xd91d_d2e1_3ce1_44af,
            0x92c6_e9ed_90d2_eb35,
            0x0680_447a_8e5f_f9a6,
        ]);

        if sqrt.square() == *self {
            return Ok(sqrt);
        } else {
            return Err(());
        }
    }

    #[inline]
    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    #[cfg(not(target_os = "zkvm"))]
    pub fn invert(&self) -> Option<Self> {
        // Exponentiate by p - 2
        let inv = self.pow_vartime(&[
            0xb9fe_ffff_ffff_aaa9,
            0x1eab_fffe_b153_ffff,
            0x6730_d2a0_f6b0_f624,
            0x6477_4b84_f385_12bf,
            0x4b1b_a7b6_434b_acd7,
            0x1a01_11ea_397f_e69a,
        ]);

        Some(inv).filter(|_| !self.is_zero())
    }

    #[cfg(target_os = "zkvm")]
    pub fn invert(&self) -> Option<Self> {
        use sp1_zkvm::io::FD_HINT;

        // Compute the inverse using the zkvm syscall
        unconstrained! {
            let mut buf = [0u8; 48];
            // Exponentiate by p - 2
            let t = self.pow_vartime(&[
                0xb9fe_ffff_ffff_aaa9,
                0x1eab_fffe_b153_ffff,
                0x6730_d2a0_f6b0_f624,
                0x6477_4b84_f385_12bf,
                0x4b1b_a7b6_434b_acd7,
                0x1a01_11ea_397f_e69a,
            ]);
            buf.copy_from_slice(&t.to_bytes_unsafe());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 48] = byte_vec.try_into().unwrap();
        unsafe {
            let inv = Fp::<C>::from_bytes_unsafe(&bytes);
            Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
        }
    }

    #[inline]
    /// Add two field elements together.
    #[cfg(not(target_os = "zkvm"))]
    pub fn add(&self, rhs: &Fp<C>) -> Fp<C> {
        use num_bigint::BigUint;

        unsafe {
            let lhs = BigUint::from_bytes_be(&self.to_bytes());
            let rhs = BigUint::from_bytes_be(&rhs.to_bytes());

            let sum =
                (lhs + rhs) % BigUint::from_slice(&transmute::<[u64; 6], [u32; 12]>(C::MODULUS));

            let mut sum_slice = sum.to_u32_digits();
            sum_slice.resize(12, 0);
            Fp::from_raw_unchecked(transmute::<[u32; 12], [u64; 6]>(
                sum_slice.try_into().unwrap(),
            ))
        }
    }

    #[cfg(target_os = "zkvm")]
    pub fn add(&self, rhs: &Fp<C>) -> Fp<C> {
        let mut result: [u32; 12] = [0; 12];
        unsafe {
            let lhs = transmute::<&[u64; 6], &[u32; 12]>(&self.0);
            let rhs = transmute::<&[u64; 6], &[u32; 12]>(&rhs.0);
            bls12381_sys_bigint(&mut result, 1, lhs, rhs);
            Fp::from_raw_unchecked(*transmute::<&mut [u32; 12], &mut [u64; 6]>(&mut result))
        }
    }

    #[inline]
    /// Returns the negation of this field element.
    pub const fn neg(&self) -> Fp<C> {
        let (d0, borrow) = sbb(C::MODULUS[0], self.0[0], 0);
        let (d1, borrow) = sbb(C::MODULUS[1], self.0[1], borrow);
        let (d2, borrow) = sbb(C::MODULUS[2], self.0[2], borrow);
        let (d3, borrow) = sbb(C::MODULUS[3], self.0[3], borrow);
        let (d4, borrow) = sbb(C::MODULUS[4], self.0[4], borrow);
        let (d5, _) = sbb(C::MODULUS[5], self.0[5], borrow);

        // Let's use a mask if `self` was zero, which would mean
        // the result of the subtraction is p.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3] | self.0[4] | self.0[5]) == 0)
            as u64)
            .wrapping_sub(1);

        Fp::from_raw_unchecked([
            d0 & mask,
            d1 & mask,
            d2 & mask,
            d3 & mask,
            d4 & mask,
            d5 & mask,
        ])
    }

    #[inline]
    /// Squares this element.
    pub fn sub(&self, rhs: &Fp<C>) -> Fp<C> {
        (&rhs.neg()).add(self)
    }

    #[inline]
    /// Multiplies two field elements
    #[cfg(not(target_os = "zkvm"))]
    pub fn mul(&self, rhs: &Fp<C>) -> Fp<C> {
        use num_bigint::BigUint;

        unsafe {
            let slice_lhs = transmute::<&[u64; 6], &[u32; 12]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs);
            let slice_rhs = transmute::<&[u64; 6], &[u32; 12]>(&rhs.0);
            let rhs = BigUint::from_slice(slice_rhs);

            let prod =
                (lhs * rhs) % BigUint::from_slice(&transmute::<[u64; 6], [u32; 12]>(C::MODULUS));

            let mut prod_slice = prod.to_u32_digits();
            prod_slice.resize(12, 0);
            Fp::from_raw_unchecked(transmute::<[u32; 12], [u64; 6]>(
                prod_slice.try_into().unwrap(),
            ))
        }
    }

    /// Multiplies two field elements
    #[cfg(target_os = "zkvm")]
    pub fn mul(&self, rhs: &Fp<C>) -> Fp<C> {
        let mut result: [u32; 12] = [0; 12];
        unsafe {
            let lhs = transmute::<&[u64; 6], &[u32; 12]>(&self.0);
            let rhs = transmute::<&[u64; 6], &[u32; 12]>(&rhs.0);
            bls12381_sys_bigint(&mut result, 0, lhs, rhs);
            Fp::from_raw_unchecked(*transmute::<&mut [u32; 12], &mut [u64; 6]>(&mut result))
        }
    }

    pub fn div(&self, rhs: &Fp<C>) -> Fp<C> {
        self * rhs.invert().unwrap()
    }

    /// Squares this element.
    pub fn square(&self) -> Self {
        self * self
    }
}

impl<C: Curve> FieldElement for Fp<C> {} // For `AffinePoint` trait

#[cfg(test)]
mod test {
    use rand::Rng;

    use crate::common::Bls12381Curve;

    use super::*;

    fn fp_rand() -> Fp<Bls12381Curve> {
        let mut rng = rand::thread_rng();
        Fp::random(&mut rng)
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp::<Bls12381Curve>::from_raw_unchecked(x.clone().try_into().unwrap());
            let b = Fp::<Bls12381Curve>::from_raw_unchecked(x.try_into().unwrap());

            assert_eq!(a, b)
        }
    }

    #[test]
    fn test_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp::<Bls12381Curve>::from_raw_unchecked(x.try_into().unwrap());
            let b = Fp::<Bls12381Curve>::from_raw_unchecked(y.try_into().unwrap());

            assert_ne!(a, b)
        }
    }

    #[test]
    fn test_addition_subtraction() {
        for _ in 0..10 {
            let a = fp_rand();
            let b = fp_rand();
            let c = fp_rand();

            // commutative
            assert_eq!(a + b, b + a);
            assert_eq!(a + (b + c), (a + b) + c);

            // additive identity
            assert_eq!(a + Fp::zero(), a); // a + 0 = a
            assert_eq!(a - Fp::zero(), a); // subtraction identity

            assert_eq!(Fp::zero() - a, -a); // 0 - a = -a
            assert_eq!(a - b, a + (-b)); // a - b = a + -b
            assert_eq!(a - b, a + (b * -Fp::one())); // a - b = a + b * -1

            assert_eq!(-a, Fp::zero() - a);
            assert_eq!(-a, a * -Fp::one());
        }
    }

    #[test]
    fn test_multiplication() {
        for _ in 0..10 {
            let a = fp_rand();
            let b = fp_rand();
            let c = fp_rand();

            // commutative
            assert_eq!(a * b, b * a);

            // associative
            assert_eq!(a * (b * c), (a * b) * c);

            // distributive
            assert_eq!(a * (b + c), a * b + a * c);
        }
    }

    #[test]
    fn test_add_equality() {
        for _ in 0..10 {
            let a = fp_rand();

            assert_eq!(a * Fp::from(0), Fp::zero());
            assert_eq!(a * Fp::zero(), Fp::zero());
            assert_eq!(a * Fp::one(), a);
            assert_eq!(a * Fp::from(1), a);
            assert_eq!(a * Fp::from(2), a + a);
            assert_eq!(a * Fp::from(3), a + a + a);
            assert_eq!(a * Fp::from(4), a + a + a + a);
        }
    }

    #[test]
    fn test_square_equality() {
        for _ in 0..10 {
            let a = fp_rand();
            assert_eq!(a.square(), a * a);
        }
    }

    #[test]
    fn test_pow_equality() {
        for _ in 0..10 {
            let a = fp_rand();
            assert_eq!(a.pow_vartime(&[1, 0, 0, 0, 0, 0]), a);
            assert_eq!(a.pow_vartime(&[2, 0, 0, 0, 0, 0]), a.square());
            assert_eq!(a.pow_vartime(&[3, 0, 0, 0, 0, 0]), a.square() * a);
            assert_eq!(a.pow_vartime(&[4, 0, 0, 0, 0, 0]), a.square().square());
        }
    }

    #[test]
    fn test_sqrt() {
        let sqr1 = Fp::<Bls12381Curve>::from_raw_unchecked([300855555557, 0, 0, 0, 0, 0])
            .sqrt()
            .unwrap();
        assert_eq!(format!("{:?}", sqr1), "0x025e51146a92917731d9d66d63f8c24ed8cae114e7c9d188e3eaa1e79bb19769f5877f9443e03723d9ed1eebbf92df98");

        assert!(
            Fp::<Bls12381Curve>::from_raw_unchecked([72057594037927816, 0, 0, 0, 0, 0])
                .sqrt()
                .is_err()
        );
    }

    #[test]
    fn test_div() {
        for _ in 0..10 {
            let a = fp_rand();

            // division by one
            assert_eq!(a / Fp::one(), a);
            assert_eq!(a / a, Fp::one());

            // division by zero
            assert_eq!(Fp::zero() / a, Fp::zero());

            // division distributivity
            let a = fp_rand();
            let b = fp_rand();
            let c = fp_rand();

            assert_eq!((a + b) / c, a / c + b / c);

            // division and multiplication equality
            let a = fp_rand();
            let b = fp_rand();
            assert_eq!(a / b, a * b.invert().unwrap());
        }
    }
}
