//! This module provides an implementation of the BLS12-381 base field `GF(p)`
//! where `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
use crate::utils::*;
use core::fmt;
use core::mem::transmute;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use num_bigint::BigUint;
use rand::RngCore;
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

// // The internal representation of this type is six 64-bit unsigned
// // integers in little-endian order. `Fp` values are always in

// #[derive(Copy, Clone)]
// /// Represents an element in the finite field Fp.
// pub struct Fp<C: Curve>(pub(crate) [u64; 6], PhantomData<C>);

// impl<C: Curve> fmt::Debug for Fp<C> {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         let tmp = self.to_bytes();
//         write!(f, "0x")?;
//         for &b in tmp.iter() {
//             write!(f, "{:02x}", b)?;
//         }
//         Ok(())
//     }
// }

// impl<C: Curve> Default for Fp<C> {
//     fn default() -> Self {
//         Fp::<C>::zero()
//     }
// }

// impl<C: Curve> From<u64> for Fp<C> {
//     fn from(value: u64) -> Self {
//         Fp::<C>::from_raw_unchecked([value, 0, 0, 0, 0, 0])
//     }
// }

// #[cfg(feature = "zeroize")]
// impl zeroize::DefaultIsZeroes for Fp {}

// impl<C: Curve> Eq for Fp<C> {}
// impl<C: Curve> PartialEq for Fp<C> {
//     #[inline]
//     fn eq(&self, other: &Self) -> bool {
//         self.0.iter().zip(other.0.iter()).all(|(a, b)| a == b)
//     }
// }

// impl<'a, C: Curve> Neg for &'a Fp<C> {
//     type Output = Fp<C>;

//     #[inline]
//     fn neg(self) -> Fp<C> {
//         self.neg()
//     }
// }

// impl<C: Curve> Neg for Fp<C> {
//     type Output = Fp<C>;

//     #[inline]
//     fn neg(self) -> Fp<C> {
//         -&self
//     }
// }

// impl<'a, 'b, C: Curve> Sub<&'b Fp<C>> for &'a Fp<C> {
//     type Output = Fp<C>;

//     #[inline]
//     fn sub(self, rhs: &'b Fp<C>) -> Fp<C> {
//         self.sub(rhs)
//     }
// }

// impl<'a, 'b, C: Curve> Add<&'b Fp<C>> for &'a Fp<C> {
//     type Output = Fp<C>;

//     #[inline]
//     fn add(self, rhs: &'b Fp<C>) -> Fp<C> {
//         self.add(rhs)
//     }
// }

// impl<'a, 'b, C: Curve> Mul<&'b Fp<C>> for &'a Fp<C> {
//     type Output = Fp<C>;

//     #[inline]
//     fn mul(self, rhs: &'b Fp<C>) -> Fp<C> {
//         self.mul(rhs)
//     }
// }

// impl<'a, 'b, C: Curve> Div<&'b Fp<C>> for &'a Fp<C> {
//     type Output = Fp<C>;

//     #[inline]
//     fn div(self, rhs: &'b Fp<C>) -> Fp<C> {
//         self.div(rhs)
//     }
// }

// cfg_if::cfg_if! {
//     if #[cfg(not(target_os = "zkvm"))] {
//         impl_binops_multiplicative!(Fp<C>, Fp<C>);
//     }
//     else {
//         impl_binops_multiplicative_mixed!(Fp<C>, Fp<C>, Fp<C>);
//         impl<C: Curve> MulAssign<Fp<C>> for Fp<C> {
//             #[inline]
//             fn mul_assign(&mut self, rhs: Fp<C>) {
//                 unsafe {
//                     let mut lhs = transmute::<[u64; 6], [u32; 12]>(self.0);
//                     let rhs = transmute::<[u64; 6], [u32; 12]>(rhs.0);
//                     syscall_bls12381_fp_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());

//                     *self = Fp::<C>::from_raw_unchecked(transmute::<[u32; 12], [u64; 6]>(lhs));
//                 }
//             }
//         }

//         impl<'b, C: Curve> MulAssign<&'b Fp<C>> for Fp<C>{
//             #[inline]
//             fn mul_assign(&mut self, rhs: &'b Fp<C>) {
//                 *self = &*self * rhs;
//             }
//         }
//     }
// }

// impl_binops_additive!(Fp<C>, Fp<C>);
// impl_binops_divisible!(Fp<C>, Fp<C>);

// impl<C: Curve> Fp<C> {
//     /// Returns zero, the additive identity.
//     #[inline]
//     pub const fn zero() -> Fp<C> {
//         Fp::from_raw_unchecked([0, 0, 0, 0, 0, 0])
//     }

//     /// Returns one, the multiplicative identity.
//     #[inline]
//     pub const fn one() -> Fp<C> {
//         Fp::from_raw_unchecked([1, 0, 0, 0, 0, 0])
//     }

//     /// Checks if this element is zero.
//     pub fn is_zero(&self) -> bool {
//         self.0.iter().all(|&e| e == 0)
//     }

//     pub fn is_one(&self) -> bool {
//         self == &Fp::one()
//     }

//     /// Attempts to convert a big-endian byte representation of
//     /// a scalar into an `Fp`, failing if the input is not canonical.
//     pub fn from_bytes(bytes: &[u8; 48]) -> Fp<C> {
//         unsafe {
//             let a = BigUint::from_bytes_be(bytes)
//                 % BigUint::from_slice(&transmute::<[u64; 6], [u32; 12]>(C::MODULUS));
//             let mut words = a.to_u64_digits();
//             words.resize(6, 0);

//             Fp::from_raw_unchecked(words.try_into().unwrap())
//         }
//     }

//     pub fn is_lexicographically_largest(&self) -> bool {
//         let lhs = self.0;
//         let rhs = (-self).0;

//         // println!("lexicographically largest");
//         // println!("lhs: {:?}", self);
//         // println!("rhs: {:?}", -self);

//         for (l, r) in lhs.iter().zip(rhs.iter()).rev() {
//             if l > r {
//                 return true;
//             } else if l < r {
//                 return false;
//             }
//         }
//         false
//     }

//     /// Converts an element of `Fp` into a byte representation in
//     /// big-endian byte order.
//     pub fn to_bytes(self) -> [u8; 48] {
//         // Turn into canonical form by computing
//         // (a.R) / R = a
//         let mut res = [0; 48];
//         res[0..8].copy_from_slice(&self.0[5].to_be_bytes());
//         res[8..16].copy_from_slice(&self.0[4].to_be_bytes());
//         res[16..24].copy_from_slice(&self.0[3].to_be_bytes());
//         res[24..32].copy_from_slice(&self.0[2].to_be_bytes());
//         res[32..40].copy_from_slice(&self.0[1].to_be_bytes());
//         res[40..48].copy_from_slice(&self.0[0].to_be_bytes());

//         res
//     }

//     pub fn from_bytes_unsafe(bytes: &[u8; 48]) -> Fp<C> {
//         unsafe { transmute::<[u8; 48], Fp<C>>(*bytes) }
//     }

//     pub fn to_bytes_unsafe(self) -> [u8; 48] {
//         unsafe { transmute::<[u64; 6], [u8; 48]>(self.0) }
//     }

//     pub fn random(mut rng: impl RngCore) -> Fp<C> {
//         Fp::from_raw_unchecked(C::MODULUS.map(|x| rng.next_u64() % x))
//     }

//     /// Constructs an element of `Fp` without checking that it is
//     /// canonical.
//     pub const fn from_raw_unchecked(v: [u64; 6]) -> Fp<C> {
//         Fp(v, PhantomData::<C>)
//     }

//     /// Although this is labeled "vartime", it is only
//     /// variable time with respect to the exponent. It
//     /// is also not exposed in the public API.
//     pub fn pow_vartime(&self, by: &[u64; 6]) -> Self {
//         let mut res = Self::one();
//         for e in by.iter().rev() {
//             for i in (0..64).rev() {
//                 res = res.square();

//                 if ((*e >> i) & 1) == 1 {
//                     res *= self;
//                 }
//             }
//         }
//         res
//     }

//     #[inline]
//     #[cfg(not(target_os = "zkvm"))]
//     /// Computes the square root of this field element.
//     pub fn sqrt(&self) -> Option<Self> {
//         // We use Shank's method, as p = 3 (mod 4). This means
//         // we only need to exponentiate by (p+1)/4. This only
//         // works for elements that are actually quadratic residue,
//         // so we check that we got the correct result at the end.

//         let sqrt = self.pow_vartime(&[
//             0xee7f_bfff_ffff_eaab,
//             0x07aa_ffff_ac54_ffff,
//             0xd9cc_34a8_3dac_3d89,
//             0xd91d_d2e1_3ce1_44af,
//             0x92c6_e9ed_90d2_eb35,
//             0x0680_447a_8e5f_f9a6,
//         ]);

//         Some(sqrt).filter(|s| s.square() == *self)
//     }

//     #[inline]
//     #[cfg(target_os = "zkvm")]
//     pub fn sqrt(&self) -> Option<Self> {
//         use sp1_zkvm::io::FD_HINT;

//         // Compute the square root using the zkvm syscall
//         unconstrained! {
//             let mut buf = [0u8; 48];
//             // Exponentiate by (p + 1) / 4
//             let t = self.pow_vartime(&[
//                 0xee7f_bfff_ffff_eaab,
//                 0x07aa_ffff_ac54_ffff,
//                 0xd9cc_34a8_3dac_3d89,
//                 0xd91d_d2e1_3ce1_44af,
//                 0x92c6_e9ed_90d2_eb35,
//                 0x0680_447a_8e5f_f9a6,
//             ]);
//             buf.copy_from_slice(&t.to_bytes_unsafe());
//             io::write(FD_HINT, &buf);
//         }

//         let byte_vec = io::read_vec();
//         let bytes: [u8; 48] = byte_vec.try_into().unwrap();
//         unsafe {
//             let sqrt = Fp::<C>::from_bytes_unsafe(&bytes);
//             Some(sqrt).filter(|s| s.square() == *self)
//         }
//     }

//     pub(crate) fn _invert(&self) -> Option<Self> {
//         // Exponentiate by p - 2
//         let inv = self.pow_vartime(&[
//             0xb9fe_ffff_ffff_aaa9,
//             0x1eab_fffe_b153_ffff,
//             0x6730_d2a0_f6b0_f624,
//             0x6477_4b84_f385_12bf,
//             0x4b1b_a7b6_434b_acd7,
//             0x1a01_11ea_397f_e69a,
//         ]);

//         Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
//     }

//     #[inline]
//     /// Computes the multiplicative inverse of this field
//     /// element, returning None in the case that this element
//     /// is zero.
//     #[cfg(not(target_os = "zkvm"))]
//     pub fn invert(&self) -> Option<Self> {
//         // Exponentiate by p - 2
//         let inv = self.pow_vartime(&[
//             0xb9fe_ffff_ffff_aaa9,
//             0x1eab_fffe_b153_ffff,
//             0x6730_d2a0_f6b0_f624,
//             0x6477_4b84_f385_12bf,
//             0x4b1b_a7b6_434b_acd7,
//             0x1a01_11ea_397f_e69a,
//         ]);

//         Some(inv).filter(|_| !self.is_zero())
//     }

//     #[cfg(target_os = "zkvm")]
//     pub fn invert(&self) -> Option<Self> {
//         use sp1_zkvm::io::FD_HINT;

//         // Compute the inverse using the zkvm syscall
//         unconstrained! {
//             let mut buf = [0u8; 48];
//             // Exponentiate by p - 2
//             let t = self.pow_vartime(&[
//                 0xb9fe_ffff_ffff_aaa9,
//                 0x1eab_fffe_b153_ffff,
//                 0x6730_d2a0_f6b0_f624,
//                 0x6477_4b84_f385_12bf,
//                 0x4b1b_a7b6_434b_acd7,
//                 0x1a01_11ea_397f_e69a,
//             ]);
//             buf.copy_from_slice(&t.to_bytes_unsafe());
//             io::write(FD_HINT, &buf);
//         }

//         let byte_vec = io::read_vec();
//         let bytes: [u8; 48] = byte_vec.try_into().unwrap();
//         unsafe {
//             let inv = Fp::<C>::from_bytes_unsafe(&bytes);
//             Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
//         }
//     }

//     #[inline]
//     /// Add two field elements together.
//     #[cfg(not(target_os = "zkvm"))]
//     pub fn add(&self, rhs: &Fp<C>) -> Fp<C> {
//         use num_bigint::BigUint;

//         unsafe {
//             let lhs = BigUint::from_bytes_le(&self.to_bytes());
//             let rhs = BigUint::from_bytes_le(&rhs.to_bytes());

//             let sum =
//                 (lhs + rhs) % BigUint::from_slice(&transmute::<[u64; 6], [u32; 12]>(C::MODULUS));

//             let mut sum_slice = sum.to_u32_digits();
//             sum_slice.resize(12, 0);
//             Fp::from_raw_unchecked(transmute::<[u32; 12], [u64; 6]>(
//                 sum_slice.try_into().unwrap(),
//             ))
//         }
//     }

//     #[cfg(target_os = "zkvm")]
//     pub fn add(&self, rhs: &Fp<C>) -> Fp<C> {
//         unsafe {
//             let mut lhs = transmute::<[u64; 6], [u32; 12]>(self.0);
//             let rhs = transmute::<[u64; 6], [u32; 12]>(rhs.0);
//             syscall_bls12381_fp_addmod(lhs.as_mut_ptr(), rhs.as_ptr());
//             Fp::from_raw_unchecked(*transmute::<&mut [u32; 12], &mut [u64; 6]>(&mut lhs))
//         }
//     }

//     #[inline]
//     #[cfg(not(target_os = "zkvm"))]
//     pub fn neg(&self) -> Fp<C> {
//         let (d0, borrow) = sbb(C::MODULUS[0], self.0[0], 0);
//         let (d1, borrow) = sbb(C::MODULUS[1], self.0[1], borrow);
//         let (d2, borrow) = sbb(C::MODULUS[2], self.0[2], borrow);
//         let (d3, borrow) = sbb(C::MODULUS[3], self.0[3], borrow);
//         let (d4, borrow) = sbb(C::MODULUS[4], self.0[4], borrow);
//         let (d5, _) = sbb(C::MODULUS[5], self.0[5], borrow);

//         // Let's use a mask if `self` was zero, which would mean
//         // the result of the subtraction is p.
//         let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3] | self.0[4] | self.0[5]) == 0)
//             as u64)
//             .wrapping_sub(1);

//         Fp::from_raw_unchecked([
//             d0 & mask,
//             d1 & mask,
//             d2 & mask,
//             d3 & mask,
//             d4 & mask,
//             d5 & mask,
//         ])
//     }

//     #[cfg(target_os = "zkvm")]
//     pub fn neg(&self) -> Fp<C> {
//         unsafe {
//             let mut lhs = transmute::<[u64; 6], [u32; 12]>(self.0);
//             let rhs = transmute::<[u64; 6], [u32; 12]>(C::MODULUS);
//             syscall_bls12381_fp_submod(lhs.as_mut_ptr(), rhs.as_ptr());
//             Fp::from_raw_unchecked(*transmute::<&mut [u32; 12], &mut [u64; 6]>(&mut lhs))
//         }
//     }

//     #[inline]
//     #[cfg(not(target_os = "zkvm"))]
//     pub fn sub(&self, rhs: &Fp<C>) -> Fp<C> {
//         (&rhs.neg()).add(self)
//     }

//     #[inline]
//     #[cfg(target_os = "zkvm")]
//     pub fn sub(&self, rhs: &Fp<C>) -> Fp<C> {
//         unsafe {
//             let mut lhs = transmute::<[u64; 6], [u32; 12]>(self.0);
//             let rhs = transmute::<[u64; 6], [u32; 12]>(rhs.0);
//             syscall_bls12381_fp_submod(lhs.as_mut_ptr(), rhs.as_ptr());
//             Fp::from_raw_unchecked(*transmute::<&mut [u32; 12], &mut [u64; 6]>(&mut lhs))
//         }
//     }
//     #[inline]
//     /// Multiplies two field elements
//     #[cfg(not(target_os = "zkvm"))]
//     pub fn mul(&self, rhs: &Fp<C>) -> Fp<C> {
//         use num_bigint::BigUint;

//         unsafe {
//             let slice_lhs = transmute::<&[u64; 6], &[u32; 12]>(&self.0);
//             let lhs = BigUint::from_slice(slice_lhs);
//             let slice_rhs = transmute::<&[u64; 6], &[u32; 12]>(&rhs.0);
//             let rhs = BigUint::from_slice(slice_rhs);

//             let prod =
//                 (lhs * rhs) % BigUint::from_slice(&transmute::<[u64; 6], [u32; 12]>(C::MODULUS));

//             let mut prod_slice = prod.to_u32_digits();
//             prod_slice.resize(12, 0);
//             Fp::from_raw_unchecked(transmute::<[u32; 12], [u64; 6]>(
//                 prod_slice.try_into().unwrap(),
//             ))
//         }
//     }

//     /// Multiplies two field elements
//     #[cfg(target_os = "zkvm")]
//     pub fn mul(&self, rhs: &Fp<C>) -> Fp<C> {
//         unsafe {
//             let mut lhs = transmute::<[u64; 6], [u32; 12]>(self.0);
//             let rhs = transmute::<[u64; 6], [u32; 12]>(rhs.0);
//             syscall_bls12381_fp_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());
//             Fp::from_raw_unchecked(*transmute::<&mut [u32; 12], &mut [u64; 6]>(&mut lhs))
//         }
//     }

//     pub fn div(&self, rhs: &Fp<C>) -> Fp<C> {
//         assert!(!rhs.is_zero(), "Division by zero");
//         self * rhs.invert().unwrap()
//     }

//     /// Squares this element.
//     pub fn square(&self) -> Self {
//         self * self
//     }
// }

// #[cfg(test)]
// mod test {
//     use num_bigint::BigUint;
//     use rand::Rng;

//     use crate::common::Bls12381Curve;

//     use super::*;

//     fn fp_rand() -> Fp<Bls12381Curve> {
//         let mut rng = rand::thread_rng();
//         Fp::random(&mut rng)
//     }

//     #[test]
//     fn test_equality() {
//         let rng = &mut rand::thread_rng();
//         for _ in 0..10 {
//             let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

//             let a = Fp::<Bls12381Curve>::from_raw_unchecked(x.clone().try_into().unwrap());
//             let b = Fp::<Bls12381Curve>::from_raw_unchecked(x.try_into().unwrap());

//             assert_eq!(a, b)
//         }
//     }

//     #[test]
//     fn test_inequality() {
//         let rng = &mut rand::thread_rng();
//         for _ in 0..10 {
//             let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
//             let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

//             let a = Fp::<Bls12381Curve>::from_raw_unchecked(x.try_into().unwrap());
//             let b = Fp::<Bls12381Curve>::from_raw_unchecked(y.try_into().unwrap());

//             assert_ne!(a, b)
//         }
//     }

//     #[test]
//     fn test_addition_subtraction() {
//         for _ in 0..10 {
//             let a = fp_rand();
//             let b = fp_rand();
//             let c = fp_rand();

//             // commutative
//             assert_eq!(a + b, b + a);
//             assert_eq!(a + (b + c), (a + b) + c);

//             // additive identity
//             assert_eq!(a + Fp::zero(), a); // a + 0 = a
//             assert_eq!(a - Fp::zero(), a); // subtraction identity

//             assert_eq!(Fp::zero() - a, -a); // 0 - a = -a
//             assert_eq!(a - b, a + (-b)); // a - b = a + -b
//             assert_eq!(a - b, a + (b * -Fp::one())); // a - b = a + b * -1

//             assert_eq!(-a, Fp::zero() - a);
//             assert_eq!(-a, a * -Fp::one());
//         }
//     }

//     #[test]
//     fn test_multiplication() {
//         for _ in 0..10 {
//             let a = fp_rand();
//             let b = fp_rand();
//             let c = fp_rand();

//             // commutative
//             assert_eq!(a * b, b * a);

//             // associative
//             assert_eq!(a * (b * c), (a * b) * c);

//             // distributive
//             assert_eq!(a * (b + c), a * b + a * c);
//         }
//     }

//     #[test]
//     fn test_add_equality() {
//         for _ in 0..10 {
//             let a = fp_rand();

//             assert_eq!(a * Fp::from(0), Fp::zero());
//             assert_eq!(a * Fp::zero(), Fp::zero());
//             assert_eq!(a * Fp::one(), a);
//             assert_eq!(a * Fp::from(1), a);
//             assert_eq!(a * Fp::from(2), a + a);
//             assert_eq!(a * Fp::from(3), a + a + a);
//             assert_eq!(a * Fp::from(4), a + a + a + a);
//         }
//     }

//     #[test]
//     fn test_square_equality() {
//         for _ in 0..10 {
//             let a = fp_rand();
//             assert_eq!(a.square(), a * a);
//         }
//     }

//     #[test]
//     fn test_pow_equality() {
//         for _ in 0..10 {
//             let a = fp_rand();
//             assert_eq!(a.pow_vartime(&[1, 0, 0, 0, 0, 0]), a);
//             assert_eq!(a.pow_vartime(&[2, 0, 0, 0, 0, 0]), a.square());
//             assert_eq!(a.pow_vartime(&[3, 0, 0, 0, 0, 0]), a.square() * a);
//             assert_eq!(a.pow_vartime(&[4, 0, 0, 0, 0, 0]), a.square().square());
//         }
//     }

//     #[test]
//     fn test_sqrt() {
//         let sqr1 = Fp::<Bls12381Curve>::from_raw_unchecked([300855555557, 0, 0, 0, 0, 0])
//             .sqrt()
//             .unwrap();
//         assert_eq!(format!("{:?}", sqr1), "0x025e51146a92917731d9d66d63f8c24ed8cae114e7c9d188e3eaa1e79bb19769f5877f9443e03723d9ed1eebbf92df98");

//         assert!(
//             Fp::<Bls12381Curve>::from_raw_unchecked([72057594037927816, 0, 0, 0, 0, 0])
//                 .sqrt()
//                 .is_none()
//         );
//     }

//     #[test]
//     fn test_div() {
//         for _ in 0..10 {
//             let a = fp_rand();

//             // division by one
//             assert_eq!(a / Fp::one(), a);
//             assert_eq!(a / a, Fp::one());

//             // division by zero
//             assert_eq!(Fp::zero() / a, Fp::zero());

//             // division distributivity
//             let a = fp_rand();
//             let b = fp_rand();
//             let c = fp_rand();

//             assert_eq!((a + b) / c, a / c + b / c);

//             // division and multiplication equality
//             let a = fp_rand();
//             let b = fp_rand();
//             assert_eq!(a / b, a * b.invert().unwrap());
//         }
//     }

//     #[test]
//     fn test_random() {
//         for _ in 0..100 {
//             let a = Fp::<Bls12381Curve>::random(&mut rand::thread_rng());
//             let b = Fp::<Bls12381Curve>::random(&mut rand::thread_rng());
//             assert_ne!(a, b);

//             assert!(a.0 < Bls12381Curve::MODULUS);
//         }
//     }

//     #[test]
//     fn test_lexicographic_largest() {
//         unsafe {
//             let modulus =
//                 BigUint::from_slice(&transmute::<[u64; 6], [u32; 12]>(Bls12381Curve::MODULUS));
//             for _ in 0..100 {
//                 let mut rng = rand::thread_rng();
//                 let a: Vec<u8> = (0..48).map(|_| rng.gen()).collect();
//                 let a = BigUint::from_bytes_le(a.as_slice()) % &modulus;
//                 let a_inv = &modulus - &a;
//                 let mut a_bytes = a.to_bytes_le();
//                 a_bytes.resize(48, 0);

//                 let a_fp = Fp::<Bls12381Curve>::from_bytes_unsafe(&a_bytes.try_into().unwrap());

//                 assert_eq!(a_fp.is_lexicographically_largest(), a > a_inv);
//             }
//         }
//     }

//     #[test]
//     fn test_inversion() {
//         for _ in 0..10 {
//             let a = fp_rand();

//             // inversion identity
//             assert_eq!(a * a.invert().unwrap(), Fp::one());
//             assert_eq!(a.invert().unwrap() * a, Fp::one());

//             // inversion of one
//             assert_eq!(Fp::<Bls12381Curve>::one().invert().unwrap(), Fp::one());

//             // inversion of zero
//             assert!(Fp::<Bls12381Curve>::zero().invert().is_none());

//             // inversion of a random element
//             let a = fp_rand();
//             assert_eq!(a * a.invert().unwrap(), Fp::one());

//             assert_eq!(a * a._invert().unwrap(), Fp::one());
//         }
//     }
// }

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
{
    const LIMBS: usize;
    const MODULUS: &'static [u64];
    fn zero() -> Self;
    fn one() -> Self;
    fn modulus() -> String;
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
    fn is_one(&self) -> bool {
        self == &Self::one()
    }
    fn _invert(&self) -> Self;
    fn invert(&self) -> Option<Self>;
    fn _sqrt(&self) -> Self;
    fn sqrt(&self) -> Option<Self>;
    fn random(rng: impl RngCore) -> Self;
    fn from_bytes(bytes: &[u8]) -> Self;
    fn to_bytes(&self) -> Vec<u8>;
    fn from_raw_unchecked(v: &[u64]) -> Self;
    fn pow_vartime(&self, by: &[u64]) -> Self;
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
    fn div(&self, rhs: &Self) -> Self {
        assert!(!rhs.is_zero(), "Division by zero");
        *self * rhs.invert().unwrap()
    }
    fn square(&self) -> Self {
        *self * *self
    }
}

impl FpElement for Bls12381 {
    const LIMBS: usize = 6;
    const MODULUS: &'static [u64] = &[
        // This should be the actual modulus values for BLS12-381
        0xb9fe_ffff_ffff_aaab,
        0x1eab_fffe_b153_ffff,
        0x6730_d2a0_f6b0_f624,
        0x6477_4b84_f385_12bf,
        0x4b1b_a7b6_434b_acd7,
        0x1a01_11ea_397f_e69a,
    ];

    fn zero() -> Self {
        Self([0; Self::LIMBS])
    }

    fn one() -> Self {
        Self([1, 0, 0, 0, 0, 0])
    }

    fn modulus() -> String {
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787".to_string()
    }

    fn _invert(&self) -> Self {
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
            buf.copy_from_slice(&self._invert().to_bytes_unsafe());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 48] = byte_vec.try_into().unwrap();
        unsafe {
            let inv = Self::from_bytes_unsafe(&bytes);
            Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
        }
    }

    fn _sqrt(&self) -> Self {
        self.pow_vartime(&[
            0xee7f_bfff_ffff_eaab,
            0x07aa_ffff_ac54_ffff,
            0xd9cc_34a8_3dac_3d89,
            0xd91d_d2e1_3ce1_44af,
            0x92c6_e9ed_90d2_eb35,
            0x0680_447a_8e5f_f9a6,
        ])
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
            buf.copy_from_slice(&self._sqrt().to_bytes_unsafe());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 48] = byte_vec.try_into().unwrap();
        unsafe {
            let sqrt = Self::from_bytes_unsafe(&bytes);
            Some(sqrt).filter(|s| s.square() == *self)
        }
    }

    fn random(mut rng: impl RngCore) -> Self {
        Self::from_raw_unchecked(
            Self::MODULUS
                .iter()
                .map(|p| rng.next_u64() % p)
                .collect::<Vec<u64>>()
                .as_slice(),
        )
    }

    fn from_bytes(bytes: &[u8]) -> Self {
        let bytes: [u8; 48] = bytes.try_into().unwrap();
        unsafe { transmute::<[u8; 48], Bls12381>(bytes) }
    }

    fn to_bytes(&self) -> Vec<u8> {
        unsafe { transmute::<Bls12381, [u8; 48]>(*self) }.to_vec()
    }

    fn from_raw_unchecked(v: &[u64]) -> Self {
        let v: [u64; 6] = v.try_into().unwrap();
        Bls12381(v)
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
            Self::from_raw_unchecked(&transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
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
            Self::from_raw_unchecked(&transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
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

        Self::from_raw_unchecked(&[
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

impl FpElement for Bn254 {
    const LIMBS: usize = 4;
    const MODULUS: &'static [u64] = &[
        0x3c208c16d87cfd47,
        0x97816a916871ca8d,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ];

    fn zero() -> Self {
        Self([0; Self::LIMBS])
    }

    fn one() -> Self {
        Self([1, 0, 0, 0])
    }

    fn modulus() -> String {
        "21888242871839275222246405745257275088696311157297823662689037894645226208583".to_string()
    }

    fn _invert(&self) -> Self {
        self.pow_vartime(&[
            0x3c208c16d87cfd45,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ])
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
            buf.copy_from_slice(&self._invert().to_bytes_unsafe());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 32] = byte_vec.try_into().unwrap();
        unsafe {
            let inv = Self::from_bytes_unsafe(&bytes);
            Some(inv).filter(|_| !self.is_zero() && self * inv == Fp::one())
        }
    }

    fn _sqrt(&self) -> Self {
        self.pow_vartime(&[
            0x4f082305b61f3f52,
            0x65e05aa45a1c72a3,
            0x6e14116da0605617,
            0xc19139cb84c680a,
        ])
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
            buf.copy_from_slice(&self._sqrt().to_bytes_unsafe());
            io::write(FD_HINT, &buf);
        }

        let byte_vec = io::read_vec();
        let bytes: [u8; 32] = byte_vec.try_into().unwrap();
        unsafe {
            let sqrt = Self::from_bytes_unsafe(&bytes);
            Some(sqrt).filter(|s| s.square() == *self)
        }
    }

    fn random(mut rng: impl RngCore) -> Self {
        Self::from_raw_unchecked(
            Self::MODULUS
                .iter()
                .map(|p| rng.next_u64() % p)
                .collect::<Vec<u64>>()
                .as_slice(),
        )
    }

    fn from_bytes(bytes: &[u8]) -> Self {
        let bytes: [u8; 32] = bytes.try_into().unwrap();
        unsafe { transmute::<[u8; 32], Bn254>(bytes) }
    }

    fn to_bytes(&self) -> Vec<u8> {
        unsafe { transmute::<Bn254, [u8; 32]>(*self) }.to_vec()
    }

    fn from_raw_unchecked(v: &[u64]) -> Self {
        let v: [u64; 4] = v.try_into().unwrap();
        Bn254(v)
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
            Self::from_raw_unchecked(&transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
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
            Self::from_raw_unchecked(&transmute::<[u32; 2 * LIMBS], [u64; LIMBS]>(
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

        Self::from_raw_unchecked(&[d0 & mask, d1 & mask, d2 & mask, d3 & mask])
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
            let a = Bls12381::from_raw_unchecked(x.as_slice());
            let b = Bls12381::from_raw_unchecked(x.as_slice());
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_bls12381_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bls12381::from_raw_unchecked(x.as_slice());
            let b = Bls12381::from_raw_unchecked(y.as_slice());
            assert_ne!(a, b);
        }
    }

    #[test]
    fn test_bn254_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..4).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bn254::from_raw_unchecked(x.clone().as_slice());
            let b = Bn254::from_raw_unchecked(x.as_slice());
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_bn254_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..4).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..4).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let a = Bn254::from_raw_unchecked(x.as_slice());
            let b = Bn254::from_raw_unchecked(y.as_slice());
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
