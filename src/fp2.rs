use crate::fp::{Bls12381, Bn254, FpElement};
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};
use rand_core::RngCore;

cfg_if::cfg_if! {
    if #[cfg(target_os = "zkvm")] {
        use sp1_zkvm::syscalls::{syscall_bls12381_fp2_addmod, syscall_bls12381_fp2_submod, syscall_bls12381_fp2_mulmod, syscall_bn254_fp2_addmod, syscall_bn254_fp2_submod, syscall_bn254_fp2_mulmod};
        use std::mem::transmute;
    }
}
pub trait Fp2Element: FpElement {
    fn from_bytes_slice(bytes: &[u8]) -> Fp2<Self>;
    fn to_bytes_vec(f: &Fp2<Self>) -> Vec<u8>;
    fn _invert(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        // (a + bu)(a - bu) = a^2 + b^2
        //
        // which holds because u^2 = -1. This can be rewritten as
        //
        // (a + bu)(a - bu)/(a^2 + b^2) = 1
        //
        // because a^2 + b^2 = 0 has no nonzero solutions for (a, b).
        // This gives that (a - bu)/(a^2 + b^2) is the inverse
        // of (a + bu). Importantly, this can be computing using
        // only a single inversion in F.

        let t = (f.c0.square() + f.c1.square())._invert();
        Some(Fp2::new(f.c0 * t, f.c1 * -t))
    }
    fn invert(f: &Fp2<Self>) -> Option<Fp2<Self>>;
    fn _sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>>;
    fn sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>>;
    fn add(lhs: &Fp2<Self>, rhs: &Fp2<Self>) -> Fp2<Self>;
    fn sub(lhs: &Fp2<Self>, rhs: &Fp2<Self>) -> Fp2<Self>;
    fn mul(lhs: &Fp2<Self>, rhs: &Fp2<Self>) -> Fp2<Self>;
}

impl Fp2Element for Bls12381 {
    fn from_bytes_slice(bytes: &[u8]) -> Fp2<Bls12381> {
        let c0_bytes: [u8; 48] = bytes[..48].try_into().expect("slice with incorrect length");
        let c1_bytes: [u8; 48] = bytes[48..].try_into().expect("slice with incorrect length");

        let c0 = Bls12381::from_bytes_unsafe(&c0_bytes);
        let c1 = Bls12381::from_bytes_unsafe(&c1_bytes);

        Fp2::new(c0, c1)
    }

    fn to_bytes_vec(f: &Fp2<Bls12381>) -> Vec<u8> {
        let mut bytes = [0u8; 96];
        bytes[..48].copy_from_slice(&Self::to_bytes_unsafe(&f.c0));
        bytes[48..].copy_from_slice(&Self::to_bytes_unsafe(&f.c1));
        bytes.to_vec()
    }

    #[cfg(not(target_os = "zkvm"))]
    fn invert(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        Fp2Element::_invert(f)
    }

    #[cfg(target_os = "zkvm")]
    fn invert(f: &Fp2<Bls12381>) -> Option<Fp2<Bls12381>> {
        use sp1_zkvm::{
            io::FD_HINT,
            lib::{io, unconstrained},
        };

        unconstrained! {
            let mut buf = [0u8; 97];
            match Fp2Element::_invert(&f) {
                Some(x) => {
                    buf[96] = 1;
                    buf[0..96].copy_from_slice(&<Bls12381 as Fp2Element>::to_bytes_vec(&x));
                }
                None => {}
            }

            io::write(FD_HINT, &buf);
        }

        let bytes: [u8; 97] = io::read_vec().try_into().unwrap();
        let is_some = bytes[96] == 1;
        let bytes = bytes[..96].try_into().unwrap();
        let out = <Bls12381 as Fp2Element>::from_bytes_slice(bytes);

        Some(out).filter(|_| is_some)
    }

    fn _sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        if f.is_zero() {
            return Some(Fp2::<Bls12381>::zero());
        }

        // a1 = self^((p - 3) / 4)
        let a1 = f.pow_vartime(&[
            0xee7f_bfff_ffff_eaaa,
            0x07aa_ffff_ac54_ffff,
            0xd9cc_34a8_3dac_3d89,
            0xd91d_d2e1_3ce1_44af,
            0x92c6_e9ed_90d2_eb35,
            0x0680_447a_8e5f_f9a6,
        ]);

        // alpha = a1^2 * self = self^((p - 3) / 2 + 1) = self^((p - 1) / 2)
        let alpha = a1.square() * *f;
        // x0 = self^((p + 1) / 4)
        let x0 = a1 * *f;

        if alpha == -Fp2::one() {
            // The element is order p - 1, so we're just trying to get the square of an element of the subfield F.
            // This is given by x0 * u, since u = sqrt(-1). Since the element x0 = a + bu has b = 0, the solution is au.
            Some(Fp2::new(-x0.c1, x0.c0))
        } else {
            // Otherwise, the correct solution is (1 + alpha)^((q - 1) // 2) * x0
            let sqrt = (alpha + Fp2::one()).pow_vartime(&[
                0xdcff_7fff_ffff_d555,
                0x0f55_ffff_58a9_ffff,
                0xb398_6950_7b58_7b12,
                0xb23b_a5c2_79c2_895f,
                0x258d_d3db_21a5_d66b,
                0x0d00_88f5_1cbf_f34d,
            ]) * x0;

            Some(sqrt).filter(|_| sqrt.square() == *f)
        }
    }

    // Computes the square root of this element.
    #[cfg(not(target_os = "zkvm"))]
    fn sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        Fp2Element::_sqrt(f)
    }

    #[cfg(target_os = "zkvm")]
    fn sqrt(f: &Fp2<Bls12381>) -> Option<Fp2<Bls12381>> {
        use sp1_zkvm::{
            io::FD_HINT,
            lib::{io, unconstrained},
        };

        unconstrained! {
            let mut buf = [0u8; 97];
            match Fp2Element::_sqrt(&f) {
                Some(x) => {
                    buf[96] = 1;
                    buf[0..96].copy_from_slice(&<Bls12381 as Fp2Element>::to_bytes_vec(&x));
                }
                None => {}
            }

            io::write(FD_HINT, &buf);
        }

        let byte_vec: [u8; 96] = io::read_vec().try_into().unwrap();
        let is_some = io::read_vec()[0] == 1;
        let out = <Bls12381 as Fp2Element>::from_bytes_slice(&byte_vec);

        Some(out).filter(|_| is_some && out.square() == *f)
    }

    /// Adds another element to this element.
    #[cfg(not(target_os = "zkvm"))]
    fn add(lhs: &Fp2<Bls12381>, rhs: &Fp2<Bls12381>) -> Fp2<Bls12381> {
        Fp2::new((&lhs.c0).add(&rhs.c0), (&lhs.c1).add(&rhs.c1))
    }

    /// Adds another element to this element.
    #[cfg(target_os = "zkvm")]
    fn add(lhs: &Fp2<Bls12381>, rhs: &Fp2<Bls12381>) -> Fp2<Bls12381> {
        unsafe {
            let mut lhs = transmute::<Fp2<Bls12381>, [u32; 24]>(*lhs);
            let rhs = transmute::<Fp2<Bls12381>, [u32; 24]>(*rhs);
            syscall_bls12381_fp2_addmod(lhs.as_mut_ptr(), rhs.as_ptr());
            *transmute::<&mut [u32; 24], &Fp2<Bls12381>>(&mut lhs)
        }
    }

    /// Subtracts another element from this element.
    #[cfg(not(target_os = "zkvm"))]
    fn sub(f: &Fp2<Bls12381>, rhs: &Fp2<Bls12381>) -> Fp2<Bls12381> {
        Fp2::new((&f.c0).sub(&rhs.c0), (&f.c1).sub(&rhs.c1))
    }

    /// Subtracts another element from this element.
    #[cfg(target_os = "zkvm")]
    fn sub(lhs: &Fp2<Bls12381>, rhs: &Fp2<Bls12381>) -> Fp2<Bls12381> {
        unsafe {
            let mut lhs = transmute::<Fp2<Bls12381>, [u32; 24]>(*lhs);
            let rhs = transmute::<Fp2<Bls12381>, [u32; 24]>(*rhs);
            syscall_bls12381_fp2_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            *transmute::<&mut [u32; 24], &Fp2<Bls12381>>(&mut lhs)
        }
    }

    #[cfg(not(target_os = "zkvm"))]
    fn mul(lhs: &Fp2<Bls12381>, rhs: &Fp2<Bls12381>) -> Fp2<Bls12381> {
        // F_{p^2} x F_{p^2} multiplication implemented with operand scanning (schoolbook)
        // computes the result as:
        //
        //   a·b = (a_0 b_0 + a_1 b_1 β) + (a_0 b_1 + a_1 b_0)i
        //
        // In BLS12-381's F_{p^2}, our β is -1, so the resulting F_{p^2} element is:
        //
        //   c_0 = a_0 b_0 - a_1 b_1
        //   c_1 = a_0 b_1 + a_1 b_0
        //
        // Each of these is a "sum of products", which we can compute efficiently.

        Fp2::new(
            lhs.c0 * rhs.c0 - lhs.c1 * rhs.c1,
            lhs.c0 * rhs.c1 + lhs.c1 * rhs.c0,
        )
    }

    /// Multiplies this element by another element.
    #[cfg(target_os = "zkvm")]
    fn mul(lhs: &Fp2<Bls12381>, rhs: &Fp2<Bls12381>) -> Fp2<Bls12381> {
        unsafe {
            let mut lhs = transmute::<Fp2<Bls12381>, [u32; 24]>(*lhs);
            let rhs = transmute::<Fp2<Bls12381>, [u32; 24]>(*rhs);
            syscall_bls12381_fp2_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());
            *transmute::<&mut [u32; 24], &Fp2<Bls12381>>(&mut lhs)
        }
    }
}

impl Fp2Element for Bn254 {
    fn from_bytes_slice(bytes: &[u8]) -> Fp2<Bn254> {
        // Split the bytes array into two 32-byte slices
        let c0_bytes: [u8; 32] = bytes[..32].try_into().expect("slice with incorrect length");
        let c1_bytes: [u8; 32] = bytes[32..].try_into().expect("slice with incorrect length");

        // Call Bn254::from_bytes_unsafe directly
        let c0 = Bn254::from_bytes_unsafe(&c0_bytes);
        let c1 = Bn254::from_bytes_unsafe(&c1_bytes);

        // Create the Fp2 element
        Fp2::new(c0, c1)
    }

    fn to_bytes_vec(f: &Fp2<Bn254>) -> Vec<u8> {
        let mut bytes = [0u8; 64];
        bytes[..32].copy_from_slice(&Self::to_bytes_unsafe(&f.c0));
        bytes[32..].copy_from_slice(&Self::to_bytes_unsafe(&f.c1));
        bytes.to_vec()
    }

    #[cfg(not(target_os = "zkvm"))]
    fn invert(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        Fp2Element::_invert(f)
    }

    #[cfg(target_os = "zkvm")]
    fn invert(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        use sp1_zkvm::{
            io::FD_HINT,
            lib::{io, unconstrained},
        };

        unconstrained! {
            let mut buf = [0u8; 65];
            match Fp2Element::_invert(f) {
                Some(x) => {
                    buf[64] = 1;
                    buf[0..64].copy_from_slice(&<Bn254 as Fp2Element>::to_bytes_vec(&x));
                }
                None => {}
            }

            io::write(FD_HINT, &buf);
        }

        let bytes: [u8; 65] = io::read_vec().try_into().unwrap();
        let is_some = bytes[64] == 1;
        let bytes: [u8; 64] = bytes[..64].try_into().unwrap();
        let out = <Bn254 as Fp2Element>::from_bytes_slice(&bytes);

        Some(out).filter(|_| is_some)
    }

    fn _sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        if f.is_zero() {
            return Some(Fp2::<Bn254>::zero());
        }

        let a1 = f.pow_vartime(&[
            0x4f082305b61f3f51,
            0x65e05aa45a1c72a3,
            0x6e14116da0605617,
            0x0c19139cb84c680a,
        ]);

        let mut alpha = a1.square() * *f;
        let x0 = (a1 * *f).frobenius_map();
        let neg1 = Fp2::new(Bn254::one().neg(), Bn254::zero());

        if x0 == neg1 {
            Some(x0)
        } else {
            let mut a1 = a1 * *f;

            if a1 == neg1 {
                a1 = a1 * Fp2::new(Bn254::zero(), Bn254::one());
            } else {
                alpha = alpha + Fp2::one();
                a1 = a1
                    * alpha.pow_vartime(&[
                        0x9e10460b6c3e7ea3,
                        0xcbc0b548b438e546,
                        0xdc2822db40c0ac2e,
                        0x183227397098d014,
                    ]);
            }
            Some(a1)
        }
    }

    #[cfg(not(target_os = "zkvm"))]
    fn sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        Fp2Element::_sqrt(f)
    }

    #[cfg(target_os = "zkvm")]
    fn sqrt(f: &Fp2<Self>) -> Option<Fp2<Self>> {
        use sp1_zkvm::{
            io::FD_HINT,
            lib::{io, unconstrained},
        };

        unconstrained! {
            let mut buf = [0u8; 97];
            match Fp2Element::_sqrt(&f) {
                Some(x) => {
                    buf[64] = 1;
                    buf[0..64].copy_from_slice(&<Bn254 as Fp2Element>::to_bytes_vec(&x));
                }
                None => {}
            }

            io::write(FD_HINT, &buf);
        }

        let bytes: [u8; 65] = io::read_vec().try_into().unwrap();
        let is_some = bytes[64] == 1;
        let bytes = bytes[..64].try_into().unwrap();
        let out = <Bn254 as Fp2Element>::from_bytes_slice(bytes);

        Some(out).filter(|_| is_some)
    }

    /// Adds another element to this element.
    #[cfg(not(target_os = "zkvm"))]
    fn add(lhs: &Fp2<Bn254>, rhs: &Fp2<Bn254>) -> Fp2<Bn254> {
        Fp2::new((&lhs.c0).add(&rhs.c0), (&lhs.c1).add(&rhs.c1))
    }

    /// Adds another element to this element.
    #[cfg(target_os = "zkvm")]
    fn add(lhs: &Fp2<Bn254>, rhs: &Fp2<Bn254>) -> Fp2<Bn254> {
        unsafe {
            let mut lhs = transmute::<Fp2<Bn254>, [u32; 16]>(*lhs);
            let rhs = transmute::<Fp2<Bn254>, [u32; 16]>(*rhs);
            syscall_bn254_fp2_addmod(lhs.as_mut_ptr(), rhs.as_ptr());
            *transmute::<&mut [u32; 16], &Fp2<Bn254>>(&mut lhs)
        }
    }

    /// Subtracts another element from this element.
    #[cfg(not(target_os = "zkvm"))]
    fn sub(f: &Fp2<Bn254>, rhs: &Fp2<Bn254>) -> Fp2<Bn254> {
        Fp2::new((&f.c0).sub(&rhs.c0), (&f.c1).sub(&rhs.c1))
    }

    /// Subtracts another element from this element.
    #[cfg(target_os = "zkvm")]
    fn sub(lhs: &Fp2<Bn254>, rhs: &Fp2<Bn254>) -> Fp2<Bn254> {
        unsafe {
            let mut lhs = transmute::<Fp2<Bn254>, [u32; 16]>(*lhs);
            let rhs = transmute::<Fp2<Bn254>, [u32; 16]>(*rhs);
            syscall_bls12381_fp2_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            *transmute::<&mut [u32; 16], &Fp2<Bn254>>(&mut lhs)
        }
    }

    #[cfg(not(target_os = "zkvm"))]
    fn mul(lhs: &Fp2<Bn254>, rhs: &Fp2<Bn254>) -> Fp2<Bn254> {
        // F_{p^2} x F_{p^2} multiplication implemented with operand scanning (schoolbook)
        // computes the result as:
        //
        //   a·b = (a_0 b_0 + a_1 b_1 β) + (a_0 b_1 + a_1 b_0)i
        //
        // In BLS12-381's F_{p^2}, our β is -1, so the resulting F_{p^2} element is:
        //
        //   c_0 = a_0 b_0 - a_1 b_1
        //   c_1 = a_0 b_1 + a_1 b_0
        //
        // Each of these is a "sum of products", which we can compute efficiently.

        Fp2::new(
            lhs.c0 * rhs.c0 - lhs.c1 * rhs.c1,
            lhs.c0 * rhs.c1 + lhs.c1 * rhs.c0,
        )
    }

    /// Multiplies this element by another element.
    #[cfg(target_os = "zkvm")]
    fn mul(lhs: &Fp2<Bn254>, rhs: &Fp2<Bn254>) -> Fp2<Bn254> {
        unsafe {
            let mut lhs = transmute::<Fp2<Bn254>, [u32; 16]>(*lhs);
            let rhs = transmute::<Fp2<Bn254>, [u32; 16]>(*rhs);
            syscall_bls12381_fp2_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());
            *transmute::<&mut [u32; 16], &Fp2<Bn254>>(&mut lhs)
        }
    }
}

#[derive(Copy, Clone)]
/// Represents an element in the field Fp2.
pub struct Fp2<F: Fp2Element> {
    /// The first component of the Fp2 element.
    pub c0: F,
    /// The second component of the Fp2 element.
    pub c1: F,
}

impl<F: Fp2Element> fmt::Debug for Fp2<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + {:?}*u", self.c0, self.c1)
    }
}

impl<F: Fp2Element> Default for Fp2<F> {
    fn default() -> Self {
        Fp2::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<F: Fp2Element> zeroize::DefaultIsZeroes for Fp2<F> {}

impl<F: Fp2Element> From<F> for Fp2<F> {
    fn from(f: F) -> Fp2<F> {
        Fp2 { c0: f, c1: f }
    }
}

impl<F: Fp2Element> From<u64> for Fp2<F> {
    fn from(f: u64) -> Fp2<F> {
        Fp2 {
            c0: F::from(f),
            c1: F::zero(),
        }
    }
}

impl<F: Fp2Element> Eq for Fp2<F> {}
impl<F: Fp2Element> PartialEq for Fp2<F> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.c0.eq(&other.c0) & self.c1.eq(&other.c1)
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

impl<'a, F: Fp2Element> Neg for &'a Fp2<F> {
    type Output = Fp2<F>;

    #[inline]
    fn neg(self) -> Fp2<F> {
        self.neg()
    }
}

impl<F: Fp2Element> Neg for Fp2<F> {
    type Output = Fp2<F>;

    #[inline]
    fn neg(self) -> Fp2<F> {
        -&self
    }
}

impl<F: Fp2Element> Sub<Fp2<F>> for Fp2<F> {
    type Output = Fp2<F>;

    /// Subtracts another element from this element.
    fn sub(self, rhs: Fp2<F>) -> Fp2<F> {
        Fp2Element::sub(&self, &rhs)
    }
}

impl<F: Fp2Element> Add<Fp2<F>> for Fp2<F> {
    type Output = Fp2<F>;

    /// Adds another element to this element.
    fn add(self, rhs: Fp2<F>) -> Fp2<F> {
        Fp2Element::add(&self, &rhs)
    }
}

impl<F: Fp2Element> Mul<Fp2<F>> for Fp2<F> {
    type Output = Fp2<F>;

    fn mul(self, rhs: Fp2<F>) -> Fp2<F> {
        Fp2Element::mul(&self, &rhs)
    }
}

impl<F: Fp2Element> Mul<F> for Fp2<F> {
    type Output = Fp2<F>;

    #[inline]
    fn mul(self, rhs: F) -> Fp2<F> {
        Fp2::new(self.c0 * rhs, self.c1 * rhs)
    }
}

impl<F: Fp2Element> Div<Fp2<F>> for Fp2<F> {
    type Output = Fp2<F>;

    #[inline]
    fn div(self, rhs: Fp2<F>) -> Fp2<F> {
        self * rhs.invert().unwrap()
    }
}

impl<F: Fp2Element> Fp2<F> {
    /// Returns the zero element of Fp2.
    #[inline]
    pub fn zero() -> Fp2<F> {
        Fp2::new(F::zero(), F::zero())
    }

    /// Returns the one element of Fp2.
    #[inline]
    pub fn one() -> Fp2<F> {
        Fp2::new(F::one(), F::zero())
    }

    pub const fn new(c0: F, c1: F) -> Fp2<F> {
        Fp2 { c0, c1 }
    }

    pub fn non_residue() -> Fp2<F> {
        Fp2::new(F::one(), F::one())
    }

    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    /// Checks if this element is zero.
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Generates a random element in Fp2.
    pub fn random(mut rng: impl RngCore) -> Fp2<F> {
        Fp2::new(F::random(&mut rng), F::random(&mut rng))
    }

    /// Raises this element to p.
    #[inline(always)]
    pub fn frobenius_map(&self) -> Self {
        // This is always just a conjugation. If you're curious why, here's
        // an article about it: https://alicebob.cryptoland.net/the-frobenius-endomorphism-with-finite-fields/
        self.conjugate()
    }

    /// Computes the conjugate of this element.
    #[inline(always)]
    pub fn conjugate(&self) -> Self {
        Fp2::new(self.c0, -self.c1)
    }

    /// Multiplies this element by the non-residue.
    #[inline(always)]
    pub fn mul_by_nonresidue(&self) -> Fp2<F> {
        // Multiply a + bu by u + 1, getting
        // au + a + bu^2 + bu
        // and because u^2 = -1, we get
        // (a - b) + (a + b)u

        // Fp2::new(self.c0 - self.c1, self.c0 + self.c1)
        *self * Fp2::non_residue()
    }

    /// Computes the square of this element.
    pub fn square(&self) -> Fp2<F> {
        *self * self.clone()
    }
    /// Negates this element.
    pub fn neg(&self) -> Fp2<F> {
        Fp2::zero() - *self
    }

    #[inline]
    pub fn lexicographically_largest(&self) -> bool {
        self.c1.is_lexicographically_largest()
            || (self.c1 == F::zero() && self.c0.is_lexicographically_largest())
    }

    pub(crate) fn invert(&self) -> Option<Self> {
        Fp2Element::invert(self)
    }

    pub(crate) fn sqrt(&self) -> Option<Fp2<F>> {
        Fp2Element::sqrt(self)
    }

    /// Although this is labeled "vartime", it is only
    /// variable time with respect to the exponent. It
    /// is also not exposed in the public API.
    pub fn pow_vartime(&self, by: &[u64]) -> Self {
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
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fp::{Bls12381, Bn254};
    use num_bigint::BigUint;
    use rand::Rng;
    use std::str::FromStr;

    fn bls12381_fp2_rand() -> Fp2<Bls12381> {
        let mut rng = rand::thread_rng();
        Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng))
    }

    fn bn254_fp2_rand() -> Fp2<Bn254> {
        let mut rng = rand::thread_rng();
        Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng))
    }

    macro_rules! fp2_tests {
        ($curve:ident, $rand_fn:ident, $curve_test: ident) => {
            mod $curve_test {
                use super::*;

                #[test]
                fn test_equality() {
                    let rng = &mut rand::thread_rng();
                    for _ in 0..10 {
                        let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
                        let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
                        let a = Fp2::<$curve>::new(
                            $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                            $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                        );
                        let b = Fp2::<$curve>::new(
                            $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                            $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                        );
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
                        assert_eq!(a + Fp2::<$curve>::zero(), a);
                        assert_eq!(a - Fp2::<$curve>::zero(), a);

                        assert_eq!(Fp2::<$curve>::zero() - a, -a);
                        assert_eq!(a - b, a + (-b));
                        assert_eq!(a - b, a + (b * -Fp2::<$curve>::one()));

                        assert_eq!(-a, Fp2::<$curve>::zero() - a);
                        assert_eq!(-a, a * -Fp2::<$curve>::one());
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
                fn test_add_equality() {
                    for _ in 0..10 {
                        let a = $rand_fn();

                        assert_eq!(a * $curve::zero(), Fp2::<$curve>::zero());
                        assert_eq!(a * Fp2::<$curve>::zero(), Fp2::<$curve>::zero());
                        assert_eq!(a * Fp2::<$curve>::one(), a);
                        assert_eq!(a * $curve::one(), a);
                        assert_eq!(a * $curve::from(2u64), a + a);
                        assert_eq!(a * $curve::from(3u64), a + a + a);
                        assert_eq!(a * $curve::from(4u64), a + a + a + a);
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
                fn test_pow_equality() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        assert_eq!(a.pow_vartime(&[1, 0, 0, 0, 0, 0]), a);
                        assert_eq!(a.pow_vartime(&[2, 0, 0, 0, 0, 0]), a.square());
                        assert_eq!(a.pow_vartime(&[3, 0, 0, 0, 0, 0]), a.square() * a);
                        assert_eq!(a.pow_vartime(&[4, 0, 0, 0, 0, 0]), a.square().square());
                    }
                }

                #[test]
                fn test_sqrt() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let a_sq = a.square();
                        let a_sqrt = a_sq.sqrt();
                        if a_sqrt.is_some().into() {
                            assert_eq!(a_sqrt.unwrap().square(), a_sq);
                        }
                    }
                }

                #[test]
                fn test_div() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = $rand_fn();
                        let c = $rand_fn();

                        // division by one
                        assert_eq!(a / Fp2::<$curve>::one(), a);
                        assert_eq!(a / a, Fp2::<$curve>::one());

                        // division by zero
                        assert_eq!(Fp2::<$curve>::zero() / a, Fp2::<$curve>::zero());

                        // division distributivity
                        assert_eq!((a + b) / c, a / c + b / c);

                        // division and multiplication equality
                        if !b.is_zero() {
                            assert_eq!(a / b, a * b.invert().unwrap());
                        }
                    }
                }

                #[test]
                fn test_inversion() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        if !a.is_zero() {
                            assert_eq!(a * a.invert().unwrap(), Fp2::<$curve>::one());
                            assert_eq!(a.invert().unwrap().invert().unwrap(), a);
                        }
                    }
                }

                #[test]
                fn test_lexicographic_largest() {
                    unsafe {
                        let modulus = BigUint::from_str($curve::modulus().as_str()).unwrap();

                        let gen_test_value = || {
                            let mut rng = rand::thread_rng();
                            let a: Vec<u8> = (0..48).map(|_| rng.gen()).collect();
                            let a = BigUint::from_bytes_le(a.as_slice()) % &modulus;
                            let a_inv = &modulus - &a;
                            let mut a_bytes = a.to_bytes_le();
                            a_bytes.resize(48, 0);

                            let a_fp = $curve::from_bytes_unsafe(&a_bytes.try_into().unwrap());
                            (a, a_inv, a_fp)
                        };

                        for _ in 0..100 {
                            let (a, a_inv, a_fp) = gen_test_value();
                            let (b, b_inv, b_fp) = gen_test_value();

                            let lhs = Fp2::new(a_fp, b_fp).lexicographically_largest();
                            let rhs = b > b_inv || (b == BigUint::ZERO && a > a_inv);

                            assert_eq!(lhs, rhs);
                        }

                        for _ in 0..100 {
                            let (a, a_inv, a_fp) = gen_test_value();
                            let b = BigUint::ZERO;
                            let b_inv = BigUint::ZERO;
                            let b_fp = $curve::zero();

                            let lhs = Fp2::new(a_fp, b_fp).lexicographically_largest();
                            let rhs = b > b_inv || (b == BigUint::ZERO && a > a_inv);

                            assert_eq!(lhs, rhs);
                        }
                    }
                }
            }
        };
    }

    fp2_tests!(Bls12381, bls12381_fp2_rand, bls12381_fp2_test);
    fp2_tests!(Bn254, bn254_fp2_rand, bn254_fp2_test);

    #[test]
    fn test_frobenius() {
        use rand::thread_rng;
        for _ in 0..10 {
            let f = Fp2::<Bls12381>::random(&mut thread_rng());
            let frob = f.frobenius_map().frobenius_map();
            assert_eq!(f, frob);
        }
    }
}
