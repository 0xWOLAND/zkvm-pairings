use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::RngCore;
use std::marker::PhantomData;

use crate::common::Curve;
use crate::fp::Fp;

#[derive(Copy, Clone)]
/// Represents an element in the field Fp2.
pub struct Fp2<C: Curve> {
    /// The first component of the Fp2 element.
    pub c0: Fp<C>,
    /// The second component of the Fp2 element.
    pub c1: Fp<C>,
    _marker: PhantomData<C>,
}

impl<C: Curve> fmt::Debug for Fp2<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + {:?}*u", self.c0, self.c1)
    }
}

impl<C: Curve> Default for Fp2<C> {
    fn default() -> Self {
        Fp2::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<C: Curve> zeroize::DefaultIsZeroes for Fp2<C> {}

impl<C: Curve> From<Fp<C>> for Fp2<C> {
    fn from(f: Fp<C>) -> Fp2<C> {
        Fp2 {
            c0: f,
            c1: Fp::zero(),
            _marker: PhantomData::<C>,
        }
    }
}

impl<C: Curve> Eq for Fp2<C> {}
impl<C: Curve> PartialEq for Fp2<C> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.c0.eq(&other.c0) & self.c1.eq(&other.c1)
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

impl<'a, C: Curve> Neg for &'a Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn neg(self) -> Fp2<C> {
        self.neg()
    }
}

impl<C: Curve> Neg for Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn neg(self) -> Fp2<C> {
        -&self
    }
}

impl<'a, 'b, C: Curve> Sub<&'b Fp2<C>> for &'a Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn sub(self, rhs: &'b Fp2<C>) -> Fp2<C> {
        self.sub(rhs)
    }
}

impl<'a, 'b, C: Curve> Add<&'b Fp2<C>> for &'a Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn add(self, rhs: &'b Fp2<C>) -> Fp2<C> {
        self.add(rhs)
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fp2<C>> for &'a Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn mul(self, rhs: &'b Fp2<C>) -> Fp2<C> {
        self.mul(rhs)
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fp<C>> for &'a Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn mul(self, rhs: &'b Fp<C>) -> Fp2<C> {
        Fp2::new(self.c0 * rhs, self.c1 * rhs)
    }
}

impl<'a, 'b, C: Curve> Div<&'b Fp2<C>> for &'a Fp2<C> {
    type Output = Fp2<C>;

    #[inline]
    fn div(self, rhs: &'b Fp2<C>) -> Fp2<C> {
        self.div(rhs)
    }
}

impl_binops_additive!(Fp2<C>, Fp2<C>);
impl_binops_multiplicative!(Fp2<C>, Fp2<C>);
impl_binops_multiplicative!(Fp2<C>, Fp<C>);
impl_binops_divisible!(Fp2<C>, Fp2<C>);

impl<C: Curve> Fp2<C> {
    /// Returns the zero element of Fp2.
    #[inline]
    pub const fn zero() -> Fp2<C> {
        Fp2::new(Fp::zero(), Fp::zero())
    }

    /// Returns the one element of Fp2.
    #[inline]
    pub const fn one() -> Fp2<C> {
        Fp2::new(Fp::one(), Fp::zero())
    }

    pub const fn new(c0: Fp<C>, c1: Fp<C>) -> Fp2<C> {
        Fp2 {
            c0,
            c1,
            _marker: PhantomData::<C>,
        }
    }

    /// Checks if this element is zero.
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Generates a random element in Fp2.
    pub fn random(mut rng: impl RngCore) -> Fp2<C> {
        Fp2::new(Fp::random(&mut rng), Fp::random(&mut rng))
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
    pub fn mul_by_nonresidue(&self) -> Fp2<C> {
        // Multiply a + bu by u + 1, getting
        // au + a + bu^2 + bu
        // and because u^2 = -1, we get
        // (a - b) + (a + b)u

        Fp2::new(self.c0 - self.c1, self.c0 + self.c1)
    }

    /// Computes the square of this element.
    pub fn square(&self) -> Fp2<C> {
        // Complex squaring:
        //
        // v0  = c0 * c1
        // c0' = (c0 + c1) * (c0 + \beta*c1) - v0 - \beta * v0
        // c1' = 2 * v0
        //
        // In BLS12-381's F_{p^2}, our \beta is -1 so we
        // can modify this formula:
        //
        // c0' = (c0 + c1) * (c0 - c1)
        // c1' = 2 * c0 * c1

        let a = (&self.c0).add(&self.c1);
        let b = (&self.c0).sub(&self.c1);
        let c = (&self.c0).add(&self.c0);

        Fp2::new((&a).mul(&b), (&c).mul(&self.c1))
    }

    /// Multiplies this element by another element.
    pub fn mul(&self, rhs: &Fp2<C>) -> Fp2<C> {
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
            self.c0 * rhs.c0 - self.c1 * rhs.c1,
            self.c0 * rhs.c1 + self.c1 * rhs.c0,
        )
    }

    pub fn div(&self, rhs: &Fp2<C>) -> Fp2<C> {
        self * rhs.invert().unwrap()
    }

    /// Adds another element to this element.
    pub fn add(&self, rhs: &Fp2<C>) -> Fp2<C> {
        Fp2::new((&self.c0).add(&rhs.c0), (&self.c1).add(&rhs.c1))
    }

    /// Subtracts another element from this element.
    pub fn sub(&self, rhs: &Fp2<C>) -> Fp2<C> {
        Fp2::new((&self.c0).sub(&rhs.c0), (&self.c1).sub(&rhs.c1))
    }

    /// Negates this element.
    pub const fn neg(&self) -> Fp2<C> {
        Fp2::new((&self.c0).neg(), (&self.c1).neg())
    }

    /// Computes the square root of this element.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Fp2::<C>::zero());
        }

        // a1 = self^((p - 3) / 4)
        let a1 = self.pow_vartime(&[
            0xee7f_bfff_ffff_eaaa,
            0x07aa_ffff_ac54_ffff,
            0xd9cc_34a8_3dac_3d89,
            0xd91d_d2e1_3ce1_44af,
            0x92c6_e9ed_90d2_eb35,
            0x0680_447a_8e5f_f9a6,
        ]);

        // alpha = a1^2 * self = self^((p - 3) / 2 + 1) = self^((p - 1) / 2)
        let alpha = a1.square() * self;
        // x0 = self^((p + 1) / 4)
        let x0 = a1 * self;

        if alpha == -Fp2::one() {
            // The element is order p - 1, so we're just trying to get the square of an element of the subfield Fp.
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

            // Only return the result if it's really the square root (and so self is actually quadratic nonresidue)
            if sqrt.square() == *self {
                Some(sqrt)
            } else {
                None
            }
        }
    }

    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    pub fn invert(&self) -> Option<Self> {
        // We wish to find the multiplicative inverse of a nonzero
        // element a + bu in Fp2. We leverage an identity
        //
        // (a + bu)(a - bu) = a^2 + b^2
        //
        // which holds because u^2 = -1. This can be rewritten as
        //
        // (a + bu)(a - bu)/(a^2 + b^2) = 1
        //
        // because a^2 + b^2 = 0 has no nonzero solutions for (a, b).
        // This gives that (a - bu)/(a^2 + b^2) is the inverse
        // of (a + bu). Importantly, this can be computing using
        // only a single inversion in Fp.

        (self.c0.square() + self.c1.square())
            .invert()
            .map(|t| Fp2::new(self.c0 * t, self.c1 * -t))
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

    pub(crate) fn pow_vartime_extended(&self, by: &[u64]) -> Self {
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
}

#[cfg(test)]
mod test {
    use rand::Rng;

    use crate::common::Bls12381Curve;

    use super::*;

    fn fp2_rand() -> Fp2<Bls12381Curve> {
        let mut rng = rand::thread_rng();
        Fp2::random(&mut rng)
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp2::<Bls12381Curve>::new(
                Fp::from_raw_unchecked(x.clone().try_into().unwrap()),
                Fp::from_raw_unchecked(y.clone().try_into().unwrap()),
            );

            let b = Fp2::<Bls12381Curve>::new(
                Fp::from_raw_unchecked(x.clone().try_into().unwrap()),
                Fp::from_raw_unchecked(y.clone().try_into().unwrap()),
            );

            assert_eq!(a, b)
        }
    }

    #[test]
    fn test_inequality() {
        for _ in 0..10 {
            let a = fp2_rand();
            let b = fp2_rand();
            assert_ne!(a, b)
        }
    }

    #[test]
    fn test_addition_subtraction() {
        for _ in 0..10 {
            let a = fp2_rand();
            let b = fp2_rand();
            let c = fp2_rand();

            // commutative
            assert_eq!(a + b, b + a);
            assert_eq!(a + (b + c), (a + b) + c);

            // additive identity
            assert_eq!(a + Fp2::zero(), a); // a + 0 = a
            assert_eq!(a - Fp2::zero(), a); // subtraction identity

            assert_eq!(Fp2::zero() - a, -a); // 0 - a = -a
            assert_eq!(a - b, a + (-b)); // a - b = a + -b
            assert_eq!(a - b, a + (b * -Fp2::one())); // a - b = a + b * -1

            assert_eq!(-a, Fp2::zero() - a);
            assert_eq!(-a, a * -Fp2::one());
        }
    }

    #[test]
    fn test_multiplication() {
        for _ in 0..10 {
            let a = fp2_rand();
            let b = fp2_rand();
            let c = fp2_rand();

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
            let a = fp2_rand();

            assert_eq!(a * Fp::from(0), Fp2::zero());
            assert_eq!(a * Fp2::zero(), Fp2::zero());
            assert_eq!(a * Fp2::one(), a);
            assert_eq!(a * Fp::from(1), a);
            assert_eq!(a * Fp::from(2), a + a);
            assert_eq!(a * Fp::from(3), a + a + a);
            assert_eq!(a * Fp::from(4), a + a + a + a);
        }
    }

    #[test]
    fn test_square_equality() {
        for _ in 0..10 {
            let a = fp2_rand();
            assert_eq!(a.square(), a * a);
        }
    }

    #[test]
    fn test_pow_equality() {
        for _ in 0..10 {
            let a = fp2_rand();
            assert_eq!(a.pow_vartime(&[1, 0, 0, 0, 0, 0]), a);
            assert_eq!(a.pow_vartime(&[2, 0, 0, 0, 0, 0]), a.square());
            assert_eq!(a.pow_vartime(&[3, 0, 0, 0, 0, 0]), a.square() * a);
            assert_eq!(a.pow_vartime(&[4, 0, 0, 0, 0, 0]), a.square().square());
        }
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..100 {
            let a = fp2_rand();
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
            let a = fp2_rand();

            // division by one
            assert_eq!(a / Fp2::one(), a);
            assert_eq!(a / a, Fp2::one());

            // division by zero
            assert_eq!(Fp2::zero() / a, Fp2::zero());

            // division distributivity
            let a = fp2_rand();
            let b = fp2_rand();
            let c = fp2_rand();

            assert_eq!((a + b) / c, a / c + b / c);

            // division and multiplication equality
            let a = fp2_rand();
            let b = fp2_rand();
            assert_eq!(a / b, a * b.invert().unwrap());
        }
    }

    #[test]
    fn test_inversion() {
        for _ in 0..10 {
            let a = fp2_rand();
            assert_eq!(a * a.invert().unwrap(), Fp2::one());
            assert_eq!(a.invert().unwrap().invert().unwrap(), a);
        }
    }
}

#[cfg(feature = "zeroize")]
#[test]
fn test_zeroize() {
    use zeroize::Zeroize;

    let mut a = Fp2::one();
    a.zeroize();
    assert!(bool::from(a.is_zero()));
}
