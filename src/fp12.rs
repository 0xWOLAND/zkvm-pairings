use crate::fp::*;
use crate::fp2::*;
use crate::fp6::*;

use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

use num_bigint::BigUint;
use rand::RngCore;
#[cfg(feature = "pairings")]
use rand_core::RngCore;

pub trait Fp12Element: Fp6Element {
    type Fp12ElementType;
    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp12ElementType;
    fn to_bytes_vec(value: &Self::Fp12ElementType) -> Vec<u8>;
    fn get_fp12_frobenius_coeffs(pow: usize) -> [Self; 2];
}

impl Fp12Element for Bls12381 {
    type Fp12ElementType = Fp12<Bls12381>;
    fn get_fp12_frobenius_coeffs(pow: usize) -> [Self; 2] {
        match pow % 12 {
            0 => [Self::one(); 2],
            1 => [
                Self::from_raw_unchecked([
                    0x8d0775ed92235fb8,
                    0xf67ea53d63e7813d,
                    0x7b2443d784bab9c4,
                    0x0fd603fd3cbd5f4f,
                    0xc231beb4202c0d1f,
                    0x1904d3bf02bb0667,
                ]),
                Self::from_raw_unchecked([
                    0x2cf78a126ddc4af3,
                    0x282d5ac14d6c7ec2,
                    0xec0c8ec971f63c5f,
                    0x54a14787b6c7b36f,
                    0x88e9e902231f9fb8,
                    0x00fc3e2b36c4e032,
                ]),
            ],
            2 => [
                Self::from_raw_unchecked([
                    0x2e01fffffffeffff,
                    0xde17d813620a0002,
                    0xddb3a93be6f89688,
                    0xba69c6076a0f77ea,
                    0x5f19672fdf76ce51,
                    0x0000000000000000,
                ]),
                Self::zero(),
            ],
            3 => [
                Self::from_raw_unchecked([
                    0xf1ee7b04121bdea2,
                    0x304466cf3e67fa0a,
                    0xef396489f61eb45e,
                    0x1c3dedd930b1cf60,
                    0xe2e9c448d77a2cd9,
                    0x135203e60180a68e,
                ]),
                Self::from_raw_unchecked([
                    0xc81084fbede3cc09,
                    0xee67992f72ec05f4,
                    0x77f76e17009241c5,
                    0x48395dabc2d3435e,
                    0x6831e36d6bd17ffe,
                    0x06af0e0437ff400b,
                ]),
            ],
            4 => [
                Self::from_raw_unchecked([
                    0x2e01fffffffefffe,
                    0xde17d813620a0002,
                    0xddb3a93be6f89688,
                    0xba69c6076a0f77ea,
                    0x5f19672fdf76ce51,
                    0x0000000000000000,
                ]),
                Self::zero(),
            ],
            5 => [
                Self::from_raw_unchecked([
                    0x1ee605167ff82995,
                    0x5871c1908bd478cd,
                    0xdb45f3536814f0bd,
                    0x70df3560e77982d0,
                    0x6bd3ad4afa99cc91,
                    0x144e4211384586c1,
                ]),
                Self::from_raw_unchecked([
                    0x9b18fae980078116,
                    0xc63a3e6e257f8732,
                    0x8beadf4d8e9c0566,
                    0xf39816240c0b8fee,
                    0xdf47fa6b48b1e045,
                    0x05b2cfd9013a5fd8,
                ]),
            ],
            6 => [
                Self::from_raw_unchecked([
                    0xb9feffffffffaaaa,
                    0x1eabfffeb153ffff,
                    0x6730d2a0f6b0f624,
                    0x64774b84f38512bf,
                    0x4b1ba7b6434bacd7,
                    0x1a0111ea397fe69a,
                ]),
                Self::zero(),
            ],
            7 => [
                Self::from_raw_unchecked([
                    0x2cf78a126ddc4af3,
                    0x282d5ac14d6c7ec2,
                    0xec0c8ec971f63c5f,
                    0x54a14787b6c7b36f,
                    0x88e9e902231f9fb8,
                    0x00fc3e2b36c4e032,
                ]),
                Self::from_raw_unchecked([
                    0x8d0775ed92235fb8,
                    0xf67ea53d63e7813d,
                    0x7b2443d784bab9c4,
                    0x0fd603fd3cbd5f4f,
                    0xc231beb4202c0d1f,
                    0x1904d3bf02bb0667,
                ]),
            ],
            8 => [
                Self::from_raw_unchecked([
                    0x8bfd00000000aaac,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
                Self::zero(),
            ],
            9 => [
                Self::from_raw_unchecked([
                    0xc81084fbede3cc09,
                    0xee67992f72ec05f4,
                    0x77f76e17009241c5,
                    0x48395dabc2d3435e,
                    0x6831e36d6bd17ffe,
                    0x06af0e0437ff400b,
                ]),
                Self::from_raw_unchecked([
                    0xf1ee7b04121bdea2,
                    0x304466cf3e67fa0a,
                    0xef396489f61eb45e,
                    0x1c3dedd930b1cf60,
                    0xe2e9c448d77a2cd9,
                    0x135203e60180a68e,
                ]),
            ],
            10 => [
                Self::from_raw_unchecked([
                    0x8bfd00000000aaad,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
                Self::zero(),
            ],
            11 => [
                Self::from_raw_unchecked([
                    0x9b18fae980078116,
                    0xc63a3e6e257f8732,
                    0x8beadf4d8e9c0566,
                    0xf39816240c0b8fee,
                    0xdf47fa6b48b1e045,
                    0x05b2cfd9013a5fd8,
                ]),
                Self::from_raw_unchecked([
                    0x1ee605167ff82995,
                    0x5871c1908bd478cd,
                    0xdb45f3536814f0bd,
                    0x70df3560e77982d0,
                    0x6bd3ad4afa99cc91,
                    0x144e4211384586c1,
                ]),
            ],
            _ => unimplemented!(),
        }
    }

    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp12ElementType {
        let c0 = <Bls12381 as Fp6Element>::from_bytes_slice(&bytes[..288]);
        let c1 = <Bls12381 as Fp6Element>::from_bytes_slice(&bytes[288..]);
        Fp12::<Bls12381>::new(c0, c1)
    }

    fn to_bytes_vec(value: &Self::Fp12ElementType) -> Vec<u8> {
        let mut res = [0u8; 576];
        let c0 = <Bls12381 as Fp6Element>::to_bytes_vec(&value.c0);
        let c1 = <Bls12381 as Fp6Element>::to_bytes_vec(&value.c1);

        res[..288].copy_from_slice(&c0);
        res[288..].copy_from_slice(&c1);

        res.to_vec()
    }
}

impl Fp12Element for Bn254 {
    type Fp12ElementType = Fp12<Bn254>;
    fn get_fp12_frobenius_coeffs(pow: usize) -> [Self; 2] {
        match pow % 12 {
            0 => [Self::one(); 2],
            1 => [
                Self::from_raw_unchecked([
                    12653890742059813127,
                    14585784200204367754,
                    1278438861261381767,
                    212598772761311868,
                ]),
                Self::from_raw_unchecked([
                    11683091849979440498,
                    14992204589386555739,
                    15866167890766973222,
                    1200023580730561873,
                ]),
            ],
            2 => [
                Self::from_raw_unchecked([
                    14595462726357228530,
                    17349508522658994025,
                    1017833795229664280,
                    299787779797702374,
                ]),
                Self::zero(),
            ],
            3 => [
                Self::from_raw_unchecked([
                    3914496794763385213,
                    790120733010914719,
                    7322192392869644725,
                    581366264293887267,
                ]),
                Self::from_raw_unchecked([
                    12817045492518885689,
                    4440270538777280383,
                    11178533038884588256,
                    2767537931541304486,
                ]),
            ],
            _ => unimplemented!(),
        }
    }

    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp12ElementType {
        let c0 = <Bn254 as Fp6Element>::from_bytes_slice(&bytes[..288]);
        let c1 = <Bn254 as Fp6Element>::from_bytes_slice(&bytes[288..]);
        Fp12::<Bn254>::new(c0, c1)
    }

    fn to_bytes_vec(value: &Self::Fp12ElementType) -> Vec<u8> {
        let mut res = [0u8; 576];
        let c0 = <Bn254 as Fp6Element>::to_bytes_vec(&value.c0);
        let c1 = <Bn254 as Fp6Element>::to_bytes_vec(&value.c1);

        res[..288].copy_from_slice(&c0);
        res[288..].copy_from_slice(&c1);

        res.to_vec()
    }
}

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
pub struct Fp12<F: Fp12Element> {
    pub c0: Fp6<F>,
    pub c1: Fp6<F>,
}

impl<F: Fp12Element> From<F> for Fp12<F> {
    fn from(f: F) -> Fp12<F> {
        Fp12 {
            c0: Fp6::<F>::from(f),
            c1: Fp6::<F>::from(f),
        }
    }
}

impl<F: Fp12Element> From<Fp2<F>> for Fp12<F> {
    fn from(f: Fp2<F>) -> Fp12<F> {
        Fp12 {
            c0: Fp6::<F>::from(f),
            c1: Fp6::<F>::zero(),
        }
    }
}

impl<F: Fp12Element> From<Fp6<F>> for Fp12<F> {
    fn from(f: Fp6<F>) -> Fp12<F> {
        Fp12 {
            c0: f,
            c1: Fp6::<F>::zero(),
        }
    }
}

impl<F: Fp12Element> Eq for Fp12<F> {}
impl<F: Fp12Element> PartialEq for Fp12<F> {
    fn eq(&self, other: &Fp12<F>) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<F: Fp12Element> Copy for Fp12<F> {}
impl<F: Fp12Element> Clone for Fp12<F> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<F: Fp12Element> Default for Fp12<F> {
    fn default() -> Self {
        Fp12::<F>::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<F: Fp12Element> zeroize::DefaultIsZeroes for Fp12 {}

impl<F: Fp12Element> fmt::Debug for Fp12<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*w", self.c0, self.c1)
    }
}

impl<F: Fp12Element> Fp12<F> {
    #[inline]
    pub fn new(c0: Fp6<F>, c1: Fp6<F>) -> Self {
        Fp12 { c0, c1 }
    }

    #[inline]
    pub fn zero() -> Self {
        Fp12::new(Fp6::zero(), Fp6::zero())
    }

    #[inline]
    pub fn one() -> Self {
        Fp12::new(Fp6::one(), Fp6::zero())
    }

    pub fn random(mut rng: impl RngCore) -> Self {
        Fp12 {
            c0: Fp6::<F>::random(&mut rng),
            c1: Fp6::<F>::random(&mut rng),
        }
    }

    pub fn mul_by_014(&self, c0: &Fp2<F>, c1: &Fp2<F>, c4: &Fp2<F>) -> Fp12<F> {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = *c1 + *c4;
        let c1 = self.c1 + self.c0;
        let c1 = c1.mul_by_01(c0, &o);
        let c1 = c1 - aa - bb;
        let c0 = bb;
        let c0 = c0.mul_by_nonresidue();
        let c0 = c0 + aa;

        Fp12 { c0, c1 }
    }

    pub fn mul_14_by_14(d0: &Fp2<F>, d1: &Fp2<F>, c0: &Fp2<F>, c1: &Fp2<F>) -> [Fp2<F>; 5] {
        let x0 = *d0 * *c0;
        let x1 = *d1 * *c1;
        let x04 = *c0 + *d0;
        let tmp = *c0 + *c1;
        let x01 = *d0 + *d1;
        let x01 = x01 * tmp;
        let tmp = x1 + x0;
        let x01 = x01 - tmp;
        let x14 = *c1 + *d1;
        let z_c0_b0 = Fp2::<F>::non_residue() + x0;

        [z_c0_b0, x01, x1, x04, x14]
    }

    fn cyclotomic_square(&self) -> Fp12<F> {
        let t0 = self.c1.c1.square();
        let t1 = self.c0.c0.square();
        let t6 = (self.c1.c1 + self.c0.c0).square();
        let t6 = t6 - t0;
        let t6 = t6 - t1;
        let t2 = &self.c0.c2.square();
        let t3 = &self.c1.c0.square();
        let t7 = (self.c0.c2 + self.c1.c0).square();
        let t7 = t7 - *t2;
        let t7 = t7 - *t3;
        let t4 = self.c1.c2.square();
        let t5 = self.c0.c1.square();
        let t8 = (self.c1.c2 + self.c0.c1).square();
        let t8 = t8 - t4;
        let t8 = t8 - t5;
        let t8 = t8.mul_by_nonresidue();
        let t0 = t0.mul_by_nonresidue();
        let t0 = t0 + t1;
        let t2 = t2.mul_by_nonresidue();
        let t2 = t2 + *t3;
        let t4 = t4.mul_by_nonresidue();
        let t4 = t4 + t5;
        let z00 = t0 - self.c0.c0;
        let z00 = z00 + z00;
        let z00 = z00 + t0;
        let z01 = t2 - self.c0.c1;
        let z01 = z01 + z01;
        let z01 = z01 + t2;
        let z02 = t4 - self.c0.c2;
        let z02 = z02 + z02;
        let z02 = z02 + t4;
        let z10 = t8 + self.c1.c0;
        let z10 = z10 + z10;
        let z10 = z10 + t8;
        let z11 = t6 + self.c1.c1;
        let z11 = z11 + z11;
        let z11 = z11 + t6;
        let z12 = t7 + self.c1.c2;
        let z12 = z12 + z12;
        let z12 = z12 + t7;
        Fp12::new(Fp6::new(z00, z01, z02), Fp6::new(z10, z11, z12))
    }

    fn n_cyclotomic_square(&self, by: u64) -> Fp12<F> {
        (0..by).fold(*self, |acc, _| acc.cyclotomic_square())
    }

    pub fn powt(&self) -> Fp12<F> {
        let a = self.cyclotomic_square();
        let a = a * *self;
        let a = a.n_cyclotomic_square(2);
        let a = a * *self;
        let a = a.n_cyclotomic_square(3);
        let a = a * *self;
        let a = a.n_cyclotomic_square(9);
        let a = a * *self;
        let a = a.n_cyclotomic_square(32);
        let a = a * *self;
        let a = a.n_cyclotomic_square(15);
        let a = a * *self;
        a.cyclotomic_square()
    }

    pub fn div(&self, rhs: &Fp12<F>) -> Fp12<F> {
        rhs.invert().unwrap() * *self
    }

    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    #[inline(always)]
    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    #[inline(always)]
    pub fn conjugate(&self) -> Self {
        Fp12::new(self.c0, -self.c1)
    }

    pub fn pow_vartime(&self, by: &[u64; 6]) -> Self {
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

    pub fn pow_vartime_extended(&self, by: &[u64]) -> Self {
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

    pub fn pow_vartime_extended_str(&self, by: &str) -> Self {
        self.pow_vartime_extended(&BigUint::from_str(by).unwrap().to_u64_digits())
    }
    /// Raises this element to p.
    #[inline(always)]
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();

        let frob_coeffs = F::get_fp12_frobenius_coeffs(1);
        let c1 = c1 * Fp6::from(Fp2::new(frob_coeffs[0], frob_coeffs[1]));

        Fp12::new(c0, c1)
    }

    #[inline(always)]
    pub(crate) fn nth_frobenius_map(&self, pow: usize) -> Self {
        let c0 = self.c0.nth_frobenius_map(pow);
        let c1 = self.c1.nth_frobenius_map(pow);

        let frob_coeffs = F::get_fp12_frobenius_coeffs(pow);
        let c1 = c1 * Fp6::from(Fp2::new(frob_coeffs[0], frob_coeffs[1]));

        Fp12::new(c0, c1)
    }

    #[inline]
    pub fn square(&self) -> Self {
        let ab = self.c0 * self.c1;
        let c0c1 = self.c0 + self.c1;
        let c0 = self.c1.mul_by_nonresidue();
        let c0 = c0 + self.c0;
        let c0 = c0 * c0c1;
        let c0 = c0 - ab;
        let c1 = ab + ab;
        let c0 = c0 - ab.mul_by_nonresidue();

        Fp12::new(c0, c1)
    }

    pub fn invert(&self) -> Option<Self> {
        (self.c0.square() - self.c1.square().mul_by_nonresidue())
            .invert()
            .map(|t| Fp12::new(self.c0 * t, self.c1 * -t))
    }
}

impl<F: Fp12Element> Mul<Fp12<F>> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn mul(self, other: Fp12<F>) -> Self::Output {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;
        let o = other.c0 + other.c1;
        let c1 = self.c1 + self.c0;
        let c1 = c1 * o;
        let c1 = c1 - aa;
        let c1 = c1 - bb;
        let c0 = bb.mul_by_nonresidue();
        let c0 = c0 + aa;

        Fp12::new(c0, c1)
    }
}

impl<F: Fp12Element> Add<Fp12<F>> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn add(self, rhs: Fp12<F>) -> Self::Output {
        Fp12::new(self.c0 + rhs.c0, self.c1 + rhs.c1)
    }
}

impl<'a, F: Fp12Element> Neg for &'a Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp12::new(-self.c0, -self.c1)
    }
}

impl<F: Fp12Element> Neg for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F: Fp12Element> Sub<Fp12<F>> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn sub(self, rhs: Fp12<F>) -> Self::Output {
        Fp12::new(self.c0 - rhs.c0, self.c1 - rhs.c1)
    }
}

impl<F: Fp12Element> Mul<F> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn mul(self, rhs: F) -> Fp12<F> {
        let rhs = Fp12::from(rhs);
        Fp12::new(self.c0 * rhs.c0, self.c1 * rhs.c1)
    }
}

impl<'a, 'b, F: Fp12Element> Div<&'b Fp12<F>> for &'a Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn div(self, rhs: &'b Fp12<F>) -> Fp12<F> {
        self.div(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    macro_rules! fp12_tests {
        ($curve:ident, $rand_fn:ident, $curve_test:ident) => {
            mod $curve_test {
                use super::*;

                #[test]
                fn test_equality() {
                    let rng = &mut rand::thread_rng();
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
                        assert_eq!(a + Fp12::<$curve>::zero(), a);
                        assert_eq!(a - Fp12::<$curve>::zero(), a);

                        assert_eq!(Fp12::<$curve>::zero() - a, -a);
                        assert_eq!(a - b, a + (-b));
                        assert_eq!(a - b, a + (b * -Fp12::<$curve>::one()));

                        assert_eq!(-a, Fp12::<$curve>::zero() - a);
                        assert_eq!(-a, a * -Fp12::<$curve>::one());
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
                fn test_square_equality() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        assert_eq!(a.square(), a * a);
                    }
                }

                #[test]
                fn test_inversion() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        if !a.is_zero() {
                            assert_eq!(a * a.invert().unwrap(), Fp12::<$curve>::one());
                            assert_eq!(a.invert().unwrap().invert().unwrap(), a);
                        }
                    }
                }

                #[test]
                fn test_frobenius() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = (0..12).fold(a, |acc, _| acc.frobenius_map());
                        assert_eq!(a, b);
                    }
                }

                #[test]
                fn test_cyclotomic_square() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        assert_eq!(a.cyclotomic_square(), a.n_cyclotomic_square(1));
                        assert_eq!(
                            a.cyclotomic_square().cyclotomic_square(),
                            a.n_cyclotomic_square(2)
                        );
                    }
                }
            }
        };
    }

    fn bls12381_fp12_rand() -> Fp12<Bls12381> {
        let mut rng = rand::thread_rng();
        Fp12::new(
            Fp6::new(
                Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
            ),
            Fp6::new(
                Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
            ),
        )
    }

    fn bn254_fp12_rand() -> Fp12<Bn254> {
        let mut rng = rand::thread_rng();
        Fp12::new(
            Fp6::new(
                Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
            ),
            Fp6::new(
                Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
            ),
        )
    }

    fp12_tests!(Bls12381, bls12381_fp12_rand, bls12381_fp12_test);
    fp12_tests!(Bn254, bn254_fp12_rand, bn254_fp12_test);
}
