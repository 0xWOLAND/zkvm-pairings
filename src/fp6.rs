use crate::fp::*;
use crate::fp2::*;

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use rand::RngCore;
#[cfg(feature = "pairings")]
use rand_core::RngCore;

pub(crate) trait Fp6Element: Fp2Element {
    type Fp6ElementType;
    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp6ElementType;
    fn to_bytes_vec(f: &Self::Fp6ElementType) -> Vec<u8>;
    fn get_fp6_frobenius_coeff(i: usize) -> (Self, Self);
}

impl Fp6Element for Bls12381 {
    type Fp6ElementType = Fp6<Bls12381>;
    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp6ElementType {
        let mut c0_bytes = [0u8; 96];
        let mut c1_bytes = [0u8; 96];
        let mut c2_bytes = [0u8; 96];

        c0_bytes.copy_from_slice(&bytes[0..96]);
        c1_bytes.copy_from_slice(&bytes[96..192]);
        c2_bytes.copy_from_slice(&bytes[192..288]);

        let c0 = <Bls12381 as Fp2Element>::from_bytes_slice(&c0_bytes);
        let c1 = <Bls12381 as Fp2Element>::from_bytes_slice(&c1_bytes);
        let c2 = <Bls12381 as Fp2Element>::from_bytes_slice(&c2_bytes);

        Fp6::<Bls12381>::new(c0, c1, c2)
    }

    fn to_bytes_vec(f: &Self::Fp6ElementType) -> Vec<u8> {
        let mut res = [0u8; 288];
        let c0_bytes = <Bls12381 as Fp2Element>::to_bytes_vec(&f.c0);
        let c1_bytes = <Bls12381 as Fp2Element>::to_bytes_vec(&f.c1);
        let c2_bytes = <Bls12381 as Fp2Element>::to_bytes_vec(&f.c2);

        res[0..96].copy_from_slice(&c0_bytes);
        res[96..192].copy_from_slice(&c1_bytes);
        res[192..288].copy_from_slice(&c2_bytes);

        res.to_vec()
    }

    fn get_fp6_frobenius_coeff(pow: usize) -> (Self, Self) {
        match pow % 6 {
            0 => (Self::one(), Self::one()),
            1 => (
                Self::zero(),
                Self::from_raw_unchecked([
                    0x8bfd00000000aaac,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
            ),
            2 => (
                Self::from_raw_unchecked([
                    0x8bfd00000000aaad,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
                Self::zero(),
            ),
            _ => panic!("Invalid power"),
        }
    }
}

impl Fp6Element for Bn254 {
    type Fp6ElementType = Fp6<Bn254>;
    fn get_fp6_frobenius_coeff(pow: usize) -> (Self, Self) {
        match pow % 6 {
            0 => (Self::one(), Self::one()),
            1 => (
                Self::from_raw_unchecked([
                    13075984984163199792,
                    3782902503040509012,
                    8791150885551868305,
                    1825854335138010348,
                ]),
                Self::from_raw_unchecked([
                    7963664994991228759,
                    12257807996192067905,
                    13179524609921305146,
                    2767831111890561987,
                ]),
            ),
            2 => (
                Self::from_raw_unchecked([
                    3697675806616062876,
                    9065277094688085689,
                    6918009208039626314,
                    2775033306905974752,
                ]),
                Self::zero(),
            ),
            _ => panic!("Invalid power"),
        }
    }

    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp6ElementType {
        let mut c0_bytes = [0u8; 96];
        let mut c1_bytes = [0u8; 96];
        let mut c2_bytes = [0u8; 96];

        c0_bytes.copy_from_slice(&bytes[0..96]);
        c1_bytes.copy_from_slice(&bytes[96..192]);
        c2_bytes.copy_from_slice(&bytes[192..288]);

        let c0 = <Bn254 as Fp2Element>::from_bytes_slice(&c0_bytes);
        let c1 = <Bn254 as Fp2Element>::from_bytes_slice(&c1_bytes);
        let c2 = <Bn254 as Fp2Element>::from_bytes_slice(&c2_bytes);

        Fp6::<Bn254>::new(c0, c1, c2)
    }

    fn to_bytes_vec(f: &Self::Fp6ElementType) -> Vec<u8> {
        let mut res = [0u8; 288];
        let c0_bytes = <Bn254 as Fp2Element>::to_bytes_vec(&f.c0);
        let c1_bytes = <Bn254 as Fp2Element>::to_bytes_vec(&f.c1);
        let c2_bytes = <Bn254 as Fp2Element>::to_bytes_vec(&f.c2);

        res[0..96].copy_from_slice(&c0_bytes);
        res[96..192].copy_from_slice(&c1_bytes);
        res[192..288].copy_from_slice(&c2_bytes);

        res.to_vec()
    }
}
/// This represents an element $c_0 + c_1 v + c_2 v^2$ of $\mathbb{F}_{p^6} = \mathbb{F}_{p^2} / v^3 - u - 1$.
pub struct Fp6<F: Fp6Element> {
    pub c0: Fp2<F>,
    pub c1: Fp2<F>,
    pub c2: Fp2<F>,
}

impl<F: Fp6Element> From<F> for Fp6<F> {
    fn from(f: F) -> Fp6<F> {
        Fp6 {
            c0: Fp2::from(f),
            c1: Fp2::from(f),
            c2: Fp2::from(f),
        }
    }
}

impl<F: Fp6Element> From<u64> for Fp6<F> {
    fn from(f: u64) -> Fp6<F> {
        Fp6 {
            c0: Fp2::from(f),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }
}

impl<F: Fp6Element> From<Fp2<F>> for Fp6<F> {
    fn from(f: Fp2<F>) -> Fp6<F> {
        Fp6 {
            c0: f,
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }
}

impl<F: Fp6Element> PartialEq for Fp6<F> {
    fn eq(&self, other: &Fp6<F>) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}

impl<F: Fp6Element> Copy for Fp6<F> {}
impl<F: Fp6Element> Clone for Fp6<F> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<F: Fp6Element> Default for Fp6<F> {
    fn default() -> Self {
        Fp6::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<F: Fp6Element> zeroize::DefaultIsZeroes for Fp6<F> {}

impl<F: Fp6Element> fmt::Debug for Fp6<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*v + ({:?})*v^2", self.c0, self.c1, self.c2)
    }
}

impl<F: Fp6Element> Eq for Fp6<F> {}

impl<F: Fp6Element> Fp6<F> {
    #[inline]
    pub fn new(c0: Fp2<F>, c1: Fp2<F>, c2: Fp2<F>) -> Self {
        Fp6 { c0, c1, c2 }
    }

    #[inline]
    pub fn zero() -> Self {
        Fp6 {
            c0: Fp2::zero(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    #[inline]
    pub fn one() -> Self {
        Fp6 {
            c0: Fp2::one(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    pub(crate) fn random(mut rng: impl RngCore) -> Fp6<F> {
        Fp6 {
            c0: Fp2::random(&mut rng),
            c1: Fp2::random(&mut rng),
            c2: Fp2::random(&mut rng),
        }
    }

    pub fn mul_by_1(&self, c1: &Fp2<F>) -> Fp6<F> {
        Fp6 {
            c0: (self.c2 * *c1).mul_by_nonresidue(),
            c1: self.c0 * *c1,
            c2: self.c1 * *c1,
        }
    }

    pub fn mul_by_01(&self, c0: &Fp2<F>, c1: &Fp2<F>) -> Fp6<F> {
        let a_a = self.c0 * *c0;
        let b_b = self.c1 * *c1;

        let t1 = (self.c2 * *c1).mul_by_nonresidue() + a_a;

        let t2 = (*c0 + *c1) * (self.c0 + self.c1) - a_a - b_b;

        let t3 = self.c2 * *c0 + b_b;

        Fp6 {
            c0: t1,
            c1: t2,
            c2: t3,
        }
    }

    /// Multiply by quadratic nonresidue v.
    pub fn mul_by_nonresidue(&self) -> Self {
        // Given a + bv + cv^2, this produces
        //     av + bv^2 + cv^3
        // but because v^3 = u + 1, we have
        //     c(u + 1) + av + v^2

        Fp6 {
            c0: self.c2.mul_by_nonresidue(),
            c1: self.c0,
            c2: self.c1,
        }
    }
    /// Raises this element to p.
    #[inline(always)]
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();
        let c2 = self.c2.frobenius_map();

        // c1 = c1 * (u + 1)^((p - 1) / 3)
        let c1_coeffs = F::get_fp6_frobenius_coeff(1);
        let c1 = c1 * Fp2::new(c1_coeffs.0, c1_coeffs.1);

        // c2 = c2 * (u + 1)^((2p - 2) / 3)
        let c2_coeffs = F::get_fp6_frobenius_coeff(2);
        let c2 = c2 * Fp2::new(c2_coeffs.0, c2_coeffs.1);

        Fp6 { c0, c1, c2 }
    }

    pub(crate) fn nth_frobenius_map(&self, pow: usize) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();
        let c2 = self.c2.frobenius_map();

        // c1 = c1 * (u + 1)^((p - 1) / 3)
        let c1_coeffs = F::get_fp6_frobenius_coeff(pow);
        let c1 = c1 * Fp2::new(c1_coeffs.0, c1_coeffs.1);

        // c2 = c2 * (u + 1)^((2p - 2) / 3)
        let c2_coeffs = F::get_fp6_frobenius_coeff(pow + 1);
        let c2 = c2 * Fp2::new(c2_coeffs.0, c2_coeffs.1);

        Fp6 { c0, c1, c2 }
    }

    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() & self.c1.is_zero() & self.c2.is_zero()
    }

    #[inline(always)]
    pub fn is_one(&self) -> bool {
        self.c0.is_one() & self.c1.is_zero() & self.c2.is_zero()
    }

    #[inline]
    pub fn mul_interleaved(&self, rhs: &Self) -> Self {
        let t0 = self.c0 * rhs.c0;
        let t1 = self.c1 * rhs.c1;
        let t2 = self.c2 * rhs.c2;
        let c0 = self.c1 + self.c2;
        let tmp = rhs.c1 + rhs.c2;
        let c0 = c0 * tmp;
        let tmp = t2 + t1;
        let c0 = c0 - tmp;
        let c0 = c0.mul_by_nonresidue();
        let c0 = c0 + t0;
        let c1 = self.c0 + self.c1;
        let tmp = rhs.c0 + rhs.c1;
        let c1 = c1 * tmp;
        let tmp = t0 + t1;
        let c1 = c1 - tmp;
        let tmp = t2.mul_by_nonresidue();
        let c1 = c1 + tmp;
        let tmp = self.c0 + self.c2;
        let c2 = rhs.c0 + rhs.c2;
        let c2 = c2 * tmp;
        let tmp = t0 + t2;
        let c2 = c2 - tmp;
        let c2 = c2 + t1;
        Fp6::new(c0, c1, c2)
    }

    #[inline]
    pub fn square(&self) -> Self {
        let s0 = self.c0.square();
        let ab = self.c0 * self.c1;
        let s1 = ab + ab;
        let s2 = (self.c0 - self.c1 + self.c2).square();
        let bc = self.c1 * self.c2;
        let s3 = bc + bc;
        let s4 = self.c2.square();

        Fp6 {
            c0: s3.mul_by_nonresidue() + s0,
            c1: s4.mul_by_nonresidue() + s1,
            c2: s1 + s2 + s3 - s0 - s4,
        }
    }

    #[inline]
    pub fn invert(&self) -> Option<Self> {
        let c0 = (self.c1 * self.c2).mul_by_nonresidue();
        let c0 = self.c0.square() - c0;

        let c1 = self.c2.square().mul_by_nonresidue();
        let c1 = c1 - (self.c0 * self.c1);

        let c2 = self.c1.square();
        let c2 = c2 - (self.c0 * self.c2);

        let tmp = ((self.c1 * c2) + (self.c2 * c1)).mul_by_nonresidue();
        let tmp = tmp + (self.c0 * c0);

        let out = tmp.invert().map(|t| Fp6 {
            c0: t * c0,
            c1: t * c1,
            c2: t * c2,
        });
        println!("out: {:?}", out);
        out
    }
}

impl<F: Fp6Element> Mul<Fp6<F>> for Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn mul(self, other: Fp6<F>) -> Self::Output {
        self.mul_interleaved(&other)
    }
}

impl<F: Fp6Element> Add<Fp6<F>> for Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn add(self, rhs: Fp6<F>) -> Self::Output {
        Fp6 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
            c2: self.c2 + rhs.c2,
        }
    }
}

impl<'a, F: Fp6Element> Neg for &'a Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp6 {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}

impl<F: Fp6Element> Neg for Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F: Fp6Element> Sub<Fp6<F>> for Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn sub(self, rhs: Fp6<F>) -> Fp6<F> {
        Fp6 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
            c2: self.c2 - rhs.c2,
        }
    }
}

impl<F: Fp6Element> Mul<F> for Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn mul(self, rhs: F) -> Fp6<F> {
        Fp6 {
            c0: self.c0 * rhs,
            c1: self.c1 * rhs,
            c2: self.c2 * rhs,
        }
    }
}

impl<F: Fp6Element> Div<Fp6<F>> for Fp6<F> {
    type Output = Fp6<F>;

    #[inline]
    fn div(self, rhs: Fp6<F>) -> Fp6<F> {
        self * rhs.invert().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    macro_rules! fp6_tests {
        ($curve:ident, $rand_fn:ident, $curve_test: ident) => {
            mod $curve_test {
                use super::*;

                #[test]
                fn test_equality() {
                    let rng = &mut rand::thread_rng();
                    for _ in 0..10 {
                        let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
                        let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
                        let a = Fp6::<$curve>::new(
                            Fp2::new(
                                $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                                $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                            ),
                            Fp2::new(
                                $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                                $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                            ),
                            Fp2::new(
                                $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                                $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                            ),
                        );
                        let b = Fp6::<$curve>::new(
                            Fp2::new(
                                $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                                $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                            ),
                            Fp2::new(
                                $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                                $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                            ),
                            Fp2::new(
                                $curve::from_raw_unchecked(x.clone().try_into().unwrap()),
                                $curve::from_raw_unchecked(y.clone().try_into().unwrap()),
                            ),
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
                        assert_eq!(a + Fp6::<$curve>::zero(), a);
                        assert_eq!(a - Fp6::<$curve>::zero(), a);

                        assert_eq!(Fp6::<$curve>::zero() - a, -a);
                        assert_eq!(a - b, a + (-b));
                        assert_eq!(a - b, a + (b * -Fp6::<$curve>::one()));

                        assert_eq!(-a, Fp6::<$curve>::zero() - a);
                        assert_eq!(-a, a * -Fp6::<$curve>::one());
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

                        assert_eq!(a * Fp6::<$curve>::zero(), Fp6::<$curve>::zero());
                        assert_eq!(a * Fp6::<$curve>::zero(), Fp6::<$curve>::zero());
                        assert_eq!(a * Fp6::<$curve>::one(), a);
                        assert_eq!(a * Fp6::<$curve>::one(), a);
                        assert_eq!(a * Fp6::<$curve>::from(2u64), a + a);
                        assert_eq!(a * Fp6::<$curve>::from(3u64), a + a + a);
                        assert_eq!(a * Fp6::<$curve>::from(4u64), a + a + a + a);
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
                        assert_eq!(a / Fp6::<$curve>::one(), a);
                        assert_eq!(a / a, Fp6::<$curve>::one());

                        // division by zero
                        assert_eq!(Fp6::<$curve>::zero() / a, Fp6::<$curve>::zero());

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
                            assert_eq!(a * a.invert().unwrap(), Fp6::<$curve>::one());
                            assert_eq!(a.invert().unwrap().invert().unwrap(), a);
                        }
                    }
                }

                #[test]
                fn test_frobenius() {
                    for _ in 0..10 {
                        let a = $rand_fn();
                        let b = (0..6).fold(a, |acc, _| acc.frobenius_map());
                        assert_eq!(a, b);
                    }
                }
            }
        };
    }

    fn bls12381_fp6_rand() -> Fp6<Bls12381> {
        let mut rng = rand::thread_rng();
        Fp6::new(
            Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
            Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
            Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
        )
    }

    fn bn254_fp6_rand() -> Fp6<Bn254> {
        let mut rng = rand::thread_rng();
        Fp6::new(
            Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
            Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
            Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
        )
    }

    fp6_tests!(Bls12381, bls12381_fp6_rand, bls12381_fp6_test);
    fp6_tests!(Bn254, bn254_fp6_rand, bn254_fp6_test);
}
