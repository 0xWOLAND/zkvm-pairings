use crate::common::Curve;
use crate::fp::*;
use crate::fp2::*;
use crate::fp6::*;

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use rand::RngCore;
#[cfg(feature = "pairings")]
use rand_core::RngCore;

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
pub struct Fp12<C: Curve> {
    pub c0: Fp6<C>,
    pub c1: Fp6<C>,
}

impl<C: Curve> From<Fp<C>> for Fp12<C> {
    fn from(f: Fp<C>) -> Fp12<C> {
        Fp12 {
            c0: Fp6::<C>::from(f),
            c1: Fp6::<C>::from(f),
        }
    }
}

impl<C: Curve> From<Fp2<C>> for Fp12<C> {
    fn from(f: Fp2<C>) -> Fp12<C> {
        Fp12 {
            c0: Fp6::<C>::from(f),
            c1: Fp6::<C>::zero(),
        }
    }
}

impl<C: Curve> From<Fp6<C>> for Fp12<C> {
    fn from(f: Fp6<C>) -> Fp12<C> {
        Fp12 {
            c0: f,
            c1: Fp6::<C>::zero(),
        }
    }
}

impl<C: Curve> Eq for Fp12<C> {}
impl<C: Curve> PartialEq for Fp12<C> {
    fn eq(&self, other: &Fp12<C>) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<C: Curve> Copy for Fp12<C> {}
impl<C: Curve> Clone for Fp12<C> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: Curve> Default for Fp12<C> {
    fn default() -> Self {
        Fp12::<C>::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<C: Curve> zeroize::DefaultIsZeroes for Fp12 {}

impl<C: Curve> fmt::Debug for Fp12<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*w", self.c0, self.c1)
    }
}

impl<C: Curve> Fp12<C> {
    #[inline]
    pub fn new(c0: Fp6<C>, c1: Fp6<C>) -> Self {
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

    pub fn from_bytes(bytes: &[u8; 576]) -> Fp12<C> {
        let c0 = Fp6::<C>::from_bytes(&bytes[..288].try_into().unwrap());
        let c1 = Fp6::<C>::from_bytes(&bytes[288..].try_into().unwrap());

        Fp12::<C>::new(c0, c1)
    }

    pub fn to_bytes(&self) -> [u8; 576] {
        let mut res = [0u8; 576];
        let c0 = self.c0.to_bytes();
        let c1 = self.c1.to_bytes();

        res[..288].copy_from_slice(&c0);
        res[288..].copy_from_slice(&c1);

        res
    }

    pub(crate) fn random(mut rng: impl RngCore) -> Self {
        Fp12 {
            c0: Fp6::<C>::random(&mut rng),
            c1: Fp6::<C>::random(&mut rng),
        }
    }

    pub fn mul_by_014(&self, c0: &Fp2<C>, c1: &Fp2<C>, c4: &Fp2<C>) -> Fp12<C> {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = c1 + c4;
        let c1 = self.c1 + self.c0;
        let c1 = c1.mul_by_01(c0, &o);
        let c1 = c1 - aa - bb;
        let c0 = bb;
        let c0 = c0.mul_by_nonresidue();
        let c0 = c0 + aa;

        Fp12 { c0, c1 }
    }

    pub fn mul_14_by_14(d0: &Fp2<C>, d1: &Fp2<C>, c0: &Fp2<C>, c1: &Fp2<C>) -> [Fp2<C>; 5] {
        let x0 = d0 * c0;
        let x1 = d1 * c1;
        let x04 = c0 + d0;
        let tmp = c0 + c1;
        let x01 = d0 + d1;
        let x01 = x01 * tmp;
        let tmp = x1 + x0;
        let x01 = x01 - tmp;
        let x14 = c1 + d1;
        let z_c0_b0 = Fp2::<C>::non_residue() + x0;

        [z_c0_b0, x01, x1, x04, x14]
    }

    fn cyclotomic_square(&self) -> Fp12<C> {
        let t0 = &self.c1.c1.square();
        let t1 = &self.c0.c0.square();
        let t6 = (&self.c1.c1 + &self.c0.c0).square();
        let t6 = t6 - t0;
        let t6 = t6 - t1;
        let t2 = &self.c0.c2.square();
        let t3 = &self.c1.c0.square();
        let t7 = (&self.c0.c2 + &self.c1.c0).square();
        let t7 = t7 - t2;
        let t7 = t7 - t3;
        let t4 = &self.c1.c2.square();
        let t5 = &self.c0.c1.square();
        let t8 = (&self.c1.c2 + &self.c0.c1).square();
        let t8 = t8 - t4;
        let t8 = t8 - t5;
        let t8 = t8.mul_by_nonresidue();
        let t0 = t0.mul_by_nonresidue();
        let t0 = t0 + t1;
        let t2 = t2.mul_by_nonresidue();
        let t2 = t2 + t3;
        let t4 = t4.mul_by_nonresidue();
        let t4 = t4 + t5;
        let z00 = t0 - &self.c0.c0;
        let z00 = z00 + z00;
        let z00 = z00 + t0;
        let z01 = t2 - &self.c0.c1;
        let z01 = z01 + z01;
        let z01 = z01 + t2;
        let z02 = t4 - &self.c0.c2;
        let z02 = z02 + z02;
        let z02 = z02 + t4;
        let z10 = t8 + &self.c1.c0;
        let z10 = z10 + z10;
        let z10 = z10 + t8;
        let z11 = t6 + &self.c1.c1;
        let z11 = z11 + z11;
        let z11 = z11 + t6;
        let z12 = t7 + &self.c1.c2;
        let z12 = z12 + z12;
        let z12 = z12 + t7;
        Fp12::new(Fp6::new(z00, z01, z02), Fp6::new(z10, z11, z12))
    }

    fn n_cyclotomic_square(&self, by: u64) -> Fp12<C> {
        (0..by).fold(*self, |acc, _| acc.cyclotomic_square())
    }

    pub fn powt(&self) -> Fp12<C> {
        let a = self.cyclotomic_square();
        let a = a * self;
        let a = a.n_cyclotomic_square(2);
        let a = a * self;
        let a = a.n_cyclotomic_square(3);
        let a = a * self;
        let a = a.n_cyclotomic_square(9);
        let a = a * self;
        let a = a.n_cyclotomic_square(32);
        let a = a * self;
        let a = a.n_cyclotomic_square(15);
        let a = a * self;
        a.cyclotomic_square()
    }

    pub fn div(&self, rhs: &Fp12<C>) -> Fp12<C> {
        rhs.invert().unwrap() * self
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
                    res *= self;
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
                    res *= self;
                }
            }
        }
        res
    }

    /// Raises this element to p.
    #[inline(always)]
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();

        // c1 = c1 * (u + 1)^((p - 1) / 6)
        let c1 = c1
            * Fp6::from(Fp2::new(
                Fp::from_raw_unchecked([
                    0x8d07_75ed_9223_5fb8,
                    0xf67e_a53d_63e7_813d,
                    0x7b24_43d7_84ba_b9c4,
                    0x0fd6_03fd_3cbd_5f4f,
                    0xc231_beb4_202c_0d1f,
                    0x1904_d3bf_02bb_0667,
                ]),
                Fp::from_raw_unchecked([
                    0x2cf7_8a12_6ddc_4af3,
                    0x282d_5ac1_4d6c_7ec2,
                    0xec0c_8ec9_71f6_3c5f,
                    0x54a1_4787_b6c7_b36f,
                    0x88e9_e902_231f_9fb8,
                    0x00fc_3e2b_36c4_e032,
                ]),
            ));

        Fp12::new(c0, c1)
        // self.pow_vartime(&MODULUS)
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

impl<'a, 'b, C: Curve> Mul<&'b Fp12<C>> for &'a Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn mul(self, other: &'b Fp12<C>) -> Self::Output {
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

impl<'a, 'b, C: Curve> Add<&'b Fp12<C>> for &'a Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn add(self, rhs: &'b Fp12<C>) -> Self::Output {
        Fp12::new(self.c0 + rhs.c0, self.c1 + rhs.c1)
    }
}

impl<'a, C: Curve> Neg for &'a Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp12::new(-self.c0, -self.c1)
    }
}

impl<C: Curve> Neg for Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, 'b, C: Curve> Sub<&'b Fp12<C>> for &'a Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn sub(self, rhs: &'b Fp12<C>) -> Self::Output {
        Fp12::new(self.c0 - rhs.c0, self.c1 - rhs.c1)
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fp<C>> for &'a Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn mul(self, rhs: &'b Fp<C>) -> Fp12<C> {
        Fp12::new(self.c0 * rhs, self.c1 * rhs)
    }
}

impl<'a, 'b, C: Curve> Div<&'b Fp12<C>> for &'a Fp12<C> {
    type Output = Fp12<C>;

    #[inline]
    fn div(self, rhs: &'b Fp12<C>) -> Fp12<C> {
        self.div(rhs)
    }
}

impl_binops_additive!(Fp12<C>, Fp12<C>);
impl_binops_multiplicative!(Fp12<C>, Fp12<C>);
impl_binops_multiplicative!(Fp12<C>, Fp<C>);
impl_binops_divisible!(Fp12<C>, Fp12<C>);
#[cfg(test)]
mod test {
    use rand::thread_rng;

    use crate::common::Bls12381Curve;

    use super::*;

    fn fp12_rand() -> Fp12<Bls12381Curve> {
        Fp12::new(
            Fp6::new(
                Fp2::random(&mut thread_rng()),
                Fp2::random(&mut thread_rng()),
                Fp2::random(&mut thread_rng()),
            ),
            Fp6::new(
                Fp2::random(&mut thread_rng()),
                Fp2::random(&mut thread_rng()),
                Fp2::random(&mut thread_rng()),
            ),
        )
    }

    #[test]
    fn test_inequality() {
        for _ in 0..10 {
            let a = fp12_rand();
            let b = fp12_rand();

            assert_ne!(a, b)
        }
    }

    #[test]
    fn test_addition_subtraction() {
        for _ in 0..10 {
            let a = fp12_rand();
            let b = fp12_rand();
            let c = fp12_rand();

            // commutative
            assert_eq!(a + b, b + a);
            assert_eq!(a + (b + c), (a + b) + c);

            // additive identity
            assert_eq!(a + Fp12::zero(), a); // a + 0 = a
            assert_eq!(a - Fp12::zero(), a); // subtraction identity

            assert_eq!(Fp12::zero() - a, -a); // 0 - a = -a
            assert_eq!(a - b, a + (-b)); // a - b = a + -b
            assert_eq!(a - b, a + (b * -Fp12::one())); // a - b = a + b * -1

            assert_eq!(-a, Fp12::zero() - a);
            assert_eq!(-a, a * -Fp12::one());
        }
    }

    #[test]
    fn test_multiplication() {
        for _ in 0..10 {
            let a = fp12_rand();
            let b = fp12_rand();
            let c = fp12_rand();

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
            let a = fp12_rand();

            assert_eq!(a * Fp::from(0), Fp12::zero());
            assert_eq!(a * Fp::zero(), Fp12::zero());
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
            let a = fp12_rand();
            assert_eq!(a.square(), a * a);
        }
    }

    #[test]
    fn test_pow_equality() {
        for _ in 0..10 {
            let a = fp12_rand();
            assert_eq!(a.pow_vartime(&[1, 0, 0, 0, 0, 0]), a);
            assert_eq!(a.pow_vartime(&[2, 0, 0, 0, 0, 0]), a.square());
            assert_eq!(a.pow_vartime(&[3, 0, 0, 0, 0, 0]), a.square() * a);
            assert_eq!(a.pow_vartime(&[4, 0, 0, 0, 0, 0]), a.square().square());
        }
    }

    #[test]
    fn test_sqrt() {
        // let sqr1 = Fp([300855555557, 0, 0, 0, 0, 0]).sqrt().unwrap();
        // assert_eq!(format!("{:?}", sqr1), "0x025e51146a92917731d9d66d63f8c24ed8cae114e7c9d188e3eaa1e79bb19769f5877f9443e03723d9ed1eebbf92df98");

        // assert!(Fp([72057594037927816, 0, 0, 0, 0, 0]).sqrt().is_err());
    }

    #[test]
    fn test_div() {
        for _ in 0..10 {
            let a = fp12_rand();

            // division by one
            assert_eq!(a / Fp12::one(), a);
            assert_eq!(a / a, Fp12::one());

            // division by zero
            assert_eq!(Fp12::zero() / a, Fp12::zero());

            // division distributivity
            let a = fp12_rand();
            let b = fp12_rand();
            let c = fp12_rand();

            assert_eq!((a + b) / c, a / c + b / c);

            // division and multiplication equality
            let a = fp12_rand();
            let b = fp12_rand();
            assert_eq!(a / b, a * b.invert().unwrap());
        }
    }

    #[test]
    fn test_arithmetic() {
        use crate::fp::*;
        use crate::fp2::*;

        let a = Fp12::new(
            Fp6::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                ),
            ),
            Fp6::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                ),
            ),
        );

        let b = Fp12::new(
            Fp6::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d272_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e348,
                    ]),
                    Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                ),
            ),
            Fp6::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd2_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a117_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                ),
            ),
        );
        // ... (previous definitions of a and b remain the same)

        let c = Fp12::new(
            Fp6::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x47f9_cb98_71b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0x7791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_133c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_40e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_1744_c040,
                    ]),
                ),
            ),
            Fp6::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d3_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_1040,
                    ]),
                ),
            ),
        );

        // because a and b and c are similar to each other and
        // I was lazy, this is just some arbitrary way to make
        // them a little more different
        let a: Fp12<Bls12381Curve> = a.square().invert().unwrap().square() + c;
        let b: Fp12<Bls12381Curve> = b.square().invert().unwrap().square() + a;
        let c: Fp12<Bls12381Curve> = c.square().invert().unwrap().square() + b;

        assert_eq!(a.square(), a * a);
        assert_eq!(b.square(), b * b);
        assert_eq!(c.square(), c * c);

        assert_eq!((a + b) * c.square(), (c * c * a) + (c * c * b));

        assert_eq!(
            a.invert().unwrap() * b.invert().unwrap(),
            (a * b).invert().unwrap()
        );
        assert_eq!(a.invert().unwrap() * a, Fp12::one());

        assert!(a != a.frobenius_map());
        assert_eq!(
            a,
            a.frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
        );
    }

    #[test]
    fn test_cyclotomic_square() {
        for _ in 0..10 {
            let a = fp12_rand();
            assert_eq!(a.cyclotomic_square(), a.n_cyclotomic_square(1));
            assert_eq!(
                a.cyclotomic_square().cyclotomic_square(),
                a.n_cyclotomic_square(2)
            );
            assert_eq!(
                a.cyclotomic_square()
                    .cyclotomic_square()
                    .cyclotomic_square(),
                a.n_cyclotomic_square(3)
            );
            assert_eq!(
                a.cyclotomic_square()
                    .cyclotomic_square()
                    .cyclotomic_square()
                    .cyclotomic_square(),
                a.n_cyclotomic_square(4)
            );
        }
    }

    #[test]
    fn test_pow_vartime() {
        for _ in 0..10 {
            let a = fp12_rand();
            let exp = (0..6).map(|_| rand::random::<u64>()).collect::<Vec<_>>();
            let lhs = a.pow_vartime(&exp.clone().try_into().unwrap());
            let rhs = a.pow_vartime_extended(exp.as_slice());

            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn test_from_bytes() {
        // let bytes = Fp::<Bls12381Curve>::one().to_bytes_unsafe();
        // println!("{:?}", bytes);
        assert_eq!(
            Fp2::<Bls12381Curve>::from_bytes(&Fp2::<Bls12381Curve>::to_bytes(
                &Fp2::<Bls12381Curve>::one()
            )),
            Fp2::<Bls12381Curve>::one()
        );
    }
}

#[cfg(feature = "zeroize")]
#[test]
fn test_zeroize() {
    use zeroize::Zeroize;

    let mut a = Fp12::one();
    a.zeroize();
    assert!(bool::from(a.is_zero()));
}
