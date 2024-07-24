use crate::fp::*;
use crate::fp2::*;
use crate::fp6::*;

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[cfg(feature = "pairings")]
use rand_core::RngCore;

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl From<Fp> for Fp12 {
    fn from(f: Fp) -> Fp12 {
        Fp12 {
            c0: Fp6::from(f),
            c1: Fp6::zero(),
        }
    }
}

impl From<Fp2> for Fp12 {
    fn from(f: Fp2) -> Fp12 {
        Fp12 {
            c0: Fp6::from(f),
            c1: Fp6::zero(),
        }
    }
}

impl From<Fp6> for Fp12 {
    fn from(f: Fp6) -> Fp12 {
        Fp12 {
            c0: f,
            c1: Fp6::zero(),
        }
    }
}

impl Eq for Fp12 {}
impl PartialEq for Fp12 {
    fn eq(&self, other: &Fp12) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl Copy for Fp12 {}
impl Clone for Fp12 {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl Default for Fp12 {
    fn default() -> Self {
        Fp12::zero()
    }
}

#[cfg(feature = "zeroize")]
impl zeroize::DefaultIsZeroes for Fp12 {}

impl fmt::Debug for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*w", self.c0, self.c1)
    }
}

impl Fp12 {
    #[inline]
    pub fn zero() -> Self {
        Fp12 {
            c0: Fp6::zero(),
            c1: Fp6::zero(),
        }
    }

    #[inline]
    pub fn one() -> Self {
        Fp12 {
            c0: Fp6::one(),
            c1: Fp6::zero(),
        }
    }

    #[cfg(feature = "pairings")]
    pub(crate) fn random(mut rng: impl RngCore) -> Self {
        Fp12 {
            c0: Fp6::random(&mut rng),
            c1: Fp6::random(&mut rng),
        }
    }

    pub fn mul_by_014(&self, c0: &Fp2, c1: &Fp2, c4: &Fp2) -> Fp12 {
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

    pub fn div(&self, rhs: &Fp12) -> Fp12 {
        rhs.invert().unwrap() * self
    }

    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    #[inline(always)]
    pub fn conjugate(&self) -> Self {
        Fp12 {
            c0: self.c0,
            c1: -self.c1,
        }
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

    /// Raises this element to p.
    #[inline(always)]
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();

        // c1 = c1 * (u + 1)^((p - 1) / 6)
        let c1 = c1
            * Fp6::from(Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x8d07_75ed_9223_5fb8,
                    0xf67e_a53d_63e7_813d,
                    0x7b24_43d7_84ba_b9c4,
                    0x0fd6_03fd_3cbd_5f4f,
                    0xc231_beb4_202c_0d1f,
                    0x1904_d3bf_02bb_0667,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x2cf7_8a12_6ddc_4af3,
                    0x282d_5ac1_4d6c_7ec2,
                    0xec0c_8ec9_71f6_3c5f,
                    0x54a1_4787_b6c7_b36f,
                    0x88e9_e902_231f_9fb8,
                    0x00fc_3e2b_36c4_e032,
                ]),
            });

        Fp12 { c0, c1 }
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

        Fp12 { c0, c1 }
    }

    pub fn invert(&self) -> Option<Self> {
        (self.c0.square() - self.c1.square().mul_by_nonresidue())
            .invert()
            .map(|t| Fp12 {
                c0: self.c0 * t,
                c1: self.c1 * -t,
            })
    }
}

impl<'a, 'b> Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn mul(self, other: &'b Fp12) -> Self::Output {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;
        let o = other.c0 + other.c1;
        let c1 = self.c1 + self.c0;
        let c1 = c1 * o;
        let c1 = c1 - aa;
        let c1 = c1 - bb;
        let c0 = bb.mul_by_nonresidue();
        let c0 = c0 + aa;

        Fp12 { c0, c1 }
    }
}

impl<'a, 'b> Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn add(self, rhs: &'b Fp12) -> Self::Output {
        Fp12 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<'a> Neg for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp12 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

impl Neg for Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn sub(self, rhs: &'b Fp12) -> Self::Output {
        Fp12 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl<'a, 'b> Mul<&'b Fp> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn mul(self, rhs: &'b Fp) -> Fp12 {
        Fp12 {
            c0: self.c0 * rhs,
            c1: self.c1 * rhs,
        }
    }
}

impl<'a, 'b> Div<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn div(self, rhs: &'b Fp12) -> Fp12 {
        self.div(rhs)
    }
}

impl_binops_additive!(Fp12, Fp12);
impl_binops_multiplicative!(Fp12, Fp12);
impl_binops_multiplicative!(Fp12, Fp);
impl_binops_divisible!(Fp12, Fp12);
#[cfg(test)]
mod test {
    use rand::{thread_rng, Rng};

    use super::*;

    fn fp12_rand() -> Fp12 {
        Fp12 {
            c0: Fp6 {
                c0: Fp2::random(&mut thread_rng()),
                c1: Fp2::random(&mut thread_rng()),
                c2: Fp2::random(&mut thread_rng()),
            },
            c1: Fp6 {
                c0: Fp2::random(&mut thread_rng()),
                c1: Fp2::random(&mut thread_rng()),
                c2: Fp2::random(&mut thread_rng()),
            },
        }
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp::from_raw_unchecked(x.clone().try_into().unwrap());
            let b = Fp::from_raw_unchecked(x.try_into().unwrap());

            assert_eq!(a, b)
        }
    }

    #[test]
    fn test_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp::from_raw_unchecked(x.try_into().unwrap());
            let b = Fp::from_raw_unchecked(y.try_into().unwrap());

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
        let sqr1 = Fp([300855555557, 0, 0, 0, 0, 0]).sqrt().unwrap();
        assert_eq!(format!("{:?}", sqr1), "0x025e51146a92917731d9d66d63f8c24ed8cae114e7c9d188e3eaa1e79bb19769f5877f9443e03723d9ed1eebbf92df98");

        assert!(Fp([72057594037927816, 0, 0, 0, 0, 0]).sqrt().is_err());
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

        let a = Fp12 {
            c0: Fp6 {
                c0: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                },
                c1: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                },
                c2: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                },
            },
            c1: Fp6 {
                c0: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                },
                c1: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                },
                c2: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                },
            },
        };

        let b = Fp12 {
            c0: Fp6 {
                c0: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d272_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                },
                c1: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e348,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                },
                c2: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                },
            },
            c1: Fp6 {
                c0: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd2_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                },
                c1: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a117_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                },
                c2: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_c040,
                    ]),
                },
            },
        };

        let c = Fp12 {
            c0: Fp6 {
                c0: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x47f9_cb98_71b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                },
                c1: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0x7791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d9_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_133c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                },
                c2: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_40e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_1744_c040,
                    ]),
                },
            },
            c1: Fp6 {
                c0: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x47f9_cb98_b1b8_2d58,
                        0x5fe9_11eb_a3aa_1d9d,
                        0x96bf_1b5f_4dd8_1db3,
                        0x8100_d27c_c925_9f5b,
                        0xafa2_0b96_7464_0eab,
                        0x09bb_cea7_d8d9_497d,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x0303_cb98_b166_2daa,
                        0xd931_10aa_0a62_1d5a,
                        0xbfa9_820c_5be4_a468,
                        0x0ba3_643e_cb05_a348,
                        0xdc35_34bb_1f1c_25a6,
                        0x06c3_05bb_19c0_e1c1,
                    ]),
                },
                c1: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x46f9_cb98_b162_d858,
                        0x0be9_109c_f7aa_1d57,
                        0xc791_bc55_fece_41d2,
                        0xf84c_5770_4e38_5ec2,
                        0xcb49_c1d3_c010_e60f,
                        0x0acd_b8e1_58bf_e3c8,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x8aef_cb98_b15f_8306,
                        0x3ea1_108f_e4f2_1d54,
                        0xcf79_f69f_a1b7_df3b,
                        0xe4f5_4aa1_d16b_1a3c,
                        0xba5e_4ef8_6105_a679,
                        0x0ed8_6c07_97be_e5cf,
                    ]),
                },
                c2: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xcee5_cb98_b15c_2db4,
                        0x7159_1082_d23a_1d51,
                        0xd762_30e9_44a1_7ca4,
                        0xd19e_3dd3_549d_d5b6,
                        0xa972_dc17_01fa_66e3,
                        0x12e3_1f2d_d6bd_e7d6,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xad2a_cb98_b173_2d9d,
                        0x2cfd_10dd_0696_1d64,
                        0x0739_6b86_c6ef_24e8,
                        0xbd76_e2fd_b1bf_c820,
                        0x6afe_a7f6_de94_d0d5,
                        0x1099_4b0c_5744_1040,
                    ]),
                },
            },
        };

        // because a and b and c are similar to each other and
        // I was lazy, this is just some arbitrary way to make
        // them a little more different
        let a = a.square().invert().unwrap().square() + c;
        let b = b.square().invert().unwrap().square() + a;
        let c = c.square().invert().unwrap().square() + b;

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
}

#[cfg(feature = "zeroize")]
#[test]
fn test_zeroize() {
    use zeroize::Zeroize;

    let mut a = Fp12::one();
    a.zeroize();
    assert!(bool::from(a.is_zero()));
}
