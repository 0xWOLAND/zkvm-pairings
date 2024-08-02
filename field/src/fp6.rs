use crate::common::Curve;
use crate::fp::*;
use crate::fp2::*;

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use rand::RngCore;
#[cfg(feature = "pairings")]
use rand_core::RngCore;

/// This represents an element $c_0 + c_1 v + c_2 v^2$ of $\mathbb{F}_{p^6} = \mathbb{F}_{p^2} / v^3 - u - 1$.
pub struct Fp6<C: Curve> {
    pub c0: Fp2<C>,
    pub c1: Fp2<C>,
    pub c2: Fp2<C>,
}

impl<C: Curve> From<Fp<C>> for Fp6<C> {
    fn from(f: Fp<C>) -> Fp6<C> {
        Fp6 {
            c0: Fp2::from(f),
            c1: Fp2::from(f),
            c2: Fp2::from(f),
        }
    }
}

impl<C: Curve> From<Fp2<C>> for Fp6<C> {
    fn from(f: Fp2<C>) -> Fp6<C> {
        Fp6 {
            c0: f,
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }
}

impl<C: Curve> PartialEq for Fp6<C> {
    fn eq(&self, other: &Fp6<C>) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}

impl<C: Curve> Copy for Fp6<C> {}
impl<C: Curve> Clone for Fp6<C> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: Curve> Default for Fp6<C> {
    fn default() -> Self {
        Fp6::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<C: Curve> zeroize::DefaultIsZeroes for Fp6<C> {}

impl<C: Curve> fmt::Debug for Fp6<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*v + ({:?})*v^2", self.c0, self.c1, self.c2)
    }
}

impl<C: Curve> Eq for Fp6<C> {}

impl<C: Curve> Fp6<C> {
    #[inline]
    pub fn new(c0: Fp2<C>, c1: Fp2<C>, c2: Fp2<C>) -> Self {
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

    pub(crate) fn random(mut rng: impl RngCore) -> Fp6<C> {
        Fp6 {
            c0: Fp2::random(&mut rng),
            c1: Fp2::random(&mut rng),
            c2: Fp2::random(&mut rng),
        }
    }

    pub fn from_bytes(bytes: &[u8; 288]) -> Fp6<C> {
        let mut c0_bytes = [0u8; 96];
        let mut c1_bytes = [0u8; 96];
        let mut c2_bytes = [0u8; 96];

        c0_bytes.copy_from_slice(&bytes[0..96]);
        c1_bytes.copy_from_slice(&bytes[96..192]);
        c2_bytes.copy_from_slice(&bytes[192..288]);

        let c0 = Fp2::from_bytes(&c0_bytes);
        let c1 = Fp2::from_bytes(&c1_bytes);
        let c2 = Fp2::from_bytes(&c2_bytes);

        Fp6::<C>::new(c0, c1, c2)
    }

    pub fn to_bytes(&self) -> [u8; 288] {
        let mut res = [0u8; 288];
        let c0_bytes = self.c0.to_bytes();
        let c1_bytes = self.c1.to_bytes();
        let c2_bytes = self.c2.to_bytes();

        res[0..96].copy_from_slice(&c0_bytes);
        res[96..192].copy_from_slice(&c1_bytes);
        res[192..288].copy_from_slice(&c2_bytes);

        res
    }

    pub fn mul_by_1(&self, c1: &Fp2<C>) -> Fp6<C> {
        Fp6 {
            c0: (self.c2 * c1).mul_by_nonresidue(),
            c1: self.c0 * c1,
            c2: self.c1 * c1,
        }
    }

    pub fn mul_by_01(&self, c0: &Fp2<C>, c1: &Fp2<C>) -> Fp6<C> {
        let a_a = self.c0 * c0;
        let b_b = self.c1 * c1;

        let t1 = (self.c2 * c1).mul_by_nonresidue() + a_a;

        let t2 = (c0 + c1) * (self.c0 + self.c1) - a_a - b_b;

        let t3 = self.c2 * c0 + b_b;

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
        let c1 = c1
            * Fp2::new(
                Fp::from_raw_unchecked([
                    0x2e01_ffff_fffe_fffe,
                    0xde17_d813_620a_0002,
                    0xddb3_a93b_e6f8_9688,
                    0xba69_c607_6a0f_77ea,
                    0x5f19_672f_df76_ce51,
                    0x0000_0000_0000_0000,
                ]),
                Fp::zero(),
            );

        // c2 = c2 * (u + 1)^((2p - 2) / 3)
        let c2 = c2
            * Fp2::new(
                Fp::from_raw_unchecked([
                    0x8bfd_0000_0000_aaac,
                    0x4094_27eb_4f49_fffd,
                    0x897d_2965_0fb8_5f9b,
                    0xaa0d_857d_8975_9ad4,
                    0xec02_4086_63d4_de85,
                    0x1a01_11ea_397f_e699,
                ]),
                Fp::zero(),
            );

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
        let t0 = &self.c0 * &rhs.c0;
        let t1 = &self.c1 * &rhs.c1;
        let t2 = &self.c2 * &rhs.c2;
        let c0 = &self.c1 + &self.c2;
        let tmp = &rhs.c1 + &rhs.c2;
        let c0 = c0 * tmp;
        let tmp = t2 + t1;
        let c0 = c0 - tmp;
        let c0 = c0.mul_by_nonresidue();
        let c0 = c0 + t0;
        let c1 = &self.c0 + &self.c1;
        let tmp = &rhs.c0 + &rhs.c1;
        let c1 = c1 * tmp;
        let tmp = t0 + t1;
        let c1 = c1 - tmp;
        let tmp = t2.mul_by_nonresidue();
        let c1 = c1 + tmp;
        let tmp = &self.c0 + &self.c2;
        let c2 = &rhs.c0 + &rhs.c2;
        let c2 = c2 * tmp;
        let tmp = t0 + t2;
        let c2 = c2 - tmp;
        let c2 = c2 + t1;
        Fp6::new(c0, c1, c2)
    }

    pub fn div(&self, rhs: &Self) -> Self {
        self * rhs.invert().unwrap()
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

        tmp.invert().map(|t| Fp6 {
            c0: t * c0,
            c1: t * c1,
            c2: t * c2,
        })
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fp6<C>> for &'a Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn mul(self, other: &'b Fp6<C>) -> Self::Output {
        self.mul_interleaved(other)
    }
}

impl<'a, 'b, C: Curve> Add<&'b Fp6<C>> for &'a Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn add(self, rhs: &'b Fp6<C>) -> Self::Output {
        Fp6 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
            c2: self.c2 + rhs.c2,
        }
    }
}

impl<'a, C: Curve> Neg for &'a Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp6 {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}

impl<C: Curve> Neg for Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, 'b, C: Curve> Sub<&'b Fp6<C>> for &'a Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn sub(self, rhs: &'b Fp6<C>) -> Fp6<C> {
        Fp6 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
            c2: self.c2 - rhs.c2,
        }
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fp<C>> for &'a Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn mul(self, rhs: &'b Fp<C>) -> Fp6<C> {
        Fp6 {
            c0: self.c0 * rhs,
            c1: self.c1 * rhs,
            c2: self.c2 * rhs,
        }
    }
}

impl<'a, 'b, C: Curve> Div<&'b Fp6<C>> for &'a Fp6<C> {
    type Output = Fp6<C>;

    #[inline]
    fn div(self, rhs: &'b Fp6<C>) -> Fp6<C> {
        self.div(rhs)
    }
}

impl_binops_additive!(Fp6<C>, Fp6<C>);
impl_binops_multiplicative!(Fp6<C>, Fp6<C>);
impl_binops_multiplicative!(Fp6<C>, Fp<C>);
impl_binops_divisible!(Fp6<C>, Fp6<C>);

#[cfg(test)]
mod test {
    use rand::{thread_rng, Rng};

    use crate::common::Bls12381Curve;

    use super::*;

    fn fp6_rand() -> Fp6<Bls12381Curve> {
        Fp6 {
            c0: Fp2::random(&mut thread_rng()),
            c1: Fp2::random(&mut thread_rng()),
            c2: Fp2::random(&mut thread_rng()),
        }
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp2::new(
                Fp::<Bls12381Curve>::from_raw_unchecked(x.clone().try_into().unwrap()),
                Fp::<Bls12381Curve>::from_raw_unchecked(y.clone().try_into().unwrap()),
            );

            let b = Fp2::new(
                Fp::<Bls12381Curve>::from_raw_unchecked(x.clone().try_into().unwrap()),
                Fp::<Bls12381Curve>::from_raw_unchecked(y.clone().try_into().unwrap()),
            );

            assert_eq!(a, b)
        }
    }

    #[test]
    fn test_inequality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x1 = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y1 = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let x2 = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y2 = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = Fp2::new(
                Fp::<Bls12381Curve>::from_raw_unchecked(x1.clone().try_into().unwrap()),
                Fp::<Bls12381Curve>::from_raw_unchecked(y1.clone().try_into().unwrap()),
            );
            let b = Fp2::new(
                Fp::<Bls12381Curve>::from_raw_unchecked(x2.try_into().unwrap()),
                Fp::<Bls12381Curve>::from_raw_unchecked(y2.try_into().unwrap()),
            );

            assert_ne!(a, b)
        }
    }

    #[test]
    fn test_addition_subtraction() {
        for _ in 0..10 {
            let a = fp6_rand();
            let b = fp6_rand();
            let c = fp6_rand();

            // commutative
            assert_eq!(a + b, b + a);
            assert_eq!(a + (b + c), (a + b) + c);

            // additive identity
            assert_eq!(a + Fp6::zero(), a); // a + 0 = a
            assert_eq!(a - Fp6::zero(), a); // subtraction identity

            assert_eq!(Fp6::zero() - a, -a); // 0 - a = -a
            assert_eq!(a - b, a + (-b)); // a - b = a + -b
            assert_eq!(a - b, a + (b * -Fp6::one())); // a - b = a + b * -1

            assert_eq!(-a, Fp6::zero() - a);
            assert_eq!(-a, a * -Fp6::one());
        }
    }

    #[test]
    fn test_multiplication() {
        for _ in 0..10 {
            let a = fp6_rand();
            let b = fp6_rand();
            let c = fp6_rand();

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
            let a = fp6_rand();

            let _a = Fp6::new(a.c0 * Fp::one(), a.c1 * Fp::one(), a.c2 * Fp::one());

            assert_eq!(a * Fp::from(0), Fp6::zero());
            assert_eq!(a * Fp::zero(), Fp6::zero());
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
            let a = fp6_rand();
            assert_eq!(a.square(), a * a);
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
                .is_none()
        );
    }

    #[test]
    fn test_div() {
        for _ in 0..10 {
            let a = fp6_rand();

            // division by one
            assert_eq!(a / Fp6::one(), a);
            assert_eq!(a / a, Fp6::one());

            // division by zero
            assert_eq!(Fp6::zero() / a, Fp6::zero());

            // division distributivity
            let a = fp6_rand();
            let b = fp6_rand();
            let c = fp6_rand();

            assert_eq!((a + b) / c, a / c + b / c);

            // division and multiplication equality
            let a = fp6_rand();
            let b = fp6_rand();
            assert_eq!(a / b, a * b.invert().unwrap());
        }
    }
    #[test]
    fn test_arithmetic() {
        use crate::fp::*;

        let a = Fp6::<Bls12381Curve>::new(
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
                    0x0ba3_643e_cc05_a348,
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
                    0x1099_4c0c_5744_c040,
                ]),
            ),
        );

        let b = Fp6::new(
            Fp2::new(
                Fp::from_raw_unchecked([
                    0xf120_cb98_b16f_d84b,
                    0x5fb5_10cf_f3de_1d61,
                    0x0f21_a5d0_69d8_c251,
                    0xaa1f_d62f_34f2_839a,
                    0x5a13_3515_7f89_913f,
                    0x14a3_fe32_9643_c247,
                ]),
                Fp::from_raw_unchecked([
                    0x3516_cb98_b16c_82f9,
                    0x926d_10c2_e126_1d5f,
                    0x1709_e01a_0cc2_5fba,
                    0x96c8_c960_b825_3f14,
                    0x4927_c234_207e_51a9,
                    0x18ae_b158_d542_c44e,
                ]),
            ),
            Fp2::new(
                Fp::from_raw_unchecked([
                    0xbf0d_cb98_b169_82fc,
                    0xa679_10b7_1d1a_1d5c,
                    0xb7c1_47c2_b8fb_06ff,
                    0x1efa_710d_47d2_e7ce,
                    0xed20_a79c_7e27_653c,
                    0x02b8_5294_dac1_dfba,
                ]),
                Fp::from_raw_unchecked([
                    0x9d52_cb98_b180_82e5,
                    0x621d_1111_5176_1d6f,
                    0xe798_8260_3b48_af43,
                    0x0ad3_1637_a4f4_da37,
                    0xaeac_737c_5ac1_cf2e,
                    0x006e_7e73_5b48_b824,
                ]),
            ),
            Fp2::new(
                Fp::from_raw_unchecked([
                    0xe148_cb98_b17d_2d93,
                    0x94d5_1104_3ebe_1d6c,
                    0xef80_bca9_de32_4cac,
                    0xf77c_0969_2827_95b1,
                    0x9dc1_009a_fbb6_8f97,
                    0x0479_3199_9a47_ba2b,
                ]),
                Fp::from_raw_unchecked([
                    0x253e_cb98_b179_d841,
                    0xc78d_10f7_2c06_1d6a,
                    0xf768_f6f3_811b_ea15,
                    0xe424_fc9a_ab5a_512b,
                    0x8cd5_8db9_9cab_5001,
                    0x0883_e4bf_d946_bc32,
                ]),
            ),
        );

        let c = Fp6::new(
            Fp2::new(
                Fp::from_raw_unchecked([
                    0x6934_cb98_b176_82ef,
                    0xfa45_10ea_194e_1d67,
                    0xff51_313d_2405_877e,
                    0xd0cd_efcc_2e8d_0ca5,
                    0x7bea_1ad8_3da0_106b,
                    0x0c8e_97e6_1845_be39,
                ]),
                Fp::from_raw_unchecked([
                    0x4779_cb98_b18d_82d8,
                    0xb5e9_1144_4daa_1d7a,
                    0x2f28_6bda_a653_2fc2,
                    0xbca6_94f6_8bae_ff0f,
                    0x3d75_e6b8_1a3a_7a5d,
                    0x0a44_c3c4_98cc_96a3,
                ]),
            ),
            Fp2::new(
                Fp::from_raw_unchecked([
                    0x8b6f_cb98_b18a_2d86,
                    0xe8a1_1137_3af2_1d77,
                    0x3710_a624_493c_cd2b,
                    0xa94f_8828_0ee1_ba89,
                    0x2c8a_73d6_bb2f_3ac7,
                    0x0e4f_76ea_d7cb_98aa,
                ]),
                Fp::from_raw_unchecked([
                    0xcf65_cb98_b186_d834,
                    0x1b59_112a_283a_1d74,
                    0x3ef8_e06d_ec26_6a95,
                    0x95f8_7b59_9214_7603,
                    0x1b9f_00f5_5c23_fb31,
                    0x125a_2a11_16ca_9ab1,
                ]),
            ),
            Fp2::new(
                Fp::from_raw_unchecked([
                    0x135b_cb98_b183_82e2,
                    0x4e11_111d_1582_1d72,
                    0x46e1_1ab7_8f10_07fe,
                    0x82a1_6e8b_1547_317d,
                    0x0ab3_8e13_fd18_bb9b,
                    0x1664_dd37_55c9_9cb8,
                ]),
                Fp::from_raw_unchecked([
                    0xce65_cb98_b131_8334,
                    0xc759_0fdb_7c3a_1d2e,
                    0x6fcb_8164_9d1c_8eb3,
                    0x0d44_004d_1727_356a,
                    0x3746_b738_a7d0_d296,
                    0x136c_144a_96b1_34fc,
                ]),
            ),
        );

        assert_eq!(a.square(), a * a);
        assert_eq!(b.square(), b * b);
        assert_eq!(c.square(), c * c);

        assert_eq!((a + b) * c.square(), (c * c * a) + (c * c * b));

        assert_eq!(
            a.invert().unwrap() * b.invert().unwrap(),
            (a * b).invert().unwrap()
        );
        assert_eq!(a.invert().unwrap() * a, Fp6::one());

        assert_eq!(
            a.frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map()
                .frobenius_map(),
            a
        )
    }

    #[cfg(feature = "zeroize")]
    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fp6::one();
        a.zeroize();
        assert!(bool::from(a.is_zero()));
    }
}
