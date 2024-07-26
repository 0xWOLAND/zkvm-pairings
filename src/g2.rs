use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::common::{AffinePoint, Curve};
use crate::fp::Fp;
use crate::{fp2::Fp2, fr::Fr};

#[derive(Clone, Copy, Debug)]
struct G2Affine<C: Curve> {
    x: Fp2<C>,
    y: Fp2<C>,
    is_infinity: bool,
    _marker: PhantomData<C>,
}

impl<C: Curve> PartialEq for G2Affine<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: Curve> AffinePoint<C> for G2Affine<C> {
    type Dtype = Fp2<C>;

    fn new(x: Self::Dtype, y: Self::Dtype, is_infinity: bool) -> Self {
        G2Affine {
            x,
            y,
            is_infinity,
            _marker: PhantomData::<C>,
        }
    }

    fn identity() -> Self {
        G2Affine {
            x: Fp2::zero(),
            y: Fp2::one(),
            is_infinity: true,
            _marker: PhantomData::<C>,
        }
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }

    fn generator() -> Self {
        G2Affine {
            x: Fp2::new(
                Fp::from_raw_unchecked(C::G2_X0),
                Fp::from_raw_unchecked(C::G2_X1),
            ),
            y: Fp2::new(
                Fp::from_raw_unchecked(C::G2_Y0),
                Fp::from_raw_unchecked(C::G2_Y1),
            ),
            is_infinity: false,
            _marker: PhantomData::<C>,
        }
    }

    fn is_valid(&self) -> Result<(), String> {
        if self.is_zero() {
            return Ok(());
        }
        if !self.is_on_curve() {
            return Err("Point is not on curve".to_string());
        }
        if !self.is_torsion_free() {
            return Err("Point is not torsion free".to_string());
        }

        Ok(())
    }

    fn random(mut rng: impl rand::Rng) -> Self {
        let x = Fp2::random(&mut rng);
        let y = Fp2::random(&mut rng);
        Self {
            x,
            y,
            is_infinity: false,
            _marker: PhantomData::<C>,
        }
    }

    fn double(&self) -> Self {
        let x = self.x;
        let y = self.y;

        if y.is_zero() {
            return Self::identity();
        }

        let three = Fp2::<C>::from(Fp::from(3));
        let two = Fp2::<C>::from(Fp::from(2));

        let slope = (three * x.square()) / (two * y);
        let xr = slope.square() - two * x;
        let yr = slope * (x - xr) - y;

        Self {
            x: xr,
            y: yr,
            is_infinity: false,
            _marker: PhantomData::<C>,
        }
    }
}

impl<C: Curve> G2Affine<C> {
    fn is_on_curve(&self) -> bool {
        let x = self.x;
        let y = self.y;

        // y^2 = x^3 + B
        y.square()
            == x.square() * x
                + Fp2::new(
                    Fp::from_raw_unchecked(C::B2_X),
                    Fp::from_raw_unchecked(C::B2_Y),
                )
    }

    fn mul_by_x(&self) -> Self {
        println!("C::X: {:?}", Fr::<C>::from(C::X));
        self * Fr::from(C::X)
    }

    fn psi(&self) -> Self {
        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2::<C>::new(
            Fp::zero(),
            Fp::from_raw_unchecked([
                0x8bfd00000000aaad,
                0x409427eb4f49fffd,
                0x897d29650fb85f9b,
                0xaa0d857d89759ad4,
                0xec02408663d4de85,
                0x1a0111ea397fe699,
            ]),
        );
        println!("psi_coeff_x: {:?}", psi_coeff_x.c1);
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fp2::<C>::new(
            Fp::from_raw_unchecked([
                0xf1ee7b04121bdea2,
                0x304466cf3e67fa0a,
                0xef396489f61eb45e,
                0x1c3dedd930b1cf60,
                0xe2e9c448d77a2cd9,
                0x135203e60180a68e,
            ]),
            Fp::from_raw_unchecked([
                0xc81084fbede3cc09,
                0xee67992f72ec05f4,
                0x77f76e17009241c5,
                0x48395dabc2d3435e,
                0x6831e36d6bd17ffe,
                0x06af0e0437ff400b,
            ]),
        );

        G2Affine::new(
            self.x.frobenius_map() * psi_coeff_x,
            self.y.frobenius_map() * psi_coeff_y,
            false,
        )
    }

    fn is_torsion_free(&self) -> bool {
        let lhs = self.psi();
        let rhs = -&self.mul_by_x();
        println!("lhs: {:?}", lhs);
        println!("rhs: {:?}", rhs);
        lhs == rhs
    }
}

impl<'a, C: Curve> Neg for &'a G2Affine<C> {
    type Output = G2Affine<C>;

    fn neg(self) -> G2Affine<C> {
        G2Affine {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
            _marker: PhantomData::<C>,
        }
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fr<C>> for &'a G2Affine<C> {
    type Output = G2Affine<C>;

    #[inline]
    fn mul(self, other: &'b Fr<C>) -> G2Affine<C> {
        let mut xself = G2Affine::<C>::identity();
        let mut acc = *self;
        let mut i = 0;

        for bit in other
            .0
            .iter()
            .flat_map(|&x| (0..64).map(move |i| (x >> i) & 1 == 1))
            .skip(1)
        {
            println!("ACC before: {:?}", acc);
            acc = acc.double();

            if bit {
                println!("({}) added: {:?}", i, acc);
                xself += &acc;
            }
            i += 1;
        }
        println!("iterations: {:?}", i);

        xself
    }
}

impl<'a, 'b, C: Curve> Add<&'b G2Affine<C>> for &'a G2Affine<C> {
    type Output = G2Affine<C>;

    #[inline]
    fn add(self, other: &'b G2Affine<C>) -> G2Affine<C> {
        if self.is_infinity {
            return *other;
        }

        if other.is_infinity {
            return *self;
        }

        let x1 = self.x;
        let y1 = self.y;
        let x2 = other.x;
        let y2 = other.y;

        if x1 == x2 && y1 == y2 {
            return self.double();
        }

        let slope = (y2 - y1) / (x2 - x1);
        let xr = slope.square() - x1 - x2;
        let yr = slope * (x1 - xr) - y1;

        G2Affine {
            x: xr,
            y: yr,
            is_infinity: false,
            _marker: PhantomData::<C>,
        }
    }
}

impl<'a, 'b, C: Curve> Sub<&'b G2Affine<C>> for &'a G2Affine<C> {
    type Output = G2Affine<C>;

    #[inline]
    fn sub(self, other: &'b G2Affine<C>) -> G2Affine<C> {
        self + (-other)
    }
}

impl_binops_multiplicative!(G2Affine<C>, Fr<C>);
impl_binops_additive!(G2Affine<C>, G2Affine<C>);

#[cfg(test)]
mod test {
    use crate::common::Bls12381Curve;

    use super::*;

    #[test]
    fn test_doubling() {
        let tmp = G2Affine::<Bls12381Curve>::identity().double();
        assert_eq!(tmp, G2Affine::<Bls12381Curve>::identity());

        let tmp = G2Affine::<Bls12381Curve>::generator().double();
        assert!(!tmp.is_zero());

        assert_eq!(
            G2Affine::from(tmp),
            G2Affine::<Bls12381Curve>::new(
                Fp2::<Bls12381Curve>::new(
                    Fp::from_raw_unchecked([
                        0xe9d9_e2da_9620_f98b,
                        0x54f1_1993_46b9_7f36,
                        0x3db3_b820_376b_ed27,
                        0xcfdb_31c9_b0b6_4f4c,
                        0x41d7_c127_8635_4493,
                        0x0571_0794_c255_c064,
                    ]),
                    Fp::from_raw_unchecked([
                        0xd6c1_d3ca_6ea0_d06e,
                        0xda0c_bd90_5595_489f,
                        0x4f53_52d4_3479_221d,
                        0x8ade_5d73_6f8c_97e0,
                        0x48cc_8433_925e_f70e,
                        0x08d7_ea71_ea91_ef81,
                    ]),
                ),
                Fp2::<Bls12381Curve>::new(
                    Fp::from_raw_unchecked([
                        0x15ba_26eb_4b0d_186f,
                        0x0d08_6d64_b7e9_e01e,
                        0xc8b8_48dd_652f_4c78,
                        0xeecf_46a6_123b_ae4f,
                        0x255e_8dd8_b6dc_812a,
                        0x1641_42af_21dc_f93f,
                    ]),
                    Fp::from_raw_unchecked([
                        0xf9b4_a1a8_9598_4db4,
                        0xd417_b114_cccf_f748,
                        0x6856_301f_c89f_086e,
                        0x41c7_7787_8931_e3da,
                        0x3556_b155_066a_2105,
                        0x00ac_f7d3_25cb_89cf,
                    ]),
                ),
                false,
            )
        );
    }

    #[test]
    fn test_torsion_free() {
        let a = G2Affine::<Bls12381Curve>::new(
            Fp2::new(
                Fp::from_raw_unchecked([
                    0x89f5_50c8_13db_6431,
                    0xa50b_e8c4_56cd_8a1a,
                    0xa45b_3741_14ca_e851,
                    0xbb61_90f5_bf7f_ff63,
                    0x970c_a02c_3ba8_0bc7,
                    0x02b8_5d24_e840_fbac,
                ]),
                Fp::from_raw_unchecked([
                    0x6888_bc53_d707_16dc,
                    0x3dea_6b41_1768_2d70,
                    0xd8f5_f930_500c_a354,
                    0x6b5e_cb65_56f5_c155,
                    0xc96b_ef04_3477_8ab0,
                    0x0508_1505_5150_06ad,
                ]),
            ),
            Fp2::new(
                Fp::from_raw_unchecked([
                    0x3cf1_ea0d_434b_0f40,
                    0x1a0d_c610_e603_e333,
                    0x7f89_9561_60c7_2fa0,
                    0x25ee_03de_cf64_31c5,
                    0xeee8_e206_ec0f_e137,
                    0x0975_92b2_26df_ef28,
                ]),
                Fp::from_raw_unchecked([
                    0x71e8_bb5f_2924_7367,
                    0xa5fe_049e_2118_31ce,
                    0x0ce6_b354_502a_3896,
                    0x93b0_1200_0997_314e,
                    0x6759_f3b6_aa5b_42ac,
                    0x1569_44c4_dfe9_2bbb,
                ]),
            ),
            false,
        );
        // assert!(!bool::from(a.is_torsion_free()));
        println!("Generator: {:?}", G2Affine::<Bls12381Curve>::generator());
        assert!(G2Affine::<Bls12381Curve>::generator().is_torsion_free());
    }
}
