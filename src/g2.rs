use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::common::{AffinePoint, Curve};
use crate::fp::Fp;
use crate::{fp2::Fp2, fr::Fr};

#[derive(Clone, Copy, Debug)]
pub struct G2Affine<C: Curve> {
    x: Fp2<C>,
    y: Fp2<C>,
    is_infinity: bool,
}

impl<C: Curve> PartialEq for G2Affine<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: Curve> AffinePoint<C> for G2Affine<C> {
    type Dtype = Fp2<C>;

    fn new(x: Self::Dtype, y: Self::Dtype, is_infinity: bool) -> Self {
        G2Affine { x, y, is_infinity }
    }

    fn identity() -> Self {
        G2Affine {
            x: Fp2::zero(),
            y: Fp2::one(),
            is_infinity: true,
        }
    }

    fn is_identity(&self) -> bool {
        self.is_infinity
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
        }
    }

    fn is_valid(&self) -> Result<(), String> {
        if self.is_infinity {
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
        let b = Fp2 {
            c0: Fp::from_raw_unchecked(C::B2_X),
            c1: Fp::from_raw_unchecked(C::B2_Y),
        };
        loop {
            let x = Fp2::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            // Obtain the corresponding y-coordinate given x as y = sqrt(x^3 + 4)
            let p = ((x.square() * x) + b).sqrt().map(|y| G2Affine {
                x,
                y: if flip_sign { -y } else { y },
                is_infinity: false,
            });

            if p.is_some().into() {
                // let p = p.unwrap().to_curve().clear_cofactor();
                let p = p.unwrap();

                if bool::from(!p.is_identity()) {
                    return p;
                }
            }
        }
    }

    fn double(&self) -> Self {
        if self.is_infinity {
            return Self::identity();
        }
        let x = self.x;
        let y = self.y;

        // Check if y is zero
        if y.is_zero() {
            return Self::identity();
        }

        let three = Fp::<C>::from(3);
        let two = Fp::<C>::from(2);

        let slope = (x.square() * three) / (y * two);
        let x_new = slope.square() - x * two;
        let y_new = slope * (x - x_new) - y;

        Self {
            x: x_new,
            y: y_new,
            is_infinity: false,
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
        }
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fr<C>> for &'a G2Affine<C> {
    type Output = G2Affine<C>;

    #[inline]
    fn mul(self, other: &'b Fr<C>) -> G2Affine<C> {
        let mut acc = G2Affine::<C>::identity();

        for bit in other
            .0
            .iter()
            .rev()
            .flat_map(|&x| (0..64).rev().map(move |i| (x >> i) & 1 == 1))
            .skip(1)
        {
            acc = acc.double();

            if bit {
                acc += self;
            }
        }

        acc
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
    use rand::Rng;

    use super::*;

    #[test]
    fn test_scalar_multiplication() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let r: u64 = rng.gen::<u64>() % 100;
            let k = Fr::<Bls12381Curve>::from(r);
            let a = G2Affine::<Bls12381Curve>::random(&mut rng);
            let lhs = &a * &k;
            let rhs = (0..r).fold(G2Affine::<Bls12381Curve>::identity(), |acc, _| acc + &a);
            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn test_affine_addition() {
        {
            let a = G2Affine::<Bls12381Curve>::identity();
            let b = G2Affine::<Bls12381Curve>::identity();
            let c = &a + &b;
            assert!(c.is_identity());
            assert!(c.is_valid().is_ok());
        }
        {
            let a = G2Affine::<Bls12381Curve>::identity();
            let b = G2Affine::<Bls12381Curve>::new(
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xd48056c8c121bdb8,
                        0x0bac0326a805bbef,
                        0xb4510b647ae3d177,
                        0xc6e47ad4fa403b02,
                        0x260805272dc51051,
                        0x024aa2b2f08f0a91,
                    ]),
                    Fp::from_raw_unchecked([
                        0xe5ac7d055d042b7e,
                        0x334cf11213945d57,
                        0xb5da61bbdc7f5049,
                        0x596bd0d09920b61a,
                        0x7dacd3a088274f65,
                        0x13e02b6052719f60,
                    ]),
                ),
                Fp2::new(
                    Fp::from_raw_unchecked([
                        0xe193548608b82801,
                        0x923ac9cc3baca289,
                        0x6d429a695160d12c,
                        0xadfd9baa8cbdd3a7,
                        0x8cc9cdc6da2e351a,
                        0x0ce5d527727d6e11,
                    ]),
                    Fp::from_raw_unchecked([
                        0xaaa9075ff05f79be,
                        0x3f370d275cec1da1,
                        0x267492ab572e99ab,
                        0xcb3e287e85a763af,
                        0x32acd2b02bc28b99,
                        0x0606c4a02ea734cc,
                    ]),
                ),
                false,
            );
            let c = a + b;
            assert!(!c.is_identity());
            assert!(c.is_on_curve());
            assert!(c == G2Affine::<Bls12381Curve>::generator());
        }
        {
            let a = G2Affine::<Bls12381Curve>::generator().double().double();
            let b = G2Affine::<Bls12381Curve>::generator().double();
            let c = a + b;

            let mut d = G2Affine::<Bls12381Curve>::generator();
            for _ in 0..5 {
                d += G2Affine::<Bls12381Curve>::generator();
            }
            assert!(!c.is_identity());
            assert!(c.is_on_curve());
            assert!(!d.is_identity());
            assert!(d.is_on_curve());
            assert_eq!(c, d);
        }
    }

    #[test]
    fn test_doubling() {
        let tmp = G2Affine::<Bls12381Curve>::identity().double();
        assert_eq!(tmp, G2Affine::<Bls12381Curve>::identity());

        let tmp = G2Affine::<Bls12381Curve>::generator().double();
        assert!(!tmp.is_zero());

        assert_eq!(
            tmp,
            G2Affine::<Bls12381Curve>::new(
                Fp2::<Bls12381Curve>::new(
                    Fp::from_raw_unchecked([
                        0xc952aacab827a053,
                        0x81f14b0bf3611b78,
                        0xe1ea1e1e4d00dbae,
                        0x3bc0b995b8825e0e,
                        0xd2370f17cc7ed586,
                        0x1638533957d540a9,
                    ]),
                    Fp::from_raw_unchecked([
                        0x6178288c47c33577,
                        0xc6c886f6b57ec72a,
                        0x728114d1031e1572,
                        0xd70662a904ba1074,
                        0x9f520e47730a124f,
                        0x0a4edef9c1ed7f72,
                    ]),
                ),
                Fp2::<Bls12381Curve>::new(
                    Fp::from_raw_unchecked([
                        0x999d95d71e4c9899,
                        0xe88dece9764bf3bd,
                        0xbfe6bd221e47aa8a,
                        0x9a66da69bf91009c,
                        0x0aeb8dca2b525678,
                        0x0468fb440d82b063,
                    ]),
                    Fp::from_raw_unchecked([
                        0xacdefd8b6e36ccf3,
                        0x422e1aa0a59c8967,
                        0x97003f7a13c308f5,
                        0xa43253d9c66c4116,
                        0x38b361543f887136,
                        0x0f6d4552fa65dd26,
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
        assert!(!a.is_torsion_free());
        assert!(G2Affine::<Bls12381Curve>::generator().is_torsion_free());
    }
}
