use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::common::AffinePoint;
use crate::fp::{Bls12381, Bn254, FpElement};
use crate::{fp2::Fp2, fr::Fr};

pub(crate) trait G2AffineType: FpElement {
    type G2AffineType;
    fn is_on_curve(p: &Self::G2AffineType) -> bool {
        let x = p.x;
        let y = p.y;

        // y^2 = x^3 + B
        y.square()
            == x.square() * x
                + Fp2::new(
                    Self::from_raw_unchecked(Self::B2_X),
                    Self::from_raw_unchecked(Self::B2_Y),
                )
    }

    fn is_valid(p: &Self::G2AffineType) -> Result<(), String>;
}

impl G2AffineType for Bls12381 {
    type G2AffineType = G2Affine<Bls12381>;

    fn is_valid(p: &Self::G2AffineType) -> Result<(), String> {
        if !Self::is_on_curve(p) {
            return Err("Point is not on curve".to_string());
        }

        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2::<Self>::new(
            Self::zero(),
            Self::from_raw_unchecked([
                0x8bfd00000000aaad,
                0x409427eb4f49fffd,
                0x897d29650fb85f9b,
                0xaa0d857d89759ad4,
                0xec02408663d4de85,
                0x1a0111ea397fe699,
            ]),
        );
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fp2::<Self>::new(
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
        );

        let lhs = G2Affine::<Self>::new(
            p.x.frobenius_map() * psi_coeff_x,
            p.y.frobenius_map() * psi_coeff_y,
            false,
        );

        let rhs = -p.mul_by_x();

        (lhs == rhs)
            .then(|| ())
            .ok_or("Point is not torsion free".to_string())
    }
}

impl G2AffineType for Bn254 {
    type G2AffineType = G2Affine<Bn254>;

    fn is_valid(p: &Self::G2AffineType) -> Result<(), String> {
        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2::<Self>::from_raw_unchecked([
            Self::zero(),
            Self::from_raw_unchecked([
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
            ]),
        ]);
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fp2::<Self>::from_raw_unchecked([
            Self::from_raw_unchecked([
                0x3e2f585da55c9ad1,
                0x4294213d86c18183,
                0x382844c88b623732,
                0x92ad2afd19103e18,
            ]),
            Self::from_raw_unchecked([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
            ]),
        ]);
        let lhs = p.mul_by_x();
        let rhs = G2Affine::<Self>::new(
            p.x.frobenius_map() * psi_coeff_x,
            p.y.frobenius_map() * psi_coeff_y,
            false,
        );

        (lhs == rhs)
            .then(|| ())
            .ok_or("Point is not torsion free".to_string())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct G2Affine<F: G2AffineType> {
    pub(crate) x: Fp2<F>,
    pub(crate) y: Fp2<F>,
    is_infinity: bool,
}

impl<F: G2AffineType> PartialEq for G2Affine<F> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<F: G2AffineType> Eq for G2Affine<F> {}

impl<F: G2AffineType> G2Affine<F> {
    fn new(x: &Self, y: &Self, is_infinity: bool) -> Self {
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
                F::from_raw_unchecked(F::G2_X0),
                F::from_raw_unchecked(F::G2_X1),
            ),
            y: Fp2::new(
                F::from_raw_unchecked(F::G2_Y0),
                F::from_raw_unchecked(F::G2_Y1),
            ),
            is_infinity: false,
        }
    }

    fn is_valid(&self) -> Result<(), String> {
        self.is_identity()
            .then(|| ())
            .ok_or_else(|| F::is_valid(&self).unwrap_err())
    }

    fn random(mut rng: impl rand::Rng) -> Self {
        let b = Fp2 {
            c0: F::from_raw_unchecked(F::B2_X),
            c1: F::from_raw_unchecked(F::B2_Y),
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

        let three = F::<F>::from(3);
        let two = F::<F>::from(2);

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

impl<F: G2AffineType> G2Affine<F> {
    fn mul_by_x(&self) -> Self {
        let mut xself = G2Affine::identity();
        let mut x = F::X >> 1;
        let mut tmp = *self;
        while x != 0 {
            tmp = tmp.double();

            if x % 2 == 1 {
                xself += tmp;
            }
            x >>= 1;
        }
        xself
    }

    pub fn from_compressed_unchecked(bytes: &[u8; 96]) -> Option<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = (bytes[0] >> 7) & 1 == 1;
        let infinity_flag_set = (bytes[0] >> 6) & 1 == 1;
        let sort_flag_set = (bytes[0] >> 5) & 1 == 1;

        // Attempt to obtain the x-coordinate
        let xc1 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[0] &= 0b0001_1111;

            Some(F::from_bytes(&tmp))
        };
        let xc0 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[48..96]);

            Some(F::from_bytes(&tmp))
        };

        xc1.and_then(|xc1| {
            xc0.and_then(|xc0| {
                let x = Fp2 { c0: xc0, c1: xc1 };

                if infinity_flag_set && compression_flag_set && !sort_flag_set && x.is_zero() {
                    // Infinity flag is set and x-coordinate is zero
                    Some(G2Affine::identity())
                } else if !infinity_flag_set && compression_flag_set {
                    // Recover a y-coordinate given x by y = sqrt(x^3 + 4)
                    let y_result = ((x.square() * x)
                        + Fp2::new(
                            F::from_raw_unchecked(F::B2_X),
                            F::from_raw_unchecked(F::B2_Y),
                        ))
                    .sqrt();

                    y_result.map(|y| {
                        // Switch to the correct y-coordinate if necessary
                        let y = if y.lexicographically_largest() ^ sort_flag_set {
                            -y
                        } else {
                            y
                        };

                        G2Affine {
                            x,
                            y,
                            is_infinity: infinity_flag_set,
                        }
                    })
                } else {
                    None
                }
            })
        })
    }
}

impl<F: G2AffineType> Neg for G2Affine<F> {
    type Output = G2Affine<F>;

    fn neg(self) -> G2Affine<F> {
        G2Affine {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

impl<'a, 'b, F: G2AffineType> Mul<&'b Fr<F>> for &'a G2Affine<F> {
    type Output = G2Affine<F>;

    #[inline]
    fn mul(self, other: &'b Fr<F>) -> G2Affine<F> {
        let mut acc = G2Affine::<F>::identity();

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

impl<'a, 'b, F: G2AffineType> Add<&'b G2Affine<F>> for &'a G2Affine<F> {
    type Output = G2Affine<F>;

    #[inline]
    fn add(self, other: &'b G2Affine<F>) -> G2Affine<F> {
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

impl<'a, 'b, F: G2AffineType> Sub<&'b G2Affine<F>> for &'a G2Affine<F> {
    type Output = G2Affine<F>;

    #[inline]
    fn sub(self, other: &'b G2Affine<F>) -> G2Affine<F> {
        if self == other {
            return G2Affine::<F>::identity();
        }
        self + (-other)
    }
}

pub struct G2Projective<F: G2AffineType> {
    pub(crate) x: Fp2<F>,
    pub(crate) y: Fp2<F>,
    pub(crate) z: Fp2<F>,
}

impl<F: G2AffineType> G2Projective<F> {
    pub fn to_affine(&self) -> G2Affine<F> {
        if self.is_identity() {
            return G2Affine::<F>::identity();
        }

        let zinv = self.z.invert().unwrap();
        let zinv2 = zinv.square();
        let zinv3 = zinv2 * zinv;

        let x = self.x * zinv2;
        let y = self.y * zinv3;

        G2Affine {
            x,
            y,
            is_infinity: false,
        }
    }

    pub fn from_affine(p: G2Affine<F>) -> Self {
        G2Projective {
            x: p.x,
            y: p.y,
            z: Fp2::one(),
        }
    }

    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    macro_rules! g2_tests {
        ($mod_name:ident, $curve:ident, $g2_affine:ident, $fr:ident, $fp2:ident) => {
            mod $mod_name {
                use super::*;

                #[test]
                fn test_scalar_multiplication() {
                    let mut rng = rand::thread_rng();
                    for _ in 0..10 {
                        let r: u64 = rng.gen::<u64>() % 100;
                        let k = $fr::<$curve>::from(r);
                        let a = $g2_affine::<$curve>::random(&mut rng);
                        let lhs = &a * &k;
                        let rhs = (0..r).fold($g2_affine::<$curve>::identity(), |acc, _| acc + &a);
                        assert_eq!(lhs, rhs);
                    }
                }

                #[test]
                fn test_affine_addition() {
                    {
                        let a = $g2_affine::<$curve>::identity();
                        let b = $g2_affine::<$curve>::identity();
                        let c = &a + &b;
                        assert!(c.is_identity());
                        assert!(c.is_valid().is_ok());
                    }
                    {
                        let a = $g2_affine::<$curve>::identity();
                        let b = $g2_affine::<$curve>::new(
                            $fp2::new(
                                $curve::from_raw_unchecked([
                                    0xd48056c8c121bdb8,
                                    0x0bac0326a805bbef,
                                    0xb4510b647ae3d177,
                                    0xc6e47ad4fa403b02,
                                    0x260805272dc51051,
                                    0x024aa2b2f08f0a91,
                                ]),
                                $curve::from_raw_unchecked([
                                    0xe5ac7d055d042b7e,
                                    0x334cf11213945d57,
                                    0xb5da61bbdc7f5049,
                                    0x596bd0d09920b61a,
                                    0x7dacd3a088274f65,
                                    0x13e02b6052719f60,
                                ]),
                            ),
                            $fp2::new(
                                $curve::from_raw_unchecked([
                                    0xe193548608b82801,
                                    0x923ac9cc3baca289,
                                    0x6d429a695160d12c,
                                    0xadfd9baa8cbdd3a7,
                                    0x8cc9cdc6da2e351a,
                                    0x0ce5d527727d6e11,
                                ]),
                                $curve::from_raw_unchecked([
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
                        assert!(c == $g2_affine::<$curve>::generator());
                    }
                    {
                        let a = $g2_affine::<$curve>::generator().double().double();
                        let b = $g2_affine::<$curve>::generator().double();
                        let c = a + b;

                        let mut d = $g2_affine::<$curve>::generator();
                        for _ in 0..5 {
                            d += $g2_affine::<$curve>::generator();
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
                    let tmp = $g2_affine::<$curve>::identity().double();
                    assert_eq!(tmp, $g2_affine::<$curve>::identity());

                    let tmp = $g2_affine::<$curve>::generator().double();
                    assert!(!tmp.is_zero());

                    assert_eq!(
                        tmp,
                        $g2_affine::<$curve>::new(
                            $fp2::<$curve>::new(
                                $curve::from_raw_unchecked([
                                    0xc952aacab827a053,
                                    0x81f14b0bf3611b78,
                                    0xe1ea1e1e4d00dbae,
                                    0x3bc0b995b8825e0e,
                                    0xd2370f17cc7ed586,
                                    0x1638533957d540a9,
                                ]),
                                $curve::from_raw_unchecked([
                                    0x6178288c47c33577,
                                    0xc6c886f6b57ec72a,
                                    0x728114d1031e1572,
                                    0xd70662a904ba1074,
                                    0x9f520e47730a124f,
                                    0x0a4edef9c1ed7f72,
                                ]),
                            ),
                            $fp2::<$curve>::new(
                                $curve::from_raw_unchecked([
                                    0x999d95d71e4c9899,
                                    0xe88dece9764bf3bd,
                                    0xbfe6bd221e47aa8a,
                                    0x9a66da69bf91009c,
                                    0x0aeb8dca2b525678,
                                    0x0468fb440d82b063,
                                ]),
                                $curve::from_raw_unchecked([
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
                    let a = $g2_affine::<$curve>::new(
                        $fp2::new(
                            $curve::from_raw_unchecked([
                                0x89f5_50c8_13db_6431,
                                0xa50b_e8c4_56cd_8a1a,
                                0xa45b_3741_14ca_e851,
                                0xbb61_90f5_bf7f_ff63,
                                0x970c_a02c_3ba8_0bc7,
                                0x02b8_5d24_e840_fbac,
                            ]),
                            $curve::from_raw_unchecked([
                                0x6888_bc53_d707_16dc,
                                0x3dea_6b41_1768_2d70,
                                0xd8f5_f930_500c_a354,
                                0x6b5e_cb65_56f5_c155,
                                0xc96b_ef04_3477_8ab0,
                                0x0508_1505_5150_06ad,
                            ]),
                        ),
                        $fp2::new(
                            $curve::from_raw_unchecked([
                                0x3cf1_ea0d_434b_0f40,
                                0x1a0d_c610_e603_e333,
                                0x7f89_9561_60c7_2fa0,
                                0x25ee_03de_cf64_31c5,
                                0xeee8_e206_ec0f_e137,
                                0x0975_92b2_26df_ef28,
                            ]),
                            $curve::from_raw_unchecked([
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
                    assert!($g2_affine::<$curve>::generator().is_torsion_free());
                }

                #[test]
                fn test_double_and_add_arithmetic() {
                    let p = $g2_affine::<$curve>::random(&mut rand::thread_rng());
                    let q = $g2_affine::<$curve>::random(&mut rand::thread_rng());
                    let r = $g2_affine::<$curve>::random(&mut rand::thread_rng());

                    let double_p_add_q = p.double() + q;
                    let p_plus_q_plus_p = (p + q) + p;

                    assert_eq!(p + q, q + p);
                    assert_eq!((p + q) + r, p + (q + r));
                    assert_eq!(double_p_add_q, p_plus_q_plus_p);
                }
            }
        };
    }

    g2_tests!(bls12381_g2_tests, Bls12381, G2Affine, Fr, Fp2);
    g2_tests!(bn254_g2_tests, Bn254, G2Affine, Fr, Fp2);
}
