use std::ops::{Add, Mul, Neg, Sub};

use crate::common::AffinePoint;
use crate::fp::{Bls12381, Bn254, FpElement};
use crate::fp2::Fp2Element;
use crate::fr::FrElement;
use crate::{fp2::Fp2, fr::Fr};

pub trait G2Element: Fp2Element + AffinePoint {
    const B: Fp2<Self>;
    fn is_on_curve(x: &Fp2<Self>, y: &Fp2<Self>) -> bool {
        // y^2 = x^3 + B
        y.square() == x.square() * *x + Self::B
    }

    fn is_valid(p: &G2Affine<Self>) -> Result<(), String>;
    fn from_compressed_unchecked(bytes: &[u8]) -> Option<G2Affine<Self>>;
    fn generator() -> G2Affine<Self>;
}

impl G2Element for Bls12381 {
    const B: Fp2<Bls12381> = Fp2::<Bls12381>::new(
        Bls12381::from_raw_unchecked([4, 0, 0, 0, 0, 0]),
        Bls12381::from_raw_unchecked([4, 0, 0, 0, 0, 0]),
    );

    fn is_valid(p: &G2Affine<Self>) -> Result<(), String> {
        if !Self::is_on_curve(&p.x, &p.y) {
            return Err("Point is not on curve".to_string());
        }

        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2::<Self>::new(
            <Self as FpElement>::zero(),
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

    fn generator() -> G2Affine<Self> {
        const G2_X0: Bls12381 = Bls12381::from_raw_unchecked([
            0xd48056c8c121bdb8,
            0x0bac0326a805bbef,
            0xb4510b647ae3d177,
            0xc6e47ad4fa403b02,
            0x260805272dc51051,
            0x024aa2b2f08f0a91,
        ]);

        const G2_X1: Bls12381 = Bls12381::from_raw_unchecked([
            0xe5ac7d055d042b7e,
            0x334cf11213945d57,
            0xb5da61bbdc7f5049,
            0x596bd0d09920b61a,
            0x7dacd3a088274f65,
            0x13e02b6052719f60,
        ]);

        const G2_Y0: Bls12381 = Bls12381::from_raw_unchecked([
            0xe193548608b82801,
            0x923ac9cc3baca289,
            0x6d429a695160d12c,
            0xadfd9baa8cbdd3a7,
            0x8cc9cdc6da2e351a,
            0x0ce5d527727d6e11,
        ]);

        const G2_Y1: Bls12381 = Bls12381::from_raw_unchecked([
            0xaaa9075ff05f79be,
            0x3f370d275cec1da1,
            0x267492ab572e99ab,
            0xcb3e287e85a763af,
            0x32acd2b02bc28b99,
            0x0606c4a02ea734cc,
        ]);
        G2Affine::<Self>::new(Fp2::new(G2_X0, G2_X1), Fp2::new(G2_Y0, G2_Y1), false)
    }

    fn from_compressed_unchecked(bytes: &[u8]) -> Option<G2Affine<Bls12381>> {
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

            Some(Bls12381::from_bytes(&tmp))
        };
        let xc0 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[48..96]);

            Some(Bls12381::from_bytes(&tmp))
        };

        xc1.and_then(|xc1| {
            xc0.and_then(|xc0| {
                let x = Fp2 { c0: xc0, c1: xc1 };

                if infinity_flag_set && compression_flag_set && !sort_flag_set && x.is_zero() {
                    // Infinity flag is set and x-coordinate is zero
                    Some(G2Affine::identity())
                } else if !infinity_flag_set && compression_flag_set {
                    // Recover a y-coordinate given x by y = sqrt(x^3 + 4)
                    let y_result = ((x.square() * x) + Self::B).sqrt();

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

impl G2Element for Bn254 {
    // The B coefficient for the BN254 curve is twist(0, 3)
    const B: Fp2<Bn254> = Fp2::<Bn254>::new(
        Bn254::from_raw_unchecked([
            0xe4a2bd0685c315d2,
            0xa74fa084e52d1852,
            0xcd2cafadeed8fdf4,
            0x9713b03af0fed4,
        ]),
        Bn254::from_raw_unchecked([
            0x3267e6dc24a138e5,
            0xb5b4c5e559dbefa3,
            0x81be18991be06ac3,
            0x2b149d40ceb8aaae,
        ]),
    );

    fn is_valid(p: &G2Affine<Self>) -> Result<(), String> {
        if !Self::is_on_curve(&p.x, &p.y) {
            return Err("Point is not on curve".to_string());
        }
        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2::<Self>::new(
            Self::zero(),
            Self::from_raw_unchecked([
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
            ]),
        );
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fp2::<Self>::new(
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
        );

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

    fn generator() -> G2Affine<Self> {
        const G2_X0: Bn254 = Bn254::from_raw_unchecked([
            0x46debd5cd992f6ed,
            0x674322d4f75edadd,
            0x426a00665e5c4479,
            0x1800deef121f1e76,
        ]);

        const G2_X1: Bn254 = Bn254::from_raw_unchecked([
            0x97e485b7aef312c2,
            0xf1aa493335a9e712,
            0x7260bfb731fb5d25,
            0x198e9393920d483a,
        ]);

        const G2_Y0: Bn254 = Bn254::from_raw_unchecked([
            0x4ce6cc0166fa7daa,
            0xe3d1e7690c43d37b,
            0x4aab71808dcb408f,
            0x12c85ea5db8c6deb,
        ]);

        const G2_Y1: Bn254 = Bn254::from_raw_unchecked([
            0x55acdadcd122975b,
            0xbc4b313370b38ef3,
            0xec9e99ad690c3395,
            0x90689d0585ff075,
        ]);

        G2Affine::<Self>::new(Fp2::new(G2_X0, G2_X1), Fp2::new(G2_Y0, G2_Y1), false)
    }

    fn from_compressed_unchecked(bytes: &[u8]) -> Option<G2Affine<Self>> {
        todo!()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct G2Affine<F: G2Element> {
    pub(crate) x: Fp2<F>,
    pub(crate) y: Fp2<F>,
    is_infinity: bool,
}

impl<F: G2Element> PartialEq for G2Affine<F> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<F: G2Element> Eq for G2Affine<F> {}

impl<F: G2Element> G2Affine<F> {
    const fn new(x: Fp2<F>, y: Fp2<F>, is_infinity: bool) -> Self {
        G2Affine { x, y, is_infinity }
    }

    fn identity() -> Self {
        G2Affine {
            x: Fp2::zero(),
            y: Fp2::one(),
            is_infinity: true,
        }
    }

    pub(crate) fn is_identity(&self) -> bool {
        self.is_infinity
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }

    pub fn generator() -> Self {
        F::generator()
    }

    pub fn is_valid(&self) -> Result<(), String> {
        self.is_identity()
            .then(|| ())
            .ok_or_else(|| F::is_valid(&self).unwrap_err())
    }

    fn random(mut rng: impl rand::Rng) -> Self {
        let b = F::B;
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

    pub fn double(&self) -> Self {
        if self.is_infinity {
            return Self::identity();
        }
        let x = self.x;
        let y = self.y;

        // Check if y is zero
        if y.is_zero() {
            return Self::identity();
        }

        let three = F::from(3);
        let two = F::from(2);

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

impl<F: G2Element> G2Affine<F> {
    fn mul_by_x(&self) -> Self {
        let mut xself = G2Affine::identity();
        let mut x = F::X >> 1;
        let mut tmp = *self;
        while x != 0 {
            tmp = tmp.double();

            if x % 2 == 1 {
                xself = xself + tmp;
            }
            x >>= 1;
        }
        xself
    }
}

impl<F: G2Element> Neg for G2Affine<F> {
    type Output = G2Affine<F>;

    fn neg(self) -> G2Affine<F> {
        G2Affine {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

impl<G: G2Element, F: FrElement> Mul<Fr<F>> for G2Affine<G> {
    type Output = G2Affine<G>;

    #[inline]
    fn mul(self, other: Fr<F>) -> G2Affine<G> {
        let mut acc = G2Affine::<G>::identity();

        for bit in other
            .0
            .iter()
            .rev()
            .flat_map(|&x| (0..64).rev().map(move |i| (x >> i) & 1 == 1))
            .skip(1)
        {
            acc = acc.double();

            if bit {
                acc = acc + self;
            }
        }

        acc
    }
}

impl<F: G2Element> Add<G2Affine<F>> for G2Affine<F> {
    type Output = G2Affine<F>;

    #[inline]
    fn add(self, other: G2Affine<F>) -> G2Affine<F> {
        if self.is_infinity {
            return other;
        }

        if other.is_infinity {
            return self;
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

impl<F: G2Element> Sub<G2Affine<F>> for G2Affine<F> {
    type Output = G2Affine<F>;

    #[inline]
    fn sub(self, other: G2Affine<F>) -> G2Affine<F> {
        if self == other {
            return G2Affine::<F>::identity();
        }
        self + (-other)
    }
}

pub struct G2Projective<F: G2Element> {
    pub(crate) x: Fp2<F>,
    pub(crate) y: Fp2<F>,
    pub(crate) z: Fp2<F>,
}

impl<F: G2Element> G2Projective<F> {
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
