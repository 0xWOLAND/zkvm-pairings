use std::ops::{Add, Mul, Neg, Sub};

use crate::common::AffinePoint;
use crate::fp::{Bls12381, Bn254, FpElement};
use crate::fr::Fr;

pub(crate) trait G1Element: FpElement {
    type G1AffineType;
    fn is_on_curve(x: &Self, y: &Self) -> bool {
        y.square() == x.square() * x + Self::from_raw_unchecked(Self::B)
    }
    fn is_valid(p: &Self::G1AffineType) -> Result<(), String>;
}

impl G1Element for Bls12381 {
    type G1AffineType = G1Affine<Bls12381>;
    fn is_valid(p: &G1Affine<Bls12381>) -> Result<(), String> {
        let x = p.x;
        let y = p.y;

        // Check if the point is on the curve
        if !Self::is_on_curve(&x, &y) {
            return Err("Point is not on curve".to_string());
        }

        // Check if the point is torsion free
        let generator = Self::from_raw_unchecked(Self::X);
        let beta = Self::from_raw_unchecked([
            0x2e01fffffffefffe,
            0xde17d813620a0002,
            0xddb3a93be6f89688,
            0xba69c6076a0f77ea,
            0x5f19672fdf76ce51,
            0x0,
        ]);
        let rhs = p.mul_by_x().mul_by_x().neg();
        let lhs = G1Affine::new(p.x * beta, p.y, false);

        (lhs == rhs)
            .then(|| ())
            .ok_or("Point is not torsion free".to_string())
    }
}

impl G1Element for Bn254 {
    type G1AffineType = G1Affine<Bn254>;

    // G1Affine Bn254 elements are always torsion free
    fn is_valid(p: &Self::G1AffineType) -> Result<(), String> {
        // Check if the point is on the curve
        let x = p.x;
        let y = p.y;

        Self::is_on_curve(&x, &y)
            .then(|| ())
            .ok_or("Point is not on curve".to_string())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct G1Affine<F: G1Element> {
    pub(crate) x: F,
    pub(crate) y: F,
    is_infinity: bool,
}

impl<F: G1Element> PartialEq for G1Affine<F> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<F: G1Element> Eq for G1Affine<F> {}

impl<F: G1Element> G1Affine<F> {
    fn new(x: F, y: F, is_infinity: bool) -> Self {
        G1Affine { x, y, is_infinity }
    }

    fn identity() -> Self {
        G1Affine {
            x: F::zero(),
            y: F::one(),
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
        G1Affine {
            x: F::from_raw_unchecked(F::G1_X),
            y: F::from_raw_unchecked(F::G1_Y),
            is_infinity: false,
        }
    }

    fn is_valid(&self) -> Result<(), String> {
        self.is_identity()
            .then(|| ())
            .ok_or_else(|| F::is_valid(&self).unwrap_err())
    }

    fn random(mut rng: impl rand::Rng) -> Self {
        let b = F::from_raw_unchecked(F::B);
        loop {
            let x = F::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            // Obtain the corresponding y-coordinate given x as y = sqrt(x^3 + 4)
            let p = ((x.square() * x) + b).sqrt().map(|y| G1Affine {
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
        let x = self.x;
        let y = self.y;

        if self.is_identity() {
            return Self::identity();
        }

        let slope = (F::from(3) * x.square()) / (F::from(2) * y);
        let xr = slope.square() - F::from(2) * x;
        let yr = slope * (x - xr) - y;

        Self {
            x: xr,
            y: yr,
            is_infinity: false,
        }
    }
}

impl<F: G1Element> G1Affine<F> {
    fn mul_by_x(&self) -> Self {
        let mut xself = G1Affine::identity();
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

    pub fn from_compressed_unchecked(bytes: &[u8; 48]) -> Option<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = (bytes[0] >> 7) & 1 == 1;
        let infinity_flag_set = (bytes[0] >> 6) & 1 == 1;
        let sort_flag_set = (bytes[0] >> 5) & 1 == 1;

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[0] &= 0b0001_1111;

            F::from_bytes(&tmp)
        };

        if infinity_flag_set && compression_flag_set && !sort_flag_set && x.is_zero() {
            // Infinity flag is set and x-coordinate is zero
            Some(G1Affine::identity())
        } else if !infinity_flag_set && compression_flag_set {
            // Recover a y-coordinate given x by y = sqrt(x^3 + B)
            let y_result = ((x.square() * x) + F::<F>::from_raw_unchecked(F::B)).sqrt();

            y_result.map(|y| {
                let y = match !(y.is_lexicographically_largest() ^ sort_flag_set) {
                    true => y,
                    false => -y,
                };

                G1Affine {
                    x,
                    y,
                    is_infinity: infinity_flag_set,
                }
            })
        } else {
            None
        }
    }
}

impl<F: G1Element> Neg for G1Affine<F> {
    type Output = G1Affine<F>;

    fn neg(self) -> G1Affine<F> {
        G1Affine {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

impl<'a, 'b, F: G1Element> Mul<&'b Fr<F>> for &'a G1Affine<F> {
    type Output = G1Affine<F>;

    #[inline]
    fn mul(self, other: &'b Fr<F>) -> G1Affine<F> {
        let mut acc = G1Affine::<F>::identity();

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

impl<'a, 'b, F: G1Element> Add<&'b G1Affine<F>> for &'a G1Affine<F> {
    type Output = G1Affine<F>;

    #[inline]
    fn add(self, other: &'b G1Affine<F>) -> G1Affine<F> {
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

        G1Affine {
            x: xr,
            y: yr,
            is_infinity: false,
        }
    }
}

impl<'a, 'b, F: G1Element> Sub<&'b G1Affine<F>> for &'a G1Affine<F> {
    type Output = G1Affine<F>;

    #[inline]
    fn sub(self, other: &'b G1Affine<F>) -> G1Affine<F> {
        if self == other {
            return G1Affine::<F>::identity();
        }
        self + -(other)
    }
}

#[cfg(test)]
mod bls12381_g1_affine_test {
    use crate::fp::Bls12381;

    use super::*;
    use rand::Rng;

    fn bls12381_g1_affine_rand() -> G1Affine<Bls12381> {
        let mut rng = rand::thread_rng();
        G1Affine::random(&mut rng)
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = G1Affine::<Bls12381>::new(
                Bls12381::from_raw_unchecked(x.clone().try_into().unwrap()),
                Bls12381::from_raw_unchecked(y.clone().try_into().unwrap()),
                false,
            );

            let b = G1Affine::<Bls12381>::new(
                Bls12381::from_raw_unchecked(x.clone().try_into().unwrap()),
                Bls12381::from_raw_unchecked(y.clone().try_into().unwrap()),
                false,
            );

            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_valid_point() {
        let a = G1Affine::<Bls12381>::new(
            Bls12381::from_raw_unchecked([
                0x1b72cc2215a57793,
                0x6263ee31e953a86d,
                0x866816fd0826ce7b,
                0x991224ec2a74beb2,
                0x041314ca93e1fee0,
                0x19cdf3807146e68e,
            ]),
            Bls12381::from_raw_unchecked([
                0xe7da608505d48616,
                0x22c9f5c77db40989,
                0x35c0752ac9742b45,
                0x41bfaf99f604d1f8,
                0xf45c6e4fc2780554,
                0x7481b1f261aabac,
            ]),
            false,
        );
        assert!(G1Affine::<Bls12381>::generator().is_valid().unwrap() == ());
    }

    #[test]
    fn test_double() {
        let a = G1Affine::<Bls12381>::new(
            Bls12381::from_raw_unchecked([
                0xfb3af00adb22c6bb,
                0x6c55e83ff97a1aef,
                0xa14e3a3f171bac58,
                0xc3688c4f9774b905,
                0x2695638c4fa9ac0f,
                0x17f1d3a73197d794,
            ]),
            Bls12381::from_raw_unchecked([
                0x0caa232946c5e7e1,
                0xd03cc744a2888ae4,
                0x00db18cb2c04b3ed,
                0xfcf5e095d5d00af6,
                0xa09e30ed741d8ae4,
                0x08b3f481e3aaa0f1,
            ]),
            false,
        );
        let a_double = G1Affine::<Bls12381>::new(
            Bls12381::from_raw_unchecked([
                0xc39a8c5529bf0f4e,
                0xe28f75bb8f1c7c42,
                0x43902d0ac358a62a,
                0x9721db3091280125,
                0x8808c8eb50a9450c,
                0x572cbea904d6746,
            ]),
            Bls12381::from_raw_unchecked([
                0xba86881979749d28,
                0x4c56d9d4cd16bd1b,
                0xf73bb9021d5fd76a,
                0x22ba3ecb8670e461,
                0x22fda673779d8e38,
                0x166a9d8cabc673a3,
            ]),
            false,
        );
        assert_eq!(a.double(), a_double);

        let mut b = G1Affine::<Bls12381>::new(
            Bls12381::from_raw_unchecked([
                0x1b72cc2215a57793,
                0x6263ee31e953a86d,
                0x866816fd0826ce7b,
                0x991224ec2a74beb2,
                0x041314ca93e1fee0,
                0x19cdf3807146e68e,
            ]),
            Bls12381::from_raw_unchecked([
                0xe7da608505d48616,
                0x22c9f5c77db40989,
                0x35c0752ac9742b45,
                0x41bfaf99f604d1f8,
                0xf45c6e4fc2780554,
                0x7481b1f261aabac,
            ]),
            false,
        );

        let b_double = G1Affine::<Bls12381>::new(
            Bls12381::from_raw_unchecked([
                0xd1e2c01839752ada,
                0x2c4c7d1e2b03e6b1,
                0x5df708511d3f6863,
                0x65f07f9a9b73f98d,
                0xb6e8189b95a60b88,
                0x1252a4ac3529f8b2,
            ]),
            Bls12381::from_raw_unchecked([
                0x77dd65f2e90f5358,
                0x15ac65f21bed27bd,
                0xa9af032870f7bbc6,
                0x18ab5c26dffca63c,
                0x1a49b9965eca3cb8,
                0x2a1bc189e36902d,
            ]),
            false,
        );

        for _ in 0..10 {
            let b = G1Affine::<Bls12381>::random(&mut rand::thread_rng());
            let b1 = b.double() + b;
            let b2 = b1 + b;
            let b3 = b * Fr::from(4);
            assert!(b2 == b.double().double());
            assert!(b2 == b3);
        }
    }

    #[test]
    fn test_double_and_add_arithmetic() {
        for _ in 0..100 {
            let p = G1Affine::<Bls12381>::random(&mut rand::thread_rng());
            let q = G1Affine::<Bls12381>::random(&mut rand::thread_rng());
            let r = G1Affine::<Bls12381>::random(&mut rand::thread_rng());

            let double_p_add_q = p.double() + q;
            let p_plus_q_plus_p = (p + q) + p;

            assert_eq!(p + q, q + p);
            assert_eq!((p + q) + r, p + (q + r));
            assert_eq!(double_p_add_q, p_plus_q_plus_p);
        }
    }

    #[test]
    fn test_scalar_multiplication() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let r: u64 = rng.gen::<u64>() % 100;
            let k = Fr::<Bls12381>::from(r);
            let a = G1Affine::<Bls12381>::random(&mut rng);
            let lhs = &a * &k;
            let rhs = (0..r).fold(G1Affine::<Bls12381>::identity(), |acc, _| acc + &a);
            assert_eq!(lhs, rhs);
        }
    }
}

#[cfg(test)]
mod bn254_g1_affine_test {
    use crate::fp::Bn254;

    use super::*;
    use rand::Rng;

    fn bn254_g1_affine_rand() -> G1Affine<Bn254> {
        let mut rng = rand::thread_rng();
        G1Affine::random(&mut rng)
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = G1Affine::<Bn254>::new(
                Bn254::from_raw_unchecked(x.clone().try_into().unwrap()),
                Bn254::from_raw_unchecked(y.clone().try_into().unwrap()),
                false,
            );

            let b = G1Affine::<Bn254>::new(
                Bn254::from_raw_unchecked(x.clone().try_into().unwrap()),
                Bn254::from_raw_unchecked(y.clone().try_into().unwrap()),
                false,
            );

            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_valid_point() {
        let a = G1Affine::<Bn254>::new(
            Bn254::from_raw_unchecked([
                0x1b72cc2215a57793,
                0x6263ee31e953a86d,
                0x866816fd0826ce7b,
                0x991224ec2a74beb2,
                0x041314ca93e1fee0,
                0x19cdf3807146e68e,
            ]),
            Bn254::from_raw_unchecked([
                0xe7da608505d48616,
                0x22c9f5c77db40989,
                0x35c0752ac9742b45,
                0x41bfaf99f604d1f8,
                0xf45c6e4fc2780554,
                0x7481b1f261aabac,
            ]),
            false,
        );
        assert!(G1Affine::<Bn254>::generator().is_valid().unwrap() == ());
    }

    #[test]
    fn test_double() {
        let a = G1Affine::<Bn254>::new(
            Bn254::from_raw_unchecked([
                0xfb3af00adb22c6bb,
                0x6c55e83ff97a1aef,
                0xa14e3a3f171bac58,
                0xc3688c4f9774b905,
                0x2695638c4fa9ac0f,
                0x17f1d3a73197d794,
            ]),
            Bn254::from_raw_unchecked([
                0x0caa232946c5e7e1,
                0xd03cc744a2888ae4,
                0x00db18cb2c04b3ed,
                0xfcf5e095d5d00af6,
                0xa09e30ed741d8ae4,
                0x08b3f481e3aaa0f1,
            ]),
            false,
        );
        let a_double = G1Affine::<Bn254>::new(
            Bn254::from_raw_unchecked([
                0xc39a8c5529bf0f4e,
                0xe28f75bb8f1c7c42,
                0x43902d0ac358a62a,
                0x9721db3091280125,
                0x8808c8eb50a9450c,
                0x572cbea904d6746,
            ]),
            Bn254::from_raw_unchecked([
                0xba86881979749d28,
                0x4c56d9d4cd16bd1b,
                0xf73bb9021d5fd76a,
                0x22ba3ecb8670e461,
                0x22fda673779d8e38,
                0x166a9d8cabc673a3,
            ]),
            false,
        );
        assert_eq!(a.double(), a_double);

        let mut b = G1Affine::<Bn254>::new(
            Bn254::from_raw_unchecked([
                0x1b72cc2215a57793,
                0x6263ee31e953a86d,
                0x866816fd0826ce7b,
                0x991224ec2a74beb2,
                0x041314ca93e1fee0,
                0x19cdf3807146e68e,
            ]),
            Bn254::from_raw_unchecked([
                0xe7da608505d48616,
                0x22c9f5c77db40989,
                0x35c0752ac9742b45,
                0x41bfaf99f604d1f8,
                0xf45c6e4fc2780554,
                0x7481b1f261aabac,
            ]),
            false,
        );

        let b_double = G1Affine::<Bn254>::new(
            Bn254::from_raw_unchecked([
                0xd1e2c01839752ada,
                0x2c4c7d1e2b03e6b1,
                0x5df708511d3f6863,
                0x65f07f9a9b73f98d,
                0xb6e8189b95a60b88,
                0x1252a4ac3529f8b2,
            ]),
            Bn254::from_raw_unchecked([
                0x77dd65f2e90f5358,
                0x15ac65f21bed27bd,
                0xa9af032870f7bbc6,
                0x18ab5c26dffca63c,
                0x1a49b9965eca3cb8,
                0x2a1bc189e36902d,
            ]),
            false,
        );

        for _ in 0..10 {
            let b = G1Affine::<Bn254>::random(&mut rand::thread_rng());
            let b1 = b.double() + b;
            let b2 = b1 + b;
            let b3 = b * Fr::from(4);
            assert!(b2 == b.double().double());
            assert!(b2 == b3);
        }
    }

    #[test]
    fn test_double_and_add_arithmetic() {
        for _ in 0..100 {
            let p = G1Affine::<Bn254>::random(&mut rand::thread_rng());
            let q = G1Affine::<Bn254>::random(&mut rand::thread_rng());
            let r = G1Affine::<Bn254>::random(&mut rand::thread_rng());

            let double_p_add_q = p.double() + q;
            let p_plus_q_plus_p = (p + q) + p;

            assert_eq!(p + q, q + p);
            assert_eq!((p + q) + r, p + (q + r));
            assert_eq!(double_p_add_q, p_plus_q_plus_p);
        }
    }

    #[test]
    fn test_scalar_multiplication() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let r: u64 = rng.gen::<u64>() % 100;
            let k = Fr::<Bn254>::from(r);
            let a = G1Affine::<Bn254>::random(&mut rng);
            let lhs = &a * &k;
            let rhs = (0..r).fold(G1Affine::<Bn254>::identity(), |acc, _| acc + &a);
            assert_eq!(lhs, rhs);
        }
    }
}
