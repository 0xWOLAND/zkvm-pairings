use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::common::{AffinePoint, Curve};
use crate::{fp::Fp, fr::Fr};

#[derive(Clone, Copy, Debug)]
pub struct G1Affine<C: Curve> {
    x: Fp<C>,
    y: Fp<C>,
    is_infinity: bool,
}

impl<C: Curve> PartialEq for G1Affine<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: Curve> AffinePoint<C> for G1Affine<C> {
    type Dtype = Fp<C>;

    fn new(x: Self::Dtype, y: Self::Dtype, is_infinity: bool) -> Self {
        G1Affine { x, y, is_infinity }
    }
    fn identity() -> Self {
        G1Affine {
            x: Fp::zero(),
            y: Fp::one(),
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
            x: Fp::from_raw_unchecked(C::G1_X),
            y: Fp::from_raw_unchecked(C::G1_Y),
            is_infinity: false,
        }
    }

    fn is_valid(&self) -> Result<(), String> {
        if self.is_identity() {
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
        let b = Fp::from_raw_unchecked(C::B);
        loop {
            let x = Fp::random(&mut rng);
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

        let slope = (Fp::from(3) * x.square()) / (Fp::from(2) * y);
        let xr = slope.square() - Fp::from(2) * x;
        let yr = slope * (x - xr) - y;

        Self {
            x: xr,
            y: yr,
            is_infinity: false,
        }
    }
}

impl<C: Curve> G1Affine<C> {
    pub(crate) fn is_on_curve(&self) -> bool {
        let x = self.x;
        let y = self.y;

        // y^2 = x^3 + B
        y.square() == x.square() * x + Fp::from_raw_unchecked(C::B)
    }

    fn endomorphism(&self) -> Self {
        G1Affine::new(self.x * Fp::from_raw_unchecked(C::BETA), self.y, false) // BETA is a nontrivial third root of unity in Fp
    }

    fn mul_by_x(&self) -> Self {
        self * Fr::from(C::X)
    }

    fn is_torsion_free(&self) -> bool {
        let lhs = self.mul_by_x().mul_by_x().neg();
        let rhs = self.endomorphism();
        lhs == rhs
    }
}

impl<'a, C: Curve> Neg for &'a G1Affine<C> {
    type Output = G1Affine<C>;

    fn neg(self) -> G1Affine<C> {
        G1Affine {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Fr<C>> for &'a G1Affine<C> {
    type Output = G1Affine<C>;

    #[inline]
    fn mul(self, other: &'b Fr<C>) -> G1Affine<C> {
        let mut xself = G1Affine::<C>::identity();
        let mut acc = *self;

        for bit in other
            .0
            .iter()
            .flat_map(|&x| (0..64).map(move |i| (x >> i) & 1 == 1))
            .skip(1)
        {
            acc = acc.double();

            if bit {
                xself += &acc;
            }
        }

        xself
    }
}

impl<'a, 'b, C: Curve> Add<&'b G1Affine<C>> for &'a G1Affine<C> {
    type Output = G1Affine<C>;

    #[inline]
    fn add(self, other: &'b G1Affine<C>) -> G1Affine<C> {
        if self.is_infinity {
            println!("self is infinity");
            return *other;
        }

        if other.is_infinity {
            println!("other is infinity");
            return *self;
        }

        let x1 = self.x;
        let y1 = self.y;
        let x2 = other.x;
        let y2 = other.y;

        if x1 == x2 && y1 == y2 {
            println!("doubling");
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

impl<'a, 'b, C: Curve> Sub<&'b G1Affine<C>> for &'a G1Affine<C> {
    type Output = G1Affine<C>;

    #[inline]
    fn sub(self, other: &'b G1Affine<C>) -> G1Affine<C> {
        self + (-other)
    }
}

impl_binops_multiplicative!(G1Affine<C>, Fr<C>);
impl_binops_additive!(G1Affine<C>, G1Affine<C>);

#[cfg(test)]
mod test {

    use crate::common::Bls12381Curve;
    use rand::Rng;

    use super::*;

    fn fp2_rand() -> G1Affine<Bls12381Curve> {
        let mut rng = rand::thread_rng();
        G1Affine::random(&mut rng)
    }

    #[test]
    fn test_equality() {
        let rng = &mut rand::thread_rng();
        for _ in 0..10 {
            let x = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();
            let y = (0..6).map(|_| rng.gen::<u64>()).collect::<Vec<_>>();

            let a = G1Affine::<Bls12381Curve>::new(
                Fp::from_raw_unchecked(x.clone().try_into().unwrap()),
                Fp::from_raw_unchecked(y.clone().try_into().unwrap()),
                false,
            );

            let b = G1Affine::<Bls12381Curve>::new(
                Fp::from_raw_unchecked(x.clone().try_into().unwrap()),
                Fp::from_raw_unchecked(y.clone().try_into().unwrap()),
                false,
            );

            assert_eq!(a, b)
        }
    }

    #[test]
    fn test_valid_point() {
        let a = G1Affine::<Bls12381Curve>::new(
            Fp::from_raw_unchecked([
                0x1b72cc2215a57793,
                0x6263ee31e953a86d,
                0x866816fd0826ce7b,
                0x991224ec2a74beb2,
                0x041314ca93e1fee0,
                0x19cdf3807146e68e,
            ]),
            Fp::from_raw_unchecked([
                0xe7da608505d48616,
                0x22c9f5c77db40989,
                0x35c0752ac9742b45,
                0x41bfaf99f604d1f8,
                0xf45c6e4fc2780554,
                0x7481b1f261aabac,
            ]),
            false,
        );
        assert!(G1Affine::<Bls12381Curve>::generator().is_valid().unwrap() == ());
    }

    #[test]
    fn test_double() {
        let a = G1Affine::<Bls12381Curve>::new(
            Fp::from_raw_unchecked([
                0xfb3af00adb22c6bb,
                0x6c55e83ff97a1aef,
                0xa14e3a3f171bac58,
                0xc3688c4f9774b905,
                0x2695638c4fa9ac0f,
                0x17f1d3a73197d794,
            ]),
            Fp::from_raw_unchecked([
                0x0caa232946c5e7e1,
                0xd03cc744a2888ae4,
                0x00db18cb2c04b3ed,
                0xfcf5e095d5d00af6,
                0xa09e30ed741d8ae4,
                0x08b3f481e3aaa0f1,
            ]),
            false,
        );
        let a_double = G1Affine::<Bls12381Curve>::new(
            Fp::from_raw_unchecked([
                0xc39a8c5529bf0f4e,
                0xe28f75bb8f1c7c42,
                0x43902d0ac358a62a,
                0x9721db3091280125,
                0x8808c8eb50a9450c,
                0x572cbea904d6746,
            ]),
            Fp::from_raw_unchecked([
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

        let mut b = G1Affine::<Bls12381Curve>::new(
            Fp::from_raw_unchecked([
                0x1b72cc2215a57793,
                0x6263ee31e953a86d,
                0x866816fd0826ce7b,
                0x991224ec2a74beb2,
                0x041314ca93e1fee0,
                0x19cdf3807146e68e,
            ]),
            Fp::from_raw_unchecked([
                0xe7da608505d48616,
                0x22c9f5c77db40989,
                0x35c0752ac9742b45,
                0x41bfaf99f604d1f8,
                0xf45c6e4fc2780554,
                0x7481b1f261aabac,
            ]),
            false,
        );

        let b_double = G1Affine::<Bls12381Curve>::new(
            Fp::from_raw_unchecked([
                0xd1e2c01839752ada,
                0x2c4c7d1e2b03e6b1,
                0x5df708511d3f6863,
                0x65f07f9a9b73f98d,
                0xb6e8189b95a60b88,
                0x1252a4ac3529f8b2,
            ]),
            Fp::from_raw_unchecked([
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
            let b = G1Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
            let b1 = b.double() + b;
            let b2 = b1 + b;
            let b3 = b * Fr::from(4);
            assert!(b2 == b.double().double());
            assert!(b2 == b3);
        }
    }

    #[test]
    fn test_double_and_add_arithmetic() {
        let p = G1Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        let q = G1Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        let r = G1Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

        let double_p_add_q = p.double() + q;
        let p_plus_q_plus_p = (p + q) + p;

        assert_eq!(p + q, q + p);
        assert_eq!((p + q) + r, p + (q + r));
        assert_eq!(double_p_add_q, p_plus_q_plus_p);
    }
}
