use std::marker::PhantomData;
use std::ops::{Mul, MulAssign, Neg};

use crate::common::Curve;
use crate::{fp::Fp, scalar::Scalar};

/// Constant representing the modulus
/// q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
const MODULUS: Scalar = Scalar([
    0xffff_ffff_0000_0001,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
]);

#[derive(Clone, Copy)]
struct G1Affine<Curve> {
    x: Fp,
    y: Fp,
    infinity: bool,
    _marker: PhantomData<Curve>,
}

impl<C: Curve> PartialEq for G1Affine<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.infinity == other.infinity
    }
}

impl<C: Curve> G1Affine<C> {
    pub fn new(x: Fp, y: Fp, infinity: bool) -> Self {
        G1Affine {
            x,
            y,
            infinity,
            _marker: PhantomData::<C>,
        }
    }

    pub const fn zero() -> Self {
        G1Affine {
            x: Fp::one(),
            y: Fp::one(),
            infinity: false,
            _marker: PhantomData::<C>,
        }
    }

    pub fn double(&self) -> Self {
        let x = self.x;
        let y = self.y;

        let slope = (Fp::from(3) * x.square() + C::A) / (Fp::from(2) * y);
        let xr = slope.square() - Fp::from(2) * x;
        let yr = slope * (x - xr) - y;

        Self {
            x: xr,
            y: yr,
            infinity: false,
            _marker: PhantomData::<C>,
        }
    }

    pub fn add(&self, rhs: &Self) -> Self {
        if self.infinity {
            return *rhs;
        }

        if rhs.infinity {
            return *self;
        }

        let x1 = self.x;
        let y1 = self.y;
        let x2 = rhs.x;
        let y2 = rhs.y;

        if x1 == x2 && y1 == y2 {
            return self.double();
        }

        let slope = (y2 - y1) / (x2 - x1);
        let xr = slope.square() - x1 - x2;
        let yr = slope * (x1 - xr) - y1;

        Self {
            x: xr,
            y: yr,
            infinity: false,
            _marker: PhantomData::<C>,
        }
    }

    pub fn mul(&self, rhs: &Scalar) -> Self {
        let mut acc = Self::zero();

        for bit in rhs
            .0
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| (byte >> i) == 1))
            .skip(1)
        {
            acc = acc.double();
            if bit {
                acc = acc.add(self);
            }
        }

        acc
    }
}

impl<C: Curve> Neg for G1Affine<C> {
    type Output = Self;

    fn neg(self) -> Self {
        G1Affine {
            x: self.x,
            y: -self.y,
            infinity: self.infinity,
            _marker: PhantomData::<C>,
        }
    }
}

impl<'a, 'b, C: Curve> Mul<&'b Scalar> for &'a G1Affine<C> {
    type Output = G1Affine<C>;

    #[inline]
    fn mul(self, rhs: &'b Scalar) -> G1Affine<C> {
        self.mul(rhs)
    }
}

impl<'a, 'b, C: Curve> Mul<&'b G1Affine<C>> for &'a Scalar {
    type Output = G1Affine<C>;

    #[inline]
    fn mul(self, rhs: &'b G1Affine<C>) -> G1Affine<C> {
        rhs.mul(self)
    }
}
