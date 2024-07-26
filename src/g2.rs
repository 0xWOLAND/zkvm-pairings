use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::common::Curve;
use crate::fp::Fp;
use crate::{fp2::Fp2, fr::Fr};

#[derive(Clone, Copy, Debug)]
struct G2Affine<C: Curve> {
    x: Fp2<C>,
    y: Fp2<C>,
    _marker: PhantomData<C>,
}

impl<C: Curve> PartialEq for G2Affine<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: Curve> G2Affine<C> {
    pub fn new(x: Fp2<C>, y: Fp2<C>) -> Self {
        G2Affine {
            x,
            y,
            _marker: PhantomData::<C>,
        }
    }

    pub const fn zero() -> Self {
        G2Affine {
            x: Fp2::zero(),
            y: Fp2::zero(),
            _marker: PhantomData::<C>,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }

    pub const fn generator() -> Self {
        G2Affine {
            x: Fp2::new(
                Fp::from_raw_unchecked(C::G2_X0),
                Fp::from_raw_unchecked(C::G2_X1),
            ),
            y: Fp2::new(
                Fp::from_raw_unchecked(C::G2_Y0),
                Fp::from_raw_unchecked(C::G2_Y1),
            ),
            _marker: PhantomData::<C>,
        }
    }
}
