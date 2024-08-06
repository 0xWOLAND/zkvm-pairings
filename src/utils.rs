/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

macro_rules! impl_add_binop_specify_output {
    ($lhs:ty, $rhs:ty, $output:ty) => {
        impl<'b, F: FpElement> Add<&'b $rhs> for $lhs {
            type Output = $output;
            #[inline]
            fn add(self, rhs: &'b $rhs) -> $output {
                <&$lhs>::add(&self, rhs)
            }
        }
        impl<'a, F: FpElement> Add<$rhs> for &'a $lhs {
            type Output = $output;
            #[inline]
            fn add(self, rhs: $rhs) -> $output {
                self.add(&rhs)
            }
        }
        impl<F: FpElement> Add<$rhs> for $lhs {
            type Output = $output;
            #[inline]
            fn add(self, rhs: $rhs) -> $output {
                <&$lhs>::add(&self, &rhs)
            }
        }
    };
}

macro_rules! impl_sub_binop_specify_output {
    ($lhs:ty, $rhs:ty, $output:ty) => {
        impl<'b, F: FpElement> Sub<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: &'b $rhs) -> $output {
                &self - rhs
            }
        }

        impl<'a, F: FpElement> Sub<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: $rhs) -> $output {
                self - &rhs
            }
        }

        impl<F: FpElement> Sub<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: $rhs) -> $output {
                &self - &rhs
            }
        }
    };
}

macro_rules! impl_binops_additive_specify_output {
    ($lhs:ty, $rhs:ty, $output:ty) => {
        impl_add_binop_specify_output!($lhs, $rhs, $output);
        impl_sub_binop_specify_output!($lhs, $rhs, $output);
    };
}

macro_rules! impl_binops_divisible_mixed {
    ($lhs:ty, $rhs:ty, $output:ty) => {
        impl<'b, F: FpElement> Div<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn div(self, rhs: &'b $rhs) -> $output {
                &self / rhs
            }
        }

        impl<'a, F: FpElement> Div<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn div(self, rhs: $rhs) -> $output {
                self / &rhs
            }
        }

        impl<F: FpElement> Div<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn div(self, rhs: $rhs) -> $output {
                &self / &rhs
            }
        }
    };
}

macro_rules! impl_binops_multiplicative_mixed {
    ($lhs:ty, $rhs:ty, $output:ty) => {
        impl<'b, F: FpElement> Mul<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: &'b $rhs) -> $output {
                &self * rhs
            }
        }

        impl<'a, F: FpElement> Mul<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: $rhs) -> $output {
                self * &rhs
            }
        }

        impl<F: FpElement> Mul<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: $rhs) -> $output {
                &self * rhs
            }
        }
    };
}

macro_rules! impl_binops_additive {
    ($lhs:ty, $rhs:ty) => {
        impl_binops_additive_specify_output!($lhs, $rhs, $lhs);

        impl<F: FpElement> SubAssign<$rhs> for $lhs {
            #[inline]
            fn sub_assign(&mut self, rhs: $rhs) {
                *self = &*self - &rhs;
            }
        }

        impl<F: FpElement> AddAssign<$rhs> for $lhs {
            #[inline]
            fn add_assign(&mut self, rhs: $rhs) {
                *self = &*self + &rhs;
            }
        }

        impl<'b, F: FpElement> SubAssign<&'b $rhs> for $lhs {
            #[inline]
            fn sub_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self - rhs;
            }
        }

        impl<'b, F: FpElement> AddAssign<&'b $rhs> for $lhs {
            #[inline]
            fn add_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self + rhs;
            }
        }
    };
}

macro_rules! impl_binops_multiplicative {
    ($lhs:ty, $rhs:ty) => {
        impl_binops_multiplicative_mixed!($lhs, $rhs, $lhs);

        impl<F: FpElement> MulAssign<$rhs> for $lhs {
            #[inline]
            fn mul_assign(&mut self, rhs: $rhs) {
                *self = &*self * &rhs;
            }
        }

        impl<'b, F: FpElement> MulAssign<&'b $rhs> for $lhs {
            #[inline]
            fn mul_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self * rhs;
            }
        }
    };
}

macro_rules! impl_binops_divisible {
    ($lhs:ty, $rhs:ty) => {
        impl_binops_divisible_mixed!($lhs, $rhs, $lhs);

        impl<F: FpElement> DivAssign<$rhs> for $lhs {
            #[inline]
            fn div_assign(&mut self, rhs: $rhs) {
                *self = &*self / &rhs;
            }
        }

        impl<'b, F: FpElement> DivAssign<&'b $rhs> for $lhs {
            #[inline]
            fn div_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self / rhs;
            }
        }
    };
}
