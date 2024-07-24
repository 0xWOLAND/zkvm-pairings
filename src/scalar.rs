use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::mem::transmute;

use crate::utils::{adc, sbb};

const MODULUS: Scalar = Scalar([
    0xffff_ffff_0000_0001,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
]);

const LEGENDRE: Scalar = Scalar([
    0x7fffffff80000000,
    0xa9ded2017fff2dff,
    0x199cec0404d0ec02,
    0x39f6d3a994cebea4,
]);

#[derive(Clone, Copy, Eq)]
pub struct Scalar(pub [u64; 4]);

impl Scalar {
    pub const fn zero() -> Self {
        Scalar([0, 0, 0, 0])
    }

    pub const fn one() -> Self {
        Scalar([1, 0, 0, 0])
    }

    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|&e| e == 0)
    }

    fn invert(&self) -> Option<Self> {
        #[inline(always)]
        fn square_assign_multi(n: &mut Scalar, num_times: usize) {
            for _ in 0..num_times {
                *n = n.square();
            }
        }
        // found using https://github.com/kwantam/addchain
        let mut t0 = self.square();
        let mut t1 = t0 * self;
        let mut t16 = t0.square();
        let mut t6 = t16.square();
        let mut t5 = t6 * t0;
        t0 = t6 * t16;
        let mut t12: Scalar = t5 * t16;
        let mut t2 = t6.square();
        let mut t7 = t5 * t6;
        let mut t15 = t0 * t5;
        let mut t17 = t12.square();
        t1 *= t17;
        let mut t3 = t7 * t2;
        let t8 = t1 * t17;
        let t4 = t8 * t2;
        let t9 = t8 * t7;
        t7 = t4 * t5;
        let t11 = t4 * t17;
        t5 = t9 * t17;
        let t14 = t7 * t15;
        let t13 = t11 * t12;
        t12 = t11 * t17;
        t15 *= &t12;
        t16 *= &t15;
        t3 *= &t16;
        t17 *= &t3;
        t0 *= &t17;
        t6 *= &t0;
        t2 *= &t6;
        square_assign_multi(&mut t0, 8);
        t0 *= &t17;
        square_assign_multi(&mut t0, 9);
        t0 *= &t16;
        square_assign_multi(&mut t0, 9);
        t0 *= &t15;
        square_assign_multi(&mut t0, 9);
        t0 *= &t15;
        square_assign_multi(&mut t0, 7);
        t0 *= &t14;
        square_assign_multi(&mut t0, 7);
        t0 *= &t13;
        square_assign_multi(&mut t0, 10);
        t0 *= &t12;
        square_assign_multi(&mut t0, 9);
        t0 *= &t11;
        square_assign_multi(&mut t0, 8);
        t0 *= &t8;
        square_assign_multi(&mut t0, 8);
        t0 *= self;
        square_assign_multi(&mut t0, 14);
        t0 *= &t9;
        square_assign_multi(&mut t0, 10);
        t0 *= &t8;
        square_assign_multi(&mut t0, 15);
        t0 *= &t7;
        square_assign_multi(&mut t0, 10);
        t0 *= &t6;
        square_assign_multi(&mut t0, 8);
        t0 *= &t5;
        square_assign_multi(&mut t0, 16);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 7);
        t0 *= &t4;
        square_assign_multi(&mut t0, 9);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 5);
        t0 *= &t1;
        square_assign_multi(&mut t0, 5);
        t0 *= &t1;

        // CtOption::new(t0, !self.ct_eq(&Self::zero()))
        Some(t0).filter(|_| *self != Self::zero())
    }

    pub fn square(&self) -> Self {
        self * *self
    }

    pub fn pow_vartime(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res *= self;
                }
            }
        }
        res
    }

    pub fn from_bytes(bytes: &[u8; 32]) -> Option<Scalar> {
        let mut tmp = Scalar([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], MODULUS.0[0], 0);
        let (_, borrow) = sbb(tmp.0[1], MODULUS.0[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], MODULUS.0[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], MODULUS.0[3], borrow);

        // If the element is smaller than MODULUS then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        Some(tmp).filter(|_| is_some == 1)
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut res = [0; 32];
        res[0..8].copy_from_slice(&self.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&self.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&self.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&self.0[3].to_le_bytes());

        res
    }

    pub fn add(&self, rhs: &Self) -> Self {
        use num_bigint::BigUint;

        unsafe {
            let modulus = BigUint::from_slice(&transmute::<[u64; 4], [u32; 8]>(MODULUS.0));
            let slice_lhs = transmute::<&[u64; 4], &[u32; 8]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs) % &modulus;
            let rhs = BigUint::from_bytes_le(&rhs.to_bytes()) % &modulus;

            let prod = (lhs + rhs) % modulus;

            let mut prod_slice = prod.to_bytes_le();
            prod_slice.resize(32, 0);
            Scalar::from_bytes(&prod_slice.try_into().unwrap()).unwrap()
        }
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self + &(-rhs)
    }
    pub fn mul(&self, rhs: &Self) -> Self {
        use num_bigint::BigUint;

        unsafe {
            let modulus = BigUint::from_slice(&transmute::<[u64; 4], [u32; 8]>(MODULUS.0));
            let slice_lhs = transmute::<&[u64; 4], &[u32; 8]>(&self.0);
            let lhs = BigUint::from_slice(slice_lhs) % &modulus;
            let rhs = BigUint::from_bytes_le(&rhs.to_bytes()) % &modulus;

            let prod = (lhs * rhs) % modulus;

            let mut prod_slice = prod.to_bytes_le();
            prod_slice.resize(32, 0);
            Scalar::from_bytes(&prod_slice.try_into().unwrap()).unwrap()
        }
    }

    pub fn div(&self, rhs: &Self) -> Self {
        self * rhs.invert().unwrap()
    }

    pub fn div_by_two(&self) -> Scalar {
        let x = self.0;
        let mut result = x;
        let mut carry = 0u64;

        // If x is odd, we add the modulus to make it even
        if x[0] & 1 == 1 {
            let mut overflow = 0u64;
            for i in 0..4 {
                let (res, new_overflow) = adc(result[i], MODULUS.0[i], overflow);
                result[i] = res;
                overflow = new_overflow;
            }
            carry = overflow;
        }

        // Perform the division by 2
        for i in (0..4).rev() {
            let new_carry = result[i] & 1;
            result[i] = (result[i] >> 1) | (carry << 63);
            carry = new_carry;
        }

        Scalar(result)
    }

    fn legendre(&self) -> Scalar {
        self.pow_vartime(&LEGENDRE.0)
    }

    pub fn sqrt(&self) -> Option<Self> {
        if !self.legendre().eq(&Self::one()) {
            return None;
        }

        let p = MODULUS; // Assuming ORDER is a wrapper around [u64; 4]
        let p_plus_one = MODULUS + Self::one();
        let p_minus_one = MODULUS - Self::one();
        let mut q = p_minus_one;
        let mut s = 0u32;

        while q.0[0] & 1 == 0 {
            q = q.div_by_two();
            s += 1;
        }

        if s == 1 {
            let exp = p_plus_one;
            return Some(self.pow_vartime(&exp.div_by_two().0));
        }

        let mut z = Self::from(2u64);
        while z.legendre() != p_minus_one {
            z += &Self::one();
        }

        let mut c = z.pow_vartime(&q.0);
        let mut r = self.pow_vartime(&p_minus_one.div_by_two().0);
        let mut t = self.pow_vartime(&q.0);

        while !(t - &Self::one()).is_zero() {
            let mut t2 = t.square();
            let mut i = 1;

            while i < s {
                if (t - &Self::one()).is_zero() {
                    break;
                }
                t2 = t2.square();
                i += 1;
            }

            let b = c.pow_vartime(
                &(Scalar::from((s - i - 1) as u64))
                    .pow_vartime(&[2, 0, 0, 0])
                    .0,
            );
            r = r.mul(&b);
            c = b.square();
            t = t.mul(&c);
            s = i;
        }

        Some(r)
    }
}

impl PartialEq for Scalar {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Neg for &Scalar {
    type Output = Scalar;

    fn neg(self) -> Scalar {
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sbb(MODULUS.0[0], self.0[0], 0);
        let (d1, borrow) = sbb(MODULUS.0[1], self.0[1], borrow);
        let (d2, borrow) = sbb(MODULUS.0[2], self.0[2], borrow);
        let (d3, _) = sbb(MODULUS.0[3], self.0[3], borrow);

        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Scalar([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }
}

impl<'a, 'b> Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn add(self, rhs: &'b Scalar) -> Scalar {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn sub(self, rhs: &'b Scalar) -> Scalar {
        self.sub(rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn mul(self, rhs: &'b Scalar) -> Scalar {
        self.mul(rhs)
    }
}

impl<'a, 'b> Div<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn div(self, rhs: &'b Scalar) -> Scalar {
        self.div(rhs)
    }
}

impl From<u64> for Scalar {
    fn from(val: u64) -> Scalar {
        Scalar([val, 0, 0, 0])
    }
}

impl_binops_additive!(Scalar, Scalar);
impl_binops_multiplicative!(Scalar, Scalar);
impl_binops_divisible!(Scalar, Scalar);
