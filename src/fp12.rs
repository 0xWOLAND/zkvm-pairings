use crate::fp::*;
use crate::fp2::*;
use crate::fp6::*;

use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

use num_bigint::BigUint;
use rand::RngCore;

#[cfg(target_os = "zkvm")]
use sp1_lib::{
    io::{self, hint_slice},
    unconstrained,
};

pub trait Fp12Element: Fp6Element {
    type Fp12ElementType;
    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp12ElementType;
    fn to_bytes_vec(value: &Self::Fp12ElementType) -> Vec<u8>;
    fn get_fp12_frobenius_coeffs(pow: usize) -> (Self, Self);
    fn _invert(f: &Fp12<Self>) -> Option<Fp12<Self>> {
        (f.c0.square() - f.c1.square().mul_by_nonresidue())
            .invert()
            .map(|t| Fp12::new(f.c0 * t, f.c1 * -t))
    }
    fn invert(f: &Fp12<Self>) -> Option<Fp12<Self>>;
}

impl Fp12Element for Bls12381 {
    type Fp12ElementType = Fp12<Bls12381>;
    fn get_fp12_frobenius_coeffs(pow: usize) -> (Self, Self) {
        match pow % 12 {
            0 => (Self::one(), Self::one()),
            1 => (
                Self::from_raw_unchecked([
                    0x8d0775ed92235fb8,
                    0xf67ea53d63e7813d,
                    0x7b2443d784bab9c4,
                    0x0fd603fd3cbd5f4f,
                    0xc231beb4202c0d1f,
                    0x1904d3bf02bb0667,
                ]),
                Self::from_raw_unchecked([
                    0x2cf78a126ddc4af3,
                    0x282d5ac14d6c7ec2,
                    0xec0c8ec971f63c5f,
                    0x54a14787b6c7b36f,
                    0x88e9e902231f9fb8,
                    0x00fc3e2b36c4e032,
                ]),
            ),
            2 => (
                Self::zero(),
                Self::from_raw_unchecked([
                    0x8bfd00000000aaac,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
            ),
            3 => (
                Self::from_raw_unchecked([
                    0xc81084fbede3cc09,
                    0xee67992f72ec05f4,
                    0x77f76e17009241c5,
                    0x48395dabc2d3435e,
                    0x6831e36d6bd17ffe,
                    0x06af0e0437ff400b,
                ]),
                Self::from_raw_unchecked([
                    0xc81084fbede3cc09,
                    0xee67992f72ec05f4,
                    0x77f76e17009241c5,
                    0x48395dabc2d3435e,
                    0x6831e36d6bd17ffe,
                    0x06af0e0437ff400b,
                ]),
            ),
            4 => (
                Self::from_raw_unchecked([
                    0x8bfd00000000aaad,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
                Self::zero(),
            ),
            5 => (
                Self::from_raw_unchecked([
                    0x9b18fae980078116,
                    0xc63a3e6e257f8732,
                    0x8beadf4d8e9c0566,
                    0xf39816240c0b8fee,
                    0xdf47fa6b48b1e045,
                    0x05b2cfd9013a5fd8,
                ]),
                Self::from_raw_unchecked([
                    0x1ee605167ff82995,
                    0x5871c1908bd478cd,
                    0xdb45f3536814f0bd,
                    0x70df3560e77982d0,
                    0x6bd3ad4afa99cc91,
                    0x144e4211384586c1,
                ]),
            ),
            6 => (
                Self::zero(),
                Self::from_raw_unchecked([
                    0xb9feffffffffaaaa,
                    0x1eabfffeb153ffff,
                    0x6730d2a0f6b0f624,
                    0x64774b84f38512bf,
                    0x4b1ba7b6434bacd7,
                    0x1a0111ea397fe69a,
                ]),
            ),
            7 => (
                Self::from_raw_unchecked([
                    0x2cf78a126ddc4af3,
                    0x282d5ac14d6c7ec2,
                    0xec0c8ec971f63c5f,
                    0x54a14787b6c7b36f,
                    0x88e9e902231f9fb8,
                    0x00fc3e2b36c4e032,
                ]),
                Self::from_raw_unchecked([
                    0x2cf78a126ddc4af3,
                    0x282d5ac14d6c7ec2,
                    0xec0c8ec971f63c5f,
                    0x54a14787b6c7b36f,
                    0x88e9e902231f9fb8,
                    0x00fc3e2b36c4e032,
                ]),
            ),
            8 => (
                Self::from_raw_unchecked([
                    0x8bfd00000000aaac,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
                Self::zero(),
            ),
            9 => (
                Self::from_raw_unchecked([
                    0x8bfd00000000aaac,
                    0x409427eb4f49fffd,
                    0x897d29650fb85f9b,
                    0xaa0d857d89759ad4,
                    0xec02408663d4de85,
                    0x1a0111ea397fe699,
                ]),
                Self::zero(),
            ),
            10 => (
                Self::from_raw_unchecked([
                    0xc81084fbede3cc09,
                    0xee67992f72ec05f4,
                    0x77f76e17009241c5,
                    0x48395dabc2d3435e,
                    0x6831e36d6bd17ffe,
                    0x06af0e0437ff400b,
                ]),
                Self::from_raw_unchecked([
                    0xf1ee7b04121bdea2,
                    0x304466cf3e67fa0a,
                    0xef396489f61eb45e,
                    0x1c3dedd930b1cf60,
                    0xe2e9c448d77a2cd9,
                    0x135203e60180a68e,
                ]),
            ),
            11 => (
                Self::zero(),
                Self::from_raw_unchecked([
                    0x2e01fffffffefffe,
                    0xde17d813620a0002,
                    0xddb3a93be6f89688,
                    0xba69c6076a0f77ea,
                    0x5f19672fdf76ce51,
                    0x0000000000000000,
                ]),
            ),
            _ => unimplemented!(),
        }
    }

    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp12ElementType {
        let c0 = <Bls12381 as Fp6Element>::from_bytes_slice(&bytes[..288]);
        let c1 = <Bls12381 as Fp6Element>::from_bytes_slice(&bytes[288..]);
        Fp12::<Bls12381>::new(c0, c1)
    }

    fn to_bytes_vec(value: &Self::Fp12ElementType) -> Vec<u8> {
        let mut res = [0u8; 576];
        let c0 = <Bls12381 as Fp6Element>::to_bytes_vec(&value.c0);
        let c1 = <Bls12381 as Fp6Element>::to_bytes_vec(&value.c1);

        res[..288].copy_from_slice(&c0);
        res[288..].copy_from_slice(&c1);

        res.to_vec()
    }

    #[cfg(not(target_os = "zkvm"))]
    fn invert(f: &Fp12<Self>) -> Option<Fp12<Self>> {
        Fp12Element::_invert(f)
    }

    #[cfg(target_os = "zkvm")]
    fn invert(f: &Fp12<Self>) -> Option<Fp12<Self>> {
        unconstrained! {
            let mut buf = [0u8; 577];
            match Fp12Element::_invert(&f) {
                Some(x) => {
                    buf[576] = 1;
                    buf[0..576].copy_from_slice(&<Self as Fp12Element>::to_bytes_vec(&x));
                }
                None => {}
            }
            hint_slice(&buf);
        }

        let bytes: [u8; 577] = io::read_vec().try_into().unwrap();
        let is_some = bytes[576] == 1;
        let bytes = bytes[..576].try_into().unwrap();
        let out = <Self as Fp12Element>::from_bytes_slice(bytes);

        Some(out).filter(|_| is_some)
    }
}

impl Fp12Element for Bn254 {
    type Fp12ElementType = Fp12<Bn254>;
    fn get_fp12_frobenius_coeffs(pow: usize) -> (Self, Self) {
        match pow % 12 {
            0 => (Self::one(), Self::zero()),
            1 => (
                Self::from_raw_unchecked([
                    0xd60b35dadcc9e470,
                    0x5c521e08292f2176,
                    0xe8b99fdd76e68b60,
                    0x1284b71c2865a7df,
                ]),
                Self::from_raw_unchecked([
                    0xca5cf05f80f362ac,
                    0x747992778eeec7e5,
                    0xa6327cfe12150b8e,
                    0x246996f3b4fae7e6,
                ]),
            ),
            2 => (
                Self::from_raw_unchecked([
                    0xe4bd44e5607cfd49,
                    0xc28f069fbb966e3d,
                    0x5e6dd9e7e0acccb0,
                    0x30644e72e131a029,
                ]),
                Self::zero(),
            ),
            3 => (
                Self::from_raw_unchecked([
                    0xe86f7d391ed4a67f,
                    0x894cb38dbe55d24a,
                    0xefe9608cd0acaa90,
                    0x19dc81cfcc82e4bb,
                ]),
                Self::from_raw_unchecked([
                    0x7694aa2bf4c0c101,
                    0x7f03a5e397d439ec,
                    0x06cbeee33576139d,
                    0xabf8b60be77d73,
                ]),
            ),
            4 => (
                Self::from_raw_unchecked([
                    0xe4bd44e5607cfd48,
                    0xc28f069fbb966e3d,
                    0x5e6dd9e7e0acccb0,
                    0x30644e72e131a029,
                ]),
                Self::zero(),
            ),
            5 => (
                Self::from_raw_unchecked([
                    0x1264475e420ac20f,
                    0x2cfa95859526b0d4,
                    0x072fc0af59c61f30,
                    0x757cab3a41d3cdc,
                ]),
                Self::from_raw_unchecked([
                    0xe85845e34c4a5b9c,
                    0xa20b7dfd71573c93,
                    0x18e9b79ba4e2606c,
                    0xca6b035381e35b6,
                ]),
            ),
            6 => (
                Self::from_raw_unchecked([
                    0x3c208c16d87cfd46,
                    0x97816a916871ca8d,
                    0xb85045b68181585d,
                    0x30644e72e131a029,
                ]),
                Self::zero(),
            ),
            7 => (
                Self::from_raw_unchecked([
                    0x6615563bfbb318d7,
                    0x3b2f4c893f42a916,
                    0xcf96a5d90a9accfd,
                    0x1ddf9756b8cbf849,
                ]),
                Self::from_raw_unchecked([
                    0x71c39bb757899a9b,
                    0x2307d819d98302a7,
                    0x121dc8b86f6c4ccf,
                    0xbfab77f2c36b843,
                ]),
            ),
            8 => (
                Self::from_raw_unchecked([
                    0x0,
                    0x5763473177fffffe,
                    0xd4f263f1acdb5c4f,
                    0x59e26bcea0d48bac,
                ]),
                Self::zero(),
            ),
            9 => (
                Self::from_raw_unchecked([
                    0x53b10eddb9a856c8,
                    0x0e34b703aa1bf842,
                    0xc866e529b0d4adcd,
                    0x1687cca314aebb6d,
                ]),
                Self::from_raw_unchecked([
                    0xc58be1eae3bc3c46,
                    0x187dc4add09d90a0,
                    0xb18456d34c0b44c0,
                    0x2fb855bcd54a22b6,
                ]),
            ),
            10 => (
                Self::from_raw_unchecked([
                    0x0,
                    0x5763473177ffffff,
                    0xd4f263f1acdb5c4f,
                    0x59e26bcea0d48bac,
                ]),
                Self::zero(),
            ),
            11 => (
                Self::from_raw_unchecked([
                    0x29bc44b896723b38,
                    0x6a86d50bd34b19b9,
                    0xb120850727bb392d,
                    0x290c83bf3d14634d,
                ]),
                Self::from_raw_unchecked([
                    0x53c846338c32a1ab,
                    0xf575ec93f71a8df9,
                    0x9f668e1adc9ef7f0,
                    0x23bd9e3da9136a73,
                ]),
            ),

            _ => unimplemented!(),
        }
    }

    fn from_bytes_slice(bytes: &[u8]) -> Self::Fp12ElementType {
        let c0 = <Bn254 as Fp6Element>::from_bytes_slice(&bytes[..288]);
        let c1 = <Bn254 as Fp6Element>::from_bytes_slice(&bytes[288..]);
        Fp12::<Bn254>::new(c0, c1)
    }

    fn to_bytes_vec(value: &Self::Fp12ElementType) -> Vec<u8> {
        let mut res = [0u8; 576];
        let c0 = <Bn254 as Fp6Element>::to_bytes_vec(&value.c0);
        let c1 = <Bn254 as Fp6Element>::to_bytes_vec(&value.c1);

        res[..288].copy_from_slice(&c0);
        res[288..].copy_from_slice(&c1);

        res.to_vec()
    }

    #[cfg(not(target_os = "zkvm"))]
    fn invert(f: &Fp12<Self>) -> Option<Fp12<Self>> {
        Fp12Element::_invert(f)
    }

    #[cfg(target_os = "zkvm")]
    fn invert(f: &Fp12<Self>) -> Option<Fp12<Self>> {
        unconstrained! {
            let mut buf = [0u8; 513];
            match Fp12Element::_invert(&f) {
                Some(x) => {
                    buf[512] = 1;
                    buf[0..512].copy_from_slice(&<Self as Fp12Element>::to_bytes_vec(&x));
                }
                None => {}
            }
            hint_slice(&buf);
        }

        let bytes: [u8; 513] = io::read_vec().try_into().unwrap();
        let is_some = bytes[512] == 1;
        let bytes = bytes[..512].try_into().unwrap();
        let out = <Self as Fp12Element>::from_bytes_slice(bytes);

        Some(out).filter(|_| is_some)
    }
}

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
pub struct Fp12<F: Fp12Element> {
    pub c0: Fp6<F>,
    pub c1: Fp6<F>,
}

impl<F: Fp12Element> From<F> for Fp12<F> {
    fn from(f: F) -> Fp12<F> {
        Fp12 {
            c0: Fp6::<F>::from(f),
            c1: Fp6::<F>::from(f),
        }
    }
}

impl<F: Fp12Element> From<Fp2<F>> for Fp12<F> {
    fn from(f: Fp2<F>) -> Fp12<F> {
        Fp12 {
            c0: Fp6::<F>::from(f),
            c1: Fp6::<F>::zero(),
        }
    }
}

impl<F: Fp12Element> From<Fp6<F>> for Fp12<F> {
    fn from(f: Fp6<F>) -> Fp12<F> {
        Fp12 {
            c0: f,
            c1: Fp6::<F>::zero(),
        }
    }
}

impl<F: Fp12Element> Eq for Fp12<F> {}
impl<F: Fp12Element> PartialEq for Fp12<F> {
    fn eq(&self, other: &Fp12<F>) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<F: Fp12Element> Copy for Fp12<F> {}
impl<F: Fp12Element> Clone for Fp12<F> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<F: Fp12Element> Default for Fp12<F> {
    fn default() -> Self {
        Fp12::<F>::zero()
    }
}

#[cfg(feature = "zeroize")]
impl<F: Fp12Element> zeroize::DefaultIsZeroes for Fp12 {}

impl<F: Fp12Element> fmt::Debug for Fp12<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*w", self.c0, self.c1)
    }
}

impl<F: Fp12Element> Fp12<F> {
    #[inline]
    pub fn new(c0: Fp6<F>, c1: Fp6<F>) -> Self {
        Fp12 { c0, c1 }
    }

    #[inline]
    pub fn zero() -> Self {
        Fp12::new(Fp6::zero(), Fp6::zero())
    }

    #[inline]
    pub fn one() -> Self {
        Fp12::new(Fp6::one(), Fp6::zero())
    }

    pub fn random(mut rng: impl RngCore) -> Self {
        Fp12 {
            c0: Fp6::<F>::random(&mut rng),
            c1: Fp6::<F>::random(&mut rng),
        }
    }

    pub fn mul_by_014(&self, c0: &Fp2<F>, c1: &Fp2<F>, c4: &Fp2<F>) -> Fp12<F> {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = *c1 + *c4;
        let c1 = self.c1 + self.c0;
        let c1 = c1.mul_by_01(c0, &o);
        let c1 = c1 - aa - bb;
        let c0 = bb;
        let c0 = c0.mul_by_nonresidue();
        let c0 = c0 + aa;

        Fp12 { c0, c1 }
    }

    pub fn mul_14_by_14(d0: &Fp2<F>, d1: &Fp2<F>, c0: &Fp2<F>, c1: &Fp2<F>) -> [Fp2<F>; 5] {
        let x0 = *d0 * *c0;
        let x1 = *d1 * *c1;
        let x04 = *c0 + *d0;
        let tmp = *c0 + *c1;
        let x01 = *d0 + *d1;
        let x01 = x01 * tmp;
        let tmp = x1 + x0;
        let x01 = x01 - tmp;
        let x14 = *c1 + *d1;
        let z_c0_b0 = Fp2::<F>::non_residue() + x0;

        [z_c0_b0, x01, x1, x04, x14]
    }

    fn cyclotomic_square(&self) -> Fp12<F> {
        let t0 = self.c1.c1.square();
        let t1 = self.c0.c0.square();
        let t6 = (self.c1.c1 + self.c0.c0).square();
        let t6 = t6 - t0;
        let t6 = t6 - t1;
        let t2 = &self.c0.c2.square();
        let t3 = &self.c1.c0.square();
        let t7 = (self.c0.c2 + self.c1.c0).square();
        let t7 = t7 - *t2;
        let t7 = t7 - *t3;
        let t4 = self.c1.c2.square();
        let t5 = self.c0.c1.square();
        let t8 = (self.c1.c2 + self.c0.c1).square();
        let t8 = t8 - t4;
        let t8 = t8 - t5;
        let t8 = t8.mul_by_nonresidue();
        let t0 = t0.mul_by_nonresidue();
        let t0 = t0 + t1;
        let t2 = t2.mul_by_nonresidue();
        let t2 = t2 + *t3;
        let t4 = t4.mul_by_nonresidue();
        let t4 = t4 + t5;
        let z00 = t0 - self.c0.c0;
        let z00 = z00 + z00;
        let z00 = z00 + t0;
        let z01 = t2 - self.c0.c1;
        let z01 = z01 + z01;
        let z01 = z01 + t2;
        let z02 = t4 - self.c0.c2;
        let z02 = z02 + z02;
        let z02 = z02 + t4;
        let z10 = t8 + self.c1.c0;
        let z10 = z10 + z10;
        let z10 = z10 + t8;
        let z11 = t6 + self.c1.c1;
        let z11 = z11 + z11;
        let z11 = z11 + t6;
        let z12 = t7 + self.c1.c2;
        let z12 = z12 + z12;
        let z12 = z12 + t7;
        Fp12::new(Fp6::new(z00, z01, z02), Fp6::new(z10, z11, z12))
    }

    fn n_cyclotomic_square(&self, by: u64) -> Fp12<F> {
        (0..by).fold(*self, |acc, _| acc.cyclotomic_square())
    }

    pub fn powt(&self) -> Fp12<F> {
        let a = self.cyclotomic_square();
        let a = a * *self;
        let a = a.n_cyclotomic_square(2);
        let a = a * *self;
        let a = a.n_cyclotomic_square(3);
        let a = a * *self;
        let a = a.n_cyclotomic_square(9);
        let a = a * *self;
        let a = a.n_cyclotomic_square(32);
        let a = a * *self;
        let a = a.n_cyclotomic_square(15);
        let a = a * *self;
        a.cyclotomic_square()
    }

    pub fn div(&self, rhs: &Fp12<F>) -> Fp12<F> {
        rhs.invert().unwrap() * *self
    }

    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    #[inline(always)]
    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    #[inline(always)]
    pub fn conjugate(&self) -> Self {
        Fp12::new(self.c0, -self.c1)
    }

    pub fn pow_vartime(&self, by: &[u64; 6]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res = res * *self;
                }
            }
        }
        res
    }

    pub fn pow_vartime_extended(&self, by: &[u64]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res = res * *self;
                }
            }
        }
        res
    }

    pub fn pow_vartime_extended_str(&self, by: &str) -> Self {
        self.pow_vartime_extended(&BigUint::from_str(by).unwrap().to_u64_digits())
    }
    /// Raises this element to p.
    #[inline(always)]
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();

        let frob_coeffs = F::get_fp12_frobenius_coeffs(1);
        let c1 = c1 * Fp6::from(Fp2::new(frob_coeffs.0, frob_coeffs.1));

        Fp12::new(c0, c1)
    }

    // #[inline(always)]
    // pub(crate) fn nth_frobenius_map(&self, pow: usize) -> Self {
    //     let c0 = self.c0.nth_frobenius_map(pow);
    //     let c1 = self.c1.nth_frobenius_map(pow);

    //     let frob_coeffs = F::get_fp12_frobenius_coeffs(pow);
    //     let c1 = c1 * Fp6::from(Fp2::new(frob_coeffs[0], frob_coeffs[1]));

    //     Fp12::new(c0, c1)
    // }

    #[inline]
    pub fn square(&self) -> Self {
        let ab = self.c0 * self.c1;
        let c0c1 = self.c0 + self.c1;
        let c0 = self.c1.mul_by_nonresidue();
        let c0 = c0 + self.c0;
        let c0 = c0 * c0c1;
        let c0 = c0 - ab;
        let c1 = ab + ab;
        let c0 = c0 - ab.mul_by_nonresidue();

        Fp12::new(c0, c1)
    }

    pub fn invert(&self) -> Option<Self> {
        (self.c0.square() - self.c1.square().mul_by_nonresidue())
            .invert()
            .map(|t| Fp12::new(self.c0 * t, self.c1 * -t))
    }
}

impl<F: Fp12Element> Mul<Fp12<F>> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn mul(self, other: Fp12<F>) -> Self::Output {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;
        let o = other.c0 + other.c1;
        let c1 = self.c1 + self.c0;
        let c1 = c1 * o;
        let c1 = c1 - aa;
        let c1 = c1 - bb;
        let c0 = bb.mul_by_nonresidue();
        let c0 = c0 + aa;

        Fp12::new(c0, c1)
    }
}

impl<F: Fp12Element> Add<Fp12<F>> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn add(self, rhs: Fp12<F>) -> Self::Output {
        Fp12::new(self.c0 + rhs.c0, self.c1 + rhs.c1)
    }
}

impl<'a, F: Fp12Element> Neg for &'a Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp12::new(-self.c0, -self.c1)
    }
}

impl<F: Fp12Element> Neg for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F: Fp12Element> Sub<Fp12<F>> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn sub(self, rhs: Fp12<F>) -> Self::Output {
        Fp12::new(self.c0 - rhs.c0, self.c1 - rhs.c1)
    }
}

impl<F: Fp12Element> Mul<F> for Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn mul(self, rhs: F) -> Fp12<F> {
        let rhs = Fp12::from(rhs);
        Fp12::new(self.c0 * rhs.c0, self.c1 * rhs.c1)
    }
}

impl<'a, 'b, F: Fp12Element> Div<&'b Fp12<F>> for &'a Fp12<F> {
    type Output = Fp12<F>;

    #[inline]
    fn div(self, rhs: &'b Fp12<F>) -> Fp12<F> {
        self.div(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod bls12381_fp12_test {
        use super::*;

        fn bls12381_fp12_rand() -> Fp12<Bls12381> {
            let mut rng = rand::thread_rng();
            Fp12::new(
                Fp6::new(
                    Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                    Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                    Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                ),
                Fp6::new(
                    Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                    Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                    Fp2::new(Bls12381::random(&mut rng), Bls12381::random(&mut rng)),
                ),
            )
        }

        #[test]
        fn test_equality() {
            let rng = &mut rand::thread_rng();
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                let b = a;
                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_inequality() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                let b = bls12381_fp12_rand();
                if a != b {
                    assert_ne!(a, b);
                }
            }
        }

        #[test]
        fn test_addition_subtraction() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                let b = bls12381_fp12_rand();
                let c = bls12381_fp12_rand();

                // commutative
                assert_eq!(a + b, b + a);
                assert_eq!(a + (b + c), (a + b) + c);

                // additive identity
                assert_eq!(a + Fp12::<Bls12381>::zero(), a);
                assert_eq!(a - Fp12::<Bls12381>::zero(), a);

                assert_eq!(Fp12::<Bls12381>::zero() - a, -a);
                assert_eq!(a - b, a + (-b));
                assert_eq!(a - b, a + (b * -Fp12::<Bls12381>::one()));

                assert_eq!(-a, Fp12::<Bls12381>::zero() - a);
                assert_eq!(-a, a * -Fp12::<Bls12381>::one());
            }
        }

        #[test]
        fn test_multiplication() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                let b = bls12381_fp12_rand();
                let c = bls12381_fp12_rand();

                // commutative
                assert_eq!(a * b, b * a);

                // associative
                assert_eq!(a * (b * c), (a * b) * c);

                // distributive
                assert_eq!(a * (b + c), a * b + a * c);
            }
        }

        #[test]
        fn test_square_equality() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                assert_eq!(a.square(), a * a);
            }
        }

        #[test]
        fn test_inversion() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                if !a.is_zero() {
                    assert_eq!(a * a.invert().unwrap(), Fp12::<Bls12381>::one());
                    assert_eq!(a.invert().unwrap().invert().unwrap(), a);
                }
            }
        }

        #[test]
        fn test_frobenius() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                {
                    let b = (0..12).fold(a, |acc, _| acc.frobenius_map());
                    assert_eq!(a, b);
                }
            }
        }

        #[test]
        fn test_cyclotomic_square() {
            for _ in 0..10 {
                let a = bls12381_fp12_rand();
                assert_eq!(a.cyclotomic_square(), a.n_cyclotomic_square(1));
                assert_eq!(
                    a.cyclotomic_square().cyclotomic_square(),
                    a.n_cyclotomic_square(2)
                );
            }
        }
    }

    mod bn254_fp12_test {
        use super::*;

        fn bn254_fp12_rand() -> Fp12<Bn254> {
            let mut rng = rand::thread_rng();
            Fp12::new(
                Fp6::new(
                    Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                    Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                    Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                ),
                Fp6::new(
                    Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                    Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                    Fp2::new(Bn254::random(&mut rng), Bn254::random(&mut rng)),
                ),
            )
        }

        #[test]
        fn test_equality() {
            let rng = &mut rand::thread_rng();
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                let b = a;
                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_inequality() {
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                let b = bn254_fp12_rand();
                if a != b {
                    assert_ne!(a, b);
                }
            }
        }

        #[test]
        fn test_addition_subtraction() {
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                let b = bn254_fp12_rand();
                let c = bn254_fp12_rand();

                // commutative
                assert_eq!(a + b, b + a);
                assert_eq!(a + (b + c), (a + b) + c);

                // additive identity
                assert_eq!(a + Fp12::<Bn254>::zero(), a);
                assert_eq!(a - Fp12::<Bn254>::zero(), a);

                assert_eq!(Fp12::<Bn254>::zero() - a, -a);
                assert_eq!(a - b, a + (-b));
                assert_eq!(a - b, a + (b * -Fp12::<Bn254>::one()));

                assert_eq!(-a, Fp12::<Bn254>::zero() - a);
                assert_eq!(-a, a * -Fp12::<Bn254>::one());
            }
        }

        #[test]
        fn test_multiplication() {
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                let b = bn254_fp12_rand();
                let c = bn254_fp12_rand();

                // commutative
                assert_eq!(a * b, b * a);

                // associative
                assert_eq!(a * (b * c), (a * b) * c);

                // distributive
                assert_eq!(a * (b + c), a * b + a * c);
            }
        }

        #[test]
        fn test_square_equality() {
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                assert_eq!(a.square(), a * a);
            }
        }

        #[test]
        fn test_inversion() {
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                if !a.is_zero() {
                    assert_eq!(a * a.invert().unwrap(), Fp12::<Bn254>::one());
                    assert_eq!(a.invert().unwrap().invert().unwrap(), a);
                }
            }
        }

        #[test]
        fn test_frobenius() {
            {
                for _ in 0..10 {
                    let a = bn254_fp12_rand();
                    let lhs = a.frobenius_map();
                    let rhs = (0..6).fold(a, |acc, _| acc.frobenius_map());
                    assert_eq!(lhs, rhs);
                }
            }
            // {
            //     for _ in 0..10 {
            //         let a = bn254_fp12_rand();
            //         // {
            //         // let b = (0..12).fold(a, |acc, _| acc.frobenius_map());
            //         // assert_eq!(a, b);
            //         // }
            //         let mut b = a;
            //         for _ in 0..12 {
            //             b = b.frobenius_map();
            //             println!("b: {:?}", b);
            //         }
            //         assert_eq!(a, b);
            //     }
            // }
        }

        #[test]
        fn test_cyclotomic_square() {
            for _ in 0..10 {
                let a = bn254_fp12_rand();
                assert_eq!(a.cyclotomic_square(), a.n_cyclotomic_square(1));
                assert_eq!(
                    a.cyclotomic_square().cyclotomic_square(),
                    a.n_cyclotomic_square(2)
                );
            }
        }

        #[test]
        fn test_bn254_frobenius() {
            for _ in 0..12 {
                let lhs = bn254_fp12_rand();
                let rhs = (0..12).fold(lhs, |acc, _| acc.frobenius_map());
                assert_eq!(lhs, rhs);
            }
        }

        #[test]
        fn test_mul_ne() {
            let a1 = bn254_fp12_rand();
            let a2 = bn254_fp12_rand();

            let b1 = bn254_fp12_rand();
            let b2 = bn254_fp12_rand();

            assert_ne!(a1 * a2, b1 * b2);
        }
    }
}
