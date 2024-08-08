use crate::fp12::Fp12Element;
use crate::fp2::Fp2Element;
use crate::g1::G1Element;
use crate::g2::G2Element;
use crate::{common::AffinePoint, fp12::Fp12, fp2::Fp2, fp6::Fp6};

use crate::{
    g1::G1Affine,
    g2::{G2Affine, G2Projective},
};

fn ell<F: Fp12Element + G1Element>(
    f: Fp12<F>,
    coeffs: &(Fp2<F>, Fp2<F>, Fp2<F>),
    p: &G1Affine<F>,
) -> Fp12<F> {
    let mut c0 = coeffs.0;
    let mut c1 = coeffs.1;

    c0.c0 *= p.y;
    c0.c1 *= p.y;

    c1.c0 *= p.x;
    c1.c1 *= p.x;

    f.mul_by_014(&coeffs.2, &c1, &c0)
}

// Adaptation of Algorithm 26, https://eprint.iacr.org/2010/354.pdf
fn doubling_step<F: Fp2Element + G2Element>(r: &mut G2Projective<F>) -> (Fp2<F>, Fp2<F>, Fp2<F>) {
    let tmp0 = r.x.square();
    let tmp1 = r.y.square();
    let tmp2 = tmp1.square();
    let tmp3 = (tmp1 + r.x).square() - tmp0 - tmp2;
    let tmp3 = tmp3 + tmp3;
    let tmp4 = tmp0 + tmp0 + tmp0;
    let tmp6 = r.x + tmp4;
    let tmp5 = tmp4.square();
    let zsquared = r.z.square();
    r.x = tmp5 - tmp3 - tmp3;
    r.z = (r.z + r.y).square() - tmp1 - zsquared;
    r.y = (tmp3 - r.x) * tmp4;
    let tmp2 = tmp2 + tmp2;
    let tmp2 = tmp2 + tmp2;
    let tmp2 = tmp2 + tmp2;
    r.y = r.y - tmp2;
    let tmp3 = tmp4 * zsquared;
    let tmp3 = tmp3 + tmp3;
    let tmp3 = -tmp3;
    let tmp6 = tmp6.square() - tmp0 - tmp5;
    let tmp1 = tmp1 + tmp1;
    let tmp1 = tmp1 + tmp1;
    let tmp6 = tmp6 - tmp1;
    let tmp0 = r.z * zsquared;
    let tmp0 = tmp0 + tmp0;

    (tmp0, tmp3, tmp6)
}

// Adaptation of Algorithm 27, https://eprint.iacr.org/2010./354.pdf
fn addition_step<F: Fp2Element + G2Element>(
    r: &mut G2Projective<F>,
    q: &G2Affine<F>,
) -> (Fp2<F>, Fp2<F>, Fp2<F>) {
    let zsquared = r.z.square();
    let ysquared = q.y.square();
    let t0 = zsquared * q.x;
    let t1 = ((q.y + r.z).square() - ysquared - zsquared) * zsquared;
    let t2 = t0 - r.x;
    let t3 = t2.square();
    let t4 = t3 + t3;
    let t4 = t4 + t4;
    let t5 = t4 * t2;
    let t6 = t1 - r.y - r.y;
    let t9 = t6 * q.x;
    let t7 = t4 * r.x;
    r.x = t6.square() - t5 - t7 - t7;
    r.z = (r.z + t2).square() - zsquared - t3;
    let t10 = q.y + r.z;
    let t8 = (t7 - r.x) * t6;
    let t0 = r.y * t5;
    let t0 = t0 + t0;
    r.y = t8 - t0;
    let t10 = t10.square() - ysquared;
    let ztsquared = r.z.square();
    let t10 = t10 - ztsquared;
    let t9 = t9 + t9 - t10;
    let t10 = r.z + r.z;
    let t6 = -t6;
    let t1 = t6 + t6;

    (t10, t1, t9)
}

fn miller_loop<F: Fp12Element + G1Element + G2Element>(
    p: &G1Affine<F>,
    q: &G2Affine<F>,
) -> Fp12<F> {
    let p = p
        .is_identity()
        .then(G1Affine::<F>::generator)
        .unwrap_or_else(|| *p);
    let q = q
        .is_identity()
        .then(<F as G2Element>::generator)
        .unwrap_or_else(|| *q);

    let mut r = G2Projective::<F>::from_affine(q);
    let mut f = Fp12::<F>::one();

    let mut found_one = false;
    (0..64)
        .rev()
        .map(|b| (((F::X >> 1) >> b) & 1) == 1)
        .for_each(|i| {
            if !found_one {
                found_one = i;
                return;
            }

            let coeffs = doubling_step(&mut r);
            f = ell(f, &coeffs, &p);

            if i {
                let coeffs = addition_step(&mut r, &q);
                f = ell(f, &coeffs, &p);
            }

            f = f.square();
        });
    let coeffs = doubling_step(&mut r);
    f = ell(f, &coeffs, &p);

    f.conjugate()
}

pub struct G2Prepared<F: Fp2Element> {
    coeffs: Vec<(Fp2<F>, Fp2<F>, Fp2<F>)>,
    is_infinity: bool,
}

impl<F: Fp12Element + G2Element> From<G2Affine<F>> for G2Prepared<F> {
    fn from(q: G2Affine<F>) -> Self {
        let is_identity = q.is_identity();
        let q = is_identity
            .then(G2Affine::<F>::generator)
            .unwrap_or_else(|| q);

        let mut coeffs = Vec::with_capacity(68);
        let mut r = G2Projective::<F>::from_affine(q);
        let mut f = Fp12::<F>::one();

        let mut found_one = false;
        for i in (0..64).rev().map(|b| (((F::X >> 1) >> b) & 1) == 1) {
            if !found_one {
                found_one = i;
                continue;
            }

            coeffs.push(doubling_step(&mut r));

            if i {
                coeffs.push(addition_step(&mut r, &q));
            }

            f = f.square();
        }
        coeffs.push(doubling_step(&mut r));

        assert_eq!(coeffs.len(), 68);

        G2Prepared {
            coeffs,
            is_infinity: is_identity,
        }
    }
}

fn multi_miller_loop<F: Fp12Element + G1Element>(
    p: &[G1Affine<F>],
    q: &[G2Prepared<F>],
) -> Fp12<F> {
    let mut f = Fp12::<F>::one();
    let mut found_one = false;
    let mut j = 0;

    for i in (0..64).rev().map(|b| (((F::X >> 1) >> b) & 1) == 1) {
        if !found_one {
            found_one = i;
            continue;
        }
        p.iter().zip(q.iter()).for_each(|(a, b)| {
            if !(a.is_identity() || b.is_infinity) {
                f = ell(f, &b.coeffs[j], a);
            }
        });
        j += 1;

        if i {
            p.iter().zip(q.iter()).for_each(|(a, b)| {
                (!(a.is_identity() || b.is_infinity)).then(|| {
                    f = ell(f, &b.coeffs[j], a);
                });
            });
            j += 1;
        }
        f = f.square();
    }

    p.iter().zip(q.iter()).for_each(|(a, b)| {
        (!(a.is_identity() || b.is_infinity)).then(|| {
            f = ell(f, &b.coeffs[j], a);
        });
    });

    f.conjugate()
}

// https://eprint.iacr.org/2009/565.pdf
pub fn final_exponentiation<F: Fp12Element + AffinePoint>(&f: &Fp12<F>) -> Fp12<F> {
    #[must_use]
    fn fp4_square<F: Fp12Element>(a: Fp2<F>, b: Fp2<F>) -> (Fp2<F>, Fp2<F>) {
        let t0 = a.square();
        let t1 = b.square();
        let mut t2 = t1.mul_by_nonresidue();
        let c0 = t2 + t0;
        t2 = a + b;
        t2 = t2.square();
        t2 = t2 - t0;
        let c1 = t2 - t1;

        (c0, c1)
    }
    #[must_use]
    fn cyclotomic_square<F: Fp12Element>(f: Fp12<F>) -> Fp12<F> {
        let mut z0 = f.c0.c0;
        let mut z4 = f.c0.c1;
        let mut z3 = f.c0.c2;
        let mut z2 = f.c1.c0;
        let mut z1 = f.c1.c1;
        let mut z5 = f.c1.c2;

        let (t0, t1) = fp4_square(z0, z1);

        // For A
        z0 = t0 - z0;
        z0 = z0 + z0 + t0;

        z1 = t1 + z1;
        z1 = z1 + z1 + t1;

        let (mut t0, t1) = fp4_square(z2, z3);
        let (t2, t3) = fp4_square(z4, z5);

        // For C
        z4 = t0 - z4;
        z4 = z4 + z4 + t0;

        z5 = t1 + z5;
        z5 = z5 + z5 + t1;

        // For B
        t0 = t3.mul_by_nonresidue();
        z2 = t0 + z2;
        z2 = z2 + z2 + t0;

        z3 = t2 - z3;
        z3 = z3 + z3 + t2;

        Fp12::new(Fp6::new(z0, z4, z3), Fp6::new(z2, z1, z5))
    }

    fn cycolotomic_exp<F: Fp12Element + AffinePoint>(f: Fp12<F>) -> Fp12<F> {
        let x = F::X;
        let mut tmp = Fp12::<F>::one();
        let mut found_one = false;
        for i in (0..64).rev().map(|b| ((x >> b) & 1) == 1) {
            if found_one {
                tmp = cyclotomic_square(tmp)
            } else {
                found_one = i;
            }

            if i {
                tmp = tmp * f;
            }
        }

        tmp.conjugate()
    }

    let mut f = f.clone();
    let mut t0 = f
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map();
    // let mut t0 = fh_frobenius_map(6);

    f.invert()
        .map(|mut t1| {
            let mut t2 = t0 * t1;
            t1 = t2;
            t2 = t2.frobenius_map().frobenius_map();
            // t2 = t2.nth_frobenius_map(2);
            t2 = t2 * t1;
            t1 = cyclotomic_square(t2).conjugate();
            let mut t3 = cycolotomic_exp(t2);
            let mut t4 = cyclotomic_square(t3);
            let mut t5 = t1 * t3;
            t1 = cycolotomic_exp(t5);
            t0 = cycolotomic_exp(t1);
            let mut t6 = cycolotomic_exp(t0);
            t6 = t6 * t4;
            t4 = cycolotomic_exp(t6);
            t5 = t5.conjugate();
            t4 = t4 * t5 * t2;
            t5 = t2.conjugate();
            t1 = t1 * t2;
            t1 = t1.frobenius_map().frobenius_map().frobenius_map();
            // t1 = t1.nth_frobenius_map(3);
            t6 = t6 * t5;
            t6 = t6.frobenius_map();
            t3 = t3 * t0;
            t3 = t3.frobenius_map().frobenius_map();
            // t3 = t3.nth_frobenius_map(2);
            t3 = t3 * t1;
            t3 = t3 * t6;
            f = t3 * t4;

            f
        })
        .unwrap()
}

pub fn verify_pairing<F: Fp12Element + G1Element + G2Element>(
    p: &[G1Affine<F>],
    q: &[G2Affine<F>],
) -> bool {
    let q = q
        .iter()
        .map(|q| G2Prepared::from(*q))
        .collect::<Vec<G2Prepared<F>>>();
    let f = multi_miller_loop(p, &q);
    final_exponentiation(&f) == Fp12::<F>::one()
}
