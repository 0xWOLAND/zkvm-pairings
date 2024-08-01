use field::{common::AffinePoint, fp::Fp, fp12::Fp12, fp2::Fp2, fp6::Fp6, fr::Fr, utils};
use std::f32::NEG_INFINITY;

use crate::{
    common::Curve,
    g1::G1Affine,
    g2::{G2Affine, G2Projective},
};

/// 6U+2 for in NAF form

pub(crate) const SIX_U_PLUS_2_NAF: [i8; 65] = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0,
    0, 1, 0, 1, 1,
];

#[derive(Clone, Copy)]
pub(crate) struct LineEvaluation<C: Curve> {
    pub(crate) x: Fp2<C>,
    pub(crate) y: Fp2<C>,
}

impl<C: Curve> LineEvaluation<C> {
    pub(crate) fn zero() -> Self {
        LineEvaluation {
            x: Fp2::<C>::zero(),
            y: Fp2::<C>::zero(),
        }
    }
}

fn ell<C: Curve>(f: Fp12<C>, coeffs: &(Fp2<C>, Fp2<C>, Fp2<C>), p: &G1Affine<C>) -> Fp12<C> {
    let mut c0 = coeffs.0;
    let mut c1 = coeffs.1;

    c0.c0 *= p.y;
    c0.c1 *= p.y;

    c1.c0 *= p.x;
    c1.c1 *= p.x;

    f.mul_by_014(&coeffs.2, &c1, &c0)
}

fn doubling_step<C: Curve>(r: &mut G2Projective<C>) -> (Fp2<C>, Fp2<C>, Fp2<C>) {
    // Adaptation of Algorithm 26, https://eprint.iacr.org/2010/354.pdf
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
    r.y -= tmp2;
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

fn addition_step<C: Curve>(r: &mut G2Projective<C>, q: &G2Affine<C>) -> (Fp2<C>, Fp2<C>, Fp2<C>) {
    // Adaptation of Algorithm 27, https://eprint.iacr.org/2010./354.pdf
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

fn miller_loop<C: Curve>(p: &G1Affine<C>, q: &G2Affine<C>) -> Fp12<C> {
    let p = p
        .is_identity()
        .then(G1Affine::<C>::generator)
        .unwrap_or_else(|| *p);
    let q = q
        .is_identity()
        .then(G2Affine::<C>::generator)
        .unwrap_or_else(|| *q);

    let mut r = G2Projective::<C>::from_affine(q);
    let mut f = Fp12::<C>::one();

    let mut found_one = false;
    (0..64)
        .rev()
        .map(|b| (((C::X >> 1) >> b) & 1) == 1)
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

pub struct G2Prepared<C: Curve> {
    coeffs: Vec<(Fp2<C>, Fp2<C>, Fp2<C>)>,
    is_infinity: bool,
}

impl<C: Curve> From<G2Affine<C>> for G2Prepared<C> {
    fn from(q: G2Affine<C>) -> Self {
        let is_identity = q.is_identity();
        let q = is_identity
            .then(G2Affine::<C>::generator)
            .unwrap_or_else(|| q);

        let mut coeffs = Vec::with_capacity(68);
        let mut r = G2Projective::<C>::from_affine(q);
        let mut f = Fp12::<C>::one();

        let mut found_one = false;
        let mut j = 0;
        println!("X: {}", C::X);
        for i in (0..64).rev().map(|b| (((C::X >> 1) >> b) & 1) == 1) {
            if !found_one {
                println!("found one at j = {}", j);
                found_one = i;
                continue;
            }

            coeffs.push(doubling_step(&mut r));

            if i {
                coeffs.push(addition_step(&mut r, &q));
            }

            f = f.square();
            j += 1;
        }
        coeffs.push(doubling_step(&mut r));

        assert_eq!(coeffs.len(), 68);

        G2Prepared {
            coeffs,
            is_infinity: is_identity,
        }
    }
}

fn multi_miller_loop<C: Curve>(p: &[G1Affine<C>], q: &[G2Prepared<C>]) -> Fp12<C> {
    let mut f = Fp12::<C>::one();
    let mut found_one = false;
    let mut j = 0;

    for i in (0..64).rev().map(|b| (((C::X >> 1) >> b) & 1) == 1) {
        if !found_one {
            found_one = i;
            continue;
        }
        p.iter().zip(q.iter()).for_each(|(a, b)| {
            (!(a.is_identity() || b.is_infinity)).then(|| {
                f = ell(f, &b.coeffs[j], a);
            });
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

fn residue_test<C: Curve>(x: &Fp12<C>) -> bool {
    let (x, scaling_factor) = C::get_root_and_scaling_factor(x);

    let t0 = x.frobenius_map();
    let t1 = x.powt();

    let lhs = t0 * t1;
    let rhs = x * scaling_factor;

    lhs == rhs
}

pub fn verify_pairing<C: Curve>(p: &[G1Affine<C>], q: &[G2Affine<C>]) -> bool {
    let q = q
        .iter()
        .map(|q| G2Prepared::from(*q))
        .collect::<Vec<G2Prepared<C>>>();
    let f = multi_miller_loop(p, &q);
    let buf = f.conjugate() / f;
    let f = buf.frobenius_map() * buf;
    residue_test(&f)
}

#[cfg(test)]
mod test {
    use field::common::Bls12381Curve;

    use super::*;

    #[test]
    fn test_miller_loop() {
        let p1 = G1Affine::<Bls12381Curve>::generator();
        let q1 = G2Affine::<Bls12381Curve>::generator();
        let p2 = G1Affine::<Bls12381Curve>::generator();
        let q2 = G2Affine::<Bls12381Curve>::generator();

        println!(
            "{:?}",
            multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)])
        );
    }
    // #[test]
    // fn test_dobule() {
    //     let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

    //     let (p1, l1) = double_step(p);
    //     let p2 = p.double();

    //     assert_eq!(p1, p2);
    // }

    // #[test]
    // fn test_add() {
    //     let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
    //     let q = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

    //     let (p1, l2) = add_step(p, q);
    //     let p2 = p + q;

    //     assert_eq!(p1, p2);
    // }

    // #[test]
    // fn test_double_and_add() {
    //     let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
    //     let q = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

    //     let (p1, l1, ll1) = double_and_add_step(&p, &q);
    //     let (_p2, l2) = double_step(p);
    //     let (p2, ll2) = add_step(_p2, q);

    //     assert_eq!(p1, p2);
    // }

    // #[test]
    // fn test_triple_step() {
    //     let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

    //     let (p1, l1, l2) = triple_step(&p);
    //     // let (_p2, ll1) = double_step(p);
    //     // let (p2, ll2) = add_step(_p2, p);
    //     let p2 = p + p + p;

    //     assert_eq!(p1, p2);
    // }
}
