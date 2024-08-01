use field::{common::AffinePoint, fp::Fp, fp12::Fp12, fp2::Fp2, fp6::Fp6, fr::Fr, utils};
use std::f32::NEG_INFINITY;

use crate::{
    common::Curve,
    g1::G1Affine,
    g2::{G2Affine, G2Projective},
};

type TwistedFieldElement<C: Curve> = Option<(Fp12<C>, Fp12<C>)>;

fn twist<C: Curve>(p: &G2Affine<C>) -> TwistedFieldElement<C> {
    let x = p.x;
    let y = p.y;

    let _x = (x.c0 - x.c1, x.c1);
    let _y = (y.c0 - y.c1, y.c1);

    let nx = Fp12::new(
        Fp6 {
            c0: Fp2::new(_x.0, Fp::zero()),
            c1: Fp2::new(Fp::zero(), Fp::zero()),
            c2: Fp2::new(Fp::zero(), Fp::zero()),
        },
        Fp6 {
            c0: Fp2::new(_x.1, Fp::zero()),
            c1: Fp2::new(Fp::zero(), Fp::zero()),
            c2: Fp2::new(Fp::zero(), Fp::zero()),
        },
    );

    let ny = Fp12::new(
        Fp6 {
            c0: Fp2::new(_y.0, Fp::zero()),
            c1: Fp2::new(Fp::zero(), Fp::zero()),
            c2: Fp2::new(Fp::zero(), Fp::zero()),
        },
        Fp6 {
            c0: Fp2::new(_y.1, Fp::zero()),
            c1: Fp2::new(Fp::zero(), Fp::zero()),
            c2: Fp2::new(Fp::zero(), Fp::zero()),
        },
    );

    let w: Fp12<C> = Fp12::new(
        Fp6 {
            c0: Fp2::new(Fp::zero(), Fp::one()),
            c1: Fp2::new(Fp::zero(), Fp::zero()),
            c2: Fp2::new(Fp::zero(), Fp::zero()),
        },
        Fp6 {
            c0: Fp2::new(Fp::zero(), Fp::zero()),
            c1: Fp2::new(Fp::zero(), Fp::zero()),
            c2: Fp2::new(Fp::zero(), Fp::zero()),
        },
    );

    Some((nx / w.square(), ny / (w * w.square())))
}

fn embed_g1_in_fp12<C: Curve>(a: &G1Affine<C>) -> Option<(Fp12<C>, Fp12<C>)> {
    Some((Fp12::from(a.x), Fp12::from(a.y)))
}

fn ell<C: Curve>(
    p1: &TwistedFieldElement<C>,
    p2: &TwistedFieldElement<C>,
    t: &TwistedFieldElement<C>,
) -> Fp12<C> {
    let (x1, y1) = p1.unwrap();
    let (x2, y2) = p2.unwrap();
    let (xt, yt) = t.unwrap();

    if x1 != x2 {
        let m = (y2 - y2) / (x2 - x1);
        m * (xt - x1) - yt + y1
    } else if y1 == y2 {
        let m = (x1.square() * Fp::from(3)) / (y1 * Fp::from(2));
        m * (xt - x1) - yt + y1
    } else {
        xt - x1
    }
}

fn double_step<C: Curve>(p: &TwistedFieldElement<C>) -> TwistedFieldElement<C> {
    if p.is_none() {
        return *p;
    }

    let (x, y) = p.unwrap();
    let m = (x.square() * Fp::from(3)) / (y * Fp::from(2));
    let x3 = m.square() - (x * Fp::from(2));
    let y3 = m * (x - x3) - y;
    Some((x3, y3))
}

fn add_step<C: Curve>(
    p: &TwistedFieldElement<C>,
    q: &TwistedFieldElement<C>,
) -> TwistedFieldElement<C> {
    if p.is_none() {
        return q.clone();
    } else if q.is_none() {
        return p.clone();
    }

    let (x1, y1) = p.unwrap();
    let (x2, y2) = q.unwrap();

    if (x1 == x2) && (y1 == y2) {
        return double_step(p);
    } else if x2 == x1 {
        return None;
    }

    let m = (y2 - y1) / (x2 - x1);
    let x3 = m.square() - x1 - x2;
    let y3 = m * (x1 - x3) - y1;
    Some((x3, y3))
}

fn miller_loop<C: Curve>(p: &G1Affine<C>, q: &G2Affine<C>) -> Fp12<C> {
    if p.is_identity() || q.is_identity() {
        return Fp12::one();
    }
    let p = embed_g1_in_fp12(p);
    let q = twist(q);
    let mut f = Fp12::<C>::one();
    let mut r = q.clone();

    C::LOOP_COUNTER.iter().rev().for_each(|&i| {
        println!("cycle-tracker-start: square * ell");
        f = f.square() * ell(&r, &r, &p);
        println!("cycle-tracker-end: square * ell");

        println!("cycle-tracker-start: double");
        r = double_step(&r);
        println!("cycle-tracker-end: double");
        // if C::LOG_ATE_LOOP_COUNT & (1 << i) != 0 {
        if i > 0 {
            f = f * ell(&r, &q, &p);
            println!("cycle-tracker-start: add");
            r = add_step(&r, &q);
            println!("cycle-tracker-end: add");
        }
    });
    // let mut r = G2Projective::<C>::from_affine(q);
    // let mut f = Fp12::<C>::one();

    // let mut found_one = false;
    // (0..64)
    //     .rev()
    //     .map(|b| (((C::X >> 1) >> b) & 1) == 1)
    //     .for_each(|i| {
    //         if !found_one {
    //             found_one = i;
    //             return;
    //         }

    //         let coeffs = doubling_step(&mut r);
    //         f = ell(f, &coeffs, &p);

    //         if i {
    //             let coeffs = addition_step(&mut r, &q);
    //             f = ell(f, &coeffs, &p);
    //         }

    //         f = f.square();
    //     });
    // let coeffs = doubling_step(&mut r);
    // f = ell(f, &coeffs, &p);

    // f.conjugate()
    f
}

// pub struct G2Prepared<C: Curve> {
//     coeffs: Vec<(Fp2<C>, Fp2<C>, Fp2<C>)>,
//     is_infinity: bool,
// }

// impl<C: Curve> From<G2Affine<C>> for G2Prepared<C> {
//     fn from(q: G2Affine<C>) -> Self {
//         let is_identity = q.is_identity();
//         let q = is_identity
//             .then(G2Affine::<C>::generator)
//             .unwrap_or_else(|| q);

//         let mut coeffs = Vec::with_capacity(68);
//         let mut r = G2Projective::<C>::from_affine(q);
//         let mut f = Fp12::<C>::one();

//         let mut found_one = false;
//         (0..64)
//             .rev()
//             .map(|b| (((C::X >> 1) >> b) & 1) == 1)
//             .for_each(|i| {
//                 if !found_one {
//                     found_one = i;
//                     return;
//                 }

//                 coeffs.push(doubling_step(&mut r));

//                 if i {
//                     coeffs.push(addition_step(&mut r, &q));
//                 }

//                 f = f.square();
//             });
//         coeffs.push(doubling_step(&mut r));

//         G2Prepared {
//             coeffs,
//             is_infinity: is_identity,
//         }
//     }
// }

// fn multi_miller_loop<C: Curve>(p: &[G1Affine<C>], q: &[G2Prepared<C>]) -> Fp12<C> {
//     let mut f = Fp12::<C>::one();
//     let mut found_one = false;
//     let mut j = 0;

//     (0..64)
//         .rev()
//         .map(|b| (((C::X >> 1) >> b) & 1) == 1)
//         .for_each(|i| {
//             if !found_one {
//                 found_one = i;
//                 return;
//             }

//             p.iter().zip(q.iter()).for_each(|(a, b)| {
//                 (!(a.is_identity() || b.is_infinity)).then(|| {
//                     f = ell(f, &b.coeffs[j], a);
//                     // println!("f = {:?}", f);
//                 });
//             });
//             j += 1;
//             // println!("j = {}", j);

//             if i {
//                 p.iter().zip(q.iter()).for_each(|(a, b)| {
//                     (!(a.is_identity() || b.is_infinity)).then(|| {
//                         f = ell(f, &b.coeffs[j], a);
//                     });
//                 });
//                 j += 1;
//                 // println!("j (add) = {}", j);
//             }

//             f = f.square();
//         });

//     f.conjugate()
// }

fn residue_test<C: Curve>(x: &Fp12<C>) -> bool {
    let (x, scaling_factor) = C::get_root_and_scaling_factor(x);

    let t0 = x.frobenius_map();
    let t1 = x.powt();

    let lhs = t0 * t1;
    let rhs = x * scaling_factor;

    lhs == rhs
}

pub fn verify_pairing<C: Curve>(p: &[G1Affine<C>], q: &[G2Affine<C>]) -> bool {
    // let q = q
    //     .iter()
    //     .map(|q| G2Prepared::from(*q))
    //     .collect::<Vec<G2Prepared<C>>>();

    // let f = multi_miller_loop(p, &q);
    p.iter().zip(q.iter()).for_each(|(p, q)| {
        let f = miller_loop(p, q);
    });
    let f = Fp12::<C>::one();
    println!("f = {:?}", f);

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
        // let p1 = G1Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        // let q1 = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        // let p2 = G1Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        // let q2 = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        println!("p1: {:?}", p1);
        println!("q1: {:?}", q1);
        println!("millier loop: {:?}", miller_loop(&p1, &q1));
        // println!("verify pairing: {:?}", verify_pairing(&[p1, p2], &[q1, q2]));
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
