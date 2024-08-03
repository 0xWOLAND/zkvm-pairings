use crate::{common::AffinePoint, fp::Fp, fp12::Fp12, fp2::Fp2, fp6::Fp6};

use crate::{
    common::Curve,
    g1::G1Affine,
    g2::{G2Affine, G2Projective},
};

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

#[derive(Debug)]
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

        println!("cycle-tracker-start: miller-loop");

        let mut found_one = false;
        for i in (0..64).rev().map(|b| (((C::X >> 1) >> b) & 1) == 1) {
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

        println!("cycle-tracker-start: miller-loop");
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
pub fn final_exponentiation<C: Curve>(&f: &Fp12<C>) -> Fp12<C> {
    println!("FINAL EXPONENTIATION");
    #[must_use]
    fn fp4_square<C: Curve>(a: Fp2<C>, b: Fp2<C>) -> (Fp2<C>, Fp2<C>) {
        let t0 = a.square();
        let t1 = b.square();
        let mut t2 = t1.mul_by_nonresidue();
        let c0 = t2 + t0;
        t2 = a + b;
        t2 = t2.square();
        t2 -= t0;
        let c1 = t2 - t1;

        (c0, c1)
    }
    #[must_use]
    fn cyclotomic_square<C: Curve>(f: Fp12<C>) -> Fp12<C> {
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
    #[must_use]
    fn cycolotomic_exp<C: Curve>(f: Fp12<C>) -> Fp12<C> {
        let x = C::X;
        let mut tmp = Fp12::<C>::one();
        let mut found_one = false;
        for i in (0..64).rev().map(|b| ((x >> b) & 1) == 1) {
            if found_one {
                tmp = cyclotomic_square(tmp)
            } else {
                found_one = i;
            }

            if i {
                tmp *= f;
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

    f.invert()
        .map(|mut t1| {
            let mut t2 = t0 * t1;
            t1 = t2;
            t2 = t2.frobenius_map().frobenius_map();
            t2 *= t1;
            t1 = cyclotomic_square(t2).conjugate();
            let mut t3 = cycolotomic_exp(t2);
            let mut t4 = cyclotomic_square(t3);
            let mut t5 = t1 * t3;
            t1 = cycolotomic_exp(t5);
            t0 = cycolotomic_exp(t1);
            let mut t6 = cycolotomic_exp(t0);
            t6 *= t4;
            t4 = cycolotomic_exp(t6);
            t5 = t5.conjugate();
            t4 *= t5 * t2;
            t5 = t2.conjugate();
            t1 *= t2;
            t1 = t1.frobenius_map().frobenius_map().frobenius_map();
            t6 *= t5;
            t6 = t6.frobenius_map();
            t3 *= t0;
            t3 = t3.frobenius_map().frobenius_map();
            t3 *= t1;
            t3 *= t6;
            f = t3 * t4;

            f
        })
        .unwrap()
}

pub fn verify_pairing<C: Curve>(p: &[G1Affine<C>], q: &[G2Affine<C>]) -> bool {
    let q = q
        .iter()
        .map(|q| G2Prepared::from(*q))
        .collect::<Vec<G2Prepared<C>>>();
    let f = multi_miller_loop(p, &q);
    final_exponentiation(&f) == Fp12::<C>::one()
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use crate::{common::Bls12381Curve, fr::Fr};
    use num_bigint::BigUint;
    use rand::thread_rng;

    use super::*;

    fn _final_exponentiation(x: &Fp12<Bls12381Curve>) -> Fp12<Bls12381Curve> {
        let h: &BigUint= &BigUint::from_str("322277361516934140462891564586510139908379969514828494218366688025288661041104682794998680497580008899973249814104447692778988208376779573819485263026159588510513834876303014016798809919343532899164848730280942609956670917565618115867287399623286813270357901731510188149934363360381614501334086825442271920079363289954510565375378443704372994881406797882676971082200626541916413184642520269678897559532260949334760604962086348898118982248842634379637598665468817769075878555493752214492790122785850202957575200176084204422751485957336465472324810982833638490904279282696134323072515220044451592646885410572234451732790590013479358343841220074174848221722017083597872017638514103174122784843925578370430843522959600095676285723737049438346544753168912974976791528535276317256904336520179281145394686565050419250614107803233314658825463117900250701199181529205942363159325765991819433914303908860460720581408201373164047773794825411011922305820065611121544561808414055302212057471395719432072209245600258134364584636810093520285711072578721435517884103526483832733289802426157301542744476740008494780363354305116978805620671467071400711358839553375340724899735460480144599782014906586543813292157922220645089192130209334926661588737007768565838519456601560804957985667880395221049249803753582637708560").unwrap();
        x.pow_vartime_extended(&h.to_u64_digits())
    }

    #[test]
    fn test_final_exponentiation() {
        for _ in 0..10 {
            let f = Fp12::<Bls12381Curve>::random(&mut thread_rng());
            // let challenge = Fp::<Bls12381Curve>::random(&mut thread_rng());
            // let f = Fp12::<Bls12381Curve>::one() * challenge;
            let lhs = final_exponentiation(&f);
            let rhs = _final_exponentiation(&f);
            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn test_random_points() {
        for _ in 0..1 {
            {
                let p1 = &G1Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());
                let q1 = &G2Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());
                let p2 = &G1Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());
                let q2 = -G2Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());

                let eq_class =
                    multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)]);
                let rhs = final_exponentiation(&eq_class) == Fp12::<Bls12381Curve>::one();
                let lhs = verify_pairing(&[p1, p2], &[q1, q2]);

                assert_eq!(lhs, rhs);
            }
            {
                let p1 = &G1Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());
                let q1 = &G2Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());
                let p2 = &G1Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());
                let q2 = -G2Affine::<Bls12381Curve>::generator()
                    * &Fr::<Bls12381Curve>::random(&mut thread_rng());

                let eq_class =
                    multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)]);
                let rhs = final_exponentiation(&eq_class) == Fp12::<Bls12381Curve>::one();
                let lhs = verify_pairing(&[p1, p2], &[q1, q2]);

                assert_eq!(lhs, rhs);
            }
        }
    }

    #[test]
    fn test_bilinearity() {
        for _ in 0..1 {
            {
                let a = Fr::<Bls12381Curve>::random(&mut thread_rng());
                let b = Fr::<Bls12381Curve>::random(&mut thread_rng());

                let p1 = &G1Affine::<Bls12381Curve>::generator() * a;
                let q1 = &G2Affine::<Bls12381Curve>::generator() * b;
                let p2 = p1.clone();
                let q2 = -q1.clone();

                let eq_class =
                    multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)]);
                let rhs = final_exponentiation(&eq_class) == Fp12::<Bls12381Curve>::one();
                let lhs = verify_pairing(&[p1, p2], &[q1, q2]);
                assert_eq!(lhs, rhs);
            }
            {
                let a = Fr::<Bls12381Curve>::random(&mut thread_rng());
                let b = Fr::<Bls12381Curve>::random(&mut thread_rng());

                let p1 = &G1Affine::<Bls12381Curve>::generator() * a;
                let q1 = &G2Affine::<Bls12381Curve>::generator() * b;
                let p2 = -p1.clone();
                let q2 = q1.clone();

                let eq_class =
                    multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)]);
                let rhs = final_exponentiation(&eq_class) == Fp12::<Bls12381Curve>::one();
                let lhs = verify_pairing(&[p1, p2], &[q1, q2]);
                assert_eq!(lhs, rhs);
            }
            {
                let a = Fr::<Bls12381Curve>::random(&mut thread_rng());
                let b = Fr::<Bls12381Curve>::random(&mut thread_rng());
                let c = Fr::<Bls12381Curve>::random(&mut thread_rng());

                let p1 = &G1Affine::<Bls12381Curve>::generator() * a;
                let q1 = &G2Affine::<Bls12381Curve>::generator() * b;
                let p2 = &G1Affine::<Bls12381Curve>::generator() * c;
                let q2 = G2Affine::<Bls12381Curve>::identity();

                let eq_class =
                    multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)]);
                let rhs = final_exponentiation(&eq_class) == Fp12::<Bls12381Curve>::one();
                let lhs = verify_pairing(&[p1, p2], &[q1, q2]);
                assert_eq!(lhs, rhs);
            }
            {
                let a = Fr::<Bls12381Curve>::random(&mut thread_rng());
                let b = Fr::<Bls12381Curve>::random(&mut thread_rng());
                let c = Fr::<Bls12381Curve>::random(&mut thread_rng());

                let p1 = &G1Affine::<Bls12381Curve>::generator() * a;
                let q1 = &G2Affine::<Bls12381Curve>::generator() * b;
                let p2 = G1Affine::<Bls12381Curve>::identity();
                let q2 = &G2Affine::<Bls12381Curve>::generator() * c;

                let eq_class =
                    multi_miller_loop(&[p1, p2], &[G2Prepared::from(q1), G2Prepared::from(q2)]);
                let rhs = final_exponentiation(&eq_class) == Fp12::<Bls12381Curve>::one();
                let lhs = verify_pairing(&[p1, p2], &[q1, q2]);
                assert_eq!(lhs, rhs);
            }
        }
    }

    #[test]
    fn test_kzg_proof() {
        // Define p_minus_y as G1Affine
        let p_minus_y: G1Affine<Bls12381Curve> = G1Affine::new(
            Fp::from_raw_unchecked([
                0xf1fbbbca6f146556,
                0xd97b05f5c8d900ac,
                0xc9dd98c56817e878,
                0x74e56183bb247c8f,
                0x7834a1246463e647,
                0x13efc82d2017e9c5,
            ]),
            Fp::from_raw_unchecked([
                0x5a32ad89774d101a,
                0x9f98f8b297c0a6ef,
                0xe6ced5347a619278,
                0x247e4e72ac6405a7,
                0x042314d67c5e1501,
                0x0d227532e5ff0df7,
            ]),
            false,
        );

        // Define G2Affine::generator
        let g = G2Affine::generator();

        // Define proof as G1Affine
        let proof = G1Affine::new(
            Fp::from_raw_unchecked([
                0x1fd63e07b7ebc348,
                0x9877db2329721299,
                0x086bc4061c7f6324,
                0x4b7ba36a0f40e2dc,
                0x71cefecd79e8274b,
                0x12c51ff81dd71dab,
            ]),
            Fp::from_raw_unchecked([
                0xf51a4d928d996704,
                0x38dec50cc439f51e,
                0x4b97bd173f0b4c01,
                0x866a52ef6cac68a2,
                0x4a1f6c1e6aa580c7,
                0x0bd8e916630efd7e,
            ]),
            false,
        );

        // Define x_minus_z as G2Affine
        let x_minus_z = G2Affine::new(
            Fp2::new(
                Fp::from_raw_unchecked([
                    0xcc42fbeea3ec7808,
                    0x3b8ef75ed5a1a1ca,
                    0x0633aea0b5b120db,
                    0xe8f1cfad5bc15230,
                    0x8fdc563271d71fb5,
                    0x0f58325a9123de31,
                ]),
                Fp::from_raw_unchecked([
                    0xe1dae420221cad3c,
                    0xd56cfc3004fc9708,
                    0x0a09e94596bb3128,
                    0x15f9874482328d1a,
                    0xeb06711c136096c8,
                    0x137758dc8afa2624,
                ]),
            ),
            Fp2::new(
                Fp::from_raw_unchecked([
                    0xd302bbe8bda77fb6,
                    0xa2b946d42ddc533b,
                    0x098cb463ef721dcd,
                    0x67e19da6e2eb73cd,
                    0xce92e1e270600b1f,
                    0x19e67ae2dd855f14,
                ]),
                Fp::from_raw_unchecked([
                    0x54f07cccf704a1ed,
                    0x3a0d28b8fc4bfb08,
                    0x37fbb971fb52744e,
                    0x0e96d351a42e2641,
                    0x46fc269fcd20eb3e,
                    0x048d9a7062ad9106,
                ]),
            ),
            false,
        );

        // Print residue output
        println!(
            "residue output: {:?}",
            verify_pairing(&[p_minus_y, proof], &[g, x_minus_z])
        );

        // Compute equivalence class using multi_miller_loop
        let eq_class = multi_miller_loop(
            &[p_minus_y, proof],
            &[G2Prepared::from(g), G2Prepared::from(x_minus_z)],
        );

        // Final exponentiation and comparison
        let rhs = final_exponentiation(&eq_class);
        println!("final exponentiation output: {:?}", rhs);
        println!(
            "final exponentiation output: {:?}",
            rhs == Fp12::<Bls12381Curve>::one()
        );
    }
}
