use std::f32::NEG_INFINITY;

use crate::{
    common::{AffinePoint, Curve},
    fp::Fp,
    fp12::Fp12,
    fp2::Fp2,
    fp6::Fp6,
    fr::Fr,
    g1::G1Affine,
    g2::G2Affine,
    precompute::{compute_lines, LineEvaluations},
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

pub(crate) fn add_step<C: Curve>(
    p1: G2Affine<C>,
    p2: G2Affine<C>,
) -> (G2Affine<C>, LineEvaluation<C>) {
    let n = p2.y - p1.y;
    let d = p2.x - p1.x;
    let l = n / d;

    let ll = l.square();
    let x_r = ll - (p1.x + p2.x);

    let y_r = l * (p1.x - x_r) - p1.y;

    let line = LineEvaluation {
        x: l,
        y: (l * p1.x) - p1.y,
    };

    let out = G2Affine::<C>::new(x_r, y_r, false);

    (out, line)
}

pub(crate) fn double_step<C: Curve>(p1: G2Affine<C>) -> (G2Affine<C>, LineEvaluation<C>) {
    let n = p1.x.square();
    let n = n + n + n;
    let d = p1.y + p1.y;
    let l = n / d;

    let x_r = l.square() - p1.x - p1.x;
    let y_r = l * (p1.x - x_r) - p1.y;

    let p = G2Affine::<C>::new(x_r, y_r, false);

    let line = LineEvaluation {
        x: l,
        y: (l - p1.x) - p1.y,
    };

    (p, line)
}

pub(crate) fn double_and_add_step<C: Curve>(
    p1: &G2Affine<C>,
    p2: &G2Affine<C>,
) -> (G2Affine<C>, LineEvaluation<C>, LineEvaluation<C>) {
    let l1 = (p1.y - p2.y) / (p1.x - p2.x);
    let line1 = LineEvaluation {
        x: l1,
        y: (l1 * p1.x) - p1.y,
    };
    let xspq = l1.square() - p1.x - p2.x;
    let l2 = -l1 - (p1.y + p1.y) / (xspq - p1.x);
    let line2 = LineEvaluation {
        x: l2,
        y: (l2 * p1.x) - p1.y,
    };
    let xr = l2.square() - p1.x - xspq;
    let yr = l2 * (p1.x - xr) - p1.y;

    let out = G2Affine::<C>::new(xr, yr, false);
    (out, line1, line2)
}

pub(crate) fn triple_step<C: Curve>(
    p1: &G2Affine<C>,
) -> (G2Affine<C>, LineEvaluation<C>, LineEvaluation<C>) {
    // λ1 = 3x²/2y
    let n = p1.x.square();
    let three = Fp::<C>::from(3);
    let n = n * three;
    let d = p1.y + p1.y;
    let l1 = n / d;

    // compute line1
    let line1 = LineEvaluation {
        x: l1,
        y: l1 * p1.x - p1.y,
    };

    // x2 = λ1² - 2x
    let x2 = l1.square() - (p1.x + p1.x);

    // compute λ2 = 2y / (x2 - x) - λ1
    let x1x2 = p1.x - x2;
    let l2 = (d / x1x2) - l1;

    // compute line2
    let line2 = LineEvaluation {
        x: l2,
        y: l2 * p1.x - p1.y,
    };

    // xr = λ2² - x2 - x
    let l2_square = l2.square();
    let x_r = l2_square - (x2 + p1.x);

    // yr = λ2 * (x - xr) - y
    let pxrx = p1.x - x_r;
    let y_r = (l2 * pxrx) - p1.y;

    let out = G2Affine::<C>::new(x_r, y_r, false);

    (out, line1, line2)
}

pub(crate) fn compute_tangent<C: Curve>(p: &G2Affine<C>) -> LineEvaluation<C> {
    let n = p.x.square();
    let three = Fp::<C>::from(3);
    let n = n * three;
    let d = p.y + p.y;
    let l = n / d;

    LineEvaluation {
        x: l,
        y: l * p.x - p.y,
    }
}

fn miller_loop_lines<C: Curve>(p: &[G1Affine<C>], lines: &[LineEvaluations<C>]) -> Fp12<C> {
    let n = p.len();
    assert!(n > 0, "Cannot perform pairing on empty slices");
    assert!(
        n == lines[0].len(),
        "Input slices must have the same length"
    );
    let y_inv = p
        .iter()
        .map(|p| p.y.invert().unwrap())
        .collect::<Vec<Fp<C>>>();
    let x_neg_over_y = p
        .iter()
        .zip(y_inv.iter())
        .map(|(p, y_inv)| -p.x * y_inv)
        .collect::<Vec<Fp<C>>>();

    let mut res = Fp12::<C>::one();

    res.c0.c0 = lines[0][0][62].y * y_inv[0];
    res.c0.c1 = lines[0][0][62].x * x_neg_over_y[0];
    res.c1.c1 = Fp2::one();

    let prod_lines = Fp12::mul_14_by_14(
        &(lines[0][1][62].x * y_inv[0]),
        &(lines[0][1][62].y * x_neg_over_y[0]),
        &res.c0.c0,
        &res.c0.c1,
    );

    let mut res = Fp12::new(
        Fp6::new(prod_lines[0], prod_lines[1], prod_lines[2]),
        Fp6::new(res.c1.c0, prod_lines[3], prod_lines[4]),
    );

    (1..n).for_each(|i| {
        res = res.mul_by_014(
            &(lines[i][0][62].y * y_inv[i]),
            &(lines[i][0][62].x * x_neg_over_y[i]),
        );
        res = res.mul_by_014(
            &(lines[i][1][62].y * y_inv[i]),
            &(lines[i][1][62].x * x_neg_over_y[i]),
        );
    });
    (0..62).rev().for_each(|i| {
        res = res.square();
        (0..n).for_each(|j| {
            if C::LOOP_COUNTER[j] == 0 {
                res = res.mul_by_014(
                    &(lines[i][0][j].y * y_inv[i]),
                    &(lines[i][0][j].x * x_neg_over_y[i]),
                );
            } else {
                res = res.mul_by_014(
                    &(lines[i][0][j].y * y_inv[i]),
                    &(lines[i][0][j].x * x_neg_over_y[i]),
                );
                res = res.mul_by_014(
                    &(lines[i][1][j].y * y_inv[i]),
                    &(lines[i][1][j].x * x_neg_over_y[i]),
                );
            }
        })
    });

    res = res.conjugate();

    res
}

pub(crate) fn miller_loop<C: Curve>(p: &[G1Affine<C>], q: &[G2Affine<C>]) -> Fp12<C> {
    assert_eq!(p.len(), q.len(), "Input slices must have the same length");
    assert!(p.len() > 0, "Cannot perform pairing on empty slices");

    let lines: Vec<LineEvaluations<C>> = q.iter().map(|q| compute_lines(q)).collect();
    miller_loop_lines(p, lines.as_slice())
}

#[cfg(test)]
mod test {
    use crate::common::Bls12381Curve;

    use super::*;

    #[test]
    fn test_dobule() {
        let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

        let (p1, l1) = double_step(p);
        let p2 = p.double();

        assert_eq!(p1, p2);
    }

    #[test]
    fn test_add() {
        let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        let q = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

        let (p1, l2) = add_step(p, q);
        let p2 = p + q;

        assert_eq!(p1, p2);
    }

    #[test]
    fn test_double_and_add() {
        let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());
        let q = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

        let (p1, l1, ll1) = double_and_add_step(&p, &q);
        let (_p2, l2) = double_step(p);
        let (p2, ll2) = add_step(_p2, q);

        assert_eq!(p1, p2);
    }

    #[test]
    fn test_triple_step() {
        let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

        let (p1, l1, l2) = triple_step(&p);
        // let (_p2, ll1) = double_step(p);
        // let (p2, ll2) = add_step(_p2, p);
        let p2 = p + p + p;

        assert_eq!(p1, p2);
    }
}
