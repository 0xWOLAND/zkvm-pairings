use crate::{
    common::{AffinePoint, Curve},
    fp::Fp,
    fp2::Fp2,
    fr::Fr,
    g1::G1Affine,
    g2::G2Affine,
};

/// 6U+2 for in NAF form

pub(crate) const SIX_U_PLUS_2_NAF: [i8; 65] = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0,
    0, 1, 0, 1, 1,
];

pub(crate) struct LineEvaluation<C: Curve> {
    x: Fp2<C>,
    y: Fp2<C>,
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
    p1: G2Affine<C>,
    p2: G2Affine<C>,
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
    p1: G2Affine<C>,
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

pub(crate) fn miller_loop<C: Curve>(p: &[G1Affine<C>], q: &[G2Affine<C>]) {
    assert_eq!(p.len(), q.len(), "Input slices must have the same length");
    assert!(p.len() > 0, "Cannot perform pairing on empty slices");

    q.iter().for_each(|qq| {})
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

        let (p1, l1, ll1) = double_and_add_step(p, q);
        let (_p2, l2) = double_step(p);
        let (p2, ll2) = add_step(_p2, q);

        assert_eq!(p1, p2);
    }

    #[test]
    fn test_triple_step() {
        let p = G2Affine::<Bls12381Curve>::random(&mut rand::thread_rng());

        let (p1, l1, l2) = triple_step(p);
        // let (_p2, ll1) = double_step(p);
        // let (p2, ll2) = add_step(_p2, p);
        let p2 = p + p + p;

        assert_eq!(p1, p2);
    }
}
