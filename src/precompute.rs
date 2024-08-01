use crate::{
    common::{AffinePoint, Curve},
    fp2::Fp2,
    g1::G1Affine,
    g2::G2Affine,
    pairings::{compute_tangent, double_and_add_step, double_step, triple_step, LineEvaluation},
};

pub type LineEvaluations<C> = [[LineEvaluation<C>; 64]; 2]; // 64 >= LOOP_COUNTER.len() - 1 based on the current curve

pub(crate) fn compute_lines<C: Curve>(q: &G2Affine<C>) -> LineEvaluations<C> {
    let n = C::LOOP_COUNTER.len();
    let mut lines: LineEvaluations<C> = [[LineEvaluation::zero(); 64]; 2];

    let (mut acc, l1, l2) = triple_step(q);
    lines[0][n - 2] = l1;
    lines[1][n - 2] = l2;

    (0..=(n - 3)).rev().for_each(|i| {
        if C::LOOP_COUNTER[i] == 0 {
            let (acc1, l1) = double_step(acc);
            acc = acc1;
            lines[0][i] = l1;
            lines[1][i] = LineEvaluation::zero();
        } else {
            let (acc1, l1, l2) = double_and_add_step(&acc, q);
            acc = acc1;
            lines[0][i] = l1;
            lines[1][i] = l2;
        }
    });

    lines[0][0] = compute_tangent(&acc);
    lines[1][0] = LineEvaluation::zero();
    lines
}
