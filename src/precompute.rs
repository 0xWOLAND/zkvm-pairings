use crate::{
    common::{AffinePoint, Curve},
    g1::G1Affine,
    g2::G2Affine,
};

pub(crate) fn precompute_lines<C: Curve>(q: &G2Affine<C>) -> Vec<G1Affine<C>> {
    let mut acc = q.clone();
    // (0..C::LOOP_COUNTER.len() - 1).rev().iter().for_each(|i| {
    //     acc = acc.double();
    // });
    todo!()
}
