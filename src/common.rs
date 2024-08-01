use field::{
    common::{Bls12381Curve, Curve as FieldCurve},
    fp12::Fp12,
};
use sp1_zkvm::{io::FD_BLS_HINT, lib::io};

pub trait Curve: FieldCurve {
    const LOG_ATE_LOOP_COUNT: usize;
    const LOOP_COUNTER: [u64; 64];
    fn get_root_and_scaling_factor(x: &Fp12<impl FieldCurve>) -> (Fp12<Self>, Fp12<Self>);
}

impl Curve for Bls12381Curve {
    fn get_root_and_scaling_factor(x: &Fp12<impl FieldCurve>) -> (Fp12<Self>, Fp12<Self>) {
        io::write(FD_BLS_HINT, x.to_bytes().as_ref());
        // let (root_buf, scaling_factor_buf) = buf.split_at(buf.len() / 2);
        // println!("root buf len: {}", root_buf.len());
        // println!("scaling factor buf len: {}", scaling_factor_buf.len());
        // let root_buf: [u8; 576] = root_buf.try_into().unwrap();
        // let scaling_factor_buf: [u8; 576] = scaling_factor_buf.try_into().unwrap();
        let root_buf: [u8; 576] = io::read_vec().try_into().unwrap();
        let scaling_factor_buf: [u8; 576] = io::read_vec().try_into().unwrap();
        let root = Fp12::<Bls12381Curve>::from_bytes(&root_buf);
        let scaling_factor = Fp12::<Bls12381Curve>::from_bytes(&scaling_factor_buf);
        (root, scaling_factor)
    }

    const LOG_ATE_LOOP_COUNT: usize = 62;
    const LOOP_COUNTER: [u64; 64] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        1, 0, 1, 1,
    ];
}
