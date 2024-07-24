use crate::{fp::Fp, scalar::Scalar};

/// Constant representing the modulus
/// q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
const MODULUS: Scalar = Scalar([
    0xffff_ffff_0000_0001,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
]);

struct G1Affine {
    x: Fp,
    y: Fp,
    infinity: bool,
}

impl PartialEq for G1Affine {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.infinity == other.infinity
    }
}

impl G1Affine {
    pub fn new(x: Fp, y: Fp, infinity: bool) -> Self {
        G1Affine { x, y, infinity }
    }

    pub const fn zero() -> Self {
        G1Affine {
            x: Fp::one(),
            y: Fp::one(),
            infinity: false,
        }
    }
}
