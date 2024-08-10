use substrate_bn::G1;

use crate::{
    fp::{Bn254, FpElement},
    fp2::Fp2,
    fr::Fr,
    g1::G1Affine,
    g2::G2Affine,
    pairings::verify_pairing,
    utils::right_pad,
};

/// Input length for the add operation.
/// `ADD` takes two uncompressed G1 points (64 bytes each).
pub const ADD_INPUT_LEN: usize = 64 + 64;

/// Input length for the multiplication operation.
/// `MUL` takes an uncompressed G1 point (64 bytes) and scalar (32 bytes).
pub const MUL_INPUT_LEN: usize = 64 + 32;

/// Pair element length.
/// `PAIR` elements are composed of an uncompressed G1 point (64 bytes) and an uncompressed G2 point
/// (128 bytes).
pub const PAIR_ELEMENT_LEN: usize = 64 + 128;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum PrecompileError {
    /// out of gas is the main error. Others are here just for completeness
    OutOfGas,
    // Blake2 errors
    Blake2WrongLength,
    Blake2WrongFinalIndicatorFlag,
    // Modexp errors
    ModexpExpOverflow,
    ModexpBaseOverflow,
    ModexpModOverflow,
    // Bn128 errors
    Bn128FieldPointNotAMember,
    Bn128AffineGFailedToCreate,
    Bn128PairLength,
    // Blob errors
    /// The input length is not exactly 192 bytes.
    BlobInvalidInputLength,
    /// The commitment does not match the versioned hash.
    BlobMismatchedVersion,
    /// The proof verification failed.
    BlobVerifyKzgProofFailed,
}

/// A precompile operation result.
///
/// Returns either `Ok((gas_used, return_bytes))` or `Err(error)`.
pub type PrecompileResult = Result<(u64, Vec<u8>), PrecompileError>;

/// Reads a single `Fq` from the input slice.
///
/// # Panics
///
/// Panics if the input is not at least 32 bytes long.
#[inline]
pub fn read_fq(input: &[u8]) -> Result<Bn254, String> {
    match Bn254::from_bytes_be(&input[..32].try_into().unwrap()) {
        Some(fq) => Ok(fq),
        None => Err("Fq is not in the field".to_string()),
    }
}

/// Reads the `x` and `y` points from the input slice.
///
/// # Panics
///
/// Panics if the input is not at least 64 bytes long.
#[inline]
pub fn read_point(input: &[u8]) -> Result<G1Affine<Bn254>, String> {
    let px = read_fq(&input[0..32])?;
    let py = read_fq(&input[32..64])?;

    G1Affine::<Bn254>::new(px, py).ok_or("Point is not on the curve".to_string())
}

pub fn run_add(input: &[u8], gas_cost: u64, gas_limit: u64) -> Vec<u8> {
    let input = right_pad::<ADD_INPUT_LEN>(input);
    let p1 = read_point(&input[..64]).expect("Failed to read point 1");
    let p2 = read_point(&input[64..]).expect("Failed to read point 2");

    let mut output = [0u8; 64];
    let sum = p1 + p2;

    let mut bytes = [0u8; 64];
    bytes[..32].copy_from_slice(&sum.x.to_bytes());
    bytes[32..].copy_from_slice(&sum.y.to_bytes());

    bytes.to_vec()
}

pub fn run_mul(input: &[u8], gas_cost: u64, gas_limit: u64) -> Vec<u8> {
    assert!(gas_cost <= gas_limit, "Gas cost exceeds gas limit");

    let input = right_pad::<MUL_INPUT_LEN>(input);
    let p = read_point(&input[..64]).unwrap();
    let fr = Fr::<Bn254>::from_bytes(&input[64..96].try_into().unwrap()).unwrap();

    let mul = p * fr;
    let mut bytes = [0u8; 64];

    bytes[..32].copy_from_slice(&mul.x.to_bytes());
    bytes[32..].copy_from_slice(&mul.y.to_bytes());
    bytes.to_vec()
}

pub fn run_pair(
    input: &[u8],
    pair_per_point_cost: u64,
    pair_base_cost: u64,
    gas_limit: u64,
) -> PrecompileResult {
    let gas_used = (input.len() / PAIR_ELEMENT_LEN) as u64 * pair_per_point_cost + pair_base_cost;
    if gas_used > gas_limit {
        return Err(PrecompileError::OutOfGas);
    }

    if input.len() % PAIR_ELEMENT_LEN != 0 {
        return Err(PrecompileError::Bn128PairLength);
    }

    let output = if input.is_empty() {
        Fr::<Bn254>::one()
    } else {
        let elements = input.len() / PAIR_ELEMENT_LEN;
        let mut vals = Vec::with_capacity(elements);

        const PEL: usize = PAIR_ELEMENT_LEN;

        for idx in 0..elements {
            let mut buf = [0u8; 32];

            buf.copy_from_slice(&input[(idx * PEL)..(idx * PEL + 32)]);
            let ax = Bn254::from_bytes_be(&buf).unwrap();
            buf.copy_from_slice(&input[(idx * PEL + 32)..(idx * PEL + 64)]);
            let ay = Bn254::from_bytes_be(&buf).unwrap();
            buf.copy_from_slice(&input[(idx * PEL + 64)..(idx * PEL + 96)]);
            let bay = Bn254::from_bytes_be(&buf).unwrap();
            buf.copy_from_slice(&input[(idx * PEL + 96)..(idx * PEL + 128)]);
            let bax = Bn254::from_bytes_be(&buf).unwrap();
            buf.copy_from_slice(&input[(idx * PEL + 128)..(idx * PEL + 160)]);
            let bby = Bn254::from_bytes_be(&buf).unwrap();
            buf.copy_from_slice(&input[(idx * PEL + 160)..(idx * PEL + 192)]);
            let bbx = Bn254::from_bytes_be(&buf).unwrap();

            let a = {
                if ax.is_zero() && ay.is_zero() {
                    G1Affine::<Bn254>::zero()
                } else {
                    G1Affine::<Bn254>::new(ax, ay).unwrap()
                }
            };
            let b = {
                let ba = Fp2::<Bn254>::new(bax, bay);
                let bb = Fp2::<Bn254>::new(bbx, bby);

                if ba.is_zero() && bb.is_zero() {
                    G2Affine::<Bn254>::zero()
                } else {
                    // G2::from(AffineG2::new(ba, bb).map_err(|_| Error::Bn128AffineGFailedToCreate)?)
                    G2Affine::<Bn254>::new(ba, bb, false).unwrap()
                }
            };
            vals.push((a, b))
        }

        match verify_pairing(&vals) {
            true => Fr::<Bn254>::one(),
            false => Fr::<Bn254>::zero(),
        }
    };

    Ok((gas_used, output.to_bytes().to_vec()))
}
