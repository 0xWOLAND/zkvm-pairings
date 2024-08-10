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

pub fn run_add(input: &[u8]) -> Vec<u8> {
    let input = right_pad::<ADD_INPUT_LEN>(input);
    let p1 = read_point(&input[..64]).expect("Failed to read point 1");
    let p2 = read_point(&input[64..]).expect("Failed to read point 2");

    let sum = p1 + p2;
    let mut bytes = [0u8; 64];

    bytes[..32].copy_from_slice(&sum.x.to_bytes());
    bytes[32..].copy_from_slice(&sum.y.to_bytes());

    bytes.to_vec()
}

pub fn run_mul(input: &[u8]) -> Vec<u8> {
    let input = right_pad::<MUL_INPUT_LEN>(input);
    let p = read_point(&input[..64]).unwrap();
    let fr = Fr::<Bn254>::from_bytes_be(&input[64..96].try_into().unwrap()).unwrap();

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_point() {
        let input = hex::decode("20b781dd0db3b7980a4b3814128c86e597e1442d0fc9eb7f932a5229494d6b7917d1cef436eb2f665670c7b34854e62c227043a7b111a539c0295518bbab3ca91e57cd7d385ce3b0d436ec73d61caaa8290ef07388228daa83627d6b2eb3e41d0000000000000000000000000000000000000000000000000000000000000000").unwrap();
        let expected = hex::decode("239f13d5597aa26f424e4afefb9d82396ae2bdb102c9ec90c68b4254c2769da90a36984eadebffe544f078bbaadb33193fda7bd1aed33e683399c13ac24fc433").unwrap();
        let out = run_mul(&input);
        assert_eq!(out, expected);
    }
}

/*
#[cfg(test)]
mod tests {
    use crate::test_utils::new_context;

    use super::*;

    #[test]
    fn test_alt_bn128_add() {
        let input = hex::decode(
            "\
             18b18acfb4c2c30276db5411368e7185b311dd124691610c5d3b74034e093dc9\
             063c909c4720840cb5134cb9f59fa749755796819658d32efc0d288198f37266\
             07c2b7f58a84bd6145f00c9c2bc0bb1a187f20ff2c92963a88019e7c6a014eed\
             06614e20c147e940f2d70da3f74c9a17df361706a4485c742bd6788478fa17d7",
        )
        .unwrap();
        let expected = hex::decode(
            "\
            2243525c5efd4b9c3d3c45ac0ca3fe4dd85e830a4ce6b65fa1eeaee202839703\
            301d1d33be6da8e509df21cc35964723180eed7532537db9ae5e7d48f195c915",
        )
        .unwrap();

        let res = Bn128Add::<Byzantium>::run(&input, 500, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // zero sum test
        let input = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();
        let expected = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();

        let res = Bn128Add::<Byzantium>::run(&input, 500, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // out of gas test
        let input = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();
        let res = Bn128Add::<Byzantium>::run(&input, 499, &new_context(), false);
        assert!(matches!(res, Err(Return::OutOfGas)));

        // no input test
        let input = [0u8; 0];
        let expected = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();

        let res = Bn128Add::<Byzantium>::run(&input, 500, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // point not on curve fail
        let input = hex::decode(
            "\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111",
        )
        .unwrap();

        let res = Bn128Add::<Byzantium>::run(&input, 500, &new_context(), false);
        assert!(matches!(
            res,
            Err(Return::Other(Cow::Borrowed("ERR_BN128_INVALID_POINT")))
        ));
    }

    #[test]
    fn test_alt_bn128_mul() {
        let input = hex::decode(
            "\
            2bd3e6d0f3b142924f5ca7b49ce5b9d54c4703d7ae5648e61d02268b1a0a9fb7\
            21611ce0a6af85915e2f1d70300909ce2e49dfad4a4619c8390cae66cefdb204\
            00000000000000000000000000000000000000000000000011138ce750fa15c2",
        )
        .unwrap();
        let expected = hex::decode(
            "\
            070a8d6a982153cae4be29d434e8faef8a47b274a053f5a4ee2a6c9c13c31e5c\
            031b8ce914eba3a9ffb989f9cdd5b0f01943074bf4f0f315690ec3cec6981afc",
        )
        .unwrap();

        let res = Bn128Mul::<Byzantium>::run(&input, 40_000, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // out of gas test
        let input = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000\
            0200000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();
        let res = Bn128Mul::<Byzantium>::run(&input, 39_999, &new_context(), false);
        assert!(matches!(res, Err(Return::OutOfGas)));

        // zero multiplication test
        let input = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000\
            0200000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();
        let expected = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();

        let res = Bn128Mul::<Byzantium>::run(&input, 40_000, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // no input test
        let input = [0u8; 0];
        let expected = hex::decode(
            "\
            0000000000000000000000000000000000000000000000000000000000000000\
            0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();

        let res = Bn128Mul::<Byzantium>::run(&input, 40_000, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // point not on curve fail
        let input = hex::decode(
            "\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            0f00000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();

        let res = Bn128Mul::<Byzantium>::run(&input, 40_000, &new_context(), false);
        assert!(matches!(
            res,
            Err(Return::Other(Cow::Borrowed("ERR_BN128_INVALID_POINT")))
        ));
    }

    #[test]
    fn test_alt_bn128_pair() {
        let input = hex::decode(
            "\
            1c76476f4def4bb94541d57ebba1193381ffa7aa76ada664dd31c16024c43f59\
            3034dd2920f673e204fee2811c678745fc819b55d3e9d294e45c9b03a76aef41\
            209dd15ebff5d46c4bd888e51a93cf99a7329636c63514396b4a452003a35bf7\
            04bf11ca01483bfa8b34b43561848d28905960114c8ac04049af4b6315a41678\
            2bb8324af6cfc93537a2ad1a445cfd0ca2a71acd7ac41fadbf933c2a51be344d\
            120a2a4cf30c1bf9845f20c6fe39e07ea2cce61f0c9bb048165fe5e4de877550\
            111e129f1cf1097710d41c4ac70fcdfa5ba2023c6ff1cbeac322de49d1b6df7c\
            2032c61a830e3c17286de9462bf242fca2883585b93870a73853face6a6bf411\
            198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2\
            1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed\
            090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b\
            12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa",
        )
        .unwrap();
        let expected =
            hex::decode("0000000000000000000000000000000000000000000000000000000000000001")
                .unwrap();

        let res = Bn128Pair::<Byzantium>::run(&input, 260_000, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // out of gas test
        let input = hex::decode(
            "\
            1c76476f4def4bb94541d57ebba1193381ffa7aa76ada664dd31c16024c43f59\
            3034dd2920f673e204fee2811c678745fc819b55d3e9d294e45c9b03a76aef41\
            209dd15ebff5d46c4bd888e51a93cf99a7329636c63514396b4a452003a35bf7\
            04bf11ca01483bfa8b34b43561848d28905960114c8ac04049af4b6315a41678\
            2bb8324af6cfc93537a2ad1a445cfd0ca2a71acd7ac41fadbf933c2a51be344d\
            120a2a4cf30c1bf9845f20c6fe39e07ea2cce61f0c9bb048165fe5e4de877550\
            111e129f1cf1097710d41c4ac70fcdfa5ba2023c6ff1cbeac322de49d1b6df7c\
            2032c61a830e3c17286de9462bf242fca2883585b93870a73853face6a6bf411\
            198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2\
            1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed\
            090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b\
            12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa",
        )
        .unwrap();
        let res = Bn128Pair::<Byzantium>::run(&input, 259_999, &new_context(), false);
        assert!(matches!(res, Err(Return::OutOfGas)));

        // no input test
        let input = [0u8; 0];
        let expected =
            hex::decode("0000000000000000000000000000000000000000000000000000000000000001")
                .unwrap();

        let res = Bn128Pair::<Byzantium>::run(&input, 260_000, &new_context(), false)
            .unwrap()
            .output;
        assert_eq!(res, expected);

        // point not on curve fail
        let input = hex::decode(
            "\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111",
        )
        .unwrap();

        let res = Bn128Pair::<Byzantium>::run(&input, 260_000, &new_context(), false);
        assert!(matches!(
            res,
            Err(Return::Other(Cow::Borrowed("ERR_BN128_INVALID_A")))
        ));

        // invalid input length
        let input = hex::decode(
            "\
            1111111111111111111111111111111111111111111111111111111111111111\
            1111111111111111111111111111111111111111111111111111111111111111\
            111111111111111111111111111111\
        ",
        )
        .unwrap();

        let res = Bn128Pair::<Byzantium>::run(&input, 260_000, &new_context(), false);
        assert!(matches!(
            res,
            Err(Return::Other(Cow::Borrowed("ERR_BN128_INVALID_LEN",)))
        ));
    }
}
*/
