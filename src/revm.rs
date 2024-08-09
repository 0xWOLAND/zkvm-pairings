pub fn run_add(input: &[u8], gas_cost: u64, gas_limit: u64) -> PrecompileResult {
    if gas_cost > gas_limit {
        return Err(Error::OutOfGas.into());
    }

    let input = right_pad::<ADD_INPUT_LEN>(input);

    let p1 = read_point(&input[..64])?;
    let p2 = read_point(&input[64..])?;

    let mut output = [0u8; 64];
    if let Some(sum) = AffineG1::from_jacobian(p1 + p2) {
        sum.x().to_big_endian(&mut output[..32]).unwrap();
        sum.y().to_big_endian(&mut output[32..]).unwrap();
    }
    Ok(PrecompileOutput::new(gas_cost, output.into()))
}

pub fn run_mul(input: &[u8], gas_cost: u64, gas_limit: u64) -> PrecompileResult {
    if gas_cost > gas_limit {
        return Err(Error::OutOfGas.into());
    }

    let input = right_pad::<MUL_INPUT_LEN>(input);

    let p = read_point(&input[..64])?;

    // `Fr::from_slice` can only fail when the length is not 32.
    let fr = bn::Fr::from_slice(&input[64..96]).unwrap();

    let mut output = [0u8; 64];
    if let Some(mul) = AffineG1::from_jacobian(p * fr) {
        mul.x().to_big_endian(&mut output[..32]).unwrap();
        mul.y().to_big_endian(&mut output[32..]).unwrap();
    }
    Ok(PrecompileOutput::new(gas_cost, output.into()))
}

pub fn run_pair(
    input: &[u8],
    pair_per_point_cost: u64,
    pair_base_cost: u64,
    gas_limit: u64,
) -> PrecompileResult {
    let gas_used = (input.len() / PAIR_ELEMENT_LEN) as u64 * pair_per_point_cost + pair_base_cost;
    if gas_used > gas_limit {
        return Err(Error::OutOfGas.into());
    }

    if input.len() % PAIR_ELEMENT_LEN != 0 {
        return Err(Error::Bn128PairLength.into());
    }

    let success = if input.is_empty() {
        true
    } else {
        let elements = input.len() / PAIR_ELEMENT_LEN;

        let mut mul = Gt::one();
        for idx in 0..elements {
            let read_fq_at = |n: usize| {
                debug_assert!(n < PAIR_ELEMENT_LEN / 32);
                let start = idx * PAIR_ELEMENT_LEN + n * 32;
                // SAFETY: We're reading `6 * 32 == PAIR_ELEMENT_LEN` bytes from `input[idx..]`
                // per iteration. This is guaranteed to be in-bounds.
                let slice = unsafe { input.get_unchecked(start..start + 32) };
                Fq::from_slice(slice).map_err(|_| Error::Bn128FieldPointNotAMember)
            };
            let ax = read_fq_at(0)?;
            let ay = read_fq_at(1)?;
            let bay = read_fq_at(2)?;
            let bax = read_fq_at(3)?;
            let bby = read_fq_at(4)?;
            let bbx = read_fq_at(5)?;

            let a = new_g1_point(ax, ay)?;
            let b = {
                let ba = Fq2::new(bax, bay);
                let bb = Fq2::new(bbx, bby);
                if ba.is_zero() && bb.is_zero() {
                    G2::zero()
                } else {
                    G2::from(AffineG2::new(ba, bb).map_err(|_| Error::Bn128AffineGFailedToCreate)?)
                }
            };

            mul = mul * bn::pairing(a, b);
        }

        mul == Gt::one()
    };
    Ok(PrecompileOutput::new(gas_used, bool_to_bytes32(success)))
}
