use std::borrow::Cow;

/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Right-pads the given slice with zeroes until `LEN`.
///
/// Returns the first `LEN` bytes if it does not need padding.
#[inline]
pub fn right_pad<const LEN: usize>(data: &[u8]) -> Cow<'_, [u8; LEN]> {
    if let Some(data) = data.get(..LEN) {
        Cow::Borrowed(data.try_into().unwrap())
    } else {
        let mut padded = [0; LEN];
        padded[..data.len()].copy_from_slice(data);
        Cow::Owned(padded)
    }
}
