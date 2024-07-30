use std::{f32::consts::LN_10, str::FromStr};

use num_bigint::BigUint;

use crate::{common::Curve, fp12::Fp12, fp6::Fp6};

fn residue_test_hint<C: Curve>(x: Fp12<C>) -> (Fp12<C>, Fp12<C>) {
    let p_factor = BigUint::from_str("5044125407647214251").unwrap();
    let exp_factor = BigUint::from_str("2366356426548243601069753987687709088104621721678962410379583120840019275952471579477684846670499039076873213559162845121989217658133790336552276567078487633052653005423051750848782286407340332979263075575489766963251914185767058009683318020965829271737924625612375201545022326908440428522712877494557944965298566001441468676802477524234094954960009227631543471415676620753242466901942121887152806837594306028649150255258504417829961387165043999299071444887652375514277477719817175923289019181393803729926249507024121957184340179467502106891835144220611408665090353102353194448552304429530104218473070114105759487413726485729058069746063140422361472585604626055492939586602274983146215294625774144156395553405525711143696689756441298365274341189385646499074862712688473936093315628166094221735056483459332831845007196600723053356837526749543765815988577005929923802636375670820616189737737304893769679803809426304143627363860243558537831172903494450556755190448279875942974830469855835666815454271389438587399739607656399812689280234103023464545891697941661992848552456326290792224091557256350095392859243101357349751064730561345062266850238821755009430903520645523345000326783803935359711318798844368754833295302563158150573540616830138810935344206231367357992991289265295323280").unwrap();

    let exp = &exp_factor * BigUint::from(27u32);
    let root = x.pow_vartime_extended(&exp.to_u64_digits());

    // pth-root inverse
    let root_p_inv = match root.is_one() {
        true => Fp12::one(),
        false => {
            let exp_inv = exp.modinv(&p_factor).unwrap();
            let exp = (&p_factor - exp_inv) % &p_factor;
            root.pow_vartime_extended(&exp.to_u64_digits())
        }
    };

    // 3rd-root inverse
    let exp = &exp_factor * &p_factor;
    let mut root = x.pow_vartime_extended(&exp.to_u64_digits());
    let mut deg_3_pow = 0u32;

    for i in 0..=3 {
        if root.is_one() {
            deg_3_pow = i;
            break;
        }
        if i != 3 {
            root = root.pow_vartime_extended(&[3])
        }
    }

    let root_27_inv: Fp12<C> = match deg_3_pow {
        0 => Fp12::one(),
        _ => {
            let deg_3 = BigUint::from(3u32).pow(deg_3_pow);
            let exp = &exp_factor * p_factor; // should have already done this
                                              // let mut root = x.pow_vartime_extended(&exp.to_u64_digits());
            let exp_inv = exp.modinv(&deg_3).unwrap();
            let exp = (&deg_3 - exp_inv) % deg_3;
            root.pow_vartime_extended(&exp.to_u64_digits())
        }
    };

    let w_full = root_p_inv * root_27_inv;
    let x = x * w_full;

    let lam = BigUint::from_str(C::LAMBDA).unwrap();
    let exp = lam.modinv(&exp_factor).unwrap();
    let w = x.pow_vartime_extended(&exp.to_u64_digits());

    (w, Fp12::new(w_full.c0, Fp6::zero()))
}
