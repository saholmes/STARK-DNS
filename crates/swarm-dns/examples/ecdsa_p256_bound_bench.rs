//! WITNESS-BINDING ECDSA-P256 verifier (the ECDSA soundness fix).
//!
//! Routes the single-row ECDSA-P256 v2 AIR through
//! `sub_air_with_trace::{prove,verify}_one_sub_air_with_trace`, using the
//! `trace_row`-gated merge + evaluator
//! (`deep_ali_merge_p256_ecdsa_v2_rowgated_streaming` /
//! `eval_ecdsa_verify_v2_rowgated_per_row`) so prove and verify agree
//! (the production merge's Lagrange `row0_indicator` gating is
//! incompatible with the generic trace-row verifier).
//!
//! DECISIVE TEST: an honest P256 signature ACCEPTS; a one-cell-tampered
//! trace REJECTS — i.e. the in-circuit ECDSA verdict binds the witness.
//!
//! Run: cargo run --release -p swarm-dns --example ecdsa_p256_bound_bench \
//!        --no-default-features --features sha3-256,mldsa-44,parallel

use std::time::Instant;

use ark_ff::Zero;
use ark_goldilocks::Goldilocks as F;
use sha2::{Digest as _, Sha256};

use p256::ecdsa::{signature::Signer, Signature as P256Signature, SigningKey, VerifyingKey};
use p256::elliptic_curve::sec1::ToEncodedPoint;

use deep_ali::{
    deep_ali_merge_p256_ecdsa_v2_rowgated_streaming,
    fri::DeepFriParams,
    p256_ecdsa::{
        reduce_digest_mod_n, verify as ecdsa_verify_native, PublicKey as EcdsaPublicKey,
        Signature as EcdsaSignature,
    },
    p256_ecdsa_air_v2::{
        build_ecdsa_verify_v2_layout, ecdsa_verify_v2_constraints,
        eval_ecdsa_verify_v2_rowgated_per_row, fill_ecdsa_verify_v2,
    },
    p256_field::NUM_LIMBS as P256_NUM_LIMBS,
    p256_group::GENERATOR as P256_GENERATOR,
    p256_scalar::ScalarElement,
    sub_air_with_trace::{
        prove_one_sub_air_with_trace_local, verify_one_sub_air_with_trace_local,
    },
};

const BLOWUP: usize = 32;
const PI_HASH: [u8; 32] = [0x33; 32];

fn scalar_to_msb_bits_256(s: &ScalarElement) -> Vec<bool> {
    let bytes = s.to_be_bytes();
    let mut bits = Vec::with_capacity(256);
    for byte in bytes.iter() {
        for shift in (0..8).rev() {
            bits.push((byte >> shift) & 1 == 1);
        }
    }
    bits
}

fn mk_params(n0: usize, r: usize, use_stir: bool, ph: [u8; 32]) -> DeepFriParams {
    DeepFriParams {
        schedule: (0..n0.trailing_zeros() as usize).map(|_| 2).collect(),
        r, seed_z: 0xDEEFu64, coeff_commit_final: true, d_final: 1,
        stir: use_stir, s0: r, public_inputs_hash: Some(ph),
    }
}

fn main() {
    let level = deep_ali::stark_level::NIST_LEVEL;
    let ext_deg = deep_ali::permutation_argument::EXT_DEGREE;
    let r = deep_ali::stark_level::NUM_QUERIES_LEVEL;
    let use_stir = matches!(std::env::var("BENCH_LDT").as_deref(), Ok("stir") | Ok("STIR"));
    let n_trace = std::env::var("NTRACE").ok().and_then(|s| s.parse().ok()).unwrap_or(64usize);
    eprintln!("=== ecdsa_p256_bound_bench: WITNESS-BINDING ECDSA-P256 v2, NIST L{level}, Fp{ext_deg}, r={r} ===");

    // ── 1. Real P256 signature (deterministic key) ──
    let sk = SigningKey::from_slice(&[0x42u8; 32]).expect("valid P256 key");
    let msg = b"STARK-DNS ECDSA-P256 witness-binding bench";
    let sig: P256Signature = sk.sign(msg);
    let vk = VerifyingKey::from(&sk);
    let ep = vk.to_encoded_point(false);
    let epb = ep.as_bytes();
    let mut qx = [0u8; 32]; qx.copy_from_slice(&epb[1..33]);
    let mut qy = [0u8; 32]; qy.copy_from_slice(&epb[33..65]);
    let sb = sig.to_bytes();
    let mut rbytes = [0u8; 32]; rbytes.copy_from_slice(&sb[0..32]);
    let mut sbytes = [0u8; 32]; sbytes.copy_from_slice(&sb[32..64]);
    let digest: [u8; 32] = Sha256::digest(msg).into();

    let pk = EcdsaPublicKey::from_be_bytes(&qx, &qy).expect("pk parse");
    let signature = EcdsaSignature::from_be_bytes(&rbytes, &sbytes).expect("sig parse");
    assert!(ecdsa_verify_native(&digest, &pk, &signature), "native ECDSA verify must hold");

    // ── 2. u1, u2 (FIPS 186-4 §6.4.2) + MSB bit decompositions ──
    let e = reduce_digest_mod_n(&digest);
    let w = signature.s.invert();
    let u_1 = e.mul(&w);
    let u_2 = signature.r.mul(&w);
    let u1_bits = scalar_to_msb_bits_256(&u_1);
    let u2_bits = scalar_to_msb_bits_256(&u_2);

    // ── 3. Layout (K=256) + single-row trace ──
    let g_x_base = 0;
    let g_y_base = P256_NUM_LIMBS;
    let g_z_base = 2 * P256_NUM_LIMBS;
    let q_x_base = 3 * P256_NUM_LIMBS;
    let q_y_base = 4 * P256_NUM_LIMBS;
    let q_z_base = 5 * P256_NUM_LIMBS;
    let start = 6 * P256_NUM_LIMBS;
    let (layout, total) = build_ecdsa_verify_v2_layout(
        start, g_x_base, g_y_base, g_z_base, q_x_base, q_y_base, q_z_base, 256,
    );
    let kk = ecdsa_verify_v2_constraints(&layout);

    let build_trace = |u1b: &[bool], u2b: &[bool], r_scalar: &ScalarElement| -> Vec<Vec<F>> {
        let mut trace: Vec<Vec<F>> = (0..total).map(|_| vec![F::zero(); n_trace]).collect();
        let mut row0: Vec<Vec<F>> = (0..total).map(|_| vec![F::zero(); 1]).collect();
        let g = *P256_GENERATOR;
        fill_ecdsa_verify_v2(&mut row0, 0, &layout, &g.x, &g.y, &pk.point.x, &pk.point.y, u1b, u2b, r_scalar);
        for c in 0..total { trace[c][0] = row0[c][0]; }
        trace
    };

    let prove_verify = |trace: &[Vec<F>]| -> (f64, f64, usize, bool) {
        // ECDSA v2 constraints are LOCAL (cur-only) -> use the local
        // witness-binding path that omits the next-row openings (half the
        // payload).  Sound: the verifier check c_eval*Z_H = Sum a_j Phi_j(cur)
        // never reads nxt for a local AIR.
        let t0 = Instant::now();
        let proof = prove_one_sub_air_with_trace_local(
            trace, n_trace, BLOWUP, PI_HASH, b"ecdsa_v2_bound", kk,
            |lde, nt, bw, cc| deep_ali_merge_p256_ecdsa_v2_rowgated_streaming(lde, cc, &layout, nt, bw).0,
            |n0, ph| mk_params(n0, r, use_stir, ph),
        );
        let prove_ms = t0.elapsed().as_secs_f64() * 1000.0;
        let t0 = Instant::now();
        let ok = verify_one_sub_air_with_trace_local(
            &proof, n_trace, BLOWUP, PI_HASH, b"ecdsa_v2_bound", total, kk,
            |cur, row| eval_ecdsa_verify_v2_rowgated_per_row(cur, row, &layout),
            |n0, ph| mk_params(n0, r, use_stir, ph),
        ).is_ok();
        (prove_ms, t0.elapsed().as_secs_f64() * 1000.0, proof.fri_proof_bytes.len(), ok)
    };

    // ── Honest: ACCEPT ──
    let honest = build_trace(&u1_bits, &u2_bits, &signature.r);
    let (p_ms, v_ms, fri_b, ok) = prove_verify(&honest);
    eprintln!("[honest]   prove {p_ms:.1} ms, verify {v_ms:.2} ms, fri {} KiB -> verify={ok}", fri_b / 1024);
    assert!(ok, "BINDING BROKEN: honest ECDSA signature must verify");

    // ── Tampered: corrupt a constrained row-0 cell (a G coordinate limb) -> REJECT ──
    let mut bad = honest.clone();
    bad[q_x_base][0] += F::from(1u64); // perturb the public-key x limb 0
    let (_, _, _, bad_ok) = prove_verify(&bad);
    eprintln!("[tampered] perturb Q_x[0] -> verify={bad_ok}");

    println!(
        "ecdsa_p256_bound level=L{level} field=Fp{ext_deg} n_trace={n_trace} r={r} \
         prove_ms={p_ms:.1} verify_ms={v_ms:.2} honest_verify={ok} tampered_verify={bad_ok}"
    );
    if !bad_ok {
        println!("=> WITNESS-BINDING WORKS for ECDSA-P256: honest accepts, tampered REJECTS");
    } else {
        println!("=> ECDSA STILL NON-BINDING (or Q_x[0] unconstrained — try another cell)");
    }
}
