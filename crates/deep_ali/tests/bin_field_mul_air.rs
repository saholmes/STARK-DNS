//! **Binary-verify AIR over deep_ali (prototype).** The prime-field outer STARK (pq-rollup
//! tests/prime_outer_stark.rs) needs the binary-field verify arithmetised as a Goldilocks AIR. This
//! builds the ATOMIC unit of that arithmetisation -- the boolean multiply-accumulate that a carryless
//! (GF(2^k)) multiply is built from -- as a real deep_ali AIR, runs it through the DEEP-ALI merge to
//! ONE composition polynomial, FRI-proves it, and verifies it. The single merged polynomial is exactly
//! what lets the QROM proof reason about one random oracle.
//!
//! AIR (width 4, 4 degree-2 constraints): columns [a, b, p, acc].
//!   a*(a-1)=0, b*(b-1)=0   (booleanity)
//!   p - a*b = 0            (AND)
//!   nxt.acc - (acc + p - 2*acc*p) = 0   (XOR-accumulate: the carryless-mul reduction step)
//!
//! A GF(2^128) carryless multiply is ~33,297 such boolean ops (pq-rollup prime_outer_stark), so one
//! binary multiply ~= a 2^15-row instance of this AIR. We measure prove/verify/size at a few sizes.
//!
//! Run: cargo test --release -p deep_ali --test bin_field_mul_air -- --nocapture

use std::time::Instant;

use ark_ff::{Field, One, Zero};
use ark_goldilocks::Goldilocks as F;

use deep_ali::cubic_ext::{ext_mul_count, reset_ext_mul_count};
use deep_ali::deep_ali_merge_per_row_no_layout;
use deep_ali::fri::DeepFriParams;
use deep_ali::sub_air_with_trace::{
    prove_one_sub_air_with_trace, serialize_proof, verify_one_sub_air_with_trace,
};

const WIDTH: usize = 4;
const NUM_CONSTRAINTS: usize = 4;

/// The AIR's per-row transition constraints (shared by prover and verifier).
fn eval_per_row(cur: &[F], nxt: &[F], _row: usize) -> Vec<F> {
    let (a, b, p, acc) = (cur[0], cur[1], cur[2], cur[3]);
    let two = F::from(2u64);
    vec![
        a * a - a,                          // booleanity(a)
        b * b - b,                          // booleanity(b)
        p - a * b,                          // AND
        nxt[3] - (acc + p - two * acc * p), // XOR-accumulate
    ]
}

/// A valid boolean multiply-accumulate trace (or a tampered one that breaks the AND constraint).
fn build_trace(n: usize, tamper: bool) -> Vec<Vec<F>> {
    let (mut a, mut b, mut p, mut acc) =
        (vec![F::zero(); n], vec![F::zero(); n], vec![F::zero(); n], vec![F::zero(); n]);
    let mut acc_bit = 0u64;
    for i in 0..n {
        let ab = ((i.wrapping_mul(2654435761)) & 1) as u64;
        let bb = ((i.wrapping_mul(40503).wrapping_add(7)) & 1) as u64;
        let pb = ab & bb;
        a[i] = F::from(ab);
        b[i] = F::from(bb);
        p[i] = F::from(pb);
        acc[i] = F::from(acc_bit);
        acc_bit ^= pb;
    }
    if tamper {
        p[n / 2] += F::one(); // break p = a*b at one row
    }
    vec![a, b, p, acc]
}

fn params(n0: usize, pi_hash: [u8; 32]) -> DeepFriParams {
    // NIST L1 UNCONDITIONAL (Johnson): binary (arity-2) folding, r=128 queries at rho=1/4 (blowup 4)
    // gives 128*1 bit = 128-bit Johnson soundness (delta=1-sqrt(1/4)=1/2 -> 1 bit/query). No conjecture.
    DeepFriParams {
        schedule: vec![2usize; n0.trailing_zeros() as usize],
        r: 128,
        seed_z: 0xDEEF_BAAD,
        coeff_commit_final: true,
        d_final: 1,
        stir: false,
        s0: 16,
        public_inputs_hash: Some(pi_hash),
    }
}

/// The DEEP-ALI merge: fold all constraints/columns into ONE composition polynomial to commit.
fn c_eval(lde: &[Vec<F>], n_trace: usize, blowup: usize, comb: &[F]) -> Vec<F> {
    deep_ali_merge_per_row_no_layout(lde, comb, F::one(), n_trace, blowup, WIDTH, NUM_CONSTRAINTS, eval_per_row).0
}

#[test]
fn bin_field_mul_air_end_to_end() {
    let pi_hash = [0x11u8; 32];
    let dsep = b"pq-rollup/bin-mul-air/v1";
    let blowup = 4usize;

    // Measured on-chain gas anchors (programs/evm-gas): base Goldilocks mul = native MULMOD (91 gas,
    // loop-inclusive); a cubic-ext (Ext) mul = 9 base muls (schoolbook) ~= 819 gas.
    const BASE_MUL_GAS: u64 = 91;
    const BASE_PER_EXT: u64 = 9;
    let gwei = 15.0f64;
    let eth = 3500.0f64;

    println!("\n===== Binary-field multiply AIR over deep_ali (Goldilocks; DEEP-ALI merge -> 1 poly) =====");
    println!("width={WIDTH}, constraints={NUM_CONSTRAINTS} (degree 2). One GF(2^128) mul ~= 2^15 rows.");
    println!("On-chain: base Goldilocks mul = {BASE_MUL_GAS} gas (MULMOD), Ext mul = {BASE_PER_EXT} base.");
    println!("{:>8} {:>9} {:>9} {:>10} {:>10} {:>11} {:>9}", "n_trace", "prove ms", "vfy ms", "proof B", "Ext muls", "verify gas", "$/verify");

    let mut samples: Vec<(f64, f64)> = Vec::new(); // (log2 n_trace, Ext muls) for extrapolation
    for &nt in &[4096usize, 16384, 32768] {
        let trace = build_trace(nt, false);
        let t0 = Instant::now();
        let proof = prove_one_sub_air_with_trace(&trace, nt, blowup, pi_hash, dsep, NUM_CONSTRAINTS, c_eval, params);
        let prove_ms = t0.elapsed().as_secs_f64() * 1e3;
        let bytes = serialize_proof(&proof).len();

        reset_ext_mul_count();
        let t1 = Instant::now();
        let ok = verify_one_sub_air_with_trace(&proof, nt, blowup, pi_hash, dsep, WIDTH, NUM_CONSTRAINTS, eval_per_row, params);
        let verify_ms = t1.elapsed().as_secs_f64() * 1e3;
        let ext_muls = ext_mul_count();
        assert!(ok.is_ok(), "honest AIR must verify at n_trace={nt}: {ok:?}");

        let verify_gas = ext_muls * BASE_PER_EXT * BASE_MUL_GAS;
        let usd = verify_gas as f64 * gwei * 1e-9 * eth;
        println!("{nt:>8} {prove_ms:>9.1} {verify_ms:>9.2} {bytes:>10} {ext_muls:>10} {verify_gas:>9} gas ${usd:>7.2}");
        samples.push(((nt as f64).log2(), ext_muls as f64));
    }

    // Extrapolate Ext-mul count to the recursion-shrunk inner-verify size (58 B256 muls ~= 7.6M rows
    // ~= 2^23). STARK verify muls grow ~linearly in log2(n) (queries x fold rounds), so fit a line.
    let m = samples.len() as f64;
    let (sx, sy) = (samples.iter().map(|p| p.0).sum::<f64>(), samples.iter().map(|p| p.1).sum::<f64>());
    let sxx = samples.iter().map(|p| p.0 * p.0).sum::<f64>();
    let sxy = samples.iter().map(|p| p.0 * p.1).sum::<f64>();
    let slope = (m * sxy - sx * sy) / (m * sxx - sx * sx);
    let intercept = (sy - slope * sx) / m;
    for &log2n in &[23.0f64] {
        let ext = (slope * log2n + intercept).max(0.0);
        let gas = (ext * BASE_PER_EXT as f64 * BASE_MUL_GAS as f64) as u64;
        let usd = gas as f64 * gwei * 1e-9 * eth;
        println!(
            "\nExtrapolated OUTER verify @ 2^{:.0} rows (~7.6M, the recursion-shrunk inner verify): {:.0} Ext muls",
            log2n, ext
        );
        println!("  -> field-op verify gas ~= {gas} gas (${usd:.2} @15gwei); binary equiv was 754M gas.");
    }

    // SOUNDNESS: a trace that breaks a constraint must NOT verify.
    let nt = 4096usize;
    let bad = build_trace(nt, true);
    let bad_proof = prove_one_sub_air_with_trace(&bad, nt, blowup, pi_hash, dsep, NUM_CONSTRAINTS, c_eval, params);
    let bad = verify_one_sub_air_with_trace(&bad_proof, nt, blowup, pi_hash, dsep, WIDTH, NUM_CONSTRAINTS, eval_per_row, params);
    assert!(bad.is_err(), "tampered trace must be REJECTED, got {bad:?}");
    println!("\ntamper (broke AND at one row): REJECTED. The whole AIR verifies via ONE merged");
    println!("composition polynomial (single commitment -> single-random-oracle QROM security).");
    println!("Extrapolate: one GF(2^128) mul ~= a 2^15-row instance; the recursion-shrunk 58-B256-mul");
    println!("inner verify ~= 58*4 = 232 such 2^15 instances (~7.6M rows), proved once off-chain.");
}
