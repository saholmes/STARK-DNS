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
    DeepFriParams {
        schedule: vec![2usize; n0.trailing_zeros() as usize],
        r: 16,
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

    println!("\n===== Binary-field multiply AIR over deep_ali (Goldilocks; DEEP-ALI merge -> 1 poly) =====");
    println!("width={WIDTH}, constraints={NUM_CONSTRAINTS} (degree 2). One GF(2^128) mul ~= 2^15 rows.");
    println!("{:>8} {:>10} {:>10} {:>12}", "n_trace", "prove ms", "verify ms", "proof bytes");

    for &nt in &[4096usize, 16384, 32768] {
        // Honest trace: prove + verify.
        let trace = build_trace(nt, false);
        let t0 = Instant::now();
        let proof = prove_one_sub_air_with_trace(&trace, nt, blowup, pi_hash, dsep, NUM_CONSTRAINTS, c_eval, params);
        let prove_ms = t0.elapsed().as_secs_f64() * 1e3;
        let bytes = serialize_proof(&proof).len();
        let t1 = Instant::now();
        let ok = verify_one_sub_air_with_trace(&proof, nt, blowup, pi_hash, dsep, WIDTH, NUM_CONSTRAINTS, eval_per_row, params);
        let verify_ms = t1.elapsed().as_secs_f64() * 1e3;
        assert!(ok.is_ok(), "honest AIR must verify at n_trace={nt}: {ok:?}");
        println!("{nt:>8} {prove_ms:>10.1} {verify_ms:>10.2} {bytes:>12}");
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
