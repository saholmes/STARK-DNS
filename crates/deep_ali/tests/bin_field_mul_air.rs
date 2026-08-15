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

/// STIR folding schedule with NON-BINARY folding factor 16 (= 2^4): STIR folds by 2^x per round with
/// no soundness loss, so log2(n0) binary rounds collapse to ~log16(n0) rounds -- the source of STIR's
/// smaller proof and cheaper verify. Returns (schedule, final_size).
fn stir_schedule(n0: usize) -> (Vec<usize>, usize) {
    let mut sched = Vec::new();
    let mut cur = n0;
    while cur % 16 == 0 && cur / 16 >= 2 {
        sched.push(16);
        cur /= 16;
    }
    while cur % 2 == 0 && cur / 2 >= 2 {
        sched.push(2);
        cur /= 2;
    }
    (sched, cur)
}

/// NIST L1 UNCONDITIONAL STIR: r=54 queries at rho=1/32 (blowup 32), Johnson bound (2.5 bits/query,
/// unconditional). Non-binary (16-ary) folding => few rounds => smaller proof + cheaper verify.
fn stir_params(n0: usize, pi_hash: [u8; 32]) -> DeepFriParams {
    let (schedule, final_size) = stir_schedule(n0);
    let d_final = (final_size / 2).max(1);
    let mut p = DeepFriParams::new(schedule, 0, 42).with_stir().with_s0(54).with_d_final(d_final);
    p.public_inputs_hash = Some(pi_hash);
    p
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

use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain};
use deep_ali::cubic_ext::{CubeExt, GoldilocksCubeConfig};
use deep_ali::fri::{deep_fri_proof_size_bytes, deep_fri_prove, deep_fri_verify, FriDomain};

type E = CubeExt<GoldilocksCubeConfig>;

/// STIR (non-binary 16-ary folding) vs FRI (binary) on the same low-degree test (the outer proof's
/// FRI layer), same rate rho=1/32 and same r=54 queries (NIST L1 unconditional Johnson). STIR folds
/// by 16 per round -> log16(n) rounds vs FRI's log2(n) -> smaller proof + fewer verify Ext muls.
#[test]
fn stir_vs_fri_outer_proof() {
    use rand::{rngs::StdRng, SeedableRng};
    let mut rng = StdRng::seed_from_u64(1234);
    const BASE_MUL_GAS: u64 = 91;
    const BASE_PER_EXT: u64 = 9;

    println!("\n===== STIR (16-ary folding) vs FRI (binary) on the outer proof, r=54 rho=1/32 (NIST L1) =====");
    println!("{:>7} {:>6} {:>7} {:>7} | {:>10} {:>10} {:>7} | {:>10} {:>10} {:>7}",
        "log2 n", "rounds", "FRI B", "STIR B", "FRIextmul", "STIRextmul", "mul x", "FRI gas", "STIR gas", "gas x");

    for &logn in &[14usize, 16, 18] {
        let n = 1usize << logn;
        let degree = n / 32 - 1; // rho = 1/32
        let dom = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let poly = DensePolynomial::<F>::rand(degree, &mut rng);
        let evals: Vec<F> = dom.fft(&poly.coeffs);
        let domain0 = FriDomain::new_radix2(n);

        // Fold both to final_size = 16 (d_final = 8). FRI: binary (logn-4 rounds). STIR: 16-ary.
        let fri_sched = vec![2usize; logn - 4];
        let mut stir_sched = Vec::new();
        let mut cur = n;
        while cur % 16 == 0 && cur / 16 >= 16 {
            stir_sched.push(16);
            cur /= 16;
        }
        while cur / 2 >= 16 {
            stir_sched.push(2);
            cur /= 2;
        }
        let d_final = (cur / 2).max(1);

        let fri_p = DeepFriParams::new(fri_sched, 54, 42).with_d_final(8);
        let stir_p = DeepFriParams::new(stir_sched.clone(), 0, 42).with_stir().with_s0(54).with_d_final(d_final);

        let fri_proof = deep_fri_prove::<E>(evals.clone(), domain0, &fri_p);
        let stir_proof = deep_fri_prove::<E>(evals.clone(), domain0, &stir_p);
        let fri_b = deep_fri_proof_size_bytes(&fri_proof, false);
        let stir_b = deep_fri_proof_size_bytes(&stir_proof, true);

        reset_ext_mul_count();
        assert!(deep_fri_verify(&fri_p, &fri_proof), "FRI verify @ 2^{logn}");
        let fri_ext = ext_mul_count();
        reset_ext_mul_count();
        assert!(deep_fri_verify(&stir_p, &stir_proof), "STIR verify @ 2^{logn}");
        let stir_ext = ext_mul_count();

        let fri_gas = fri_ext * BASE_PER_EXT * BASE_MUL_GAS + fri_b as u64 * 16;
        let stir_gas = stir_ext * BASE_PER_EXT * BASE_MUL_GAS + stir_b as u64 * 16;
        println!("{logn:>7} {:>6} {fri_b:>7} {stir_b:>7} | {fri_ext:>10} {stir_ext:>10} {:>6.1}x | {fri_gas:>8} g {stir_gas:>8} g {:>5.1}x",
            stir_sched.len(), fri_ext as f64 / stir_ext.max(1) as f64, fri_gas as f64 / stir_gas.max(1) as f64);
    }
    println!("\nSTIR's non-binary (16-ary) folding cuts rounds (log16 vs log2), shrinking BOTH the proof");
    println!("(calldata/chunking bottleneck) AND the verify Ext-mul count -- unconditionally (Johnson).");
}
