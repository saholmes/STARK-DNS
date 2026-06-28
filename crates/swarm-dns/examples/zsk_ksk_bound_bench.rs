//! WITNESS-BINDING ZSK→KSK Ed25519 bound bench (zone-pipeline path).
//!
//! Exercises the newly-wired `prove_zsk_ksk_binding_v2_bound` /
//! `verify_zsk_ksk_binding_v2_bound` — the SOUND replacement for
//! `prove_zsk_ksk_binding_v2` (a bare `deep_fri_prove(c_eval)` proof,
//! sound only for low-degreeness).  The bound path routes the v16
//! Ed25519 verify AIR through `prove_one_sub_air_with_trace`: the trace
//! LDE is Merkle-committed (root folded into the FS pi_hash) and the
//! verifier re-checks `c_eval(x)·Z_H(x) = Σ αⱼ Φⱼ(trace[x])` at every
//! authenticated query opening.
//!
//! Reports the COMPLETE sound proof size (serialize_proof: LDT + 2rw
//! trace openings + Merkle paths).  Demonstrates honest-ACCEPT,
//! cross-signature REJECT, and (constraint-level witness-binding is the
//! property of `prove_one_sub_air_with_trace`, independently proven by
//! `deep_ali --example ed25519_verify_bound_bench`, flip-cell → reject).
//!
//! Run (L1): cargo run --release -p swarm-dns --example zsk_ksk_bound_bench
//!  (L3: --no-default-features --features sha3-384,mldsa-65,parallel ;
//!   L5: --no-default-features --features sha3-512,mldsa-87,parallel)

use std::time::Instant;

use ed25519_dalek::{Signer, SigningKey};
use rand::rngs::StdRng;
use rand::SeedableRng;

use swarm_dns::prover::{
    prove_zsk_ksk_binding_v2_bound, slack_johnson_queries,
    verify_zsk_ksk_binding_v2_bound, LdtMode,
};

const FS_BINDING:  [u8; 32] = [0xCA; 32];
const MERKLE_ROOT: [u8; 32] = [0x77; 32];

/// Deterministic real Ed25519 (pubkey, signature, message) triple.
fn synth_ed25519(seed: u64, msg: &[u8]) -> ([u8; 32], [u8; 64]) {
    let mut rng = StdRng::seed_from_u64(seed);
    let sk = SigningKey::generate(&mut rng);
    let pk = sk.verifying_key().to_bytes();
    let sig = sk.sign(msg).to_bytes();
    (pk, sig)
}

fn main() {
    eprintln!("=== zsk_ksk_bound_bench: WITNESS-BINDING Ed25519 ZSK→KSK (zone path) ===");
    let ldt = LdtMode::Fri; // deployed constant-rate ÷4 fold
    let k_scalar = 256usize;

    // Target soundness bits for the active NIST level (feature-gated).
    let level = deep_ali::stark_level::NIST_LEVEL;
    let lambda = match level { 1 => 128, 3 => 192, _ => 256 };

    // blowup ↔ r are COUPLED via the slack Johnson bound: per-query yield
    // is 0.5·log2(blowup) − 0.07, so a smaller blowup needs MORE queries.
    // blowup=32 (r≈55 at L1) is the compact config but its LDE for the wide
    // v16 trace is ~21 GB; BENCH_BLOWUP picks a memory-feasible point and r
    // is recomputed to stay sound.  Default 8 (~5 GB LDE) — bigger proof,
    // fits commodity RAM.
    let blowup: usize = std::env::var("BENCH_BLOWUP").ok()
        .and_then(|s| s.parse().ok()).unwrap_or(8usize);
    let r = slack_johnson_queries(blowup, lambda);
    let bits_per_q = 0.5_f64 * (blowup as f64).log2() - 0.07_f64;
    eprintln!("    NIST L{level} (λ={lambda}) | blowup={blowup} (rate 1/{blowup}) | \
               slack Johnson {bits_per_q:.2} b/q × r={r} = {:.0} bits", bits_per_q * r as f64);

    // ── Signature A (the one we prove) and B (a distinct cross test). ──
    let msg_a: &[u8] = b"DNSKEY ZSK RRset, zone example., RFC 8080 alg 15";
    let msg_b: &[u8] = b"DNSKEY ZSK RRset, attacker.example., RFC 8080 alg 15";
    let (pk_a, sig_a) = synth_ed25519(0xA11CE, msg_a);
    let (pk_b, sig_b) = synth_ed25519(0xB0B, msg_b);

    // ── Prove A through the BOUND path. ──
    let t = Instant::now();
    let out = prove_zsk_ksk_binding_v2_bound(
        &pk_a, &sig_a, msg_a, &FS_BINDING, &MERKLE_ROOT, k_scalar, blowup, r, ldt,
    );
    let prove_ms = t.elapsed().as_secs_f64() * 1e3;
    let full_mib =
        deep_ali::sub_air_with_trace::serialize_proof(&out.proof).len() as f64 / 1048576.0;
    let fri_kib = out.proof.fri_proof_bytes.len() as f64 / 1024.0;

    // ── Honest verify → must ACCEPT. ──
    let t = Instant::now();
    let honest = verify_zsk_ksk_binding_v2_bound(
        &out, &pk_a, &sig_a, msg_a, &FS_BINDING, &MERKLE_ROOT, ldt,
    );
    let verify_ms = t.elapsed().as_secs_f64() * 1e3;
    eprintln!("[honest]   prove {prove_ms:.1} ms, verify {verify_ms:.2} ms, \
               fri {fri_kib:.1} KiB, FULL sound proof {full_mib:.2} MiB -> verify={honest}");
    assert!(honest, "BINDING BROKEN: honest ZSK→KSK signature must verify");

    // ── Cross-signature: A's proof under B's public inputs → must REJECT. ──
    let cross = verify_zsk_ksk_binding_v2_bound(
        &out, &pk_b, &sig_b, msg_b, &FS_BINDING, &MERKLE_ROOT, ldt,
    );
    eprintln!("[cross]    A's proof under B's (pk,sig,msg) -> verify={cross}");

    println!(
        "zsk_ksk_v2_bound L{level} k_scalar={k_scalar} n_trace={} blowup={blowup} r={r} \
         prove_ms={prove_ms:.1} verify_ms={verify_ms:.2} fri_kib={fri_kib:.1} \
         proof_mib={full_mib:.2} honest_verify={honest} cross_verify={cross}", out.n_trace);
    if honest && !cross {
        println!("=> SOUND PATH WIRED: honest accepts, cross-signature REJECTS \
                  (zone ZSK→KSK now witness-binding)");
    } else {
        println!("=> FAILED: honest={honest} cross={cross}");
    }
}
