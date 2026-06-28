# Width-efficient witness-binding for P256 ECDSA — multi-row AIR adaptation plan

Status: planned (resumption point). Author: anonymised for review.

## Problem

The witness-binding construction (`deep_ali::sub_air_with_trace::{prove,verify}_one_sub_air_with_trace`)
commits the trace LDE by rows and opens the queried rows, so its cost is
`r × width` cells (r = FRI query count). This is **sound** — it is the same
construction that rejects tampering for RSA, Ed25519 and the NSEC3 chain — but
its memory/payload scales with the trace **width**.

Both current P256 ECDSA AIRs (`p256_ecdsa_air` (demo) and `p256_ecdsa_air_v2`)
are **single-row, fully-unrolled**: the entire 256-step double-scalar-mult
`u1·G + u2·Q` is packed into **one row** of tens of thousands of columns
(`p256_ecdsa_air.rs:559` — *"Pads the single-row witness to n_trace"*;
`make_trace_row` is `width × 1`). Binding such a layout opens
`r × (entire computation)` per query → multi-GB. The bench
`ecdsa_p256_bound_bench` confirms this: it is functionally correct (native
verify passes) but is SIGKILL'd at ~3 GB while building the openings.

**Key point:** this is a *layout artifact*, not ECDSA's inherent cost. The old
non-binding path looked "fast" (≈95 s) partly because it never materialises the
trace openings at all (only the width-collapsed `c_eval`); binding makes the
single-row layout's latent width visible.

Three width-efficiency shortcuts were each ruled out on soundness grounds:
- **Column-subset binding** — unsound for a densely-constrained AIR: every
  constraint-referenced column must be opened, and the ECDSA scalar-mult
  references essentially all columns. Skipping any is a soundness hole.
- **OOD-batched DEEP** — does not reduce the openings: `r × width` row openings
  are inherent to FRI-STARKs over a wide trace; OOD only *adds* per-column
  evaluations.
- **Recursion** — needs an already-*bound* inner. Recursing the non-binding
  ECDSA `c_eval` adds no witness-binding; ML-DSA only recurses cheaply because
  it is *already* internally decomposed into narrow sub-AIRs.

## The fix: a StarkWare-style multi-row P256 ECDSA AIR

StarkWare's ECDSA (Cairo `common/ec.cairo`, the Stone prover's Ec-op / ECDSA
builtin) is cheap because it is a **multi-row** AIR: **one scalar-mult step per
row**, a *narrow* width (the running accumulator point + the current scalar
bits + the doubled/added points), ~256 rows, with the accumulator threaded
**row-to-row by transition constraints**. This matches the paper's existing
soundness model (*"a boundary constraint at row K−1 binds the chain output"*).

A narrow trace makes `r × width` small, so the **existing, audited** binding
prover binds it cheaply *and* soundly — **no new cryptography**, just the right
layout. StarkWare's AIR is open-source and well-vetted, so it is a sound
starting point to adapt rather than design from scratch.

## Design

1. **Layout (`p256_ecdsa_multirow_air.rs`):** `n_trace` = power-of-two ≥ ~256
   (K steps + setup/teardown). Narrow per-row columns (P256 field limbs):
   running accumulator `(x, y[, z])`, the current `u1`/`u2` bits, the doubled
   point, the conditionally-added point, partial sums — order *hundreds* of
   cells/row instead of tens of thousands in one row.
2. **Transition constraints (row r → r+1):** the double-and-add step for
   `u1·G + u2·Q` (Shamir's trick, or two chains): `acc_{r+1} = 2·acc_r ⊕
   bit_r·base`, expressed as short-Weierstrass (a = −3) point doubling +
   conditional add over the P256 base field (mod p), reusing
   `p256_field`/`p256_group`/`p256_scalar`. Plus a row-0 bit-decomposition
   boundary and a final-row boundary `R.x ≡ r (mod n)` binding the output to the
   public signature `r`.
3. **Fill:** standard double-scalar-mult recording each step's accumulator.
4. **Soundness model:** identical to the RSA/Ed25519 in-circuit verifiers — the
   transition constraints verify the EC arithmetic in-circuit; the final
   boundary binds `R.x` to the verifier-derivable public `r`. Fully in-circuit;
   no native pre-check.

## Binding (reuse the audited construction)

The multi-row evaluator reads **`cur` and `nxt`** (transition) → it is **not**
local, so use the **full** binding (`prove/verify_one_sub_air_with_trace`, with
next-row openings), *not* the `_local` variant.

- Register the AIR as an `AirType` (e.g. `AirType::P256EcdsaMultirow`) so
  `deep_ali_merge_general` + `air_workloads::evaluate_constraints` +
  `build_execution_trace` work generically — the existing
  `airtype_bound_bench` then validates honest-accept / tampered-reject
  immediately.
- Expected memory: LDE `width × n_trace × blowup` + openings `r × width × 2`.
  For width ≈ a few hundred, `n_trace = 256`, `blowup = 32`: LDE ≈ a few million
  cells (tens of MB), openings ≈ 10⁵ cells — fits comfortably.

## Steps

1. Port StarkWare's EC double-and-add field-limb constraints (P256 short
   Weierstrass) into `p256_ecdsa_multirow_air.rs`, reusing the existing P256
   primitives.
2. Implement `build_layout` (narrow), `fill` (per-row double-scalar-mult),
   `eval_per_row` (transition double-and-add + `R.x ≡ r` boundary), constraint
   count.
3. Register `AirType::P256EcdsaMultirow` (merge + eval + trace builder dispatch).
4. Bound bench: real P256 signature (`p256` crate) → fill → the existing
   `prove/verify_one_sub_air_with_trace` → confirm honest accepts, a tampered
   step/`R.x` rejects in-circuit, and measure prove/verify/proof at L1/L3/L5.
5. Wire into the product ECDSA path (`se_zone_demo` / the epoch pipeline),
   replacing the single-row v2 verifier.

## Validation

- Honest signature accepts; a corrupted accumulator step or `R.x` rejects
  **in-circuit** (not via a native pre-check).
- Cross-check against `ecdsa_verify_native` on NIST/RFC test vectors.
- Re-measure the sound cost at NIST Levels 1/3/5 (slack `r = 55/81/108`,
  Fp⁶/Fp⁶/Fp⁸) and add to the paper's level table.

## References

- Current single-row AIRs: `crates/deep_ali/src/p256_ecdsa_air.rs`,
  `p256_ecdsa_air_v2.rs`; merges `deep_ali_merge_p256_ecdsa_streaming` /
  `_v2_streaming` / `_v2_rowgated_streaming` (lib.rs).
- Audited binding: `crates/deep_ali/src/sub_air_with_trace.rs`
  (`prove/verify_one_sub_air_with_trace[_local]`).
- StarkWare ECDSA: Cairo `common/ec.cairo`; Stone prover Ec-op / ECDSA builtin
  AIR (open-source) — the EC double-and-add constraints to adapt.
