# AWS blowup-vs-RAM benchmark sweep plan

Plan for measuring the **final, sound** proof sizes of the witness-binding +
public-input-bound signature/DNS AIRs on AWS, by sweeping the LDE blowup to
find the smallest sound proof that fits available RAM.

Script: [`scripts/bench-sweep-blowup.sh`](../scripts/bench-sweep-blowup.sh)
(RSS-capturing, portable macOS/Linux).

## Why a sweep (the blowup↔r coupling)

The FRI query count `r` is derived from the LDE blowup via the slack-Johnson
bound (`stark_level::num_queries_for_blowup`):

```
bits/query = −log2(√(1/blowup) · 1.05)      # η-included slack Johnson
r          = ⌈(λ + log2(M+2) + 1) / bits_per_query⌉   # λ=128/192/256, M=8
```

So `r` shrinks as blowup grows: **55/81/108 at blowup=32, but 309/457/606 at
blowup=2** (L1/L3/L5). Reusing the blowup=32 floor (55) at a lower blowup is
**under-secured** (≈24 bits at b=2, ≈51 at b=4) — the benches now compute the
correct `r` per blowup (commit that introduced `num_queries_for_blowup` in the
benches).

Two consequences:

1. **Proof size ≈ linear in `r`.** Higher blowup → fewer queries → smaller
   sound proof.
2. **Prover RAM grows with the LDE (∝ blowup)** but is partly offset because
   the `2·r·w` trace openings shrink as `r` falls. Net peak RSS grows
   **sublinearly** in blowup.

The sweep finds, per wide AIR, the **largest blowup whose measured peak RSS ≤
~80% of instance RAM** → the smallest publishable sound proof.

## M4 calibration (L1, measured — `/usr/bin/time -l`)

| AIR | blowup | r | proof | peak RSS | prove |
|-----|-------:|--:|------:|---------:|------:|
| RSA-2048   | 32 | 55  | 9.82 MiB | 0.26 GiB | 0.80 s |
| Ed25519    | 2  | 309 | 199 MiB  | 4.4 GiB  | 6.8 s  |
| Ed25519    | 4  | 143 | 92 MiB   | 6.9 GiB  | 9.6 s  |
| Ed25519    | 8  | 93  | 60 MiB   | 8.4 GiB  | 20.9 s |
| ECDSA-P256 | 4  | 143 | 430 MiB  | 8.4 GiB  | 31 s   |
| ML-DSA-44  | 4  | 143 | 23.9 MiB | 0.70 GiB | 3.5 s  |
| ML-DSA-44  | 8  | 93  | 20.9 MiB | 1.1 GiB  | 6.5 s  |
| DS-KSK     | 32 | —   | 1.24 MiB | 0.08 GiB | 0.20 s |

**Key observation:** Ed25519 RSS rose only 1.9× (4.4→8.4 GiB) for a 4× blowup
increase — sublinear, because `r` fell 309→93. RAM is **not** the binding
constraint; prove time (≈linear in blowup) is.

### Calibrated extrapolation to blowup=32

| AIR | peak RSS @ b=32 | sound proof @ b=32 (L1/L3/L5) | prove @ b=32 (≈) |
|-----|----------------:|------------------------------:|-----------------:|
| RSA-2048   | <0.3 GiB | 9.8 / 14.8 / 20.2 MiB  | ~1 s |
| Ed25519    | ~18–24 GiB | ~35 / ~52 / ~70 MiB  | ~100 s |
| ECDSA-P256 | ~28–32 GiB | ~165 / ~243 / ~324 MiB | ~5–6 min/proof |
| ML-DSA     | ~3.5 GiB | ~17 / ~32 / ~55 MiB    | ~30 s |
| NSEC3/DS   | <0.1 GiB | sub-MiB / 1.24 MiB     | <0.2 s |

At blowup=32 the sound sizes return to the original small headline figures —
now correctly parameterized (`r=55/81/108` is sound at b=32). M4's 16 GiB
forced the low blowups that inflated the figures to 199/430 MiB.

## Instance recommendation

- **Primary: `r8g.4xlarge`** — Graviton4 (ARM, matches the deployment
  narrative), 16 vCPU, **128 GiB**, DDR5-5600. Real peak need is ~32 GiB, so
  this has comfortable headroom and good memory bandwidth (the prover is
  bandwidth-bound). Spot pricing.
- **Sufficient: `r8g.2xlarge`** — 8 vCPU, **64 GiB** — fits b=32 for RAM, fewer
  cores (slower, but the workload saturates ~4 cores anyway).
- The 256 GiB+ instances are **not** needed.

Set `RAYON_NUM_THREADS = <vCPU>`.

## Sweep matrix

RSA, NSEC3, DS stay at blowup=32 (tiny). Sweep the three memory-bound AIRs:

| AIR | blowups | levels |
|-----|---------|--------|
| Ed25519    | 16, 32 (32 expected best) | L1/L3/L5 |
| ECDSA-P256 | 16, 32 (optionally 64)    | L1/L3/L5 |
| ML-DSA     | 16, 32                    | L1/L3/L5 |

`b=64` only buys ~16% size (ECDSA 165→138 MiB) for ~2× prove time — not worth
it. **Settle on b=32.**

## Run

```bash
export STARKDNS_DIR=/path/to/STARK-DNS   # (script also derives root from its own path)
cd "$STARKDNS_DIR"
LEVELS="1 3 5" \
BLOWUPS_ED="16 32" BLOWUPS_ECDSA="16 32" BLOWUPS_MLDSA="16 32" \
  ./scripts/bench-sweep-blowup.sh
# results -> scripts/results/sweep_blowup.txt   (RSS-annotated, with soundness checks)
```

Full b=16+32 × L1/L3/L5 ≈ ~1 hour (dominated by ECDSA b=32, ~5–6 min/proof).

## Decision criterion & deliverable

For each wide AIR, pick the **largest blowup with measured peak RSS ≤ ~80% of
instance RAM** (expected: b=32 for all on a 64–128 GiB box). The deliverable is
the RSS-annotated `(AIR × level)` table at the chosen blowup — the final sound
proof sizes for the paper, each with its `(blowup, r)` and the
honest-accept / tamper-reject / cross-signature-reject confirmations. These
feed the deferred paper projection-recompute.
