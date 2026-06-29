#!/usr/bin/env bash
# bench-sweep-blowup.sh — sound proof-size vs. blowup sweep with peak-RSS capture.
#
# For the witness-binding + public-input-bound signature/DNS AIRs, the FRI
# query count r is derived from the LDE blowup via the slack-Johnson bound
# (stark_level::num_queries_for_blowup): higher blowup -> fewer queries ->
# SMALLER sound proof, but a larger LDE.  This script sweeps the blowup per
# AIR and records, for each (AIR, NIST level, blowup):
#     prove ms | verify ms | complete sound proof size | PEAK RSS
# plus the soundness checks (honest-accept / tamper-reject / cross-sig-reject)
# emitted by each bench.  Peak RSS is the AWS instance-sizing decision variable.
#
# Portable: macOS (`/usr/bin/time -l`, bytes) and Linux (`/usr/bin/time -v`,
# kbytes).  RSS grows SUBLINEARLY in blowup (as blowup rises, r falls, so the
# 2*r*w trace openings shrink and offset the LDE growth) — see
# docs/aws-benchmark-sweep-plan.md.
#
# Output:
#   scripts/results/sweep_blowup.txt   — RSS-annotated results (one line per run)
#   $TMPDIR/sweep_<tag>.log            — full per-run log (stdout+stderr+time)
#
# Usage (defaults shown):
#   LEVELS="1 3 5" \
#   BLOWUPS_ED="2 4 8 16 32" BLOWUPS_ECDSA="4 8 16 32" BLOWUPS_MLDSA="4 8 16 32" \
#   ./scripts/bench-sweep-blowup.sh
#
# Calibrated guidance (see plan): all AIRs fit blowup=32 in ~32 GiB peak RSS;
# the small sound sizes (RSA 9.8 / Ed25519 ~35 / ECDSA ~165 / ML-DSA ~17 MiB at
# L1) are recovered at blowup=32.  An r8g.4xlarge (128 GiB) has ample headroom.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p scripts/results
OUT="${OUT:-scripts/results/sweep_blowup.txt}"
TMPD="${TMPDIR:-/tmp}"
: > "$OUT"
log() { printf '%s\n' "$*" | tee -a "$OUT"; }

OS="$(uname)"
peak_rss_mib() { # $1 = time-wrapped log file -> echo peak RSS in MiB (or '?')
  local f="$1" v
  if [ "$OS" = "Darwin" ]; then          # /usr/bin/time -l : bytes
    v=$(grep -a "maximum resident set size" "$f" | awk '{print $1}' | tail -1)
    [ -n "$v" ] && echo $(( v / 1048576 )) || echo '?'
  else                                    # GNU time -v : kbytes
    v=$(grep -a "Maximum resident set size" "$f" | awk -F': ' '{print $2}' | tail -1)
    [ -n "$v" ] && echo $(( v / 1024 )) || echo '?'
  fi
}
time_wrap() { # peak-RSS timer wrapper (platform-specific flag)
  if [ "$OS" = "Darwin" ]; then /usr/bin/time -l "$@"; else /usr/bin/time -v "$@"; fi
}
run_one() { # $1=tag  $2=result-grep-ERE  $3...=shell command
  local tag="$1" pat="$2"; shift 2
  local lf="$TMPD/sweep_$(echo "$tag" | tr -c 'A-Za-z0-9' '_').log"
  time_wrap bash -c "$*" > "$lf" 2>&1
  local rss res
  rss=$(peak_rss_mib "$lf")
  res=$(grep -aE "$pat" "$lf" | grep -avE "QUERY POS" | tail -1)
  log "[$tag]  RSS=${rss} MiB | ${res:-<no result line — see $lf>}"
}

LEVELS="${LEVELS:-1 3 5}"
BLOWUPS_ED="${BLOWUPS_ED:-2 4 8 16 32}"
BLOWUPS_ECDSA="${BLOWUPS_ECDSA:-4 8 16 32}"
BLOWUPS_MLDSA="${BLOWUPS_MLDSA:-4 8 16 32}"

for LVL in $LEVELS; do
  case "$LVL" in
    1) FEAT="sha3-256,mldsa-44,parallel" ;;
    3) FEAT="sha3-384,mldsa-65,parallel" ;;
    5) FEAT="sha3-512,mldsa-87,parallel" ;;
    *) echo "unknown level $LVL"; exit 1 ;;
  esac
  log "==================== NIST L$LVL ($FEAT) ===================="

  # RSA-2048 — fixed blowup=32 (tiny trace), single point.
  run_one "RSA_L${LVL}_b32" "level=L.*proof_mib=" \
    "cargo run --release -q -p deep_ali --example rsa2048_exp_bound_bench --no-default-features --features '$FEAT'"

  for B in $BLOWUPS_ED; do
    run_one "Ed25519_L${LVL}_b${B}" "ed25519_verify_bound.*proof_mib=" \
      "BENCH_BLOWUP=$B cargo run --release -q -p deep_ali --example ed25519_verify_bound_bench --no-default-features --features '$FEAT'"
  done

  for B in $BLOWUPS_ECDSA; do
    run_one "ECDSA_L${LVL}_b${B}" "ecdsa_verify_multirow_pub.*proof_mib=" \
      "BENCH_BLOWUP=$B cargo run --release -q -p deep_ali --example ecdsa_verify_multirow_bound_bench --no-default-features --features '$FEAT'"
  done

  for B in $BLOWUPS_MLDSA; do
    run_one "MLDSA_L${LVL}_b${B}" "v2_bench level=.*proof_kib=" \
      "BENCH_BLOWUP=$B cargo test --release -q -p deep_ali --no-default-features --features '$FEAT' v2_bench -- --include-ignored --nocapture --test-threads=1"
  done

  # DNS denial-of-existence + DS digest (fixed blowup=32, tiny).
  run_one "DNS_L${LVL}_b32" "NODATA|OPTOUT|WILDCARD|DS-KSK" \
    "cargo test --release -q -p swarm-dns --no-default-features --features '$FEAT' nsec3_bound_sizes_l1 -- --ignored --nocapture --test-threads=1"
done
log "==================== DONE ===================="
echo "results -> $OUT"
