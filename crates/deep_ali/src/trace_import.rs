// src/trace_import.rs

use ark_goldilocks::Goldilocks as F;
use ark_ff::Zero;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain as Domain};
use std::collections::HashMap;

/// Four evaluation vectors over the FRI domain, derived from a real
/// execution trace rather than random sampling.
///
/// Each vector contains evaluations of a polynomial with degree < n0/rate_inv
/// on an n0-point domain.  This gives FRI the same algebraic structure
/// (bounded-degree polynomials evaluated on a larger domain) as a real STARK,
/// which random vectors do NOT have — random vectors are full-rank, so they
/// don't exercise the degree-testing logic that FRI actually performs.
pub struct RealTraceInputs {
    pub a_eval: Vec<F>,
    pub s_eval: Vec<F>,
    pub e_eval: Vec<F>,
    pub t_eval: Vec<F>,
}

/// Build 4 Fibonacci trace columns over Goldilocks, interpolate, and
/// LDE-evaluate on an n0-point domain.
///
/// `n0`:       FRI domain size (power of 2)
/// `rate_inv`: blowup factor (2, 4, 8, 16…).  Polynomials have degree < n0/rate_inv.
pub fn real_trace_inputs(n0: usize, rate_inv: usize) -> RealTraceInputs {
    assert!(n0.is_power_of_two());
    assert!(rate_inv >= 2 && rate_inv.is_power_of_two());
    let trace_len = n0 / rate_inv;
    assert!(trace_len >= 2, "trace too short");

    // Four Fibonacci-like columns with different seeds.
    // Each column satisfies col[i] = col[i-1] + col[i-2],
    // so it's a valid execution trace of a degree-1 transition constraint.
    let seeds: [(u64, u64); 4] = [
        (1, 1),
        (2, 3),
        (5, 8),
        (13, 21),
    ];

    let trace_dom = Domain::<F>::new(trace_len).unwrap();
    let lde_dom   = Domain::<F>::new(n0).unwrap();

    let mut evals = Vec::with_capacity(4);

    for &(s0, s1) in &seeds {
        // 1. Build trace column
        let mut col = Vec::with_capacity(trace_len);
        col.push(F::from(s0));
        col.push(F::from(s1));
        for i in 2..trace_len {
            col.push(col[i - 1] + col[i - 2]);
        }

        // 2. Interpolate: IFFT over trace domain → coefficients
        //    Polynomial has degree trace_len - 1 = n0/rate_inv - 1
        let coeffs = trace_dom.ifft(&col);

        // 3. LDE: pad coefficients to n0 (zeros for high degrees),
        //    then FFT over the larger domain
        let mut padded = coeffs;
        padded.resize(n0, F::zero());
        evals.push(lde_dom.fft(&padded));
    }

    RealTraceInputs {
        a_eval: evals.remove(0),
        s_eval: evals.remove(0),
        e_eval: evals.remove(0),
        t_eval: evals.remove(0),
    }
}

/// Convert an arbitrary execution trace (produced by an AIR workload)
/// into `RealTraceInputs` by interpolating each column and LDE-evaluating
/// on the extended domain.
///
/// `trace_columns`: each inner Vec is one column of length `n0 / blowup`.
/// `n0`:            FRI / extended-evaluation domain size (power of 2).
/// `blowup`:        rate inverse (typically 4).
///
/// The function maps the first four columns to `a_eval … t_eval`.
/// If the trace has fewer than four columns, columns are reused with
/// wraparound (same strategy as `import_winterfell_trace`).
/// If the trace has more than four columns, the extra columns are ignored.
pub fn trace_inputs_from_air(
    trace_columns: Vec<Vec<F>>,
    n0: usize,
    blowup: usize,
) -> RealTraceInputs {
    let num_cols = trace_columns.len();
    assert!(num_cols >= 1, "need at least 1 trace column");
    assert!(n0.is_power_of_two());
    assert!(blowup >= 2 && blowup.is_power_of_two());

    let trace_len = n0 / blowup;
    assert!(trace_len >= 2, "trace too short");

    // Sanity-check that every column has the expected length
    for (i, col) in trace_columns.iter().enumerate() {
        assert_eq!(
            col.len(),
            trace_len,
            "column {} has length {} but expected {}",
            i,
            col.len(),
            trace_len
        );
    }

    let trace_dom = Domain::<F>::new(trace_len).unwrap();
    let lde_dom   = Domain::<F>::new(n0).unwrap();

    let lde = |col: &[F]| -> Vec<F> {
        let coeffs = trace_dom.ifft(col);
        let mut padded = coeffs;
        padded.resize(n0, F::zero());
        lde_dom.fft(&padded)
    };

    // Map columns to the four required vectors with wraparound
    let a_eval = lde(&trace_columns[0]);
    let s_eval = lde(&trace_columns[1 % num_cols]);
    let e_eval = lde(&trace_columns[2 % num_cols]);
    let t_eval = lde(&trace_columns[3 % num_cols]);

    RealTraceInputs { a_eval, s_eval, e_eval, t_eval }
}

/// Same as above but reads trace columns from a binary file exported
/// by Winterfell's FibSmall example (f64 = Goldilocks).
///
/// File format (produced by the export binary in Path B):
///   Header line:  "TRACE <trace_len> <num_cols> <field_bits>\n"
///   Body:         column-major, each element as u64 little-endian (8 bytes)
pub fn import_winterfell_trace(path: &str, n0: usize) -> RealTraceInputs {
    use std::io::{BufRead, BufReader, Read};
    use std::fs::File;

    let file = File::open(path).expect("cannot open trace file");
    let mut reader = BufReader::new(file);

    // Parse header
    let mut header = String::new();
    reader.read_line(&mut header).unwrap();
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    assert_eq!(parts[0], "TRACE");
    let trace_len: usize = parts[1].parse().unwrap();
    let num_cols: usize  = parts[2].parse().unwrap();
    assert!(num_cols >= 2, "need at least 2 trace columns");

    // Read columns (u64 LE → Goldilocks)
    let mut columns: Vec<Vec<F>> = Vec::with_capacity(num_cols);
    let mut buf = [0u8; 8];

    for _col in 0..num_cols {
        let mut column = Vec::with_capacity(trace_len);
        for _row in 0..trace_len {
            reader.read_exact(&mut buf).unwrap();
            let val = u64::from_le_bytes(buf);
            column.push(F::from(val));
        }
        columns.push(column);
    }

    // Interpolate and LDE, same as above
    let trace_dom = Domain::<F>::new(trace_len).unwrap();
    let lde_dom   = Domain::<F>::new(n0).unwrap();

    let lde = |col: &[F]| -> Vec<F> {
        let coeffs = trace_dom.ifft(col);
        let mut padded = coeffs;
        padded.resize(n0, F::zero());
        lde_dom.fft(&padded)
    };

    // Map columns to the four vectors.
    // With 2 trace columns we duplicate; with 4+ we use the first 4.
    let a_eval = lde(&columns[0]);
    let s_eval = lde(&columns[1 % num_cols]);
    let e_eval = lde(&columns[2 % num_cols]);
    let t_eval = lde(&columns[3 % num_cols]);

    RealTraceInputs { a_eval, s_eval, e_eval, t_eval }
}
// ═══════════════════════════════════════════════════════════════════
//  StarkWare column-major JSON trace import
// ═══════════════════════════════════════════════════════════════════

/// Import a StarkWare column-major JSON string.
///
/// Column values are u64 Goldilocks field elements.
/// Column ordering follows `column_order` when provided, otherwise sorted by name.
/// Returns raw trace columns as `Vec<Vec<F>>` (not LDE-evaluated).
/// Call `lde_trace_columns` on the result before `deep_ali_merge_general`.
pub fn import_starkware_json(
    json_str: &str,
    column_order: Option<&[&str]>,
) -> Result<Vec<Vec<F>>, String> {
    let v: serde_json::Value = serde_json::from_str(json_str)
        .map_err(|e| format!("JSON parse error: {e}"))?;

    let length = v["length"].as_u64()
        .ok_or("missing 'length' field")? as usize;
    let cols_val = v["columns"].as_object()
        .ok_or("'columns' must be an object")?;

    let extract = |name: &str| -> Result<Vec<F>, String> {
        let arr = cols_val.get(name)
            .ok_or_else(|| format!("column '{name}' not found"))?
            .as_array()
            .ok_or_else(|| format!("column '{name}' must be an array"))?;
        if arr.len() != length {
            return Err(format!("column '{name}': expected {length} rows, got {}", arr.len()));
        }
        arr.iter().map(|v| {
            v.as_u64()
                .ok_or_else(|| format!("column '{name}' has non-u64 value: {v}"))
                .map(F::from)
        }).collect()
    };

    let names: Vec<String> = if let Some(order) = column_order {
        order.iter().map(|s| s.to_string()).collect()
    } else {
        let mut names: Vec<String> = cols_val.keys().cloned().collect();
        names.sort();
        names
    };

    names.iter().map(|name| extract(name)).collect()
}

/// Import a StarkWare JSON trace file from disk.
pub fn import_starkware_json_file(
    path: &str,
    column_order: Option<&[&str]>,
) -> Result<Vec<Vec<F>>, String> {
    let s = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read '{path}': {e}"))?;
    import_starkware_json(&s, column_order)
}

/// LDE-evaluate all columns of a raw trace.
///
/// Each column must have `trace_len` elements.  Returns columns evaluated
/// on an `n0 = trace_len * blowup`-point domain, ready for `deep_ali_merge_general`.
pub fn lde_trace_columns(
    columns: &[Vec<F>],
    trace_len: usize,
    blowup: usize,
) -> Result<Vec<Vec<F>>, String> {
    if columns.is_empty() {
        return Err("no columns provided".into());
    }
    for (i, col) in columns.iter().enumerate() {
        if col.len() != trace_len {
            return Err(format!("column {i}: expected {trace_len} rows, got {}", col.len()));
        }
    }
    if !trace_len.is_power_of_two() {
        return Err(format!("trace_len {trace_len} must be a power of 2"));
    }
    if blowup < 2 || !blowup.is_power_of_two() {
        return Err(format!("blowup {blowup} must be a power-of-2 >= 2"));
    }

    let n0 = trace_len * blowup;
    let trace_dom = Domain::<F>::new(trace_len).unwrap();
    let lde_dom = Domain::<F>::new(n0).unwrap();

    columns.iter().map(|col| {
        let coeffs = trace_dom.ifft(col);
        let mut padded = coeffs;
        padded.resize(n0, F::zero());
        Ok(lde_dom.fft(&padded))
    }).collect()
}

/// Deterministic splitmix64 PRG (no external rand dependency).
/// Used only to draw the ZK trace-mask coefficients.
#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Zero-knowledge (hidden-trace) variant of [`lde_trace_columns`].
///
/// For every column index in `mask_cols`, a random low-degree multiple of the
/// trace-domain vanishing polynomial `Z_H(X) = X^{trace_len} - 1` is added to
/// the column's interpolant **before** the LDE:
///
/// ```text
///   p_masked(X) = p(X) + r(X) * Z_H(X),   deg(r) < mask_deg
/// ```
///
/// `Z_H` is zero on every trace row (the order-`trace_len` subgroup), so
/// `p_masked` agrees with `p` on the trace and the AIR is unaffected; off the
/// trace domain `p_masked` is randomized, so the LDE values opened by FRI reveal
/// nothing about the hidden column beyond what the (public) constraints force.
///
/// Soundness for a LINEAR constraint `C` is preserved: `C(p + r*Z_H) =
/// C(p) + r'*Z_H` is still divisible by `Z_H`, so the quotient stays
/// low-degree. (For non-linear constraints the degree budget must account for
/// the mask; this helper is intended for the linear conservation AIR.)
///
/// Hiding is statistical and governed by `mask_deg`: with `mask_deg` at least
/// the number of FRI query openings, every opened point of a masked column is
/// blinded by an independent uniform term. `seed` selects the mask; distinct
/// seeds give independent masks (used by the two-witness hiding gate).
///
/// This is purely additive — existing callers of [`lde_trace_columns`] are
/// unchanged (it is the `mask_cols = []` case, bit-for-bit).
pub fn lde_trace_columns_masked(
    columns: &[Vec<F>],
    trace_len: usize,
    blowup: usize,
    mask_cols: &[usize],
    mask_deg: usize,
    seed: u64,
) -> Result<Vec<Vec<F>>, String> {
    if columns.is_empty() {
        return Err("no columns provided".into());
    }
    for (i, col) in columns.iter().enumerate() {
        if col.len() != trace_len {
            return Err(format!("column {i}: expected {trace_len} rows, got {}", col.len()));
        }
    }
    if !trace_len.is_power_of_two() {
        return Err(format!("trace_len {trace_len} must be a power of 2"));
    }
    if blowup < 2 || !blowup.is_power_of_two() {
        return Err(format!("blowup {blowup} must be a power-of-2 >= 2"));
    }
    for &c in mask_cols {
        if c >= columns.len() {
            return Err(format!("mask column {c} out of range ({} columns)", columns.len()));
        }
    }

    let n0 = trace_len * blowup;
    // The mask r*Z_H occupies coefficients [0, trace_len + mask_deg); it must fit
    // in the n0-length LDE coefficient buffer or the high terms would alias.
    let eff_mask_deg = mask_deg.min(n0.saturating_sub(trace_len));
    let trace_dom = Domain::<F>::new(trace_len).unwrap();
    let lde_dom = Domain::<F>::new(n0).unwrap();

    columns.iter().enumerate().map(|(ci, col)| {
        let coeffs = trace_dom.ifft(col);
        let mut padded = coeffs;
        padded.resize(n0, F::zero());
        if mask_cols.contains(&ci) {
            // Draw a per-column mask stream, domain-separated by column index.
            let mut st = seed ^ ((ci as u64).wrapping_mul(0xD1B5_4A32_D192_ED03));
            for i in 0..eff_mask_deg {
                let r = F::from(splitmix64(&mut st));
                // r * Z_H = r * (X^{trace_len} - 1):  +r at X^{trace_len+i}, -r at X^i.
                padded[i] -= r;
                padded[trace_len + i] += r;
            }
        }
        Ok(lde_dom.fft(&padded))
    }).collect()
}

/// Coset variant of [`lde_trace_columns_masked`] — evaluates the masked interpolants on the COSET
/// `coset_offset * H_{n0}` instead of the subgroup `H_{n0}`. With `coset_offset` a multiplicative
/// generator the coset is disjoint from the trace subgroup, so no committed/opened LDE point is ever a
/// trace row and `Z_H` never vanishes on it. Pair with a GLOBAL (ungated) coset merge + coset verify
/// using the same `coset_offset`: hiding only holds when the composition is a genuine global polynomial
/// (see `docs/goldilocks-confidential-amounts-scope.md` §4ter — a per-row-GATED composition breaks on a
/// coset because the gating no longer controls the interpolant on the trace subgroup).
#[allow(clippy::too_many_arguments)]
pub fn lde_trace_columns_masked_coset(
    columns: &[Vec<F>],
    trace_len: usize,
    blowup: usize,
    mask_cols: &[usize],
    mask_deg: usize,
    seed: u64,
    coset_offset: F,
) -> Result<Vec<Vec<F>>, String> {
    if columns.is_empty() {
        return Err("no columns provided".into());
    }
    for (i, col) in columns.iter().enumerate() {
        if col.len() != trace_len {
            return Err(format!("column {i}: expected {trace_len} rows, got {}", col.len()));
        }
    }
    if !trace_len.is_power_of_two() {
        return Err(format!("trace_len {trace_len} must be a power of 2"));
    }
    if blowup < 2 || !blowup.is_power_of_two() {
        return Err(format!("blowup {blowup} must be a power-of-2 >= 2"));
    }
    for &c in mask_cols {
        if c >= columns.len() {
            return Err(format!("mask column {c} out of range ({} columns)", columns.len()));
        }
    }

    let n0 = trace_len * blowup;
    let eff_mask_deg = mask_deg.min(n0.saturating_sub(trace_len));
    let trace_dom = Domain::<F>::new(trace_len).unwrap();
    let lde_coset = Domain::<F>::new(n0)
        .unwrap()
        .get_coset(coset_offset)
        .ok_or("could not form LDE coset")?;

    columns.iter().enumerate().map(|(ci, col)| {
        let coeffs = trace_dom.ifft(col);
        let mut padded = coeffs;
        padded.resize(n0, F::zero());
        if mask_cols.contains(&ci) {
            let mut st = seed ^ ((ci as u64).wrapping_mul(0xD1B5_4A32_D192_ED03));
            for i in 0..eff_mask_deg {
                let r = F::from(splitmix64(&mut st));
                padded[i] -= r;
                padded[trace_len + i] += r;
            }
        }
        Ok(lde_coset.fft(&padded))
    }).collect()
}

#[cfg(test)]
mod zk_mask_tests {
    use super::*;
    use ark_ff::Field;

    #[test]
    fn masked_lde_vanishes_on_trace_rows_and_matches_default_when_unmasked() {
        let trace_len = 8usize;
        let blowup = 4usize;
        let col0: Vec<F> = (0..trace_len).map(|i| F::from((100 + i) as u64)).collect();
        let col1: Vec<F> = (0..trace_len).map(|i| F::from((7 * i + 3) as u64)).collect();
        let columns = vec![col0.clone(), col1.clone()];

        let plain = lde_trace_columns(&columns, trace_len, blowup).unwrap();
        let masked = lde_trace_columns_masked(&columns, trace_len, blowup, &[0], 6, 0xABCD_1234).unwrap();

        for k in 0..trace_len {
            assert_eq!(masked[0][k * blowup], col0[k], "masked col0 altered trace row {k}");
        }
        let mut differs = false;
        for j in 0..(trace_len * blowup) {
            if j % blowup != 0 && masked[0][j] != plain[0][j] {
                differs = true;
                break;
            }
        }
        assert!(differs, "mask had no off-domain effect (not hiding)");
        assert_eq!(masked[1], plain[1], "unmasked column diverged from default LDE");
    }

    #[test]
    fn distinct_seeds_give_distinct_masks_same_trace() {
        let trace_len = 8usize;
        let blowup = 4usize;
        let columns = vec![(0..trace_len).map(|i| F::from((42 + i) as u64)).collect::<Vec<F>>()];

        let a = lde_trace_columns_masked(&columns, trace_len, blowup, &[0], 6, 1).unwrap();
        let b = lde_trace_columns_masked(&columns, trace_len, blowup, &[0], 6, 2).unwrap();

        for k in 0..trace_len {
            assert_eq!(a[0][k * blowup], b[0][k * blowup]);
        }
        assert_ne!(a[0], b[0], "distinct seeds produced identical masked LDE");
    }

    #[test]
    fn empty_mask_is_identity() {
        let trace_len = 16usize;
        let blowup = 2usize;
        let columns: Vec<Vec<F>> = (0..3)
            .map(|c| (0..trace_len).map(|i| F::from((c * 17 + i) as u64)).collect())
            .collect();
        let plain = lde_trace_columns(&columns, trace_len, blowup).unwrap();
        let masked = lde_trace_columns_masked(&columns, trace_len, blowup, &[], 8, 999).unwrap();
        assert_eq!(plain, masked, "empty mask_cols must equal the default LDE exactly");
        let _ = F::ONE;
    }
}

#[cfg(test)]
mod starkware_tests {
    use super::*;

    #[test]
    fn import_starkware_json_basic() {
        let json = r#"{
            "format": "starkware-v1",
            "width": 3,
            "length": 4,
            "columns": {
                "pc": [0, 1, 2, 3],
                "ap": [100, 101, 102, 103],
                "fp": [100, 100, 100, 100]
            }
        }"#;
        let cols = import_starkware_json(json, None).unwrap();
        assert_eq!(cols.len(), 3);
        assert_eq!(cols[0].len(), 4);
    }

    #[test]
    fn import_starkware_json_ordered() {
        let json = r#"{
            "format": "starkware-v1",
            "width": 3,
            "length": 4,
            "columns": {
                "pc": [10, 11, 12, 13],
                "ap": [100, 101, 102, 103],
                "fp": [200, 200, 200, 200]
            }
        }"#;
        let order = ["pc", "ap", "fp"];
        let cols = import_starkware_json(json, Some(&order)).unwrap();
        assert_eq!(cols.len(), 3);
        assert_eq!(cols[0][0], F::from(10u64));
        assert_eq!(cols[1][0], F::from(100u64));
    }

    #[test]
    fn lde_trace_columns_produces_correct_size() {
        let cols: Vec<Vec<F>> = vec![
            (0u64..8).map(F::from).collect(),
            (0u64..8).map(F::from).collect(),
        ];
        let lde = lde_trace_columns(&cols, 8, 4).unwrap();
        assert_eq!(lde.len(), 2);
        assert_eq!(lde[0].len(), 32);
    }
}
