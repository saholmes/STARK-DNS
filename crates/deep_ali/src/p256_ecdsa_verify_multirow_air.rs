// p256_ecdsa_verify_multirow_air.rs — END-TO-END witness-binding
// multi-row P256 ECDSA *verification* AIR.
//
// Extends the narrow multi-row double-scalar-mult kernel
// (`p256_ecdsa_double_multirow_air`) with a small TAIL that finishes a
// full ECDSA-P256 verify in-circuit:
//
//   rows 0..K-1   : the two scalar-mult chains (one step/row)
//                   A = u1·G, B = u2·Q, with the projective outputs
//                   bound to the `r_a_proj` / `r_b_proj` columns at the
//                   last chain row (K-1).
//   row  K        : the TAIL —
//                     1. R = R_a + R_b              (group_add gadget)
//                     2. r_plus_n = r + n  (mod p)  (add + freeze)
//                     3. r·R.Z and (r+n)·R.Z        (two mul gadgets)
//                     4. selected = sel ? (r+n)·R.Z : r·R.Z   (select)
//                     5. R.X == selected            (10-limb equality)
//                     6. n column == n constant     (10-limb bind)
//   rows K+1..    : padding (all constraints zeroed).
//
// The cross-multiply (step 4/5) is the inverse-free affine-x check:
// the affine x-coordinate x1 = R.X · R.Z^{-1} ∈ [0,p) satisfies
// `x1 mod n == r`  iff  `x1 ∈ {r, r+n}`  iff
// `R.X ≡ r·R.Z (mod p)`  OR  `R.X ≡ (r+n)·R.Z (mod p)`.
// A boolean selector `sel` (inside the select gadget) picks the branch.
// This avoids the ~677k-cell Fermat inverse that would OOM.
//
// SOUNDNESS CHAIN: the scalar-mult chains bind R_a, R_b to the witness
// (u1/u2 bits, G, Q) via the kernel's transition constraints; the
// r_proj columns are constant (column-constancy) so the value bound to
// the chain output at row K-1 is exactly what the tail's group_add
// consumes at row K; r_plus_n is bound to r+n; and the final equality
// ties the projective R.X to the signature's r.  Tampering r (or any
// witness cell) makes some gadget constraint fire → FRI rejects.

#![allow(non_snake_case, non_upper_case_globals, dead_code)]

use ark_ff::{PrimeField, Zero};
use ark_goldilocks::Goldilocks as F;

use crate::p256_ecdsa_double_multirow_air::{
    build_ecdsa_double_multirow_layout, ecdsa_double_multirow_local_constraints,
    fill_ecdsa_double_multirow, EcdsaDoubleMultirowLayout,
};
use crate::p256_field::{FieldElement, NUM_LIMBS};
use crate::p256_field_air::{
    eval_add_gadget, eval_freeze_gadget, eval_mul_gadget, eval_select_gadget,
    fill_add_gadget, fill_freeze_gadget, fill_mul_gadget, fill_select_gadget,
    AddGadgetLayout, FreezeGadgetLayout, MulGadgetLayout, SelectGadgetLayout,
    ADD_GADGET_CONSTRAINTS, ELEMENT_BIT_CELLS, FREEZE_GADGET_CONSTRAINTS,
    MUL_GADGET_CONSTRAINTS, SELECT_GADGET_CONSTRAINTS,
};
use crate::p256_group_air::{
    alloc_add_layout, alloc_freeze_layout, alloc_mul_layout, build_group_add_layout,
    eval_group_add_gadget, fill_group_add_gadget, group_add_gadget_constraints,
    GroupAddGadgetLayout,
};
use crate::p256_scalar::N_LIMBS_TIGHT;

// ── per-row constraint-slot counts of the kernel sub-blocks ──
const DSM_BOUNDARY: usize = 6 * NUM_LIMBS; // r_proj boundary at row K-1
const DSM_ACC_TRANSITION: usize = 2 * 3 * NUM_LIMBS; // acc[r+1]=select[r]
const DSM_RPROJ_CONSTANCY: usize = 6 * NUM_LIMBS; // r_proj column-constancy

/// END-TO-END multi-row ECDSA-verify layout.
#[derive(Clone, Debug)]
pub struct EcdsaVerifyMultirowLayout {
    /// The reused double-scalar-mult kernel (chains A,B + r_proj cols).
    pub dsm: EcdsaDoubleMultirowLayout,
    /// R = R_a + R_b over the (constant) r_a_proj / r_b_proj columns.
    pub group_add: GroupAddGadgetLayout,
    /// Public column: signature `r` (as a mod-p field element, 10 limbs).
    pub r_base: usize,
    /// Public column: curve order `n` (10 limbs, bound to the constant).
    pub n_const_base: usize,
    /// r + n (integer, < 2p).
    pub add_rn: AddGadgetLayout,
    /// canonical (r + n) mod p.
    pub freeze_rn: FreezeGadgetLayout,
    /// r · R.Z (mod p).
    pub mul_r: MulGadgetLayout,
    /// (r+n) · R.Z (mod p).
    pub mul_rn: MulGadgetLayout,
    /// selected = sel ? mul_rn : mul_r.
    pub select: SelectGadgetLayout,
    /// Number of scalar-mult steps (= number of MSB bits, e.g. 256).
    pub k_steps: usize,
    pub width: usize,
}

pub fn build_ecdsa_verify_multirow_layout(
    start: usize,
    k_steps: usize,
) -> (EcdsaVerifyMultirowLayout, usize) {
    let (dsm, dsm_end) = build_ecdsa_double_multirow_layout(start);
    let mut cursor = dsm_end;

    // R = R_a + R_b over the constant r_proj columns.
    let (group_add, ga_end) = build_group_add_layout(
        cursor,
        dsm.r_a_proj_x_base,
        dsm.r_a_proj_y_base,
        dsm.r_a_proj_z_base,
        dsm.r_b_proj_x_base,
        dsm.r_b_proj_y_base,
        dsm.r_b_proj_z_base,
    );
    cursor = ga_end;

    // Public columns: r and n (limbs only — gadget inputs need no bits).
    let r_base = cursor;
    cursor += NUM_LIMBS;
    let n_const_base = cursor;
    cursor += NUM_LIMBS;

    // r_plus_n = freeze(r + n).
    let add_rn = alloc_add_layout(&mut cursor, r_base, n_const_base);
    let freeze_rn = alloc_freeze_layout(&mut cursor, add_rn.c_limbs_base);
    let r_plus_n_base = freeze_rn.c_limbs_base;

    // Cross-multiply muls against R.Z.
    let z3 = group_add.result_z3_limbs_base;
    let mul_r = alloc_mul_layout(&mut cursor, r_base, z3);
    let mul_rn = alloc_mul_layout(&mut cursor, r_plus_n_base, z3);

    // select: sel ? mul_rn : mul_r.
    let sel_cell = cursor;
    cursor += 1;
    let sel_c_limbs = cursor;
    cursor += NUM_LIMBS;
    let sel_c_bits = cursor;
    cursor += ELEMENT_BIT_CELLS;
    let select = SelectGadgetLayout {
        a_limbs_base: mul_rn.c_limbs_base,
        b_limbs_base: mul_r.c_limbs_base,
        c_limbs_base: sel_c_limbs,
        c_bits_base: sel_c_bits,
        sel_cell,
    };

    let layout = EcdsaVerifyMultirowLayout {
        dsm,
        group_add,
        r_base,
        n_const_base,
        add_rn,
        freeze_rn,
        mul_r,
        mul_rn,
        select,
        k_steps,
        width: cursor,
    };
    (layout, cursor)
}

/// Number of tail (row-K) constraint slots.
pub fn ecdsa_verify_tail_constraints(layout: &EcdsaVerifyMultirowLayout) -> usize {
    group_add_gadget_constraints(&layout.group_add)
        + ADD_GADGET_CONSTRAINTS
        + FREEZE_GADGET_CONSTRAINTS
        + 2 * MUL_GADGET_CONSTRAINTS
        + SELECT_GADGET_CONSTRAINTS
        + NUM_LIMBS // R.X == selected
        + NUM_LIMBS // n column == n constant
}

pub fn ecdsa_verify_multirow_constraints(layout: &EcdsaVerifyMultirowLayout) -> usize {
    ecdsa_double_multirow_local_constraints(&layout.dsm)
        + DSM_BOUNDARY
        + DSM_ACC_TRANSITION
        + DSM_RPROJ_CONSTANCY
        + ecdsa_verify_tail_constraints(layout)
}

#[inline]
fn read_fe_row(trace: &[Vec<F>], base: usize, row: usize) -> FieldElement {
    let mut limbs = [0i64; NUM_LIMBS];
    for i in 0..NUM_LIMBS {
        let bi = trace[base + i][row].into_bigint();
        limbs[i] = bi.as_ref()[0] as i64;
    }
    FieldElement { limbs }
}

#[inline]
fn place_limbs(trace: &mut [Vec<F>], base: usize, row: usize, fe: &FieldElement) {
    for i in 0..NUM_LIMBS {
        trace[base + i][row] = F::from(fe.limbs[i] as u64);
    }
}

/// `n` (curve order) as a tight-form mod-p FieldElement.
pub fn order_n_field() -> FieldElement {
    FieldElement { limbs: *N_LIMBS_TIGHT }
}

/// Fill the full verify trace.
///
/// `n_trace >= k_steps + 1` (power of two).  The kernel runs `k_steps`
/// scalar-mult steps on rows `0..k_steps`; the tail is filled at row
/// `k_steps`.  Two-pass: the caller first fills with placeholder r_proj
/// to capture the chain outputs at row `k_steps-1`, then refills with
/// the captured projective outputs.  This function performs BOTH passes
/// internally and additionally fills the tail.
///
/// `r_scalar_fe` is the signature `r` reduced into a mod-p FieldElement.
#[allow(clippy::too_many_arguments)]
pub fn fill_ecdsa_verify_multirow(
    trace: &mut [Vec<F>],
    layout: &EcdsaVerifyMultirowLayout,
    n_trace: usize,
    a_initial: (&FieldElement, &FieldElement, &FieldElement),
    a_base: (&FieldElement, &FieldElement, &FieldElement),
    a_bits: &[bool],
    b_initial: (&FieldElement, &FieldElement, &FieldElement),
    b_base: (&FieldElement, &FieldElement, &FieldElement),
    b_bits: &[bool],
    r_scalar_fe: &FieldElement,
) {
    let k = layout.k_steps;
    assert!(n_trace.is_power_of_two());
    assert!(n_trace >= k + 1, "need a tail row at index k_steps");
    assert_eq!(a_bits.len(), k);
    assert_eq!(b_bits.len(), k);

    let zfe = FieldElement::zero();

    // ── Pass 1: r_proj = 0 to capture chain outputs at row k-1. ──
    fill_ecdsa_double_multirow(
        trace, &layout.dsm, n_trace, k, k,
        a_initial.0, a_initial.1, a_initial.2,
        a_base.0, a_base.1, a_base.2, a_bits,
        b_initial.0, b_initial.1, b_initial.2,
        b_base.0, b_base.1, b_base.2, b_bits,
        &zfe, &zfe, &zfe, &zfe, &zfe, &zfe,
    );

    let last = k - 1;
    let r_a_x = read_fe_row(trace, layout.dsm.step_a.select_x.c_limbs_base, last);
    let r_a_y = read_fe_row(trace, layout.dsm.step_a.select_y.c_limbs_base, last);
    let r_a_z = read_fe_row(trace, layout.dsm.step_a.select_z.c_limbs_base, last);
    let r_b_x = read_fe_row(trace, layout.dsm.step_b.select_x.c_limbs_base, last);
    let r_b_y = read_fe_row(trace, layout.dsm.step_b.select_y.c_limbs_base, last);
    let r_b_z = read_fe_row(trace, layout.dsm.step_b.select_z.c_limbs_base, last);

    // ── Pass 2: refill kernel with the captured r_proj values. ──
    fill_ecdsa_double_multirow(
        trace, &layout.dsm, n_trace, k, k,
        a_initial.0, a_initial.1, a_initial.2,
        a_base.0, a_base.1, a_base.2, a_bits,
        b_initial.0, b_initial.1, b_initial.2,
        b_base.0, b_base.1, b_base.2, b_bits,
        &r_a_x, &r_a_y, &r_a_z, &r_b_x, &r_b_y, &r_b_z,
    );

    // ── Tail at row k. ──
    let row = k;
    let n_fe = order_n_field();

    // R = R_a + R_b.
    fill_group_add_gadget(
        trace, row, &layout.group_add,
        &r_a_x, &r_a_y, &r_a_z, &r_b_x, &r_b_y, &r_b_z,
    );
    let r_x3 = read_fe_row(trace, layout.group_add.result_x3_limbs_base, row);
    let r_z3 = read_fe_row(trace, layout.group_add.result_z3_limbs_base, row);

    // Affine x1 = R.X · R.Z^{-1}; choose sel.
    let mut x1 = r_x3.mul(&r_z3.invert());
    x1.freeze();
    let mut r_can = *r_scalar_fe;
    r_can.freeze();
    let mut rn_native = r_scalar_fe.add(&n_fe);
    rn_native.freeze();
    let sel = if x1.ct_eq(&rn_native) {
        true
    } else {
        debug_assert!(x1.ct_eq(&r_can), "honest fill: x1 not in {{r, r+n}}");
        false
    };

    // r, n public columns.
    place_limbs(trace, layout.r_base, row, r_scalar_fe);
    place_limbs(trace, layout.n_const_base, row, &n_fe);

    // r_plus_n = freeze(r + n).
    fill_add_gadget(trace, row, &layout.add_rn, r_scalar_fe, &n_fe);
    let rn_raw = read_fe_row(trace, layout.add_rn.c_limbs_base, row);
    fill_freeze_gadget(trace, row, &layout.freeze_rn, &rn_raw);
    let r_plus_n = read_fe_row(trace, layout.freeze_rn.c_limbs_base, row);

    // Cross-multiply.
    fill_mul_gadget(trace, row, &layout.mul_r, r_scalar_fe, &r_z3);
    fill_mul_gadget(trace, row, &layout.mul_rn, &r_plus_n, &r_z3);
    let mul_r_c = read_fe_row(trace, layout.mul_r.c_limbs_base, row);
    let mul_rn_c = read_fe_row(trace, layout.mul_rn.c_limbs_base, row);

    // select: sel ? mul_rn : mul_r.
    fill_select_gadget(trace, row, &layout.select, &mul_rn_c, &mul_r_c, sel);
}

/// Per-row evaluator (transition-aware).  Order MUST match
/// `ecdsa_verify_multirow_constraints`.
pub fn eval_ecdsa_verify_multirow_per_row(
    cur: &[F],
    nxt: &[F],
    trace_row: usize,
    _n_trace: usize,
    layout: &EcdsaVerifyMultirowLayout,
) -> Vec<F> {
    use crate::p256_scalar_mul_air::eval_scalar_mul_step_gadget;

    let k = layout.k_steps;
    let total = ecdsa_verify_multirow_constraints(layout);
    let mut out = Vec::with_capacity(total);

    let in_chain = trace_row < k; // rows 0..k-1 run the scalar mult
    let is_chain_last = trace_row + 1 == k; // row k-1: bind r_proj
    let is_tail = trace_row == k; // row k: the verify tail

    // (1) DSM local: both step gadgets, gated to chain rows.
    let dsm_local = ecdsa_double_multirow_local_constraints(&layout.dsm);
    if in_chain {
        out.extend(eval_scalar_mul_step_gadget(cur, &layout.dsm.step_a));
        out.extend(eval_scalar_mul_step_gadget(cur, &layout.dsm.step_b));
        debug_assert_eq!(out.len(), dsm_local);
    } else {
        out.resize(dsm_local, F::zero());
    }

    // (2) DSM boundary: bind chain outputs to r_proj cols at row k-1.
    let push_boundary = |out: &mut Vec<F>, chain_base: usize, proj_base: usize| {
        for i in 0..NUM_LIMBS {
            if is_chain_last {
                out.push(cur[chain_base + i] - cur[proj_base + i]);
            } else {
                out.push(F::zero());
            }
        }
    };
    push_boundary(&mut out, layout.dsm.step_a.select_x.c_limbs_base, layout.dsm.r_a_proj_x_base);
    push_boundary(&mut out, layout.dsm.step_a.select_y.c_limbs_base, layout.dsm.r_a_proj_y_base);
    push_boundary(&mut out, layout.dsm.step_a.select_z.c_limbs_base, layout.dsm.r_a_proj_z_base);
    push_boundary(&mut out, layout.dsm.step_b.select_x.c_limbs_base, layout.dsm.r_b_proj_x_base);
    push_boundary(&mut out, layout.dsm.step_b.select_y.c_limbs_base, layout.dsm.r_b_proj_y_base);
    push_boundary(&mut out, layout.dsm.step_b.select_z.c_limbs_base, layout.dsm.r_b_proj_z_base);

    // (3) DSM acc transition: acc[r+1] = select[r], for r in 0..k-2.
    let acc_link = trace_row + 1 < k;
    if acc_link {
        for (acc, sel) in [
            (layout.dsm.step_a.acc_x_base, layout.dsm.step_a.select_x.c_limbs_base),
            (layout.dsm.step_a.acc_y_base, layout.dsm.step_a.select_y.c_limbs_base),
            (layout.dsm.step_a.acc_z_base, layout.dsm.step_a.select_z.c_limbs_base),
            (layout.dsm.step_b.acc_x_base, layout.dsm.step_b.select_x.c_limbs_base),
            (layout.dsm.step_b.acc_y_base, layout.dsm.step_b.select_y.c_limbs_base),
            (layout.dsm.step_b.acc_z_base, layout.dsm.step_b.select_z.c_limbs_base),
        ] {
            for i in 0..NUM_LIMBS {
                out.push(nxt[acc + i] - cur[sel + i]);
            }
        }
    } else {
        let base = out.len();
        out.resize(base + DSM_ACC_TRANSITION, F::zero());
    }

    // (4) DSM r_proj column-constancy, for r in 0..k-1 (links k-1 -> k).
    let rproj_link = trace_row < k;
    if rproj_link {
        for base in &[
            layout.dsm.r_a_proj_x_base, layout.dsm.r_a_proj_y_base, layout.dsm.r_a_proj_z_base,
            layout.dsm.r_b_proj_x_base, layout.dsm.r_b_proj_y_base, layout.dsm.r_b_proj_z_base,
        ] {
            for i in 0..NUM_LIMBS {
                out.push(nxt[*base + i] - cur[*base + i]);
            }
        }
    } else {
        let base = out.len();
        out.resize(base + DSM_RPROJ_CONSTANCY, F::zero());
    }

    // (5) TAIL at row k.
    let tail_count = ecdsa_verify_tail_constraints(layout);
    if is_tail {
        let tail_start = out.len();
        out.extend(eval_group_add_gadget(cur, &layout.group_add));
        out.extend(eval_add_gadget(cur, &layout.add_rn));
        out.extend(eval_freeze_gadget(cur, &layout.freeze_rn));
        out.extend(eval_mul_gadget(cur, &layout.mul_r));
        out.extend(eval_mul_gadget(cur, &layout.mul_rn));
        out.extend(eval_select_gadget(cur, &layout.select));
        // R.X == selected.
        for i in 0..NUM_LIMBS {
            out.push(cur[layout.group_add.result_x3_limbs_base + i] - cur[layout.select.c_limbs_base + i]);
        }
        // n column == n constant.
        let n_fe = order_n_field();
        for i in 0..NUM_LIMBS {
            out.push(cur[layout.n_const_base + i] - F::from(n_fe.limbs[i] as u64));
        }
        debug_assert_eq!(out.len() - tail_start, tail_count);
    } else {
        let base = out.len();
        out.resize(base + tail_count, F::zero());
    }

    debug_assert_eq!(out.len(), total);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::p256_group::GENERATOR;

    fn z_one() -> FieldElement {
        let mut t = FieldElement::zero();
        t.limbs[0] = 1;
        t
    }
    fn identity() -> (FieldElement, FieldElement, FieldElement) {
        // Projective point at infinity (0:1:0).
        let mut y = FieldElement::zero();
        y.limbs[0] = 1;
        (FieldElement::zero(), y, FieldElement::zero())
    }

    // Small smoke test: K=4 scalar mult + tail, native cross-check via
    // a constructed (R, r) that we KNOW lands on the equality.
    #[test]
    fn verify_multirow_k4_tail_consistent() {
        let k = 4usize;
        let (layout, total) = build_ecdsa_verify_multirow_layout(0, k);
        let n_trace = (k + 1).next_power_of_two(); // 8

        let g = *GENERATOR;
        let q = g.double();
        let zo = z_one();
        let (ix, iy, iz) = identity();

        // Drive both chains with bits; the resulting R is some real
        // point.  We must set r = R.x mod n for the tail equality.
        let a_bits = vec![true, false, true, true];
        let b_bits = vec![false, true, true, false];

        // Compute the chains natively-ish by filling pass 1 to read R_a,R_b,
        // then group_add, then derive r = affine_x(R) mod ... we just read
        // it back after a trial fill below.
        let mut trace = vec![vec![F::zero(); n_trace]; total];

        // First fill with a placeholder r to obtain R, then recompute r.
        // We do a throwaway pass to read R, choosing r afterwards.
        let placeholder_r = FieldElement::zero();
        fill_ecdsa_verify_multirow(
            &mut trace, &layout, n_trace,
            (&ix, &iy, &iz), (&g.x, &g.y, &zo), &a_bits,
            (&ix, &iy, &iz), (&q.x, &q.y, &zo), &b_bits,
            &placeholder_r,
        );
        // Read affine x1 of R from the (correct) group_add output at row k.
        let r_x3 = read_fe_row(&trace, layout.group_add.result_x3_limbs_base, k);
        let r_z3 = read_fe_row(&trace, layout.group_add.result_z3_limbs_base, k);
        let mut x1 = r_x3.mul(&r_z3.invert());
        x1.freeze();
        // Use r := x1 (as a mod-p element; sel=0 branch).  This is a
        // self-consistent in-circuit check of the cross-multiply.
        let r_fe = x1;

        // Refill with the real r.
        let mut trace = vec![vec![F::zero(); n_trace]; total];
        fill_ecdsa_verify_multirow(
            &mut trace, &layout, n_trace,
            (&ix, &iy, &iz), (&g.x, &g.y, &zo), &a_bits,
            (&ix, &iy, &iz), (&q.x, &q.y, &zo), &b_bits,
            &r_fe,
        );

        let mut failures = 0usize;
        for r in 0..n_trace {
            let cur: Vec<F> = (0..total).map(|c| trace[c][r]).collect();
            let nxt: Vec<F> = (0..total).map(|c| trace[c][(r + 1) % n_trace]).collect();
            let cons = eval_ecdsa_verify_multirow_per_row(&cur, &nxt, r, n_trace, &layout);
            failures += cons.iter().filter(|v| !v.is_zero()).count();
        }
        assert_eq!(failures, 0, "verify-multirow K=4 had {failures} non-zero constraints");
    }

    #[test]
    fn verify_multirow_tampered_r_violates() {
        let k = 4usize;
        let (layout, total) = build_ecdsa_verify_multirow_layout(0, k);
        let n_trace = (k + 1).next_power_of_two();
        let g = *GENERATOR;
        let q = g.double();
        let zo = z_one();
        let (ix, iy, iz) = identity();
        let a_bits = vec![true, false, true, true];
        let b_bits = vec![false, true, true, false];

        let mut trace = vec![vec![F::zero(); n_trace]; total];
        fill_ecdsa_verify_multirow(
            &mut trace, &layout, n_trace,
            (&ix, &iy, &iz), (&g.x, &g.y, &zo), &a_bits,
            (&ix, &iy, &iz), (&q.x, &q.y, &zo), &b_bits,
            &FieldElement::zero(),
        );
        let r_x3 = read_fe_row(&trace, layout.group_add.result_x3_limbs_base, k);
        let r_z3 = read_fe_row(&trace, layout.group_add.result_z3_limbs_base, k);
        let mut x1 = r_x3.mul(&r_z3.invert());
        x1.freeze();

        let mut trace = vec![vec![F::zero(); n_trace]; total];
        fill_ecdsa_verify_multirow(
            &mut trace, &layout, n_trace,
            (&ix, &iy, &iz), (&g.x, &g.y, &zo), &a_bits,
            (&ix, &iy, &iz), (&q.x, &q.y, &zo), &b_bits,
            &x1,
        );
        // Tamper the r public column -> mul_r position identity fires.
        trace[layout.r_base][k] += F::from(1u64);

        let mut failures = 0usize;
        for r in 0..n_trace {
            let cur: Vec<F> = (0..total).map(|c| trace[c][r]).collect();
            let nxt: Vec<F> = (0..total).map(|c| trace[c][(r + 1) % n_trace]).collect();
            let cons = eval_ecdsa_verify_multirow_per_row(&cur, &nxt, r, n_trace, &layout);
            failures += cons.iter().filter(|v| !v.is_zero()).count();
        }
        assert!(failures >= 1, "tampered r must violate >=1 constraint");
    }
}
