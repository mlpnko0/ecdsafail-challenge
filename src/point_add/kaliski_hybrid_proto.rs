//! Main-circuit profiling prototypes for the hybrid Kaliski-jump moonshot.
//!
//! This file started with a lower-bound question:
//!
//! > If the 3-step bulk prefix were already known exactly up front, how cheap
//! > would a branch-free fused implementation be?
//!
//! But the more important follow-up is to stress-test the hidden assumptions:
//!
//! - **Assumption A (derivability)**: can `cmp1, cmp2` be known before running
//!   the first two micro-steps?
//! - **Assumption B (cleanup)**: once compare bits are used for selection,
//!   should they be kept live across the core or recomputed from the output
//!   state for selector cleanup?
//!
//! The branch-free `exact3` profile is therefore just a lower bound. The more
//! realistic staged profiles below quantify what happens once those assumptions
//! are checked against the real builder.

use std::collections::BTreeMap;

use super::{
    add_nbit_qq_fast, cmp_lt_into_fast, cswap, kaliski_iteration, mod_double_no_corr,
    sub_nbit_qq_fast, with_gt, B, N, OperationType, QubitId, SECP256K1_P,
};
use super::kaliski_jump::{KCase, Sampler, observe_window};
use super::test_timeout::{check_deadline, two_min_deadline};

#[derive(Debug, Clone, Copy, Default)]
pub struct ProtoCost {
    pub ccx: u64,
    pub cliff: u64,
    pub other: u64,
    pub ops: usize,
    pub peak_qubits: u32,
}

#[derive(Debug, Clone)]
pub struct PrefixCostRow {
    pub prefix: [KCase; 3],
    pub windows: usize,
    pub exact3: ProtoCost,
    pub staged_exact3: ProtoCost,
    pub recompute_cleanup: ProtoCost,
}

#[derive(Debug, Clone)]
pub struct HybridProtoSummary {
    pub full_windows: usize,
    pub distinct_prefixes: usize,

    pub baseline3: ProtoCost,
    pub baseline4: ProtoCost,
    pub baseline_step3: ProtoCost,

    pub weighted_exact3_ccx: f64,
    pub weighted_exact3_cliff: f64,
    pub min_exact3_ccx: u64,
    pub max_exact3_ccx: u64,

    pub weighted_hybrid4_lower_bound_ccx: f64,
    pub weighted_hybrid4_lower_bound_cliff: f64,
    pub lower_bound_ccx_savings_vs_baseline4: f64,
    pub lower_bound_ccx_savings_pct_vs_baseline4: f64,

    pub gt_compare: ProtoCost,

    pub weighted_staged_exact3_ccx: f64,
    pub weighted_staged_exact3_cliff: f64,
    pub weighted_recompute_cleanup_ccx: f64,
    pub weighted_recompute_cleanup_cliff: f64,
    pub total_keep_live_ccx: f64,
    pub total_keep_live_cliff: f64,
    pub total_recompute_ccx: f64,
    pub total_recompute_cliff: f64,

    pub real_bulk3_forward: ProtoCost,
    pub real_bulk3_specialized: ProtoCost,

    pub top_prefix_rows: Vec<PrefixCostRow>,
}

fn summarize_builder(b: &B) -> ProtoCost {
    let mut out = ProtoCost::default();
    out.ops = b.ops.len();
    out.peak_qubits = b.peak_qubits;
    for op in &b.ops {
        match op.kind {
            OperationType::CCX | OperationType::CCZ => out.ccx += 1,
            OperationType::CX
            | OperationType::CZ
            | OperationType::Swap
            | OperationType::Hmr
            | OperationType::R => out.cliff += 1,
            _ => out.other += 1,
        }
    }
    out
}

fn shift_right_known_even(b: &mut B, v: &[QubitId]) {
    for i in 0..(v.len() - 1) {
        b.swap(v[i], v[i + 1]);
    }
}

fn prefix_to_string(prefix: [KCase; 3]) -> String {
    let mut s = String::new();
    for (i, kc) in prefix.into_iter().enumerate() {
        if i > 0 { s.push('-'); }
        s.push_str(match kc {
            KCase::UEven => "UE",
            KCase::VEven => "VE",
            KCase::UGtV => "UG",
            KCase::VGtU => "VG",
        });
    }
    s
}

/// Branch-free realization of a *known* Kaliski case, intended only as a
/// lower-bound prototype for profiling. This is not yet a full reversible
/// hybrid primitive: there is no selector, no cleanup, and no tail handling.
fn c_shift_left_no_corr(b: &mut B, v: &[QubitId], ctrl: QubitId) {
    for i in (0..(v.len() - 1)).rev() {
        cswap(b, ctrl, v[i], v[i + 1]);
    }
}

fn cucc_add_ctrl_fast(b: &mut B, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = a.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    add_nbit_qq_fast(b, &tmp, acc);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    b.free_vec(&tmp);
}

fn cucc_sub_ctrl_fast(b: &mut B, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = a.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    sub_nbit_qq_fast(b, &tmp, acc);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    b.free_vec(&tmp);
}

fn exact_known_case_step(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    iter_idx: usize,
    kc: KCase,
) {
    let rs_add_width = (iter_idx + 2).min(N);
    match kc {
        KCase::UEven => {
            // (u, v) -> (u/2, v), (r, s) -> (r, 2s)
            shift_right_known_even(b, u);
            mod_double_no_corr(b, s);
        }
        KCase::VEven => {
            // (u, v) -> (u, v/2), (r, s) -> (2r, s)
            shift_right_known_even(b, v);
            mod_double_no_corr(b, r);
        }
        KCase::UGtV => {
            // (u, v) -> ((u-v)/2, v), (r, s) -> (r+s, 2s)
            sub_nbit_qq_fast(b, v, u);
            shift_right_known_even(b, u);
            let s_slice = s[..rs_add_width].to_vec();
            let r_slice = r[..rs_add_width].to_vec();
            add_nbit_qq_fast(b, &s_slice, &r_slice);
            mod_double_no_corr(b, s);
        }
        KCase::VGtU => {
            // (u, v) -> (u, (v-u)/2), (r, s) -> (2r, r+s)
            sub_nbit_qq_fast(b, u, v);
            shift_right_known_even(b, v);
            let r_slice = r[..rs_add_width].to_vec();
            let s_slice = s[..rs_add_width].to_vec();
            add_nbit_qq_fast(b, &r_slice, &s_slice);
            mod_double_no_corr(b, r);
        }
    }
}

fn exact_known_case_step_inverse(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    iter_idx: usize,
    kc: KCase,
) {
    let rs_add_width = (iter_idx + 2).min(N);
    match kc {
        KCase::UEven => {
            mod_double_no_corr(b, u);
            shift_right_known_even(b, s);
        }
        KCase::VEven => {
            mod_double_no_corr(b, v);
            shift_right_known_even(b, r);
        }
        KCase::UGtV => {
            shift_right_known_even(b, s);
            let s_slice = s[..rs_add_width].to_vec();
            let r_slice = r[..rs_add_width].to_vec();
            sub_nbit_qq_fast(b, &s_slice, &r_slice);
            mod_double_no_corr(b, u);
            add_nbit_qq_fast(b, v, u);
        }
        KCase::VGtU => {
            shift_right_known_even(b, r);
            let r_slice = r[..rs_add_width].to_vec();
            let s_slice = s[..rs_add_width].to_vec();
            sub_nbit_qq_fast(b, &r_slice, &s_slice);
            mod_double_no_corr(b, v);
            add_nbit_qq_fast(b, u, v);
        }
    }
}

fn emit_gt_compare(b: &mut B, u: &[QubitId], v: &[QubitId]) {
    let flag = b.alloc_qubit();
    with_gt(b, u, v, flag, |_b| {});
    b.free(flag);
}

fn profile_exact_prefix(prefix: [KCase; 3]) -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    for (i, kc) in prefix.into_iter().enumerate() {
        b.set_phase("hybrid_exact3_step");
        exact_known_case_step(&mut b, &u, &v, &r, &s, i, kc);
    }
    summarize_builder(&b)
}

fn profile_baseline_iters(n_iters: usize) -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m_hist = b.alloc_qubits(n_iters);
    let f = b.alloc_qubit();
    for i in 0..n_iters {
        b.set_phase("hybrid_baseline_iter");
        kaliski_iteration(&mut b, SECP256K1_P, &u, &v, &r, &s, m_hist[i], f, i);
    }
    summarize_builder(&b)
}

fn profile_single_baseline_iter(iter_idx: usize) -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    let f = b.alloc_qubit();
    b.set_phase("hybrid_baseline_iter");
    kaliski_iteration(&mut b, SECP256K1_P, &u, &v, &r, &s, m, f, iter_idx);
    summarize_builder(&b)
}

fn profile_staged_exact3(prefix: [KCase; 3]) -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    for (i, kc) in prefix.into_iter().enumerate() {
        b.set_phase("hybrid_staged_cmp");
        emit_gt_compare(&mut b, &u, &v);
        b.set_phase("hybrid_staged_step");
        exact_known_case_step(&mut b, &u, &v, &r, &s, i, kc);
    }
    summarize_builder(&b)
}

fn profile_recompute_cleanup(prefix: [KCase; 3]) -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);

    // Recompute compare bits from the output state of the exact 3-step core by
    // running the known prefix backward, comparing at the recovered states, and
    // then restoring the output state.
    for (i, kc) in prefix.into_iter().enumerate().rev() {
        b.set_phase("hybrid_cleanup_inv_step");
        exact_known_case_step_inverse(&mut b, &u, &v, &r, &s, i, kc);
        b.set_phase("hybrid_cleanup_cmp");
        emit_gt_compare(&mut b, &u, &v);
    }
    for (i, kc) in prefix.into_iter().enumerate() {
        b.set_phase("hybrid_cleanup_fwd_restore");
        exact_known_case_step(&mut b, &u, &v, &r, &s, i, kc);
    }
    summarize_builder(&b)
}

fn profile_gt_compare() -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    emit_gt_compare(&mut b, &u, &v);
    summarize_builder(&b)
}

#[derive(Debug)]
struct BulkStepLive {
    u_even: QubitId,
    v_even: QubitId,
    both_odd: QubitId,
    cmp: QubitId,
}

fn bulk_step_forward_keep_live(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    iter_idx: usize,
) -> BulkStepLive {
    let rs_add_width = (iter_idx + 2).min(N);

    let u_even = b.alloc_qubit();
    b.x(u_even);
    b.cx(u[0], u_even); // u_even = !u0

    let v_even = b.alloc_qubit();
    b.x(v[0]);
    b.ccx(u[0], v[0], v_even); // v_even = u0 & !v0
    b.x(v[0]);

    let both_odd = b.alloc_qubit();
    b.ccx(u[0], v[0], both_odd);

    // Even cases.
    for i in 0..(u.len() - 1) { cswap(b, u_even, u[i], u[i + 1]); }
    c_shift_left_no_corr(b, s, u_even);

    for i in 0..(v.len() - 1) { cswap(b, v_even, v[i], v[i + 1]); }
    c_shift_left_no_corr(b, r, v_even);

    // Always compute compare for now; derivability tests showed cmp1/cmp2 are
    // not available up front. Keep cmp live for cleanup.
    let cmp = b.alloc_qubit();
    cmp_lt_into_fast(b, v, u, cmp); // cmp = (u > v)

    let ug = b.alloc_qubit();
    b.ccx(both_odd, cmp, ug);
    let vg = b.alloc_qubit();
    b.x(cmp);
    b.ccx(both_odd, cmp, vg);
    b.x(cmp);

    cucc_sub_ctrl_fast(b, v, u, ug);
    for i in 0..(u.len() - 1) { cswap(b, ug, u[i], u[i + 1]); }
    let s_slice = s[..rs_add_width].to_vec();
    let r_slice = r[..rs_add_width].to_vec();
    cucc_add_ctrl_fast(b, &s_slice, &r_slice, ug);
    c_shift_left_no_corr(b, s, ug);

    cucc_sub_ctrl_fast(b, u, v, vg);
    for i in 0..(v.len() - 1) { cswap(b, vg, v[i], v[i + 1]); }
    let r2_slice = r[..rs_add_width].to_vec();
    let s2_slice = s[..rs_add_width].to_vec();
    cucc_add_ctrl_fast(b, &r2_slice, &s2_slice, vg);
    c_shift_left_no_corr(b, r, vg);

    // ug/vg are derivable from cmp + original parities, so drop their direct
    // conjunctions and keep the base flags live.
    b.ccx(both_odd, cmp, ug);
    b.x(cmp);
    b.ccx(both_odd, cmp, vg);
    b.x(cmp);
    b.free(ug);
    b.free(vg);

    BulkStepLive { u_even, v_even, both_odd, cmp }
}

fn kaliski_bulk_iteration_specialized(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    iter_idx: usize,
) {
    let a_f = b.alloc_qubit();
    let b_f = b.alloc_qubit();
    let add_f = b.alloc_qubit();
    let m_i = b.alloc_qubit();

    // Specialized nonterminal, f=1 version of STEP 1.
    b.x(a_f);
    b.cx(u[0], a_f); // a_f = !u0
    b.x(v[0]);
    b.ccx(u[0], v[0], m_i); // m_i = u0 & !v0
    b.x(v[0]);
    b.cx(a_f, b_f);
    b.cx(m_i, b_f); // b_f = a_f xor m_i

    // Specialized STEP 2: use compare only for odd/odd cases.
    let l_gt = b.alloc_qubit();
    with_gt(b, u, v, l_gt, |b| {
        b.x(b_f);
        let t = b.alloc_qubit();
        b.ccx(l_gt, b_f, t);
        b.cx(t, a_f);
        b.cx(t, m_i);
        b.ccx(l_gt, b_f, t);
        b.free(t);
        b.x(b_f);
    });
    b.free(l_gt);

    // STEP 3 cswaps.
    for j in 0..u.len() { cswap(b, a_f, u[j], v[j]); }
    let rs_width_step3 = (iter_idx + 1).min(N);
    for j in 0..rs_width_step3 { cswap(b, a_f, r[j], s[j]); }

    // STEP 4 with add_f = !b_f.
    b.x(add_f);
    b.cx(b_f, add_f);
    {
        let tmp = b.alloc_qubits(N);
        for i in 0..N { b.ccx(add_f, u[i], tmp[i]); }
        sub_nbit_qq_fast(b, &tmp, v);
        let transform_width = (iter_idx + 1).min(N);
        for i in 0..transform_width { b.cx(r[i], u[i]); }
        for i in 0..transform_width { b.ccx(add_f, u[i], tmp[i]); }
        for i in 0..transform_width { b.cx(r[i], u[i]); }
        let add_width = (iter_idx + 2).min(N);
        let mut tmp_slice: Vec<QubitId> = tmp[..transform_width].to_vec();
        let tmp_pad = if add_width > transform_width {
            let q = b.alloc_qubit();
            tmp_slice.push(q);
            Some(q)
        } else {
            None
        };
        let s_slice: Vec<QubitId> = s[..add_width].to_vec();
        add_nbit_qq_fast(b, &tmp_slice, &s_slice);
        if let Some(q) = tmp_pad { b.free(q); }
        for i in 0..N {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            if i < transform_width {
                b.cz_if(add_f, r[i], m);
            } else {
                b.cz_if(add_f, u[i], m);
            }
        }
        b.free_vec(&tmp);
    }

    // STEP 5 uncompute add_f,b_f.
    b.cx(b_f, add_f);
    b.x(add_f);
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);

    // STEP 6-8 unconditional shift/double.
    for i in 0..(N - 1) { b.swap(v[i], v[i + 1]); }
    mod_double_no_corr(b, r);

    // STEP 9 cswaps again.
    for j in 0..u.len() { cswap(b, a_f, u[j], v[j]); }
    let rs_width_step9 = (iter_idx + 2).min(N);
    for j in 0..rs_width_step9 { cswap(b, a_f, r[j], s[j]); }

    // STEP 10 uncompute a.
    b.x(s[0]);
    b.cx(s[0], a_f);
    b.x(s[0]);

    b.free(m_i);
    b.free(add_f);
    b.free(b_f);
    b.free(a_f);
}

fn profile_real_bulk3_forward() -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let _live0 = bulk_step_forward_keep_live(&mut b, &u, &v, &r, &s, 0);
    let _live1 = bulk_step_forward_keep_live(&mut b, &u, &v, &r, &s, 1);
    let _live2 = bulk_step_forward_keep_live(&mut b, &u, &v, &r, &s, 2);
    summarize_builder(&b)
}

fn profile_real_bulk3_specialized() -> ProtoCost {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    kaliski_bulk_iteration_specialized(&mut b, &u, &v, &r, &s, 0);
    kaliski_bulk_iteration_specialized(&mut b, &u, &v, &r, &s, 1);
    kaliski_bulk_iteration_specialized(&mut b, &u, &v, &r, &s, 2);
    summarize_builder(&b)
}

fn collect_full_prefix_counts(seed: &[u8], n_inputs: usize, w: usize) -> BTreeMap<[KCase; 3], usize> {
    let deadline = two_min_deadline();
    let mut sampler = Sampler::new(seed, SECP256K1_P);
    let mut out = BTreeMap::new();
    for input_idx in 0..n_inputs {
        if (input_idx & 31) == 0 { check_deadline(deadline, "kaliski_hybrid_proto::collect_full_prefix_counts"); }
        let mut u = SECP256K1_P;
        let mut v = sampler.next();
        for _ in 0..742 {
            if v.is_zero() { break; }
            let (_u4, _v4, obs) = observe_window(u, v, w, 4);
            if obs.cases.len() == 4 {
                let prefix = [obs.cases[0], obs.cases[1], obs.cases[2]];
                *out.entry(prefix).or_default() += 1;
            }
            let (u1, v1, _kc) = super::kaliski_jump::kaliski_step_uv(u, v);
            u = u1;
            v = v1;
        }
    }
    out
}

pub fn profile_exact3_bulk_core(seed: &[u8], n_inputs: usize, w: usize) -> HybridProtoSummary {
    let prefix_counts = collect_full_prefix_counts(seed, n_inputs, w);
    let full_windows: usize = prefix_counts.values().sum();
    let baseline3 = profile_baseline_iters(3);
    let baseline4 = profile_baseline_iters(4);
    let baseline_step3 = profile_single_baseline_iter(3);

    let mut rows = Vec::new();
    let mut weighted_exact3_ccx = 0.0;
    let mut weighted_exact3_cliff = 0.0;
    let mut min_exact3_ccx = u64::MAX;
    let mut max_exact3_ccx = 0u64;

    for (&prefix, &windows) in &prefix_counts {
        let exact3 = profile_exact_prefix(prefix);
        weighted_exact3_ccx += exact3.ccx as f64 * windows as f64;
        weighted_exact3_cliff += exact3.cliff as f64 * windows as f64;
        min_exact3_ccx = min_exact3_ccx.min(exact3.ccx);
        max_exact3_ccx = max_exact3_ccx.max(exact3.ccx);
        rows.push(PrefixCostRow {
            prefix,
            windows,
            exact3,
            staged_exact3: ProtoCost::default(),
            recompute_cleanup: ProtoCost::default(),
        });
    }
    weighted_exact3_ccx /= full_windows as f64;
    weighted_exact3_cliff /= full_windows as f64;

    let weighted_hybrid4_lower_bound_ccx = weighted_exact3_ccx + baseline_step3.ccx as f64;
    let weighted_hybrid4_lower_bound_cliff = weighted_exact3_cliff + baseline_step3.cliff as f64;
    let lower_bound_ccx_savings_vs_baseline4 = baseline4.ccx as f64 - weighted_hybrid4_lower_bound_ccx;
    let lower_bound_ccx_savings_pct_vs_baseline4 = 100.0 * lower_bound_ccx_savings_vs_baseline4 / baseline4.ccx as f64;

    let gt_compare = profile_gt_compare();
    let real_bulk3_forward = profile_real_bulk3_forward();
    let real_bulk3_specialized = profile_real_bulk3_specialized();
    let mut weighted_staged_exact3_ccx = 0.0;
    let mut weighted_staged_exact3_cliff = 0.0;
    let mut weighted_recompute_cleanup_ccx = 0.0;
    let mut weighted_recompute_cleanup_cliff = 0.0;

    rows.clear();
    for (&prefix, &windows) in &prefix_counts {
        let exact3 = profile_exact_prefix(prefix);
        let staged_exact3 = profile_staged_exact3(prefix);
        let recompute_cleanup = profile_recompute_cleanup(prefix);
        weighted_staged_exact3_ccx += staged_exact3.ccx as f64 * windows as f64;
        weighted_staged_exact3_cliff += staged_exact3.cliff as f64 * windows as f64;
        weighted_recompute_cleanup_ccx += recompute_cleanup.ccx as f64 * windows as f64;
        weighted_recompute_cleanup_cliff += recompute_cleanup.cliff as f64 * windows as f64;
        rows.push(PrefixCostRow { prefix, windows, exact3, staged_exact3, recompute_cleanup });
    }
    weighted_staged_exact3_ccx /= full_windows as f64;
    weighted_staged_exact3_cliff /= full_windows as f64;
    weighted_recompute_cleanup_ccx /= full_windows as f64;
    weighted_recompute_cleanup_cliff /= full_windows as f64;
    let total_keep_live_ccx = weighted_staged_exact3_ccx;
    let total_keep_live_cliff = weighted_staged_exact3_cliff;
    let total_recompute_ccx = weighted_staged_exact3_ccx + weighted_recompute_cleanup_ccx;
    let total_recompute_cliff = weighted_staged_exact3_cliff + weighted_recompute_cleanup_cliff;

    rows.sort_by(|a, b| b.windows.cmp(&a.windows).then_with(|| a.prefix.cmp(&b.prefix)));
    rows.truncate(12);

    HybridProtoSummary {
        full_windows,
        distinct_prefixes: prefix_counts.len(),
        baseline3,
        baseline4,
        baseline_step3,
        weighted_exact3_ccx,
        weighted_exact3_cliff,
        min_exact3_ccx,
        max_exact3_ccx,
        weighted_hybrid4_lower_bound_ccx,
        weighted_hybrid4_lower_bound_cliff,
        lower_bound_ccx_savings_vs_baseline4,
        lower_bound_ccx_savings_pct_vs_baseline4,
        gt_compare,
        weighted_staged_exact3_ccx,
        weighted_staged_exact3_cliff,
        weighted_recompute_cleanup_ccx,
        weighted_recompute_cleanup_cliff,
        total_keep_live_ccx,
        total_keep_live_cliff,
        total_recompute_ccx,
        total_recompute_cliff,
        real_bulk3_forward,
        real_bulk3_specialized,
        top_prefix_rows: rows,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact3_bulk_core_profile() {
        let s = profile_exact3_bulk_core(b"kaliski-hybrid-proto-seed-v1", 10_000, 8);
        eprintln!("=== Hybrid exact-3 bulk-core prototype profile (w=8) ===");
        eprintln!("full windows observed                  : {}", s.full_windows);
        eprintln!("distinct full 3-step prefixes         : {}", s.distinct_prefixes);
        eprintln!("baseline 3 iters                      : ccx={} cliff={} peak_qubits={}", s.baseline3.ccx, s.baseline3.cliff, s.baseline3.peak_qubits);
        eprintln!("baseline 4 iters                      : ccx={} cliff={} peak_qubits={}", s.baseline4.ccx, s.baseline4.cliff, s.baseline4.peak_qubits);
        eprintln!("baseline iter #3 only                 : ccx={} cliff={} peak_qubits={}", s.baseline_step3.ccx, s.baseline_step3.cliff, s.baseline_step3.peak_qubits);
        eprintln!("weighted exact-3 core                 : ccx={:.1} cliff={:.1}", s.weighted_exact3_ccx, s.weighted_exact3_cliff);
        eprintln!("exact-3 core ccx range                : min={} max={}", s.min_exact3_ccx, s.max_exact3_ccx);
        eprintln!("hybrid 4-step lower bound (no selector): ccx={:.1} cliff={:.1}", s.weighted_hybrid4_lower_bound_ccx, s.weighted_hybrid4_lower_bound_cliff);
        eprintln!("lower-bound ccx savings vs 4 iters    : {:.1} ({:.2}%)", s.lower_bound_ccx_savings_vs_baseline4, s.lower_bound_ccx_savings_pct_vs_baseline4);
        eprintln!("one 256-bit gt comparator             : ccx={} cliff={} peak_qubits={}", s.gt_compare.ccx, s.gt_compare.cliff, s.gt_compare.peak_qubits);
        eprintln!("weighted staged exact-3 (A false)     : ccx={:.1} cliff={:.1}", s.weighted_staged_exact3_ccx, s.weighted_staged_exact3_cliff);
        eprintln!("weighted cleanup recompute overhead   : ccx={:.1} cliff={:.1}", s.weighted_recompute_cleanup_ccx, s.weighted_recompute_cleanup_cliff);
        eprintln!("keep-live strategy total              : ccx={:.1} cliff={:.1}", s.total_keep_live_ccx, s.total_keep_live_cliff);
        eprintln!("recompute strategy total              : ccx={:.1} cliff={:.1}", s.total_recompute_ccx, s.total_recompute_cliff);
        eprintln!("real forward bulk3 keep-live proto    : ccx={} cliff={} peak_qubits={}", s.real_bulk3_forward.ccx, s.real_bulk3_forward.cliff, s.real_bulk3_forward.peak_qubits);
        eprintln!("real bulk3 specialized primitive      : ccx={} cliff={} peak_qubits={}", s.real_bulk3_specialized.ccx, s.real_bulk3_specialized.cliff, s.real_bulk3_specialized.peak_qubits);
        eprintln!("top full-window prefixes:");
        for row in &s.top_prefix_rows {
            eprintln!(
                "  {:>8} windows : {:<11}  exact3_ccx={:<6} staged3_ccx={:<6} cleanup_recmp_ccx={:<6}",
                row.windows,
                prefix_to_string(row.prefix),
                row.exact3.ccx,
                row.staged_exact3.ccx,
                row.recompute_cleanup.ccx,
            );
        }
        eprintln!("========================================================");

        assert_eq!(s.distinct_prefixes, 36);
        assert!(s.full_windows > 3_500_000);
        assert!(s.min_exact3_ccx < s.baseline3.ccx);
        assert!(s.max_exact3_ccx < s.baseline3.ccx);
        assert!(s.weighted_exact3_ccx < s.baseline3.ccx as f64);
        assert!(s.weighted_hybrid4_lower_bound_ccx < s.baseline4.ccx as f64);
        assert!(s.weighted_staged_exact3_ccx > s.weighted_exact3_ccx);
        assert!(s.weighted_recompute_cleanup_ccx > 0.0);
        assert!(s.total_recompute_ccx > s.total_keep_live_ccx);
        assert!(s.real_bulk3_forward.ccx > s.weighted_staged_exact3_ccx as u64);
        assert!(s.real_bulk3_specialized.ccx < s.baseline3.ccx);
    }
}
