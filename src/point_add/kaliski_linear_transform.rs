//! Ground-up structural probe: use Kaliski's coefficient update as a linear
//! transform on the *data* y-register instead of treating it as disposable
//! ancilla.
//!
//! This is analysis-only (`#[cfg(test)]` module imported from `mod.rs`). It
//! tests a possible 600-scratch architecture:
//!
//! - keep `tx = dx` as the preserved x-difference,
//! - use `ty` as Kaliski's coefficient register `s`, initialized to `dy`,
//! - run a canonical-mod-p coefficient version of Kaliski.
//!
//! If this worked naively, the forward Kaliski would turn `ty=dy` into
//! `s=0` while `r = raw_inv(dx) * dy`, i.e. the scaled slope. Then Kaliski's
//! backward coefficient transform might be used to write the final `Ry` into
//! `ty` without a second inversion. The tests below verify the linear algebra
//! and isolate the remaining obstruction.

#![cfg(test)]
#![allow(dead_code)]

use alloy_primitives::U256;
use sha3::{digest::{ExtendableOutput, Update, XofReader}, Shake128};

use super::SECP256K1_P;

const ITERS: usize = 407;

fn random_element(seed: u64) -> U256 {
    let mut h = Shake128::default();
    h.update(&seed.to_le_bytes());
    let mut reader = h.finalize_xof();
    loop {
        let mut buf = [0u8; 32];
        reader.read(&mut buf);
        let v = U256::from_be_bytes(buf);
        if v != U256::ZERO && v < SECP256K1_P {
            return v;
        }
    }
}

#[inline]
fn sub_mod(a: U256, b: U256, p: U256) -> U256 {
    let (r, borrow) = a.overflowing_sub(b);
    if borrow { r.wrapping_add(p) } else { r }
}

#[inline]
fn neg_mod(a: U256, p: U256) -> U256 {
    if a.is_zero() { a } else { p.wrapping_sub(a) }
}

#[inline]
fn add_mod(a: U256, b: U256, p: U256) -> U256 {
    a.add_mod(b, p)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Branch {
    a_swap: bool,
    add: bool,
}

#[derive(Clone, Copy, Debug)]
struct LinState {
    u: U256,
    v: U256,
    r: U256,
    s: U256,
    f: u8,
}

fn limbs(x: U256) -> [u64; 4] {
    *x.as_limbs()
}

/// The branch sequence depends only on `(u,v,f)`, not on the coefficient
/// values, so it can be separated from the coefficient linear transform.
fn branch_sequence(dx: U256, iters: usize) -> Vec<Branch> {
    let p = SECP256K1_P;
    let mut u = p;
    let mut v = dx;
    let mut f = 1u8;
    let mut out = Vec::with_capacity(iters);
    for _ in 0..iters {
        let mut m = 0u8;
        if f == 1 && v == U256::ZERO { m ^= 1; }
        f ^= m;

        let u0 = if u.bit(0) { 1u8 } else { 0u8 };
        let v0 = if v.bit(0) { 1u8 } else { 0u8 };
        let mut a = 0u8;
        if f == 1 && u0 == 0 { a ^= 1; }
        if f == 1 && u0 == 1 && v0 == 0 { m ^= 1; }
        let b = a ^ m;
        let gt = if u > v { 1u8 } else { 0u8 };
        let delta = (f & gt) & (1 ^ b);
        a ^= delta;
        m ^= delta;
        let add = (f & (1 ^ b)) == 1;
        let a_swap = a == 1;
        out.push(Branch { a_swap, add });

        if a_swap { core::mem::swap(&mut u, &mut v); }
        if add { v = v.wrapping_sub(u); }
        v >>= 1;
        if a_swap { core::mem::swap(&mut u, &mut v); }
        let _ = m;
    }
    out
}

/// Apply the coefficient-side transform with canonical mod-p arithmetic.
/// This is *not* exactly the current circuit's noncanonical `s=p` sentinel;
/// it is the modified architecture needed if `s` is a data register like `dy`.
fn apply_coeffs(seq: &[Branch], mut r: U256, mut s: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    for br in seq {
        if br.a_swap { core::mem::swap(&mut r, &mut s); }
        if br.add { s = add_mod(s, r, p); }
        r = add_mod(r, r, p);
        if br.a_swap { core::mem::swap(&mut r, &mut s); }
    }
    (r, s)
}

fn pow2_mod(e: usize) -> U256 {
    let mut r = U256::from(1u64);
    for _ in 0..e {
        r = add_mod(r, r, SECP256K1_P);
    }
    r
}

fn step_linear_canonical(st: &mut LinState) -> Branch {
    step_linear_canonical_with_flags(st).0
}

fn step_linear_canonical_with_flags(st: &mut LinState) -> (Branch, u8, u8) {
    let mut m = 0u8;
    if st.f == 1 && st.v == U256::ZERO { m ^= 1; }
    st.f ^= m;

    let u0 = if st.u.bit(0) { 1u8 } else { 0u8 };
    let v0 = if st.v.bit(0) { 1u8 } else { 0u8 };
    let mut a = 0u8;
    if st.f == 1 && u0 == 0 { a ^= 1; }
    if st.f == 1 && u0 == 1 && v0 == 0 { m ^= 1; }
    let b = a ^ m;
    let gt = if st.u > st.v { 1u8 } else { 0u8 };
    let delta = (st.f & gt) & (1 ^ b);
    a ^= delta;
    m ^= delta;
    let br = Branch { a_swap: a == 1, add: (st.f & (1 ^ b)) == 1 };

    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
    }
    if br.add {
        st.v = st.v.wrapping_sub(st.u);
        st.s = add_mod(st.s, st.r, SECP256K1_P);
    }
    st.v >>= 1;
    st.r = add_mod(st.r, st.r, SECP256K1_P);
    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
    }
    (br, a, m)
}

#[test]
fn coefficient_transform_shape() {
    let p = SECP256K1_P;
    let scale = pow2_mod(ITERS);
    for seed in 1..50u64 {
        let dx = random_element(seed);
        let seq = branch_sequence(dx, ITERS);
        let (a, c) = apply_coeffs(&seq, U256::from(1u64), U256::ZERO);
        let (k, d) = apply_coeffs(&seq, U256::ZERO, U256::from(1u64));

        // Empirical theorem for the canonical coefficient transform T(dx):
        //      T = [[a(dx), k(dx)], [dx, 0]]
        // with k(dx) * dx = -2^ITERS mod p.
        assert_eq!(c, dx, "lower-left coefficient is exactly dx");
        assert_eq!(d, U256::ZERO, "lower-right coefficient is zero");
        assert_eq!(k.mul_mod(dx, p), neg_mod(scale, p), "k is the raw inverse scale");
        assert_eq!(k.mul_mod(c, p), neg_mod(scale, p), "determinant relation");
        let _ = a;
    }
}

#[test]
fn single_coefficient_pair_cannot_preserve_x_and_expose_quotient_by_constant_tag() {
    // Try the most tempting one-pair DIV rescue.  Set r0=ρ (nonzero constant)
    // so the lower output s=ρ*x preserves the denominator while seed
    // s0=y+β.  The upper output is
    //     r = k*y + (ρ*a + β*k).
    // If the parenthesized contaminant were a known constant, one coefficient
    // pair would simultaneously expose y/x and keep x, fitting the ~600q
    // target.  This requires an affine relation ρ*a(x)+β*k(x)=C across all x.
    // Three sampled transforms already make (a,k,1) non-collinear, killing all
    // constant-tag/constant-r0 variants of this rescue.
    let p = SECP256K1_P;
    let mut pts = Vec::new();
    for seed in 1..=3u64 {
        let x = random_element(seed);
        let seq = branch_sequence(x, ITERS);
        let (a, lower) = apply_coeffs(&seq, U256::from(1u64), U256::ZERO);
        let (k, zero) = apply_coeffs(&seq, U256::ZERO, U256::from(1u64));
        assert_eq!(lower, x);
        assert_eq!(zero, U256::ZERO);
        pts.push((a, k));
    }
    let (a0, k0) = pts[0];
    let (a1, k1) = pts[1];
    let (a2, k2) = pts[2];
    let da10 = sub_mod(a1, a0, p);
    let dk10 = sub_mod(k1, k0, p);
    let da20 = sub_mod(a2, a0, p);
    let dk20 = sub_mod(k2, k0, p);
    let det = sub_mod(da10.mul_mod(dk20, p), da20.mul_mod(dk10, p), p);
    eprintln!("constant-tag coefficient-pair relation determinant = {det:#x}");
    assert!(!det.is_zero(), "sampled (a,k) were affine-collinear; constant-tag DIV rescue may exist");
}

#[test]
fn dx_tagged_seed_recovers_division_with_negligible_exception() {
    // Approximate tolerance reopens the self-cleaning DIV route. Seed the
    // coefficient with (y + x) instead of y. Then
    //   T(x)*(0, y+x) = (k*y + k*x, 0) = (k*y - 2^ITERS, 0)
    // because k*x = -2^ITERS. Adding the known scale recovers k*y, and a
    // known rescale gives y/x. The only zero-coefficient exceptional set is
    // y = -x, probability ≈ 1/p for random field inputs.
    let p = SECP256K1_P;
    let scale = pow2_mod(ITERS);
    let scale_inv = scale.inv_mod(p).unwrap();
    for seed in 1..100u64 {
        let x = random_element(seed);
        let y = random_element(seed + 10_000);
        let tagged = add_mod(y, x, p);
        assert_ne!(tagged, U256::ZERO, "random sample hit y=-x exceptional set");
        let seq = branch_sequence(x, ITERS);
        let (r_tagged, s_out) = apply_coeffs(&seq, U256::ZERO, tagged);
        assert_eq!(s_out, U256::ZERO);
        let k_y = add_mod(r_tagged, scale, p); // r + 2^ITERS = k*y
        let quotient = neg_mod(k_y, p).mul_mod(scale_inv, p);
        assert_eq!(quotient, y.mul_mod(x.inv_mod(p).unwrap(), p));
    }
}

#[test]
fn stored_a_and_m_bits_recover_branch_pair() {
    // If we abandon qrisp's full inverse coefficient `(r,s)` sentinel, one
    // plausible branch-only cleanup stores the final swap bit `a` in addition
    // to the existing `m_hist`. The per-step add bit is then not independent:
    // for active steps, add = !(a xor m); after termination f=0 forces add=0.
    // This does not solve the 600-scratch target by itself (it still stores
    // history), but it validates the next branch-only circuit scaffold.
    for seed in 1..200u64 {
        let mut st = LinState {
            u: SECP256K1_P,
            v: random_element(seed),
            r: U256::ZERO,
            s: add_mod(random_element(seed + 10_000), random_element(seed), SECP256K1_P),
            f: 1,
        };
        for _ in 0..ITERS {
            let (br, a, m) = step_linear_canonical_with_flags(&mut st);
            assert_eq!(br.a_swap, a == 1);
            let recovered_add = st.f == 1 && ((a ^ m) == 0);
            assert_eq!(br.add, recovered_add, "add should be recoverable from stored a,m and post f");
        }
    }
}

#[test]
fn dy_seeded_forward_computes_scaled_slope_and_zeroes_s() {
    let p = SECP256K1_P;
    let scale = pow2_mod(ITERS);
    for seed in 1..50u64 {
        let dx = random_element(seed);
        let dy = random_element(seed + 10_000);
        let seq = branch_sequence(dx, ITERS);
        let (r, s) = apply_coeffs(&seq, U256::ZERO, dy);
        let expect = neg_mod(scale, p)
            .mul_mod(dy, p)
            .mul_mod(dx.inv_mod(p).unwrap(), p);
        assert_eq!(r, expect, "r = raw_inv(dx) * dy = scaled slope");
        assert_eq!(s, U256::ZERO, "s/ty is consumed to zero in canonical form");
    }
}

#[test]
fn end_state_needs_coefficient_registers_to_recover_branch() {
    // A forward-only low-qubit DIV would like to run Kaliski without storing
    // m_hist. That requires each iteration's branch bit to be uncomputed from
    // the updated live state. This diagnostic separates two facts:
    //   1. denominator state alone (u,v,f) is NOT enough; many collisions occur.
    //   2. full coefficient state (u,v,r,s,f) WAS enough on this sample set.
    // So a self-cleaning DIV, if it exists, must use the coefficient registers
    // in the branch-recovery predicate; a tiny parity/comparator fingerprint is
    // not enough.
    use std::collections::HashMap;

    let mut denom_seen: HashMap<([u64; 4], [u64; 4], u8), Branch> = HashMap::new();
    let mut full_seen: HashMap<([u64; 4], [u64; 4], [u64; 4], [u64; 4], u8), Branch> = HashMap::new();
    let mut denom_conflicts = 0usize;
    let mut full_conflicts = 0usize;

    for seed in 1..=200u64 {
        let mut st = LinState {
            u: SECP256K1_P,
            v: random_element(seed),
            r: U256::ZERO,
            s: random_element(seed + 10_000),
            f: 1,
        };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let dk = (limbs(st.u), limbs(st.v), st.f);
            if let Some(prev) = denom_seen.insert(dk, br) {
                if prev != br { denom_conflicts += 1; }
            }
            let fk = (limbs(st.u), limbs(st.v), limbs(st.r), limbs(st.s), st.f);
            if let Some(prev) = full_seen.insert(fk, br) {
                if prev != br { full_conflicts += 1; }
            }
        }
    }

    assert!(denom_conflicts > 0, "denominator-only end-state unexpectedly recovered branches");
    assert_eq!(full_conflicts, 0, "full end-state branch recovery collided in samples");
}

#[test]
fn bilinear_invariant_does_not_recover_inverse_branch() {
    // The obvious algebraic invariant of the coefficient transform is
    //     r*v + s*u = 0 (mod p)
    // starting from (u,v,r,s)=(p,x,0,tag). Unfortunately it is preserved by
    // almost all locally valid inverse candidates, so it does not provide the
    // cheap self-cleaning branch predicate we need.
    let p = SECP256K1_P;
    let inv2 = U256::from(2u64).inv_mod(p).unwrap();
    let mut ambiguous = 0usize;
    let mut total = 0usize;

    for seed in 1..=200u64 {
        let x = random_element(seed);
        let y = random_element(seed + 10_000);
        let mut st = LinState { u: p, v: x, r: U256::ZERO, s: add_mod(y, x, p), f: 1 };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            if st.f == 0 { continue; }
            let mut survivors = 0usize;
            let candidates = [
                // (case_is_true, pre_u, pre_v, pre_r, pre_s)
                (!br.a_swap && !br.add, st.u, st.v << 1, st.r.mul_mod(inv2, p), st.s),
                ( br.a_swap && !br.add, st.u << 1usize, st.v, st.r, st.s.mul_mod(inv2, p)),
                ( br.a_swap &&  br.add, (st.u << 1usize).wrapping_add(st.v), st.v, sub_mod(st.r, st.s.mul_mod(inv2, p), p), st.s.mul_mod(inv2, p)),
                (!br.a_swap &&  br.add, st.u, (st.v << 1usize).wrapping_add(st.u), st.r.mul_mod(inv2, p), sub_mod(st.s, st.r.mul_mod(inv2, p), p)),
            ];
            for (_is_true, pu, pv, pr, ps) in candidates {
                let branch_valid = if pu.bit(0) == false {
                    // U-even candidate.
                    true
                } else if pv.bit(0) == false {
                    // V-even candidate.
                    true
                } else {
                    // Odd/odd candidate; either ordering is locally valid.
                    true
                };
                let invariant = add_mod(pr.mul_mod(pv % p, p), ps.mul_mod(pu % p, p), p) == U256::ZERO;
                if branch_valid && invariant {
                    survivors += 1;
                }
            }
            if survivors > 1 { ambiguous += 1; }
            total += 1;
        }
    }
    let frac = ambiguous as f64 / total as f64;
    assert!(frac > 0.90, "bilinear invariant unexpectedly disambiguated branches: ambiguous={frac}");
}

#[test]
fn low_bit_end_state_branch_classifier_is_not_approx_good_enough() {
    // Approximate incorrectness reopens rare exceptional sets, but it does not
    // make a crude local branch predicate viable. Train a best-majority lookup
    // table from low bits of the end-state registers, then test on disjoint
    // samples. Even with coefficient registers included, the error is huge.
    use std::collections::HashMap;

    type Key = (u16, u16, u16, u16, u8);
    const LOW_BITS: u32 = 3;
    let mask = (1u64 << LOW_BITS) - 1;
    let key_of = |st: &LinState| -> Key {
        (
            (st.u.as_limbs()[0] & mask) as u16,
            (st.v.as_limbs()[0] & mask) as u16,
            (st.r.as_limbs()[0] & mask) as u16,
            (st.s.as_limbs()[0] & mask) as u16,
            st.f,
        )
    };

    let mut counts: HashMap<Key, [usize; 4]> = HashMap::new();
    let idx = |br: Branch| -> usize { (br.a_swap as usize) * 2 + (br.add as usize) };

    for seed in 1..=120u64 {
        let mut st = LinState { u: SECP256K1_P, v: random_element(seed), r: U256::ZERO, s: random_element(seed + 10_000), f: 1 };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let k = key_of(&st);
            counts.entry(k).or_insert([0; 4])[idx(br)] += 1;
        }
    }

    let mut table: HashMap<Key, usize> = HashMap::new();
    for (k, c) in counts {
        let mut best_i = 0usize;
        let mut best_c = 0usize;
        for (i, &v) in c.iter().enumerate() {
            if v > best_c { best_c = v; best_i = i; }
        }
        table.insert(k, best_i);
    }

    let mut wrong = 0usize;
    let mut total = 0usize;
    for seed in 10_001..=10_120u64 {
        let mut st = LinState { u: SECP256K1_P, v: random_element(seed), r: U256::ZERO, s: random_element(seed + 10_000), f: 1 };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let k = key_of(&st);
            // All 3-bit keys are present in the train set; fallback is arbitrary.
            let pred = table.get(&k).copied().unwrap_or(0);
            if pred != idx(br) { wrong += 1; }
            total += 1;
        }
    }
    let err_rate = wrong as f64 / total as f64;
    assert!(err_rate > 0.50, "low-bit branch classifier unexpectedly good: err={err_rate}");
}

#[test]
fn zero_coefficient_seed_loses_branch_information() {
    // Exact DIV must also handle y=0 (or any value making the coefficient
    // channel uninformative). With r=s=0, full state collapses to the
    // denominator state, and branch recovery collides. Therefore any
    // self-cleaning forward-only Kaliski needs either an additional nonzero
    // tag mixed into the coefficient state or a branch predicate independent
    // of the coefficient scalar.
    use std::collections::HashMap;

    let mut seen: HashMap<([u64; 4], [u64; 4], [u64; 4], [u64; 4], u8), Branch> = HashMap::new();
    let mut conflicts = 0usize;
    for seed in 1..=200u64 {
        let mut st = LinState {
            u: SECP256K1_P,
            v: random_element(seed),
            r: U256::ZERO,
            s: U256::ZERO,
            f: 1,
        };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let key = (limbs(st.u), limbs(st.v), limbs(st.r), limbs(st.s), st.f);
            if let Some(prev) = seen.insert(key, br) {
                if prev != br { conflicts += 1; }
            }
        }
    }
    assert!(conflicts > 0, "zero coefficient seed unexpectedly preserved branch information");
}

#[test]
fn backward_write_condition_for_ry() {
    // If the coefficient transform is T=[[a,k],[dx,0]], then to have the
    // backward pass finish with `(r_initial=0, s_initial=Ry)`, the final
    // coefficient pair before backward MUST be T*(0,Ry) = (k*Ry, 0).
    // Starting from dy-seeded forward gives (k*dy, 0). So the structural
    // task is exactly to add k*(Ry-dy) into r, while s remains zero.
    // This test records the identity on random field values. It is not a
    // proof of impossibility; it is the crisp algebraic subproblem.
    let p = SECP256K1_P;
    for seed in 1..50u64 {
        let dx = random_element(seed);
        let dy = random_element(seed + 10_000);
        let ry = random_element(seed + 20_000);
        let seq = branch_sequence(dx, ITERS);
        let (k, _) = apply_coeffs(&seq, U256::ZERO, U256::from(1u64));
        let (r_dy, s_dy) = apply_coeffs(&seq, U256::ZERO, dy);
        let (r_ry, s_ry) = apply_coeffs(&seq, U256::ZERO, ry);
        assert_eq!(s_dy, U256::ZERO);
        assert_eq!(s_ry, U256::ZERO);
        assert_eq!(sub_mod(r_ry, r_dy, p), k.mul_mod(sub_mod(ry, dy, p), p));
    }
}
