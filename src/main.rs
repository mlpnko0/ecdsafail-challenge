//! Experiment harness: builds the point-addition circuit defined in
//! `point_add.rs`, runs it against the zenodo `Simulator` with random
//! test shots, and reports Toffoli / Clifford / qubit counts.
//!
//! Research-loop contract: ONLY `point_add.rs` is edited by the loop.
//! This file, `builder.rs`, `circuit.rs`, `sim.rs`, and
//! `weierstrass_elliptic_curve.rs` are harness and must not be touched.
//!
//! Attribution: `circuit.rs`, `sim.rs`, and `weierstrass_elliptic_curve.rs`
//! are reused verbatim from the `zkp_ecc` Zenodo project under CC BY 4.0.
//! See `NOTICE` at the repository root for details.

#[allow(dead_code)]
mod circuit;
#[allow(dead_code)]
mod sim;
#[allow(dead_code)]
mod weierstrass_elliptic_curve;
#[allow(dead_code)]
mod builder;
mod point_add;

use alloy_primitives::U256;
use builder::{Builder, Layout};
use circuit::{Op, OperationType, QubitOrBit, analyze_ops};
use sha3::{digest::{ExtendableOutput, Update, XofReader}, Shake256};
use sim::Simulator;
use std::fs::OpenOptions;
use std::io::Write;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};
use weierstrass_elliptic_curve::WeierstrassEllipticCurve;

// ─── secp256k1 parameters ──────────────────────────────────────────────────

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap(),
        a: U256::from(0),
        b: U256::from(7),
        gx: U256::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap(),
        gy: U256::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap(),
        order: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap(),
    }
}

// ─── Test runner ───────────────────────────────────────────────────────────

const NUM_TESTS: usize = 64;

/// Build the inverse op list for a forward-then-reverse reversibility check.
///
/// All Clifford+Toffoli gates used by our circuit (X, Z, CX, CZ, CCX, CCZ,
/// Swap) are self-inverse, so reversing just means reversing the op order.
/// Metadata ops (Register / AppendToRegister / DebugPrint) are skipped.
/// `R` is treated as an assertion (identity on zero qubits) and also
/// skipped — `strict_apply` already enforces its precondition on the
/// forward pass.
/// `PushCondition` and `PopCondition` are mirrored (push↔pop) so a
/// conditional block in the forward pass is still scoped correctly in
/// the reverse pass — not used by current point_add.rs, but kept for
/// future correctness.
/// `Hmr`, `BitInvert`, `BitStore0`, `BitStore1`, `Neg` are errors: they
/// are either non-reversible or should not appear in a reversible
/// circuit; panic to fail loud.
fn invert_ops(ops: &[Op]) -> Vec<Op> {
    let mut out = Vec::with_capacity(ops.len());
    for op in ops.iter().rev() {
        match op.kind {
            OperationType::X
            | OperationType::Z
            | OperationType::CX
            | OperationType::CZ
            | OperationType::CCX
            | OperationType::CCZ
            | OperationType::Swap => out.push(*op),
            OperationType::PushCondition => {
                let mut o = *op;
                o.kind = OperationType::PopCondition;
                out.push(o);
            }
            OperationType::PopCondition => {
                let mut o = *op;
                o.kind = OperationType::PushCondition;
                out.push(o);
            }
            OperationType::R
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {
                // skip: assertions and metadata
            }
            OperationType::Hmr
            | OperationType::BitInvert
            | OperationType::BitStore0
            | OperationType::BitStore1
            | OperationType::Neg => {
                panic!(
                    "invert_ops: op kind {:?} is not allowed in a reversible circuit",
                    op.kind
                );
            }
        }
    }
    out
}

/// Apply `ops` to `sim` while enforcing that every `R` op targets a qubit
/// whose value is 0 on every live shot. This turns `R` back into a
/// reversibility assertion (the `Builder::free_qubit` contract) rather
/// than a free "forget this ancilla" primitive that happens to work
/// because the sim zeros unconditionally.
///
/// `push_condition` / `pop_condition` are not used in our circuit, so
/// we check the raw qubit value across all 64 shots. If ever needed,
/// this would need to thread the live-condition mask.
fn strict_apply<'a, R: sha3::digest::XofReader>(
    sim: &mut Simulator<'a, R>,
    ops: &[Op],
) -> Result<(), String> {
    let mut start = 0usize;
    for (i, op) in ops.iter().enumerate() {
        if op.kind == OperationType::R {
            if start < i {
                sim.apply(&ops[start..i]);
            }
            let v = sim.qubit(op.q_target);
            if v != 0 {
                return Err(format!(
                    "REVERSIBILITY VIOLATION at op #{i}: R targets qubit {:?} \
                     but its value across 64 shots is {v:#018x} (nonzero). \
                     Builder::free_qubit requires the ancilla to already be |0⟩ — \
                     proper reversible uncomputation is missing.",
                    op.q_target
                ));
            }
            sim.apply(std::slice::from_ref(op));
            start = i + 1;
        }
    }
    if start < ops.len() {
        sim.apply(&ops[start..]);
    }
    Ok(())
}

fn run_tests(ops: &[Op], layout_regs: &[Vec<QubitOrBit>], total_qubits: u32, num_bits: u32)
    -> (bool, f64, f64, u64, u64, usize)
{
    let curve = secp256k1();
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-baseline-seed-v1");
    let mut xof = hasher.finalize_xof();

    // Generate random target/offset points as k*G.
    let mut targets = Vec::with_capacity(NUM_TESTS);
    let mut offsets = Vec::with_capacity(NUM_TESTS);
    let mut expected = Vec::with_capacity(NUM_TESTS);
    for _ in 0..NUM_TESTS {
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let t = curve.mul(curve.gx, curve.gy, k1);
        let o = curve.mul(curve.gx, curve.gy, k2);
        // Avoid the doubling / inverse-pair cases the baseline doesn't handle.
        if t.0 == o.0 { continue; }
        if t.0.is_zero() && t.1.is_zero() { continue; }
        if o.0.is_zero() && o.1.is_zero() { continue; }
        let e = curve.add(t.0, t.1, o.0, o.1);
        targets.push(t);
        offsets.push(o);
        expected.push(e);
    }
    let n = targets.len();

    let mut rng_hasher = Shake256::default();
    rng_hasher.update(b"quantum_ecc-sim-rng-v1");
    let mut rng = rng_hasher.finalize_xof();

    let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut rng);
    let mut ok = true;

    let mut got = vec![(U256::ZERO, U256::ZERO); n];

    // Precompute inverse circuit for the reversibility identity test.
    let inv_ops = invert_ops(ops);

    const BATCH: usize = 64;
    let num_batches = (n + BATCH - 1) / BATCH;
    for batch in 0..num_batches {
        let bs = BATCH.min(n - batch * BATCH);
        sim.clear_for_shot();
        for shot in 0..bs {
            let i = batch * BATCH + shot;
            sim.set_register(&layout_regs[0], targets[i].0, shot);
            sim.set_register(&layout_regs[1], targets[i].1, shot);
            sim.set_register(&layout_regs[2], offsets[i].0, shot);
            sim.set_register(&layout_regs[3], offsets[i].1, shot);
        }

        // Snapshot every qubit state before the forward pass so we can
        // verify forward∘reverse = identity at the end.
        let snapshot: Vec<u64> = (0..total_qubits)
            .map(|q| sim.qubit(circuit::QubitId(q)))
            .collect();

        if let Err(e) = strict_apply(&mut sim, ops) {
            eprintln!("\n!! {e}");
            return (false, 0.0, 0.0, 0, 0, n);
        }
        for shot in 0..bs {
            let i = batch * BATCH + shot;
            let gx = sim.get_register(&layout_regs[0], shot);
            let gy = sim.get_register(&layout_regs[1], shot);
            got[i] = (gx, gy);
            if gx != expected[i].0 || gy != expected[i].1 {
                ok = false;
            }
        }

        // Reversibility identity check: apply the inverse op list and
        // verify every qubit returns to its pre-forward snapshot. This
        // catches any form of information loss (dirty R, bad uncompute,
        // missing inverse, etc.) regardless of whether the forward
        // result happened to look right.
        sim.apply(&inv_ops);
        for q in 0..total_qubits {
            let before = snapshot[q as usize];
            let after = sim.qubit(circuit::QubitId(q));
            if before != after {
                eprintln!(
                    "\n!! REVERSIBILITY FAILURE: qubit {} differs after forward∘reverse\n   before: {:#018x}\n   after : {:#018x}",
                    q, before, after
                );
                ok = false;
                break;
            }
        }
        if !ok { break; }
    }

    println!("  test points:");
    for i in 0..n {
        let mark = if got[i] == expected[i] { "OK  " } else { "FAIL" };
        println!("    [{i:02}] {mark}");
        println!("         T   =({:#x}, {:#x})", targets[i].0, targets[i].1);
        println!("         O   =({:#x}, {:#x})", offsets[i].0, offsets[i].1);
        println!("         got =({:#x}, {:#x})", got[i].0, got[i].1);
        println!("         exp =({:#x}, {:#x})", expected[i].0, expected[i].1);
    }

    let denom = n.max(1) as f64;
    let avg_cliff = sim.stats.clifford_gates as f64 / denom;
    let avg_tof = sim.stats.toffoli_gates as f64 / denom;
    (ok, avg_cliff, avg_tof, sim.stats.toffoli_gates, sim.stats.clifford_gates, n)
}

fn parse_note() -> String {
    let mut args = std::env::args().skip(1);
    let mut note = String::new();
    while let Some(a) = args.next() {
        if a == "--note" {
            if let Some(v) = args.next() { note = v; }
        } else if let Some(rest) = a.strip_prefix("--note=") {
            note = rest.to_string();
        }
    }
    note.replace('\t', " ").replace('\n', " ")
}

fn git_commit_short() -> String {
    Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| if o.status.success() { Some(String::from_utf8_lossy(&o.stdout).trim().to_string()) } else { None })
        .unwrap_or_else(|| "nogit".to_string())
}

fn append_results_row(
    correct: &str,
    avg_tof: f64,
    avg_cliff: f64,
    qubits: u32,
    ops_len: usize,
    note: &str,
) {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let commit = git_commit_short();
    let row = format!(
        "{ts}\t{commit}\t{avg_tof:.3}\t{avg_cliff:.3}\t{qubits}\t{ops_len}\t{correct}\t{note}\n"
    );
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/results.tsv");
    match OpenOptions::new().create(true).append(true).open(path) {
        Ok(mut f) => {
            if let Err(e) = f.write_all(row.as_bytes()) {
                eprintln!("warning: failed to write results.tsv: {e}");
            }
        }
        Err(e) => eprintln!("warning: failed to open results.tsv: {e}"),
    }
}

fn main() {
    let note = parse_note();
    println!("=== quantum_ecc: secp256k1 point addition baseline ===\n");
    let curve = secp256k1();

    println!("-- building circuit --");
    let mut builder = Builder::new();
    let layout = point_add::build(&mut builder);
    let ops = builder.ops.clone();
    let peak = builder.peak_qubits();
    let total_alloc = builder.total_qubits();

    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter().copied());

    // Sanity-check layout matches zenodo's program interface.
    assert!(regs.len() == 4, "expected 4 registers (target_x, target_y, offset_x, offset_y); got {}", regs.len());
    for (i, r) in regs.iter().enumerate() {
        assert_eq!(r.len(), 256, "register {i} should be 256 wide, got {}", r.len());
    }
    for q in &regs[0] { assert!(matches!(q, QubitOrBit::Qubit(_)), "register 0 must be qubits"); }
    for q in &regs[1] { assert!(matches!(q, QubitOrBit::Qubit(_)), "register 1 must be qubits"); }
    for q in &regs[2] { assert!(matches!(q, QubitOrBit::Bit(_)),   "register 2 must be bits"); }
    for q in &regs[3] { assert!(matches!(q, QubitOrBit::Bit(_)),   "register 3 must be bits"); }
    let _ = layout; // layout IDs are implicit in ordering

    println!("  total ops : {}", ops.len());
    println!("  qubits    : {} (peak live {}, total alloc {})", total_qubits, peak, total_alloc);
    println!("  bits      : {}", num_bits);

    println!("\n-- running correctness tests --");
    let (ok, avg_cliff, avg_tof, tot_tof, tot_cliff, n_shots) = run_tests(&ops, &regs, total_qubits, num_bits);
    if !ok {
        println!("\n!! correctness FAILED");
        append_results_row("FAIL", avg_tof, avg_cliff, total_qubits, ops.len(), &note);
        std::process::exit(1);
    }
    println!("  all {} shots OK", NUM_TESTS);

    println!("\n=== circuit metrics (secp256k1, n=256) ===");
    println!("  avg executed Toffoli  : {:.3}", avg_tof);
    println!("  avg executed Clifford : {:.3}", avg_cliff);
    println!("  total Toffoli (sum)   : {} over {} shots", tot_tof, n_shots);
    println!("  total Clifford (sum)  : {}", tot_cliff);
    println!("  emitted ops           : {}", ops.len());
    println!("  qubits                : {}", total_qubits);

    append_results_row("OK", avg_tof, avg_cliff, total_qubits, ops.len(), &note);

    println!("\n=== experiment OK ===");
}
