# Research notes — inversion moonshots inside `src/point_add/`

Session: 2026-04-22 (continued, moonshot-only work).

This file keeps all moonshot literature / classical-analysis work under
`src/point_add/`, per the current scope rules.

## Deliverable 1 (classical B-Y on secp256k1) — confirmed

Implemented classical `divstep2` reference and modular-inverse recovery in
`src/point_add/by.rs`, then ran a 10,000-input secp256k1 survey.

Results:

| metric | value |
|---|---|
| theoretical bound `⌈(49·256 + 57)/17⌉` | 742 |
| observed minimum iters | 502 |
| observed maximum iters | 567 |
| observed mean iters | 531.01 |
| max `|δ|` observed | 20 |
| modinv matches (vs Fermat) | 10,000 / 10,000 |

Interpretation:
- The BY safegcd upper bound is pessimistic by ~24% on secp256k1 inputs.
- However, this is **not enough** to save plain B-Y: the per-iter reversible
  cost is still too high relative to Kaliski.

## Deliverable 2 (algorithm-space survey) — corrected final version

### 1. Kaliski almost-inverse (baseline)
- Classical ref: Burton S. Kaliski Jr., “The Montgomery inverse and its
  applications,” IEEE Trans. Computers 44(8), 1995.
- Quantum / reversible refs:
  - Roetteler–Naehrig–Svore–Lauter 2017, arXiv:1706.06752.
  - Häner–Roetteler–Soeken 2020, arXiv:2001.09580 / ePrint 2020/077.
- Iterations in our tuned circuit: 399.
- Measured per-iter reversible cost: ~2180 CCX.
- Per-pass cost: ~1.81M CCX.

### 2. Bernstein–Yang divstep2 (w = 1)
- Ref: Bernstein–Yang 2019, ePrint 2019/266.
- Reversible implementation: unpublished / would be novel.
- Empirical iterations on secp256k1: max 567, mean 531.
- Per-iter reversible estimate: 10–12n CCX.
- Conclusion: still worse than Kaliski.

### 3. Bernstein–Yang jumpdivsteps2 (w > 1)
- Ref: Bernstein–Yang 2019, Figure 10.2 / §10.
- Reversible implementation: unpublished / would be novel.

#### 3a. Corrected matrix-growth result
A previous version of the jump survey undercounted the scaled transition
matrix. After fixing it, the 100,000-sample survey now shows the **full
scaled** transition matrices do hit the theoretical `2^w` growth.

Corrected survey over 100,000 random low-word states:

| w | max observed `|entry|` | max log2 | mean log2 | theoretical max log2 |
|---|---:|---:|---:|---:|
| 4  | 16    | 4.00  | 2.03 | 4  |
| 8  | 256   | 8.00  | 4.28 | 8  |
| 12 | 4096  | 12.00 | 6.34 | 12 |
| 16 | 65536 | 16.00 | 8.19 | 16 |

Interpretation:
- The **maximum** entry size really does hit the full `2^w` growth.
- So a faithful reversible matrix-apply must still handle `w`-bit classical
  coefficients.
- That restores the pessimistic reversible cost model: batching by `w` does
  not automatically beat Kaliski.

#### 3b. Exact matrix-family compression result
Even if entries hit `2^w`, a quantum QROM implementation might still benefit
if the number of **distinct** transition matrices is tiny compared to the raw
state space. I measured this exactly for all low-word states with
`delta ∈ [-20, 20]`, odd `f_low`, and arbitrary `g_low`.

Results:

| w | total states | distinct matrices | compression factor |
|---|---:|---:|---:|
| 4 | 5,248 | 656 | 8× |
| 6 | 83,968 | 2,624 | 32× |
| 8 | 1,343,488 | 10,496 | 128× |

Pattern:
- compression factor = `2^(w−1)` exactly on the observed range.
- equivalently, distinct matrix count appears to scale like `2^(w+2)`.

This does **not** rescue full jumped B-Y by itself, but it is a strong sign
that *compressed local transition classes* are real and exploitable.

#### 3c. Updated verdict on jumped B-Y
Full jumped B-Y still looks too expensive as a drop-in replacement, because:
- matrix entries hit the full `2^w` growth,
- full coefficient tracking would still need to carry those `w`-bit entries,
- cleanup is all-new machinery.

But the compression result changes the local-batching story.

### 4. Montgomery inverse (Savaş–Koç)
- Classical ref: Savaş–Koç 2000, “The Montgomery modular inverse revisited.”
- Quantum / reversible refs: effectively same family as RNSL/HRSL Kaliski.
- Conclusion: not a distinct win over Kaliski in our setting.

### 5. Lehmer-style GCDs
- Classical refs: Lehmer 1938; Jebelean 1993.
- Reversible implementation: unpublished / novel.
- Main issue: runtime matrix selection depends on quantum data, so a faithful
  reversible implementation needs a QROM keyed by top bits. No concrete,
  literature-backed reversible cost win established yet.
- Still potentially interesting as novel research, but now less grounded than
  a compressed Kaliski-local batching route, because we have exact empirical
  class-compression data for the latter.

### 6. Fermat / addition-chain inversion
- Standard classical method; discussed in cryptographic resource estimates.
- Prime-field reversible cost is far too large (hundreds of multiplications).
- Not competitive.

### 7. Itoh–Tsujii
- Only for GF(2^n), not GF(p).
- Not applicable to secp256k1.

## Stronger result: coefficient-side compression matches (u, v) compression

A remaining risk in the hybrid Kaliski-jump idea was that even if the `(u, v)`
window transition family compressed well, the coefficient-side `(r, s)`
transforms might explode and ruin the QROM story.

I derived the per-case coefficient matrices directly from the implemented
`kaliski_iteration` logic:

- UEven: `(r, s) -> (r, 2s)`
- VEven: `(r, s) -> (2r, s)`
- UGtV : `(r, s) -> (r+s, 2s)`
- VGtU : `(r, s) -> (2r, r+s)`

Then I ran the same exact 10,000-input window survey for those coefficient-side
matrices.

**Result:** the `(r, s)` side compresses **identically** to the `(u, v)` side.

| w | t | distinct `(u,v)` mats | distinct `(r,s)` mats | max `|entry|` | mean mats/class |
|---|---:|---:|---:|---:|---:|
| 6 | 4 | 125 | 125 | 16 | 4.506 |
| 8 | 4 | 125 | 125 | 16 | 4.493 |
| 8 | 6 | 1133 | 1133 | 64 | 9.461 |

This removed the biggest remaining objection to the hybrid Kaliski-jump
moonshot.

## Strongest result so far: the **joint** transition family also stays tiny

I pushed the classical analysis one step further and measured the *joint* local
transition object that a reversible batched primitive would actually need to
know: the pair `(uv_mat, rs_mat)`, not just each side separately.

Result on the same 10,000 secp256k1 trajectories:

| w | t | distinct `(u,v)` mats | distinct `(r,s)` mats | distinct joint pairs |
|---|---:|---:|---:|---:|
| 6 | 4 | 125 | 125 | **125** |
| 8 | 4 | 125 | 125 | **125** |
| 8 | 6 | 1133 | 1133 | **1133** |

This is the strongest empirical result in the project so far.

Interpretation:
- The coefficient-side transform is not merely similarly compressible — in the
  sampled data it is effectively **functionally locked** to the `(u, v)` side.
- So a hybrid batched primitive may need only **one compressed lookup** for the
  whole local Kaliski window.

## Strongest result so far, refined again: modest side information collapses ambiguity

The remaining practical question is whether the raw key `(u mod 2^w, v mod 2^w)`
is already enough to select the local transition class, or whether we need extra
metadata (which would cost qubits / logic in the eventual quantum version).

I added `src/point_add/kaliski_jump_extra.rs` and measured how much the branch-
sequence ambiguity drops as we augment the key.

For `w = 8`, `t = 4` on 10,000 real secp256k1 trajectories:

| key | mean sequences/class | max sequences/class | singleton classes |
|---|---:|---:|---:|
| `low = (u_low, v_low)` | 4.492 | 16 | 4,102 |
| `low + cmp0` | 2.570 | 8 | 28,731 |
| `low + cmp0 + cmp1` | 1.742 | 4 | 78,817 |
| `low + cmp0 + cmp1 + low1` | 1.696 | 4 | 163,675 |

Interpretation:
- Just adding the **initial compare bit** nearly halves the ambiguity.
- Adding the **compare bit after the first micro-step** cuts the average class
  ambiguity to ~1.74 and the maximum to 4.
- Even the strongest tested key only gets down to ~1.70 average, so there is
  still some residual ambiguity. But it is *tiny*.

This is a huge deal:
- it suggests a practical hybrid batched primitive does **not** need a full
  branch history or a massive QROM key,
- and that a small amount of dynamically-computed side information may be enough
  to select from a very small family of local transition classes.

## New result: brute-force feature search still ranks `(cmp0, cmp1, cmp2)` best

I added `src/point_add/kaliski_key_search.rs` and brute-forced feature subsets
of size up to 4 over a reasonable feature family built from:
- compare bits `cmp0, cmp1, cmp2`,
- a few low bits of `(u1, v1)` and `(u2, v2)`.

On a 300-trajectory heuristic sample (~108k 4-step windows), the best key is still:

> **`(u_low, v_low, cmp0, cmp1, cmp2)`**

with max ambiguity 2.

Important correction:
- that 300-input search is useful for **ranking candidate features**,
- but the absolute mean ambiguity from that small sample was overly optimistic
  as a global estimate.

The real breakthrough came from looking at the *shape* of the residual
ambiguity on the full 10,000-trajectory dataset.

## New strongest result: the 125 four-step classes split into 108 bulk classes + 17 tiny tail classes

I added `src/point_add/kaliski_window_decomp.rs` and decomposed the actual
`w = 8`, `t = 4` class family over 10,000 real secp256k1 trajectories.

### Window-length census
Out of 3,619,614 observed 4-step windows:

| window type | count | fraction |
|---|---:|---:|
| full 4-step windows | 3,589,614 | 99.1712% |
| short windows (last 3 cleanup windows / trajectory) | 30,000 | 0.8288% |

Exactly one length-1, one length-2, and one length-3 tail window appears per
trajectory, i.e. the short windows are a tiny deterministic end effect.

### Distinct class counts
The observed 125 classes decompose as:

| class family | count |
|---|---:|
| distinct length-4 sequences | 108 |
| distinct length-3 sequences | 12 |
| distinct length-2 sequences | 4 |
| distinct length-1 sequences | 1 |
| total | **125** |

So the scary-looking 125-class family is really:
- **108 bulk 4-step classes**, plus
- only **17 tail classes** near termination.

That is much better news for a real reversible design.

## New strongest result: three compare bits determine the **3-step bulk core exactly**

Using the same key

> **`(u_low, v_low, cmp0, cmp1, cmp2)`**

I measured the residual ambiguity structure of the 4-step family.

### Full 10,000-trajectory results for the 4-step family
For key `(u_low, v_low, cmp0, cmp1, cmp2)`:

| metric | value |
|---|---:|
| key classes observed | 261,870 |
| mean sequences / key | 1.275 |
| max sequences / key | 2 |
| ambiguous key classes | 71,936 |
| ambiguous windows | 1,120,661 |
| ambiguous window fraction | 30.9608% |

At first sight, that looks less miraculous than the early 300-input search.
But the **structure** of the ambiguity is the real moonshot result.

### Residual ambiguity structure
Among those 71,936 ambiguous key classes:
- **71,920** are exactly a pair of full 4-step sequences that share the same
  first 3 steps and differ only in the final odd/odd direction:
  - `...-UG` vs `...-VG`
- only **16** ambiguous key classes involve tail windows at all.
- only **4** ambiguous key classes fail to share a common 3-step prefix, and
  those are tiny end-of-algorithm tail effects.

Representative high-frequency residual pairs are exactly of the form:
- `VE-UG-UG-UG  <->  VE-UG-UG-VG`
- `UE-VG-VG-UG  <->  UE-VG-VG-VG`
- `UG-UE-UG-UG  <->  UG-UE-UG-VG`
- `VG-VE-VG-UG  <->  VG-VE-VG-VG`

So for almost the entire inversion trajectory, the key `(low, cmp0, cmp1, cmp2)`
does **not** leave a complicated residual search; it leaves only the **last
odd/odd branch bit** unresolved.

### Exact 3-step bulk core
The most actionable number is this:

| object | distinct classes |
|---|---:|
| all observed 3-step prefixes | 41 |
| full-window 3-step prefixes only | **36** |

And on **full 4-step windows**, the key `(u_low, v_low, cmp0, cmp1, cmp2)`
determines the 3-step prefix **exactly**.

This is the clearest prototype path found so far.

## Updated current best moonshot conclusion

**Conclusion: `hybrid Kaliski-jump` is still the bet, but the best first
prototype target is no longer “exact 4-step lookup.” It is an _exact 3-step
bulk core_ plus a cheap residual step / tail fallback.**

### Why full B-Y replacement is still not the best bet
Full BY jumpdivsteps2 still has two major problems:
1. matrix entries hit the full `2^w` growth;
2. coefficient tracking and cleanup are all-new machinery.

So a *full* B-Y replacement remains very high-risk.

### Why the new decomposition result matters
The old story was:
- 125 four-step classes,
- three compare bits almost collapse them.

The new, much stronger story is:
- the 125 classes split into **108 bulk + 17 tail**,
- 99.17% of actual windows are bulk 4-step windows,
- and on those bulk windows, `(u_low, v_low, cmp0, cmp1, cmp2)` identifies the
  **first 3-step transform exactly**.

That means a practical hybrid primitive can plausibly look like:
1. compute `(u_low, v_low, cmp0, cmp1, cmp2)`,
2. lookup one of only **36** exact bulk 3-step transforms,
3. apply that exact 3-step batched transform to both `(u, v)` and `(r, s)`,
4. do one ordinary final Kaliski micro-step,
5. use a tiny separate fallback for the last 3 windows near `v = 0`.

That is a much cleaner reversible interface than a monolithic exact 4-step
selector.

## New classical proposal: exact 3-step hybrid Kaliski core

### Model
Standard Kaliski / binary almost-inverse update on `(u, v)` has four branch
cases:

```text
if u even:                   (u, v) ← (u/2, v)
elif v even:                 (u, v) ← (u, v/2)
elif u > v:                  (u, v) ← ((u-v)/2, v)
else:                        (u, v) ← (u, (v-u)/2)
```

Each step is a linear map with a shared `1/2` factor. Over `t` steps we get
an integer 2×2 matrix `P_t` with

```text
(u_t, v_t)^T = (1 / 2^t) · P_t · (u_0, v_0)^T.
```

### Best current empirical lead
For `w = 8`:
- the full `t = 4` family has only **125** joint classes,
- but that family decomposes into **108 bulk classes + 17 tail classes**,
- the full-window 3-step prefix family has only **36** classes,
- and `(u_low, v_low, cmp0, cmp1, cmp2)` selects that 3-step prefix exactly on
  99.17% of actual windows.

This now looks like the most actionable structural lead toward reducing the 81%
inversion budget.

## New strongest result: the exact bulk key shrank to 9 bits, but the branch-free 385 CCX core does **not** survive stress-testing

I added `src/point_add/kaliski_prefix_key.rs` and extended
`src/point_add/kaliski_hybrid_proto.rs` to test two hidden assumptions behind
the optimistic exact-3 bulk-core story.

## Exact bulk key width: now much smaller than expected
For **full 4-step windows only**, the exact 3-step bulk prefix is determined by

> **`(u mod 8, v mod 8, cmp0, cmp1, cmp2)`**

That is only a **9-bit key**, not the earlier 19-bit `(u mod 256, v mod 256,
cmp0, cmp1, cmp2)` key.

Measured results:

### secp256k1 trajectories (10,000 inputs)
| width `w` | keys seen | mean prefixes/key | max | exact? |
|---|---:|---:|---:|---:|
| 1 | 20 | 3.200 | 4 | no |
| 2 | 64 | 2.000 | 2 | no |
| 3 | **256** | **1.000** | **1** | **yes** |

### generic exhaustive full-window survey (all `1 ≤ u,v ≤ 512`)
| width `w` | keys seen | mean prefixes/key | max | exact? |
|---|---:|---:|---:|---:|
| 1 | 25 | 3.120 | 4 | no |
| 2 | 78 | 1.923 | 2 | no |
| 3 | **312** | **1.000** | **1** | **yes** |

So the *selector key size* story got materially better.

## Assumption A (derivability) — refuted
Question tested:

> Are `cmp1` and `cmp2` uniquely determined by `(u mod 8, v mod 8, cmp0)`
> alone?

If yes, the branch-free weighted **385.3 CCX** exact-3 core would be a real
forward primitive.

It is **not**.

### secp256k1 result
Using full windows on 10,000 trajectories, for base key
`(u mod 8, v mod 8, cmp0)`:

| metric | value |
|---|---:|
| base classes | 96 |
| mean `(cmp1,cmp2)` pairs per base | **2.667** |
| max `(cmp1,cmp2)` pairs per base | **4** |
| ambiguous base classes | **80** |
| singleton base classes | 16 |

### generic exhaustive result
For all full windows with `1 ≤ u,v ≤ 512`:

| metric | value |
|---|---:|
| base classes | 128 |
| mean `(cmp1,cmp2)` pairs per base | **2.438** |
| max `(cmp1,cmp2)` pairs per base | **4** |
| ambiguous base classes | **96** |
| singleton base classes | 32 |

Interpretation:
- `cmp1` and `cmp2` are **not derivable** from the initial 6 low bits plus
  `cmp0`.
- Therefore the weighted **385.3 CCX** branch-free exact-3 core is a **lower
  bound only**, not a directly implementable forward primitive.
- In a real circuit, `cmp1` and `cmp2` must be produced by actually reaching
  the `u1,v1` and `u2,v2` states (or by an equivalent construction of equal
  cost/complexity).

This is the most important correction in the project so far.

## Assumption A fallout in the real builder
I extended `src/point_add/kaliski_hybrid_proto.rs` to profile a more realistic
builder-level staged exact-3 skeleton:

```text
compare -> exact step 0 -> compare -> exact step 1 -> compare -> exact step 2
```

Measured costs in the real builder:

| object | CCX | Clifford-ish |
|---|---:|---:|
| baseline 3 Kaliski iterations | 4,647 | 23,199 |
| weighted exact-3 lower bound | 385.3 | 5,006.7 |
| one 256-bit `gt` comparator | 256 | 2,308 |
| weighted staged exact-3 (`A` false) | **1,153.3** | **11,930.7** |

So the correct forward-bulk number is no longer 385 CCX; it is **1,153.3 CCX**
for the staged exact-3 skeleton, still much better than 3 ordinary Kaliski
iterations, but not magical.

## Assumption B (cleanup) — keeping compare bits live is much cheaper
Question tested:

> After selection, should the original compare bits be kept live across the
> exact-3 bulk core, or recomputed from the output state when cleaning up the
> selector?

I profiled the recompute strategy directly in the builder by:
- running the known 3-step prefix **backward** from the output state,
- re-computing the three `u > v` compares at the recovered states,
- then restoring the output state with the forward exact-3 prefix.

### Measured cleanup costs
| strategy | CCX | Clifford-ish |
|---|---:|---:|
| keep compare bits live | **0 extra** | **0 extra** |
| recompute cleanup overhead | **1,538.6** | **16,937.4** |

Corresponding total staged exact-3 bulk costs:

| strategy | CCX | Clifford-ish |
|---|---:|---:|
| keep-live total | **1,153.3** | **11,930.7** |
| recompute total | **2,691.9** | **28,868.1** |

Interpretation:
- The **keep-live** strategy is decisively cheaper.
- In this framework, the compare bits are classical measurement results, so
  keeping them live costs classical bit storage, not extra Toffoli.
- Recomputing from output state is far too expensive to be the default plan.

## New strongest result: a naive “real” exact-step controlization is terrible, but a specialized nonterminal bulk primitive is genuinely cheaper

I pushed one step beyond cost-modeling and built two more **real builder-level**
prototypes in `src/point_add/kaliski_hybrid_proto.rs`.

### 1. Naive real forward bulk-3 keep-live prototype
This prototype literally controlizes the case-specific exact operations using
runtime case flags and kept-live compare bits.

Measured cost:

| object | CCX | Clifford-ish | peak qubits |
|---|---:|---:|---:|
| naive real forward bulk-3 keep-live prototype | **11,556** | **34,635** | **1,550** |

This is **worse** than 3 ordinary Kaliski iterations (`4,647 CCX`).

Interpretation:
- directly turning the exact branch formulas into controlled sub/add/shift
  networks is a dead end,
- because controlled shifts and controlled arithmetic destroy the lower-bound
  savings.

So the moonshot cannot be realized as a naive “selector + fully-controlled exact
branch bodies” construction.

### 2. Real specialized nonterminal bulk-3 primitive
I then derived a more faithful prototype directly from the existing
`kaliski_iteration` structure:
- assume the early-window bulk precondition (`f = 1`, `v != 0`),
- remove STEP 0 / eqzero handling,
- specialize the flag logic for `f = 1`,
- keep the rest close to the existing reversible structure.

Measured cost:

| object | CCX | Clifford-ish | peak qubits |
|---|---:|---:|---:|
| real specialized bulk-3 primitive | **3,873** | **20,877** | **1,540** |
| baseline 3 Kaliski iterations | **4,647** | **23,199** | **1,543** |

So this real specialized primitive is an actual win:
- **774 CCX saved** over 3 ordinary iterations,
- about **16.7%** lower Toffoli,
- with slightly lower peak qubits as well.

This is the first genuinely cheaper **real builder-level** bulk primitive in the
project.

### Attempted direct main-path integration: not yet correct
I also tried swapping the specialized bulk primitive into the live Kaliski
forward path for the first 3 iterations. That attempt **failed full point-add
correctness**, so it was backed out immediately.

Current diagnosis:
- the forward primitive is promising,
- but its compatibility with the stored `m_hist` history / generic backward path
  is not yet nailed down.

So this remains a real integrated **prototype**, not yet an enabled circuit
replacement.

## Revised state of the moonshot
The good news:
- the exact bulk key shrank all the way to **9 bits**,
- and the staged exact-3 skeleton is still dramatically cheaper than three
  ordinary Kaliski iterations.

The bad news:
- the optimistic **385 CCX** branch-free forward core is not physically real as
  stated, because `cmp1` and `cmp2` are not derivable from the initial 6 low
  bits plus `cmp0`.

So the right concrete target is now:

> a **staged exact-3 bulk primitive** keyed by the 9-bit exact key, with the
> compare bits generated sequentially and then **kept live** across the core for
> selector cleanup.

I am deliberately pausing all new point-add-level projections until the staged
selector / cleanup path is worked out more concretely.

## Proposed next sessions

### P1. Make the real specialized bulk-3 primitive backward-compatible
The highest-value next step is now:
- pin down why the direct forward-path swap broke correctness,
- determine the exact history-state (`m_hist`, flag semantics, etc.) mismatch,
- and either derive a matching backward for the specialized primitive or prove
  how to make its persistent state identical to the generic step.

This is the shortest path from prototype to an actual circuit replacement.

### P2. Enumerate the exact 36 bulk 3-step classes
For the full-window bulk family, produce:
- canonical representative branch sequences,
- the exact `(uv_mat, rs_mat)` pair,
- the low-bit / compare-bit conditions under which each occurs.

This remains the cleanest classical-to-reversible handoff point for the more
selector-driven variants.

### P3. Build a reversible cost model for the staged exact 3-step core
Estimate the real cost of:
- forming `cmp0, cmp1, cmp2` **sequentially**,
- using the exact 9-bit bulk key `(u mod 8, v mod 8, cmp0, cmp1, cmp2)`,
- keeping the compare bits live across the core,
- applying the corresponding staged `(uv, rs)` transforms,
- then doing one ordinary residual Kaliski step.

This should be compared directly against 3 ordinary Kaliski micro-steps.

### P4. Design the tiny tail fallback
Because every trajectory has exactly three short terminal windows, we can likely
handle them with a separate, tiny cleanup path rather than bloating the bulk
primitive.

## Bottom line

The strongest current research judgement is:

> The best moonshot is still **hybrid Kaliski-jump batching**, but there are now
> two concrete implementation tracks:
> 1. a **staged exact 3-step bulk primitive** using the exact 9-bit key
>    `(u mod 8, v mod 8, cmp0, cmp1, cmp2)` with sequential compare generation
>    and keep-live cleanup, and
> 2. a more conservative **specialized nonterminal bulk-3 primitive** derived
>    directly from the existing Kaliski step, which already shows a real
>    16.7% builder-level win over 3 generic iterations but is not yet wired
>    into the live circuit because backward compatibility is unresolved.

That is still novel research, but it is now tied to a very concrete empirical
structure in the 81%-of-budget hot path, rather than just a vague hope that a
4-step lookup will be small enough.
