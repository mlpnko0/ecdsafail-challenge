# Hard truth on the current phase bug

## What has now been checked and survived
The specialized bulk-prefix replacement matches the generic path on reachable
state probes much more strongly than before.

In particular, targeted tests now pass for:
- isolated specialized step vs generic step,
- deep reachable-state samples,
- classical transition agreement,
- 3-step forward agreement,
- 3-step forward+backward agreement,
- isolated inverse identity (`with_kal_inv_raw(..., body=[])`) for small tested
  prefix lengths,
- and full prefix forward / forward+backward agreement for a grid of prefix
  lengths:
  `k = 4, 8, 16, 32, 40, 64, 96, 128`.

That means the specialized prefix itself, including its `m_hist` bits and local
phase behavior, looks equivalent to the generic prefix on all tested reachable
samples.

## What still fails
The full `main.rs` harness still reports phase-garbage failures for many `k`
values, even though some larger values pass.

So the remaining bug is now very unlikely to be:
- a local specialized-step phase error,
- a local backward mismatch,
- or a simple `m_hist` incompatibility inside the prefix itself.

## Best current hypothesis
The phase defect is now most likely introduced by a **global interaction**
between the specialized prefix replacement and the outer circuit-seeded,
measurement-based scaffold.

The fact that narrow batch probes gave inconsistent localization results is
consistent with that: the bug does not behave like a simple local wrong gate,
but like a delicate full-scaffold phase effect.
