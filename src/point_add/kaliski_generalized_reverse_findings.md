# Generalized reverse finding

A direct test was run on the best strict current setting:
- `KAL_BULK3_EXPERIMENT=1`
- `KAL_BULK3_ITERS=255`

Two modes were compared:
1. current explicit specialized backward,
2. a generalized reverse mode that replays the inverse of the whole forward
   block shape while skipping `Hmr`/`R` in reverse, analogous to `emit_inverse`.

## Result
- explicit backward: **passes** cleanly
- generalized reverse mode: catastrophic failure
  - classical mismatches: `9024`
  - phase-garbage batches: `141`
  - ancilla-garbage batches: `141`

## Interpretation
A naive generalized reverse is not the right fix.

The specialized path is not just missing a convenient inverse wrapper; the whole
measurement-based uncompute protocol is too intertwined with the forward block
for a simple `emit_inverse`-style reversal to work.

So the remaining path is not “just use the exact reverse of the specialized
forward block”. It still requires a deliberately designed reverse / cleanup that
matches the measurement protocol.
