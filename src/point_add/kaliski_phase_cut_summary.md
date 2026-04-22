# Phase cut summary for strict failing k=4 batch

Using a main-like replay of the circuit-seeded sampling/randomness stream for
`KAL_BULK3_ITERS = 4`, the first strict phase-failing batch occurs at batch
index **10**.

For that exact batch, the experimental circuit's phase mask is already:

- `0x0000040000000000`

at every late top-level cut that was probed:
- after pair1
- after `tx <- Rx - Qx`
- after `mul3_between_pair`
- before `lam` free

## Interpretation
This means the strict phase bug is already present by the end of pair1 on the
actual failing batch. The later scaffold operations tested so far preserve that
same phase mask rather than introducing it.

So the bug hunt should move earlier, not later:
- deeper inside pair1's inverse-carrying path,
- or at the pair1/body boundary,
- rather than around `lam` free or the pair2 composition.
