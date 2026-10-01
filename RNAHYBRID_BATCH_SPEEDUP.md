# Batched RNAhybrid for Cas13design

`scripts/RfxCas13d_GuideScoring.R` scores each candidate guide's hybridization
energy by shelling out to RNAhybrid. The original Install.txt workflow spawns
one RNAhybrid process per guide/window pair, which dominates runtime once a
transcript has more than a few hundred candidate guides.

`scripts/rnahybrid_paired.c` replaces that per-pair loop with one persistent
RNAhybrid 2.1.2 process that reads all target/query pairs for a transcript at
once. It links the original RNAhybrid 2.1.2 dynamic-programming and energy
code unchanged — it only changes how input/output is batched, not the energy
calculation. `scripts/RNAhyb.sh` is the entry point `RfxCas13d_GuideScoring.R`
calls; it execs the compiled `scripts/rnahybrid_paired` binary.

## Building it

```
cd scripts
tar -xzf RNAhybrid-2.1.2.tar.gz     # per Install.txt
./build_rnahybrid_paired.sh
```

This produces `scripts/rnahybrid_paired`, which `RNAhyb.sh` requires.

## Correctness

Verified against vanilla `RNAhybrid -c -s 3utr_human` with zero MFE mismatches
on:
- 5,636 real 9-mer and 12-mer target/query pairs
- six explicit lowercase/U/N/homopolymer/9-mer/12-mer edge cases
- a 982-nt full Cas13design run (960 guide records) — identical SHA-256 output
- an 8,389-nt full Cas13design run (8,367 guide records) — identical SHA-256 output

## Measured performance

Measured 2026-08-05 in the track-cas13 workspace against the human v50
Cas13design workload:

| Workload | Original (per-pair) | Batched | Speedup |
|---|---:|---:|---:|
| 1,460 RNAhybrid pairs | 18.95 s | 0.04 s | 474x |
| Full Cas13design run, 982 nt | 36.35 s | 6.73 s | 5.4x |
| Full Cas13design run, 8,389 nt | 233.97 s | 18.40 s | 12.7x |

For the 8,389-nt run, the hybridization step itself fell from ~215 s to
under 1 s; remaining runtime is RNAfold, model loading/prediction, and CSV
generation, none of which this change touches.

At observed production concurrency (one worker 571.9 bp/s vs. an observed
active-run range of 2.0–43.3 bp/s), this projects to roughly **53x–1,144x**
depending on load. In one concrete example, an average 4.4 Mb chunk fell from
about 28 hours to about 32 minutes.

A 100,000-pair stress test completed in 2.38 s with 45 MB peak RSS.

## RfxCas13d_GuideScoring.R changes that go with this

This copy of `RfxCas13d_GuideScoring.R` also differs from the original
Install.txt version in two smaller ways, both aimed at per-transcript
overhead rather than the hybridization step above:

- RNAplfold is invoked with `-u 23 -c 1` instead of `-u 50`: it now computes
  unpaired-probability windows only up to the actual guide length (23 nt)
  instead of 50, and `-c 1` suppresses the dot-plot postscript files RNAplfold
  otherwise writes per transcript.
- An optional 4th CLI argument lets many transcripts append to one shared
  output CSV instead of each transcript creating and closing its own file —
  avoiding one file-create per transcript on a shared filesystem.

Neither of these has an isolated benchmark the way the RNAhybrid batching
does above; they reduce RNAplfold output volume and filesystem I/O rather
than changing an algorithm, so no multiplier is claimed for them here.
