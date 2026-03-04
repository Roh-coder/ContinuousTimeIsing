# Algorithm README (for QMC collaborators)

This note explains the method in this repository for someone comfortable with condensed-matter QMC ideas, but less familiar with this classical cluster formulation.

## What we are simulating

The code samples the equilibrium distribution of the transverse-field Ising chain in a **continuous-time worldline representation**.

- Spatial coupling: `K`
- Transverse field rate: `Γ`
- Imaginary-time extent: `Δt` (called `layers` in the CLI, used as inverse-temperature extent)

You can think of this as the anisotropic classical `(space × imaginary-time)` Ising representation, handled directly in continuous time (no Trotter slicing grid).

## QMC analogy in one table

If your mental model is SSE / worldline QMC:

- **Operator string / kink worldlines** → here: per-site flip times from a Poisson process
- **Imaginary-time segments** → here: rails (piecewise-constant spin segments)
- **Loop/cluster updates** → here: FK graph clusters on overlapping space-time segments
- **Autocorrelation control via nonlocal updates** → here: whole-cluster flips with probability 1/2

So, algorithmically, this behaves like a nonlocal worldline update scheme, but implemented as a classical FK connectivity problem on continuous-time segments.

## Core update cycle

Each Monte Carlo sweep does:

1. **Worldline representation (`makeRails`)**
   - At each site, flip events are Poisson with mean `Γ Δt`.
   - These events define alternating spin segments along imaginary time.

2. **Phantom cuts (`splitRails`)**
   - Additional cuts are sampled (also Poisson) inside segments.
   - These are helper breakpoints used to construct FK bonds; they do not necessarily survive.

3. **FK graph + cluster build (`make_cluster_with_graph`)**
   - Segments at neighboring sites are considered if their time intervals overlap.
   - For aligned spins, activate bond with
     - `p = 1 - exp(-2 K Δ_overlap)`
   - Connected components of active bonds define clusters.

4. **Cluster flip**
   - Flip each cluster with probability `1/2`.
   - This is the nonlocal decorrelating move (analog of loop updates).

5. **Compression (`clearStaleFlips`)**
   - Remove redundant cuts and merge adjacent equal-spin segments.

## Why this works well

- Continuous time avoids finite-Trotter-step discretization artifacts.
- Cluster flips change large correlated regions in one move.
- Near criticality, this strongly suppresses critical slowing down compared to local updates.

Empirically in this repo, integrated autocorrelation times remain around order-unity configs across the scanned regime, including near the transition.

## Measurements and outputs

Main observable used for crossings:

- **Binder cumulant**
  - `U = 1 - <m^4> / (3 <m^2>^2)`

Data products:

- Time series CSVs (`bench/*.csv`) with per-configuration values
- Summary tables with
  - mean `U`
  - variance
  - standard error
  - integrated autocorrelation time `τ_int`
  - configs required for target precision

## How to explain it quickly to a QMC person

“We sample continuous-time worldlines, then do a Swendsen–Wang/FK cluster update on overlapping space-time segments instead of local spin updates. It is the same spirit as loop updates in QMC: nonlocal collective moves on a worldline representation, yielding very short autocorrelation times.”

## Practical files

- Main implementation: `ConTimeIsing.cpp`
- Build config: `CMakeLists.txt`
- Analysis pipeline used in this workspace: `analyze_timeseries.py`
- Extended performance notes: `ALGORITHM_PERFORMANCE.md`
