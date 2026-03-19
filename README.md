# Modeling Self-Organized Criticality — The Bak–Tang–Wiesenfeld Sandpile Model

> Simulation and analysis of the BTW Abelian Sandpile under SAS IML — power-law avalanches, 1/f noise, and fractal geometry

---

## Overview

This project implements and analyzes the **Bak–Tang–Wiesenfeld (BTW) sandpile model** (1987), a cellular automaton that spontaneously evolves toward a critical state — without any external parameter tuning. This phenomenon is called **Self-Organized Criticality (SOC)**.

The core mechanism:
```
Slow grain deposition  →  Local threshold rule (≥ 4 grains → topple)
       →  Avalanche cascade  →  Power-law size & duration distributions
       →  1/f noise in the dissipation time series
```

Two driving schemes are implemented and compared:

| Model | Driving | Grid | Grains |
|---|---|---|---|
| Center-driven (intuitive) | One grain at center per step | 50×50 | 4,000 |
| Center-driven (optimized) | All grains at once, `loc()` vectorization | 201×201 | up to 1,000,000 |
| Corner-driven variant | One grain at each of 4 corners per step | 201×201 | 40,000,000 |
| Random-driven | Uniform random deposition | 50×50 | 200,000 |

---

## Model Rules

On an **L × L** integer lattice, each site `(r, c)` holds a grain count `z[r,c] ≥ 0`.

At each step:
1. **Drive** — add 1 grain at the selected site
2. **Topple** — if `z[r,c] ≥ 4`:
```
z[r,c]   -= 4
z[r±1,c] += 1   (if neighbor exists)
z[r,c±1] += 1   (if neighbor exists)
```
3. **Repeat** until all sites are stable (`max(z) < 4`)
4. **Boundary** — grains falling outside the grid are lost (open boundary = dissipation)

The toppling rule is the **discrete Laplacian**. Because topplings commute regardless of order, the model is called **Abelian** — the final stable configuration is independent of the toppling sequence.

---

## Key Results

### Center-driven model — Fractal geometry

The Abelian structure produces striking **fractal spatial patterns** visible in the heatmaps. Colors represent grain counts mod 4 (0 to 3). As grain count grows, the sandpile expands in a self-similar, rotationally symmetric pattern:

| Configuration | Observation |
|---|---|
| 201×201 — 4,000 grains | Small localized diamond pattern at center |
| 201×201 — 40,000 grains | Expanding fractal oval emerges |
| 201×201 — 400,000 grains | Dense fractal fills most of the grid |
| 151×151 — 1,000,000 grains | Full fractal coverage, fine-grained structure |

Corner-driven variant (201×201, 40,000,000 grains) produces a **4-fold symmetric fractal** with hyperbolic arms extending from each corner.

---

### Random-driven model — Power-law (1/f) statistics

On a 50×50 grid with 200,000 random grain depositions, the system reaches a **stationary critical state**. Avalanche sizes and durations were tracked for every grain addition.

**Avalanche size distribution — log-log regression:**

| Parameter | Value |
|---|---|
| Estimated slope `τ̂` | **−1.048** |
| R² | **0.99** |
| p-value | < 0.0001 |
| Observations | 148 |

> The exponent `τ̂ ≈ 1.04` falls within the canonical BTW range `1 ≲ τ ≲ 1.3` for 2D lattices, confirming SOC under random drive.

**Interpretation:**
- Most avalanches are small and localized
- Large, system-spanning avalanches are rare but follow the **same statistical law** — no characteristic scale
- The log-log linearity confirms `P(s) ∼ s⁻τ` and `P(T) ∼ T⁻α`
- Slight curvature for large `s` is a finite-size effect (grid boundary limits avalanche spread)

---

## Three Implementations

### 1 — Intuitive version (center-driven)
Sequential grain deposition with nested loops. Conceptually clear, computationally slow for large grids.
```sas
/* Core toppling loop */
if Tas[r,c] >= 4 then do;
   Tas[r,c] = Tas[r,c] - 4;
   if r > 1  then Tas[r-1,c] = Tas[r-1,c] + 1;
   if r < 50 then Tas[r+1,c] = Tas[r+1,c] + 1;
   if c > 1  then Tas[r,c-1] = Tas[r,c-1] + 1;
   if c < 50 then Tas[r,c+1] = Tas[r,c+1] + 1;
end;
```

### 2 — Optimized version (center-driven)
All grains deposited at once. Vectorized toppling using `loc()` + auxiliary `Add` matrix. Suitable for large grids (201×201, 400,000+ grains).
```sas
/* Vectorized toppling */
idx = loc(Tas >= 4);
rr  = ceil(idx / n);
cc  = idx - (rr - 1)*n;
/* Simultaneous removal + redistribution via Add matrix */
Tas = Tas + Add;
```

### 3 — Random-driven model
Uniform random deposition across the grid. Tracks avalanche `Size` (number of cells toppled) and `Time` (number of toppling rounds) per grain. Exports log-log data for power-law regression.

---

## Statistical Analysis Pipeline
```
Random deposition (200,000 grains)
    └── Track Size[i] and Time[i] per grain
            └── Filter non-zero events
                    └── call tabulate() → frequency tables
                            └── log-log transformation
                                    └── PROC SGPLOT → visual power-law
                                            └── PROC REG → estimate τ̂
```


