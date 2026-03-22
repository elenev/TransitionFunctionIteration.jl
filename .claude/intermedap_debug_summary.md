# intermedap.jl Debugging Session Summary

**Date:** 2026-03-22
**Branch:** rev202603
**Files:** `examples/intermedap.jl`, `examples/intermedap.typ`

---

## Problem

`calc_at_guess(mdl, xopt)` should return all-zero residuals at the steady state, but several residuals were non-zero.

---

## Bugs Found and Status

### Bug 1 — Household bond demand BH (FIXED by user)

**Location:** `calc_eqm`, line computing `BH`

**Error:** The stationarized bond repayment from last period is `(1/G)*b̂_old`, but the code used `B` (= `p*b̂_old`, the bond *value*) directly as if it were the repayment, giving a factor-of-G error.

**Correct formula:**
```julia
BH = (B/(G*p) + D + w - C) / p
```

**Effect:** MC_B was -0.22. After fix, reduces to -0.006.

---

### Bug 2 — λ̃ formula in steadystate_guess (FIXED by user in code; needs .typ fix)

**Location:** `steadystate_guess`, line 67; and `intermedap.typ` Step 2 (lines 267–272)

**Error:** The normalized leverage multiplier was derived as `λ̃ = 1 - v^F/μ̂`. The correct steady-state bond Euler gives:
```
p = λ̃ + p*(v^F/μ̂)   ⟹   λ̃ = p * (1 - v^F/μ̂)
```
The `.typ` file's Step 2 implicitly assumes `M̂_ss * exp(-g) = 1` (valid only if `β = G^γ`).

**Fix needed in `.typ`:** Replace Step 2 text and formula with:
```
E[M̂ v^F exp(-g)] / μ̂ = p · (v^F/μ̂)   at steady state
⟹ λ̃(x) = p · (1 - (v^F/μ̂)(x))
```

**Effect:** EE_BI residual. Also changes RK, k, K, N — x* must be re-solved (user did this with Roots.jl, getting xopt ≈ 0.005576).

---

### Bug 3 — BI and LC use un-stationarized K' (PENDING)

**Location:** `calc_eqm`, BI formula and LC formula

**Error:** State `K` = K̂_{t-1} (in Z_{t-1} units). The capital choice `K' = (1-δ+i)*K` is still in Z_{t-1} units. At steady state `K' = G*K`. The stationarized budget and leverage constraint need `K̂_t = K'/G` (in Z_t units), which equals `K` at steady state. Using `K'` instead of `K'/G` introduces a factor of G.

**Magnitude:** `MC_B = -K*(G-1)/p ≈ -0.02*K/p ≈ -0.006` ✓

**Fix:**
```julia
# BI formula
BI = (q*K′/G - ((1-η) + x - ξ/2 * x^2) * N) / p

# Leverage constraint
LC = fisher_burmeister(λ, θ*q*K′/G - p*BI)

# Net worth transition (IterateMode)
N′ = @. ((K′/G) * np.continuations.RK - BI) / np.continuations.G
```

**Effect:** Fixes MC_B = -0.006 and LC = -0.015.

---

### Bug 4 — Wage formula in .typ (correctness issue, does not affect residuals)

**Location:** `intermedap.typ` line 329

**Error:** The .typ writes `ŵ = (1-α)*K̂^α`. From the labor FOC `(1-α)*l^{-α} = ŵ` with labor market clearing `l = G/K̂` at steady state:
```
ŵ = (1-α)*(G/K̂)^{-α} = (1-α)*(K̂/G)^α
```
The G^{-α} factor is missing.

**Fix in .typ line 329:** `hat(w) = (1-alpha)(hat(K)/G)^alpha`

Similarly, stationarized output is `Ŷ = K̂^α * G^{1-α}` (not `K̂^α`).

**Effect:** The w error cancels in equilibrium residuals (same formula in both `steadystate_guess` and `calc_eqm`), so fixing this does not change any residuals. It does fix the economic interpretation of Y, w, C, and the household labor wedge.

---

## What Remains To Do

1. Fix `.typ` Step 2 (λ̃ derivation) — text edit only
2. Fix `.typ` line 329 (ŵ and Ŷ) — text edit only
3. Fix `calc_eqm` BI formula: `K′/G` instead of `K′`
4. Fix `calc_eqm` LC formula: `K′/G` instead of `K′`
5. Fix `calc_eqm` IterateMode N′ formula: use `K′/G` as capital base
6. After fixes: run `calc_at_guess(mdl, xopt)` to verify all residuals ≈ 0

---

## Key Mathematical Relationships

```
State:  K = K̂_{t-1}  (stationarized by Z_{t-1})
        L = N̂_t / (q_t * K̂_{t-1})   [balance sheet fraction]
        N = L*q*K = N̂_t              (in Z_t units)

Capital accumulation:
        K′ = (1-δ+i)*K               (still in Z_{t-1} units)
        K̂_t = K′/G                  (stationarized by Z_t)

At steady state: K′ = G*K,  K′/G = K

Bond market at steady state:
        b̂ = θ*K/p                   (leverage constraint binds)
        B = q*K - N = p*b̂           (balance sheet identity)

Correct λ̃:  λ̃ = p*(1 - v^F/μ̂)    (from bond Euler, p = β*exp(-γ*g))
Correct ŵ:  ŵ = (1-α)*(K̂/G)^α    (from labor FOC with l = G/K̂)
```
