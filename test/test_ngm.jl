using LinearAlgebra

# Runs the full NGM: main_setup, main_execute!, simulate, describe, calc_error_stats.
# Top-level bindings (ngm, obj, stv, simseries, errs, stats, errstats) land in global
# scope via `include`, making them available to test_linearize.jl.
include("../examples/ngm.jl")

# ── Klein (2000) perturbation solver ──────────────────────────────────────────

"""
Solves the linear expectations model y' = W y + M z + Π η, z' = Φ z + ϵ
using the Klein (2000) ordered real Schur decomposition. Returns (Ψ_x, Ψ_z).
"""
function klein(W::AbstractMatrix, M::AbstractVecOrMat, Φ::Union{Real,AbstractMatrix}, nx::Integer)
    F_schur  = schur(W)
    F_ordered = ordschur(F_schur, abs.(F_schur.values) .< 1.0)

    T = F_ordered.T
    Z = F_ordered.Z

    Z11 = Z[1:nx, 1:nx];        Z12 = Z[1:nx, nx+1:end]
    Z21 = Z[nx+1:end, 1:nx];    Z22 = Z[nx+1:end, nx+1:end]

    Psi_x = Z21 / Z11

    M_rot = Z' * M
    M2    = M_rot[nx+1:end, :]
    T22   = T[nx+1:end, nx+1:end]

    F = if length(T22) == 1 && length(Φ) == 1
        -M2 / (T22[1] - Φ)
    else
        sylvester(-T22, Φ, -M2)
    end

    Psi_z = (Z22 - Psi_x * Z12) * F
    return Psi_x, Psi_z
end

function klein(; XX, XC, CX, CC, XZ, CZ, Φ)
    return klein(vcat(hcat(XX,XC), hcat(CX,CC)), vcat(XZ,CZ), Φ, size(XX,1))
end

# ── Analytic linearization for the NGM ────────────────────────────────────────

function solve_linearized(model::NeoclassicalGrowthModel)
    stv  = steadystate(model)
    Kbar = stv.endog.K;  Cbar = stv.values.C;  Abar = 1.0
    (; α, β, γ, ρ) = model

    η = Abar * Kbar^(α-1)   # steady-state MPK
    κ = α * β * η / γ

    KK = 1/β;          KC = -Cbar/Kbar;              KA = η
    CK = κ*(α-1)/β;    CC = 1 - κ*(α-1)*Cbar/Kbar;   CA = κ*((α-1)*η + ρ)

    Ck, Ca = klein(XX=KK, XC=KC, CX=CK, CC=CC, XZ=KA, CZ=CA, Φ=ρ)
    Ck = length(Ck) == 1 ? Ck[1] : Ck
    Ca = length(Ca) == 1 ? Ca[1] : Ca

    Kk = KK + KC * Ck
    Ka = KA + KC * Ca

    vcoefs = TransitionFunctionIteration.MomentArray([0.0 Ca Ck],      (:C,),  (:constant, :A, :K))
    tcoefs = TransitionFunctionIteration.MomentArray([0.0 Ka Kk 0.0],  (:K′,), (:constant, :A, :K, :A′))
    return vcoefs, tcoefs
end

# ── Tests ─────────────────────────────────────────────────────────────────────

@testset "NGM — perturbation vs simulation linearization" begin
    Vlin, Tlin = solve_linearized(ngm)
    Vsim, Tsim = linearize(simseries, stv)

    @test Vlin isa TFI.MomentArray
    @test Tlin isa TFI.MomentArray
    @test Vsim isa TFI.MomentArray
    @test Tsim isa TFI.MomentArray

    # Exclude the constant: stochastic SS differs slightly from deterministic SS,
    # so we don't expect the intercepts to match.
    @test isapprox(Vlin[:, 2:end], Vsim[:, 2:end], atol=5e-4)
    @test isapprox(Tlin[:, 2:end], Tsim[:, 2:end], atol=5e-4)
end
