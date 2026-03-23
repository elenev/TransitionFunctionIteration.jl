# Convenience constructors for exogenous environments

import QuantEcon.rouwenhorst

# Discretize AR(1) using Rouwenhorst method
function ExogenvAR1(;
    ρ::Real = 0.0,
    σ::Real = 1.0,
    μ::Real = 0.0,
    n::Int = 3,
    jensen::Bool = false,
    shock_volatility::Bool = false,
    transform::Function = jensen ? exp : identity,
    shock_name::Symbol = :z)

    if ~shock_volatility
        # If the supplied volatility is the volatility of the state rather than the shock,
        # we need to adjust the parameters
        σ = σ * sqrt(1 - ρ^2)
    end

    if jensen
        # Jensen correction to the mean
        # such that expectation is exp(μ)
        # rather than exp(μ + σ^2/2)
        μ -= σ^2/2
    end


    # Markov chain constructor from QuantEcon
    mc = rouwenhorst(n, ρ, σ, μ*(1-ρ))

    state_values = mc.state_values
    mtrans = mc.p
    if n <= 10
        state_values = SVector{n}(state_values)
        mtrans = SMatrix{n, n}(mtrans)
    end

    pts = StructArray((; shock_name => transform.(state_values)))

    # Steady state index is middle point
    ssidx = (n + 1) ÷ 2

    #@info "Created Exogenv with ρ=$(ρ), σ=$(σ), μ=$(μ), n=$(n), jensen=$(jensen)"
    #return mc, pts, ssidx
    return Exogenv(; pts, mtrans, ssidx)

end

# Blend multiple independent Markov chains into a product-space chain
function Exogenv(chains::Exogenv...)
    # Validate: all field names must be distinct across chains
    all_names = reduce(vcat, [collect(propertynames(c.pts)) for c in chains])
    if length(all_names) != length(unique(all_names))
        error("Field names must be distinct across blended chains; got: $all_names")
    end

    # Transition matrix: kron product, reversed so ordering matches Iterators.product
    # (first index varies fastest in both kron(P_last,...,P_1) and product(1:N1,...,1:N_last))
    P_blended = reduce(kron, reverse([transition_matrix(c) for c in chains]))

    # Enumerate all state combinations — first chain index varies fastest
    ranges = tuple((eachindex(c.pts) for c in chains)...)
    combos = vec(collect(Iterators.product(ranges...)))

    # Build StructArray by merging each chain's NamedTuple at the relevant state index
    pts_blended = StructArray([
        reduce(merge, (chains[d].pts[combo[d]] for d in eachindex(chains)))
        for combo in combos
    ])

    # Steady-state index: linear position of the tuple of per-chain steady states
    ss_combo = tuple((ssidx(c) for c in chains)...)
    ssidx_blended = findfirst(==(ss_combo), combos)

    return Exogenv(pts = pts_blended, mtrans = P_blended, ssidx = ssidx_blended)
end

# --- Ergodic statistics ---

"Compute the ergodic (stationary) distribution of the Markov chain."
function ergodic_distro(exogenv::AbstractExogenv)
    # Stationary distribution = right eigenvector of P' at eigenvalue 1
    F = eigen(Matrix(transition_matrix(exogenv))')
    idx = argmin(i -> abs(real(F.values[i]) - 1), eachindex(F.values)) #abs.(real.(F.values) .- 1))
    pvec = real(F.vectors[:, idx])
    return pvec ./ sum(pvec)
end

"Ergodic mean of each state variable. Returns a NamedTuple."
function Statistics.mean(exogenv::AbstractExogenv)
    pvec = ergodic_distro(exogenv)
    return map(col -> sum(p * x for (p, x) in zip(pvec, col)),
               StructArrays.components(points(exogenv)))
end

"Ergodic variance of each state variable. Returns a NamedTuple."
function Statistics.var(exogenv::AbstractExogenv)
    pvec = ergodic_distro(exogenv)
    μ = mean(exogenv)
    return map((col, m) -> sum(p * x^2 for (p, x) in zip(pvec, col)) - m^2,
               StructArrays.components(points(exogenv)), μ)
end

"Ergodic standard deviation of each state variable. Returns a NamedTuple."
Statistics.std(exogenv::AbstractExogenv) = map(sqrt, var(exogenv))

"State with highest ergodic probability. Returns `(idx, state_NamedTuple)`."
function mode(exogenv::AbstractExogenv)
    pvec = ergodic_distro(exogenv)
    idx = argmax(pvec)
    return idx, points(exogenv)[idx]
end

"Middle state by index. Returns `(idx, state_NamedTuple)`, or `nothing` if Nstates is even."
function Statistics.median(exogenv::AbstractExogenv)
    N = exnpt(exogenv)
    if iseven(N)
        @warn "Even number of states ($N). Cannot compute single median point."
        return nothing, nothing
    end
    idx = (N + 1) ÷ 2
    return idx, points(exogenv)[idx]
end
