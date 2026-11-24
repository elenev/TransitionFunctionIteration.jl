# Ingredients of provenance records: the TFI version and solution ids. The package itself does
# not write solutions to disk.

"TFI version string, recorded in provenance"
tfi_version() = string(pkgversion(@__MODULE__))

"""
    solution_id(obj::ModelSolution) -> String

Identifier of a solution's content: a SHA-1 hash of the grid points and of the values of Pfct,
Vfct, Tfct and Ffct (if any) at the grid nodes. Solutions with identical values have identical
ids, including a solution restored from a file of these values. Simulations and IRFs record it
in their provenance.
"""
function solution_id(obj::ModelSolution)
    values = Float64[]
    foreach(g -> _append_numbers!(values, g), unigrids(obj.ssg))
    for fct in (obj.Pfct, obj.Vfct, obj.Tfct, obj.Ffct)
        isnothing(fct) || _append_numbers!(values, coefficients(fct))
    end
    bytes2hex(SHA.sha1(reinterpret(UInt8, values)))
end

# every number in a (nested) container, in iteration order: grid axes, then the nodes' SVectors
# (Pfct, Vfct, Ffct) or SMatrices (Tfct)
_append_numbers!(out, x::Number) = push!(out, x)
_append_numbers!(out, x) = (foreach(y -> _append_numbers!(out, y), x); out)
