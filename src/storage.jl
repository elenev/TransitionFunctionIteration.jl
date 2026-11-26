# Saving and comparing results with OutputVault. The methods of these functions come with the
# package extension TransitionFunctionIterationOutputVaultExt, loaded with `using OutputVault`;
# `assemble_moments` for results in memory is defined in moments.jl.

"""
    save_result(rf, result; start = :ss, dirty_repo_check = :warn) -> String

Save the output of [`solve_and_simulate`](@ref) (or a NamedTuple with some of its properties) in
the OutputVault run `rf`. With `obj`, the run starts over: `save_solution(rf, (; obj, stv))`, then
the parameters at `params/params` ([`collect_parameters`](@ref)). Then, for each property present:
`simseries` and `errs` at `simulation/simseries` and `simulation/errs`, `errstats` at
`errstats/errstats`, `stats` at `stats/<sample>`, the IRF `irf` at `irf/<start>/stats` and
`irf/<start>/paths`, and the interpolants of `recursions` (e.g. `claims`, `welfare`) at
`recursions/<name>`. Without `obj` the run must exist, and only the other properties are saved.
Returns the run's folder. Needs `using OutputVault`.
"""
function save_result end

"""
    load_result(rf, item::Symbol; model = nothing, native = false)
    load_result(rf; model = nothing, native = false) -> NamedTuple

Load a saved item of the OutputVault run `rf` as TFI's types: `:solution` (`(; obj, stv)`),
`:params` and `:errstats` (`MomentArray`s), `:simulation` (`(; simseries, errs)`), `:stats` (moments by
sample), `:irf` (a NamedTuple of `IRFResult`s by starting point) or `:recursions` (a NamedTuple of
`Interpolant`s). Without `item`, every item the run holds, in one NamedTuple with the properties
of `solve_and_simulate`'s output.

The portable solution files store the model's metadata, not the model: `model` is a model (its
metadata are compared with the stored ones, and differences are a warning) or a model type (the
model is rebuilt with [`model_from_metadata`](@ref)). `native = true` reads the solution's native
copy instead, which needs no model but the code that saved it. Needs `using OutputVault`.
"""
function load_result end

"""
    migrate_irf!(rf)
    migrate_irf!(pattern; move = false)

Move the IRF files of runs saved with OutputVault 0.3 (`irf/irf` and `irf/irf_paths`) to the
current paths `irf/ss/stats` and `irf/ss/paths`, so that queries find them. For a `RunPattern`,
every run found is listed with what would move, and the files move only with `move = true`.
Needs `using OutputVault`.
"""
function migrate_irf! end

"""
    result_shape(kind::Symbol)

The OutputVault `QueryShape` of TFI's saved tables: `:stats` (`"stats/{sample}"`, key `variable`,
column `stat`), `:irf` (`"irf/{start}/stats"`, keys `__shock`, `__stat` and `__period`, column
`variable`) or `:params` (`"params/params"`, key `parameter`, column `stat`). Needs
`using OutputVault`.
"""
function result_shape end
