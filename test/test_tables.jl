# Tables.jl interfaces, DataAPI.jl metadata and table constructors (MomentArray, simulation and
# error panels, IRF results), including round trips through a table that keeps metadata. Real
# files (parquet) are tested in OutputVault. Uses obj, stv, simseries, errs from test_ngm.jl.
using Tables, DataAPI

# A table with table-level DataAPI metadata, standing in for a file format that keeps metadata
# (e.g. parquet), so that the constructors can be tested without one.
struct MetaTable{C<:NamedTuple}
    columns::C
    metadata::Dict{String,String}
end
Tables.istable(::Type{<:MetaTable}) = true
Tables.columnaccess(::Type{<:MetaTable}) = true
Tables.columns(t::MetaTable) = t.columns
DataAPI.metadatasupport(::Type{<:MetaTable}) = (read = true, write = false)
DataAPI.metadatakeys(t::MetaTable) = keys(t.metadata)
DataAPI.metadata(t::MetaTable, key::AbstractString; style::Bool = false) =
    style ? (t.metadata[key], :default) : t.metadata[key]
DataAPI.metadata(t::MetaTable, key::AbstractString, default; style::Bool = false) =
    (v = get(t.metadata, key, default); style ? (v, :default) : v)

# what a writer and a reader do together: copy the columns into plain vectors, keep the metadata
# as strings
roundtrip(x) = MetaTable(map(collect, Tables.columntable(x)),
                         Dict(String(k) => string(DataAPI.metadata(x, k)) for k in DataAPI.metadatakeys(x)))

@testset "Tables — MomentArray" begin
    ma = MomentArray([1.0 2.0; 3.0 4.0; 5.0 6.0], (:a, :b, :c), (:mean, :sd))
    @test Tables.istable(ma)
    cols = Tables.columns(ma)
    @test keys(cols) == (:__row, :mean, :sd)
    @test cols.__row == ["a", "b", "c"] && cols.sd == [2.0, 4.0, 6.0]
    Tables.columns(ma)
    @test (@allocated Tables.columns(ma)) == 0
    @test DataAPI.metadata(ma, "tfi.rowkey") == "__row"
    @test DataAPI.metadata(ma)["tfi.kind"] == "MomentArray"

    # restore without names (run-time type), with names (checked), with names and eltype (inferred)
    nt = Tables.columntable(ma)
    @test MomentArray(nt) == ma && row_names(MomentArray(nt)) == (:a, :b, :c)
    @test col_names(MomentArray(nt)) == (:mean, :sd)
    @test MomentArray{(:a, :b, :c), (:mean, :sd)}(nt) == ma
    @test_throws ArgumentError MomentArray{(:a, :b), (:mean, :sd)}(nt)
    restore(t) = MomentArray{(:a, :b, :c), (:mean, :sd), Float64}(t)
    @test (@inferred restore(nt)) == ma
    @test MomentArray(nt, (:a, :b, :c), (:mean, :sd)) == ma
    # a renamed row key, recorded nowhere: supplied by keyword
    nt2 = (variable = ["a", "b", "c"], mean = [1.0, 3.0, 5.0], sd = [2.0, 4.0, 6.0])
    @test MomentArray(nt2; rowkey = :variable) == ma

    # the row key and kind travel in the metadata
    ds = roundtrip(ma)
    @test DataAPI.metadata(ds, "tfi.kind") == "MomentArray"
    @test MomentArray(ds) == ma
    @test MomentArray{(:a, :b, :c), (:mean, :sd)}(ds) == ma
    @test_throws ArgumentError SimulationPanel(ds)            # wrong kind
end

@testset "Tables — simulation and error panels" begin
    T, R = size(simseries)
    prov = simseries.provenance
    @test prov isa Provenance{SimulationSettings}
    @test prov.solution == solution_id(obj) && prov.tfi_version == TFI.tfi_version()
    @test prov.settings.seed == 100 && prov.settings.dynamic_align == :drop && prov.settings.exog_start isa String
    # one provenance object, shared by the error panel and by slices
    @test errs.provenance === prov
    @test simseries[:, 1:2].provenance === prov && view(errs, :, 1).provenance === prov

    cols = Tables.columns(simseries)
    @test keys(cols)[1:2] == (:__run, :__period)
    @test length(cols.__run) == T * R
    @test cols.__period[1:3] == [1, 2, 3] && cols.__period[T + 1] == 1
    @test cols.__run[T] == 1 && cols.__run[T + 1] == 2
    @test cols.K == vec(simseries.K)
    @test eltype(cols.__run) == Int64
    @test DataAPI.metadata(simseries, "tfi.endog") == "[\"K\"]"
    # provenance is one JSON entry, parsed back from the struct types
    @test TFI.JSON.parse(DataAPI.metadata(simseries, "tfi.provenance"), Provenance{SimulationSettings}) == prov

    # round trip: data, endogenous states, provenance
    sp = SimulationPanel(roundtrip(simseries))
    @test sp.data == simseries.data
    @test endog_states(sp) == endog_states(simseries)
    @test sp.provenance == prov
    ep = SimulationErrorPanel(roundtrip(errs))
    @test ep.data == errs.data && ep.provenance == prov
    @test_throws ArgumentError SimulationPanel(roundtrip(errs))   # wrong kind

    # sources without metadata need the endogenous names and have no provenance; row order does not matter
    nt = Tables.columntable(simseries)
    @test_throws ArgumentError SimulationPanel(nt)
    @test SimulationPanel(nt, (:K,)).data == simseries.data
    @test isnothing(SimulationPanel(nt, (:K,)).provenance)
    rev = map(reverse, nt)
    @test SimulationPanel(rev, (:K,)).data == simseries.data
    @test SimulationErrorPanel(Tables.columntable(errs)).data == errs.data
end

@testset "simulate — supplied exogenous path" begin
    exog = fill(TFI.ssidx(obj.exogenv), 30, 2)
    s2, e2 = simulate(obj, exog; endog_start = stv.endog)
    @test size(s2) == (29, 2) == size(e2)           # 30 rows: 29 output periods and the extra one
    @test isnothing(s2.provenance.settings.seed) && isnothing(s2.provenance.settings.exog_start)
    @test e2.provenance === s2.provenance
    @test SimulationPanel(roundtrip(s2)).provenance == s2.provenance
    # the old keywords that conflicted with seed and exog_start are gone
    @test_throws MethodError simulate(obj, stv; exog_states = exog)
    @test_throws MethodError simulate(obj, stv; shocks = rand(30, 2))
end

@testset "Tables — IRF results" begin
    irf = simulate_irf(obj, simseries; shocks = (pos = 5, neg = 1), path_vars = (:K, :C),
                       other_stats = (; p90 = x -> quantile(x, 0.9)), T = 8, N = 40)
    cols = Tables.columns(irf)
    @test keys(cols)[1:3] == (:__shock, :__stat, :__period)
    @test Set(cols.__stat) == Set(["mean", "std", "p90"])

    function same(a, b)
        a.means == b.means && a.stdevs == b.stdevs && a.other_stats == b.other_stats &&
            all(k -> a.paths[k] == b.paths[k], keys(a.paths))
    end
    @test same(IRFResult(Tables.columntable(irf), Tables.columntable(irf_paths(irf))), irf)

    # provenance: settings the statistics do not record
    prov = irf.provenance
    @test prov isa Provenance{IRFSettings} && prov.solution == solution_id(obj)
    @test prov.settings.seed == 100 && prov.settings.N == 40
    @test prov.settings.shocks == repr((pos = 5, neg = 1))

    dir = mktempdir()
    ds_stats = roundtrip(irf)
    ds_paths = roundtrip(irf_paths(irf))
    @test DataAPI.metadata(ds_paths, "tfi.kind") == "IRFPaths"
    @test DataAPI.metadata(ds_paths, "tfi.provenance") == DataAPI.metadata(ds_stats, "tfi.provenance")
    @test same(IRFResult(ds_stats, ds_paths), irf)
    @test IRFResult(ds_stats, ds_paths).provenance == prov
    @test isnothing(IRFResult(Tables.columntable(irf)).provenance)       # no metadata, no provenance

    # no saved paths: empty paths table, statistics only
    irf0 = simulate_irf(obj, simseries; shocks = (pos = 5,), T = 8, N = 40)
    p0 = Tables.columntable(irf_paths(irf0))
    @test isempty(p0.__shock)
    r0 = IRFResult(roundtrip(irf0), p0)
    @test r0.means == irf0.means && all(isempty, r0.paths)
    @test IRFResult(roundtrip(irf0)).stdevs == irf0.stdevs
    @test r0.provenance == irf0.provenance
end

# Allocation checks (allocs_after_warmup from test_performance.jl): wrappers allocate no more
# than the arrays they wrap, writing builds columns without copying data, and reading copies
# nothing beyond what the source's own reads allocate.
read_lazy_indices(cols) = sum(cols.__run) + sum(cols.__period)
# restores with the target type fixed at compile time, as in user code
const ALLOC_ROWS = ntuple(i -> Symbol("r$i"), 12)
const ALLOC_COLS = ntuple(j -> Symbol("c$j"), 6)
restore_static(table) = MomentArray{ALLOC_ROWS, ALLOC_COLS, Float64, SMatrix{12, 6, Float64}}(table)
restore_dense(table) = MomentArray{ALLOC_ROWS, ALLOC_COLS, Float64}(table)

@testset "Allocations — MomentArray and panel operations" begin
    ma = MomentArray([1.0 2.0; 3.0 4.0; 5.0 6.0], (:a, :b, :c), (:mean, :sd))
    ms = MomentArray(SMatrix{3,2}(ma.data), (:a, :b, :c), (:mean, :sd))
    for m in (ma, ms)
        @test allocs_after_warmup(m -> TFI.getmoment(m, :a, :mean), m) == 0
        @test allocs_after_warmup(m -> m(:b, :sd), m) == 0
        @test allocs_after_warmup(m -> TFI.getmoment(m, (:a, :c), :mean), m) == 0
        @test allocs_after_warmup(m -> view(m, 1:2, :), m) <= allocs_after_warmup(d -> view(d, 1:2, :), m.data)
        @test allocs_after_warmup(sum, m) <= allocs_after_warmup(sum, m.data)
        # broadcasting runs on the wrapped matrix: same result type, no extra allocation
        @test allocs_after_warmup(m -> m .* 2.0, m) <= allocs_after_warmup(d -> d .* 2.0, m.data)
        @test allocs_after_warmup(m -> m .* 2.0 .+ m, m) <= allocs_after_warmup(d -> d .* 2.0 .+ d, m.data)
    end

    errfield = first(propertynames(errs))
    for (panel, field) in ((simseries, :K), (errs, errfield))
        for f in (p -> p[3, 1], p -> view(p, :, 1), p -> p[1:5, :], p -> getproperty(p, field), size)
            @test allocs_after_warmup(f, panel) <= allocs_after_warmup(f, panel.data)
        end
    end
    @test allocs_after_warmup(endog_states, simseries) == 0
end

@testset "Allocations — writing tables" begin
    ma = MomentArray([1.0 2.0; 3.0 4.0; 5.0 6.0], (:a, :b, :c), (:mean, :sd))
    @test allocs_after_warmup(Tables.columns, ma) == 0
    # panels: no more than `vec` of the stored arrays (a new array header per field, no data;
    # linear views of a Matrix allocate the same), and independent of the panel's size
    for panel in (simseries, errs)
        @test allocs_after_warmup(Tables.columns, panel) <= allocs_after_warmup(d -> map(vec, components(d)), panel.data)
        @test allocs_after_warmup(Tables.columns, panel[1:10, 1:1]) == allocs_after_warmup(Tables.columns, panel)
    end
    # the lazy index columns are read without allocating (row names allocate one String per access)
    @test allocs_after_warmup(read_lazy_indices, Tables.columns(simseries)) == 0

    # IRF tables look values up in the stored lines and panels: their cost does not grow with
    # the number of draws
    irf_sizes = map(N -> simulate_irf(obj, simseries; shocks = (pos = 5, neg = 1), path_vars = (:K, :C), T = 8, N), (40, 80))
    @test allocs_after_warmup(Tables.columns, irf_sizes[1]) == allocs_after_warmup(Tables.columns, irf_sizes[2])
    paths_columns(irf) = Tables.columns(irf_paths(irf))
    @test allocs_after_warmup(paths_columns, irf_sizes[1]) == allocs_after_warmup(paths_columns, irf_sizes[2])
    @test length(paths_columns(irf_sizes[2]).__run) == 2 * length(paths_columns(irf_sizes[1]).__run)
end

@testset "Allocations — restoring a MomentArray" begin
    R, C = length(ALLOC_ROWS), length(ALLOC_COLS)
    ma = MomentArray(reshape(collect(1.0:R*C), R, C), ALLOC_ROWS, ALLOC_COLS)

    # the static target is filled from the columns; restored values match
    @test restore_static(Tables.columntable(ma)).data isa SMatrix{R, C, Float64}
    @test restore_static(Tables.columntable(ma)) == ma
    @test restore_dense(Tables.columntable(ma)).data isa Matrix{Float64}

    # in-memory source: nothing (static) or exactly one matrix (dense)
    nt = (__row = collect(map(String, ALLOC_ROWS)), map(j -> ma.data[:, j], NamedTuple{ALLOC_COLS}(1:C))...)
    @test allocs_after_warmup(restore_static, nt) == 0
    @test allocs_after_warmup(restore_dense, nt) == allocs_after_warmup((r, c) -> Matrix{Float64}(undef, r, c), R, C)
    # a source with metadata: the same restore, values checked (parquet: OutputVault's tests)
    @test restore_static(roundtrip(ma)) == ma && restore_dense(roundtrip(ma)) == ma
end
