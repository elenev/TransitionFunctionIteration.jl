# Completely unimplemented model for testing interface stubs.
# Must be defined at module scope (not inside @testset) so the struct is valid.
struct UnimplModel <: AbstractModel end

@testset "Interface stubs — throw MethodNotImplementedError" begin

    mdl = UnimplModel()

    # define_endog_grids: second arg must be AbstractSteadyState (not Nothing)
    dummy_ss = SteadyState(model=mdl, guess=(;), endog=(;), exog=(;), values=(;))
    @test_throws TFI.MethodNotImplementedError TFI.define_endog_grids(mdl, dummy_ss)

    @test_throws TFI.MethodNotImplementedError TFI.steadystate(mdl, AnalyticalSteadyStateMode())
    @test_throws TFI.MethodNotImplementedError TFI.steadystate_system((;), mdl, SolveMode())

    # calc_eqm: fourth arg must be NextPeriod (not Nothing)
    dummy_np = NextPeriod(expectations=(;), continuations=(;))
    @test_throws TFI.MethodNotImplementedError TFI.calc_eqm(nothing, nothing, mdl, dummy_np, SolveMode())

    @test_throws TFI.MethodNotImplementedError TFI.compute_next_period_values(nothing, mdl, SolveMode())
    @test_throws TFI.MethodNotImplementedError Exogenv(mdl)

end

@testset "MethodNotImplementedError — construction and display" begin

    err = TFI.MethodNotImplementedError(:my_func, UnimplModel)
    @test err isa TFI.MethodNotImplementedError

    # showerror produces a human-readable string
    io = IOBuffer()
    showerror(io, err)
    msg = String(take!(io))
    @test occursin("my_func", msg)
    @test occursin("UnimplModel", msg)

    # Constructor with custom message
    err2 = TFI.MethodNotImplementedError(:other_func, UnimplModel, " for some context")
    io2 = IOBuffer()
    showerror(io2, err2)
    msg2 = String(take!(io2))
    @test occursin("some context", msg2)

end
