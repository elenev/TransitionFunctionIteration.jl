using Logging: with_logger, @logmsg

@testset "ExecutionLogger — terminal mode" begin

    buf = IOBuffer()
    logger = TerminalLogger(ITERATION; io=buf, varname_width=10, gridpoint_widths=(exst=1, K=5, Lev=5))
    @test TFI.Logging.min_enabled_level(logger) == ITERATION
    @test logger isa TerminalLogger

    with_logger(logger) do
        @logmsg ITERATION "" iter=1 max_vdist=-3.5 max_tdist=-4.2 worst_vvar=:K worst_tvar=:Lev worst_vpoint=(exst=5, K=1.234, Lev=0.567) worst_tpoint=(exst=2, K=1.234, Lev=0.567) worst_texst=2
        @logmsg COMPLETE  "" msg="Converged."
    end

    out = String(take!(buf))
    @test contains(out, "Max V Dist")
    @test contains(out, "Max T Dist")
    @test contains(out, ":K")
    @test contains(out, "(:Lev, 2)")
    @test contains(out, "COMPLETE")
    @test contains(out, "Converged.")
    @test count(==("Max V Dist"), [m.match for m in eachmatch(r"Max V Dist", out)]) == 1

end

@testset "ExecutionLogger — CSV mode" begin

    tmp = tempname() * ".csv"
    logger = CSVLogger(tmp; min_level=ITERATION)
    @test logger isa CSVLogger

    with_logger(logger) do
        @logmsg ITERATION "" iter=1 max_vdist=-3.0 max_tdist=-4.0 worst_vvar=:C worst_tvar=:K worst_vpoint=(K=1.0,) worst_tpoint=(K=1.0,) worst_texst=1
        @logmsg ITERATION "" iter=2 max_vdist=-4.5 max_tdist=-5.1
        @logmsg COMPLETE  "" iter=2 max_vdist=-4.5 max_tdist=-5.1 msg="Converged."
    end
    close(logger)

    lines = readlines(tmp)
    @test lines[1] == TFI.CSV_HEADER
    @test count(l -> l == TFI.CSV_HEADER, lines) == 1
    @test length(lines) == 4
    @test startswith(lines[2], "ITERATION,")
    @test startswith(lines[4], "COMPLETE,")
    @test contains(lines[2], "-3")
    @test !contains(lines[2], "-3.00000")
    @test contains(lines[4], "Converged.")

end

@testset "ExecutionLogger — FAILEDPOINT in CSV mode" begin

    tmp = tempname() * ".csv"
    logger = CSVLogger(tmp; min_level=FAILEDPOINT)

    with_logger(logger) do
        @logmsg FAILEDPOINT "" fp_state="(exst=1, x=0.9)" fp_retcode="Failure"
    end
    close(logger)

    lines = readlines(tmp)
    @test lines[1] == TFI.CSV_HEADER
    @test startswith(lines[2], "FAILEDPOINT,")
    @test contains(lines[2], "Failure")

end

@testset "ExecutionLogger — min_level filtering" begin

    buf = IOBuffer()
    logger = TerminalLogger(COMPLETE; io=buf)

    with_logger(logger) do
        @logmsg ITERATION "" iter=1 max_vdist=-3.0 max_tdist=-4.0 worst_vvar=:K worst_tvar=:K worst_vpoint=(K=1.0,) worst_tpoint=(K=1.0,) worst_texst=nothing
        @logmsg COMPLETE  "" msg="done"
    end

    out = String(take!(buf))
    @test !contains(out, "Max V Dist")
    @test contains(out, "COMPLETE")
    @test contains(out, "done")

end

@testset "ExecutionLogger — shouldlog" begin

    # When min_level=ITERATION, levels above should pass, below should not
    logger = TerminalLogger(ITERATION; io=devnull)

    @test TFI.Logging.shouldlog(logger, COMPLETE,        nothing, nothing, nothing) == true
    @test TFI.Logging.shouldlog(logger, ITERATION,       nothing, nothing, nothing) == true
    @test TFI.Logging.shouldlog(logger, ITERATIONFAILED, nothing, nothing, nothing) == false
    @test TFI.Logging.shouldlog(logger, FAILEDPOINT,     nothing, nothing, nothing) == false

end

@testset "ExecutionLogger — ITERATIONFAILED message" begin

    buf = IOBuffer()
    logger = TerminalLogger(ITERATIONFAILED; io=buf)

    with_logger(logger) do
        @logmsg ITERATIONFAILED "" pass=1 n_failed=3 n_total=50
    end

    out = String(take!(buf))
    @test contains(out, "3")
    @test contains(out, "50")

    # Pass > 1 format
    buf2 = IOBuffer()
    logger2 = TerminalLogger(ITERATIONFAILED; io=buf2)
    with_logger(logger2) do
        @logmsg ITERATIONFAILED "" pass=2 n_solved=2 n_unsolved=1
    end
    out2 = String(take!(buf2))
    @test contains(out2, "Pass")

end

@testset "ExecutionLogger — TerminalLogger(level, obj) constructor" begin

    # The second TerminalLogger overload infers widths from a ModelSolution.
    # Uses the DummyModel/obj fixture defined in runtests.jl.
    logger = TerminalLogger(ITERATION, obj; io=devnull)
    @test logger isa TerminalLogger
    @test TFI.Logging.min_enabled_level(logger) == ITERATION
    # varname_width should be positive (inferred from Vfct/Tfct variable names)
    @test logger.varname_width >= 1

end

@testset "ExecutionLogger — close" begin

    # close(TerminalLogger) is a no-op; should not error
    logger_t = TerminalLogger(ITERATION; io=devnull)
    @test close(logger_t) === nothing

    # close(CSVLogger) closes the underlying IO
    tmp = tempname() * ".csv"
    logger_c = CSVLogger(tmp; min_level=ITERATION)
    @test isopen(logger_c.io) == true
    close(logger_c)
    @test isopen(logger_c.io) == false

end

@testset "_fmt_width — scalar methods" begin
    @test TFI._fmt_width(5)    == 1
    @test TFI._fmt_width(123)  == 3
    @test TFI._fmt_width(0)    == 1
    @test TFI._fmt_width(2.5)  == 5      # 1 integer digit + 4
    @test TFI._fmt_width(12.3) == 6
    @test TFI._fmt_width(0.1)  == 5      # integer part floored at 1 digit
    @test TFI._fmt_width(-2.5) == 6      # extra character for the sign

    # Scalar and one-element-array methods agree
    for x in (7, 123, 2.5, -2.5, 12.3)
        @test TFI._fmt_width(x) == TFI._fmt_width([x])
    end
end

@testset "_gridpoint_widths and TerminalLogger — singleton grid axis" begin
    ssg_s = StateSpaceGrid(; exst=1:3, x=collect(1.0:3.0), z=2.5)   # z is a scalar singleton
    @test singletons(ssg_s) == (:z,)

    widths = TFI._gridpoint_widths(ssg_s)    # used to throw: no _fmt_width method for scalars
    @test widths == (exst = 1, x = 5, z = 5)

    buf = IOBuffer()
    logger = TerminalLogger(ITERATION; io=buf, varname_width=10, gridpoint_widths=widths)
    with_logger(logger) do
        @logmsg ITERATION "" iter=1 max_vdist=-3.5 max_tdist=-4.2 worst_vvar=:x worst_tvar=:x worst_vpoint=(exst=2, x=1.5, z=2.5) worst_tpoint=(exst=2, x=1.5, z=2.5) worst_texst=2
    end
    @test occursin("2.500", String(take!(buf)))
end
