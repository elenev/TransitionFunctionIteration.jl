using Logging: with_logger
import TransitionFunctionIteration: @iteration, @complete, @iterationfailed, @failedpoint

@testset "ExecutionLogger — terminal mode" begin

    buf = IOBuffer()
    logger = TerminalLogger(ITERATION; io=buf, varname_width=10, gridpoint_widths=(exst=1, K=5, Lev=5))
    @test TFI.Logging.min_enabled_level(logger) == ITERATION
    @test logger isa TerminalLogger

    with_logger(logger) do
        @iteration iter=1 max_vdist=-3.5 max_tdist=-4.2 worst_vvar=:K worst_tvar=:Lev worst_vpoint=(exst=5, K=1.234, Lev=0.567) worst_tpoint=(exst=2, K=1.234, Lev=0.567) worst_texst=2
        @complete  msg="Converged."
    end

    out = String(take!(buf))
    @test contains(out, "Max V Dist")
    @test contains(out, "Max T Dist")
    @test contains(out, ":K")
    @test contains(out, "(:Lev, 2)")
    @test contains(out, "COMPLETE")
    @test contains(out, "Converged.")
    # COMPLETE no longer prints distance lines
    @test count(==("Max V Dist"), [m.match for m in eachmatch(r"Max V Dist", out)]) == 1

end

@testset "ExecutionLogger — CSV mode" begin

    tmp = tempname() * ".csv"
    logger = CSVLogger(tmp; min_level=ITERATION)
    @test logger isa CSVLogger

    with_logger(logger) do
        @iteration iter=1 max_vdist=-3.0 max_tdist=-4.0 worst_vvar=:C worst_tvar=:K worst_vpoint=(K=1.0,) worst_tpoint=(K=1.0,) worst_texst=1
        @iteration iter=2 max_vdist=-4.5 max_tdist=-5.1
        @complete iter=2 max_vdist=-4.5 max_tdist=-5.1 msg="Converged."
    end
    close(logger)

    lines = readlines(tmp)

    # Header written exactly once
    @test lines[1] == TFI.CSV_HEADER
    @test count(l -> l == TFI.CSV_HEADER, lines) == 1

    # Correct number of data rows (2 ITERATION + 1 COMPLETE)
    @test length(lines) == 4

    # First data row is ITERATION
    @test startswith(lines[2], "ITERATION,")

    # Last data row is COMPLETE
    @test startswith(lines[4], "COMPLETE,")

    # 4-sig-fig float formatting: -3.0 should appear as "-3" or "-3.0" (not "-3.00000")
    @test contains(lines[2], "-3")
    @test !contains(lines[2], "-3.00000")

    # msg field is populated for COMPLETE row
    @test contains(lines[4], "Converged.")

end

@testset "ExecutionLogger — FAILEDPOINT in CSV mode" begin

    tmp = tempname() * ".csv"
    logger = CSVLogger(tmp; min_level=FAILEDPOINT)

    with_logger(logger) do
        @failedpoint fp_state="(exst=1, x=0.9)" fp_retcode="Failure"
    end
    close(logger)

    lines = readlines(tmp)
    @test lines[1] == TFI.CSV_HEADER
    @test startswith(lines[2], "FAILEDPOINT,")
    @test contains(lines[2], "Failure")

end

@testset "ExecutionLogger — min_level filtering" begin

    buf = IOBuffer()
    logger = TerminalLogger(COMPLETE; io=buf)  # only COMPLETE and above

    with_logger(logger) do
        @iteration iter=1 max_vdist=-3.0 max_tdist=-4.0 worst_vvar=:K worst_tvar=:K worst_vpoint=(K=1.0,) worst_tpoint=(K=1.0,) worst_texst=nothing   # suppressed
        @complete  msg="done"
    end

    out = String(take!(buf))
    @test !contains(out, "Max V Dist")   # @iteration suppressed
    @test contains(out, "COMPLETE")
    @test contains(out, "done")

end
