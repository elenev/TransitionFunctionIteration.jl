using Logging: with_logger

@testset "ExecutionLogger — terminal mode" begin

    # Colored mode writes to the provided IO
    buf = IOBuffer()
    logger = ExecutionLogger(ITERATION; io=buf)
    @test TFI.Logging.min_enabled_level(logger) == ITERATION
    @test !logger.csv_mode

    with_logger(logger) do
        @iteration iter=1 max_vdist=-3.0
        @complete iter=1 max_vdist=-3.0 msg="Converged."
    end

    out = String(take!(buf))
    @test contains(out, "ITERATION")
    @test contains(out, "COMPLETE")

end

@testset "ExecutionLogger — CSV mode" begin

    tmp = tempname() * ".csv"
    logger = ExecutionLogger(tmp; min_level=ITERATION)
    @test logger.csv_mode

    with_logger(logger) do
        @iteration iter=1 max_vdist=-3.0 max_tdist=-4.0 worst_vvar=:C worst_tvar=:K worst_vpoint="(K=1.0,)" worst_tpoint="(K=1.0,)" worst_texst="1"
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
    logger = ExecutionLogger(tmp; min_level=FAILEDPOINT)

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
    logger = ExecutionLogger(COMPLETE; io=buf)  # only COMPLETE and above

    with_logger(logger) do
        @iteration iter=1 max_vdist=-3.0   # should be suppressed
        @complete  iter=1 max_vdist=-3.0 msg="done"
    end

    out = String(take!(buf))
    @test !contains(out, "ITERATION")
    @test contains(out, "COMPLETE")

end
