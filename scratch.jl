using TransitionFunctionIteration
using StructArrays, StatsBase

function create_test_simvector()
    X = rand(10_000)
    Y = rand(10_000)
    simvector = StructVector((X=X, Y=Y))
    return simvector
end

function create_test_structarray()
    sv1 = create_test_simvector()
    sv2 = create_test_simvector()
    sa = hcat(sv1, sv2)  # This creates a StructArray of StructVectors
    return permutedims(sa,(2,1))
end

teststatfcn = (mean = mean, std = std, ac = TransitionFunctionIteration.ac1);
tss = TransitionFunctionIteration.StatisticsSet(teststatfcn);
mss = TransitionFunctionIteration.StatisticsSet((:corY,), (TransitionFunctionIteration.MultivariateStatFun{:Y}(cor),));
sm = create_test_structarray();
mask = trues(size(sm));

function run_many_times(f, N)
    tmp = ntuple(i -> f(), N)
    return nothing
end

ds() = TransitionFunctionIteration._describe_sample(sm, tss, mss, mask);

runstats = Matrix{Float64}(undef,2,4);
sv = sm[1,:];
dr!() = TransitionFunctionIteration._describe_run!(runstats, sv, tss, mss)