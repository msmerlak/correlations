using DrWatson;
@quickactivate;
begin
    using Plots
    gr(margin=5Plots.mm)
    using ThreadsX
    using Revise
    using DataFrames, StatsPlots

    includet(srcdir("correlations.jl"))
    includet(scriptsdir("distribution.jl"))

end


function overlaps_outliers(dist, Q; nmodes=10, repetitions=10)
    N = dist.dim
    T = Int(floor(Q * N))
    return Dict(
        Symbol(method) =>
            ThreadsX.collect(mean(dot(
                corr_estimator(rand(dist, T)', method).evectors[:, end-k],
                C.evectors[:, end-k]
            )^2 for _ in 1:repetitions)
                             for k in 0:(nmodes-1)) for method in METHODS) |> DataFrame
end

@time df = overlaps_outliers(D(3), 1.5; repetitions = 100)


@df df plot(cols(), xlabel = "Mode", ylabel = "Squared overlap")
savefig(plotsdir("overlaps.png"))