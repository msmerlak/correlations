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



Q = 1.5
dist = D(3)
begin
    plts = []
    for rie in (true, false)
        p = plot(ylabel="Estimated eigenvalues", xlabel="True eigenvalues")
        for method in METHODS
            estimators = ThreadsX.collect(corr_estimator(rand(dist, Int(floor(Q * N)))', method; rie=rie) for _ in 1:3)
            mean_evals = mean([estimator.evalues for estimator in estimators])
            se_evals = std([estimator.evalues for estimator in estimators]) / length(estimators)
            plot!(p, true_evalues, mean_evals, ribbon=(se_evals, se_evals), marker=:auto, label=string(Symbol(method)))
        end
        plot!(p, x -> x, ls=:dash, color=:black, label=false)
        push!(plts, p)
    end
    plot(plts..., size=(1000, 500), margin=5Plots.mm)
end
savefig(plotsdir("eigenvalues.png"))