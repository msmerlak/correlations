using DrWatson
using Random
using Distributions
include(srcdir("correlations.jl"))

Random.seed!(123)

N = 200
M = 10

vectors = [randn(N) for _ in 1:M]
vectors = [v / norm(v) for v in vectors]
values = 1:10
scale = Symmetric(I + sum(λ * v * v' for (λ, v) in zip(values, vectors)))

D(ν) = MvTDist(ν, Matrix(scale))


true_corr = Bigsimr.cov2cor(cov(D(5)))
true_evectors = eigen(true_corr).vectors
true_evalues = eigen(true_corr).values

C = TrueCorrelation(true_corr, true_evectors, true_evalues)
histogram(C.evalues, xlabel = "True eigenvalues", ylabel = "Density", label = false)
savefig(plotsdir("true_eigenvalues.png"))