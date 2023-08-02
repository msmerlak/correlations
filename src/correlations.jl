using LinearAlgebra
using Bigsimr
using Distances

resolvent(z, M) = inv(z * I - M)
stieljes(z, evalues) = mean([1 / (z - l) for l in evalues])
ledoit_peche(ξ, q, evalues) = ξ / abs(1 - q + q * ξ * stieljes(ξ - im / sqrt(size(M, 1)), evalues))^2

abstract type CorrelationMatrix end

struct TrueCorrelation <: CorrelationMatrix
    matrix
    evectors
    evalues
end

struct EmpiricalCorrelation <: CorrelationMatrix
    data
    q
    evectors
    evalues
    matrix
end

struct PearsonClipped end
struct SpearmanUncorrected end
struct KendallUncorrected end

METHODS = (Pearson, PearsonClipped, Kendall, Spearman, SpearmanUncorrected, KendallUncorrected)

function corr_estimator(
    X::AbstractMatrix, method::DataType;
    rie=true,
    cor_near=false
)
    q = size(X, 2) / size(X, 1)
    if method == Pearson
        corr = cor(X, Pearson)
    elseif method in (Spearman, Kendall)
        corr = cor_convert(cor(X, method), method, Pearson)
        cor_near && (corr = cor_nearPD(corr))
    elseif method == SpearmanUncorrected
        corr = cor(X, Spearman)
    elseif method == KendallUncorrected
        corr = cor(X, Kendall)
    elseif method == PearsonClipped
        flat = collect(Iterators.flatten(X))
        m = quantile(flat, 0.1)
        M = quantile(flat, 0.9)
        corr = cor(clamp.(X, m, M), Pearson)
    end
    spectrum = eigen(Symmetric(corr))
    L = spectrum.values
    U = spectrum.vectors
    if rie
        L = ledoit_peche.(L, Ref(q), Ref(L))
        corr = U * Diagonal(L) * U'
    end
    return EmpiricalCorrelation(X, q, U, L, corr)
end

function oracle(E::EmpiricalCorrelation, C::TrueCorrelation)
    U = E.evectors
    return U * Diagonal([u' * C.matrix * u for u in eachcol(U)]) * U'
end

function outliers(evalues)
    N = length(evalues)
    D = pairwise(Euclidean(), evalues)
    D[diagind(D)] .= Inf
    indices = eachindex(evalues)
    bulk_condition = map(minimum, eachcol(D)) .< 1 / sqrt(N)
    bulk = indices[bulk_condition]
    outliers = indices[.!bulk_condition]
    return outliers
end



function outliers_overlaps(A::CorrelationMatrix, B::CorrelationMatrix)
    O = intersect(outliers(A.evalues), outliers(B.evalues))
    return [dot(A.evectors[:, i], B.evectors[:, i])^2 for i in reverse(O)]
end

function risk(E::EmpiricalCorrelation, C::TrueCorrelation)
    E_inv = inv(E.matrix)
    true_inv = inv(true_corr)
    return Dict(
        :in => 1 / tr(estimated_inv),
        :true => 1 / tr(true_inv),
        :out => tr(estimated_inv * true_corr * estimated_inv) / tr(estimated_inv)^2
    )
end
