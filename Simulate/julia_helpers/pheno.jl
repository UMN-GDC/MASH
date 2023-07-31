using Distributions
using Missings
using LinearAlgebra
using Statistics

function simPheno(Xg, h2)
  nSubjects, nSNPs = size(Xg)
  freq = sum(Xg, dims = 1) / (2 * nSubjects)
  Z = (Xg .- mean(Xg, dims = 2)) ./ .√(2 .* freq .* (1 .- freq))
  
  Beta = rand(Normal(0, sqrt(h2 / nSNPs)), nSNPs)
  eps = rand(Normal(0,sqrt(1-h2)), nSubjects)
  Y = Z * Beta +  eps
  return Y  
end

