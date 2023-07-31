using Distributions
using Missings
using LinearAlgebra
using Statistics

function sim_genos(nSubjects= 1000, nSNPs = 10000)
  X = Binomial(2, 0.5)
  Xg = missings(Float64, nSubjects, nSNPs)
  for i in 1:nSubjects
    Xg[i,:] = rand(X, nSNPs)
  end
  return Xg
end

