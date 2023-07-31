using DataFrames

struct phenoSimulation
  genotypes::Vector{Int64}
  df::DataFrame
  nsubjects::Int64
  simSettings::Dict{String,Any}
end
Genes = Array{Int64}([1,2,3])
df = DataFrame(Genes = Genes, A = [1,2,3], B = [4,5,6])

sampleSize(p::phenoSimulation) = p.nsubjects

settings = Dict([("A", 1), ("B", 2)])
X = phenoSimulation(Genes, df, 30, settings)

