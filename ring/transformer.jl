include("./transform/cyclic.jl")
include("./transform/polybarrett.jl")
include("./transform/cyclotomic.jl")

"""
Transformer is an abstract type which supports number theoretic transformations (NTT).
"""
const Transformer = Union{CyclicTransformer, CyclotomicTransformer}