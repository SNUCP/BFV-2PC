include("./transform/cyclic.jl")
include("./transform/subring.jl")
include("./transform/polybarrett.jl")
include("./transform/cyclotomic.jl")

"""
Transformer is an abstract type which supports number theoretic transformations (NTT).
"""
const Transformer = Union{CyclicTransformer,SubringTransformer,CyclotomicTransformer}

(::Type{Transformer})(param::RingParam, Q::Modulus) = begin
    if typeof(param) == CyclicParam
        CyclicTransformer(param.m, Q)     
    elseif typeof(param) == CyclotomicParam
        CyclotomicTransformer(param.m, Q)
    elseif typeof(param) == SubringParam
        SubringTransformer(param.m, param.d, Q)
    end
end