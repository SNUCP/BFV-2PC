"""
RLWEkey is a struct for the RLWE secret key s ∈ R.
"""
struct RLWEkey
    coeffs::Vector{Int64}
    N::Int64
end

"""
Outputs secret key for RLWE.
If hw = 0, the secret is sampled from a uniform binary distribution.
Otherwise, it outputs a secret key with hamming weight at most hw. 
"""
binary_ringkey(us::UniformSampler, N::Int64, hw::Int64 = 0) = RLWEkey(uniform_binary(us, N, hw), N)

"""
Outputs secret key for RLWE.
If hw = 0, the secret is sampled from a uniform ternary distribution.
Otherwise, it outputs a secret key with hamming weight at most hw. 
"""
ternary_ringkey(us::UniformSampler, N::Int64, hw::Int64 = 0) = RLWEkey(uniform_ternary(us, N, hw), N)

"""
RLWEkeyQ is a struct for the embedding of the RLWE secret key s ∈ Rₚ.
"""
const RLWEkeyQ = ModPoly

RLWEkeyQ(key::RLWEkey, moduli::Moduli) = begin
    coeffs = Array{UInt64, 2}(undef, key.N, length(moduli))
    for j = eachindex(moduli)
        @simd for i = 1 : key.N
            coeffs[i, j] = key.coeffs[i] ≥ 0 ? key.coeffs[i] : key.coeffs[i] + moduli[j].Q
        end
    end
    ModPoly(coeffs, false, false)
end