mutable struct CyclotomicTransformer_pow2
    const Q::Modulus
    const m::Int64
    const N::Int64
    const N⁻¹::UInt64
    const Ψ::Vector{UInt64}
    const Ψinv::Vector{UInt64}

    function CyclotomicTransformer_pow2(N::Int64, Q::Modulus)
        @assert Q.Q & (2N-1) == 1 "ϕ($(Q.Q)) does not have enough 2 factor to perform NTT."

        ξ = primitive_root_finder(Q)
        ζ = powermod(ξ, (Q.Q - 1) ÷ 2N, Q)
        ζinv = powermod(ζ, 2N - 1, Q)
        Ψ, Ψinv = Vector{UInt64}(undef, N), Vector{UInt64}(undef, N)
        Ψ[1], Ψinv[1], Ψ[2], Ψinv[2] = Mform(1, Q), Mform(1, Q), ζ, ζinv
        @inbounds for i = 3:N
            Ψ[i], Ψinv[i] = Mmul(ζ, Ψ[i-1], Q), Mmul(ζinv, Ψinv[i-1], Q)
        end
        scramble!(Ψ, 2)
        scramble!(Ψinv, 2)

        new(Q, 2N, N, Mform(invmod(N, Q), Q), Ψ, Ψinv)
    end
end

@views ntt!(a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_pow2) = begin
    ntt_2a3b5c7d!(a, ntter.Ψ, ntter.Ψ, ntter.Ψ, ntter.Ψ, ntter.Q)

    @inbounds @simd for i = eachindex(a)
        a[i] ≥ ntter.Q.Q && (a[i] -= ntter.Q.Q)
    end
end

@views ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_pow2) = begin
    @inbounds @simd for i = eachindex(res)
        res[i] = a[i]
    end
    ntt!(res, ntter)
end

@views intt!(a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_pow2) = begin
    intt_2a3b5c7d!(a, ntter.Ψinv, ntter.Ψinv, ntter.Ψinv, ntter.Ψinv, ntter.Q)

    @inbounds for i = eachindex(a)
        a[i] = lazy_Mmul(a[i], ntter.N⁻¹, ntter.Q)
        a[i] ≥ ntter.Q.Q && (a[i] -= ntter.Q.Q)
    end
end

@views intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_pow2) = begin
    @inbounds @simd for i = eachindex(res)
        res[i] = a[i]
    end
    intt!(res, ntter)
end

# We use Bluestein NTT for the sake of the low space complexity. 
# This approach is somewhat similar to HElib, although we simply apply Bluestein NTT of degree m.
# In the original cuHE and polynomial Barrett reduction paper, the authors make use of power-of-two NTT instead.
# Although the time complexity is at least half, the space complexity of their approach is as twice as ours.
mutable struct CyclotomicTransformer_arb
    const Q::Modulus
    const m::Int64
    const N::Int64
    const ntter::CyclicTransformer_bluestein
    const rdtor::Reductor
    const buff::Vector{UInt64}
    const nttidx::Vector{Int64}
    const autidxset::Dict{Int64, Vector{Int64}}
    const dims::Vector{Int64}
    const gens::Vector{Int64}

    function CyclotomicTransformer_arb(m::Int64, Q::Modulus)
        @assert isodd(m) "Even cyclotomic degree is not supported, unless power of two."
        
        dims, gens = find_generators_mod_m(m)
        @assert length(dims) < 5 "Cyclotomic degree with factors more than five is not supported."

        ntter = CyclicTransformer_bluestein(m, Q)

        rdtor = Reductor(m, Q)

        buff = Vector{UInt64}(undef, m)

        nttidx = Array{Int64}(undef, dims...)
        autidxset = Dict{Int64, Vector{Int64}}()
        for i = Iterators.product([1:dimi for dimi = dims]...)
            nttidx[i...] = prod([powermod(gens[j], i[j]-1, m) for j = eachindex(gens)]) % m + 1
            autidxset[nttidx[i...]-1] = [i...]
        end
        nttidx = reshape(nttidx, rdtor.N)

        new(Q, m, rdtor.N, ntter, rdtor, buff, nttidx, autidxset, dims, gens)
    end
end

# Bluestein NTT
@views function ntt!(a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_arb)
    @. ntter.buff = zero(UInt64)
    @inbounds @simd for i = 1:ntter.N
        ntter.buff[i] = a[i]
    end

    ntt!(ntter.buff, ntter.ntter)
    @inbounds @simd for i = 1:ntter.N
        a[i] = ntter.buff[ntter.nttidx[i]]
    end
end

@views ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_arb) = begin
    @inbounds @simd for i = eachindex(res)
        res[i] = a[i]
    end
    ntt!(res, ntter)
end

# Bluestein iNTT
@views function intt!(a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_arb, islazy::Bool=false)
    @. ntter.buff = 0
    for i = 1:ntter.N
        ntter.buff[ntter.nttidx[i]] = a[i]
    end

    intt!(ntter.buff, ntter.ntter)

    if islazy
        resize!(a, ntter.m)
        @. a = ntter.buff
    else
        Barrett(ntter.buff, ntter.rdtor)
        @inbounds @simd for i = 1:ntter.N
            a[i] = ntter.buff[i]
        end
    end
end

@views intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicTransformer_arb, islazy::Bool=false) = begin
    @inbounds @simd for i = eachindex(res)
        res[i] = a[i]
    end
    intt!(res, ntter, islazy)
end

const CyclotomicTransformer = Union{CyclotomicTransformer_pow2,CyclotomicTransformer_arb}

CyclotomicTransformer(m::Int64, Q::Modulus) =
    ispow2(m) ? CyclotomicTransformer_pow2(m ÷ 2, Q) : CyclotomicTransformer_arb(m, Q)