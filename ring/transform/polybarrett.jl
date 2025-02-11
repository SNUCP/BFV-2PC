# We use Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
# In this implementation, Qₛₚ = (Xᵐ-1)/(X^(m/p)-1) for the smallest prime factor p.
# In the original paper, the authors used polynomials with smaller degree than ours.
# However, for a generalised code, we shall stick to these polynomials.
mutable struct Reductor
    const p::Int64
    const moverp::Int64
    const iseasy::Bool
    const m::Int64
    const N::Int64
    const Q::Modulus
    const α::Int64
    const A::Int64
    const ñ::Int64
    const ntterñ::Union{CyclicTransformer_2a3b5c7d,Missing}
    const ntterA::Union{CyclicTransformer_2a3b5c7d,Missing}
    const ϕm::Vector{UInt64}
    const Xm_over_ϕm::Vector{UInt64}
    const fbuff::Vector{UInt64}
    const rbuff::Vector{UInt64}

    function Reductor(m::Int64, Q::Modulus)
        N = totient(m)
        factors = factor(Vector, m)
        p = factors[1]
        moverp = m ÷ p
        deg = m - moverp

        if deg == N
            new(p, moverp, true, m, N, Q, 0, 0, 0, missing, missing, [], [], [], [])
        else
            ñ = next2a3b5c7d(m)
            α = deg - N
            A = next2a3b5c7d(2α + 1)

            ntterñ = CyclicTransformer_2a3b5c7d(ñ, Q)
            ntterA = CyclicTransformer_2a3b5c7d(A, Q)

            ϕm = cyclotomic_finder(m)
            Xm_over_ϕm = division(vcat(zeros(Int64, N + α), 1), ϕm)

            ϕm = Mform.(zeropadto(ϕm, ñ), Ref(Q))
            Xm_over_ϕm = Mform.(zeropadto(Xm_over_ϕm, A), Ref(Q))

            ntt!(ϕm, ntterñ)
            ntt!(Xm_over_ϕm, ntterA)

            fbuff = Vector{UInt64}(undef, A)
            rbuff = Vector{UInt64}(undef, ñ)

            new(p, moverp, false, m, N, Q, α, A, ñ, ntterñ, ntterA, ϕm, Xm_over_ϕm, fbuff, rbuff)
        end
    end
end

"""
Barrett computes a (mod Φₘ).
"""
function Barrett(a::Vector{UInt64}, rdtor::Reductor)
    p, moverp, N, m, α, ñ, Q = rdtor.p, rdtor.moverp, rdtor.N, rdtor.m, rdtor.α, rdtor.ñ, rdtor.Q

    # reduction by Qₛₚ
    @inbounds for j = 1:moverp, i = 0:p-2
        a[i*moverp+j] = sub(a[i*moverp+j], a[m-moverp+j], Q)
    end

    if !rdtor.iseasy
        # Compute f = ⌊c/Xⁿ⌋
        @inbounds @simd for i = 1:α
            rdtor.fbuff[i] = a[N+i]
        end
        @. rdtor.fbuff[α+1:end] = zero(UInt64)

        # Compute f = f × ⌊Xⁿ⁺ᵅ/ϕₘ⌋
        ntt!(rdtor.fbuff, rdtor.ntterA)
        Mmul_to!(rdtor.fbuff, rdtor.fbuff, rdtor.Xm_over_ϕm, Q)
        intt!(rdtor.fbuff, rdtor.ntterA)

        # Compute r = ⌊f/Xᵅ⌋
        @inbounds @simd for i = 1:α
            rdtor.rbuff[i] = rdtor.fbuff[α+i]
        end
        @. rdtor.rbuff[α+1:end] = zero(UInt64)

        # Compute r = r × ϕₘ (mod Xⁿ̃ + 1)
        ntt!(rdtor.rbuff, rdtor.ntterñ)
        Mmul_to!(rdtor.rbuff, rdtor.rbuff, rdtor.ϕm, Q)
        intt!(rdtor.rbuff, rdtor.ntterñ)

        # Compute a = a % (Xⁿ̃ + 1)
        @inbounds for i = 1:N+α-ñ
            a[i] = sub(a[i], a[ñ+i], Q)
        end

        # Compute a -= r
        @inbounds for i = 1:min(length(a), ñ)
            a[i] = sub(a[i], rdtor.rbuff[i], Q)
        end
    end
end