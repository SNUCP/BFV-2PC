"""
BasisExtender provides basis extension from Q to P.
"""
mutable struct BasisExtender
    const Q::Vector{Modulus}
    const P::Vector{Modulus}
    const Qtilde::Vector{UInt64}
    const Qstar::Array{UInt64,2}
    const QmodP::Vector{UInt64}
    const oneoverQ::Vector{UInt128}
    const buff::Array{UInt64,2}
    const v128buff::Vector{UInt128}
    const v64buff::Vector{UInt64}

    function BasisExtender(Q::Moduli, P::Moduli, N::Int64)
        Qlen, Plen = length(Q), length(P)

        Qtilde = Vector{UInt64}(undef, Qlen)
        Qstar = Array{UInt64,2}(undef, Qlen, Plen)
        QmodP = Vector{UInt64}(undef, Plen)
        oneoverQ = Vector{UInt128}(undef, Qlen)
        buff = Array{UInt64,2}(undef, N, Qlen)
        v128buff = Vector{UInt128}(undef, N)
        v64buff = Vector{UInt64}(undef, N)

        Qbig = prod(BigInt[Qi.Q for Qi = Q])

        @inbounds for i = 1:Qlen
            Qtilde[i] = invmod(Qbig ÷ Q[i].Q, Q[i].Q) % UInt64
        end

        @inbounds for j = 1:Plen
            for i = 1:Qlen
                Qstar[i, j] = Qbig ÷ Q[i].Q % P[j].Q % UInt64
            end
        end

        @inbounds for i = 1:Plen
            QmodP[i] = Qbig % P[i].Q % UInt64
        end

        setprecision(256)
        @inbounds @simd for i = 1:Qlen
            oneoverQ[i] = round(UInt128, (Int128(1) << fixed_prec) / Q[i].Q)
        end

        new(collect(Q), collect(P), Qtilde, Qstar, QmodP, oneoverQ, buff, v128buff, v64buff)
    end
end

"""
basis_extend! performs the HPS18 Basis Extension algorithm. The input should not be in Montgomery domain.
"""
@views function basis_extend!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, be::BasisExtender)
    Q, P = be.Q, be.P
    Qlen, Plen = length(Q), length(P)

    Qtilde, Qstar, QmodP, oneoverQ, buff = be.Qtilde, be.Qstar, be.QmodP, be.oneoverQ, be.buff

    @inbounds for i = 1:Qlen
        buff[i] = Bmul(a[i], Qtilde[i], Q[i])
    end

    v128 = zero(UInt128)
    @inbounds for i = 1:Qlen
        v128 += mult_and_round(buff[i], oneoverQ[i])
    end
    v64 = round_to_uint64(v128)

    @inbounds for j = 1:Plen
        res[j] = lazy_Bmul(neg(QmodP[j], P[j]), v64, P[j])
        for i = 1:Qlen-1
            res[j] = lazy_Bred(widemul(buff[i], Qstar[i, j]) + res[j], P[j])
        end
        res[j] = Bred(widemul(buff[Qlen], Qstar[Qlen, j]) + res[j], P[j])
    end
end

@views function basis_extend!(res::AbstractArray{UInt64,2}, a::AbstractArray{UInt64,2}, be::BasisExtender)
    Q, P = be.Q, be.P
    Qlen, Plen = length(Q), length(P)

    Qtilde, Qstar, QmodP, oneoverQ, buff, v128, v64 = be.Qtilde, be.Qstar, be.QmodP, be.oneoverQ, be.buff, be.v128buff, be.v64buff

    @inbounds for j = 1:Qlen, i = eachindex(v64)
        buff[i, j] = Bmul(a[i, j], Qtilde[j], Q[j])
    end

    @. v128 = 0
    @inbounds for i = 1:Qlen
        @. v128 += mult_and_round(buff[:, i], oneoverQ[i])
    end
    @. v64 = round_to_uint64(v128)

    @inbounds for j = 1:Plen
        # Decompose more efficiently.
        for i = 1:Qlen
            if P[j].Q == Q[i].Q
                @. res[:, j] = a[:, i]
                @goto endloop
            end
        end
        
        QmodPneg = neg(QmodP[j], P[j])
        for k = eachindex(v64)
            res[k, j] = lazy_Bmul(QmodPneg, v64[k], P[j])
        end
        for i = 1:Qlen-1, k = eachindex(v64)
            res[k, j] = lazy_Bred(widemul(buff[k, i], Qstar[i, j]) + res[k, j], P[j])
        end
        for k = eachindex(v64)
            res[k, j] = Bred(widemul(buff[k, Qlen], Qstar[Qlen, j]) + res[k, j], P[j])
        end

        @label endloop
    end
end

"""
SimpleScaler provides simple scaling from Q to T.
"""
mutable struct SimpleScaler
    const Q::Vector{Modulus}
    const T::Vector{Modulus}
    const ω::Array{UInt64,2}
    const θ::Vector{UInt128}
    const v128buff::Vector{UInt128}

    function SimpleScaler(Q::Moduli, T::Moduli, N::Int64)
        Qlen, Tlen = length(Q), length(T)

        ω = Array{UInt64}(undef, Qlen, Tlen)
        θ = Vector{UInt128}(undef, Qlen)
        v128buff = Vector{UInt128}(undef, N)

        Tbig = prod(BigInt[Ti.Q for Ti = T])
        Qbig = prod(BigInt[Qi.Q for Qi = Q])

        setprecision(64 * (Tlen + Qlen + 4))
        # Compute ω and θ
        @inbounds for i = 1:Qlen
            tmp = invmod(Qbig ÷ Q[i].Q, Q[i].Q) * Tbig / Q[i].Q
            ωi = floor(BigInt, tmp)
            for j = 1:Tlen
                ω[i, j] = (ωi % T[j].Q) % UInt64
            end
            θ[i] = round(UInt128, (Int128(1) << fixed_prec) * (tmp - ωi))
        end

        new(collect(Q), collect(T), ω, θ, v128buff)
    end
end

"""
simple_scale! performs the HPS18 simple scaling algorithm. The input should not be in Montgomery form.
"""
@views function simple_scale!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ss::SimpleScaler)
    Q, T = ss.Q, ss.T
    Qlen, Tlen = length(Q), length(T)

    ω, θ = ss.ω, ss.θ

    # Compute w.
    @. res = 0
    @inbounds for i = 1:Qlen, j = 1:Tlen
        res[j] = lazy_Bred(widemul(a[i], ω[i, j]) + res[j], T[j])
    end

    v128 = zero(UInt128)
    @inbounds for i = 1:Qlen
        v128 += mult_and_round(a[i], θ[i])
    end
    v128 = round_to_uint128(v128)

    @inbounds for i = 1:Tlen
        res[i] = Bred(v128 + res[i], T[i])
    end
end

"""
simple_scale! performs the HPS18 simple scaling algorithm. The input should not be in Montgomery form.
"""
@views function simple_scale!(res::AbstractArray{UInt64,2}, a::AbstractArray{UInt64,2}, ss::SimpleScaler)
    Q, T = ss.Q, ss.T
    Qlen, Tlen = length(Q), length(T)
    v128 = ss.v128buff

    ω, θ = ss.ω, ss.θ

    # Compute w.
    @. res = 0
    @inbounds for i = 1:Qlen, j = 1:Tlen, k = eachindex(v128)
        res[k, j] = lazy_Bred(widemul(a[k, i], ω[i, j]) + res[k, j], T[j])
    end

    @. v128 = 0
    @inbounds for i = 1:Qlen
        @. v128 += mult_and_round(a[:, i], θ[i])
    end
    @. v128 = round_to_uint128(v128)

    @inbounds for i = 1:Tlen, j = eachindex(v128)
        res[j, i] = Bred(v128[j] + res[j, i], T[i])
    end
end

"""
ComplexScaler provides rounding from PQ to Q, with a factor T/Q.
"""
mutable struct ComplexScaler
    const Q::Vector{Modulus}
    const P::Vector{Modulus}
    const T::Vector{Modulus}
    const PQtilde::Vector{UInt64}
    const oneoverPQ::Vector{UInt128}
    const ω::Array{UInt64,2}
    const θ::Vector{UInt128}
    const ζ::Array{UInt64,2}
    const λ::Vector{UInt64}
    const buffPQ::Array{UInt64,2}
    const v128buff::Vector{UInt128}
    const v64buff::Vector{UInt64}

    function ComplexScaler(Q::Moduli, P::Moduli, T::Moduli, N::Int64)
        Qlen, Plen, Tlen = length(Q), length(P), length(T)

        PQtilde = Vector{UInt64}(undef, Qlen + Plen)
        oneoverPQ = Vector{UInt128}(undef, Qlen + Plen)
        ω = Array{UInt64}(undef, Qlen, Qlen)
        θ = Vector{UInt128}(undef, Qlen)
        ζ = Array{UInt64}(undef, Plen, Qlen)
        λ = Vector{UInt64}(undef, Qlen)
        buffPQ = Array{UInt64,2}(undef, N, Qlen + Plen)
        v128buff = Vector{UInt128}(undef, N)
        v64buff = Vector{UInt64}(undef, N)

        Pbig = prod(BigInt[Pi.Q for Pi = P])
        Tbig = prod(BigInt[Ti.Q for Ti = T])
        Qbig = prod(BigInt[Qi.Q for Qi = Q])

        setprecision(256)
        @inbounds for i = 1:Qlen
            PQtilde[i] = invmod(Qbig * Pbig ÷ Q[i].Q, Q[i].Q) % UInt64
            oneoverPQ[i] = round(UInt128, (Int128(1) << fixed_prec) / Q[i].Q)
        end
        @inbounds for i = 1:Plen
            PQtilde[Qlen+i] = invmod(Qbig * Pbig ÷ P[i].Q, P[i].Q) % UInt64
            oneoverPQ[Qlen+i] = round(UInt128, (Int128(1) << fixed_prec) / P[i].Q)
        end

        setprecision(64 * (Tlen + Plen + Qlen + 4))
        @inbounds for i = 1:Qlen
            tmp = Pbig * Tbig / Q[i].Q
            ωi = floor(BigInt, tmp)
            for j = 1:Qlen
                ω[i, j] = (ωi % Q[j].Q) % UInt64
            end
            θ[i] = round(UInt128, (Int128(1) << fixed_prec) * (tmp - ωi))
            λ[i] = (Pbig * Tbig % Q[i].Q) % UInt64
        end

        @inbounds for j = 1:Plen
            tmp = Pbig * Tbig ÷ P[j].Q
            for i = 1:Qlen
                ζ[j, i] = (tmp % Q[i].Q) % UInt64
            end
        end

        new(collect(Q), collect(P), collect(T), PQtilde, oneoverPQ, ω, θ, ζ, λ, buffPQ, v128buff, v64buff)
    end
end

"""
complex_scale! performs the HPS18 complex scaling algorithm. The input should not be in Montgomery form.
"""
@views function complex_scale!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, cs::ComplexScaler)
    Q, P = cs.Q, cs.P
    Qlen, Plen = length(Q), length(P)

    PQtilde, oneoverPQ, ω, θ, ζ, λ = cs.PQtilde, cs.oneoverPQ, cs.ω, cs.θ, cs.ζ, cs.λ
    buffPQ = cs.buffPQ

    v128 = UInt128(0)
    @inbounds for i = 1:Qlen
        buffPQ[i] = Bmul(a[i], PQtilde[i], Q[i])
        v128 += mult_and_round(buffPQ[i], oneoverPQ[i])
    end
    @inbounds for i = 1:Plen
        buffPQ[Qlen+i] = Bmul(a[Qlen+i], PQtilde[Qlen+i], P[i])
        v128 += mult_and_round(buffPQ[Qlen+i], oneoverPQ[Qlen+i])
    end
    v64 = round_to_uint64(v128)

    v128 = UInt128(0)
    @inbounds for i = 1:Qlen
        v128 += mult_and_round(buffPQ[i], θ[i])
    end
    v128 = round_to_uint128(v128)

    @inbounds for i = 1:Qlen
        res[i] = Bred(v128 + widemul(Q[i].Q - v64, λ[i]), Q[i])
        for j = 1:Qlen
            res[i] = Bred(widemul(buffPQ[j], ω[j, i]) + res[i], Q[i])
        end
        for j = 1:Plen
            res[i] = Bred(widemul(buffPQ[Qlen+j], ζ[j, i]) + res[i], Q[i])
        end
    end
end

"""
complex_scale! performs the HPS18 complex scaling algorithm. The input should not be in Montgomery form.
"""
@views function complex_scale!(res::AbstractArray{UInt64,2}, a::AbstractArray{UInt64,2}, cs::ComplexScaler)
    Q, P = cs.Q, cs.P
    Qlen, Plen = length(Q), length(P)

    PQtilde, oneoverPQ, ω, θ, ζ, λ = cs.PQtilde, cs.oneoverPQ, cs.ω, cs.θ, cs.ζ, cs.λ
    buffPQ, v128, v64 = cs.buffPQ, cs.v128buff, cs.v64buff

    @. v128 = 0
    @inbounds for j = 1:Qlen, i = eachindex(v128)
        buffPQ[i, j] = Bmul(a[i, j], PQtilde[j], Q[j])
        v128[i] += mult_and_round(buffPQ[i, j], oneoverPQ[j])
    end
    @inbounds for j = 1:Plen, i = eachindex(v128)
        buffPQ[i, Qlen+j] = Bmul(a[i, Qlen+j], PQtilde[Qlen+j], P[j])
        v128[i] += mult_and_round(buffPQ[i, Qlen+j], oneoverPQ[Qlen+j])
    end
    @. v64 = round_to_uint64(v128)

    @. v128 = 0
    @inbounds for j = 1:Qlen, i = eachindex(v128)
        v128[i] += mult_and_round(buffPQ[i, j], θ[j])
    end
    @. v128 = round_to_uint128(v128)

    @inbounds for j = 1:Qlen
        for i = eachindex(v128)
            res[i, j] = lazy_Bred(v128[i] + widemul(Q[j].Q - v64[i], λ[j]), Q[j])
        end
        for k = 1:Qlen
            for i = eachindex(v128)
                res[i, j] = lazy_Bred(widemul(buffPQ[i, k], ω[k, j]) + res[i, j], Q[j])
            end
        end
        for k = 1:Plen
            for i = eachindex(v128)
                res[i, j] = Bred(widemul(buffPQ[i, Qlen+k], ζ[k, j]) + res[i, j], Q[j])
            end
        end
    end
end