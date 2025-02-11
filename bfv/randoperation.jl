mutable struct Randomizer
    us::UniformSampler
    rgs::RGSampler
    entor::PKencryptor
    buffRLWE::Vector{RLWE}
    buffQ::Vector{ModPoly}
    buffT::Vector{ModPoly}
    DiQ::Vector{Float64x2}
    DiPQ::Vector{Float64x2}
    Δ::ModScalar

    @views function Randomizer(oper::BFVoperator, entor::PKencryptor)
        us = UniformSampler()
        rgs = RGSampler()

        N, len = oper.buffQ[1].N, oper.buffQ[1].len

        buffRLWE = [RLWE(N, len) for _ = 1:2]
        buffQ = [ModPoly(N, len) for _ = 1:2]
        buffT = [ModPoly(N, 1) for _ = 1:2]

        PQ, Q, ΔQ = oper.evalPQ.moduli, oper.operQ.eval.moduli, oper.ecder.Δ

        Qlen, Qbig = length(Q), prod([big(Qi.Q) for Qi = Q])
        DiQ = Vector{Float64x2}(undef, Qlen)
        @inbounds for j = 1:Qlen
            Qtilde = Qbig ÷ Q[j].Q
            DiQ[j] = Float64x2(invmod(Qtilde, Q[j].Q) * Qtilde / Qbig)
        end

        PQlen, PQbig = length(PQ), prod([big(PQi.Q) for PQi = PQ])

        Δval = to_big(ΔQ, oper.operQ.eval)
        DiPQ = Vector{Float64x2}(undef, PQlen)
        @inbounds for j = 1:PQlen
            PQtilde = PQbig ÷ PQ[j].Q
            DiPQ[j] = Float64x2(invmod(PQtilde, PQ[j].Q) * PQtilde / Δval)
        end

        Δ = ModScalar(length(oper.evalPQ.moduli))
        copy!(Δ.vals[1:Qlen], ΔQ.vals)
        basis_extend!(Δ.vals[Qlen+1:end], ΔQ.vals, oper.beQ2P)

        new(us, rgs, entor, buffRLWE, buffQ, buffT, DiQ, DiPQ, Δ)
    end
end

rand_mul(x::ModPoly, y::RLWE, dgs::Vector{CDTSampler}, τ::Float64, rand::Randomizer, oper::BFVoperator) = begin
    res = deepcopy(y)
    rand_mul_to!(res, x, y, dgs, τ, rand, oper)
    res
end

@views function rand_mul_to!(res::RLWE, x::ModPoly, y::RLWE, dgs::Vector{CDTSampler}, τ::Float64, rand::Randomizer, oper::BFVoperator)
    evalQ, evalT, ecder, buffQ, buffT = oper.operQ.eval, oper.evalT, oper.ecder, oper.buffQ[1], oper.buffT
    @assert x.len == buffT.len == 1 "The input polynomial should be in the plaintext space."

    N, T = buffT.N, evalT.moduli[1]

    copy!(buffT, x)
    buffT.isntt[] && intt!(buffT, evalT)
    buffT.isMform[] && iMform!(buffT, evalT)
    @inbounds for i = 1:N
        I = sample(dgs[T.Q-buffT.coeffs[i, 1]]) * T.Q
        buffQ.coeffs[i, 1] += Bred(I, T)
    end

    extend_to!(buffQ, buffT, ecder)
    Mform!(buffQ, evalQ)
    ntt!(buffQ, evalQ)

    copy!(res, y)
    !res.b.isMform[] && Mform!(res.b, evalQ)
    !res.b.isntt[] && ntt!(res.b, evalQ)
    !res.a.isMform[] && Mform!(res.a, evalQ)
    !res.a.isntt[] && ntt!(res.a, evalQ)

    mul_to!(res.b, res.b, buffQ, evalQ)
    mul_to!(res.a, res.a, buffQ, evalQ)

    buffQ.isntt[] = false
    buffQ.isMform[] = false
    Q = evalQ.moduli
    if τ * 6 * √(2π) < 2.0^63
        @inbounds for i = 1:N
            ei = sample(rand.rgs, τ)
            @simd for j = eachindex(Q)
                buffQ.coeffs[i, j] = Bred(ei, Q[j])
            end
        end
    else
        setprecision(ceil(Int64, log2(τ)))
        @inbounds for i = 1:N
            ei = sample128(rand.rgs, τ)
            @simd for j = eachindex(Q)
                buffQ.coeffs[i, j] = Bred(ei, Q[j])
            end
        end
    end
    ntt!(buffQ, evalQ)
    add_to!(res.b, res.b, buffQ, evalQ)
end

rand_mul(x::RLWE, y::RLWE, rlk::RLEV, dgs::Vector{CDTSampler}, rand::Randomizer, oper::BFVoperator) = begin
    res = deepcopy(x)
    rand_mul_to!(res, x, y, rlk, dgs, rand, oper)
    res
end

@views function rand_mul_to!(res::RLWE, x::RLWE, y::RLWE, rlk::RLEV, dgslift::VarCenSampler, dgsswitch::VarCenSampler, τ_tensor::Float64, τ_modswitch::Float64, rand::Randomizer, oper::BFVoperator)
    # Lift the ciphertexts to R.
    rand_lift_to!(oper.mulbuffPQ[1:2], x, dgslift, rand, oper)
    rand_lift_to!(oper.mulbuffPQ[3:4], y, dgslift, rand, oper)

    # Tensoring.
    tensor_to!(oper.mulbuffPQ[5:7], oper.mulbuffPQ[1:2], oper.mulbuffPQ[3:4], oper)

    # Scale by 1/Delta.
    buffQ = oper.buffQ
    rand_scale_by_delta_to!(buffQ[1:3], oper.mulbuffPQ[5:7], dgsswitch, τ_tensor, τ_modswitch, rand, oper)
    copy!(res.b, buffQ[1])
    copy!(res.a, buffQ[2])

    # Key-switching
    buffRLWE, operQ, evalQ = oper.operQ.buffRLWE, oper.operQ, oper.operQ.eval
    rlwe_sample_to!(buffRLWE, rand.entor)
    ntt!(buffQ[3], evalQ)
    add_to!(buffQ[3], buffQ[3], buffRLWE.a, evalQ)
    ntt!(res, operQ)
    add_to!(res.a, res.a, buffRLWE.b, evalQ)
    gadgetprod_to!(buffRLWE, buffQ[3], rlk, operQ)
    add_to!(res, res, buffRLWE, operQ)
end

@views function to_float(x::AbstractVector{UInt64}, DiQ::Vector{Float64x2})
    res = Float64x2(.0)
    @inbounds @simd for i = eachindex(x)
        res += x[i] * DiQ[i]
    end
    (res - round(res._limbs[1]) - round(res._limbs[2]))._limbs[1]
end

@views function rand_lift_to!(res::AbstractVector{ModPoly}, x::RLWE, dgslift::VarCenSampler, rand::Randomizer, oper::BFVoperator)
    evalPQ, evalQ, beQ2P, buffQ = oper.evalPQ, oper.operQ.eval, oper.beQ2P, oper.buffQ[1]
    N, Qlen, Plen = x.b.N, length(beQ2P.Q), length(beQ2P.P)

    b, a = res[1], res[2]

    # Lift the ciphertexts to R.
    copy!(buffQ, x.b)
    copy!(b.coeffs[:, 1:Qlen], x.b.coeffs)
    buffQ.isntt[] && intt!(buffQ, evalQ)
    buffQ.isMform[] && iMform!(buffQ, evalQ)
    @inbounds for i = 1:N
        basis_extend!(b.coeffs[i, Qlen+1:end], buffQ.coeffs[i, :], beQ2P)
        bi = to_float(buffQ.coeffs[i, :], rand.DiQ)
        I = sample(dgslift, -bi)
        for j = 1:Plen
            b.coeffs[i, Qlen+j] += Bred(widemul(beQ2P.QmodP[j], I), evalPQ.moduli[Qlen+j]) 
        end
    end
    @inbounds for i = 1:Qlen
        (evalQ.ismod[i] || x.b.isMform[]) && continue
        Mform!(b.coeffs[:, i], evalQ.moduli[i])
        x.b.isntt[] && continue
        ntt!(b.coeffs[:, i], evalQ.ntter[i])
    end
    @inbounds for i = 1:Plen
        Mform!(b.coeffs[:, Qlen+i], evalPQ.moduli[Qlen+i])
        ntt!(b.coeffs[:, Qlen+i], evalPQ.ntter[Qlen+i])
    end
    b.isMform[] = true
    b.isntt[] = true

    copy!(buffQ, x.a)
    copy!(a.coeffs[:, 1:Qlen], x.a.coeffs)
    buffQ.isntt[] && intt!(buffQ, evalQ)
    buffQ.isMform[] && iMform!(buffQ, evalQ)
    @inbounds for i = 1:N
        basis_extend!(a.coeffs[i, Qlen+1:end], buffQ.coeffs[i, :], beQ2P)
        ai = to_float(buffQ.coeffs[i, :], rand.DiQ)
        I = sample(dgslift, -ai)
        for j = 1:Plen
            a.coeffs[i, Qlen+j] += Bred(widemul(beQ2P.QmodP[j], I), evalPQ.moduli[Qlen+j]) 
        end
    end
    @inbounds for i = 1:Qlen
        (evalQ.ismod[i] || x.a.isMform[]) && continue
        Mform!(a.coeffs[:, i], evalQ.moduli[i])
        x.a.isntt[] && continue
        ntt!(a.coeffs[:, i], evalQ.ntter[i])
    end
    @inbounds for i = 1:Plen
        Mform!(a.coeffs[:, Qlen+i], evalPQ.moduli[Qlen+i])
        ntt!(a.coeffs[:, Qlen+i], evalPQ.ntter[Qlen+i])
    end
    a.isMform[] = true
    a.isntt[] = true
end

@views function rand_scale_by_delta_to!(res::AbstractVector{ModPoly}, xvec::AbstractVector{ModPoly}, dgsswitch::VarCenSampler, τ_tensor::Float64, τ_modswitch::Float64, rand::Randomizer, oper::BFVoperator)
    evalPQ, cscaler, N = oper.evalPQ, oper.cscaler, res[1].N

    d0, d1, d2 = xvec[1], xvec[2], xvec[3]
    d̃0, d̃1, d̃2 = res[1], res[2], res[3]

    intt!(d0, evalPQ)
    intt!(d1, evalPQ)
    intt!(d2, evalPQ)

    iMform!(d0, evalPQ)
    iMform!(d1, evalPQ)
    iMform!(d2, evalPQ)

    d̃0.isntt[] = false
    d̃0.isMform[] = false
    d̃1.isntt[] = false
    d̃1.isMform[] = false
    d̃2.isntt[] = false
    d̃2.isMform[] = false

    PQ = evalPQ.moduli
    @inbounds for i = 1:N
        complex_scale!(d̃0.coeffs[i, :], d0.coeffs[i, :], cscaler)
        complex_scale!(d̃1.coeffs[i, :], d1.coeffs[i, :], cscaler)
        complex_scale!(d̃2.coeffs[i, :], d2.coeffs[i, :], cscaler)
        
        dd0, dd1, dd2 = to_float(d0.coeffs[i, :], rand.DiPQ), to_float(d1.coeffs[i, :], rand.DiPQ), to_float(d2.coeffs[i, :], rand.DiPQ)
        ed0, ed1, ed2 = sample(dgsswitch, -dd0), sample(dgsswitch, -dd1), sample(dgsswitch, -dd2)
        ẽi = sample(rand.rgs, τ_modswitch) + sample128(rand.rgs, τ_tensor)
        for j = 1:d̃0.len
            d̃0.coeffs[i, j] = add(d̃0.coeffs[i, j], Bred(ẽi + ed0, PQ[j]), PQ[j])
            d̃1.coeffs[i, j] = add(d̃1.coeffs[i, j], Bred(ed1, PQ[j]), PQ[j])
            d̃2.coeffs[i, j] = add(d̃2.coeffs[i, j], Bred(ed2, PQ[j]), PQ[j])
        end
    end
end

mask_randomize(x::RLWE, rand::Randomizer, oper::BFVoperator) = begin
    res = deepcopy(x)
    mask_randomize!(res, rand, oper)
    res
end

mask_randomize!(x::RLWE, rand::Randomizer, oper::BFVoperator) = mask_randomize_to!(x, x, rand, oper)

mask_randomize_to!(res::RLWE, x::RLWE, rand::Randomizer, oper::BFVoperator) = begin
    buff = oper.operQ.buffRLWE
    rlwe_sample_to!(buff, rand.entor)
    add_to!(res, x, buff, oper.operQ)
end

rand_automorphism(x::RLWE, idx::Int64, rand::Randomizer, oper::BFVoperator) = begin
    res = deepcopy(x)
    rand_automorphism!(res, idx, atk, rand, oper)
    res
end

rand_automorphism!(x::RLWE, idx::Int64, atk::RLEV, rand::Randomizer, oper::BFVoperator) = rand_automorphism_to!(x, x, idx, atk, rand, oper)

rand_automorphism_to!(res::RLWE, x::RLWE, idx::Int64, atk::RLEV, rand::Randomizer, oper::BFVoperator) = begin
    mask_randomize_to!(res, x, rand, oper)
    automorphism!(res, idx, atk, oper.operQ)
end