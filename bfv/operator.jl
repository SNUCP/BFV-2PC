"""
BFVoperator is a struct for arithmetic operations over BFV ciphertexts.
"""
mutable struct BFVoperator
    operQ::Operator
    evalPQ::PolyEvaluator
    evalT::PolyEvaluator
    ecder::BFVencoder
    beQ2P::BasisExtender
    cscaler::ComplexScaler
    mulbuffPQ::Vector{ModPoly}
    buffQ::Vector{ModPoly}
    buffT::ModPoly

    function BFVoperator(evalQ::PolyEvaluator, evalT::PolyEvaluator, decer::Decomposer)
        m, N, Q, T = evalQ.ntter[1].m, evalQ.ntter[1].N, evalQ.moduli, evalT.moduli

        # Setting P.
        Pprimes = collect(find_prime_cyclotomic(62, m, 2length(Q) + 1))
        filter!(x -> x ∉ [Qi.Q for Qi ∈ Q], Pprimes)
        Pidx, bigP, bigQ = 1, big(1), prod(BigInt[Qi.Q for Qi ∈ Q])
        while true
            bigP *= Pprimes[Pidx]
            if bigP > bigQ * N
                break
            end
            Pidx += 1
        end
        P = Modulus.(Pprimes[1:Pidx])

        lenQ, lenT, lenPQ = length(Q), length(T), length(Q) + length(P)

        PQ = vcat(Q, P)
        PQntter = vcat(evalQ.ntter, [CyclotomicTransformer(m, Pi) for Pi ∈ P])
        
        operQ = Operator(evalQ, decer)
        evalPQ = PolyEvaluator(PQ, PQntter)
        ecder = BFVencoder(evalQ, evalT)
        beQ2P = BasisExtender(Q, P, N)
        cscaler = ComplexScaler(Q, P, T, N)
        mulbuffPQ = [ModPoly(N, lenPQ) for _ = 1:7]
        buffQ = [ModPoly(N, lenQ) for _ = 1:3]
        buffT = ModPoly(N, lenT)

        new(operQ, evalPQ, evalT, ecder, beQ2P, cscaler, mulbuffPQ, buffQ, buffT)
    end
end

#TODO RLWE + ModScalar, RLWE - ModScalar, RLWE * ModScalar

add(x::RLWE, y::Union{ModPoly,RLWE}, oper::BFVoperator) = begin
    res = deepcopy(x)
    add_to!(res, x, y, oper)
    res
end

#PAdd
function add_to!(res::RLWE, x::RLWE, y::ModPoly, oper::BFVoperator)
    evalQ, evalT, ecder, buffQ, buffT = oper.operQ.eval, oper.evalT, oper.ecder, oper.buffQ[1], oper.buffT

    if y.len == buffT.len
        copy!(buffT, y)
        buffT.isntt[] && intt!(buffT, evalT)
        buffT.isMform[] && iMform!(buffT, evalT)
        extend_to!(buffQ, buffT, ecder)
        encode_to!(buffQ, buffQ, ecder)
    elseif y.len == buffQ.len
        encode_to!(buffQ, y, ecder)
    else
        throw(DimensionMismatch("The dimension of the input plaintext does not match the parameters."))
    end

    if x.b.isntt[]
        !buffQ.isntt[] && ntt!(buffQ, evalQ)
        !x.b.isMform[] && iMform!(buffQ, evalQ)
    elseif x.b.isMform[]
        buffQ.isntt[] && intt!(buffQ, evalQ)
    else
        buffQ.isntt[] && intt!(buffQ, evalQ)
        buffQ.isMform[] && iMform!(buffQ, evalQ)
    end

    res ≠ x && copy!(res.a, x.a)
    add_to!(res.b, x.b, buffQ, evalQ)
end

#CAdd
add_to!(res::RLWE, x::RLWE, y::RLWE, oper::BFVoperator) = add_to!(res, x, y, oper.operQ)

sub(x::RLWE, y::Union{ModPoly,RLWE}, oper::BFVoperator) = begin
    res = deepcopy(x)
    sub_to!(res, x, y, oper)
    res
end

#PSub
function sub_to!(res::RLWE, x::RLWE, y::ModPoly, oper::BFVoperator)
    evalQ, evalT, ecder, buffQ, buffT = oper.operQ.eval, oper.evalT, oper.ecder, oper.buffQ[1], oper.buffT

    if y.len == buffT.len
        copy!(buffT, y)
        buffT.isntt[] && intt!(buffT, evalT)
        buffT.isMform[] && iMform!(buffT, evalT)
        extend_to!(buffQ, buffT, ecder)
        encode_to!(buffQ, buffQ, ecder)
    elseif y.len == buffQ.len
        encode_to!(buffQ, y, ecder)
    else
        throw(DimensionMismatch("The dimension of the input plaintext does not match the parameters."))
    end

    if x.b.isntt[]
        !buffQ.isntt[] && ntt!(buffQ, evalQ)
        !x.b.isMform[] && iMform!(buffQ, evalQ)
    elseif x.b.isMform[]
        buffQ.isntt[] && intt!(buffQ, evalQ)
    else
        buffQ.isntt[] && intt!(buffQ, evalQ)
        buffQ.isMform[] && iMform!(buffQ, evalQ)
    end

    res ≠ x && copy!(res.a, x.a)
    sub_to!(res.b, x.b, buffQ, evalQ)
end

#CSub
sub_to!(res::RLWE, x::RLWE, y::RLWE, oper::BFVoperator) = sub_to!(res, x, y, oper.operQ)

mul(x::ModPoly, y::RLWE, oper::BFVoperator) = begin
    res = deepcopy(y)
    mul_to!(res, x, y, oper)
    res
end

#Pmul
function mul_to!(res::RLWE, x::ModPoly, y::RLWE, oper::BFVoperator)
    evalQ, evalT, ecder, buffQ, buffT = oper.operQ.eval, oper.evalT, oper.ecder, oper.buffQ[1], oper.buffT

    if x.len == buffT.len
        copy!(buffT, x)
        buffT.isntt[] && intt!(buffT, evalT)
        buffT.isMform[] && iMform!(buffT, evalT)
        extend_to!(buffQ, buffT, ecder)
    elseif x.len == buffQ.len
        copy!(buffQ, x)
    else
        throw(DimensionMismatch("The dimension of the input plaintext does not match the parameters."))
    end

    !buffQ.isntt[] && ntt!(buffQ, evalQ)
    !buffQ.isMform[] && Mform!(buffQ, evalQ)

    copy!(res, y)
    !res.b.isMform[] && Mform!(res.b, evalQ)
    !res.b.isntt[] && ntt!(res.b, evalQ)
    !res.a.isMform[] && Mform!(res.a, evalQ)
    !res.a.isntt[] && ntt!(res.a, evalQ)

    mul_to!(res.b, res.b, buffQ, evalQ)
    mul_to!(res.a, res.a, buffQ, evalQ)
end

mul(x::RLWE, y::RLWE, rlk::RLEV, oper::BFVoperator) = begin
    res = deepcopy(x)
    mul_to!(res, x, y, rlk, oper)
    res
end

#Cmul
@views function mul_to!(res::RLWE, x::RLWE, y::RLWE, rlk::RLEV, oper::BFVoperator)
    # Lift the ciphertexts to R.
    lift_to!(oper.mulbuffPQ[1:2], x, oper)
    lift_to!(oper.mulbuffPQ[3:4], y, oper)

    # Tensoring.
    tensor_to!(oper.mulbuffPQ[5:7], oper.mulbuffPQ[1:2], oper.mulbuffPQ[3:4], oper)

    # Scale by 1/Delta.
    buffQ = oper.buffQ
    scale_by_delta_to!(buffQ[1:3], oper.mulbuffPQ[5:7], oper)
    copy!(res.b, buffQ[1])
    copy!(res.a, buffQ[2])

    # Key-switching
    buffRLWE = oper.operQ.buffRLWE
    gadgetprod_to!(buffRLWE, buffQ[3], rlk, oper.operQ)
    ntt!(res, oper.operQ)
    add_to!(res, res, buffRLWE, oper)
end

@views function lift_to!(res::AbstractVector{ModPoly}, x::RLWE, oper::BFVoperator)
    evalPQ, evalQ, beQ2P, buffQ = oper.evalPQ, oper.operQ.eval, oper.beQ2P, oper.buffQ[1]
    N, Qlen, Plen = x.b.N, length(beQ2P.Q), length(beQ2P.P)

    b, a = res[1], res[2]

    # Lift the ciphertexts to R.
    copy!(buffQ, x.b)
    copy!(b.coeffs[:, 1:Qlen], x.b.coeffs)
    buffQ.isntt[] && intt!(buffQ, evalQ)
    buffQ.isMform[] && iMform!(buffQ, evalQ)
    basis_extend!(b.coeffs[:, Qlen+1:end], buffQ.coeffs, beQ2P)
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
    basis_extend!(a.coeffs[:, Qlen+1:end], buffQ.coeffs, beQ2P)
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

@views function tensor_to!(res::AbstractVector{ModPoly}, xvec::AbstractVector{ModPoly}, yvec::AbstractVector{ModPoly}, oper::BFVoperator)
    evalPQ = oper.evalPQ

    b0, a0, b1, a1 = xvec[1], xvec[2], yvec[1], yvec[2]
    d0, d1, d2 = res[1], res[2], res[3]

    mul_to!(d0, b0, b1, evalPQ)
    mul_to!(d1, b0, a1, evalPQ)
    muladd_to!(d1, a0, b1, evalPQ)
    mul_to!(d2, a0, a1, evalPQ)
end

@views function scale_by_delta_to!(res::AbstractVector{ModPoly}, xvec::AbstractVector{ModPoly}, oper::BFVoperator)
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

    complex_scale!(d̃0.coeffs, d0.coeffs, cscaler)
    complex_scale!(d̃1.coeffs, d1.coeffs, cscaler)
    complex_scale!(d̃2.coeffs, d2.coeffs, cscaler)
end

automorphism(x::RLWE, idx::Int64, oper::BFVoperator) = automorphism(x, idx, oper.operQ)

automorphism!(x::RLWE, idx::Int64, atk::RLEV, oper::BFVoperator) = automorphism!(x, idx, atk, oper.operQ)

automorphism_to!(res::RLWE, x::RLWE, idx::Int64, atk::RLEV, oper::BFVoperator) = automorphism_to!(res, x, idx, atk, oper.operQ)