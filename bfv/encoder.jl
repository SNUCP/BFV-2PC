mutable struct BFVencoder
    const Δ::ModScalar
    const be::BasisExtender
    const ss::SimpleScaler

    function BFVencoder(evalQ::PolyEvaluator, evalT::PolyEvaluator)
        Q, T = evalQ.moduli, evalT.moduli
        N = evalQ.ntter[1].N

        Qbig = big(1)
        @inbounds for i = eachindex(Q)
            Qbig *= Q[i].Q
        end

        Tbig = big(1)
        @inbounds for i = eachindex(T)
            Tbig *= T[i].Q
        end

        setprecision(64 * (length(Q) + length(T) + 1))
        Δ = ModScalar(round(BigInt, Qbig / Tbig), evalQ)

        be = BasisExtender(T, Q, N)
        ss = SimpleScaler(Q, T, N)

        new(Δ, be, ss)
    end
end

function extend(msg::ModPoly, ecder::BFVencoder)
    res = ModPoly(msg.N, length(ecder.ss.Q))
    extend_to!(res, msg, ecder)
    res
end

@views function extend_to!(res::ModPoly, msg::ModPoly, ecder::BFVencoder)
    @assert !msg.isMform[] && !msg.isntt[] "The input message should not be in Montgomery or evaluation form."
    @assert res.N == msg.N "The dimension of input and output polynomials should match."

    basis_extend!(res.coeffs, msg.coeffs, ecder.be)

    res.isMform[] = false
    res.isntt[] = false
end

function encode(msg::ModPoly, ecder::BFVencoder)
    res = ModPoly(msg.N, length(ecder.ss.Q))
    encode_to!(res, msg, ecder)
    res
end

@views function encode_to!(res::ModPoly, msg::ModPoly, ecder::BFVencoder)
    @assert res.len == msg.len "The length of input and output polynomials should match."
    @assert res.N == msg.N "The dimension of input and output polynomials should match."

    Δ, be = ecder.Δ, ecder.be

    @inbounds for i = 1:res.len
        Bmul_to!(res.coeffs[:, i], Δ.vals[i], msg.coeffs[:, i], be.P[i])
    end

    res.isMform[] = msg.isMform[]
    res.isntt[] = msg.isntt[]
end

function extend(msg::ModScalar, ecder::BFVencoder)
    res = ModScalar(length(ecder.ss.Q))
    extend_to!(res, msg, ecder)
    res
end

@views function extend_to!(res::ModScalar, msg::ModScalar, ecder::BFVencoder)
    @assert !msg.isMform[] "The input message should not be in Montgomery form."

    basis_extend!(res.vals, msg.vals, ecder.be)
    res.isMform[] = false
end

function encode(msg::ModScalar, ecder::BFVencoder)
    res = ModScalar(length(ecder.ss.Q))
    encode_to!(res, msg, ecder)
    res
end

@views function encode_to!(res::ModScalar, msg::ModScalar, ecder::BFVencoder)
    @assert res.len == msg.len "The length of input and output scalars should match."

    Δ, be = ecder.Δ, ecder.be

    @inbounds for i = 1:res.len
        res.vals[i] = Bmul(msg.vals[i], Δ.vals[i], be.P[i])
    end

    res.isMform[] = msg.isMform[]
end

function decode(pt::ModPoly, ecder::BFVencoder)
    res = ModPoly(pt.N, length(ecder.ss.T))
    decode_to!(res, pt, ecder)
    res
end

@views function decode_to!(res::ModPoly, pt::ModPoly, ecder::BFVencoder)
    @assert !pt.isMform[] && !pt.isntt[] "The input plaintext should not be in Montgomery or evaluation form."
    @assert res.N == pt.N "The dimension of input and output polynomials should match."

    simple_scale!(res.coeffs, pt.coeffs, ecder.ss)

    res.isntt[] = false
    res.isMform[] = false
end

function decode(pt::ModScalar, ecder::BFVencoder)
    res = ModScalar(length(ecder.ss.T))
    decode_to!(res, pt, ecder)
    res
end

@views function decode_to!(res::ModScalar, pt::ModScalar, ecder::BFVencoder)
    @assert !pt.isMform[] "The input plaintext should not be in Montgomery form."

    simple_scale!(res.vals, pt.vals, ecder.ss)

    res.isMform[] = false
end

function error(pt::ModPoly, ecder::BFVencoder)
    res = ModPoly(pt.N, length(ecder.ss.Q))
    error_to!(res, pt, ecder)
    res
end

@views function error_to!(res::ModPoly, pt::ModPoly, ecder::BFVencoder)
    @assert !pt.isMform[] && !pt.isntt[] "The input plaintext should not be in Montgomery or evaluation form."
    @assert res.N == pt.N "The dimension of input and output polynomials should match."

    ss, be, Δ = ecder.ss, ecder.be, ecder.Δ

    Qlen, Tlen = length(ss.Q), length(ss.T)
    simple_scale!(res.coeffs[:, 1:Tlen], pt.coeffs, ss)
    basis_extend!(res.coeffs, res.coeffs[:, 1:Tlen], be)
    @inbounds for i = 1:res.N
        for j = 1:Qlen
            res.coeffs[i, j] = sub(Bmul(res.coeffs[i, j], Δ.vals[j], be.P[j]), pt.coeffs[i, j], be.P[j])
        end
    end

    res.isntt[] = false
    res.isMform[] = false
end