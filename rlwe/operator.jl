struct Operator
    param::RingParam
    evalP::Union{Missing,PolyEvaluator}
    evalQ::PolyEvaluator
    auxeval::_PolyEvaluatorWord
    decer::Decomposer
    buffRLWE::Vector{RLWE}
    buffTensor::Tensor{3}
    buffpoly::Vector{ModPoly}

    function Operator(param::RLWEParameters)
        ring_param, P, Q, dlen = param.ring_param, param.P, param.Q, param.dlen
        N = ring_param.N

        if ismissing(P)
            Qlen = length(Q)
            Qmoduli = Modulus.(Q)
            evalP, evalQ = missing, PolyEvaluator(ring_param, Qmoduli)
            auxeval = _PolyEvaluatorWord(ring_param, Modulus(1 << 62))
            decer = Decomposer(Qmoduli, dlen)
            buffRLWE = [RLWE(N, Qlen, isPQ=false) for _ = 1:2]
            buffTensor = Tensor(N, Qlen, isPQ=true)
            buffpoly = [ModPoly(N, Qlen) for _ = 1:decer.glen+2]
        else
            Plen, Qlen = length(P), length(Q)
            Pmoduli, Qmoduli = Modulus.(P), Modulus.(Q)
            evalP, evalQ = PolyEvaluator(ring_param, Pmoduli), PolyEvaluator(ring_param, Qmoduli)
            auxeval = _PolyEvaluatorWord(ring_param, Modulus(1 << 62))
            decer = Decomposer(Pmoduli, Qmoduli, dlen)
            buffRLWE = [RLWE(N, Plen + Qlen, isPQ=true) for _ = 1:2]
            buffTensor = Tensor(N, Plen + Qlen, isPQ=true)
            buffpoly = [ModPoly(N, Plen + Qlen) for _ = 1:decer.glen+2]
        end

        new(ring_param, evalP, evalQ, auxeval, decer, buffRLWE, buffTensor, buffpoly)
    end
end

function _geteval_at(len::Int64, oper::Operator; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
    evalP, evalQ = oper.evalP, oper.evalQ
    if isPQ
        @assert !ismissing(evalP) "Special modulus is not defined for the operator."
        Plen, Qlen = length(evalP), length(evalQ)
        @assert Qlen + Plen ≥ len > Plen
        eval = auxQ == 0 ? vcat(evalP, evalQ[1:len-Plen]) : vcat(evalP, evalQ[1:len-Plen-1], _PolyEvaluatorWord(oper.auxeval, Modulus(auxQ)))
    else
        Qlen = length(evalQ)
        @assert len ≤ Qlen
        eval = auxQ == 0 ? evalQ[1:len] : vcat(evalQ[1:len-1], _PolyEvaluatorWord(oper.auxeval, Modulus(auxQ)))
    end

    eval
end

#============================================================================================#

add(x::RLWE, y::RLWE, oper::Operator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

function add_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = _geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    add_to!(res.b, x.b, y.b, eval)
    add_to!(res.a, x.a, y.a, eval)

    res.auxQ[] = x.auxQ[]
    res.isPQ[] = x.isPQ[]
end

sub(x::RLWE, y::RLWE, oper::Operator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = _geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    sub_to!(res.b, x.b, y.b, eval)
    sub_to!(res.a, x.a, y.a, eval)

    res.auxQ[] = x.auxQ[]
    res.isPQ[] = x.isPQ[]
end

rescale_by_P(x::RLWE, oper::Operator) = begin
    res = similar(x)
    rescale_by_P!(res, x, oper)
    res
end

function rescale_by_P!(res::RLWE, x::RLWE, oper::Operator)
    @assert x.isPQ[] "The input ciphertext should be in PQ."

    len, Plen, auxQ = length(x) - length(oper.evalP), length(oper.evalP), x.auxQ[]

    # Define Moduli.
    evalP, evalQ = oper.evalP, _geteval_at(len, oper, auxQ=auxQ)
    buff = oper.buffpoly[1][1:Plen+len]
    buffP, buffQ = buff[1:Plen], buff[Plen+1:end]

    P, Q = evalP.moduli, evalQ.moduli
    ss = SimpleScaler(P, Q)

    # Compute Pinv and Qinv.
    Pinv = ModScalar(1, evalQ)
    Qinv = ModScalar(1, evalP)

    for i = eachindex(P), j = eachindex(Q)
        Pinv.vals[j] = _Bmul(Pinv.vals[j], invmod(P[i].Q, Q[j].Q), Q[j])
        Qinv.vals[i] = _Bmul(Qinv.vals[i], invmod(Q[j].Q, P[i].Q), P[i])
    end

    # Rescale b.
    copy!(buff, x.b)
    buffP.isntt[] = x.b.isntt[]
    buffQ.isntt[] = x.b.isntt[]

    mul_to!(buffP, Qinv, buffP, evalP)
    mul_to!(buffQ, Pinv, buffQ, evalQ)
    buffP.isntt[] && intt!(buffP, evalP)

    resize!(res.b, len)
    simple_scale!(res.b.coeffs, buffP.coeffs, ss)
    res.b.isntt[] = false

    buffQ.isntt[] && ntt!(res.b, evalQ)
    add_to!(res.b, res.b, buffQ, evalQ)

    # Rescale a.
    copy!(buff, x.a)
    buffP.isntt[] = x.a.isntt[]
    buffQ.isntt[] = x.a.isntt[]

    mul_to!(buffP, Qinv, buffP, evalP)
    mul_to!(buffQ, Pinv, buffQ, evalQ)
    buffP.isntt[] && intt!(buffP, evalP)

    resize!(res.a, len)
    simple_scale!(res.a.coeffs, buffP.coeffs, ss)
    res.a.isntt[] = false

    buffQ.isntt[] && ntt!(res.a, evalQ)
    add_to!(res.a, res.a, buffQ, evalQ)

    res.isPQ[] = false
    res.auxQ[] = auxQ
end

rescale_to(x::RLWE, oper::Operator, newlen::Int64, newauxQ::UInt64=UInt64(0)) = begin
    res = similar(x)
    rescale_to!(res, x, oper, newlen, newauxQ)
    res
end

function rescale_to!(res::RLWE, x::RLWE, oper::Operator, newlen::Int64, newauxQ::UInt64=UInt64(0))
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(x), x.auxQ[]

    if len == newlen && auxQ == newauxQ
        copy!(res, x)
        return
    end

    evalQ = _geteval_at(len, oper, auxQ=auxQ)
    Q = evalQ.moduli

    if newauxQ == 0
        # Define Moduli.
        # RT -> R
        evalR, evalT = evalQ[1:newlen], evalQ[newlen+1:end]
        buff = oper.buffpoly[1][1:len]
        buffR, buffT = buff[1:newlen], buff[newlen+1:end]

        R, T = Q[1:newlen], Q[newlen+1:end]
        ss = SimpleScaler(T, R)

        # Compute Rinv and Tinv.
        Rinv = ModScalar(1, evalT)
        Tinv = ModScalar(1, evalR)

        @inbounds for i = eachindex(R), j = eachindex(T)
            Rinv.vals[j] = _Bmul(Rinv.vals[j], invmod(R[i].Q, T[j].Q), T[j])
            Tinv.vals[i] = _Bmul(Tinv.vals[i], invmod(T[j].Q, R[i].Q), R[i])
        end

        # Rescale b.
        copy!(buff, x.b)
        buffR.isntt[] = x.b.isntt[]
        buffT.isntt[] = x.b.isntt[]
        resize!(res.b, newlen)

        mul_to!(buffR, Tinv, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.b.coeffs, buffT.coeffs, ss)
        res.b.isntt[] = false

        buffR.isntt[] && ntt!(res.b, evalR)
        add_to!(res.b, res.b, buffR, evalR)

        # Rescale a.
        copy!(buff, x.a)
        buffR.isntt[] = x.a.isntt[]
        buffT.isntt[] = x.a.isntt[]
        resize!(res.a, newlen)

        mul_to!(buffR, Tinv, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.a.coeffs, buffT.coeffs, ss)
        res.a.isntt[] = false

        buffR.isntt[] && ntt!(res.a, evalR)
        add_to!(res.a, res.a, buffR, evalR)
    elseif newlen > 1
        # Define Moduli.
        # RT -> RS
        evalR, evalT = evalQ[1:newlen-1], evalQ[newlen:end]
        evalRS = _geteval_at(newlen, oper, auxQ=newauxQ)
        buff = oper.buffpoly[1][1:len]
        buffR, buffT = buff[1:newlen-1], buff[newlen:end]
        buffRS = buff[1:newlen]

        R, T, RS = Q[1:newlen-1], Q[newlen:end], evalRS.moduli
        ss = SimpleScaler(T, RS)

        # Compute Rinv and TinvS.
        Rinv = ModScalar(1, evalT)
        TinvS = ModScalar(newauxQ, evalR)

        @inbounds for i = eachindex(R), j = eachindex(T)
            Rinv.vals[j] = _Bmul(Rinv.vals[j], invmod(R[i].Q, T[j].Q), T[j])
            TinvS.vals[i] = _Bmul(TinvS.vals[i], invmod(T[j].Q, R[i].Q), R[i])
        end

        # Rescale b.
        copy!(buff, x.b)
        buffR.isntt[] = x.b.isntt[]
        buffT.isntt[] = x.b.isntt[]
        buffRS.isntt[] = x.b.isntt[]
        resize!(res.b, newlen)

        mul_to!(buffR, TinvS, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.b.coeffs, buffT.coeffs, ss)
        res.b.isntt[] = false

        buffRS.isntt[] && ntt!(res.b, evalRS)
        @. buffRS.coeffs[end] = 0
        add_to!(res.b, res.b, buffRS, evalRS)

        # Rescale a.
        copy!(buff, x.a)
        buffR.isntt[] = x.a.isntt[]
        buffT.isntt[] = x.a.isntt[]
        buffRS.isntt[] = x.a.isntt[]
        resize!(res.a, newlen)

        mul_to!(buffR, TinvS, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.a.coeffs, buffT.coeffs, ss)
        res.a.isntt[] = false

        buffRS.isntt[] && ntt!(res.a, evalRS)
        @. buffRS.coeffs[end] = 0
        add_to!(res.a, res.a, buffRS, evalRS)
    else
        # Define Moduli.
        # Q -> newQ
        ss = SimpleScaler(Q, [Modulus(newauxQ)])
        buff = oper.buffpoly[1][1:len]

        copy!(buff, x.b)
        buff.isntt[] && intt!(buff, evalQ)
        simple_scale!((@view res.b.coeffs[1:newlen]), buff.coeffs, ss)
        res.b.isntt[] = false

        copy!(buff, x.a)
        buff.isntt[] && intt!(buff, evalQ)
        simple_scale!((@view res.a.coeffs[1:newlen]), buff.coeffs, ss)
        res.a.isntt[] = false

        resize!(res, newlen)
    end

    res.isPQ[] = x.isPQ[]
    res.auxQ[] = newauxQ
end

rational_rescale(x::RLWE, Δ::Real, oper::Operator) = begin
    res = similar(x)
    rational_rescale_to!(res, x, Δ, oper)
    res
end

"""
    Compute res = ⌊x (mod Q) / Δ⌉ (mod ⌊Q / Δ⌉).
"""
function rational_rescale_to!(res::RLWE, x::RLWE, Δ::Real, oper::Operator)
    @assert Δ ≥ 1 "The scaling factor should be greater than or equal to 1."
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(x), x.auxQ[]

    # Set Q.
    if auxQ == 0
        Q = oper.evalQ.moduli[1:len]
        auxQ = Q[end].Q
    else
        Q = vcat(oper.evalQ.moduli[1:len-1], Modulus(auxQ))
    end

    # Find the next Modulus.
    newlen = len
    setprecision(BigFloat, 64 * len)
    Δ = BigFloat(Δ)
    while Δ > auxQ
        newlen -= 1
        Δ /= Q[newlen].Q
    end

    newauxQ = round(UInt64, auxQ / Δ)

    if newauxQ == 1
        newlen -= 1
        newauxQ = UInt64(0)
    end

    rescale_to!(res, x, oper, newlen, newauxQ)
end

#============================================================================================#

add(x::Tensor, y::Tensor, oper::Operator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

function add_to!(res::Tensor, x::Tensor, y::Tensor, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = _geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    @inbounds for i = 1:N
        add_to!(res.val[i], x.val[i], y.val[i], eval)
    end

    res.auxQ[] = x.auxQ[]
    res.isPQ[] = x.isPQ[]
end

sub(x::Tensor, y::Tensor, oper::Operator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::Tensor, x::Tensor, y::Tensor, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = _geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    @inbounds for i = 1:N
        sub_to!(res.val[i], x.val[i], y.val[i], eval)
    end

    res.auxQ[] = x.auxQ[]
    res.isPQ[] = x.isPQ[]
end

#=============================================================================================#

add(x::RLEV, y::RLEV, oper::Operator) = begin
    res = similar(x)

    @assert x.glen == y.glen "The length of input and output ciphertext should match."
    RLEV((@. add(x.stack, y.stack, oper.eval)))
end

add_to!(res::RLEV, x::RLEV, y::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen == y.glen "The length of input and output ciphertext should match."

    @inbounds for i = 1:res.len
        add_to!(res.stack[i], x.stack[i], y.stack[i], oper.eval)
    end
end

sub(x::RLEV, y::RLEV, oper::Operator) = begin
    @assert x.glen == y.glen "The length of input and output ciphertext should match."
    RLEV((@. sub(x.stack, y.stack, oper.eval)))
end

sub_to!(res::RLEV, x::RLEV, y::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen == y.glen "The length of input and output ciphertext should match."
    @inbounds for i = 1:res.len
        sub_to!(res.stack[i], x.stack[i], y.stack[i], oper.eval)
    end
end

#================================================================================================#

add(x::RGSW, y::RGSW, oper::Operator) =
    RGSW(add(x.basketb, y.basketb, oper.eval), add(x.basketa, y.basketa, oper.eval))

add_to!(res::RGSW, x::RGSW, y::RGSW, oper::Operator) = begin
    add_to!(res.basketb, x.basketb, y.basketb, oper.eval)
    add_to!(res.basketa, x.basketa, y.basketa, oper.eval)
end

sub(x::RGSW, y::RGSW, oper::Operator) =
    RGSW(sub(x.basketb, y.basketb, oper.eval), sub(x.basketa, y.basketa, oper.eval))

sub_to!(res::RGSW, x::RGSW, y::RGSW, oper::Operator) = begin
    sub_to!(res.basketb, x.basketb, y.basketb, oper.eval)
    sub_to!(res.basketa, x.basketa, y.basketa, oper.eval)
end

#============================================================================================#

decompose(x::ModPoly, oper::Operator; auxQ::UInt64=UInt64(0)) = begin
    len, decer = length(x), oper.decer

    if ismissing(oper.evalP)
        dlen = ceil(Int64, len / decer.dlen)
        decx = [oper.buffpoly[i][1:len] for i = 1:dlen]

        decompose_to!(decx, x, oper, auxQ=auxQ)
    else
        dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.evalP)
        decx = [oper.buffpoly[i][1:len+Plen] for i = 1:dlen]

        decompose_to!(decx, x, oper, auxQ=auxQ)
    end

    decx
end

function decompose_to!(decx::Vector{ModPoly}, x::ModPoly, oper::Operator; auxQ::UInt64=UInt64(0))
    len, decer = length(x), oper.decer
    evalQ = _geteval_at(len, oper)
    buff = oper.buffpoly[end][1:len]

    if auxQ ≠ 0
        evalauxQ = _geteval_at(len, oper, auxQ=auxQ)
        ss = SimpleScaler(evalauxQ.moduli, evalQ.moduli)

        copy!(buff, x)
        buff.isntt[] && intt!(buff, evalauxQ)

        simple_scale!(buff.coeffs, buff.coeffs, ss)
    else
        copy!(buff, x)
        buff.isntt[] && intt!(buff, evalQ)
    end

    _decompose_to!(decx, buff, decer)
end

function _hoisted_gadgetprod_to!(res::RLWE, decx::Vector{ModPoly}, ct::RLEV, oper::Operator; auxQ::UInt64=UInt64(0))
    len = length(decx[1])
    !ismissing(oper.evalP) && (len -= length(oper.evalP))
    evalQ = _geteval_at(len, oper)

    if ismissing(oper.evalP)
        @assert length(res) == len "The ciphertext length should match the parameters."
        buffRLWE = oper.buffRLWE[1][1:len]

        initialise!(buffRLWE, isntt=true, isPQ=false)

        for i = eachindex(decx)
            !decx[i].isntt[] && ntt!(decx[i], evalQ)
            muladd_to!(buffRLWE.b, decx[i], ct.stack[i].b[1:len], evalQ)
            muladd_to!(buffRLWE.a, decx[i], ct.stack[i].a[1:len], evalQ)
        end

        copy!(res, buffRLWE)
    else
        Plen = length(oper.evalP)
        evalPQ = _geteval_at(len + Plen, oper, isPQ=true)
        buffRLWE = oper.buffRLWE[1][1:len+Plen]

        initialise!(buffRLWE, isntt=true, isPQ=true)

        for i = eachindex(decx)
            !decx[i].isntt[] && ntt!(decx[i], evalPQ)
            muladd_to!(buffRLWE.b, decx[i], ct.stack[i].b[1:Plen+len], evalPQ)
            muladd_to!(buffRLWE.a, decx[i], ct.stack[i].a[1:Plen+len], evalPQ)
        end

        if auxQ ≠ 0
            intt!(buffRLWE.b, evalPQ)
            intt!(buffRLWE.a, evalPQ)
        end

        rescale_by_P!(res, buffRLWE, oper)
    end

    if auxQ ≠ 0
        rational_rescale_to!(res, res, evalQ.moduli[end].Q // auxQ, oper)
        evalQ = _geteval_at(len, oper, auxQ=auxQ)
    end

    !res.b.isntt[] && ntt!(res.b, evalQ)
    !res.a.isntt[] && ntt!(res.a, evalQ)
end

# We don't check if x is in PQ or Q. 
function _gadgetprod_to!(res::RLWE, x::ModPoly, ct::RLEV, oper::Operator; auxQ::UInt64=UInt64(0))
    len, decer = length(x), oper.decer

    if ismissing(oper.evalP)
        dlen = ceil(Int64, len / decer.dlen)
        decbuff = [oper.buffpoly[i][1:len] for i = 1:dlen]
    else
        dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.evalP)
        decbuff = [oper.buffpoly[i][1:len+Plen] for i = 1:dlen]
    end

    decompose_to!(decbuff, x, oper, auxQ=auxQ)

    _hoisted_gadgetprod_to!(res, decbuff, ct, oper, auxQ=auxQ)
end

relinearise(ct::Tensor{3}, rlk::RLEV, oper::Operator) = begin
    res = RLWE(similar(ct.vals[1]), similar(ct.vals[1]))
    relinearise_to!(res, ct, rlk, oper)
    res
end

function relinearise_to!(res::RLWE, ct::Tensor{3}, rlk::RLEV, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = _geteval_at(len, oper, auxQ=auxQ)
    buff = oper.buffpoly[end][1:len]

    _gadgetprod_to!(res, ct.vals[3], rlk, oper, auxQ=auxQ)

    copy!(buff, ct.vals[1])
    !buff.isntt[] && ntt!(buff, evalQ)
    add_to!(res.b, res.b, buff, evalQ)

    copy!(buff, ct.vals[2])
    !buff.isntt[] && ntt!(buff, evalQ)
    add_to!(res.a, res.a, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]
end

keyswitch(ct::RLWE, ksk::RLEV, oper::Operator) = begin
    res = similar(ct)
    keyswitch_to!(res, ct, ksk, oper)
    res
end

function keyswitch_to!(res::RLWE, ct::RLWE, ksk::RLEV, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = _geteval_at(len, oper, auxQ=auxQ)
    buff = oper.buffpoly[end-1][1:len]

    copy!(buff, ct.b)
    !buff.isntt[] && ntt!(buff, evalQ)

    _gadgetprod_to!(res, ct.a, ksk, oper, auxQ=auxQ)

    add_to!(res.b, res.b, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]
end

hoisted_keyswitch(adec::Vector{ModPoly}, ct::RLWE, ksk::RLEV, oper::Operator) = begin
    res = similar(ct)
    hoisted_keyswitch_to!(res, adec, ct, ksk, oper)
    res
end

function hoisted_keyswitch_to!(res::RLWE, adec::Vector{ModPoly}, ct::RLWE, ksk::RLEV, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = _geteval_at(len, oper, auxQ=auxQ)
    buff = oper.buffpoly[end-1][1:len]

    copy!(buff, ct.b)
    !buff.isntt[] && ntt!(buff, evalQ)

    _hoisted_gadgetprod_to!(res, adec, ksk, oper, auxQ=auxQ)

    add_to!(res.b, res.b, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]
end

automorphism(x::RLWE, idx::Integer, atk::RLEV, oper::Operator) = begin
    res = deepcopy(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

automorphism!(x::RLWE, idx::Integer, atk::RLEV, oper::Operator) = automorphism_to!(x, x, idx, atk, oper)

function automorphism_to!(res::RLWE, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)
    eval = _geteval_at(length(ct), oper, isPQ=ct.isPQ[], auxQ=ct.auxQ[])

    keyswitch_to!(res, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)
end

function hoisted_automorphism_to!(res::RLWE, adec::AbstractVector{ModPoly}, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)
    eval = _geteval_at(length(ct), oper, isPQ=ct.isPQ[], auxQ=ct.auxQ[])

    hoisted_keyswitch_to!(res, adec, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)
end

#======================================================================================#

function extprod(ct::RLWE, rgsw::RGSW, oper::Operator)
    res = similar(ct)
    extprod_to!(res, ct, rgsw, oper)
    res
end

function extprod_to!(res::RLWE, ct::RLWE, rgsw::RGSW, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    buff = oper.buffpoly[end-1][1:len]
    buffRLWE = oper.buffRLWE[1][1:len]

    copy!(buff, ct.a)

    _gadgetprod_to!(res, ct.b, rgsw.basketb, oper, auxQ=auxQ)
    _gadgetprod_to!(buffRLWE, buff, rgsw.basketa, oper, auxQ=auxQ)
    add_to!(res, res, buffRLWE, oper)
end

@views function hoisted_extprod_to!(res::RLWE, ctdec::Vector{ModPoly}, rgsw::RGSW, oper::Operator)
    buffRLWE, glen = oper.buffRLWE, length(ctdec) >> 1

    _hoisted_gadgetprod_to!(res, ctdec[1:glen], rgsw.basketb, oper, auxQ=auxQ)
    _hoisted_gadgetprod_to!(buffRLWE, ctdec[glen+1:end], rgsw.basketa, oper, auxQ=auxQ)
    add_to!(res, res, buffRLWE, oper)
end