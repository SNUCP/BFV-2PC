struct Operator
    eval::PolyEvaluator
    decer::Decomposer
    buffRLWE::RLWE
    buffpoly::Vector{ModPoly}

    function Operator(eval::PolyEvaluator, decer::Decomposer)
        N, len = eval.ntter[1].N, length(eval.moduli)

        buffRLWE = RLWE(N, len)
        buffpoly = [ModPoly(N, len) for _ = 1:decer.glen+1]

        new(eval, decer, buffRLWE, buffpoly)
    end
end

add(x::RLWE, y::RLWE, oper::Operator) =
    RLWE(add(x.b, y.b, oper.eval), add(x.a, y.a, oper.eval))

add_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator) = begin
    add_to!(res.b, x.b, y.b, oper.eval)
    add_to!(res.a, x.a, y.a, oper.eval)
end

sub(x::RLWE, y::RLWE, oper::Operator) =
    RLWE(sub(x.b, y.b, oper.eval), sub(x.a, y.a, oper.eval))

sub_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator) = begin
    sub_to!(res.b, x.b, y.b, oper.eval)
    sub_to!(res.a, x.a, y.a, oper.eval)
end

mul(x::ModPoly, y::RLWE, oper::Operator) =
    RLWE(mul(x, y.b, oper.eval), mul(x, y.a, oper.eval))

mul_to!(res::RLWE, x::ModPoly, y::RLWE, oper::Operator) = begin
    mul_to!(res.b, x, y.b, oper.eval)
    mul_to!(res.a, x, y.a, oper.eval)
end

muladd_to!(res::RLWE, x::ModPoly, y::RLWE, oper::Operator) = begin
    muladd_to!(res.b, x, y.b, oper.eval)
    muladd_to!(res.a, x, y.a, oper.eval)
end

mulsub_to!(res::RLWE, x::ModPoly, y::RLWE, oper::Operator) = begin
    mulsub_to!(res.b, x, y.b, oper.eval)
    mulsub_to!(res.a, x, y.a, oper.eval)
end

initialise!(x::RLWE, isMform::Bool=true, isntt::Bool=true) = begin
    initialise!(x.b, isMform, isntt)
    initialise!(x.a, isMform, isntt)
end

ntt(ct::RLWE, oper::Operator) =
    RLWE(ntt(ct.b, oper.eval), ntt(ct.a, oper.eval))

ntt!(ct::RLWE, oper::Operator) =
    ntt_to!(ct, ct, oper)

ntt_to!(res::RLWE, ct::RLWE, oper::Operator) = begin
    ntt_to!(res.b, ct.b, oper.eval)
    ntt_to!(res.a, ct.a, oper.eval)
end

intt(ct::RLWE, oper::Operator) =
    RLWE(intt(ct.b, oper.eval), intt(ct.a, oper.eval))

intt!(ct::RLWE, oper::Operator) =
    intt_to!(ct, ct, oper)

intt_to!(res::RLWE, ct::RLWE, oper::Operator) = begin
    intt_to!(res.b, ct.b, oper.eval)
    intt_to!(res.a, ct.a, oper.eval)
end

#============================================================================================#

add(x::RLEV, y::RLEV, oper::Operator) = begin
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

rotate!(idx::Int64, x::RLEV) = begin
    @inbounds for i = 1:x.len
        rotate!(idx, x.stack[i].a)
        rotate!(idx, x.stack[i].b)
    end
end

initialise!(x::RLEV, isMform::Bool=true, isntt::Bool=true) = begin
    @inbounds for i = 1:x.len
        initialise!(x.stack[i], isMform, isntt)
    end
end

ntt(ct::RLEV, oper::Operator) = RLEV((ntt.(ct.stack, Ref(oper.eval))))

ntt!(ct::RLEV, oper::Operator) = ntt_to!(ct, ct, oper)

ntt_to!(res::RLEV, ct::RLEV, oper::Operator) = begin
    @assert res.glen == ct.glen "The length of input and output ciphertext should match."

    @inbounds for i = 1:res.len
        ntt_to!(res.stack[i], ct.stack[i], oper.eval)
    end
end

intt(ct::RLEV, oper::Operator) = RLEV((intt.(ct.stack, Ref(oper.eval))))

intt!(ct::RLEV, oper::Operator) = intt_to!(ct, ct, oper)

intt_to!(res::RLEV, ct::RLEV, oper::Operator) = begin
    @assert res.glen == ct.glen "The length of input and output ciphertext should match."
    @inbounds for i = 1:res.len
        intt_to!(res.stack[i], ct.stack[i], oper.eval)
    end
end

function gadgetprod(x::ModPoly, gadgetct::RLEV, oper::Operator)
    res = deepcopy(gadgetct.stack[1])
    gadgetprod_to!(res, x, gadgetct, oper)
    res
end

@views function gadgetprod_to!(res::RLWE, x::ModPoly, gadgetct::RLEV, oper::Operator)
    buffpoly, a, decer, eval = oper.buffpoly, oper.buffRLWE.a, oper.decer, oper.eval
    glen = decer.glen

    @assert gadgetct.glen == glen "The length of input ciphertext should match the parameters."
    @assert res.a.N == x.N == gadgetct.stack[1].a.N "The dimension of the input and output ciphertexts should match."

    copy!(a, x)
    if a.isntt[]
        !a.isMform[] && Mform!(a, eval)
        intt!(a, eval)
    end
    a.isMform[] && iMform!(a, eval)

    decomposeto!(buffpoly[1:glen], a, decer)

    initialise!(res, true, true)
    for i = 1:glen
        ntt!(buffpoly[i], eval)
        muladd_to!(res, buffpoly[i], gadgetct.stack[i], oper)
    end
end

@views function hoisted_gadgetprod_to!(res::RLWE, xdec::AbstractVector{ModPoly}, gadgetct::RLEV, oper::Operator)
    decer, eval = oper.decer, oper.eval
    glen = decer.glen

    @assert gadgetct.glen == glen == length(xdec) "The length of input ciphertext should match the parameters."
    @assert res.a.N == xdec[1].N == gadgetct.stack[1].a.N "The dimension of the input and output ciphertexts should match."

    initialise!(res, true, true)
    for i = 1:glen
        !xdec[i].isntt[] && ntt!(xdec[i], eval)
        muladd_to!(res, xdec[i], gadgetct.stack[i], oper)
    end
end

keyswitch(x::RLWE, ksk::RLEV, oper::Operator) = begin
    res = deepcopy(x)
    keyswitch_to!(res, x, ksk, oper)
    res
end

function keyswitch_to!(res::RLWE, ct::RLWE, ksk::RLEV, oper::Operator)
    buff, eval = oper.buffpoly[end], oper.eval
    copy!(buff, ct.b)
    if !buff.isntt[]
        !buff.isMform[] && Mform!(buff, eval)
        ntt!(buff, eval)
    elseif !ct.b.isMform[]
        Mform!(buff, eval)
    end

    gadgetprod_to!(res, ct.a, ksk, oper)
    add_to!(res.b, res.b, buff, eval)
end

function hoisted_keyswitch_to!(res::RLWE, adec::AbstractVector{ModPoly}, ct::RLWE, ksk::RLEV, oper::Operator)
    buff, eval = oper.buffpoly[end], oper.eval
    copy!(buff, ct.b)
    if !buff.isntt[]
        !buff.isMform[] && Mform!(buff, eval)
        ntt!(buff, eval)
    elseif !ct.b.isMform[]
        Mform!(buff, eval)
    end

    hoisted_gadgetprod_to!(res, adec, ksk, oper)
    add_to!(res.b, res.b, buff, eval)
end

automorphism(x::RLWE, idx::Int64, atk::RLEV, oper::Operator) = begin
    res = deepcopy(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

automorphism!(x::RLWE, idx::Int64, atk::RLEV, oper::Operator) = automorphism_to!(x, x, idx, atk, oper)

function automorphism_to!(res::RLWE, ct::RLWE, idx::Int64, atk::RLEV, oper::Operator)
    eval = oper.eval

    keyswitch_to!(res, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)
end

function hoisted_automorphism_to!(res::RLWE, adec::AbstractVector{ModPoly}, ct::RLWE, idx::Int64, atk::RLEV, oper::Operator)
    eval = oper.eval

    hoisted_keyswitch_to!(res, adec, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)
end

#======================================================================================#

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

initialise!(x::RGSW, isMform::Bool=true, isntt::Bool=true) = begin
    initialise!(x.basketb, isMform, isntt)
    initialise!(x.basketa, isMform, isntt)
end

ntt(ct::RGSW, oper::Operator) =
    RGSW(ntt(ct.basketb, oper.eval), ntt(ct.basketa, oper.eval))

ntt_to!(res::RGSW, ct::RGSW, oper::Operator) = begin
    ntt_to!(res.basketb, ct.basketb, oper.eval)
    ntt_to!(res.basketa, ct.basketa, oper.eval)
end

intt(ct::RGSW, oper::Operator) =
    RGSW(intt(ct.basketb, oper.eval), intt(ct.basketa, oper.eval))

intt_to!(res::RGSW, ct::RGSW, oper::Operator) = begin
    intt_to!(res.basketb, ct.basketb, oper.eval)
    intt_to!(res.basketa, ct.basketa, oper.eval)
end

function extprod_to!(res::RLWE, ct::RLWE, gswct::RGSW, oper::Operator)
    buffRLWE = oper.buffRLWE

    gadgetprod_to!(res, ct.b, gswct.basketb, oper)
    gadgetprod_to!(buffRLWE, ct.a, gswct.basketa, oper)
    add_to!(res, res, buffRLWE, oper)
end

@views function hoisted_extprod_to!(res::RLWE, ctdec::Vector{ModPoly}, gswct::RGSW, oper::Operator)
    buffRLWE, glen = oper.buffRLWE, length(ctdec) >> 1

    hoisted_gadgetprod_to!(res, ctdec[1:glen], gswct.basketb, oper)
    hoisted_gadgetprod_to!(buffRLWE, ctdec[glen+1:end], gswct.basketa, oper)
    add_to!(res, res, buffRLWE, oper)
end