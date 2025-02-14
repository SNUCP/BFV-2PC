"""
SKEncryptor is a struct for secret-key encryption.
"""
struct SKEncryptor
    key::RLWEkey
    keyPQ::RLWEkeyPQ
    usampler::UniformSampler
    gsampler::CDTSampler
    oper::Operator

    SKEncryptor(key::RLWEkey, σ::Real, oper::Operator) = begin
        evalP, evalQ = oper.evalP, oper.evalQ

        if ismissing(evalP)
            keyPQ = RLWEkeyPQ(key, evalQ)
            ntt!(keyPQ, evalQ)
        else
            evalPQ = vcat(evalP, evalQ)

            keyPQ = RLWEkeyPQ(key, evalPQ)
            ntt!(keyPQ, evalPQ)
        end

        new(key, keyPQ, UniformSampler(), CDTSampler(0.0, σ), oper)
    end
end

_getkey_at(len::Int64, entor::SKEncryptor; isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    evalP, evalQ = entor.oper.evalP, entor.oper.evalQ

    if isPQ
        @assert !ismissing(evalP) "The key is not defined over R_PQ."
        Plen, Qlen = length(evalP), length(evalQ)
        @assert 0 < len ≤ Plen + Qlen "The length of input RLWE ciphertext does not match the parameters."
        if auxQ == 0
            key = entor.keyPQ[1:len]
        else
            key = entor.keyPQ[1:len]
            key[len] = entor.oper.buffpoly[1][1]
            Q = Modulus(auxQ)
            @inbounds for i = 1:key.N
                key[len][i] = _Bred(entor.key.coeffs[i], Q)
            end
        end
    else
        Qlen = length(evalQ)
        @assert 0 < len ≤ Qlen "The length of input RLWE ciphertext does not match the parameters."
        if auxQ == 0
            if ismissing(evalP)
                key = entor.keyPQ[1:len]
            else
                Plen = length(evalP)
                key = entor.keyPQ[Plen+1:Plen+len]
            end
        else
            if ismissing(evalP)
                key = entor.keyPQ[1:len]
                key[len] = entor.oper.buffpoly[1][1]
                Q = Modulus(auxQ)
                @inbounds for i = 1:key.N
                    key[len][i] = _Bred(entor.key.coeffs[i], Q)
                end
            else
                Plen = length(evalP)
                key = entor.keyPQ[Plen+1:Plen+len]
                key[len] = entor.oper.buffpoly[1][1]
                Q = Modulus(auxQ)
                @inbounds for i = 1:key.N
                    key[len][i] = _Bred(entor.key.coeffs[i], Q)
                end
            end
        end
    end

    key
end

function rlwe_sample(entor::SKEncryptor, Qlen::Int64=typemax(Int64); isPQ::Bool=false, auxQ::UInt64=UInt64(0))
    oper, N = entor.oper, entor.keyPQ.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_sample_to!(res, entor)
    res
end

function rlwe_sample_to!(res::RLWE, entor::SKEncryptor)
    us, gs = entor.usampler, entor.gsampler
    N, len, isPQ, auxQ = res.a.N, length(res.a), res.isPQ[], res.auxQ[]

    key = _getkey_at(len, entor, isPQ=isPQ, auxQ=auxQ)
    eval = _geteval_at(len, entor.oper, isPQ=isPQ, auxQ=auxQ)

    b, a = res.b, res.a
    b.isntt[] = false
    a.isntt[] = true

    uniform_random_to!(us, a, eval)

    @inbounds for j = 1:N
        ei = sample(gs)
        for i = 1:len
            b.coeffs[i][j] = _Bred(ei, eval[i])
        end
    end

    ntt!(b, eval)
    mulsub_to!(b, a, key, eval)
end

public_keygen(entor::SKEncryptor) = rlwe_sample(entor, isPQ=!ismissing(entor.oper.evalP))

"""
PKEncryptor is a struct for public key encryptions.
"""
struct PKEncryptor
    pk::RLWE
    gsampler::CDTSampler
    rgsampler::RGSampler
    oper::Operator

    PKEncryptor(pk::RLWE, σ::Real, τ::Real, oper::Operator) = begin
        evalP, evalQ = oper.evalP, oper.evalQ

        if ismissing(evalP)
            len = length(evalQ)

            @assert !pk.isPQ[] && pk.a.N == evalQ.param.N && length(pk.a) == len "The public key length does not match the parameters."
        else
            len = length(evalP) + length(evalQ)

            @assert pk.isPQ[] && pk.a.N == evalP.param.N == evalQ.param.N && length(pk.a) == len "The public key length does not match the parameters."
        end

        new(pk, CDTSampler(0.0, σ), RGSampler(τ), oper)
    end
end

function rlwe_sample(entor::PKEncryptor, Qlen::Int64=typemax(Int64); isPQ::Bool=false, auxQ::UInt64=UInt64(0))
    oper, N = entor.oper, entor.oper.param.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_sample_to!(res, entor)
    res
end

function rlwe_sample_to!(res::RLWE, entor::PKEncryptor)
    pk, gs, rgs, evalP = entor.pk, entor.gsampler, entor.rgsampler, entor.oper.evalP
    N, len, isPQ, auxQ = res.a.N, length(res.a), res.isPQ[], res.auxQ[]

    if isPQ
        pk = entor.pk[1:len]
    else
        if ismissing(evalP)
            pk = entor.pk[1:len]
        else
            Plen = length(evalP)
            pk = entor.pk[Plen+1:Plen+len]
        end
    end
    eval = _geteval_at(len, entor.oper, isPQ=isPQ)

    buff = entor.oper.buffpoly[1][1:len]
    b, a = res.b, res.a

    buff.isntt[] = false
    @inbounds for j = 1:N
        rj = sample(gs)
        for i = 1:len
            buff.coeffs[i][j] = _Bred(rj, eval[i])
        end
    end

    ntt!(buff, eval)
    mul_to!(b, buff, pk.b, eval)
    mul_to!(a, buff, pk.a, eval)

    buff.isntt[] = false
    @inbounds for j = 1:N
        e1j = sample(gs)
        for i = 1:len
            buff.coeffs[i][j] = _Bred(e1j, eval[i])
        end
    end

    ntt!(buff, eval)
    add_to!(a, buff, a, eval)

    buff.isntt[] = false
    @inbounds for j = 1:N
        e0j = sample(rgs)
        for i = 1:len
            buff.coeffs[i][j] = _Bred(e0j, eval[i])
        end
    end

    ntt!(buff, eval)
    add_to!(b, buff, b, eval)

    if auxQ ≠ 0
        intt!(b, eval)
        intt!(a, eval)

        oldQ = eval.moduli
        newQ = vcat(oldQ[1:end-1], Modulus(auxQ))

        ss = SimpleScaler(oldQ, newQ)
        eval = _geteval_at(len, entor.oper, isPQ=isPQ, auxQ=auxQ)

        simple_scale!(res.b.coeffs, res.b.coeffs, ss)
        res.b.isntt[] = false

        simple_scale!(res.a.coeffs, res.a.coeffs, ss)
        res.a.isntt[] = false
    end
end

"""
Encryptor is a struct for encryption and decryption of RLWE-based ciphertexts.
"""
Encryptor = Union{SKEncryptor,PKEncryptor}

(::Type{Encryptor})(key::RLWEkey, σ::Real, oper::Operator) = SKEncryptor(key, σ, oper)
(::Type{Encryptor})(pk::RLWE, σ::Real, τ::Real, oper::Operator) = PKEncryptor(pk, σ, τ, oper)

"""
rlwe_encrypt is a function to encrypt a plaintext into RLWE ciphertext.
"""
rlwe_encrypt(m::Union{ModScalar,ModPoly}, entor::Encryptor, Qlen::Int64=typemax(Int64); isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    oper, N = entor.oper, entor.oper.param.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_encrypt_to!(res, m, entor)
    res
end

rlwe_encrypt_to!(res::RLWE, m::ModScalar, entor::Encryptor) = begin
    @assert length(m) == length(res.b) == length(res.a) "The length of input plaintext does not match the parameters."

    rlwe_sample_to!(res, entor)
    eval = _geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = deepcopy(m)

    add_to!(res.b, res.b, buff, eval)
    res
end

rlwe_encrypt_to!(res::RLWE, m::ModPoly, entor::Encryptor) = begin
    @assert length(m) == length(res.b) == length(res.a) "The length of input plaintext does not match the parameters."

    rlwe_sample_to!(res, entor)
    eval = _geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = entor.oper.buffpoly[1][1:len]

    copy!(buff, m)
    !buff.isntt[] && ntt!(buff, eval)

    add_to!(res.b, res.b, buff, entor.eval)
    res
end

rlwe_encrypt_a(m::Union{ModScalar,ModPoly}, entor::Encryptor, Qlen::Int64=typemax(Int64); isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    oper, N = entor.oper, entor.oper.param.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_encrypt_a_to!(res, m, entor)
    res
end

rlwe_encrypt_a_to!(res::RLWE, m::ModScalar, entor::Encryptor) = begin
    @assert length(m) == length(res.b) == length(res.a) "The length of input plaintext does not match the parameters."

    rlwe_sample_to!(res, entor)
    eval = _geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = deepcopy(m)

    add_to!(res.a, res.a, buff, eval)
    res
end

rlwe_encrypt_a_to!(res::RLWE, m::ModPoly, entor::Encryptor) = begin
    @assert length(m) == length(res.b) == length(res.a) "The length of input plaintext does not match the parameters."

    rlwe_sample_to!(res, entor)
    eval = _geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = entor.oper.buffpoly[1][1:len]

    copy!(buff, m)
    !buff.isntt[] && ntt!(buff, eval)

    add_to!(res.a, res.a, buff, entor.eval)
    res
end

function phase(ct::RLWE, entor::SKEncryptor)
    res = similar(ct.a)
    phase_to!(res, ct, entor)
    res
end

function phase_to!(res::ModPoly, ct::RLWE, entor::SKEncryptor)
    @assert res.N == ct.a.N && length(res) == length(ct.a) "The parameters of the output plaintext and input ciphertext should match."
    len, isPQ, auxQ = length(ct.a), ct.isPQ[], ct.auxQ[]
    buff = entor.oper.buffpoly[end][1:len]

    key = _getkey_at(len, entor, isPQ=isPQ, auxQ=auxQ)
    eval = _geteval_at(len, entor.oper, isPQ=isPQ, auxQ=auxQ)

    copy!(res, ct.b)
    copy!(buff, ct.a)

    !res.isntt[] && ntt!(res, eval)
    !buff.isntt[] && ntt!(buff, eval)

    muladd_to!(res, buff, key, eval)

    intt!(res, eval)
end

#=====================================================================================================#

rlev_encrypt(m::Union{ModScalar,ModPoly}, entor::Encryptor) = begin
    N, evalP, evalQ, glen = entor.oper.param.N, entor.oper.evalP, entor.oper.evalQ, entor.oper.decer.glen

    if ismissing(evalP)
        len = length(evalQ)
    else
        len = length(evalP) + length(evalQ)
    end

    res = RLEV(N, len, glen)
    rlev_encrypt_to!(res, m, entor)
    res
end

function rlev_encrypt_to!(res::RLEV, m::ModScalar, entor::Encryptor)
    @assert res.glen == entor.oper.decer.glen "The length of input ciphertext does not match the parameters."

    evalQ, evalP, decer = entor.oper.evalQ, entor.oper.evalP, entor.oper.decer

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)
        @assert length(m) == len "The length of input plaintext does not match the parameters."

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].b, decer.gvec[i], m, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)
        @assert length(m) == length(eval) "The length of input plaintext does not match the parameters."

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], length(eval))
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].b, decer.gvec[i], m, eval)
        end
    end
end

function rlev_encrypt_to!(res::RLEV, m::ModPoly, entor::Encryptor)
    @assert res.glen == entor.oper.decer.glen "The length of input ciphertext does not match the parameters."

    evalQ, evalP, decer = entor.oper.evalQ, entor.oper.evalP, entor.oper.decer

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)
        @assert length(m) == len "The length of input plaintext does not match the parameters."

        buff = entor.oper.buffpoly[2]
        copy!(buff, m)
        !buff.isntt[] && ntt!(buff, eval)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].b, decer.gvec[i], buff, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)
        @assert length(m) == length(eval) "The length of input plaintext does not match the parameters."

        buff = entor.oper.buffpoly[2]
        copy!(buff, m)
        !buff.isntt[] && ntt!(buff, eval)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].b, decer.gvec[i], buff, eval)
        end
    end
end

rlev_encrypt_a(m::Union{ModScalar,ModPoly}, entor::Encryptor) = begin
    N, evalP, evalQ, glen = entor.oper.param.N, entor.oper.evalP, entor.oper.evalQ, entor.oper.decer.glen

    if ismissing(evalP)
        len = length(evalQ)
    else
        len = length(evalP) + length(evalQ)
    end

    res = RLEV(N, len, glen)
    rlev_encrypt_a_to!(res, m, entor)
    res
end

function rlev_encrypt_a_to!(res::RLEV, m::ModScalar, entor::Encryptor)
    @assert res.glen == entor.oper.decer.glen "The length of input ciphertext does not match the parameters."

    evalQ, evalP, decer = entor.oper.evalQ, entor.oper.evalP, entor.oper.decer

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)
        @assert length(m) == len "The length of input plaintext does not match the parameters."

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], m, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)
        @assert length(m) == length(eval) "The length of input plaintext does not match the parameters."

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], m, eval)
        end
    end
end

function rlev_encrypt_a_to!(res::RLEV, m::ModPoly, entor::Encryptor)
    @assert res.glen == entor.oper.decer.glen "The length of input ciphertext does not match the parameters."

    buff, evalQ, evalP, decer = entor.oper.buffpoly[2], entor.oper.evalQ, entor.oper.evalP, entor.oper.decer

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)
        @assert length(m) == len "The length of input plaintext does not match the parameters."

        copy!(buff, m)
        !buff.isntt[] && ntt!(buff, eval)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], buff, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)
        @assert length(m) == length(eval) "The length of input plaintext does not match the parameters."

        copy!(buff, m)
        !buff.isntt[] && ntt!(buff, eval)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], buff, eval)
        end
    end
end

#=====================================================================================================#

function rgsw_encrypt(m::ModScalar, entor::Encryptor)
    basketb = rlev_encrypt(m, entor)
    basketa = rlev_encrypt_a(m, entor)
    RGSW(basketb, basketa)
end

function rgsw_encrypt(m::ModPoly, entor::Encryptor)
    basketb = rlev_encrypt(m, entor)
    basketa = rlev_encrypt_a(m, entor)
    RGSW(basketb, basketa)
end

#=====================================================================================================#

function relin_keygen(entor::SKEncryptor)
    evalQ, evalP, decer = entor.oper.evalQ, entor.oper.evalP, entor.oper.decer

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)

        res = RLEV(entor.oper.param.N, len, decer.glen)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], entor.keyPQ, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)

        res = RLEV(entor.oper.param.N, len, decer.glen)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], entor.keyPQ, eval)
        end
    end

    res
end

function automorphism_keygen(idx::Int64, entor::SKEncryptor)
    buff, evalQ, evalP, decer, param = entor.oper.buffpoly[2], entor.oper.evalQ, entor.oper.evalP, entor.oper.decer, entor.oper.param

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)

        res = RLEV(param.N, len, decer.glen)

        copy!(buff, entor.keyPQ)
        automorphism!(entor.keyPQ, invmod(idx, param.m), eval)
        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].b, decer.gvec[i], buff, eval)
        end
        copy!(entor.keyPQ, buff)
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)

        res = RLEV(param.N, len, decer.glen)

        copy!(buff, entor.keyPQ)
        automorphism!(entor.keyPQ, invmod(idx, param.m), eval)
        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].b, decer.gvec[i], buff, eval)
        end
        copy!(entor.keyPQ, buff)
    end

    res
end