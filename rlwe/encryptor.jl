"""
SKencryptor is a struct for secret-key encryption.
"""
struct SKencryptor
    key::RLWEkeyQ
    usampler::UniformSampler
    gsampler::CDTSampler
    eval::PolyEvaluator
    buff::ModPoly

    SKencryptor(key::RLWEkey, σ::Real, eval::PolyEvaluator) = begin
        newkey = RLWEkeyQ(key, eval.moduli)
        buff = ModPoly(newkey.N, newkey.len)

        Mform!(newkey, eval)
        ntt!(newkey, eval)

        new(newkey, UniformSampler(), CDTSampler(0.0, σ), eval, buff)
    end
end

function rlwe_sample(entor::SKencryptor)
    N, len = entor.key.N, entor.key.len
    res = RLWE(N, len)
    rlwe_sample_to!(res, entor)
    res
end

@views function rlwe_sample_to!(res::RLWE, entor::SKencryptor)
    key, us, gs, eval = entor.key, entor.usampler, entor.gsampler, entor.eval
    N, len = key.N, key.len

    b, a = res.b, res.a
    b.isMform[] = true
    b.isntt[] = false
    a.isMform[] = true
    a.isntt[] = true

    @inbounds for i = 1:len
        uniform_random_to!(us, a.coeffs[:, i], eval.moduli[i])
    end

    @inbounds for i = 1:N
        ei = sample(gs)
        for j = 1:len
            if eval.ismod[j]
                b.coeffs[i, j] = Mform(ei, eval.moduli[j])
            else
                b.coeffs[i, j] = Bred(ei, eval.moduli[j])
            end
        end
    end

    ntt!(b, eval)
    mulsub_to!(b, a, key, eval)
end

"""
PKencryptor is a struct for public key encryptions.
"""
struct PKencryptor
    pk::RLWE
    gsampler::CDTSampler
    rgsampler::RGSampler
    eval::PolyEvaluator
    buff::ModPoly

    PKencryptor(pk::RLWE, σ::Real, τ::Real, eval::PolyEvaluator) = begin
        buff = ModPoly(pk.b.N, pk.b.len)

        new(pk, CDTSampler(0.0, σ), RGSampler(τ), eval, buff)
    end
end

function rlwe_sample(entor::PKencryptor)
    N, len = entor.pk.a.N, entor.pk.a.len
    res = RLWE(N, len)
    rlwe_sample_to!(res, entor)
    res
end

@views function rlwe_sample_to!(res::RLWE, entor::PKencryptor)
    pk, gs, rgs, eval, buff = entor.pk, entor.gsampler, entor.rgsampler, entor.eval, entor.buff
    N, len = pk.b.N, pk.b.len

    b, a = res.b, res.a

    buff.isntt[] = false
    buff.isMform[] = true
    @inbounds for i = 1:N
        ri = sample(gs)
        for j = 1:len
            if eval.ismod[j]
                buff.coeffs[i, j] = Mform(ri, eval.moduli[j])
            else
                buff.coeffs[i, j] = Bred(ri, eval.moduli[j])
            end
        end
    end

    ntt!(buff, eval)
    mul_to!(b, buff, pk.b, eval)
    mul_to!(a, buff, pk.a, eval)

    buff.isntt[] = false
    buff.isMform[] = true
    @inbounds for i = 1:N
        e1i = sample(gs)
        for j = 1:len
            if eval.ismod[j]
                buff.coeffs[i, j] = Mform(e1i, eval.moduli[j])
            else
                buff.coeffs[i, j] = Bred(e1i, eval.moduli[j])
            end
        end
    end

    ntt!(buff, eval)
    add_to!(a, buff, a, eval)

    buff.isntt[] = false
    buff.isMform[] = true
    @inbounds for i = 1:N
        e0i = sample(rgs)
        for j = 1:len
            if eval.ismod[j]
                buff.coeffs[i, j] = Mform(e0i, eval.moduli[j])
            else
                buff.coeffs[i, j] = Bred(e0i, eval.moduli[j])
            end
        end
    end

    ntt!(buff, eval)
    add_to!(b, buff, b, eval)
end

"""
Encryptor is a struct for encryption and decryption of RLWE-based ciphertexts.
"""
Encryptor = Union{SKencryptor,PKencryptor}

"""
Input should be in Montgomery form.
"""
rlwe_encrypt(m::ModScalar, entor::Encryptor) = begin
    res = rlwe_sample(entor)
    add_to!(res.b, res.b, m, eval)
    res
end

rlwe_encrypt_a(m::ModScalar, entor::Encryptor) = begin
    res = rlwe_sample(entor)
    add_to!(res.a, res.a, m, eval)
    res
end

rlwe_encrypt(m::ModPoly, entor::Encryptor) = begin
    res = rlwe_sample(entor)
    buff, eval = entor.buff, entor.eval

    copy!(buff, m)
    !buff.isMform[] && Mform!(buff, eval)
    !buff.isntt[] && ntt!(buff, eval)

    add_to!(res.b, res.b, buff, entor.eval)

    res
end

rlwe_encrypt_a(m::ModPoly, entor::Encryptor) = begin
    res = rlwe_sample(entor)
    buff, eval = entor.buff, entor.eval

    copy!(buff, m)
    !buff.isMform[] && Mform!(buff, eval)
    !buff.isntt[] && ntt!(buff, eval)

    add_to!(res.a, res.a, buff, entor.eval)

    res
end

function phase(ct::RLWE, entor::SKencryptor)
    res = ModPoly(ct.a.N, ct.a.len)
    phase_to!(res, ct, entor)
    res
end

function phase_to!(res::ModPoly, ct::RLWE, entor::SKencryptor)
    @assert res.N == ct.a.N && res.len == ct.a.len "The parameters of the output plaintext and input ciphertext should match."

    key, buff, eval = entor.key, entor.buff, entor.eval

    copy!(res, ct.b)
    copy!(buff, ct.a)

    !res.isntt[] && ntt!(res, eval)
    !buff.isntt[] && ntt!(buff, eval)

    muladd_to!(res, buff, key, eval)

    intt!(res, eval)
    iMform!(res, eval)
end

#=====================================================================================================#

function rlev_encrypt(m::ModScalar, entor::Encryptor, decer::Decomposer)
    stack = [rlwe_sample(entor) for _ = 1:decer.glen]
    @inbounds for i = 1:decer.glen
        muladd_to!(stack[i].b, decer.gvec[i], m, entor.eval)
    end
    RLEV(stack)
end

function rlev_encrypt_a(m::ModScalar, entor::Encryptor, decer::Decomposer)
    stack = [rlwe_sample(entor) for _ = 1:decer.glen]
    @inbounds for i = 1:decer.glen
        muladd_to!(stack[i].a, decer.gvec[i], m, entor.eval)
    end
    RLEV(stack)
end

function rlev_encrypt(m::ModPoly, entor::Encryptor, decer::Decomposer)
    stack = [rlwe_sample(entor) for _ = 1:decer.glen]
    buff, eval = entor.buff, entor.eval

    @inbounds for i = 1:decer.glen
        copy!(buff, m)
        !buff.isMform[] && Mform!(buff, eval)
        !buff.isntt[] && ntt!(buff, eval)
        mul_to!(buff, decer.gvec[i], buff, eval)
        add_to!(stack[i].b, stack[i].b, buff, eval)
    end

    RLEV(stack)
end

function rlev_encrypt_a(m::ModPoly, entor::Encryptor, decer::Decomposer)
    stack = [rlwe_sample(entor) for _ = 1:decer.glen]
    buff, eval = entor.buff, entor.eval

    @inbounds for i = 1:decer.glen
        copy!(buff, m)
        !buff.isMform[] && Mform!(buff, eval)
        !buff.isntt[] && ntt!(buff, eval)
        mul_to!(buff, decer.gvec[i], buff, eval)
        add_to!(stack[i].a, stack[i].a, buff, eval)
    end

    RLEV(stack)
end

#=====================================================================================================#

function rgsw_encrypt(m::ModScalar, entor::Encryptor, decer::Decomposer)
    basketb = rlev_encrypt(m, entor, decer)
    basketa = rlev_encrypt_a(m, entor, decer)
    RGSW(basketb, basketa)
end

function rgsw_encrypt(m::ModPoly, entor::Encryptor, decer::Decomposer)
    basketb = rlev_encrypt(m, entor, decer)
    basketa = rlev_encrypt_a(m, entor, decer)
    RGSW(basketb, basketa)
end

#=====================================================================================================#

function relin_keygen(entor::SKencryptor, decer::Decomposer)
    stack = [rlwe_sample(entor) for _ = 1:decer.glen]
    buff, eval = entor.buff, entor.eval

    mul_to!(buff, entor.key, entor.key, eval)
    @inbounds for i = 1:decer.glen
        muladd_to!(stack[i].b, decer.gvec[i], buff, eval)
    end

    RLEV(stack)
end

function aut_keygen(idx::Int64, entor::SKencryptor, decer::Decomposer)
    buff, eval, m = entor.buff, entor.eval, entor.eval.ntter[1].m

    copy!(buff, entor.key)
    automorphism!(entor.key, invmod(idx, m), eval)
    stack = [rlwe_sample(entor) for _ = 1:decer.glen]
    copy!(entor.key, buff)

    @inbounds for i = 1:decer.glen
        muladd_to!(stack[i].b, decer.gvec[i], entor.key, eval)
    end

    RLEV(stack)
end