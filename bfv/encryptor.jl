struct BFVEncryptor
    ptxt_modulus::Modulus
    entor::Encryptor
    packer::Union{Missing,IntPacker}
    Qatlevel::Vector{Int64}

    function BFVEncryptor(key::RLWEkey, σ::Real, oper::BFVOperator)
        t = oper.ptxt_modulus
        entor = SKEncryptor(key, σ, oper.operQ)
        packer = oper.packer
        Qatlevel = oper.Qatlevel

        new(t, entor, packer, Qatlevel)
    end

    function BFVEncryptor(pk::RLWE, σ::Real, τ::Real, oper::BFVOperator)
        t = oper.ptxt_modulus
        entor = PKEncryptor(pk, σ, τ, oper.operQ)
        packer = oper.packer
        Qatlevel = oper.Qatlevel

        new(t, entor, packer, Qatlevel)
    end
end

public_keygen(entor::BFVEncryptor) = begin
    @assert typeof(entor.entor) == SKEncryptor "The Encryptor should be a secret key encryptor."
    public_keygen(entor.entor)
end

function bfv_encrypt(m::UInt64, entor::BFVEncryptor; level::Int64=typemax(Int64))
    t, entor, Qatlevel = entor.ptxt_modulus.Q, entor.entor, entor.Qatlevel

    level == typemax(Int64) && (level = length(Qatlevel) - 1)
    tar_Qlen = Qatlevel[level+1]

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    val = rlwe_sample(entor, tar_Qlen)
    evalQ = _geteval_at(tar_Qlen, entor.oper)

    # Compute (b + Δ⋅m, a) for Δ = ⌊Q/t⌉
    Δ = round(BigInt, prod(evalQ.moduli) // t)
    @inbounds for i = 1:tar_Qlen
        tmp = _mul(m, (Δ % evalQ[i].Q) % UInt64, evalQ[i])
        _add_to!(val.b.coeffs[i], val.b.coeffs[i], tmp, val.b.isntt[], evalQ[i])
    end

    return BFV(val, level)
end

function bfv_encrypt(m::Vector{UInt64}, entor::BFVEncryptor; level::Int64=typemax(Int64), ispacking::Bool=true)
    t, entor, packer, Qatlevel = entor.ptxt_modulus.Q, entor.entor, entor.packer, entor.Qatlevel

    level == typemax(Int64) && (level = length(Qatlevel) - 1)
    tar_Qlen = Qatlevel[level+1]

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    val = rlwe_sample(entor, tar_Qlen)
    evalQ = _geteval_at(tar_Qlen, entor.oper)

    # Pack the input message.
    buff_ntt = entor.oper.buffpoly[1].coeffs[1]
    buff_pack = entor.oper.buffpoly[2].coeffs[1]

    if ispacking && !ismissing(packer)
        pack_to!(buff_pack, m, packer)
    else
        @assert length(m) == length(buff_pack) "The length of the plaintext should match the ring size."
        @. buff_pack = m
    end

    # Compute (b + Δ⋅m, a) for Δ = ⌊Q/t⌉.
    Δ = round(BigInt, prod(evalQ.moduli) // t)
    @inbounds for i = 1:tar_Qlen
        _Bred_to!(buff_ntt, buff_pack, evalQ[i])
        _ntt!(buff_ntt, evalQ[i])
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        _muladd_to!(val.b.coeffs[i], Δi, buff_ntt, evalQ[i])
    end

    return BFV(val, level)
end

function bfv_decrypt(x::BFV, entor::BFVEncryptor; ispacking::Bool=true)
    @assert typeof(entor.entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    t, entor, packer, Qatlevel = entor.ptxt_modulus, entor.entor, entor.packer, entor.Qatlevel

    level = x.level[]
    now_Qlen = Qatlevel[level+1]
    nowQ = _geteval_at(now_Qlen, entor.oper).moduli

    # Obtain the phase of the ciphertext.
    phase_x = phase(x.val, entor)
    buff = entor.oper.buffpoly[1].coeffs[1:1]

    # Compute plaintext.
    ss = SimpleScaler(nowQ, [t])
    simple_scale!(buff, phase_x.coeffs, ss)
    buff = entor.oper.buffpoly[1].coeffs[1]

    # Unpack the plaintext.
    if ispacking && !ismissing(packer) "The Encryptor does not support packing."
        res = Vector{UInt64}(undef, packer.k)
        unpack_to!(res, buff, packer)
    else
        res = deepcopy(buff)
    end

    res
end

function error(x::BFV, entor::BFVEncryptor)
    t, entor, Qatlevel = entor.ptxt_modulus, entor.entor, entor.Qatlevel

    level = x.level[]
    Qlen = Qatlevel[level+1]
    evalQ = _geteval_at(Qlen, entor.oper)

    # Obtain the phase of the ciphertext.
    phase_x = phase(x.val, entor)
    buff = entor.oper.buffpoly[1].coeffs[1:1]

    # Compute plaintext mod t.
    ss = SimpleScaler(evalQ.moduli, [t])
    simple_scale!(buff, phase_x.coeffs, ss)
    buff = entor.oper.buffpoly[1].coeffs[1]

    # Compute the error.
    err = phase_x
    Δ = round(BigInt, prod(evalQ.moduli) // t.Q)
    @inbounds for i = eachindex(evalQ)
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        for j = eachindex(err.coeffs[i])
            err.coeffs[i][j] = _sub(err.coeffs[i][j], _mul(Δi, buff[j], evalQ[i]), evalQ[i])
        end
    end

    to_big(err, evalQ)
end

relin_keygen(entor::BFVEncryptor) = relin_keygen(entor.entor)

function rotate_keygen(idx::NTuple{N, Int64}, entor::BFVEncryptor) where N
    packer = entor.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_keygen(autidx, entor)
end

automorphism_keygen(idx::Int64, entor::BFVEncryptor) = automorphism_keygen(idx, entor.entor)