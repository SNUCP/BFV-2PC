struct BGVEncryptor
    ptxt_modulus::Modulus
    entor::Encryptor
    packer::Union{Missing,IntPacker}
    Qatlevel::Vector{Tuple{Int64,UInt64}}

    function BGVEncryptor(key::RLWEkey, σ::Real, oper::BGVOperator)
        t = oper.ptxt_modulus
        entor = SKEncryptor(key, σ, oper.operQ)
        packer = oper.packer
        Qatlevel = oper.Qatlevel

        new(t, entor, packer, Qatlevel)
    end

    function BGVEncryptor(pk::RLWE, σ::Real, τ::Real, oper::BGVOperator)
        t = oper.ptxt_modulus
        entor = PKEncryptor(pk, σ, τ, oper.operQ)
        packer = oper.packer
        Qatlevel = oper.Qatlevel

        new(t, entor, packer, Qatlevel)
    end
end

public_keygen(entor::BGVEncryptor) = begin
    @assert typeof(entor.entor) == SKEncryptor "The Encryptor should be a secret key encryptor."
    public_keygen(entor.entor)
end

function bgv_encrypt(m::UInt64, entor::BGVEncryptor; level::Int64=typemax(Int64))
    t, entor, Qatlevel = entor.ptxt_modulus.Q, entor.entor, entor.Qatlevel

    level == typemax(Int64) && (level = length(Qatlevel) - 1)
    tar_Qlen, tar_auxQ = Qatlevel[level+1]

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    val = rlwe_sample(entor, tar_Qlen, auxQ=tar_auxQ)
    evalQ = _geteval_at(tar_Qlen, entor.oper, auxQ=tar_auxQ)

    # Compute (tb + m, ta), which has the phase m + te mod Q.
    @inbounds for i = 1:tar_Qlen
        _mul_to!(val.b.coeffs[i], t, val.b.coeffs[i], evalQ[i])
        _add_to!(val.b.coeffs[i], val.b.coeffs[i], m, val.b.isntt[], evalQ[i])
        _mul_to!(val.a.coeffs[i], t, val.a.coeffs[i], evalQ[i])
    end

    return BGV(val, level)
end

function bgv_encrypt(m::Vector{UInt64}, entor::BGVEncryptor; level::Int64=typemax(Int64), ispacking::Bool=true)
    t, entor, packer, Qatlevel = entor.ptxt_modulus.Q, entor.entor, entor.packer, entor.Qatlevel

    level == typemax(Int64) && (level = length(Qatlevel) - 1)
    tar_Qlen, tar_auxQ = Qatlevel[level+1]

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    val = rlwe_sample(entor, tar_Qlen, auxQ=tar_auxQ)
    evalQ = _geteval_at(tar_Qlen, entor.oper, auxQ=tar_auxQ)

    # Pack the input message.
    buff_ntt = entor.oper.buffpoly[1].coeffs[1]
    buff_pack = entor.oper.buffpoly[2].coeffs[1]

    if ispacking && !ismissing(packer)
        pack_to!(buff_pack, m, packer)
    else
        @assert length(m) == length(buff_pack) "The length of the plaintext should match the ring size."
        @. buff_pack = m
    end

    # Compute (tb + m, ta), which has the phase m + te mod Q.
    @inbounds for i = 1:tar_Qlen
        _mul_to!(val.b.coeffs[i], t, val.b.coeffs[i], evalQ[i])
        _mul_to!(val.a.coeffs[i], t, val.a.coeffs[i], evalQ[i])
        
        _Bred_to!(buff_ntt, buff_pack, evalQ[i])
        _ntt!(buff_ntt, evalQ[i])
        _add_to!(val.b.coeffs[i], val.b.coeffs[i], buff_ntt, evalQ[i])
    end

    return BGV(val, level)
end

function bgv_decrypt(x::BGV, entor::BGVEncryptor; ispacking::Bool=true)
    @assert typeof(entor.entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    t, entor, packer, Qatlevel = entor.ptxt_modulus, entor.entor, entor.packer, entor.Qatlevel

    level = x.level[]
    now_Qlen, now_auxQ = Qatlevel[level+1]
    nowQ = _geteval_at(now_Qlen, entor.oper, auxQ=now_auxQ).moduli

    # Obtain the phase of the ciphertext.
    phase_x = phase(x.val, entor)
    buff = entor.oper.buffpoly[1].coeffs[1:1]

    # Compute plaintext mod t.
    be = BasisExtender(nowQ, [t])
    basis_extend!(buff, phase_x.coeffs, be)
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

function error(x::BGV, entor::BGVEncryptor)
    t, entor, Qatlevel = entor.ptxt_modulus, entor.entor, entor.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]
    evalQ = _geteval_at(Qlen, entor.oper, auxQ=auxQ)

    # Obtain the phase of the ciphertext.
    phase_x = phase(x.val, entor)
    buff = entor.oper.buffpoly[1].coeffs[1:1]

    # Compute plaintext mod t.
    be = BasisExtender(evalQ.moduli, [t])
    basis_extend!(buff, phase_x.coeffs, be)
    buff = entor.oper.buffpoly[1].coeffs[1]

    # Compute the error.
    err = phase_x
    @inbounds for i = eachindex(evalQ)
        Qinv = invmod(t.Q, evalQ[i].Q)
        for j = eachindex(err.coeffs[i])
            err.coeffs[i][j] = _Bred(widemul(Qinv, err.coeffs[i][j] - buff[j] + evalQ[i].Q.Q), evalQ[i])
        end
    end

    to_big(err, evalQ)
end

relin_keygen(entor::BGVEncryptor) = relin_keygen(entor.entor)

function rotate_keygen(idx::NTuple{N, Int64}, entor::BGVEncryptor) where N
    packer = entor.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_keygen(autidx, entor)
end

automorphism_keygen(idx::Int64, entor::BGVEncryptor) = automorphism_keygen(idx, entor.entor)