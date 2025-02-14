"""
BGVOperator is a struct for arithmetic operations over BGV ciphertexts.
"""
struct BGVOperator
    ptxt_modulus::Modulus
    operQ::Operator
    packer::Union{IntPacker,Missing}
    bgv_buff::Vector{BGV}
    tensor_buff::Tensor{3}
    Qatlevel::Vector{Tuple{Int64,UInt64}}

    function BGVOperator(param::BGVParameters)
        ring_param, P, Q, dlen, t, ispacking = param.ring_param, param.P, param.Q, param.dlen, param.ptxt_modulus, param.ispacking

        # define operator.
        rlwe_param = RLWEParameters(ring_param, P, Q, dlen)
        operQ = Operator(rlwe_param)

        # define packer
        packer = ispacking ? IntPacker(t, ring_param) : missing

        # define buff.
        bgv_buff = BGV[BGV(ring_param.N, length(Q), 0) for _ = 1:3]
        tensor_buff = Tensor(ring_param.N, length(Q))

        # Compute the length of modulus chain, and auxQ at each level.
        sfac = ring_param.m * big(√12) * t
        minQ = Float64(log2(2 * t * 6sfac))
        Qatlevel = Vector{Tuple{Int64,UInt64}}(undef, 0)

        Qlen = length(Q)
        auxQ = UInt64(0)
        while true
            pushfirst!(Qatlevel, (Qlen, auxQ))

            nowQ = 0.0
            @inbounds for i = 1:Qlen-1
                nowQ += log2(Q[i])
            end
            nowQ += auxQ == 0 ? log2(Q[Qlen]) : log2(auxQ)
            nowQ < minQ && break

            auxQ == 0 && (auxQ = Q[Qlen])
            if sfac > auxQ
                Qlen -= 1
                auxQ = round(UInt64, 1 / (t * sfac) * Q[Qlen] * auxQ) * t + 1
                auxQ == 1 && (auxQ = UInt64(0))
            else
                auxQ = round(UInt64, 1 / (t * sfac) * auxQ) * t + 1
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end

        new(Modulus(t), operQ, packer, bgv_buff, tensor_buff, Qatlevel)
    end

    function BGVOperator(ptxt_modulus::Modulus, operQ::Operator, packer::Union{IntPacker,Missing}, bgv_buff::Vector{BGV}, tensor_buff::Tensor{3}, Qatlevel::Vector{Tuple{Int64,UInt64}})
        new(ptxt_modulus, operQ, packer, bgv_buff, tensor_buff, Qatlevel)
    end
end

#=====================================================================================#

drop_level!(x::BGV, targetlvl::Integer, oper::BGVOperator) = drop_level_to!(x, x, targetlvl, oper)

function drop_level_to!(res::BGV, x::BGV, targetlvl::Integer, oper::BGVOperator)
    currentlvl = x.level[]

    @assert targetlvl ≥ 0 "The target level should be greater than 0."
    @assert targetlvl ≤ currentlvl "The target level should be less than or equal to the current level."

    if targetlvl == currentlvl
        resize!(res.val, length(x.val))
        copy!(res, x)
        return
    end

    # define operator.
    t, operQ, Qatlevel = oper.ptxt_modulus.Q, oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert now_Qlen == length(x.val) && now_auxQ == x.val.auxQ[] "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer, while dropping the unnecessary moduli for faster arithmetic.
    buff = oper.operQ.buffRLWE[1][1:tar_Qlen]
    copy!(buff, x.val[1:tar_Qlen])

    # Convert the ciphertext into a BFV ciphertext.
    # More precisely, we compute buff = (1-Q)/t * buff = [t⁻¹]_Q * buff.
    evalQ = tar_Qlen == now_Qlen ? _geteval_at(tar_Qlen, operQ, auxQ=now_auxQ) : _geteval_at(tar_Qlen, operQ)
    @inbounds for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(buff.b[i], tinvQi, buff.b[i], evalQ[i])
        _mul_to!(buff.a[i], tinvQi, buff.a[i], evalQ[i])
    end

    # Rational rescale to the new modulus.
    resize!(res.val, tar_Qlen)
    rescale_to!(res.val, buff, operQ, tar_Qlen, tar_auxQ)

    # Convert the ciphertext into BGV again.
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    @inbounds for i = eachindex(evalQ)
        _mul_to!(res.val.b[i], t, res.val.b[i], evalQ[i])
        _mul_to!(res.val.a[i], t, res.val.a[i], evalQ[i])
    end

    # Update the level of the ciphertext.
    res.level[] = targetlvl
end

rescale(x::BGV, oper::BGVOperator) = begin
    res = similar(x)
    rescale_to!(res, x, oper)
    res
end

function rescale_to!(res::BGV, x::BGV, oper::BGVOperator)
    currentlvl = x.level[]

    @assert currentlvl > 0 "The level of ciphertext should be greater than 0."

    # define operator.
    t, operQ, Qatlevel = oper.ptxt_modulus.Q, oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[currentlvl]

    # Sanity check
    @assert now_Qlen == length(x.val) && now_auxQ == x.val.auxQ[] "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer.
    buff = oper.bgv_buff[end][1:now_Qlen]
    copy!(buff, x)

    # Convert the ciphertext into a BFV ciphertext.
    # More precisely, we compute buff = (1-Q)/t * buff = [t⁻¹]_Q * buff.
    now_evalQ = _geteval_at(now_Qlen, operQ, auxQ=now_auxQ)
    tar_evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    @inbounds for i = eachindex(now_evalQ)
        tinvQi = invmod(t, now_evalQ[i].Q.Q)
        _mul_to!(buff.val.b[i], tinvQi, buff.val.b[i], now_evalQ[i])
        _mul_to!(buff.val.a[i], tinvQi, buff.val.a[i], now_evalQ[i])
    end

    # Rational rescale to the new modulus.
    rescale_to!(res.val, buff.val, operQ, tar_Qlen, tar_auxQ)

    # Convert the ciphertext into BGV again.
    @inbounds for i = eachindex(tar_evalQ)
        _mul_to!(res.val.b[i], t, res.val.b[i], tar_evalQ[i])
        _mul_to!(res.val.a[i], t, res.val.a[i], tar_evalQ[i])
    end

    # Update the level of the ciphertext.
    res.level[] = currentlvl - 1
end

function neg(x::BGV, oper::BGVOperator)
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::BGV, x::BGV, oper::BGVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Compute res = -x.
    resize!(res.val, tar_Qlen)
    @inbounds for i = 1:tar_Qlen
        _neg_to!(res.val.b[i], x.val.b[i], evalQ[i])
        _neg_to!(res.val.a[i], x.val.a[i], evalQ[i])
    end
    res.level[] = level
end

add(x::BGV, y::AbstractVector{UInt64}, oper::BGVOperator; ispacking=true) = begin
    res = similar(x)
    add_to!(res, x, y, oper, ispacking=ispacking)
    res
end

add(x::AbstractVector{UInt64}, y::BGV, oper::BGVOperator; ispacking=true) = add(y, x, oper, ispacking=ispacking)

function add_to!(res::BGV, x::BGV, y::AbstractVector{UInt64}, oper::BGVOperator; ispacking=true)
    operQ, packer, Qatlevel = oper.operQ, oper.packer, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Pack the input message.
    buff_ntt = operQ.buffpoly[1].coeffs[1]
    buff_pack = operQ.buffpoly[2].coeffs[1]
    if ispacking && !ismissing(packer)
        pack_to!(buff_pack, y, packer)
    else
        @assert length(y) == length(buff_pack) "The length of the plaintext should match the ring size."
        @. buff_pack = y
    end

    # Compute x + m.
    resize!(res.val, tar_Qlen)
    copy!(res, x)
    @inbounds for i = 1:tar_Qlen
        _Bred_to!(buff_ntt, buff_pack, evalQ[i])
        x.val.b.isntt[] && _ntt!(buff_ntt, evalQ[i])
        _add_to!(res.val.b[i], res.val.b[i], buff_ntt, evalQ[i])
    end
end

add_to!(res::BGV, x::AbstractVector{UInt64}, y::BGV, oper::BGVOperator; ispacking=true) = add_to!(res, y, x, oper, ispacking=ispacking)

add(x::BGV, y::UInt64, oper::BGVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::UInt64, y::BGV, oper::BGVOperator) = add(y, x, oper)

function add_to!(res::BGV, x::BGV, y::UInt64, oper::BGVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Compute x + m.
    resize!(res.val, tar_Qlen)
    copy!(res, x)
    y = _Bred(y, oper.ptxt_modulus)
    @inbounds for i = 1:tar_Qlen
        tmp = _Bred(y, evalQ[i])
        _add_to!(res.val.b[i], res.val.b[i], tmp, res.val.b.isntt[], evalQ[i])
    end
end

add_to!(res::BGV, x::UInt64, y::BGV, oper::BGVOperator) = add_to!(res, y, x, oper)

add(x::BGV, y::BGV, oper::BGVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

"""
Add two BGV ciphertexts.
"""
function add_to!(res::BGV, x::BGV, y::BGV, oper::BGVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    tar_Qlen, _ = oper.Qatlevel[targetlvl+1]
    tmpx, tmpy = oper.bgv_buff[1][1:tar_Qlen], oper.bgv_buff[2][1:tar_Qlen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    resize!(res.val, tar_Qlen)
    add_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

sub(x::BGV, y::AbstractVector{UInt64}, oper::BGVOperator; ispacking=true) = begin
    res = similar(x)
    sub_to!(res, x, y, oper, ispacking=ispacking)
    res
end

@views function sub_to!(res::BGV, x::BGV, y::AbstractVector{UInt64}, oper::BGVOperator; ispacking=true)
    t = oper.ptxt_modulus
    buff = oper.bgv_buff[3].val.b.coeffs[1][1:length(y)]
    _Bred_to!(buff, y, t)
    _neg_to!(buff, buff, t)
    add_to!(res, x, buff, oper, ispacking=ispacking)
end

sub(x::AbstractVector{UInt64}, y::BGV, oper::BGVOperator; ispacking=true) = begin
    res = similar(y)
    sub_to!(res, x, y, oper, ispacking=ispacking)
    res
end

function sub_to!(res::BGV, x::AbstractVector{UInt64}, y::BGV, oper::BGVOperator; ispacking=true)
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BGV, y::UInt64, oper::BGVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BGV, x::BGV, y::UInt64, oper::BGVOperator) = begin
    tmp = _neg(_Bred(y, oper.ptxt_modulus), oper.ptxt_modulus)
    add_to!(res, x, tmp, oper)
end

sub(x::UInt64, y::BGV, oper::BGVOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BGV, x::UInt64, y::BGV, oper::BGVOperator) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BGV, y::BGV, oper::BGVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

"""
Subtract two BGV ciphertexts.
"""
function sub_to!(res::BGV, x::BGV, y::BGV, oper::BGVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.bgv_buff[1][1:xlen], oper.bgv_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    sub_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

mul(x::BGV, y::AbstractVector{UInt64}, oper::BGVOperator; ispacking=true) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, ispacking=ispacking)
    res
end

mul(x::AbstractVector{UInt64}, y::BGV, oper::BGVOperator; ispacking=true) = mul(y, x, oper, ispacking=ispacking)

function mul_to!(res::BGV, x::BGV, y::AbstractVector{UInt64}, oper::BGVOperator; ispacking=true)
    operQ, packer, Qatlevel = oper.operQ, oper.packer, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Pack the input message.
    buff_ntt = operQ.buffpoly[1].coeffs[1]
    buff_pack = operQ.buffpoly[2].coeffs[1]
    if ispacking && !ismissing(packer)
        pack_to!(buff_pack, y, packer)
    else
        @assert length(y) == length(buff_pack) "The length of the plaintext should match the ring size."
        _Bred_to!(buff_pack, y, oper.ptxt_modulus)
    end

    # Compute x + m.
    resize!(res.val, tar_Qlen)
    copy!(res, x)
    !res.val.b.isntt[] && ntt!(res.val.b, evalQ)
    !res.val.a.isntt[] && ntt!(res.val.a, evalQ)
    @inbounds for i = 1:tar_Qlen
        copy!(buff_ntt, buff_pack)
        _ntt!(buff_ntt, evalQ[i])
        _mul_to!(res.val.b[i], res.val.b[i], buff_ntt, evalQ[i])
        _mul_to!(res.val.a[i], res.val.a[i], buff_ntt, evalQ[i])
    end

    # Rescale the ciphertext.
    rescale_to!(res, res, oper)
end

mul_to!(res::BGV, x::AbstractVector{UInt64}, y::BGV, oper::BGVOperator; ispacking=true) = mul_to!(res, y, x, oper, ispacking=ispacking)

mul(x::BGV, y::UInt64, oper::BGVOperator) = begin
    res = similar(x)
    mul_to!(res, x, y, oper)
    res
end

mul(x::UInt64, y::BGV, oper::BGVOperator) = mul(y, x, oper)

function mul_to!(res::BGV, x::BGV, y::UInt64, oper::BGVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Compute x + y.
    resize!(res.val, tar_Qlen)
    copy!(res, x)
    tmp = _Bred(y, oper.ptxt_modulus)
    @inbounds for i = 1:tar_Qlen
        _mul_to!(res.val.b[i], tmp, res.val.b[i], evalQ[i])
        _mul_to!(res.val.a[i], tmp, res.val.a[i], evalQ[i])
    end

    # Rescale the ciphertext.
    rescale_to!(res, res, oper)
end

mul_to!(res::BGV, x::UInt64, y::BGV, oper::BGVOperator) = mul_to!(res, y, x, oper)

function mul(x::BGV, y::BGV, rlk::RLEV, oper::BGVOperator)
    res = similar(x)
    mul_to!(res, x, y, rlk, oper)
    res
end

function mul_to!(res::BGV, x::BGV, y::BGV, rlk::RLEV, oper::BGVOperator)
    tar_Qlen = min(length(x.val), length(y.val))

    # Tensor the input ciphertexts.
    buff = oper.tensor_buff[1:tar_Qlen]
    _tensor_to!(buff, x, y, oper)

    # Relinearisation.
    resize!(res.val, tar_Qlen)
    _relinearise_to!(res.val, buff, rlk, oper)
    res.level[] = min(x.level[], y.level[])

    # Rescale the ciphertext.
    rescale_to!(res, res, oper)
end

function _tensor_to!(res::Tensor{3}, x::BGV, y::BGV, oper::BGVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    @assert targetlvl > 0 "The input ciphertexts should be at least at level 1."

    # define operator
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the target level.
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(res) == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Drop the unnecessary levels of input ciphertexts.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.bgv_buff[1][1:xlen], oper.bgv_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    # Tensor.
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    !tmpx.val.b.isntt[] && ntt!(tmpx.val.b, evalQ)
    !tmpx.val.a.isntt[] && ntt!(tmpx.val.a, evalQ)
    !tmpy.val.b.isntt[] && ntt!(tmpy.val.b, evalQ)
    !tmpy.val.a.isntt[] && ntt!(tmpy.val.a, evalQ)

    mul_to!(res.vals[1], tmpx.val.b, tmpy.val.b, evalQ)
    mul_to!(res.vals[2], tmpx.val.a, tmpy.val.b, evalQ)
    muladd_to!(res.vals[2], tmpx.val.b, tmpy.val.a, evalQ)
    mul_to!(res.vals[3], tmpx.val.a, tmpy.val.a, evalQ)

    res.auxQ[] = tar_auxQ
end

function _relinearise_to!(res::RLWE, x::Tensor{3}, rlk::RLEV, oper::BGVOperator)
    # define operator.
    t, operQ = oper.ptxt_modulus.Q, oper.operQ

    # Get the length of the modulus chain and the auxiliary modulus.
    Qlen, auxQ = length(x), x.auxQ[]

    # Copy the ciphertext to the buffer.
    buff = oper.tensor_buff[1:Qlen]
    copy!(buff, x)

    # Convert the ciphertext into a BFV ciphertext.
    evalQ = _geteval_at(Qlen, operQ, auxQ=auxQ)
    @inbounds for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(buff.vals[1][i], tinvQi, buff.vals[1][i], evalQ[i])
        _mul_to!(buff.vals[2][i], tinvQi, buff.vals[2][i], evalQ[i])
        _mul_to!(buff.vals[3][i], tinvQi, buff.vals[3][i], evalQ[i])
    end

    # Relinearise.
    resize!(res, Qlen)
    relinearise_to!(res, buff, rlk, operQ)

    # Convert back to BGV format.
    @inbounds for i = eachindex(evalQ)
        _mul_to!(res.b.coeffs[i], t, res.b.coeffs[i], evalQ[i])
        _mul_to!(res.a.coeffs[i], t, res.a.coeffs[i], evalQ[i])
    end
end

rotate(x::BGV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BGVOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, oper)
    res
end

function rotate_to!(res::BGV, x::BGV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BGVOperator) where {N}
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, oper)
end

automorphism(x::BGV, idx::Int64, atk::RLEV, oper::BGVOperator) = begin
    res = similar(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

function automorphism_to!(res::BGV, x::BGV, idx::Int64, atk::RLEV, oper::BGVOperator)
    targetlvl = x.level[]

    # define operator.
    t, operQ, Qatlevel = oper.ptxt_modulus.Q, oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."

    # Copy the ciphertext to the buffer.
    buff = oper.bgv_buff[1][1:tar_Qlen]
    copy!(buff, x)

    # Convert the ciphertext into a BFV ciphertext.
    evalQ = _geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    @inbounds for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(buff.val.b[i], tinvQi, buff.val.b[i], evalQ[i])
        _mul_to!(buff.val.a[i], tinvQi, buff.val.a[i], evalQ[i])
    end

    # Automorphism.
    resize!(res.val, tar_Qlen)
    automorphism_to!(res.val, buff.val, idx, atk, oper.operQ)

    # Convert back to BGV format.
    @inbounds for i = eachindex(evalQ)
        _mul_to!(res.val.b.coeffs[i], t, res.val.b.coeffs[i], evalQ[i])
        _mul_to!(res.val.a.coeffs[i], t, res.val.a.coeffs[i], evalQ[i])
    end

    res.level[] = x.level[]
end