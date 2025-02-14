"""
BFVOperator is a struct for arithmetic operations over BFV ciphertexts.
"""
struct BFVOperator
    ptxt_modulus::Modulus
    operQ::Operator
    evalR::PolyEvaluator
    packer::Union{IntPacker,Missing}
    bfv_buff::Vector{BFV}
    tensor_buff::Tensor{3}
    Qatlevel::Vector{Int64}
    Ratlevel::Vector{Int64}

    function BFVOperator(param::BFVParameters)
        ring_param, P, Q, dlen, t, ispacking, islevelled = param.ring_param, param.P, param.Q, param.dlen, param.ptxt_modulus, param.ispacking, param.islevelled

        # define operator.
        rlwe_param = RLWEParameters(ring_param, P, Q, dlen)
        operQ = Operator(rlwe_param)

        # define packer
        packer = ispacking ? IntPacker(t, ring_param) : missing

        # Setting R.
        Rprimes = collect(find_prime(ring_param, 62, 2length(Q) + 1))
        filter!(x -> x ∉ Q, Rprimes)
        Ridx, logR, logQ = 1, 0.0, log2(operQ.evalQ.moduli)
        while true
            logR += log2(Rprimes[Ridx])
            if logR > logQ + log2(ring_param.m)
                break
            end
            Ridx += 1
        end
        R = Modulus.(Rprimes[1:Ridx])
        evalR = PolyEvaluator(ring_param, R)

        # define buff.
        Qlen, Rlen = length(Q), length(R)
        bfv_buff = BFV[BFV(ring_param.N, Qlen + Rlen, 0) for _ = 1:3]
        tensor_buff = Tensor(ring_param.N, Qlen + Rlen)

        # Compute the length of modulus chain.
        sfac = ring_param.m * big(√12) * t
        minQ = Float64(log2(2 * t * 6sfac))

        if islevelled
            Qatlevel = Vector{Int64}(undef, 0)
            Qlen = length(Q)
            auxQ = UInt64(0)
            while true
                pushfirst!(Qatlevel, Qlen)

                nowQ = 0.0
                @inbounds for i = 1:Qlen-1
                    nowQ += log2(Q[i])
                end
                nowQ += auxQ == 0 ? log2(Q[Qlen]) : log2(auxQ)
                nowQ < minQ && break

                auxQ == 0 && (auxQ = Q[Qlen])
                if sfac > auxQ
                    Qlen -= 1
                    auxQ = round(UInt64, 1 / sfac * Q[Qlen] * auxQ)
                    auxQ == 1 && (auxQ = UInt64(0))
                else
                    auxQ = round(UInt64, 1 / sfac * auxQ)
                    if auxQ == 1
                        auxQ = UInt64(0)
                        Qlen -= 1
                    end
                end
            end
        else
            maxlevel = floor(Int64, log2(operQ.evalQ.moduli) / log2(sfac))
            Qatlevel = fill(length(Q), maxlevel)
        end

        Ratlevel = similar(Qatlevel)
        @views for j = 1:length(Qatlevel)
            if islevelled
                Rlen, logR, logQ = 1, 0.0, log2(operQ.evalQ.moduli[1:Qatlevel[j]])
                while true
                    logR += log2(Rprimes[Rlen])
                    if logR > logQ + log2(ring_param.m)
                        break
                    end
                    Rlen += 1
                end
            else
                Rlen = length(R)
            end

            Ratlevel[j] = Rlen
        end

        new(Modulus(t), operQ, evalR, packer, bfv_buff, tensor_buff, Qatlevel, Ratlevel)
    end

    function BFVOperator(ptxt_modulus::Modulus, operQ::Operator, evalR::PolyEvaluator, packer::Union{IntPacker,Missing}, bfv_buff::Vector{BFV}, tensor_buff::Tensor{3}, Qatlevel::Vector{Int64}, Ratlevel::Vector{Int64})
        new(ptxt_modulus, operQ, evalR, packer, bfv_buff, tensor_buff, Qatlevel, Ratlevel)
    end
end

drop_level!(x::BFV, level::Integer, oper::BFVOperator) = drop_level_to!(x, x, level, oper)

function drop_level_to!(res::BFV, x::BFV, targetlvl::Integer, oper::BFVOperator)
    currentlvl = x.level[]

    @assert targetlvl ≥ 0 "The target level should be greater than 0."
    @assert targetlvl ≤ currentlvl "The target level should be less than or equal to the current level."

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain at the current level.
    now_Qlen, tar_Qlen = Qatlevel[currentlvl+1], Qatlevel[targetlvl+1]

    if now_Qlen == tar_Qlen
        resize!(res.val, length(x.val))
        copy!(res, x)
        return
    end

    # Sanity check
    @assert now_Qlen == length(x.val) "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer, while dropping the unnecessary moduli for faster arithmetic.
    buff = oper.operQ.buffRLWE[1][1:now_Qlen]
    copy!(buff, x.val)

    # Rational rescale to the new modulus.
    resize!(res.val, tar_Qlen)
    rescale_to!(res.val, buff, operQ, tar_Qlen)

    # Update the level of the ciphertext.
    res.level[] = targetlvl
end

rescale(x::BFV, oper::BFVOperator) = begin
    res = similar(x)
    rescale_to!(res, x, oper)
    res
end

rescale_to!(res::BFV, x::BFV, oper::BFVOperator) = drop_level_to!(res, x, x.level[] - 1, oper)

function neg(x::BFV, oper::BFVOperator)
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::BFV, x::BFV, oper::BFVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ)

    # Compute res = -x.
    resize!(res.val, tar_Qlen)
    @inbounds for i = 1:tar_Qlen
        _neg_to!(res.val.b[i], x.val.b[i], evalQ[i])
        _neg_to!(res.val.a[i], x.val.a[i], evalQ[i])
    end
    res.level[] = level
end

add(x::BFV, y::AbstractVector{UInt64}, oper::BFVOperator; ispacking=true) = begin
    res = similar(x)
    add_to!(res, x, y, oper, ispacking=ispacking)
    res
end

add(x::AbstractVector{UInt64}, y::BFV, oper::BFVOperator; ispacking=true) = add(y, x, oper, ispacking=ispacking)

function add_to!(res::BFV, x::BFV, y::AbstractVector{UInt64}, oper::BFVOperator; ispacking=true)
    operQ, packer, Qatlevel = oper.operQ, oper.packer, oper.Qatlevel

    level = x.level[]
    tar_Qlen = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ)

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
    Δ = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q)
    @inbounds for i = 1:tar_Qlen
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        _Bred_to!(buff_ntt, buff_pack, evalQ[i])
        x.val.b.isntt[] && _ntt!(buff_ntt, evalQ[i])
        _muladd_to!(res.val.b[i], Δi, buff_ntt, evalQ[i])
    end
end

add_to!(res::BFV, x::AbstractVector{UInt64}, y::BFV, oper::BFVOperator; ispacking=true) = add_to!(res, y, x, oper, ispacking=ispacking)

add(x::BFV, y::UInt64, oper::BFVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::UInt64, y::BFV, oper::BFVOperator) = add(y, x, oper)

function add_to!(res::BFV, x::BFV, y::UInt64, oper::BFVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ)

    # Compute x + m.
    resize!(res.val, tar_Qlen)
    copy!(res, x)
    Δy = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q * y)
    @inbounds for i = 1:tar_Qlen
        tmp = (Δy % evalQ[i].Q.Q) % UInt64
        _add_to!(res.val.b[i], res.val.b[i], tmp, res.val.b.isntt[], evalQ[i])
    end
end

add_to!(res::BFV, x::UInt64, y::BFV, oper::BFVOperator) = add_to!(res, y, x, oper)

add(x::BFV, y::BFV, oper::BFVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

"""
Add two BFV ciphertexts.
"""
function add_to!(res::BFV, x::BFV, y::BFV, oper::BFVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    tar_Qlen = oper.Qatlevel[targetlvl+1]
    tmpx, tmpy = oper.bfv_buff[1][1:tar_Qlen], oper.bfv_buff[2][1:tar_Qlen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    resize!(res.val, tar_Qlen)
    add_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

sub(x::BFV, y::AbstractVector{UInt64}, oper::BFVOperator; ispacking=true) = begin
    res = similar(x)
    sub_to!(res, x, y, oper, ispacking=ispacking)
    res
end

@views function sub_to!(res::BFV, x::BFV, y::AbstractVector{UInt64}, oper::BFVOperator; ispacking=true)
    t = oper.ptxt_modulus
    buff = oper.bfv_buff[3].val.b.coeffs[1][1:length(y)]
    _Bred_to!(buff, y, t)
    _neg_to!(buff, buff, t)
    add_to!(res, x, buff, oper, ispacking=ispacking)
end

sub(x::AbstractVector{UInt64}, y::BFV, oper::BFVOperator; ispacking=true) = begin
    res = similar(y)
    sub_to!(res, x, y, oper, ispacking=ispacking)
    res
end

sub_to!(res::BFV, x::AbstractVector{UInt64}, y::BFV, oper::BFVOperator; ispacking=true) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper, ispacking=ispacking)
end

sub(x::BFV, y::UInt64, oper::BFVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::BFV, x::BFV, y::UInt64, oper::BFVOperator)
    tmp = _neg(_Bred(y, oper.ptxt_modulus), oper.ptxt_modulus)
    add_to!(res, x, tmp, oper)
end

sub(x::UInt64, y::BFV, oper::BFVOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BFV, x::UInt64, y::BFV, oper::BFVOperator) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BFV, y::BFV, oper::BFVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

"""
Add sub BFV ciphertexts.
"""
function sub_to!(res::BFV, x::BFV, y::BFV, oper::BFVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    tar_Qlen = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(res.val) == tar_Qlen "The length of the ciphertext should match the length of the modulus chain."

    # Match the levels.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.bfv_buff[1][1:xlen], oper.bfv_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    sub_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

mul(x::BFV, y::AbstractVector{UInt64}, oper::BFVOperator; ispacking=true) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, ispacking=ispacking)
    res
end

mul(x::AbstractVector{UInt64}, y::BFV, oper::BFVOperator; ispacking=true) = mul(y, x, oper, ispacking=ispacking)

function mul_to!(res::BFV, x::BFV, y::AbstractVector{UInt64}, oper::BFVOperator; ispacking=true)
    operQ, packer, Qatlevel = oper.operQ, oper.packer, oper.Qatlevel

    level = x.level[]
    tar_Qlen = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ)

    # Pack the input message.
    buff_ntt = operQ.buffpoly[1].coeffs[1]
    buff_pack = operQ.buffpoly[2].coeffs[1]
    if ispacking && !ismissing(packer)
        pack_to!(buff_pack, y, packer)
    else
        @assert length(y) == length(buff_pack) "The length of the plaintext should match the ring size."
        _Bred_to!(buff_pack, y, oper.ptxt_modulus)
    end

    # Compute x * m.
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

mul_to!(res::BFV, x::AbstractVector{UInt64}, y::BFV, oper::BFVOperator; ispacking=true) = mul_to!(res, y, x, oper, ispacking=ispacking)

mul(x::BFV, y::UInt64, oper::BFVOperator) = begin
    res = similar(x)
    mul_to!(res, x, y, oper)
    res
end

mul(x::UInt64, y::BFV, oper::BFVOperator) = mul(y, x, oper)

function mul_to!(res::BFV, x::BFV, y::UInt64, oper::BFVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen = Qatlevel[level+1]
    evalQ = _geteval_at(tar_Qlen, operQ)

    # Compute x * y.
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

mul_to!(res::BFV, x::UInt64, y::BFV, oper::BFVOperator) = mul_to!(res, y, x, oper)

mul(x::BFV, y::BFV, rlk::RLEV, oper::BFVOperator) = begin
    res = similar(x)
    mul_to!(res, x, y, rlk, oper)
    res
end

function mul_to!(res::BFV, x::BFV, y::BFV, rlk::RLEV, oper::BFVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    @assert targetlvl > 0 "The input ciphertexts should be at least at level 1."

    # Drop the unnecessary levels of input ciphertexts.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.bfv_buff[1][1:xlen], oper.bfv_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    # Lift.
    Qlen, Rlen = oper.Qatlevel[targetlvl+1], oper.Ratlevel[targetlvl+1]
    liftlen = Qlen + Rlen
    liftx, lifty = oper.bfv_buff[1][1:liftlen], oper.bfv_buff[2][1:liftlen]

    _lift_to!(liftx, tmpx, oper)
    _lift_to!(lifty, tmpy, oper)

    # Tensoring.
    buffQR = oper.tensor_buff[1:liftlen]
    _tensor_to!(buffQR, liftx, lifty, oper)

    # Division by Δ.
    evalQ, evalR = oper.operQ.evalQ[1:Qlen], oper.evalR[1:Rlen]
    evalQR = vcat(evalQ, evalR)
    intt_to!(buffQR, buffQR, evalQR)

    Q, R = evalQ.moduli, evalR.moduli
    cs = ComplexScaler(vcat(Q, R), Q, oper.ptxt_modulus.Q // prod(Q))

    buffQ = buffQR[1:Qlen]
    complex_scale!(buffQ.vals[1].coeffs, buffQR.vals[1].coeffs, cs)
    complex_scale!(buffQ.vals[2].coeffs, buffQR.vals[2].coeffs, cs)
    complex_scale!(buffQ.vals[3].coeffs, buffQR.vals[3].coeffs, cs)

    # Relinearisation.
    resize!(res.val, Qlen)
    relinearise_to!(res.val, buffQ, rlk, oper.operQ)
    res.level[] = targetlvl

    # Rescale the ciphertext.
    rescale_to!(res, res, oper)
end

function _lift_to!(res::BFV, x::BFV, oper::BFVOperator)
    xlevel = x.level[]

    # define operator.
    operQ = oper.operQ
    Qlen, Rlen = oper.Qatlevel[xlevel+1], oper.Ratlevel[xlevel+1]

    # Sanity check
    @assert length(x.val) == Qlen && length(res.val) == Qlen + Rlen "The length of the ciphertext should match the length of the modulus chain."

    # define evaluators.
    evalQ = operQ.evalQ[1:Qlen]
    evalR = oper.evalR[1:Rlen]

    # define BasisExtender.
    be = BasisExtender(evalQ.moduli, evalR.moduli)

    # Define buffer.
    buffQ = oper.bfv_buff[3].val[1:Qlen]
    buffR = oper.bfv_buff[3].val[Qlen+1:Qlen+Rlen]

    # Lifting.
    if x.val.b.isntt[]
        for i = 1:Qlen
            _intt_to!(buffQ.b.coeffs[i], x.val.b.coeffs[i], evalQ[i])
            _intt_to!(buffQ.a.coeffs[i], x.val.a.coeffs[i], evalQ[i])
            copy!(res.val.b.coeffs[i], x.val.b.coeffs[i])
            copy!(res.val.a.coeffs[i], x.val.a.coeffs[i])
        end
    else
        for i = 1:Qlen
            copy!(buffQ.b.coeffs[i], x.val.b.coeffs[i])
            copy!(buffQ.a.coeffs[i], x.val.a.coeffs[i])
            _ntt_to!(res.val.b.coeffs[i], x.val.b.coeffs[i], evalQ[i])
            _ntt_to!(res.val.a.coeffs[i], x.val.a.coeffs[i], evalQ[i])
        end
    end

    basis_extend!(buffR.b.coeffs, buffQ.b.coeffs, be)
    basis_extend!(buffR.a.coeffs, buffQ.a.coeffs, be)

    # Copy the result to the output in NTT form.
    for i = 1:Rlen
        _ntt_to!(res.val.b.coeffs[Qlen+i], buffR.b.coeffs[i], evalR[i])
        _ntt_to!(res.val.a.coeffs[Qlen+i], buffR.a.coeffs[i], evalR[i])
    end

    res.level[] = xlevel
    res.val.auxQ[] = 0
    res.val.isPQ[] = false
    res.val.b.isntt[] = true
    res.val.a.isntt[] = true
end

function _tensor_to!(res::Tensor{3}, x::BFV, y::BFV, oper::BFVOperator)
    # Sanity check
    @assert x.level[] == y.level[] "The level of input ciphertexts should match."

    # define operator.
    lvl = x.level[]
    Qlen, Rlen = oper.Qatlevel[lvl+1], oper.Ratlevel[lvl+1]

    # Sanity check
    @assert length(res) == length(x.val) == length(y.val) == Qlen + Rlen "The length of the ciphertext should match the parameters."

    # Tensor.
    evalQR = vcat(oper.operQ.evalQ[1:Qlen], oper.evalR[1:Rlen])

    mul_to!(res.vals[1], x.val.b, y.val.b, evalQR)
    mul_to!(res.vals[2], x.val.a, y.val.b, evalQR)
    muladd_to!(res.vals[2], x.val.b, y.val.a, evalQR)
    mul_to!(res.vals[3], x.val.a, y.val.a, evalQR)

    res.auxQ[] = 0
    res.isPQ[] = false
end

rotate(x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BFVOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, oper)
    res
end

function rotate_to!(res::BFV, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BFVOperator) where {N}
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, oper)
end

automorphism(x::BFV, idx::Int64, rtk::RLEV, oper::BFVOperator) = begin
    res = similar(x)
    automorphism_to!(res, x, idx, rtk, oper)
    res
end

function automorphism_to!(res::BFV, x::BFV, idx::Int64, atk::RLEV, oper::BFVOperator)
    @assert length(res.val) == length(x.val) "The length of input and output ciphertexts should match."

    # Sanity check
    tar_Qlen = oper.Qatlevel[x.level[]+1]
    @assert length(x.val) == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Automorphism.
    automorphism_to!(res.val, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
end