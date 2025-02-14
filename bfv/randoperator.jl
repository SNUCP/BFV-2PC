const Randomiser = Union{Vector{CDTSampler},TwinCDTSampler}

struct BFVRandOperator
    oper::BFVOperator
    bfv_buff::Vector{BFV}
    poly_buff::Vector{Vector{UInt64}}
    buff_rand::Vector{Int64}
    evalT::_PolyEvaluatorWord
    us::UniformSampler
    rgs::RGSampler
    entor::PKEncryptor

    @views function BFVRandOperator(oper::BFVOperator, entor::PKEncryptor)
        bfv_buff = BFV[similar(oper.bfv_buff[1]) for _ = 1:2]
        poly_buff = [Vector{UInt64}(undef, oper.operQ.param.N) for _ = 1:2]
        buff_rand = Vector{Int64}(undef, oper.operQ.param.N)
        evalT = _PolyEvaluatorWord(oper.operQ.auxeval, oper.ptxt_modulus)
        us = UniformSampler()
        rgs = RGSampler()

        new(oper, bfv_buff, poly_buff, buff_rand, evalT, us, rgs, entor)
    end
end

neg(x::BFV, rand::BFVRandOperator) = neg(x, rand.oper)
neg_to!(res::BFV, x::BFV, rand::BFVRandOperator) = neg_to!(res, x, rand.oper)
add(x::BFV, y::AbstractVector{UInt64}, rand::BFVRandOperator; ispacking=true) = add(x, y, rand.oper, ispacking=ispacking)
add(x::AbstractVector{UInt64}, y::BFV, rand::BFVRandOperator; ispacking=true) = add(x, y, rand.oper, ispacking=ispacking)
add_to!(res::BFV, x::BFV, y::AbstractVector{UInt64}, rand::BFVRandOperator; ispacking=true) = add_to!(res, x, y, rand.oper, ispacking=ispacking)
add_to!(res::BFV, x::AbstractVector{UInt64}, y::BFV, rand::BFVRandOperator; ispacking=true) = add_to!(res, x, y, rand.oper, ispacking=ispacking)
add(x::BFV, y::UInt64, rand::BFVRandOperator) = add(x, y, rand.oper)
add(x::UInt64, y::BFV, rand::BFVRandOperator) = add(x, y, rand.oper)
add_to!(res::BFV, x::BFV, y::UInt64, rand::BFVRandOperator) = add_to!(res, x, y, rand.oper)
add_to!(res::BFV, x::UInt64, y::BFV, rand::BFVRandOperator) = add_to!(res, x, y, rand.oper)
add(x::BFV, y::BFV, rand::BFVRandOperator) = add(x, y, rand.oper)
add_to!(res::BFV, x::BFV, y::BFV, rand::BFVRandOperator) = add_to!(res, x, y, rand.oper)

sub(x::BFV, y::AbstractVector{UInt64}, rand::BFVRandOperator; ispacking=true) = sub(x, y, rand.oper, ispacking=ispacking)
sub(x::AbstractVector{UInt64}, y::BFV, rand::BFVRandOperator; ispacking=true) = sub(x, y, rand.oper, ispacking=ispacking)
sub_to!(res::BFV, x::BFV, y::AbstractVector{UInt64}, rand::BFVRandOperator; ispacking=true) = sub_to!(res, x, y, rand.oper, ispacking=ispacking)
sub_to!(res::BFV, x::AbstractVector{UInt64}, y::BFV, rand::BFVRandOperator; ispacking=true) = sub_to!(res, x, y, rand.oper, ispacking=ispacking)
sub(x::BFV, y::UInt64, rand::BFVRandOperator) = sub(x, y, rand.oper)
sub(x::UInt64, y::BFV, rand::BFVRandOperator) = sub(x, y, rand.oper)
sub_to!(res::BFV, x::BFV, y::UInt64, rand::BFVRandOperator) = sub_to!(res, x, y, rand.oper)
sub_to!(res::BFV, x::UInt64, y::BFV, rand::BFVRandOperator) = sub_to!(res, x, y, rand.oper)
sub(x::BFV, y::BFV, rand::BFVRandOperator) = sub(x, y, rand.oper)
sub_to!(res::BFV, x::BFV, y::BFV, rand::BFVRandOperator) = sub_to!(res, x, y, rand.oper)

mul(x::BFV, y::AbstractVector{UInt64}, dgs::Randomiser, τ::Float64, rand::BFVRandOperator; ispacking=true) = begin
    res = similar(y)
    mul_to!(res, x, y, dgs, τ, rand, ispacking=ispacking)
    res
end

function mul_to!(res::BFV, x::BFV, y::AbstractVector{UInt64}, dgs::Randomiser, τ::Float64, rand::BFVRandOperator; ispacking::Bool=true)
    operQ, packer = rand.oper.operQ, rand.oper.packer

    @assert length(res.val) == length(x.val) "The input and output ciphertext length should match."
    @assert x.val.b.isntt[] && x.val.a.isntt[] "The input ciphertext should be in NTT domain."

    evalQ, Qlen, t = operQ.evalQ, length(x.val), rand.oper.ptxt_modulus

    # Pack the input message.
    buff_pack = operQ.buffpoly[1].coeffs[1]
    buff_rand = rand.buff_rand
    if ispacking && !ismissing(packer)
        pack_to!(buff_pack, y, packer)
    else
        @assert length(y) == length(buff_pack) "The length of the plaintext should match the ring size."
        _Bred_to!(buff_pack, y, t)
    end

    # Randomise the input message.
    @. buff_rand = buff_pack
    if typeof(dgs) == TwinCDTSampler
        @inbounds for i = eachindex(buff_pack)
            I = sample((t.Q - buff_pack[i]) / t.Q, dgs)
            buff_rand[i] += I * Int64(t.Q)
        end
    elseif typeof(dgs) == COSACSampler
        @inbounds for i = eachindex(buff_pack)
            I = sample((t.Q - buff_pack[i]) / t.Q, dgs)
            buff_rand[i] += I * Int64(t.Q)
        end
    else
        @inbounds for i = eachindex(buff_pack)
            I = sample(dgs[t.Q-buff_pack[i]])
            buff_rand[i] += I * Int64(t.Q)
        end
    end

    # Compute x * y
    resize!(res.val, Qlen)
    copy!(res, x)
    !res.val.b.isntt[] && ntt!(res.val.b, evalQ)
    !res.val.a.isntt[] && ntt!(res.val.a, evalQ)
    @inbounds for i = 1:Qlen
        _Bred_to!(buff_pack, buff_rand, t)
        _ntt!(buff_pack, evalQ[i])
        _mul_to!(res.val.b[i], res.val.b[i], buff_pack, evalQ[i])
        _mul_to!(res.val.a[i], res.val.a[i], buff_pack, evalQ[i])
    end

    # Add noise.
    buffQ = operQ.buffpoly[1][1:Qlen]
    buffQ.isntt[] = false

    @inbounds for j = 1:buffQ.N
        ej = sample(τ, rand.rgs)
        @simd for i = eachindex(evalQ)
            buffQ.coeffs[i][j] = _Bred(ej, evalQ[i])
        end
    end
    ntt!(buffQ, evalQ)
    add_to!(res.val.b, res.val.b, buffQ, evalQ)
end

mul(x::AbstractVector{UInt64}, y::BFV, dgs::Randomiser, τ::Float64, rand::BFVRandOperator; ispacking=true) = mul(y, x, dgs, τ, rand, ispacking=ispacking)

mul_to!(res::BFV, x::AbstractVector{UInt64}, y::BFV, dgs::Randomiser, τ::Float64, rand::BFVRandOperator; ispacking=true) = mul_to!(res, y, x, dgs, τ, rand, ispacking=ispacking)

mul(x::RLWE, y::RLWE, rlk::RLEV, dgs::Randomiser, τx::Float64, τy::Float64, rand::BFVRandOperator) = begin
    res = similar(x)
    mul_to!(res, x, y, rlk, dgs, τx, τy, rand)
    res
end

#Cmul
function mul_to!(res::BFV, x::BFV, y::BFV, rlk::RLEV, dgs::Randomiser, τx::Float64, τy::Float64, rand::BFVRandOperator)
    bfvoper, evalT, us = rand.oper, rand.evalT, rand.us
    t = evalT.Q

    # sample r1, r2 
    r1, r2 = rand.poly_buff[1], rand.poly_buff[2]
    uniform_random_to!(us, r1, t)
    uniform_random_to!(us, r2, t)

    # Compute res = (x + r1)(y + r2)
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = rand.bfv_buff[1][1:xlen], rand.bfv_buff[2][1:ylen]
    add_to!(tmpx, x, r1, bfvoper; ispacking=false)
    add_to!(tmpy, y, r2, bfvoper; ispacking=false)

    # Randomize
    buff = rand.oper.bfv_buff[end][1:length(tmpx.val)]
    rlwe_sample_to!(buff.val, rand.entor)
    buff.level[] = tmpx.level[]
    add_to!(tmpx, tmpx, buff, bfvoper)

    buff = rand.oper.bfv_buff[end][1:length(tmpy.val)]
    rlwe_sample_to!(buff.val, rand.entor)
    buff.level[] = tmpy.level[]
    add_to!(tmpy, tmpy, buff, bfvoper)

    mul_to!(res, tmpx, tmpy, rlk, bfvoper)

    # Compute res -= r1 * (y + r2) + r2 * (x + r1)
    mul_to!(tmpy, r1, tmpy, dgs, τy, rand, ispacking=false)
    sub_to!(res, res, tmpy, bfvoper)

    mul_to!(tmpx, r2, tmpx, dgs, τx, rand, ispacking=false)
    sub_to!(res, res, tmpx, bfvoper)

    # Compute res += r1 * r2
    _mul_to!(r1, r1, r2, evalT)
    add_to!(res, res, r1, bfvoper, ispacking=false)
end

#Cmul
function square_to!(res::BFV, x::BFV, rlk::RLEV, dgs::Randomiser, τx::Float64, rand::BFVRandOperator)
    bfvoper, evalT, us = rand.oper, rand.evalT, rand.us
    t = evalT.Q

    # sample r
    r = rand.poly_buff[1]
    uniform_random_to!(us, r, t)

    # Compute res = (x + r)^2
    Qlen = length(x.val)
    tmpx = rand.bfv_buff[1][1:Qlen]
    add_to!(tmpx, x, r, bfvoper; ispacking=false)

    # Randomize
    buff = rand.oper.bfv_buff[end][1:length(tmpx.val)]
    rlwe_sample_to!(buff.val, rand.entor)
    buff.level[] = tmpx.level[]
    add_to!(tmpx, tmpx, buff, bfvoper)

    mul_to!(res, tmpx, tmpx, rlk, bfvoper)

    # Compute res -= 2r * (x + r)
    r2 = rand.poly_buff[2]
    @. r2 = 2 * r
    mul_to!(tmpx, r2, tmpx, dgs, τx, rand, ispacking=false)
    sub_to!(res, res, tmpx, bfvoper)

    # Compute res += r^2
    _mul_to!(r2, r, r, evalT)
    add_to!(res, res, r2, bfvoper, ispacking=false)
end

mask_randomize(x::BFV, rand::BFVRandOperator) = begin
    res = similar(x)
    mask_randomize!(res, rand)
    res
end

mask_randomize!(x::BFV, rand::BFVRandOperator) = mask_randomize_to!(x, x, rand)

mask_randomize_to!(res::BFV, x::BFV, rand::BFVRandOperator) = begin
    buff = rand.bfv_buff[1][1:length(x.val)]
    buff.level[] = x.level[]
    rlwe_sample_to!(buff.val, rand.entor)
    add_to!(res, x, buff, rand.oper)
end

rotate(x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, rand::BFVRandOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, rand)
    res
end

function rotate_to!(res::BFV, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, rand::BFVRandOperator) where {N}
    packer = rand.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, rand)
end

automorphism(x::BFV, idx::Int64, atk::RLEV, rand::BFVRandOperator) = begin
    res = deepcopy(x)
    automorphism_to!(res, x, idx, atk, rand)
    res
end

automorphism_to!(res::BFV, x::BFV, idx::Int64, atk::RLEV, rand::BFVRandOperator) = begin
    mask_randomize_to!(res, x, rand)
    automorphism_to!(res, res, idx, atk, rand.oper)
end