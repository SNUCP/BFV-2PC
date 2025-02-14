using Test: @test
include("../HIENAA.jl")

function skenc_Q()
    m, hw, logQ, σ = 3^4 * 5, 8, 62 * 4, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, 0, logQ, 1)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, auxQ = UInt64(1 << 32))
    res = phase(ct, entor)
    eval = vcat(oper.evalQ[1:length(ct.a)-1], _PolyEvaluatorWord(oper.auxeval, Modulus(ct.auxQ[])))

    @test all(abs.(to_big(res, eval)) .< 6σ)
end

function skenc_PQ()
    m, hw, logP, logQ, σ = 3^4 * 7, 8, 150, 62 * 4, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ, 1)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, isPQ=true, auxQ = UInt64(1 << 32))
    res = phase(ct, entor)
    eval = vcat(vcat(oper.evalP, oper.evalQ)[1:length(ct.a)-1], _PolyEvaluatorWord(oper.auxeval, Modulus(ct.auxQ[])))

    @test all(abs.(to_big(res, eval)) .< 6σ)
end

function pkenc_Q()
    m, hw, logQ, σ = 3 * 7 * 11, 25, 62 * 4, 3.2
    τ = σ * √m
    
    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, 0, logQ, 1)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    pk = rlwe_sample(entor)
    entorpk = PKEncryptor(pk, σ, τ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entorpk, level, auxQ = UInt64(1 << 32))
    res = phase(ct, entor)
    eval = vcat(oper.evalQ[1:length(ct.a)-1], _PolyEvaluatorWord(oper.auxeval, Modulus(ct.auxQ[])))
    @test all(abs.(to_big(res, eval)) .< 6σ * hw + 6τ)
end

function pkenc_PQ()
    m, hw, logP, logQ, σ = 3 * 7 * 17, 25, 150, 62 * 4, 3.2
    τ = σ * √m

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ, 1)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    pk = rlwe_sample(entor, isPQ=true)
    entorpk = PKEncryptor(pk, σ, τ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entorpk, level, isPQ=true, auxQ = UInt64(1 << 32))
    res = phase(ct, entor)
    eval = vcat(vcat(oper.evalP, oper.evalQ)[1:length(ct.a)-1], _PolyEvaluatorWord(oper.auxeval, Modulus(ct.auxQ[])))

    @test all(abs.(to_big(res, eval)) .< 6σ * hw + 6τ)
end

println("TEST... Encryption Algorithms")
skenc_Q()
skenc_PQ()
pkenc_Q()
pkenc_PQ()
println("TEST COMPLETED... WITHOUT ERROR.")