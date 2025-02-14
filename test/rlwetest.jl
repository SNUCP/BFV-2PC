include("../HIENAA.jl")
using Test: @test

function test_rescale()
    m, d, hw, logQ, σ = 65537, 32, 192, 62 * 4, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, 0, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level)
    rescale_to!(ct, ct, oper, level - 2)
    res = phase(ct, entor)
    eval = _geteval_at(length(ct), oper; auxQ=ct.auxQ[])

    @test all(abs.(to_big(res, eval)) .< m * hw + 6σ)
end

function test_rational_rescale()
    m, hw, logQ, σ = 17 * 127, 16, 62 * 4, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, 0, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 2
    ct = rlwe_sample(entor, level)
    rational_rescale_to!(ct, ct, big(3)^3, oper)
    res = phase(ct, entor)
    eval = _geteval_at(length(ct), oper; auxQ=ct.auxQ[])

    @test all(abs.(to_big(res, eval)) .< m * hw + 6σ)
end

function test_rescale_P()
    m, hw, logP, logQ, σ = 3^4 * 5, 16, 150, 62 * 4, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, isPQ=true, auxQ=UInt64(1 << 32))
    rescale_by_P!(ct, ct, oper)
    res = phase(ct, entor)
    eval = _geteval_at(length(ct), oper; auxQ=ct.auxQ[])

    @test all(abs.(to_big(res, eval)) .< m * hw + 6σ)
end

function test_relin()
    m, d, hw, logP, logQ, σ = 65537, 64, 32, 80, 250, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 2
    ct = rlwe_sample(entor, level, auxQ=auxQ)
    ct2 = rlwe_sample(entor, level, auxQ=auxQ)
    evalQ = _geteval_at(level, oper; auxQ=auxQ)

    vals = [mul(ct.b, ct2.b, evalQ), add(mul(ct.b, ct2.a, evalQ), mul(ct.a, ct2.b, evalQ), evalQ), mul(ct.a, ct2.a, evalQ)]
    ct3 = Tensor(Tuple(vals), auxQ=auxQ)

    rlk = relin_keygen(entor)
    relinearise_to!(ct, ct3, rlk, oper)
    res = phase(ct, entor)

    @test all(abs.(to_big(res, evalQ) .< m * hw + 6σ))
end

function test_keyswitch()
    m, d, hw, logP, logQ, σ = 174763, 19, 32, 120, 250, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    key2 = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)
    entor2 = SKEncryptor(key2, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, auxQ=auxQ)
    evalQ = _geteval_at(level, oper, auxQ=auxQ)

    ksk = rlev_encrypt(entor.keyPQ, entor2)
    ksk2 = rlev_encrypt(entor2.keyPQ, entor)

    keyswitch_to!(ct, ct, ksk, oper)
    adec = decompose(ct.a, oper, auxQ=auxQ)
    ct2 = hoisted_keyswitch(adec, ct, ksk2, oper)

    res1 = phase(ct, entor2)
    res2 = phase(ct2, entor)

    @test all(abs.(to_big(res1, evalQ)) .< m * hw + 6σ) && all(abs.(to_big(res2, evalQ)) .< m * hw + 6σ)
end

function test_automorphism()
    m, hw, logP, logQ, σ = 3 * 5 * 7 * 11, 32, 120, 250, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, auxQ=auxQ)
    evalQ = _geteval_at(level, oper, auxQ=auxQ)

    idx = 4
    atk = automorphism_keygen(idx, entor)
    automorphism_to!(ct, ct, idx, atk, oper)
    res = phase(ct, entor)
    @test all(abs.(to_big(res, evalQ) .< m * hw + 6σ))
end

function test_external()
    m, hw, logP, logQ, σ = 127 * 17, 32, 120, 250, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 1
    evalPQ = vcat(oper.evalP, oper.evalQ)
    evalQ = _geteval_at(level, oper, auxQ=auxQ)

    message = ModScalar(prod(evalQ.moduli) >> 4, evalQ)
    mu = ModScalar(2, evalPQ)

    ct = rlwe_encrypt(message, entor, level, auxQ=auxQ)
    rgsw = rgsw_encrypt(mu, entor)
    extprod_to!(ct, ct, rgsw, oper)
    res = to_big(phase(ct, entor), evalQ)
    @test @views round(Int64, res[1] / prod(evalQ.moduli) * 16) == 2 && all(round.(Int64, res[2:end] / prod(evalQ.moduli) * 16) .== 0)
end

println("TEST... RLWE Algorithms")
test_rescale()
test_rational_rescale()
test_rescale_P()
test_relin()
test_keyswitch()
test_automorphism()
test_external()
println("TEST COMPLETED... WITHOUT ERROR.")