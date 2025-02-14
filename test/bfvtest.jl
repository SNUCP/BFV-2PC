include("../HIENAA.jl")
using Test: @test

function drop_level_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 100, 100, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct = bfv_encrypt(msg, entor)

    drop_level!(ct, 1, oper)
    res = bfv_decrypt(ct, entor)

    @test all(msg .== res)
end

function rescale_test()
    m, hw, logP, logQ, σ = 127, 32, 100, 100, 3.2

    ring_param = SubringParam(m, 1)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct = bfv_encrypt(msg, entor)
    
    rescale_to!(ct, ct, oper)
    res = bfv_decrypt(ct, entor)
    
    @test all(msg .== res)
end

function Padd_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 40, 40, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg1 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    msg2 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct = bfv_encrypt(msg1, entor)

    add_to!(ct, ct, msg2, oper)
    res = bfv_decrypt(ct, entor)
    
    @test all(mod.(msg1 + msg2, entor.ptxt_modulus.Q) .== res)
end

function Cadd_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 40, 40, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg1 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    msg2 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)

    ct1 = bfv_encrypt(msg1, entor)
    ct2 = bfv_encrypt(msg2, entor)
    
    add_to!(ct1, ct1, ct2, oper)
    res = bfv_decrypt(ct1, entor)

    @test all(mod.(msg1 + msg2, entor.ptxt_modulus.Q) .== res)
end

function Psub_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 40, 40, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg1 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    msg2 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct = bfv_encrypt(msg1, entor)

    sub_to!(ct, ct, msg2, oper)
    res = bfv_decrypt(ct, entor)

    @test all(mod.(msg1 - msg2, entor.ptxt_modulus.Q) .== res)
end

function Csub_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 40, 40, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg1 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    msg2 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct1 = bfv_encrypt(msg1, entor)
    ct2 = bfv_encrypt(msg2, entor)

    sub_to!(ct1, ct1, ct2, oper)
    res = bfv_decrypt(ct1, entor)
    
    @test all(mod.(msg1 - msg2, entor.ptxt_modulus.Q) .== res)
end

function Pmul_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 80, 160, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg1 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    msg2 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct = bfv_encrypt(msg1, entor)

    mul_to!(ct, ct, msg2, oper)
    res = bfv_decrypt(ct, entor)
    
    @test all(mod.(msg1 .* msg2, entor.ptxt_modulus.Q) .== res)
end

function Cmul_test()
    m, hw, logP, logQ, σ, τ = 127, 32, 30, 70, 5, 1024

    ring_param = SubringParam(m, 1)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    pk = public_keygen(entor)
    entorpk = BFVEncryptor(pk, σ, τ, oper)

    msg1 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    msg2 = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct1 = bfv_encrypt(msg1, entorpk)
    ct2 = bfv_encrypt(msg2, entorpk)

    rlk = relin_keygen(entor)

    mul_to!(ct1, ct1, ct2, rlk, oper)
    res = bfv_decrypt(ct1, entor)
    
    @test all(mod.(msg1 .* msg2, entor.ptxt_modulus.Q) .== res)
end

function rotate_test()
    m, hw, logP, logQ, σ = 17 * 31, 32, 80, 140, 3.2

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    param = BFVParameters(sketch)

    oper = BFVOperator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = BFVEncryptor(key, σ, oper)

    msg = uniform_random(us, entor.packer.k, entor.ptxt_modulus)
    ct = bfv_encrypt(msg, entor)

    idx = (1, 1)
    rtk = rotate_keygen(idx, entor)
    rotate_to!(ct, ct, idx, rtk, oper)
    res = bfv_decrypt(ct, entor)

    @test @views all(msg[1:2:end] .== vcat(res[4:2:end], res[2])) && all(msg[2:2:end] .== vcat(res[3:2:end], res[1]))
end

println("TEST... BFV FHE SCHEME.")
drop_level_test()
rescale_test()
Padd_test()
Cadd_test()
Psub_test()
Csub_test()
Pmul_test()
Cmul_test()
rotate_test()
println("TEST COMPLETED... WITHOUT ERROR.")