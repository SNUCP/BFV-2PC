using Test: @test
include("../HIENAA.jl")

function test_intpacker_arb()
    m, p, r = 13107, 2, 4
    param = CyclotomicParam(m)

    packer = IntPacker(p^r, param)
    eval = _PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    _mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

function test_intpacker_ntt()
    m, r = 3 * 5 * 7, 2
    param = CyclotomicParam(m)

    p = find_prime(param, 1)[1]
    packer = IntPacker(p^r, param)
    eval = _PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    _mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

function test_intpacker_subring()
    m, p, r = 8191, 2, 32
    param = SubringParam(m, ord(p, m) << 1)

    packer = IntPacker(p^r, param)
    eval = _PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    _mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

function test_intpacker_nttpow2()
    m, p, r = 1 << 6, 17, 1
    param = CyclotomicParam(m)

    packer = IntPacker(p^r, param)
    eval = _PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    _mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

function test_complexpacker_pow2()
    m = 1 << 6
    packer = ComplexPackerPow2(m)

    res = Vector{Int128}(undef, packer.N)
    msg = rand(ComplexDF64, packer.N >> 1)
    out = similar(msg)
    Δ = 1 << 40

    encode_to!(res, msg, Δ, packer)
    decode_to!(out, res, Δ, packer)

    @test all(isapprox.(msg, out, atol=packer.N / Δ))
end

println("TEST... Packing Algorithms")
test_intpacker_arb()
test_intpacker_ntt()
test_intpacker_subring()
test_intpacker_nttpow2()
test_complexpacker_pow2()
println("TEST COMPLETED... WITHOUT ERROR.")
