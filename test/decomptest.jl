using Test: @test
include("../HIENAA.jl")

function test_decomp()
    m, Qlen, Plen, dlen = 16, 10, 2, 3

    param = CyclotomicParam(m)
    N = param.N

    Q = Modulus.(find_prime(param, 30, Qlen))
    P = Modulus.(find_prime(param, 40, Plen))
    PQ = vcat(P, Q)

    evalQ = PolyEvaluator(param, Q)
    evalPQ = PolyEvaluator(param, PQ)

    decer = Decomposer(P, Q, dlen)
    us = UniformSampler()

    a = ModPoly(N, Qlen, isntt=false)
    uniform_random_to!(us, a.coeffs, Q)

    decers = _decompose(a, decer)
    res = ModPoly(N, Plen + Qlen, isntt=false)

    for i = eachindex(decers)
        muladd_to!(res, decer.gvec[i], decers[i], evalPQ)
    end

    @test all(to_big(a, evalQ) .== to_big(res, evalPQ) .÷ prod(P))
end

println("TEST... Decomposition Algorithms")
test_decomp()
println("TEST COMPLETED... WITHOUT ERROR.")