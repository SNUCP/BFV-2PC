using Test: @test
include("../HIENAA.jl")

function test_cyclic()
    m = (2 * 3 * 5 * 7)^2

    param = CyclicParam(m)
    Q = Modulus(find_prime(param, 15)[1])
    us = UniformSampler()

    ffter = Transformer(param, Q)

    a = uniform_random(us, ffter.N, Q)
    b = deepcopy(a)

    _ntt!(a, ffter)
    _intt!(a, ffter)

    @assert all(a .== b)
end

function test_cyclotomic()
    m = 3^2 * 5 * 7 * 11
    param = CyclotomicParam(m)

    Q = Modulus(find_prime(param, 10)[1])
    us = UniformSampler()

    ffter = Transformer(param, Q)

    a = uniform_random(us, param.N, Q)
    b = deepcopy(a)

    _ntt!(a, ffter)
    _intt!(a, ffter)

    @assert all(a .== b)
end

function test_subring()
    m, d = 174763, 73
    param = SubringParam(m, d)

    Q = Modulus(find_prime(param, 3)[1])
    us = UniformSampler()

    ffter = Transformer(param, Q)

    a = uniform_random(us, ffter.N, Q)
    b = deepcopy(a)

    _ntt!(a, ffter)
    _intt!(a, ffter)

    @assert all(a .== b)
end

println("TEST... FFT Algorithms")
test_cyclic()
test_cyclotomic()
test_subring()
println("TEST COMPLETED... WITHOUT ERROR.")