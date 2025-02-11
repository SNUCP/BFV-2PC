using ChaChaCiphers: ChaChaStream, ChaCha12Stream
using Primes: totient, isprime, prevprime, nextprime, factor
using Base.Threads: @threads, @spawn, @sync
using LinearAlgebra: mul! as matmul!
using MultiFloats: Float64x2

const RefBool = Base.RefValue{Bool}

include("math/modular.jl")
include("math/basischange.jl")
include("math/arithmetic.jl")
include("math/sampler.jl")

include("ring/transformer.jl")
include("ring/modpoly.jl")
include("ring/automorphism.jl")

include("rlwe/key.jl")
include("rlwe/ciphertext.jl")
include("rlwe/decomposition.jl")
include("rlwe/encryptor.jl")
include("rlwe/operator.jl")

include("bfv/encoder.jl")
include("bfv/operator.jl")
include("bfv/randoperation.jl")
include("bfv/params.jl")
include("bfv/scheme.jl")