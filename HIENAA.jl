using ChaChaCiphers: ChaChaStream, ChaCha12Stream
using Primes: totient, isprime, prevprime, nextprime, factor
using Nemo: finite_field, polynomial_ring, defining_polynomial, ZZ, coeff, lift, zzModPolyRingElem, residue_ring, cyclotomic
using Base.Threads: @threads, @spawn, @sync
using NormalForms: snf
using LinearAlgebra: mul! as matmul!, inv
using DoubleFloats: Double64, ComplexDF64

const RefBool = Base.RefValue{Bool}
const RefInt = Base.RefValue{UInt64}
const RefDouble = Base.RefValue{Double64}

include("math/modular.jl")
include("math/basischange.jl")
include("math/arithmetic.jl")
include("math/sampler.jl")
include("math/polynomials.jl")

include("ring/parameters.jl")
include("ring/transformer.jl")
include("ring/automorphism.jl")
include("ring/modpoly.jl")
include("ring/packing.jl")

include("rlwe/key.jl")
include("rlwe/ciphertext.jl")
include("rlwe/decomposition.jl")
include("rlwe/parameters.jl")
include("rlwe/operator.jl")
include("rlwe/encryptor.jl")

include("bgv/ciphertext.jl")
include("bgv/parameters.jl")
include("bgv/operator.jl")
include("bgv/encryptor.jl")

include("bfv/ciphertext.jl")
include("bfv/parameters.jl")
include("bfv/operator.jl")
include("bfv/encryptor.jl")
include("bfv/randoperator.jl")

include("algorithm/common.jl")
include("algorithm/polynomial.jl")