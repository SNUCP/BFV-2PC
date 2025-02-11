include("OCE.jl")
using Printf, BenchmarkTools

function deterministic(x::RLWE, y::ModPoly, rlk::RLEV, atk::RLEV, oper::BFVoperator)
    ctout = deepcopy(x)

    @time begin
        mul_to!(ctout, y, x, oper)
        automorphism!(ctout, 5, atk, oper)
        mul_to!(ctout, ctout, x, rlk, oper)
    end

    ctout
end

function randomised(x::RLWE, y::ModPoly, rlk::RLEV, atk::RLEV, oper::BFVoperator, pkentor::PKencryptor)
    # Parameters
    σ = 48.373546489791295
    τ_pmul, τ_tensor, τ_modswitch = 3.0514693905e+08 / √(2π), 2.1840567572e+28 / √(2π), 1.0144668279e+08 / √(2π)

    # Generate samplers and randomizer.
    t = oper.evalT.moduli[1].Q
    dgs_pmul = [CDTSampler(i / t, σ) for i = vcat(collect(1:t-1), 0)]
    dgs_cmul = VarCenSampler(σ)
    rand = Randomizer(oper, pkentor)

    ctout = deepcopy(x)

    @time begin
        rand_mul_to!(ctout, y, x, dgs_pmul, τ_pmul, rand, oper)
        rand_automorphism!(ctout, 5, atk, rand, oper)
        rand_mul_to!(ctout, ctout, x, rlk, dgs_cmul, dgs_cmul, τ_tensor, τ_modswitch, rand, oper)
        mask_randomize!(ctout, rand, oper)
    end

    ctout
end

function deterministic_flood(x::RLWE, y::ModPoly, rlk::RLEV, atk::RLEV, oper::BFVoperator, pkentor::PKencryptor)
    ctout = deepcopy(x)

    @time begin
        mul_to!(ctout, y, x, oper)
        automorphism!(ctout, 5, atk, oper)
        mul_to!(ctout, ctout, x, rlk, oper)
        flood_to!(ctout, oper, pkentor)
    end

    ctout
end

function flood_to!(x::RLWE, oper::BFVoperator, pkentor::PKencryptor)
    # Parameters
    bits = 202

    evalQ, buff = oper.operQ.eval, oper.operQ.buffRLWE

    # mask randomize
    rlwe_sample_to!(buff, pkentor)
    add_to!(x, x, buff, oper)

    # error randomize
    rng = pkentor.rgsampler.rng
    bound = big(1) << bits
    error = ModPoly(rand(rng, -bound:bound, oper.buffQ[1].N), evalQ, false)
    ntt!(error, evalQ)
    add_to!(x.b, x.b, error, evalQ)
end

function get_meanstdmax(err::Vector{BigInt})
    mean = sum(err) / length(err)
    std = sqrt(sum([(e - mean)^2 for e = err]) / length(err))
    max = maximum(abs.(err))

    mean, std, max
end

function main()
    println("IN PREPARATION...")

    # Scheme parameters
    N = 1 << 12 # Ring Dimension
    ntterQ = Modulus.(UInt64[0x000001f0a3d88001, 0x000001f0a3d94001, 0x3e147ae147af6001])   # Modulus for Number Theoretic Transform
    Q = Modulus.(UInt64[0x000001f0a3d88001, 0x000001f0a3d94001, 0x0000000002000000])   # The actual ciphertext modulus.
    T = [Modulus(256)] # The plaintext modulus.

    # Generate the basic structs.
    evalQ = PolyEvaluator(Q, [CyclotomicTransformer(2N, Qi) for Qi = ntterQ])
    evalT = PolyEvaluator(T, [CyclotomicTransformer(2N, ntterQ[end])])
    decer = Decomposer(evalQ, 1)    # Gadget decomposition parameter
    ecder = BFVencoder(evalQ, evalT)    # Structure for BFV encoding.
    oper = BFVoperator(evalQ, evalT, decer) # Structure for BFV operations.

    # Generate the secret and evaluation keys.
    sk = ternary_ringkey(UniformSampler(), N)
    skentor = SKencryptor(sk, 3.2, evalQ)
    rlk = relin_keygen(skentor, decer)
    atk = aut_keygen(5, skentor, decer)

    # Generate the public key.
    pk = rlwe_sample(skentor)
    pkentor = PKencryptor(pk, 91.455871, 7365.320243 / √(2π), evalQ)

    # Generate the encryptions.
    msg = ModPoly([i for i = 1:N], evalT, false)
    pt = extend(msg, ecder)
    encode_to!(pt, pt, ecder)
    x = rlwe_encrypt(pt, skentor)

    #===========================================#

    # DETERMINISTIC COMPUTATION

    print("DETERMINISTIC COMPUTATION : ")

    ct_deterministic = deterministic(x, msg, rlk, atk, oper)
    res = phase(ct_deterministic, skentor)
    out = error(res, ecder)
    mean, std, max = get_meanstdmax(to_big(out, evalQ))

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))

    #===========================================#

    # RANDOMISED COMPUTATION

    print("RANDOMISED COMPUTATION : ")

    ct_randomised = randomised(x, msg, rlk, atk, oper, pkentor)
    res = phase(ct_randomised, skentor)
    out = error(res, ecder)
    mean, std, max = get_meanstdmax(to_big(out, evalQ))

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))
end

function main2()
    println("IN PREPARATION...")

    # Scheme parameters
    N = 1 << 13 # Ring Dimension
    ntterQ = Modulus.(UInt64[0x00000000f85a8001, 0x00000000f85bc001, 0x00000000f85cc001, 0x00000000f85e0001, 0x00000000f8608001, 0x00000000f8610001, 0x3e147ae147b90001])   # Modulus for Number Theoretic Transform
    Q = Modulus.(UInt64[0x00000000f85a8001, 0x00000000f85bc001, 0x00000000f85cc001, 0x00000000f85e0001, 0x00000000f8608001, 0x00000000f8610001, 0x0000000001000000])   # The actual ciphertext modulus.
    T = [Modulus(256)] # The plaintext modulus.

    # Generate the basic structs.
    evalQ = PolyEvaluator(Q, [CyclotomicTransformer(2N, Qi) for Qi = ntterQ])
    evalT = PolyEvaluator(T, [CyclotomicTransformer(2N, ntterQ[end])])
    decer = Decomposer(evalQ, 1)    # Gadget decomposition parameter
    ecder = BFVencoder(evalQ, evalT)    # Structure for BFV encoding.
    oper = BFVoperator(evalQ, evalT, decer) # Structure for BFV operations.

    # Generate the secret and evaluation keys.
    sk = ternary_ringkey(UniformSampler(), N)
    skentor = SKencryptor(sk, 3.2, evalQ)
    rlk = relin_keygen(skentor, decer)
    atk = aut_keygen(5, skentor, decer)

    # Generate the public key.
    pk = rlwe_sample(skentor)
    pkentor = PKencryptor(pk, 91.455871, 7365.320243 / √(2π), evalQ)

    # Generate the encryptions.
    msg = ModPoly([i for i = 1:N], evalT, false)
    pt = extend(msg, ecder)
    encode_to!(pt, pt, ecder)
    x = rlwe_encrypt(pt, skentor)

    #===========================================#

    # DETERMINISTIC COMPUTATION

    print("DETERMINISTIC COMPUTATION : ")

    ct_flood = deterministic_flood(x, msg, rlk, atk, oper, pkentor)
    res = phase(ct_flood, skentor)
    out = error(res, ecder)
    mean, std, max = get_meanstdmax(to_big(out, evalQ))

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))
end

main()
main2()