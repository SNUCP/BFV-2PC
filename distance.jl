include("HIENAA.jl")
using Printf

function deterministic(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, oper::BFVOperator)
    ctout = similar(x)
    tmp = similar(x)

    @time begin
        # (x - y)^2
        sub_to!(ctout, x, y, oper)
        mul_to!(ctout, ctout, ctout, rlk, oper)

        # Rotsum, assuming that the rotation group is a subgroup of <5>.
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], oper)
            add_to!(ctout, ctout, tmp, oper)
        end
    end

    ctout
end

function randomised(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, oper::BFVOperator, pkentor::PKEncryptor)
    # Parameters
    σ = 48.373546489791295
    τ_cmul = 1.992990407684452e13

    # Generate samplers and randomizer.
    dgs_pmul = TwinCDTSampler(σ)
    rand = BFVRandOperator(oper, pkentor)

    ctout = similar(x)
    tmp = similar(x)

    @time begin
        # (x - y)^2
        sub_to!(ctout, x, y, rand)
        mul_to!(ctout, ctout, ctout, rlk, dgs_pmul, τ_cmul, τ_cmul, rand)

        # rotsum
        # rotation group <5>
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], rand)
            add_to!(ctout, ctout, tmp, rand)
        end

        mask_randomize!(ctout, rand)
    end

    ctout
end

function optimised(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, oper::BFVOperator, pkentor::PKEncryptor)
    # Parameters
    σ = 48.373546489791295
    τ_cmul = 1.992990407684452e13

    # Generate samplers and randomizer.
    dgs_pmul = TwinCDTSampler(σ)
    rand = BFVRandOperator(oper, pkentor)

    ctout = similar(x)
    tmp = similar(x)

    @time begin
        ylen = length(y)
        evalT, oper = rand.evalT, rand.oper
        @views buff = tmp.val.b.coeffs[1][1:ylen]

        # -2xy
        _mul_to!(buff, evalT.Q.Q - 2, y, evalT)
        mul_to!(ctout, x, buff, dgs_pmul, τ_cmul, rand)

        # tmp = x^2
        mul_to!(tmp, x, x, rlk, oper)

        # ctout += tmp
        add_to!(ctout, ctout, tmp, rand)

        # ctout += y^2
        _Bmul_to!(buff, y, y, evalT.Q)
        add_to!(ctout, ctout, buff, rand)

        # rotsum
        # rotation group <5>
        mask_randomize!(ctout, rand)
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], oper)
            add_to!(ctout, ctout, tmp, oper)
        end

        mask_randomize!(ctout, rand)
    end

    ctout
end

function deterministic_flood(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, oper::BFVOperator, pkentor::PKEncryptor)
    ctout = similar(x)
    tmp = similar(x)

    @time begin
        # (x - y)^2
        sub_to!(ctout, x, y, oper)
        mul_to!(ctout, ctout, ctout, rlk, oper)

        # Rotsum, assuming that the rotation group is a subgroup of <5>.
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], oper)
            add_to!(ctout, ctout, tmp, oper)
        end
        flood_to!(ctout, oper, pkentor)
    end

    ctout
end

function flood_to!(x::BFV, oper::BFVOperator, pkentor::PKEncryptor)
    # Parameters
    bits = 154

    evalQ, buff = oper.operQ.evalQ, oper.bfv_buff[end][1:length(x.val)]

    # mask randomize
    rlwe_sample_to!(buff.val, pkentor)
    buff.level[] = x.level[]
    add_to!(x, x, buff, oper)

    # error randomize
    rng = pkentor.rgsampler.rng
    bound = big(1) << bits
    error = ModPoly(rand(rng, -bound:bound, oper.operQ.param.N), evalQ)
    ntt!(error, evalQ)
    add_to!(x.val.b, x.val.b, error, evalQ)
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
    m = 1 << 13
    P, Q = missing, UInt64[0x0040000000006001, 8404993*2147565569]
    dlen, t, ispacking, islevelled = 1, 8404993, true, false
    packlen = 128

    # Generate the basic structs.
    ring_param = CyclotomicParam(m)
    param = BFVParameters(ring_param, P, Q, dlen, t, ispacking, islevelled)
    oper = BFVOperator(param)

    # Generate the secret and evaluation keys.
    sk = ternary_ringkey(UniformSampler(), ring_param.N)
    skentor = BFVEncryptor(sk, 3.2, oper)
    rlk = relin_keygen(skentor)
    rtk = [rotate_keygen((1 << i, 0), skentor) for i = 0:trailing_zeros(packlen)-1]

    # Generate the public key.
    pk = public_keygen(skentor)
    pkentor = PKEncryptor(pk, 91.13921251821257, 7383.043789707919, oper.operQ)

    # Generate the encryptions.
    msg = UInt64[i % 256 for i = 1:packlen]
    packmsg = repeat(msg, (m >> 1) ÷ packlen)   # Sparse packing.
    x = bfv_encrypt(packmsg, skentor)
    y = reverse(packmsg)

    #===========================================#

    # DETERMINISTIC COMPUTATION

    print("DETERMINISTIC COMPUTATION : ")

    ct_deterministic = deterministic(x, y, rlk, rtk, oper)
    out = error(ct_deterministic, skentor)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))

    #===========================================#

    # RANDOMISED COMPUTATION

    print("RANDOMISED COMPUTATION : ")

    ct_randomised = randomised(x, y, rlk, rtk, oper, pkentor)
    out = error(ct_randomised, skentor)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))

    #===========================================#

    # OPTIMISED COMPUTATION

    print("OPTIMISED COMPUTATION : ")

    ct_optimised = optimised(x, y, rlk, rtk, oper, pkentor)
    out = error(ct_optimised, skentor)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))
end

function main2()
    println("IN PREPARATION...")

    # Scheme parameters
    m = 1 << 14
    P, Q = missing, UInt64[0x0000800000020001, 0x000080000008c001, 0x00008000000dc001, 0x00008000000f4001]
    dlen, t, ispacking, islevelled = 1, 8404993, true, false
    packlen = 128

    # Generate the basic structs.
    ring_param = CyclotomicParam(m)
    param = BFVParameters(ring_param, P, Q, dlen, t, ispacking, islevelled)
    oper = BFVOperator(param)

    # Generate the secret and evaluation keys.
    sk = ternary_ringkey(UniformSampler(), ring_param.N)
    skentor = BFVEncryptor(sk, 3.2, oper)
    rlk = relin_keygen(skentor)
    rtk = [rotate_keygen((1 << i, 0), skentor) for i = 0:trailing_zeros(packlen)-1]

    # Generate the public key.
    pk = public_keygen(skentor)
    pkentor = PKEncryptor(pk, 91.13921251821257, 7383.043789707919, oper.operQ)

    # Generate the encryptions.
    msg = UInt64[i % 256 for i = 1:packlen]
    packmsg = repeat(msg, (m >> 1) ÷ packlen)   # Sparse packing.
    x = bfv_encrypt(packmsg, skentor)
    y = reverse(packmsg)

    #===========================================#

    # DETERMINISTIC COMPUTATION

    print("DETERMINISTIC COMPUTATION : ")

    ct_flood = deterministic_flood(x, y, rlk, rtk, oper, pkentor)
    out = error(ct_flood, skentor)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))
end

main()
main2()