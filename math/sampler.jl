struct UniformSampler
    rng::ChaChaStream

    UniformSampler() = new(ChaCha12Stream())
end

"""
uniform_binary(N, hw = h) outputs a uniform binary Int64 vector with Hamming weight h.
If hw is not given or set to zero, it outputs a truly uniform binary vector.
"""
uniform_binary(us::UniformSampler, N::Int64, hw::Int64=0) = begin
    if hw == 0
        rand(us.rng, [0, 1], N)
    else
        res = zeros(Int64, N)
        @inbounds for _ = 1:hw
            res[ceil(Int64, rand(us.rng) * N)] = 1
        end
        res
    end
end

"""
uniform_ternary(N, hw = h) outputs a uniform ternary Int64 vector with Hamming weight h.
If hw is not given or set to zero, it outputs a truly uniform ternary vector.
"""
uniform_ternary(us::UniformSampler, N::Int64, hw::Int64=0) = begin
    if hw == 0
        rand(us.rng, [-1, 0, 1], N)
    else
        res = zeros(Int64, N)
        tmp = [-1, 1]
        @inbounds for _ = 1:hw
            res[ceil(Int64, rand(us.rng) * N)] = rand(us.rng, tmp)
        end
        res
    end
end

"""
uniform_random outputs a uniform random distribution.
"""
uniform_random(us::UniformSampler, N::Int64, Q::Modulus) = rand(us.rng, 0:Q.Q-1, N)

uniform_random(us::UniformSampler, Q::Modulus) = rand(us.rng, 0:Q.Q-1)

uniform_random_to!(us::UniformSampler, x::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(x)
        x[i] = rand(us.rng, 0:Q.Q-1)
    end
end

#=============================================================================================#

samplebit(rng::ChaChaStream) = rand(rng, UInt8) & one(UInt64)

const global tailcut = 6

struct CDTSampler
    rng::ChaChaStream       # Base sampler.
    table::Vector{UInt64}   # Probability Table.
    center::Float64         # Distribution center
    sigma::Float64          # Standard Deviation
    taillow::Int64
    tailhigh::Int64
    cint::Int64
    cfrac::Float64

    function CDTSampler(center::Float64, sigma::Float64)
        normFactor = typemax(UInt64)
        cInt = floor(Int64, center)
        cFrac = center - cInt

        tailLow = round(Int64, center - tailcut * sigma)
        tailHigh = round(Int64, center + tailcut * sigma)
        tailCount = tailHigh - tailLow + 1

        table = Vector{UInt64}(undef, tailCount)
        cdf = 0.0

        for i = 0:tailHigh-tailLow
            xf = Float64(i + tailLow)
            rho = exp(-π * (xf - cFrac)^2 / sigma^2) / sigma
            cdf += rho
            table[i+1] = cdf ≥ 1 ? typemax(UInt64) : round(UInt64, cdf * normFactor)
        end

        new(ChaCha12Stream(), table, center, sigma, tailLow, tailHigh, cInt, cFrac)
    end
end

sample(s::CDTSampler) = searchsortedlast(s.table, rand(s.rng, UInt64)) + s.cint + s.taillow

const global baselog = 10
const global sampledepth = 3
const global hipreclog = baselog * sampledepth
const global lowpreclog = 53 - hipreclog

struct VarCenSampler
    basesamplers::Vector{CDTSampler}

    function VarCenSampler(sigma::Float64)
        bk = 0.0
        for i = 0:sampledepth-1
            bk += Float64(1 << baselog)^(-2i)
        end
        sigma /= bk

        basesamplers = Vector{CDTSampler}(undef, 1 << baselog)
        for i = 0:1<<baselog-1
            c = i / (1 << baselog)
            basesamplers[i+1] = CDTSampler(c, sigma)
        end

        new(basesamplers)
    end
end

function sample(s::VarCenSampler, center::Float64)
    cInt = floor(Int64, center)
    cFrac = center - cInt
    cFrac64 = unsafe_trunc(UInt64, cFrac * (1 << 53))

    cFrac64Hi = (cFrac64 >> lowpreclog) % Int64
    r = rand(s.basesamplers[1].rng, UInt64)
    @inbounds for i = reverse(0:lowpreclog-1)
        b = (r >> i) & 1
        cFracBit = (cFrac64 >> i) & 1
        b > cFracBit && (return sampleC(s, cFrac64Hi) + cInt)
        b < cFracBit && (return sampleC(s, cFrac64Hi + 1) + cInt)
    end

    sampleC(s, cFrac64Hi + 1) + cInt
end

function sampleC(s::VarCenSampler, c::Int64)
    mask = (1 << baselog) - 1
    
    @inbounds for i = 0:sampledepth-1
        r = sample(s.basesamplers[c&mask+1])
        (c & mask > 0 && c < 0) && (r -= 1)
        c >>= baselog
        c += r
    end

    c
end

struct RGSampler
    rng::ChaChaStream
    τ::Float64

    function RGSampler()
        new(ChaCha12Stream(), 3.2)
    end

    function RGSampler(τ::Float64)
        new(ChaCha12Stream(), τ)
    end
end

sample(s::RGSampler) = round(Int64, randn(s.rng) * s.τ / √(2π))

sample(s::RGSampler, τ::Float64) = round(Int64, randn(s.rng) * τ / √(2π))

sample128(s::RGSampler, τ::Float64) = round(Int128, randn(s.rng) * τ / √(2π))