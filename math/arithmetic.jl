"""
    find_prime_cyclic(b::Real, m::Int64, n::Int64=1)

Return `n` primes of `b` bits, for cyclic NTT with degree `m`.
"""
function find_prime_cyclic(b::Real, m::Int64, n::Int64=1)
    @assert b ≤ 62 "We do not support primes with bit length bigger than 62-bit"

    primes = Vector{UInt64}(undef, n)

    if is2a3b5c7d(m)
        initialn = round(UInt64, 2^b / 2m) * 2m + 1
        currentn = nextprime(initialn, interval=2m)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + 2m, interval=2m)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=2m)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - 2m, interval=2m)
                cnt += 1
            end
        end
    else
        factor = next2a3b5c7d(2m - 1)
        factor2 = iseven(m) ? lcm(factor, 2m) : lcm(factor, m)
        initialn = round(UInt64, 2^b / factor2) * factor2 + 1
        currentn = nextprime(initialn, interval=factor2)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + factor2, interval=factor2)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=factor2)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - factor2, interval=factor2)
                cnt += 1
            end
        end
    end

    Tuple(primes)
end

"""
    find_prime_cyclotomic(b::Real, m::Int64, n::Int64=1)

Return `m` primes of `b`-bits, for cyclotomic NTT with degree `m`.
"""
function find_prime_cyclotomic(b::Real, m::Int64, n::Int64=1)
    @assert b ≤ 62 "We do not support primes with bit length bigger than 62-bit"

    primes = Vector{UInt64}(undef, n)

    if ispow2(m)
        initialn = round(UInt64, 2^b / m) * m + 1
        currentn = nextprime(initialn, interval=m)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + m, interval=m)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=m)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - m, interval=m)
                cnt += 1
            end
        end
    else
        bluestein = next2a3b5c7d(2m - 1)

        N = totient(m)
        factors = factor(Vector, m)
        p = filter(isodd, factors)[1]
        moverp = iseven(m) ? m ÷ 2p : m ÷ p
        deg = iseven(m) ? m ÷ 2 - moverp : m - moverp
        if deg == N
            fact = 2lcm(bluestein, m)
        else
            ñ = next2a3b5c7d(m)
            A = next2a3b5c7d(2(deg - N) + 1)
            fact = 2lcm(bluestein, m, ñ, A)
        end

        initialn = round(UInt64, 2^b / fact) * fact + 1
        currentn = nextprime(initialn, interval=fact)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + fact, interval=fact)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=fact)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - fact, interval=fact)
                cnt += 1
            end
        end
    end

    Tuple(primes)
end

"""
    ord(p::Integer, m::Integer)

Returns the multiplicative order of p modulo m.
"""
function ord(p::Integer, m::Integer)::Int64
    res = 1
    acc = p % m
    while acc != 1
        acc = (acc * p) % m
        res += 1
    end

    res
end

@inline zeropadto(a::Vector{T}, n::Int64) where {T} = begin
    @assert n ≥ length(a)
    vcat(a, zeros(T, n - length(a)))
end

@inline ispow3(n::Int64)::Bool = n == 3^round(Int64, log(n) / log(3))

@inline ispow5(n::Int64)::Bool = n == 5^round(Int64, log(n) / log(5))

@inline ispow7(n::Int64)::Bool = n == 7^round(Int64, log(n) / log(7))

@noinline is2a3b5c7d(n::Int64)::Bool = keys(factor(Dict, n)) ⊆ [2, 3, 5, 7]

@noinline next2a3b5c7d(n::Int64)::Int64 = begin
    res = n + 1
    while !is2a3b5c7d(res)
        res += 1
    end
    res
end

@inline factor2357(n::Int64)::NTuple{4,Int64} = begin
    n2, n3, n5, n7 = 1, 1, 1, 1

    while n % 2 == 0
        n >>= 1
        n2 <<= 1
    end

    while n % 3 == 0
        n ÷= 3
        n3 *= 3
    end

    while n % 5 == 0
        n ÷= 5
        n5 *= 5
    end

    while n % 7 == 0
        n ÷= 7
        n7 *= 7
    end

    n2, n3, n5, n7
end

# This parameter guarantees at most -52 b of fixed point error (with 64 levels). 
# The decryption failure will happen in a very unlikely situation...
const fixed_prec = 124
const float_prec = 62
const fixed_mask = UInt64(1) << float_prec - 1
const round_mask = UInt64(1) << (fixed_prec - float_prec) - 1

mult_and_round(a::UInt64, b::UInt128) = a * (b >> float_prec) + ((a * (b & fixed_mask)) >> float_prec)

round_to_uint64(a::UInt128) = round_to_uint128(a) % UInt64

round_to_uint128(a::UInt128) = (a >> (fixed_prec - float_prec) + (a & round_mask) >> (fixed_prec - float_prec - 1))

"""
scramble!(a, r) computes an in-place r-radix reversal algorithm. 
The length of the input vector should be a power-of-r.
"""
function scramble!(a::AbstractVector{UInt64}, r::Int64)
    if r == 2
        j = 0
        N = length(a)
        @inbounds for i = 1:N-1
            bit = N >> 1
            while j ≥ bit
                j -= bit
                bit >>= 1
            end
            j += bit
            if i < j
                a[i+1], a[j+1] = a[j+1], a[i+1]
            end
        end
    else
        logN = floor(Int64, log(length(a)) / log(r))
        N = r^logN

        @inbounds for i = 1:N-1
            idx, j = i, 0
            @simd for _ = 1:logN
                j = j * r + (idx % r)
                idx ÷= r
            end
            if j > i
                a[i+1], a[j+1] = a[j+1], a[i+1]
            end
        end
    end
end

"""
primitive_root_finder(Q) finds the primitive root of Q.
"""
function primitive_root_finder(Q::Modulus)::UInt64
    #TODO This code can be improved.

    test = (Q.Q - 1) .÷ collect(keys(factor(Dict, Q.Q - 1)))

    one = Mform(1, Q)
    g = one
    while true
        if all(powermod.(g, test, Ref(Q)) .!= one)
            break
        end
        g += 1

        if g > Q.Q
            throw("$Q.Q may not be a prime nber.")
        end
    end

    g
end

function primitive_root_finder(Q::Int64)::Int64
    test = (Q - 1) .÷ collect(keys(factor(Dict, Q - 1)))

    g = 2
    while true
        if all(powermod.(g, test, Q) .!= 1)
            break
        end
        g += 1

        if g > Q
            throw("$Q may not be a prime nber.")
        end
    end

    g
end

function find_generators_mod_m(m::Int64)
    fact_dict = factor(Dict, m)
    facs = collect(fact_dict)
    sort!(facs, by=x -> x[1])

    dims = Memory{Int64}(undef, length(facs))
    gens = Memory{Int64}(zeros(Int64, length(facs)))

    @inbounds for i = eachindex(facs)
        p, e = facs[i].first, facs[i].second

        for j = eachindex(facs)
            i == j && continue
            pj, ej = facs[j].first, facs[j].second
            gens[i] += invmod(m ÷ pj^ej, pj^ej) * m ÷ pj^ej
        end

        dims[i] = p^e - p^(e - 1)

        for j = 2:p^e-1
            j % p == 0 && continue
            if ord(j, p^e) == dims[i]
                gens[i] += j * invmod(m ÷ p^e, p^e) * m ÷ p^e
                break
            end
        end

        gens[i] %= m
    end

    dims, gens
end

"""
division(a, b) performs the long division of poynomials.
This function affects the vector a, so be cautious when using.
"""
function division(a::Vector{Int64}, b::Vector{Int64})::Vector{Int64}
    @assert length(a) ≥ length(b)

    Q = zeros(Int64, length(a) - length(b) + 1)

    @inbounds for i = 0:length(a)-length(b)
        if a[end-i] != 0
            Q[end-i] = a[end-i] ÷ b[end]

            @simd for j = 0:length(b)-1
                a[end-i-j] -= b[end-j] * Q[end-i]
            end
        end
    end

    Q
end

"""
cyclotomic_finder(m) returns the m-th cyclotomic polynomial. 
"""
function cyclotomic_finder(m::Int64)::Vector{Int64}
    # Find the factors of m. Sorting is necessary to enhance the efficiency.
    factors = factor(Dict, m)
    primes = collect(keys(factors))
    sort!(primes)

    if primes[1] == 2
        is2divm = true
        popfirst!(primes)
    else
        is2divm = false
    end

    # Make lazy_Mmul arrays for less allocations.
    ϕm = zeros(Int64, m + 1)
    ϕm[1] = -1
    ϕm[2] = 1
    tmp = Vector{Int64}(undef, m + 1)
    tmp2 = Vector{Int64}(undef, m + 1)

    # Perform naive polynomial long divisions.
    η = 1
    ηprime = 1
    @inbounds for p = primes
        @. tmp = ϕm
        @. tmp2 = 0
        @. ϕm = 0
        @simd for i = 0:η
            tmp2[i*p+1] = tmp[i+1]
        end

        ηtimesp = η * p
        ηprime = ηtimesp - η

        for i = 0:(p-1)*η
            if tmp2[ηtimesp+1-i] != 0
                ϕm[ηprime+1-i] = tmp2[ηtimesp+1-i] ÷ tmp[η+1]

                @simd for j = 0:η
                    tmp2[ηtimesp+1-i-j] -= tmp[η+1-j] * ϕm[ηprime+1-i]
                end
            end
        end

        η = ηprime
    end

    s = m ÷ reduce(*, primes)

    if is2divm
        @inbounds @simd for i = 1:2:η
            ϕm[i+1] = -ϕm[i+1]
        end

        s >>= 1
    end

    if s > 1
        @. tmp = ϕm
        @. ϕm = 0
        @inbounds @simd for i = 0:η
            ϕm[s*i+1] = tmp[i+1]
        end
    end

    resize!(ϕm, totient(m) + 1)
    ϕm
end