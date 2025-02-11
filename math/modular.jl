"""
Modulus supports a modular arithmetic given modulus Q.
"""
struct Modulus
    Q::UInt64
    Q⁻¹::UInt64       # Used for Montgomery reduction.
    R⁻¹::UInt64       # Used for Montgomery reduction.
    r0::UInt64        # Used for Barrett reduction.
    r1::UInt64        # Used for Barrett reduction.
    logQ::Int64       # Used for gadget decomposition.
    halfQ::Int64      # Used for balanced representation.

    function Modulus(Q::Integer)
        logQ = ceil(Int64, log2(Q))
        @assert logQ ≤ 62 "Modulus should be smaller than 2⁶²."

        if isodd(Q)
            Q⁻¹ = UInt64(invmod(-Q, 0x00000000000000010000000000000000))
            R⁻¹ = UInt64(invmod(0x00000000000000010000000000000000, Q))
        else
            Q⁻¹ = UInt64(0)
            R⁻¹ = UInt128(0)
        end

        r = floor(UInt128, (big(1) << 128) / Q)
        r0, r1 = r % UInt64, (r >> 64) % UInt64
        halfQ = (Q - 1) ÷ 2
        new(Q, Q⁻¹, R⁻¹, r0, r1, logQ, halfQ)
    end
end

"""
Bred(x, Q) returns x % Q using barret reduction.
"""
Bred(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) >> 64) % UInt64
    res = x - Q.Q * t
    res ≥ Q.Q ? res - Q.Q : res
end

Bred(x::Int64, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return Bred(UInt64(x), Q)
    else
        res = Bred(UInt64(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

Bred(x::UInt128, Q::Modulus)::UInt64 = begin
    x0 = x % UInt64
    x1 = (x >> 64) % UInt64
    t = (widemul(Q.r1, x1) + (widemul(Q.r1, x0) + widemul(Q.r0, x1) + widemul(Q.r0, x0) >> 64) >> 64) % UInt64
    res = x0 - t * Q.Q
    res ≥ Q.Q ? res - Q.Q : res
end

Bred(x::Int128, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return Bred(UInt128(x), Q)
    else
        res = Bred(UInt128(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

lazy_Bred(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) >> 64) % UInt64
    x - Q.Q * t
end

lazy_Bred(x::UInt128, Q::Modulus)::UInt64 = begin
    x0 = x % UInt64
    x1 = (x >> 64) % UInt64
    t = (widemul(Q.r1, x1) + (widemul(Q.r1, x0) + widemul(Q.r0, x1) + widemul(Q.r0, x0) >> 64) >> 64) % UInt64
    x0 - t * Q.Q
end

# Taken from the invmod.
Base.:invmod(x::Integer, Q::Modulus)::UInt64 = begin
    n, m = UInt64(x), Q.Q
    g, x, _ = gcdx(n, m)
    g != 1 && throw(DomainError((n, m), LazyString("Greatest common divisor is ", g, ".")))
    x > typemax(UInt64) >> 1 && (x += m)
    Bred(x, Q)
end

Mform(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) + widemul(Q.r0, x) >> 64) % UInt64
    res = -t * Q.Q
    res ≥ Q.Q ? res - Q.Q : res
end

Mform(x::Int64, Q::Modulus)::UInt64 = Mform(Bred(x, Q), Q)

iMform(x::UInt64, Q::Modulus)::UInt64 = Bred(widemul(x, Q.R⁻¹), Q)

Mform!(x::AbstractVector{UInt64}, Q::Modulus) = Mform_to!(x, x, Q)

Mform_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(x)
        res[i] = Mform(x[i], Q)
    end
end

iMform!(x::AbstractVector{UInt64}, Q::Modulus) = iMform_to!(x, x, Q)

iMform_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(x)
        res[i] = iMform(x[i], Q)
    end
end

"""
Convert UInt64 to a balanced representation.
"""
Base.:signed(x::UInt64, Q::Modulus)::Int64 = signed(x > Q.halfQ ? x - Q.Q : x)

"""
Convert (x mod P) mod Q in balanced representation.
"""
embed_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, P::Modulus, Q::Modulus) = begin
    @inbounds for i = eachindex(x)
        res[i] = Bred(signed(x[i], P), Q)
    end
end

add(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    res = a + b
    res ≥ Q.Q ? res - Q.Q : res
end

add_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = add(a[i], b[i], Q)
    end
end

add_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::UInt64, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = add(a[i], b, Q)
    end
end

neg(a::UInt64, Q::Modulus)::UInt64 = a == 0 ? a : Q.Q - a

neg_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = neg(a[i], Q)
    end
end

sub(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    res = a - b
    res ≥ Q.Q ? res + Q.Q : res
end

sub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = sub(a[i], b[i], Q)
    end
end

sub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::UInt64, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = sub(a[i], b, Q)
    end
end

Mmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    q, q⁻¹ = Q.Q, Q.Q⁻¹

    ab = widemul(a, b)
    w = ((widemul(q, (ab % UInt64) * q⁻¹) + ab) >> 64) % UInt64
    w ≥ q ? w - q : w
end

lazy_Mmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    q, q⁻¹ = Q.Q, Q.Q⁻¹

    ab = widemul(a, b)
    ((widemul(q, (ab % UInt64) * q⁻¹) + ab) >> 64) % UInt64
end

Bmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = Bred(widemul(a, b), Q)

lazy_Bmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = lazy_Bred(widemul(a, b), Q)

Mmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = Mmul(a[i], b[i], Q)
    end
end

Mmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(b)
        res[i] = Mmul(a, b[i], Q)
    end
end

Bmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(b)
        res[i] = Bmul(a, b[i], Q)
    end
end

Mmuladd_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = add(res[i], Mmul(a[i], b[i], Q), Q)
    end
end

Mmuladd_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(b)
        res[i] = add(res[i], Mmul(a, b[i], Q), Q)
    end
end

Bmuladd_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(b)
        res[i] = Bred(widemul(a[i], b[i]) + res[i], Q)
    end
end

Bmuladd_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(b)
        res[i] = Bred(widemul(a, b[i]) + res[i], Q)
    end
end

Mmulsub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = sub(res[i], Mmul(a[i], b[i], Q), Q)
    end
end

Mmulsub_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(b)
        res[i] = sub(res[i], Mmul(a, b[i], Q), Q)
    end
end

Bmulsub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @inbounds for i = eachindex(a)
        res[i] = Bred(widemul(Q.Q - a[i], b[i]) + res[i], Q)
    end
end

Bmulsub_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    neg_a = Q.Q - a
    @inbounds for i = eachindex(a)
        res[i] = Bred(widemul(neg_a, b[i]) + res[i], Q)
    end
end

Base.:sum(a::AbstractVector{UInt64}, Q::Modulus)::UInt64 = begin
    res = zero(UInt64)
    @inbounds for i = eachindex(a)
        res = add(res, a[i], Q)
    end
    res
end

# Input is in Mform, output is in Mform.
Base.:powermod(a::UInt64, p::Integer, Q::Modulus)::UInt64 = begin
    @assert p ≥ 0
    p == 0 && return Mform(1, Q)

    t = prevpow(2, p)
    r = Mform(1, Q)
    while true
        if p >= t
            r = Mmul(r, a, Q)
            p -= t
        end
        t >>>= 1
        t <= 0 && break
        r = Mmul(r, r, Q)
    end
    return r
end

function reduce!(a::AbstractVector{UInt64}, Q::Modulus)
    @inbounds for i = eachindex(a)
        a[i] = Bred(a[i], Q)
    end
end

const Moduli = AbstractVector{Modulus}