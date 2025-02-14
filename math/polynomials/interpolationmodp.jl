function ν(m::Int64, p::Int64)
    res = 0
    while m % p == 0
        m ÷= p
        res += 1
    end
    res
end

function μ(p::Int64, e::Int64)
    i, cnt = p, 0
    while cnt < e
        cnt += ν(i, p)
        i += p
    end
    i - p
end

function finite_diff(x::Vector{UInt64}, y::Vector{UInt64}, pe::Modulus, p::Int64)
    len = length(y)

    res = zeros(UInt64, len)
    tmp1 = UInt64.(y)
    tmp2 = zeros(UInt64, len)

    @inbounds for i = 2 : len
        res[i-1] = tmp1[1]
        for j = 1 : len-i+1
            resj = sub(tmp1[j+1], tmp1[j], pe)
            denom = sub(x[j+i-1], x[j], pe)

            while denom % p == 0
                denom ÷= p
                @assert resj % p == 0
                resj ÷= p
            end

            tmp2[j] = Bmul(resj, invmod(denom, pe.Q), pe)
        end
        @. tmp1 = tmp2
    end
    res[end] = tmp1[1]

    res
end

function interpolate(xi::Vector{<:Integer}, yi::Vector{<:Integer}, p::Int64, e::Int64)
    @assert length(xi) == length(yi)
    @assert length(xi) ≤ μ(p, e)

    pe = Modulus(p^e)
    x = Bred.(xi, Ref(pe))
    y = Bred.(yi, Ref(pe))

    diff = finite_diff(x, y, pe, p)

    len = length(xi)
    res = zeros(UInt64, len)

    basis = zeros(UInt64, len)
    basis[1] = 1
    @inbounds for i = 1 : len
        Bmuladd_to!(res, diff[i], basis, pe)

        if i < len
            @. basis[2:i+1] = basis[1:i]
            basis[1] = 0
            @. (@view basis[1:i]) += (pe.Q - x[i]) * basis[2:i+1]
            Bred.((@view basis[1:i]), Ref(pe))
        end
    end

    res
end

function horner(x::UInt64, nc::Vector{UInt64}, pe::Modulus)
    res = nc[end]
    @inbounds @simd for k = length(nc)-1:-1:1
        res = Bred(nc[k] + widemul(x, res), pe)
    end
    res
end