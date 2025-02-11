"""
PolyEvaluator is a struct for the modulo polynomial arithmetic over a CRT domain.
"""
mutable struct PolyEvaluator
    moduli::Vector{Modulus}
    ntter::Vector{<:Transformer}
    ismod::Vector{Bool}
    buff::Array{UInt64,2}
    autbuff::Vector{UInt64}

    function PolyEvaluator(ntter::AbstractVector{<:Transformer})
        moduli = [ntteri.Q for ntteri = ntter]
        ismod = fill(true, length(ntter))
        buff = Array{UInt64,2}(undef, ntter[1].N, 2)
        
        if typeof(ntter[1]) == CyclotomicTransformer_pow2
            autbuff = Vector{UInt64}(undef, ntter[1].N)
        elseif typeof(ntter[1]) == CyclotomicTransformer_arb
            autbuff = Vector{UInt64}(undef, ntter[1].m)
        end
        
        new(moduli, collect(ntter), ismod, buff, autbuff)
    end

    function PolyEvaluator(moduli::Moduli, ntter::AbstractVector{<:Transformer})
        @assert length(moduli) == length(ntter) "The length of moduli and transformers should match."

        ismod = [moduli[i].Q == ntter[i].Q.Q for i = eachindex(moduli)]
        buff = Array{UInt64,2}(undef, ntter[1].N, 2)
       
        if typeof(ntter[1]) == CyclotomicTransformer_pow2
            autbuff = Vector{UInt64}(undef, ntter[1].N)
        elseif typeof(ntter[1]) == CyclotomicTransformer_arb
            autbuff = Vector{UInt64}(undef, ntter[1].m)
        end

        new(collect(moduli), collect(ntter), ismod, buff, autbuff)
    end
end

#=================================================================================================================#

"""
ModScalar is a struct to store scalar in RNS form.
"""
struct ModScalar
    vals::Vector{UInt64}
    len::Int64
    isMform::RefBool

    ModScalar(vals::AbstractVector{UInt64}, isMform::Bool) =
        new(collect(vals), length(vals), Ref(isMform))

    ModScalar(len::Int64, isMform::Bool=false) =
        new(zeros(UInt64, len), len, Ref(isMform))

    ModScalar(val::Union{Int64,UInt64,Int128,UInt128}, eval::PolyEvaluator, isMform::Bool=false) = begin
        moduli, ismod = eval.moduli, eval.ismod

        vals = Vector{UInt64}(undef, length(moduli))
        @inbounds for i = eachindex(vals)
            if ismod[i] && isMform
                vals[i] = Mform(val, moduli[i]) 
            else
                vals[i] = Bred(val, moduli[i])
            end
        end

        new(vals, length(moduli), Ref(isMform))
    end

    ModScalar(val::BigInt, eval::PolyEvaluator, isMform::Bool=false) = begin
        moduli, ismod = eval.moduli, eval.ismod

        vals = Vector{UInt64}(undef, length(moduli))
        @inbounds for i = eachindex(vals)
            if ismod[i] && isMform
                vals[i] = Mform((val % moduli[i].Q) % UInt64, moduli[i]) 
            else
                vals[i] = (val % moduli[i].Q) % UInt64
            end
        end

        new(vals, length(moduli), Ref(isMform))
    end
end

Base.copy!(dst::ModScalar, src::ModScalar) = begin
    @assert dst.len == src.len "The length of the destination scalar should be same to the input scalar."
    copy!(dst.vals, src.vals)
    dst.isMform[] = src.isMform[]
end

initialise!(s::ModScalar, isMform::Bool=true) = begin
    @. s.vals = zero(UInt64)
    s.isMform[] = isMform
end

Mform!(x::ModScalar, eval::PolyEvaluator) = Mform_to!(x, x, eval)

@views Mform_to!(res::ModScalar, x::ModScalar, eval::PolyEvaluator) = begin
    @assert !x.isMform[] "Scalar is already in Montgomery form."
    @assert res.len == x.len "the length of input and output scalars should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:res.len
        !eval.ismod[i] && continue
        res.vals[i] = Mform(x.vals[i], eval.moduli[i])
    end

    res.isMform[] = true
end

iMform!(x::ModScalar, eval::PolyEvaluator) = iMform_to!(x, x, eval)

@views iMform_to!(res::ModScalar, x::ModScalar, eval::PolyEvaluator) = begin
    @assert !x.isMform[] "Scalar is already in Montgomery form."
    @assert res.len == x.len "the length of input and output scalars should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:res.len
        !eval.ismod[i] && continue
        res.vals[i] = iMform(x.vals[i], eval.moduli[i])
    end

    res.isMform[] = true
end

neg(x::ModScalar, eval::PolyEvaluator) = begin
    res = ModScalar(x.len)
    neg_to!(res, x, eval)
end

@views neg_to!(res::ModScalar, x::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == x.len "the length of input and output scalars should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:x.len
        res.vals[i] = neg(x.vals[i], eval.moduli[i])
    end

    res.isMform[] = x.isMform[]
end

add(x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    res = ModScalar(x.len)
    add_to!(res, x, y, eval)
end

add_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == x.len == y.len "the length of input and output scalars should be the same."
    @assert x.isMform[] == y.isMform[] "Input scalars should be both in or not in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:x.len
        res.vals[i] = add(x.vals[i], y.vals[i], eval.moduli[i])
    end

    res.isMform[] = x.isMform[]
end

sub(x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    res = ModScalar(x.len)
    sub_to!(res, x, y, eval)
end

sub_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == x.len == y.len "The length of input and output scalars should be the same."
    @assert x.isMform[] == y.isMform[] "Input scalars should be both in or not in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:x.len
        res.vals[i] = sub(x.vals[i], y.vals[i], eval.moduli[i])
    end

    res.isMform[] = x.isMform[]
end

mul(x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    res = ModScalar(x.N, x.len)
    mul_to!(res, x, y, eval)
end

mul_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == x.len == y.len "The length of input and output scalars should be the same."
    @assert x.isMform[] && y.isMform[] "Input scalars should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            res.vals[i] = Mmul(x.vals[i], y.vals[i], eval.moduli[i])
        else
            res.vals[i] = Bmul(x.vals[i], y.vals[i], eval.moduli[i])
        end
    end

    res.isMform[] = true
end

muladd_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == x.len == y.len "The length of input and output scalars should be the same."
    @assert x.isMform[] && y.isMform[] "Input scalars should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            res.vals[i] = add(res.vals[i], Mmul(x.vals[i], y.vals[i], eval.moduli[i]), eval.moduli[i])
        else
            res.vals[i] = Bred(widemul(x.vals[i], y.vals[i]) + res.vals[i], eval.moduli[i])
        end
    end

    res.isMform[] = true
end

mulsub_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == x.len == y.len "The length of input and output scalars should be the same."
    @assert x.isMform[] && y.isMform[] "Input scalars should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            res.vals[i] = sub(res.vals[i], Mmul(x.vals[i], y.vals[i], eval.moduli[i]), eval.moduli[i])
        else
            res.vals[i] = Bred(widemul(x.vals[i], y.vals[i]) + eval.moduli[i].Q - res.vals[i], eval.moduli[i])
        end
    end

    res.isMform[] = true
end

reduce!(x::ModScalar, eval::PolyEvaluator) = begin
    @inbounds for i = 1:x.len
        x.vals[i] = Bred(x.vals[i], eval.moduli[i])
    end
end

# There can be more optimisations.
function to_big(x::ModScalar, eval::PolyEvaluator)
    ismod, moduli = eval.ismod, eval.moduli

    @assert x.len == length(moduli) "The lengths of input scalar and moduli do not match."

    Qlen = x.len

    res = big(0)
    Qbig = prod([big(Qi.Q) for Qi = moduli])

    @inbounds for i = 1:Qlen
        Qtilde = Qbig ÷ moduli[i].Q
        if ismod[i] && x.isMform[]
            res += (iMform(x.vals[i], moduli[i]) * invmod(Qtilde, moduli[i].Q) % moduli[i].Q) * Qtilde
        else
            res += (x.vals[i] * invmod(Qtilde, moduli[i].Q) % moduli[i].Q) * Qtilde
        end
    end

    res %= Qbig
    if res > (Qbig >> 1)
        res -= Qbig
    end

    res
end

#=================================================================================================================#

"""
ModPoly is a struct to store the polynomial in RNS form.
"""
struct ModPoly
    coeffs::Array{UInt64,2}
    N::Int64
    len::Int64
    isntt::RefBool
    isMform::RefBool

    ModPoly(coeffs::AbstractArray{UInt64,2}, isntt::Bool, isMform::Bool) =
        new(collect(coeffs), size(coeffs)[1], size(coeffs)[2], Ref(isntt), Ref(isMform))

    ModPoly(N::Int64, len::Int64, isntt::Bool=true, isMform::Bool=true) =
        new(zeros(UInt64, N, len), N, len, Ref(isntt), Ref(isMform))

    ModPoly(coeff::Vector{<:Union{Int64,UInt64}}, eval::PolyEvaluator, isMform::Bool=false) = begin
        moduli, ismod = eval.moduli, eval.ismod

        coeffs = Array{UInt64}(undef, length(coeff), length(moduli))
        @inbounds for j = eachindex(moduli)
            if ismod[j] && isMform
                for i = eachindex(coeff)
                    coeffs[i, j] = Mform(coeff[i], moduli[j])
                end
            else
                for i = eachindex(coeff)
                    coeffs[i, j] = Bred(coeff[i], moduli[j])
                end
            end
        end

        new(coeffs, length(coeff), length(moduli), Ref(false), Ref(isMform))
    end

    ModPoly(coeff::Vector{BigInt}, eval::PolyEvaluator, isMform::Bool=false) = begin
        moduli, ismod = eval.moduli, eval.ismod

        coeffs = Array{UInt64}(undef, length(coeff), length(moduli))
        @inbounds for j = eachindex(moduli)
            if ismod[j] && isMform
                for i = eachindex(coeff)
                    coeffs[i, j] = Mform((coeff[i] % moduli[j].Q) % UInt64, moduli[j])
                end
            else
                for i = eachindex(coeff)
                    coeffs[i, j] = (coeff[i] % moduli[j].Q) % UInt64
                end
            end
        end

        new(coeffs, length(coeff), length(moduli), Ref(false), Ref(isMform))
    end
end

Base.copy!(dst::ModPoly, src::ModPoly) = begin
    @assert dst.N == src.N "The ring dimension of the destination polynomial should be same to the input polynomial."
    @assert dst.len == src.len "The length of the destination polynomial should be same to the input polynomial."
    copy!(dst.coeffs, src.coeffs)
    dst.isMform[] = src.isMform[]
    dst.isntt[] = src.isntt[]
end

initialise!(p::ModPoly, isMform::Bool=true, isntt::Bool=true) = begin
    @. p.coeffs = zero(UInt64)
    p.isMform[] = isMform
    p.isntt[] = isntt
end

Mform!(x::ModPoly, eval::PolyEvaluator) = Mform_to!(x, x, eval)

@views Mform_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator) = begin
    @assert !x.isMform[] "Polynomial is already in Montgomery form."
    @assert res.N == x.N "The ring dimension of the input and output polynomials should be the same."
    @assert res.len == x.len "the length of input and output polynomials should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:res.len
        !eval.ismod[i] && continue
        Mform_to!(res.coeffs[:, i], x.coeffs[:, i], eval.moduli[i])
    end

    res.isMform[] = true
end

iMform!(x::ModPoly, eval::PolyEvaluator) = iMform_to!(x, x, eval)

@views iMform_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator) = begin
    @assert x.isMform[] "Polynomial is not in Montgomery form."
    @assert res.N == x.N "The ring dimension of the input and output polynomials should be the same."
    @assert res.len == x.len "the length of input and output polynomials should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:res.len
        !eval.ismod[i] && continue
        iMform_to!(res.coeffs[:, i], x.coeffs[:, i], eval.moduli[i])
    end

    res.isMform[] = false
end

neg(x::ModPoly, eval::PolyEvaluator) = begin
    res = ModPoly(x.N, x.len)
    neg_to!(res, x, eval)
end

@views neg_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator) = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    @assert res.len == x.len "the length of input and output polynomials should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        neg_to!(res.coeffs[:, i], x.coeffs[:, i], eval.moduli[i])
    end

    res.isntt[] = x.isntt[]
    res.isMform[] = x.isMform[]
end

add(x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    res = ModPoly(x.N, x.len)
    add_to!(res, x, y, eval)
end

@views add_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert res.N == x.N == y.N "The polynomials should have the same ring degree."
    @assert res.len == x.len == y.len "the length of input and output polynomials should be the same."
    @assert x.isntt[] == y.isntt[] "Input polynomials should be both in coeff form or NTT form."
    @assert x.isMform[] == y.isMform[] "Input polynomials should be both in or not in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        add_to!(res.coeffs[:, i], x.coeffs[:, i], y.coeffs[:, i], eval.moduli[i])
    end

    res.isntt[] = x.isntt[]
    res.isMform[] = x.isMform[]
end

@views add_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    @assert res.len == x.len == y.len "The length of input and output polynomials should be the same."
    @assert res.isMform[] == y.isMform[] "Input polynomial and scalar should be both in or not in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if x.isntt[] && eval.ismod[i]
            add_to!(res.coeffs[:, i], x.coeffs[:, i], y.vals[i], eval.moduli[i])
        else
            res.coeffs[1, i] = add(x.coeffs[1, i], y.vals[i], eval.moduli[i])
        end
    end

    res.isntt[] = x.isntt[]
    res.isMform[] = x.isMform[]
end

sub(x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    res = ModPoly(x.N, x.len)
    sub_to!(res, x, y, eval)
end

@views sub_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert res.N == x.N == y.N "The polynomials should have the same ring degree."
    @assert res.len == x.len == y.len "The length of input and output polynomials should be the same."
    @assert x.isntt[] == y.isntt[] "Input polynomials should be both in coeff form or NTT form."
    @assert x.isMform[] == y.isMform[] "Input polynomials should be both in or not in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        sub_to!(res.coeffs[:, i], x.coeffs[:, i], y.coeffs[:, i], eval.moduli[i])
    end

    res.isntt[] = x.isntt[]
    res.isMform[] = x.isMform[]
end

@views sub_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    @assert res.len == x.len == y.len "The length of input and output polynomials should be the same."
    @assert res.isMform[] == y.isMform[] "Input polynomial and scalar should be both in or not in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if x.isntt[] && eval.ismod[i]
            sub_to!(res.coeffs[:, i], x.coeffs[:, i], y.vals[i], eval.moduli[i])
        else
            res.coeffs[1, i] = sub(x.coeffs[1, i], y.vals[i], eval.moduli[i])
        end
    end

    res.isntt[] = x.isntt[]
    res.isMform[] = x.isMform[]
end

mul(x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    res = ModPoly(x.N, x.len)
    mul_to!(res, x, y, eval)
end

@views mul_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert x.N == y.N == res.N "The input and output polynomials should have the same length."
    @assert res.len == x.len == y.len "The length of input and output polynomials should be the same."
    @assert x.isntt[] && y.isntt[] "Input polynomials should be both in NTT form."
    @assert x.isMform[] && y.isMform[] "Input polynomials should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            Mmul_to!(res.coeffs[:, i], x.coeffs[:, i], y.coeffs[:, i], eval.moduli[i])
        else
            Mform_to!(eval.buff[:, 1], x.coeffs[:, i], eval.ntter[i].Q)
            Mform_to!(eval.buff[:, 2], y.coeffs[:, i], eval.ntter[i].Q)
            ntt!(eval.buff[:, 1], eval.ntter[i])
            ntt!(eval.buff[:, 2], eval.ntter[i])
            Mmul_to!(eval.buff[:, 1], eval.buff[:, 1], eval.buff[:, 2], eval.ntter[i].Q)
            intt!(eval.buff[:, 1], eval.ntter[i])
            iMform!(eval.buff[:, 1], eval.ntter[i].Q)
            embed_to!(res.coeffs[:, i], eval.buff[:, 1], eval.ntter[i].Q, eval.moduli[i])
        end
    end

    res.isMform[] = true
    res.isntt[] = true
end

@views mul_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    @assert y.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert res.len == y.len == x.len "The length of input and output polynomials should be the same."
    @assert x.isMform[] && y.isMform[] "Input polynomials should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:y.len
        if eval.ismod[i]
            Mmul_to!(res.coeffs[:, i], x.vals[i], y.coeffs[:, i], eval.moduli[i])
        else
            Bmul_to!(res.coeffs[:, i], x.vals[i], y.coeffs[:, i], eval.moduli[i])
        end
    end

    res.isMform[] = true
    res.isntt[] = y.isntt[]
end

@views muladd_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert x.N == y.N == res.N "The input and output polynomials should have the same length."
    @assert res.len == x.len == y.len "The length of input and output polynomials should be the same."
    @assert x.isntt[] && y.isntt[] "Input polynomials should be both in NTT form."
    @assert x.isMform[] && y.isMform[] "Input polynomials should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            Mmuladd_to!(res.coeffs[:, i], x.coeffs[:, i], y.coeffs[:, i], eval.moduli[i])
        else
            Mform_to!(eval.buff[:, 1], x.coeffs[:, i], eval.ntter[i].Q)
            Mform_to!(eval.buff[:, 2], y.coeffs[:, i], eval.ntter[i].Q)
            ntt!(eval.buff[:, 1], eval.ntter[i])
            ntt!(eval.buff[:, 2], eval.ntter[i])
            Mmul_to!(eval.buff[:, 1], eval.buff[:, 1], eval.buff[:, 2], eval.ntter[i].Q)
            intt!(eval.buff[:, 1], eval.ntter[i])
            iMform!(eval.buff[:, 1], eval.ntter[i].Q)
            embed_to!(eval.buff[:, 1], eval.buff[:, 1], eval.ntter[i].Q, eval.moduli[i])
            add_to!(res.coeffs[:, i], res.coeffs[:, i], eval.buff[:, 1], eval.moduli[i])
        end
    end

    res.isntt[] = true
    res.isMform[] = true
end

@views muladd_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    @assert y.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert res.len == y.len == x.len "The length of input and output polynomials should be the same."
    @assert x.isMform[] && y.isMform[] "Input polynomials should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            Mmuladd_to!(res.coeffs[:, i], x.vals[i], y.coeffs[:, i], eval.moduli[i])
        else
            Bmuladd_to!(res.coeffs[:, i], x.vals[i], y.coeffs[:, i], eval.moduli[i])
        end
    end

    res.isntt[] = y.isntt[]
    res.isMform[] = true
end

@views muladd_to!(res::ModPoly, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == y.len == x.len "The length of input and output polynomials should be the same."
    @assert res.isMform[] && x.isMform[] && y.isMform[] "Input scalars and output polynomial should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for j = 1:x.len
        if eval.ismod[j]
            tmp = Mmul(x.vals[j], y.vals[j], eval.moduli[j])
            if res.isntt[]
                for i = 1:res.N
                    res.coeffs[i, j] = add(res.coeffs[i, j], tmp, eval.moduli[j])
                end
            else
                res.coeffs[1, j] = add(res.coeffs[1, j], tmp, eval.moduli[j])
            end
        else
            tmp = Bmul(x.vals[j], y.vals[j], eval.moduli[j])
            res.coeffs[1, j] = add(res.coeffs[1, j], tmp, eval.moduli[j])
        end
    end
end

@views mulsub_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert x.N == y.N == res.N "The input and output polynomials should have the same length."
    @assert res.len == x.len == y.len "The length of input and output polynomials should be the same."
    @assert x.isntt[] && y.isntt[] "Input polynomials should be both in NTT form."
    @assert x.isMform[] && y.isMform[] "Input polynomials should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            Mmulsub_to!(res.coeffs[:, i], x.coeffs[:, i], y.coeffs[:, i], eval.moduli[i])
        else
            Mform_to!(eval.buff[:, 1], x.coeffs[:, i], eval.ntter[i].Q)
            Mform_to!(eval.buff[:, 2], y.coeffs[:, i], eval.ntter[i].Q)
            ntt!(eval.buff[:, 1], eval.ntter[i])
            ntt!(eval.buff[:, 2], eval.ntter[i])
            Mmul_to!(eval.buff[:, 1], eval.buff[:, 1], eval.buff[:, 2], eval.ntter[i].Q)
            intt!(eval.buff[:, 1], eval.ntter[i])
            iMform!(eval.buff[:, 1], eval.ntter[i].Q)
            embed_to!(eval.buff[:, 1], eval.buff[:, 1], eval.ntter[i].Q, eval.moduli[i])
            sub_to!(res.coeffs[:, i], res.coeffs[:, i], eval.buff[:, 1], eval.moduli[i])
        end
    end

    res.isntt[] = true
    res.isMform[] = true
end

@views mulsub_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    @assert y.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert res.len == y.len == x.len "The length of input and output polynomials should be the same."
    @assert x.isMform[] && y.isMform[] "Input polynomials should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            Mmulsub_to!(res.coeffs[:, i], x.vals[i], y.coeffs[:, i], eval.moduli[i])
        else
            Bmulsub_to!(res.coeffs[:, i], x.vals[i], y.coeffs[:, i], eval.moduli[i])
        end
    end

    res.isntt[] = y.isntt[]
    res.isMform[] = true
end

@views mulsub_to!(res::ModPoly, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.len == y.len == x.len "The length of input and output polynomials should be the same."
    @assert res.isMform[] && x.isMform[] && y.isMform[] "Input scalars and output polynomial should be both in Montgomery form."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for j = 1:x.len
        if eval.ismod[j]
            tmp = Mmul(x.vals[j], y.vals[j], eval.moduli[j])
            if res.isntt[]
                for i = 1:res.N
                    res.coeffs[i, j] = sub(res.coeffs[i, j], tmp, eval.moduli[j])
                end
            else
                res.coeffs[1, j] = sub(res.coeffs[1, j], tmp, eval.moduli[j])
            end
        else
            tmp = Bmul(x.vals[j], y.vals[j], eval.moduli[j])
            res.coeffs[1, j] = sub(res.coeffs[1, j], tmp, eval.moduli[j])
        end
    end
end

@views reduce!(x::ModPoly, eval::PolyEvaluator) = begin
    @inbounds for i = 1:x.len
        reduce!(x.coeffs[:, i], eval.moduli[i])
    end
end

ntt(x::ModPoly, eval::PolyEvaluator) = begin
    res = ModPoly(x.N, x.len)
    ntt_to!(res, x, eval)
end

ntt!(x::ModPoly, eval::PolyEvaluator) = ntt_to!(x, x, eval)

@views function ntt_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator)
    @assert !x.isntt[] "Polynomial is already in NTT form."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert x.len == res.len "the length of input and output polynomials should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    if x.isMform[]
        @inbounds for i = 1:x.len
            !eval.ismod[i] && continue
            ntt_to!(res.coeffs[:, i], x.coeffs[:, i], eval.ntter[i])
        end
    else
        @inbounds for i = 1:x.len
            !eval.ismod[i] && continue
            !x.isMform[] && Mform_to!(res.coeffs[:, i], x.coeffs[:, i], eval.moduli[i])
            ntt!(res.coeffs[:, i], eval.ntter[i])
        end    
    end
    
    res.isMform[] = true
    res.isntt[] = true
end

intt!(x::ModPoly, eval::PolyEvaluator) = intt_to!(x, x, eval)

@views function intt_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator)
    @assert x.isntt[] "Polynomial is already in coefficient form."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert x.len == res.len "The length of input and output polynomials should be the same."
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    if x.isMform[]
        @inbounds for i = 1:x.len
            !eval.ismod[i] && continue
            intt_to!(res.coeffs[:, i], x.coeffs[:, i], eval.ntter[i])
        end
    else
        @inbounds for i = 1:x.len
            !eval.ismod[i] && continue
            !x.isMform[] && Mform_to!(res.coeffs[:, i], x.coeffs[:, i], eval.moduli[i])
            intt!(res.coeffs[:, i], eval.ntter[i])
        end
    end

    res.isMform[] = true
    res.isntt[] = false
end

function automorphism(x::ModPoly, idx::Int64, eval::PolyEvaluator)
    res = deepcopy(x)
    automorphism!(res, idx, eval)
end

function automorphism_to!(res::ModPoly, x::ModPoly, idx::Int64, eval::PolyEvaluator)
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert x.len == res.len "The length of input and output polynomials should be the same."

    copy!(res, x)
    automorphism!(res, idx, eval)
end

@views function automorphism!(x::ModPoly, idx::Int64, eval::PolyEvaluator)
    @assert x.len == length(eval.moduli) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = 1:x.len
        if eval.ismod[i]
            automorphism!(idx, x.coeffs[:, i], x.isntt[], eval.autbuff, eval.ntter[i])
        else
            automorphism!(idx, x.coeffs[:, i], false, eval.autbuff, eval.ntter[i])
            embed_to!(x.coeffs[:, i], x.coeffs[:, i], eval.ntter[i].Q, eval.moduli[i])
        end
    end
end

# There can be more optimisations.
function to_big(x::ModPoly, eval::PolyEvaluator)
    ismod, moduli = eval.ismod, eval.moduli
    @assert x.len == length(moduli) "The lengths of input scalar and moduli do not match."

    N, Qlen = x.N, x.len

    res = zeros(BigInt, N)
    Qbig = prod([big(Qi.Q) for Qi = moduli])

    @inbounds for j = 1:Qlen
        Qtilde = Qbig ÷ moduli[j].Q
        Qtildeinv = invmod(UInt64(Qtilde % moduli[j].Q), moduli[j])
        if ismod[j] && x.isMform[]
            for i = 1:N
                res[i] += Bmul(iMform(x.coeffs[i, j], moduli[j]), Qtildeinv, moduli[j]) * Qtilde
            end
        else
            for i = 1:N
                res[i] += Bmul(x.coeffs[i, j], Qtildeinv, moduli[j]) * Qtilde
            end
        end
    end

    @. res %= Qbig
    for i = 1:N
        if res[i] > (Qbig >> 1)
            res[i] -= Qbig
        end
    end

    res
end