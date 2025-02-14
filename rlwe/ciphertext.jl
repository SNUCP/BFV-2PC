"""
RLWE is a struct for RLWE ciphertext.
"""
struct RLWE
    b::ModPoly
    a::ModPoly
    isPQ::RefBool
    auxQ::RefInt

    function RLWE(b::ModPoly, a::ModPoly; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        @assert b.N == a.N && length(b) == length(a) && b.isntt[] == a.isntt[] "The mask and body should have the same parameters."
        new(b, a, Ref(isPQ), Ref(auxQ))
    end

    function RLWE(b::ModPoly; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        a = similar(b)
        new(b, a, Ref(isPQ), Ref(auxQ))
    end

    function RLWE(N::Int64, len::Int64; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(ModPoly(N, len, isntt=isntt), ModPoly(N, len, isntt=isntt), Ref(isPQ), Ref(auxQ))
    end
end

Base.copy(x::RLWE) = RLWE(copy(x.b), copy(x.a), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.copy!(dst::RLWE, src::RLWE) = begin
    @assert length(dst.b) == length(src.b) "The length of input and output ciphertexts should match."
    copy!(dst.b, src.b)
    copy!(dst.a, src.a)
    dst.isPQ[] = src.isPQ[]
    dst.auxQ[] = src.auxQ[]
end

Base.:length(x::RLWE) = begin
    @assert length(x.a) == length(x.b) "Something is wrong with the ciphertext."
    length(x.b)
end

Base.:getindex(x::RLWE, idx::UnitRange) = RLWE(x.b[idx], x.a[idx], isPQ=x.isPQ[], auxQ=length(x) ∈ idx ? x.auxQ[] : UInt64(0))
Base.:similar(x::RLWE) = RLWE(similar(x.b), similar(x.a), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:resize!(x::RLWE, len::Int64) = begin
    resize!(x.b, len)
    resize!(x.a, len)
end

initialise!(x::RLWE; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    initialise!(x.b, isntt=isntt)
    initialise!(x.a, isntt=isntt)

    x.isPQ[] = isPQ
    x.auxQ[] = auxQ
end

ntt(ct::RLWE, eval::PolyEvaluator) =
    RLWE(ntt(ct.b, eval), ntt(ct.a, eval))

ntt!(ct::RLWE, eval::PolyEvaluator) = ntt_to!(ct, ct, eval)

ntt_to!(res::RLWE, ct::RLWE, eval::PolyEvaluator) = begin
    ntt_to!(res.b, ct.b, eval)
    ntt_to!(res.a, ct.a, eval)
end

intt(ct::RLWE, eval::PolyEvaluator) =
    RLWE(intt(ct.b, eval), intt(ct.a, eval))

intt!(ct::RLWE, eval::PolyEvaluator) = intt_to!(ct, ct, eval)

intt_to!(res::RLWE, ct::RLWE, eval::PolyEvaluator) = begin
    intt_to!(res.b, ct.b, eval)
    intt_to!(res.a, ct.a, eval)
end

#=================================================================================================#

"""
Tensor is a struct for Tensored RLWE ciphertext.
"""
struct Tensor{N}
    vals::NTuple{N,ModPoly}
    isPQ::RefBool
    auxQ::RefInt

    function Tensor{N}(vals::NTuple{N,ModPoly}; isPQ::Bool=false, auxQ::UInt64=UInt64(0)) where {N}
        new{N}(vals, Ref(isPQ), Ref(auxQ))
    end

    function Tensor(vals::NTuple{N,ModPoly}; isPQ::Bool=false, auxQ::UInt64=UInt64(0)) where {N}
        new{N}(vals, Ref(isPQ), Ref(auxQ))
    end

    function Tensor(val::ModPoly, degree::Int64; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        vals = Vector{ModPoly}(undef, degree)

        vals[1] = val
        @inbounds for i = 2:degree
            vals[i] = similar(val)
        end

        new{degree}(Tuple(vals), Ref(isPQ), Ref(auxQ))
    end

    function Tensor(N::Int64, len::Int64, degree::Int64=3; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new{degree}(Tuple(ModPoly(N, len, isntt=isntt) for _ = 1:degree), Ref(isPQ), Ref(auxQ))
    end
end

Base.:getindex(ct::Tensor{N}, idx::UnitRange) where {N} = Tensor{N}(NTuple{N,ModPoly}(ct.vals[i][idx] for i = 1:N), isPQ=ct.isPQ[], auxQ=ct.auxQ[])

degree(ct::Tensor) = length(ct.vals)

Base.:length(ct::Tensor{N}) where {N} = begin
    len = length(ct.vals[1])
    @inbounds for i = 2:N
        @assert len == length(ct.vals[i]) "Something is wrong with the ciphertext."
    end
    len
end

Base.copy(src::Tensor{N}) where {N} = Tensor{N}(Tuple(copy(vali) for vali = src.vals), isPQ=src.isPQ[], auxQ=src.auxQ[])

Base.copy!(dst::Tensor{N}, src::Tensor{N}) where {N} = begin
    for i = 1:N
        copy!(dst.vals[i], src.vals[i])
    end
    dst.isPQ[] = src.isPQ[]
    dst.auxQ[] = src.auxQ[]
end

Base.:similar(x::Tensor) = Tensor([similar(vali) for vali = x.val], isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:resize!(x::Tensor, len::Int64) = begin
    @inbounds for vali = x.vals
        resize!(vali, len)
    end
end

initialise!(x::Tensor{N}; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) where {N} = begin
    @inbounds for i = 1:N
        initialise!(x.vals[i], isntt=isntt)
    end

    x.isPQ[] = isPQ
    x.auxQ[] = auxQ
end

ntt(ct::Tensor{N}, eval::PolyEvaluator) where {N} = Tensor([ntt(ct.vals[i], eval) for i = 1:N], isPQ=ct.isPQ[], auxQ=ct.auxQ[])

ntt!(ct::Tensor{N}, eval::PolyEvaluator) where {N} = ntt_to!(ct, ct, eval)

ntt_to!(res::Tensor{N}, ct::Tensor{N}, eval::PolyEvaluator) where {N} = begin
    @inbounds for i = 1:N
        ntt_to!(res.vals[i], ct.vals[i], eval)
    end
end

intt(ct::Tensor{N}, eval::PolyEvaluator) where {N} = Tensor([intt(ct.vals[i], eval) for i = 1:N], isPQ=ct.isPQ[], auxQ=ct.auxQ[])

intt!(ct::Tensor{N}, eval::PolyEvaluator) where {N} = intt_to!(ct, ct, eval)

intt_to!(res::Tensor{N}, ct::Tensor{N}, eval::PolyEvaluator) where {N} = begin
    @inbounds for i = 1:N
        intt_to!(res.vals[i], ct.vals[i], eval)
    end
end

#=================================================================================================#

"""
RLEV is a struct for RLEV ciphertext.
"""
struct RLEV
    glen::Int64
    stack::Vector{RLWE}

    function RLEV(stack::Vector{RLWE})
        @inbounds for rlwe = stack
            @assert rlwe.auxQ[] == 0 "The auxiliary modulus should be zero for RLEV encryption."
        end
        new(length(stack), stack)
    end

    function RLEV(N::Int64, len::Int64, glen::Int64; isntt::Bool=true, isPQ::Bool=false)
        new(glen, [RLWE(N, len, isntt=isntt, isPQ=isPQ) for _ = 1:glen])
    end
end

Base.copy(x::RLEV) = RLEV(copy(x.stack))

Base.copy!(dst::RLEV, src::RLEV) = begin
    @assert dst.glen == src.glen "The length of input and output ciphertexts should match."
    copy!(dst.stack, src.stack)
end

Base.:similar(ct::RLEV) = RLEV(ct.glen, [similar(stacki) for stacki = ct.stack])

Base.:resize!(ct::RLEV, len::Int64) = begin
    @inbounds for stacki = ct.stack
        resize!(stacki, len)
    end
end

initialise!(x::RLEV; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    @inbounds for i = 1:x.len
        initialise!(x.stack[i], isntt=isntt, isPQ=isPQ, auxQ=auxQ)
    end
end

ntt(ct::RLEV, eval::PolyEvaluator) = RLEV((ntt.(ct.stack, Ref(eval))))

ntt!(ct::RLEV, eval::PolyEvaluator) = ntt_to!(ct, ct, eval)

ntt_to!(res::RLEV, ct::RLEV, eval::PolyEvaluator) = begin
    @assert res.glen == ct.glen "The length of input and output ciphertext should match."

    @inbounds for i = 1:res.len
        ntt_to!(res.stack[i], ct.stack[i], eval)
    end
end

intt(ct::RLEV, eval::PolyEvaluator) = RLEV((intt.(ct.stack, Ref(eval))))

intt!(ct::RLEV, eval::PolyEvaluator) = intt_to!(ct, ct, eval)

intt_to!(res::RLEV, ct::RLEV, eval::PolyEvaluator) = begin
    @assert res.glen == ct.glen "The length of input and output ciphertext should match."
    @inbounds for i = 1:res.len
        intt_to!(res.stack[i], ct.stack[i], eval)
    end
end

#=================================================================================================#

"""
RGSW is a struct for RGSW ciphertext.
"""
struct RGSW
    basketb::RLEV
    basketa::RLEV

    function RGSW(basketb::RLEV, basketa::RLEV)
        @assert basketb.glen == basketa.glen "The length of input and output ciphertexts should match."
        new(basketb, basketa)
    end

    function RGSW(N::Int64, len::Int64, glen::Int64; isntt::Bool=true, isPQ::Bool=false)
        new(RLEV(N, len, glen, isntt=isntt, isPQ=isPQ), RLEV(N, len, glen, isntt=isntt, isPQ=isPQ))
    end
end

Base.copy(x::RGSW) = RGSW(copy(x.basketb), copy(x.basketa))
Base.copy!(dst::RGSW, src::RGSW) = begin
    copy!(dst.basketb, src.basketb)
    copy!(dst.basketa, src.basketa)
end

Base.:similar(ct::RGSW) = RGSW(similar(ct.basketb), similar(ct.basketa))

Base.:resize!(ct::RGSW, len::Int64) = begin
    resize!(ct.basketb, len)
    resize!(ct.basketa, len)
end

initialise!(x::RGSW; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    initialise!(x.basketb, isntt=isntt, isPQ=isPQ, auxQ=auxQ)
    initialise!(x.basketa, isntt=isntt, isPQ=isPQ, auxQ=auxQ)
end

ntt(ct::RGSW, eval::PolyEvaluator) =
    RGSW(ntt(ct.basketb, eval), ntt(ct.basketa, eval))

ntt_to!(res::RGSW, ct::RGSW, eval::PolyEvaluator) = begin
    ntt_to!(res.basketb, ct.basketb, eval)
    ntt_to!(res.basketa, ct.basketa, eval)
end

intt(ct::RGSW, eval::PolyEvaluator) =
    RGSW(intt(ct.basketb, eval), intt(ct.basketa, eval))

intt_to!(res::RGSW, ct::RGSW, eval::PolyEvaluator) = begin
    intt_to!(res.basketb, ct.basketb, eval)
    intt_to!(res.basketa, ct.basketa, eval)
end

const Ciphertext = Union{RLWE,Tensor,RLEV,RGSW}