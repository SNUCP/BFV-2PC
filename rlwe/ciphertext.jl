"""
RLWE is a struct for RLWE ciphertext.
"""
struct RLWE
    b::ModPoly
    a::ModPoly

    function RLWE(b::ModPoly, a::ModPoly)
        @assert b.N == a.N && b.len == a.len && b.isMform[] == a.isMform[] && b.isntt[] == a.isntt[] "The mask and body should have the same parameters."
        new(b, a)
    end

    function RLWE(b::ModPoly)
        a = ModPoly(b.N, b.len, b.isntt[], b.isMform[])
        new(b, a)
    end

    function RLWE(N::Int64, len::Int64, isntt::Bool=true, isMform::Bool=true)
        new(ModPoly(N, len, isntt, isMform), ModPoly(N, len, isntt, isMform))
    end
end

Base.copy!(dst::RLWE, src::RLWE) = begin
    @assert dst.b.len == src.b.len "The length of input and output ciphertexts should match."
    copy!(dst.b, src.b)
    copy!(dst.a, src.a)
end

"""
RLEV is a struct for RLEV ciphertext.
"""
struct RLEV
    glen::Int64
    stack::Vector{RLWE}

    function RLEV(stack::Vector{RLWE})
        new(length(stack), stack)
    end

    function RLEV(N::Int64, len::Int64, glen::Int64, isntt::Bool=true, isMform::Bool=true)
        new(glen, [RLWE(N, len, isntt, isMform) for _ = 1:glen])
    end
end

Base.copy!(dst::RLEV, src::RLEV) = begin
    @assert dst.glen == src.glen "The length of input and output ciphertexts should match."
    copy!(dst.stack, src.stack)
end

"""
RGSW is a struct for RGSW ciphertext.
"""
struct RGSW
    basketb::RLEV
    basketa::RLEV

    function RGSW(basketb::RLEV, basketa::RLEV)
        new(basketb, basketa)
    end

    function RGSW(N::Int64, len::Int64, glen::Int64, isntt::Bool=true, isMform::Bool=true)
        new(RLEV(N, len, glen, isntt, isMform), RLEV(N, len, glen, isntt, isMform))
    end
end

Base.copy!(dst::RGSW, src::RGSW) = begin
    copy!(dst.basketb, src.basketb)
    copy!(dst.basketa, src.basketa)
end