"""
Decomposer is a struct for the gadget decomposition. It decomposes the polynomial using RNS decomposition, from bottom to top, with digit length `dlen`.
"""
mutable struct Decomposer
    const moduli::Vector{Modulus}
    const be::Vector{BasisExtender}
    const gvec::Vector{ModScalar}
    const dlen::Int64
    const glen::Int64

    function Decomposer(eval::PolyEvaluator, dlen::Int64)
        ismod, moduli = eval.ismod, eval.moduli
        lenQ = length(moduli)
        glen = ceil(Int64, lenQ / dlen)
        N = eval.ntter[1].N

        be = Vector{BasisExtender}(undef, glen)
        gvec = Vector{ModScalar}(undef, glen)
        @views @inbounds for i = 1:glen
            lo, hi = 1 + (i - 1) * dlen, min(i * dlen, lenQ)
            be[i] = BasisExtender(moduli[lo:hi], moduli, N)

            tildeQ, Di = big(1), big(1)
            @simd for j = 1:lenQ
                if lo ≤ j ≤ hi
                    Di *= moduli[j].Q
                else
                    tildeQ *= moduli[j].Q
                end
            end
            tildeQ *= invmod(tildeQ, Di)
            gveci = Vector{UInt64}(undef, lenQ)
            for j = 1:lenQ
                if ismod[j]
                    gveci[j] = Mform((tildeQ % moduli[j].Q) % UInt64, moduli[j])
                else
                    gveci[j] = (tildeQ % moduli[j].Q) % UInt64
                end
            end

            gvec[i] = ModScalar(gveci, true)
        end

        new(collect(moduli), be, gvec, dlen, glen)
    end
end

function decomp(x::ModPoly, decer::Decomposer)
    res = [ModPoly(x.N, x.len) for _ = 1:decer.glen]
    decomposeto!(res, x, decer)
    res
end

@views function decomposeto!(res::AbstractVector{ModPoly}, x::ModPoly, decer::Decomposer)
    lenQ, dlen, glen = length(decer.moduli), decer.dlen, decer.glen

    @assert !x.isMform[] && !x.isntt[] "The input polynomial should not be in Montgomery or evaluation form."
    @assert lenQ == x.len "The decomposition parameter does not match the input polynomial."
    @assert glen == length(res) "The decomposition parameter does not match the output polynomials."

    @inbounds for i = eachindex(res)
        idx1, idx2 = 1 + (i - 1) * dlen, min(i * dlen, lenQ)
        basis_extend!(res[i].coeffs, x.coeffs[:, idx1:idx2], decer.be[i])
        res[i].isMform[] = false
        res[i].isntt[] = false
    end
end
