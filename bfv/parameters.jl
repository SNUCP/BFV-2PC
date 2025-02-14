"""
BFVParamSketch is a struct for a very vague parameter setting.
"""
struct BFVParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64
    ptxt_modulus::Int64
    ispacking::Bool
    islevelled::Bool

    BFVParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, ptxt_modulus::Int64; dlen::Int64=0, ispacking::Bool=true, islevelled::Bool=true) =
        new(ring_param, logP, logQ, dlen, ptxt_modulus, ispacking, islevelled)
end

struct BFVParameters
    ring_param::RingParam
    P::Union{Missing,Vector{UInt64}}
    Q::Vector{UInt64}
    dlen::Int64
    ptxt_modulus::UInt64
    ispacking::Bool
    islevelled::Bool

    BFVParameters(ring_param::RingParam, P::Union{Missing, Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, ptxt_modulus::Int64, ispacking::Bool, islevelled::Bool) = 
        new(ring_param, P, Q, dlen, ptxt_modulus, ispacking, islevelled)

    function BFVParameters(sketch::BFVParamSketch)
        ring_param, logP, logQ, dlen, t, ispacking, islevelled = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen, sketch.ptxt_modulus, sketch.ispacking, sketch.islevelled
    
        if logP == 0
            P = missing
    
            Qlen = ceil(Int64, logQ / 62)
            Qbits = round(Int64, logQ / Qlen)
            Q0bit = logQ - Qbits * (Qlen - 1)
    
            while Q0bit > 62
                Q0bit -= Qbits
                Qlen += 1
            end
    
            Qprimes = find_prime(ring_param, Qbits, Qlen + 1)
            Q0 = find_prime(ring_param, Q0bit)[1]
            filter!(x -> x ≠ Q0, Qprimes)
            Q = vcat(Q0, Qprimes[1:Qlen-1])
    
            if dlen == 0
                dlen = 1
            end
        else
            Plen = ceil(Int64, logP / 62)
            Pbits = logP / Plen
            P = find_prime(ring_param, Pbits, Plen)
    
            maxbits = min(62, logP)
            Qlen = ceil(Int64, logQ / maxbits)
            Qbits = round(Int64, logQ / Qlen)
            Q0bit = logQ - Qbits * (Qlen - 1)
    
            while Q0bit > maxbits
                Q0bit -= Qbits
                Qlen += 1
            end
    
            Qprimes = find_prime(ring_param, Qbits, Plen + Qlen + 1)
            Q0 = find_prime(ring_param, Q0bit)[1]
            filter!(x -> x ∉ P && x ≠ Q0, Qprimes)
            Q = vcat(Q0, Qprimes[1:Qlen-1])
    
            if dlen == 0
                dlen = floor(Int64, logP / max(Q0bit, Qbits))
            end
        end
    
        new(ring_param, P, Q, dlen, t, ispacking, islevelled)
    end
end