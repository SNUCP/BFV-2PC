# function mul_BSGS(x::PlainText, M::PlainMatrix)
#     n = x.param.n
#     n₁ = 2^ceil(Int64, log(n)/2)
#     n₂ = n ÷ n₁

#     res = zeroplaintxt(x.param)
#     rotx = [rot(x, i) for i = 0 : n₁ - 1]

#     for j = 0 : n₂ - 1
#         c = zeroplaintxt(x.param)
#         for i = 1 : n₁
#             tmpv = rot(M.ptxts[n₁*j + i], -n₁ * j)
#             c += PlainText(mul_fft(rotx[i].coeffs, tmpv.coeffs, x.param.ζ), x.param)
#         end
#         res += rot(c, n₁ * j)
#     end

#     res / M.Δ
# end

struct PackedMatrix
    val::Vector{Vector{UInt64}}
    n1::Int64
    n2::Int64
end

function mul_BSGS()

end