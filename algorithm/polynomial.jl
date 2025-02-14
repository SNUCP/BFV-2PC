# Baby step.
function polynomial_evaluate(coeffs::Vector{UInt64}, x::HECiphertext, rlk::RLEV, oper::HEOperator)
    deg = length(coeffs) - 1
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    halfdeg = 1 << (maxlevel - 1)
    first_half = Vector{typeof(x)}(undef, halfdeg)

    # Compute the monomials for the first half.
    first_half[1] = x
    @inbounds for i = 1:maxlevel-2
        # Compute the power-of-two monomial.
        first_half[1<<i] = mul(first_half[1<<(i-1)], first_half[1<<(i-1)], rlk, oper)
        for j = 1:(1<<i)-1
            # Checks if the polynomial is sparse.
            checkj = coeffs[(1<<i)+j+1] == 0
            for idx = i+1:maxlevel-2
                checkj = checkj && coeffs[(1<<idx)+(1<<i)+j+1] == 0
            end
            if j + halfdeg ≤ deg
                checkj = checkj && coeffs[j+halfdeg+1] == 0
            end
            checkj && continue

            # Compute the monomial if needed in the future.
            first_half[(1<<i)+j] = mul(first_half[j], first_half[1<<i], rlk, oper)
        end
    end
    first_half[1<<(maxlevel-1)] = mul(first_half[1<<(maxlevel-2)], first_half[1<<(maxlevel-2)], rlk, oper)

    _poly_eval_from_monomial(coeffs, first_half, x, rlk, oper)
end

function polynomial_evaluate_many(coeffs::Vector{Vector{UInt64}}, x::HECiphertext, rlk::RLEV, oper::HEOperator)
    maxdeg = maximum([length(coeffsi) for coeffsi = coeffs]) - 1
    maxlevel = trailing_zeros(Base._nextpow2(maxdeg + 1))
    halfdeg = 1 << (maxlevel - 1)
    first_half = Vector{typeof(x)}(undef, halfdeg)

    # Compute the monomials for the first half.
    first_half[1] = x
    @inbounds for i = 1:maxlevel-2
        # Compute the power-of-two monomial.
        first_half[1<<i] = mul(first_half[1<<(i-1)], first_half[1<<(i-1)], rlk, oper)
        for j = 1:(1<<i)-1
            # Checks if the polynomial is sparse.
            checkj = true
            for idx1 = eachindex(coeffs)
                (1 << i) + j + 1 > length(coeffs[idx1]) && continue
                checkj = checkj && coeffs[idx1][(1<<i)+j+1] == 0
                for idx2 = i+1:maxlevel-2
                    (1 << idx2) + (1 << i) + j + 1 > length(coeffs[idx1]) && continue
                    checkj = checkj && coeffs[idx1][(1<<idx2)+(1<<i)+j+1] == 0
                end
                if j + halfdeg ≤ deg
                    checkj = checkj && coeffs[idx1][j+halfdeg+1] == 0
                end
            end
            checkj && continue

            # Compute the monomial if needed in the future.
            first_half[(1<<i)+j] = mul(first_half[j], first_half[1<<i], rlk, oper)
        end
    end
    first_half[1<<(maxlevel-1)] = mul(first_half[1<<(maxlevel-2)], first_half[1<<(maxlevel-2)], rlk, oper)

    res = Vector{typeof(x)}(undef, length(coeffs))
    for i = eachindex(coeffs)
        deg = length(coeffs[i]) - 1
        maxleveli = trailing_zeros(Base._nextpow2(deg + 1))
        halfdegi = 1 << (maxleveli - 1)

        @views res[i] = _poly_eval_from_monomial(coeffs[i], first_half[1:halfdegi], x, rlk, oper)
    end

    res
end

function _poly_eval_from_monomial(coeffs::Vector{UInt64}, first_half::AbstractVector{T}, x::T, rlk::RLEV, oper::HEOperator) where {T<:HECiphertext}
    deg = length(coeffs) - 1
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    halfdeg = 1 << (maxlevel - 1)
    @assert length(first_half) == halfdeg "first_half must have length 2^(maxlevel-1)."

    res = similar(x)
    tmp = similar(x)

    # Compute the monomials for the last half and sum, using the binary tree.
    halfhalfdeg = halfdeg >> 1
    @inbounds for i = 1:deg-halfdeg
        coeffs[i+halfdeg+1] == 0 && continue

        tmpi = tmp[1:length(tmp.val)]

        if i ≤ halfhalfdeg
            mul_to!(tmpi, coeffs[i+halfdeg+1], first_half[i], oper)
        else
            checki = false
            for j = 0:maxlevel-2
                if (i >> j) & 1 == 1
                    if !checki
                        mul_to!(tmpi, first_half[1<<j], coeffs[1+i+halfdeg], oper)
                        checki = true
                    else
                        mul_to!(tmpi, tmpi, first_half[1<<j], rlk, oper)
                    end
                end
            end
        end

        add_to!(res, tmpi, res, oper)
    end

    # Multiply x^halfdeg for a better computation complexity.
    mul_to!(res, res, first_half[halfdeg], rlk, oper)

    # Compute the linear sum.
    @inbounds for i = 1:halfdeg
        coeffs[i+1] == 0 && continue

        mul_to!(tmp, coeffs[i+1], first_half[i], oper)
        add_to!(res, tmp, res, oper)
    end

    add_to!(res, coeffs[1], res, oper)

    res
end