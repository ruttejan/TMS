
using LinearAlgebra
include("maxplus.jl")

# Max-Power function using classic algebra multiplication
# fails!!!
function max_power_normal_mult(C::Matrix, r::Vector)
    n = length(r)
    p = 1
    V = r
    vp = copy(r)
    c = -1
    q = 0

    while true

        vp_new = C * vp
        p += 1
        for (i, x) in enumerate(eachcol(V))
            tmp = vp_new .- x
            if all(y -> isapprox(y, tmp[1]), tmp) && tmp[1] >= 0
                c = tmp[1]
                q = i
                break
            end
        end

        vp = vp_new
        V = hcat(V, vp_new)

        if c >= 0
            break
        end
    end

    λ = c / (p - q)

    max_v = repeat([-Inf], n)

    for i in 1:(p-q)
        tmp = ((p - q - i) * λ) .+ V[:, q + i - 1]
        max_v = max.(max_v, tmp)
    end
    
    return λ, max_v
end

function max_power_normal_mult(C::Matrix)
    n = size(C)[1]
    r = ones(n) ./ n
    return max_power_normal_mult(C, r)
end




# =====================




# Max-Power function using Max-Plus algebra multiplication

function max_power(C::Matrix, r::Vector)
    n = length(r)
    p = 1
    V = r
    vp = copy(r)
    c = -1
    q = 0

    while true

        vp_new = maxplus_mul(copy(C'), vp)
        p += 1
        for (i, x) in enumerate(eachcol(V))
            tmp = vp_new .- x
            if all(y -> isapprox(y, tmp[1]), tmp) && tmp[1] >= 0
                c = tmp[1]
                q = i
                break
            end
        end

        vp = vp_new
        V = hcat(V, vp_new)

        if c >= 0
            break
        end
    end

    λ = c / (p - q)

    max_v = repeat([-Inf], n)

    for i in 1:(p-q+1)
        tmp = ((p - q - i) * λ) .+ V[:, q + i - 1]
        max_v = max.(max_v, tmp)
    end
    
    return λ, max_v
end

function max_power(C::Matrix)
    n = size(C)[1]
    r = zeros(n)
    return max_power(C, r)
end

# Power Method for a matrix - not its transposition
function max_power_not_transposed(C::Matrix, r::Vector)
    n = length(r)
    p = 1
    V = r
    vp = copy(r)
    c = -1
    q = 0

    while true

        vp_new = maxplus_mul(C, vp)
        p += 1
        for (i, x) in enumerate(eachcol(V))
            tmp = vp_new .- x
            if all(y -> isapprox(y, tmp[1]), tmp) && tmp[1] >= 0
                c = tmp[1]
                q = i
                break
            end
        end

        vp = vp_new
        V = hcat(V, vp_new)

        if c >= 0
            break
        end
    end

    λ = c / (p - q)

    max_v = repeat([-Inf], n)

    for i in 1:(p-q+1)
        tmp = ((p - q - i) * λ) .+ V[:, q + i - 1]
        max_v = max.(max_v, tmp)
    end
    
    return λ, max_v
end

function max_power_not_transposed(C::Matrix)
    n = size(C)[1]
    r = zeros(n)
    return max_power(C, r)
end


