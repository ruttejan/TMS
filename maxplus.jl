using LinearAlgebra

# =================
# Some essential MAX-PLUS operations

function maxplus_mul(C::Adjoint, v::Vector)
    n = size(C)[1]
    v_new = zeros(n)

    for i in 1:n
        v_new[i] = maxplus_mul(C[i, :], v)
    end

    return v_new
end

function maxplus_mul(C::Matrix, v::Vector)
    n = size(C)[1]
    v_new = zeros(n)

    for i in 1:n
        v_new[i] = maxplus_mul(C[i, :], v)
    end

    return v_new
end

function maxplus_mul(A::Matrix, B::Matrix)
    n, m = size(A)
    o, p = size(B)

    C = zeros(n, p)

    if m != o
        error("Dimensions of matrices doesn't match for multiplication: ", m ,"̸=", o)
    end

    for i in 1:n
        for j in 1:p
            C[i, j] = maxplus_mul(A[i, :],  B[:, j])
        end
    end
    return C
end

function maxplus_mul(v1::Vector, v2::Vector)
    if length(v1) != length(v2)
        error("Vectors must have same lengt")
    end
    n = length(v1)
    v_new = zeros(n)
    for i in 1:n
        v_new[i] = maxplus_mul(v1[i], v2[i])
    end
    output = maximum(v_new)
    return output
end

function maxplus_mul(v::Vector, a::Real)
    n = length(v)
    v_new = zeros(n)
    for i in 1:n
        v_new[i] = maxplus_mul(v[i], a)
    end
    
    return v_new
end

function maxplus_mul(a::Real, b::Real)
    if a == -Inf && b == -Inf 
        return -Inf
    elseif  a == -Inf && b == Inf || a == Inf && b == -Inf
        return - Inf
    else 
        return a + b
    end
end

function maxplus_pow(C::Matrix, p::Integer)
    if p <= 0
        error("Power has to be non-negative Integer.")
    end

    n,m = size(C)

    if n != m
        error("Matrix is not a square matrix")
    end

    if p == 0 
        A = Matrix(I, n, m)
        A[findall(x -> x == 0, A)] .= -Inf
        return A
    end

    for i in 1:p-1
        C = maxplus_mul(C, C)
    end

    return C
end

function maxplus_tr(C::Matrix)
    return maximum(diag(C))
end

function calculate_eigenvalue(C::Matrix)
    n, m = size(C)
    if n != m
        error("The input matrix is not square matrix")
    end

    λ = -Inf
    for i in 1:n

        tmp = maxplus_tr(maxplus_pow(C, i)) / i
        if tmp > λ
            λ = tmp
        end
    end
    return λ
end