
using Pkg

using LinearAlgebra
using Test

function eigentrust_matmul(C::Matrix)
    
    # check if rows sums to 1
    row_sums = sum(C, dims=2)
    # if !all(x -> isapprox(x, 1), row_sums)
    #     error("Rows does not add to one")
    # end

    m, n = size(C)
    # check if x is a square matrix
    if m != n
        error("Dimension of given matrix are wrong. (n ≠ m)")
    end

    A = C' - I
    A = vcat(A, ones(m)')
    b = zeros(n + 1)
    b[end] = 1

    # sA = size(A)
    # sB = size(b)
    # println("Dims of A: $sA; Dims of b: $sB")

    v = A \ b
    return v
end

function eigentrust_iteration(C::Matrix, trust_v::Vector, ϵ = 1e-10 :: Real, max_iter = 200 :: Integer)

    # check epsilon
    if ϵ < 0
        error("Epsilon is a negative number.")
    end

    # check if rows sums to 1
    # row_sums = sum(C, dims=2)
    # if !all(x -> isapprox(x, 1), row_sums)
    #     error("Rows does not add to one")
    # end

    m, n = size(C)
    # check if x is a square matrix
    if m != n
        error("Dimension of given matrix are wrong. (n ≠ m)")
    end
   
    δ = Inf
    for i in 1:max_iter
        trust_v_new = C' * trust_v
        δ = norm(trust_v_new - trust_v)
        trust_v = trust_v_new
        if δ < ϵ
            break
        end
    end
    
    return trust_v

end

function eigentrust_iteration(C::Matrix, ϵ = 1e-10 :: Real, max_iter = 200 :: Integer)
    m, n = size(C)
    trust_v = ones(m) ./ m
    eigentrust_iteration(C, trust_v, ϵ, max_iter)
end

function basic_eigentrust(C::Matrix, pretrusted::Vector, α::Real, ϵ = 1e-10 :: Real)

    # check epsilon
    if ϵ < 0
        error("Epsilon is a negative number.")
    end

    num_of_pretrusted = length(pretrusted)

    if num_of_pretrusted == 0 || α == 0.0
        return eigentrust_iteration(C)
    end

    C = normalize_trust_matrix(C, pretrusted, α)

    
    m, n = size(C)
    ρ = zeros(n)

    ρ[pretrusted] .= 1/num_of_pretrusted



    # check if x is a square matrix
    if m != n
        error("Dimension of given matrix are wrong. (n ≠ m)")
    end

    trust_v = ones(m) ./ m
    δ = Inf
    while δ > ϵ
        trust_v_new = C' * trust_v
        trust_v_new = (1 - α) .* trust_v_new + α .* ρ
        δ = norm(trust_v_new - trust_v)
        trust_v = trust_v_new
    end
    
    return trust_v

end

eigentrust_iteration([[0 7/11 4/11]
[1 0 0]
[0 0 0]], [0.4437410968480305,
0.35398293836943506,
0.20227596478253435])

eigentrust_matmul([[0 7/11 4/11]
[1 0 0]
[0 0 0]])

# @testset "Eigentrust and Eigentrust_matmul returns same solution:" begin
#     @test  begin              
#         n = 5
#         C = rand(n, n)
#         C ./= sum(C, dims=2)

#         sol1 = eigentrust_matmul(C)
#         sol2 = eigentrust_iteration(C)

#         isapprox(sol1, sol2)
#     end

#     @test  begin              
#         n = 6
#         C = rand(n, n)
#         C ./= sum(C, dims=2)

#         sol1 = eigentrust_matmul(C)
#         sol2 = eigentrust_iteration(C)

#         isapprox(sol1, sol2)
#     end

#     @test  begin              
#         n = 10
#         C = rand(n, n)
#         C ./= sum(C, dims=2)

#         sol1 = eigentrust_matmul(C)
#         sol2 = eigentrust_iteration(C)

#         isapprox(sol1, sol2)
#     end

#     @test  begin              
#         n = 1000
#         C = rand(n, n)
#         C ./= sum(C, dims=2)

#         sol1 = eigentrust_matmul(C)
#         sol2 = eigentrust_iteration(C)

#         isapprox(sol1, sol2)
#     end

#     @test begin
#         C = [[0.1 0.9 0.0]
#              [0.2 0.8 0.0]
#              [0.3 0.7 0.0]]

#         sol1 = eigentrust_matmul(C)
#         sol2 = eigentrust_iteration(C)

#         isapprox(sol1, sol2)
#     end

#     @test begin
#         C = [[0.0 1.0 0.0]
#              [0.0 0.0 1.0]
#              [1.0 0.0 0.0]]
             
#         sol1 = eigentrust_matmul(C)
#         sol2 = eigentrust_iteration(C)

#         isapprox(sol1, sol2)
#     end
# end

A = [[0 7/11 4/11]
[1 0 0]
[0 0 0]]