
using Pkg

using LinearAlgebra
using Test

"""
eigentrust_matmul(C)

Computes the eigenvector which coordinates sums to 1. 
Used only as a check for the dominant eigenvector. 


"""
function eigentrust_matmul(C::Matrix)

    m, n = size(C)
    # check if x is a square matrix
    if m != n
        error("Dimension of given matrix are wrong. (n ≠ m)")
    end

    A = C' - I
    A = vcat(A, ones(m)')
    b = zeros(n + 1)
    b[end] = 1

    v = A \ b
    return v
end


"""
EigenTrust implementation without pre-trusted peers.
It iterates until the difference is less than ϵ.

Can be run with a given initial trust vector or without it. (Then it uses a unitary distribution vector)
"""
function eigentrust_iteration(C::Matrix, trust_v::Vector, ϵ = 1e-10 :: Real, max_iter = 200 :: Integer)

    # check epsilon
    if ϵ < 0
        error("Epsilon is a negative number.")
    end

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


"""
Basic EigenTrust - implementation with pre-trusted peers

C : trust matrix
pretrusted : vector of node numbers of pre-trusted peers
α : coeficient of how much weight we want to give the pre-trusted peers
"""
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
