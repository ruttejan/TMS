# using Pkg

using LinearAlgebra
include("normal_form.jl")
include("maxpower.jl")
include("maxplus.jl")
using Test



# Creates start and end positions of each block matrix in the normal form
function prepare_matrices_dimensions(classes::Vector, n::Int)
    starts = zeros(Int, length(classes))
    ends = zeros(Int, length(classes))
    starts[1] = 1
    ends[length(classes)] = n

    for i in 1 : length(classes) - 1
        cur_start = starts[i]
        class_len = length(classes[i])

        ends[i] = cur_start + class_len - 1
        starts[i + 1] = ends[i] + 1
    end
    return starts, ends
end


function maxtrust(C::Matrix, w::Vector, T::Real)

    # variables set-up
    n = size(C, 1) 
    lambdas = zeros(n) # λ vector (eigenvalues)
    v = zeros(Float64, n) # eigenvector 
    e = zeros(n) # ξ vector

    # D ≣ Cᵀ
    D, perm, inessential_classes, essential_classes = get_normal_form(copy(C'))


    classes = vcat(inessential_classes, essential_classes)

    
    starts, ends = prepare_matrices_dimensions(classes, n)
    
    num_of_classes = length(classes)

    # matrix C is irreducible
    if num_of_classes == 1
        lambda, v = max_power(C, w)
        lambdas = zeros(n) .+ lambda
        return v, lambdas
    end

    Dnn = D[starts[num_of_classes] : ends[num_of_classes], starts[num_of_classes] : ends[num_of_classes]]
    λn = 0
    vn = zeros(length(classes[num_of_classes]))

    # check if the diagonal block is a matrix or scalar
    dim1, dim2 = size(Dnn)
    if dim1 == 1 && dim1 == dim2
        λn = Dnn[1, 1]
    else
        wn = w[starts[num_of_classes] : ends[num_of_classes]]
        λn, vn = max_power_not_transposed(Dnn, wn)
    end


    lambdas[starts[num_of_classes] : ends[num_of_classes]] .= λn
    e[starts[num_of_classes] : ends[num_of_classes]] .= λn
    v[starts[num_of_classes] : ends[num_of_classes]] = vn

    # go through all diagonal block from last to first
    for j in (num_of_classes - 1) : -1 : 1
      
        
        Djj = D[starts[j] : ends[j], starts[j] : ends[j]]
        dim1, dim2 = size(Djj)
        
        # check if the diagonal block is a matrix or scalar
        if dim1 == dim2 && dim1 == 1
            λj = Djj[1, 1]
        else
            wj = w[starts[j] : ends[j]]
            λj, _ = max_power_not_transposed(Djj, wj)
        end
        
        lambdas[starts[j] : ends[j]] .= λj
        vj = zeros(length(classes[j])) .+ -Inf

        # BIG-O-PLUS (Maximum) on lines 10 and 13
        for k in j:num_of_classes
            Djk = D[starts[j] : ends[j], starts[k] : ends[k]]
            wk = w[starts[k] : ends[k]]
            # vk = v[starts[k] : ends[k]]
            if λj == -Inf && j == 1
                curr_vj = maxplus_mul(Djk, wk)
                # curr_vj = maxplus_mul(Djk, vk)
            else
                curr_vj = maxplus_mul(Djk, wk) .+ (λj * (j - 1))
                # curr_vj = maxplus_mul(maxplus_mul(Djk, vk) , (λj * (j - 1)))
            end
            vj = max.(curr_vj, vj)
        end

        # finished line 10
        if λj > e[starts[j + 1]]
            e[starts[j] : ends[j]] .= λj
        else # finished line 13
            e[starts[j] : ends[j]] .= lambdas[starts[j + 1]]
            ej = e[starts[j]]
            # vj = vj .+ ej ^ (-1)                      # Which one is correct?
            vj = maxplus_mul(vj, -ej)
        end

        v[starts[j] : ends[j]] = vj
    end

    t = v + e .* T

    # returning the vector to the initial permutation (1,2,3,....,n)
    t_final = zeros(n)
    e_final = zeros(n)

    for i in 1:n
        t_final[perm[i]] = t[i]
        e_final[perm[i]] = e[i]
    end

    return t_final, e_final
    
end

# Cᵀ ⊗ v = λ ⊗ v
# (Cᵀ ⊗ v) - v = λ


# Irreducible example

# A = [[-Inf 7/11 4/11]
# [1 -Inf -Inf]
# 0 0 -Inf]

# v = maxtrust(A, [0, 0, 0], 1)

# maxplus_mul(copy(A'), v) .- v



e = -Inf

m = [[e e 1 e]
    [0.4 e 0.6 0]    #11
    [1 e e e]
    [e 0.6 0.4 e]]

v, h = maxtrust(m, [0.25, 0.25, 0.25, 0.25], 1)

maxplus_mul(copy(m'),  v)

h + v




# A = [[e 0.8 0.2 e e]
#      [0.7 e 0.3 e e]
#      [e e e 0 0]
#      [e e e e 1]
#      [e e e 1 e]]

# v = maxtrust(A, [0, 0, 0, 0, 0], 1)

# maxplus_mul(copy(A'), v) .- v
e = -Inf

C = [[e e e e e e e e e e]
    [0.318182 e 0.272727 0.409091 e e e e e e]
    [e e e 1.0 e e e e e e]
    [e 0.533333 0.466667 e e e e e e e]
    [e e 0.0 e e 1.0 0.0 e e e]
    [e e e 0.0 1.0 e e e 0.0 e]
    [e e e e e e e 0.363636 0.363636 0.272727]
    [e e e e e e e e e 1.0]
    [e e e e e e 1.0 e e e]
    [e e e e e e e e 1.0 e]]

# display(get_normal_form(copy(C'))[1])

get_normal_form(copy(C'))[3:4]

A = get_normal_form(copy(C'))[1]

w = zeros(10)
v, h = maxtrust(C, w, 1)

left = maxplus_mul(A, v)

right = h + v

left - right

# maxplus_mul(copy(C'), v) .- v

# max_power_not_transposed([[e 1] #[5, 6]
#                         [1 e]])

# max_power_not_transposed([-Inf        -Inf     0.533333 
#                         0.272727  -Inf     0.466667 # [2, 3, 4]
#                         0.409091    1.0  -Inf])

# max_power_not_transposed([[-Inf        -Inf     1.0  -Inf ]
#                         [0.363636  -Inf   -Inf   -Inf ]     # 7, 8, 9, 10
#                         [0.363636  -Inf   -Inf     1.0 ]
#                         [0.272727    1.0  -Inf   -Inf ] ] )   
                        
