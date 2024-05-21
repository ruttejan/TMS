using Pkg

using LinearAlgebra
using Combinatorics
using Test



function set_undef_to_inf(edges)
    undef_indices = findall.(x -> x == UndefInitializer(), edges)
    for i in eachindex(edges)
        edges[i][undef_indices[i]] .= Inf
    end
    return edges
end

function get_games(C::Matrix)
    n = size(C, 1)

    vector_C = [C[:, i] for i in 1:n]
    tmp = [i for i in 1:n]

    indices_collection = collect(powerset(tmp, 0,n))
    indices_complement = [setdiff(tmp, x) for x in indices_collection]
    vectors_collection = collect(powerset(vector_C, 0, n))

    m = length(indices_collection)
    games = Dict{Array, Float64}()
    games[indices_collection[1]] = 0
    for i in 2:m

        tmp = 0
        # first part of the game equation (sum of all edges)

        partial_matrix = C[indices_collection[i], indices_collection[i]]
        # here we can substitute undefined values with 0
        partial_matrix[findall(x -> x == UndefInitializer(), partial_matrix)] .= 0 
        tmp += sum(partial_matrix)


        # second part of the game equation (minimum from edges)

        # get edges that are coming into the group
        edges = [vectors_collection[i][j][indices_complement[i]] for j in 1:length(vectors_collection[i])]  
        if !all(x -> isempty(x), edges)
            edges = set_undef_to_inf(edges) # set the undef values to infinity
            min = minimum.(edges)
            # if some vector had only undef values, its minimum is infinity, so we set it to 0 (there was no minimum)
            min[findall(x -> x == Inf, min)] .= 0 
            tmp += sum(min)    
        end
        
        games[indices_collection[i]] = tmp
    end

    return games, n
end


function get_other_indices(indices, i)
    new_indices = copy(indices)
    deleteat!(new_indices, findall(x -> i in x, indices))
    return new_indices
end

# calculates the shapley value for obtained games
function calculate_shapley(games::Dict, n::Int)    
    shapley = zeros(n)

    indices = collect(keys(games))

    for i in 1:n
        tmp = 0
        new_indices = get_other_indices(indices, i)
        for A in new_indices
            l = length(A)
            set_with_i = vcat(A, i)
            sort!(set_with_i)
            marginal_contribution = games[set_with_i] - games[A]
            pA = (factorial(l) * factorial(n - l - 1)) / factorial(n)
            tmp += pA * marginal_contribution
        end

        shapley[i] = tmp
    end

    return shapley
end

# function calculate_banzhaf(games::Dict, n::Int)
#     shapley = zeros(n)

#     indices = collect(keys(games))

#     pA = 1/(2^(n - 1))

#     for i in 1:n
#         tmp = 0
#         new_indices = get_other_indices(indices, i)
#         for A in new_indices
#             l = length(A)
#             set_with_i = vcat(A, i)
#             sort!(set_with_i)
#             marginal_contribution = games[set_with_i] - games[A]
#             tmp += pA * marginal_contribution
#         end

#         shapley[i] = tmp
#     end

#     return shapley
# end


function shapetrust(C::Matrix)
    games, n = get_games(C)
    return calculate_shapley(games, n)
end

