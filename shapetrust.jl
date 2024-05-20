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
    # indices_complement = reverse(indices_collection)
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

function calculate_banzhaf(games::Dict, n::Int)
    shapley = zeros(n)

    indices = collect(keys(games))

    pA = 1/(2^(n - 1))

    for i in 1:n
        tmp = 0
        new_indices = get_other_indices(indices, i)
        for A in new_indices
            l = length(A)
            set_with_i = vcat(A, i)
            sort!(set_with_i)
            marginal_contribution = games[set_with_i] - games[A]
            tmp += pA * marginal_contribution
        end

        shapley[i] = tmp
    end

    return shapley
end


function shapetrust(C::Matrix)
    games, n = get_games(C)
    return calculate_shapley(games, n)
end


A = [[undef 7/11 4/11]
[1 undef undef]
0 0 undef]

v = shapetrust(A)

# get_games(A)

# sum(v)



### TESTS
# @testset "tests" begin
#     @testset "get_games tests" begin
#         @test begin
#             C = [
#                 [undef undef 1 1]
#                 [1 undef undef undef]
#                 [undef 1 undef 1]
#                 [1 1 undef undef]
#             ]

#             n = 4
#             tmp = [i for i in 1:n]
#             indices_collection = collect(powerset(tmp, 0,n))

#             games = [0, 1, 1, 1, 1, 3, 2, 4, 3, 3, 3, 5, 6, 5, 5, 7]
#             result = get_games(C)[1]

#             value = true

#             for i in eachindex(indices_collection)
#                 if !isapprox(result[indices_collection[i]], games[i])
#                     value = false
#                 end
#             end
            
#             value

#         end

#         @test begin
#             C = [
#                 [undef 0.2 0.3]
#                 [0.1 undef 0.4]
#                 [0.5 0.7 undef]
#             ]
#             n = 3
#             tmp = [i for i in 1:n]
#             indices_collection = collect(powerset(tmp, 0,n))

#             games = [0, 0.1, 0.2, 0.3, 1.5, 1.3, 1.6, 2.2]
#             result = get_games(C)[1]
            

#             value = true

#             for i in eachindex(indices_collection)
#                 if !isapprox(result[indices_collection[i]], games[i])
#                     value = false
#                 end
#             end
            
#             value
#         end
#     end


#     @testset "calculate_shapley test" begin
#         @test begin

#             C = [
#                 [undef  undef 1 1]
#                 [1 undef undef undef]
#                 [undef 1 undef 1]
#                 [1 1 undef undef]
#             ]

#             games, n = get_games(C)

#             shapley = [1.833333333, 1.833333333, 1.333333333, 2]
#             result = calculate_shapley(games, n)

#             isapprox(shapley, result)
#         end

#         @test begin


#             C = [
#                 [undef 0.2 0.3]
#                 [0.1 undef 0.4]
#                 [0.5 0.7 undef]
#             ]

#             games, n = get_games(C)

#             shapley = [0.61666666, 0.81666666, 0.76666666]
#             result = calculate_shapley(games, n)

#             isapprox(shapley, result)
#         end
#     end

#     @testset "shapetrust test" begin
#         @test begin

#             C = [
#                 [undef undef 1 1]
#                 [1 undef undef undef]
#                 [undef 1 undef 1]
#                 [1 1 undef undef]
#             ]

#             shapley = [1.833333333, 1.833333333, 1.333333333, 2]
#             result = shapetrust(C)

#             isapprox(shapley, result)
#         end

#         @test begin

#             C = [
#                 [undef 0.2 0.3]
#                 [0.1 undef 0.4]
#                 [0.5 0.7 undef]
#             ]

#             shapley = [0.61666666, 0.81666666, 0.76666666]
#             result = shapetrust(C)

#             isapprox(shapley, result)
#         end
#     end
# end


# @testset "examples 11-15" begin

#     #11
#     @test begin

#         C = [
#             [undef 0.3 0.1]
#             [0.2 undef 0.7]
#             [undef undef undef]
#         ]

#         shapley = [0.316666666, 0.4166666666, 0.56666666666]
#         result = shapetrust(C)

#         isapprox(shapley, result)
#     end

#     #12
#     @test begin

#         C = [
#             [undef 0.2 undef]
#             [0.8 undef 0.9]
#             [undef 0.9 undef]
#         ]

#         shapley = [0.95, 0.8, 1.05]
#         result = shapetrust(C)

#         isapprox(shapley, result)
#     end

#     #13
#     @test begin

#         C = [
#             [undef 0.1 undef]
#             [0.9 undef 0.1]
#             [1 undef undef]
#         ]

#         shapley = [1.5666666666, 0.26666666666, 0.26666666666]
#         result = shapetrust(C)

#         isapprox(shapley, result)
#     end

#     #14
#     @test begin

#         C = [
#             [undef undef 0.9]
#             [0.1 undef 0.1]
#             [0.9 undef undef]
#         ]

#         shapley = [0.85, 0.3, 0.85]
#         result = shapetrust(C)

#         isapprox(shapley, result)
#     end

#     #15
#     @test begin

#         C = [
#             [undef 0.9 0.9]
#             [0.1 undef 1]
#             [undef 0.1 undef]
#         ]

#         shapley = [0.4166666666, 0.8666666666, 1.71666666666]
#         result = shapetrust(C)

#         isapprox(shapley, result)
#     end
# end



