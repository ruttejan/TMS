using LinearAlgebra
# ======
# This file provides the EigenTrust NORMALIZATION of a trust matrix
# and convertion of a reducible matrix to its NORMAL FORM
#
# Source form NORMAL FORM is (Seneta, 2006) as in MAXTRUST paper Definition 24
# ======


# =============================
# NORMALIZATION

# Matrix normalization with PRE-TRUSTED peers
function normalize_trust_matrix(C::Matrix, pretrusted::Vector, α::Real)

    indexes = findall(x -> x == -Inf, C)

    if α < 0 || α > 1
        error("Alpha should be between 0 and 1")
    end

    n = size(C)[1]
    num_of_pretrusted = length(pretrusted)
    ρ = zeros(n)

    if num_of_pretrusted != 0
        ρ[pretrusted] .= 1/num_of_pretrusted
    end

    max_applied = max.(C, 0)

    row_sums = sum(max_applied, dims=2)

    final = zeros(n, n)

    for i in 1:n
        if row_sums[i] == 0
            final[i, :] = ρ
        else
            final[i, :] = max_applied[i, :] ./ row_sums[i]
        end
        
    end
    final = final .* (1 - α)
    final = copy((final' .+ (ρ .* α))')

    final[indexes] .= -Inf

    return final
end

# Matrix normalization

function normalize_trust_matrix(C::Matrix)

    indexes = findall(x -> x == -Inf, C)
    n = size(C)[1]

    max_applied = max.(C, 0)
    row_sums = sum(max_applied, dims=2)

    final = zeros(n, n)

    for i in 1:n
        if row_sums[i] == 0
            final[i, :] .= 0
        else
            final[i, :] = max_applied[i, :] ./ row_sums[i]
        end
    end
    final[indexes] .= -Inf
    return final
end

# ===============

# NORMAL FORM

# creates a dictionary of reachable nodes (node_number -> reachable nodes)
function get_indices_scopes(communication_matrix::Matrix)
    n = size(communication_matrix, 1)

    can_reach_dict = Dict()

    visited = falses(n)

    # dfs function for finding all reachable indices from an index
    function dfs(node, can_reach)
        visited[node] = true
        push!(can_reach, node)
        
        for i in 1:n
            if communication_matrix[node, i] >= 0 && !visited[i]
                dfs(i, can_reach)
            end
        end
    end

    # find all reachable indices for every index
    # creates dictionary: can_reach_dict
    for i in 1:n
        visited = falses(n)
        can_reach = []
        visited[i] = true
        for j in 1:n
            if communication_matrix[i, j] >= 0 && !visited[j]
                dfs(j, can_reach)
            end
        end

        if isempty(can_reach) && communication_matrix[i, i] >= 0
            push!(can_reach, i)
        end

        can_reach_dict[i] = can_reach
    end

    return can_reach_dict
end

# separates nodes into two groups - essential and inessential
function separate_indices(communication_matrix::Matrix, can_reach_dict::Dict)
    
    n = size(communication_matrix, 1)

    essential_indices = []
    inessential_indices = []

    for i in 1:n
        essential = true

        if isempty(can_reach_dict[i])
            essential = false
        end

        for x in can_reach_dict[i]
            if !(i in can_reach_dict[x])
                essential = false
                break
            end
        end

        if essential
            push!(essential_indices, i)
        else
            push!(inessential_indices, i)
        end
    end

    return essential_indices, inessential_indices
end

# removing used nodes from a dictionary
# creates a new dictionary with filtered nodes
function filter_dict(indices::Vector, can_reach_dict::Dict)
    copy_of_dict = copy(can_reach_dict)
    for x in indices
        delete!(copy_of_dict, x)
    end

    for key in keys(copy_of_dict)
        setdiff!(copy_of_dict[key], indices)
    end

    return copy_of_dict
end

# removing used nodes from a dictionary
# filters the input dictionary
function filter_dict_inplace(indices::Vector, can_reach_dict::Dict)

    for x in indices
        delete!(can_reach_dict, x)
    end

    for key in keys(can_reach_dict)
        setdiff!(can_reach_dict[key], collect(Iterators.Flatten(indices)))
    end
end

# creates Essential classes (graph components) from essential nodes
function get_essential_classes(can_reach_dict::Dict)
    essential_classes = []
    while true
        dict_keys = collect(keys(can_reach_dict))

        if isempty(dict_keys)
            return essential_classes
        end
        idx = dict_keys[1]
        class = can_reach_dict[idx]
        
        if class[1] == idx
            push!(essential_classes, class)
        else
            push!(class, idx)
            push!(essential_classes, class)
        end
        filter_dict_inplace(class, can_reach_dict)
        
    end

    return essential_classes
end

# creates Inessential classes (graph components) from inessential nodes
function get_inessential_classes(can_reach_dict::Dict)
    inessential_classes = []
    while true
        dict_keys = collect(keys(can_reach_dict))

        if isempty(dict_keys)
            return inessential_classes
        end

        idx = dict_keys[1]
        can_reach = can_reach_dict[idx]
        
        if isempty(can_reach)
            filter_dict_inplace([idx], can_reach_dict)
            push!(inessential_classes, [idx])
        else
            class = []
            for x in can_reach
                if idx in can_reach_dict[x]
                    push!(class, x)
                end
            end
            push!(class, idx)
            filter_dict_inplace(class, can_reach_dict)
            push!(inessential_classes, class)
        end
        
    end

    return inessential_classes
end

# creates a dictionary of reduced graph (class -> class)
function get_reduced_dict(classes, from_dict)
    reduced_dict = Dict()

    for class in classes
        reachable_from_indices = []
        for x in class
            reachable_from_indices = vcat(reachable_from_indices, from_dict[x])
        end
        setdiff!(reachable_from_indices, class)
        reduced_dict[class] = reachable_from_indices
    end

    return reduced_dict
end

# sorts inessential classes
function sort_inessential_classes(inessential_classes, reduced_dict)
    sorted_classes = []

    while true
        dict_keys = collect(keys(reduced_dict))

        if isempty(dict_keys)
            return sorted_classes
        end

        currently_pushed = []

        for key in dict_keys
            if isempty(reduced_dict[key])
                push!(sorted_classes, key)
                push!(currently_pushed, key)
            end
        end
        
        filter_dict_inplace(currently_pushed, reduced_dict)

    end
    return sorted_classes
end

# gets inessential and essential classes
function get_classes(communication_matrix::Matrix)
    can_reach_dict = get_indices_scopes(communication_matrix)


    essential_indices, inessential_indices = separate_indices(communication_matrix, can_reach_dict)

    essential_dict = filter_dict(inessential_indices, can_reach_dict)
    inessential_dict = filter_dict(essential_indices, can_reach_dict)

    
    # essential classes can be sorted easily, because there are no classes, that can be reached from them
    essential_classes = sort.(get_essential_classes(essential_dict))

    inessential_classes = sort.(get_inessential_classes(copy(inessential_dict)))

    is_reachable_from_dict = get_indices_scopes(copy(communication_matrix'))
    inessential_from_dict = filter_dict(essential_indices, is_reachable_from_dict)
    reduced_inessential_from_dict = get_reduced_dict(inessential_classes, inessential_from_dict)
    inessential_classes_sorted = sort_inessential_classes(inessential_classes, reduced_inessential_from_dict)


    return inessential_classes_sorted, essential_classes
end

# creates the normal form through row-column permutation
# return normal form, permutation vector and all classes (essential and inessential)
function get_normal_form(communication_matrix::Matrix)
    inessential_classes, essential_classes  = get_classes(communication_matrix)
    first_part = collect(Iterators.Flatten(inessential_classes))
    second_part = collect(Iterators.Flatten(essential_classes))
    perm = vcat(first_part, second_part)
    P = I(length(perm))[:, perm]
    return transpose(P) * communication_matrix * P, perm, inessential_classes, essential_classes
end


# Example usage:

# communication_matrix = [1 1 0 0 0 0 0 0 0;
#                         1 1 1 0 0 0 1 0 0;
#                         0 0 0 0 0 0 1 0 0;
#                         0 0 0 1 0 0 0 0 1;
#                         0 0 0 0 1 0 0 0 0;
#                         0 0 1 0 0 1 0 0 0;
#                         0 0 1 0 0 0 0 0 0;
#                         0 1 0 0 0 1 0 1 0;
#                         0 0 0 1 0 0 0 0 1]

# normal_form, perm, inessential_classes, essential_classes = get_normal_form(communication_matrix)







