module TRP

using DataFrames, CSV

mutable struct Cities
    xy :: Array{Float64,2}
    prime :: Vector{Bool}
    nprime :: Vector{Float64}
end

"""
    read_cities_and_subm(cities_fn, submission_fn)

Get the cities submission file as a cities object 
"""
function read_cities_and_subm(cities_fn, submission_fn)
    cities_csv = CSV.read(cities_fn);
    subm_df = CSV.read(submission_fn);
    xy_cities   = zeros(size(cities_csv)[1],2)
    xy_cities[:,1] = cities_csv[:X]
    xy_cities[:,2] = cities_csv[:Y]
    cities = Cities(xy_cities, cities_csv[:prime], cities_csv[:nprime])

    subm_path = collect(skipmissing(subm_df[:Path]));
    subm_path .+= 1
    return cities, subm_path
end

"""
    calc_score_with_extra(cities::Cities, list_path, tenth)

Calculate the santa score given a Cities object the path (or a subpath) and the correponding array with a one/true every tenth step.
Additonally the current extra costs are returned and an approximation of the minimum amount of extra costs possible.
Use the `calc_score` if you are not interested in the extra information

Attention: tenth might need to be shifted for a subpath. Assumes that list_path has cities 1 indexed so not as in a normal submission file
"""
function calc_score_with_extra(cities::Cities, list_path, tenth)
    @views xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    @views xy_path   = xy_cities[list_path,:]
    @views @inbounds dist_path = sqrt.(sum((xy_path[1:end-1,:] .- xy_path[2:end,:]).^2; dims=2))[:,1]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views is_path_from_non_prime   = cities.nprime[list_path][1:end-1]   
    # If both flags are true, *1.1, else * 1.0
    extra = sum(dist_path .* 0.1 .* is_path_from_non_prime .* tenth)
    min_extra = 0
    sum_tenth = convert(Int,sum(tenth))
    nof_primes = convert(Int,len_path-1-sum(is_path_from_non_prime))
    if nof_primes < sum_tenth
        k_idx = partialsortperm(dist_path, 1:sum_tenth-nof_primes)
        min_extra = sum(0.1 .* dist_path[k_idx])
    end
    return sum(dist_path)+extra, extra, min_extra
end

function calc_score_dists_and_extras(cities::Cities, list_path, tenth)
    @views xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    xy_path   = xy_cities[list_path,:]
    @inbounds dist_path = sqrt.(sum((xy_path[1:end-1,:] .- xy_path[2:end,:]).^2; dims=2))[:,1]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views is_path_from_non_prime   = cities.nprime[list_path][1:end-1]   
    # If both flags are true, *1.1, else * 1.0
    extra = dist_path .* 0.1 .* is_path_from_non_prime .* tenth
    return dist_path, extra
end

function calc_score_extras(cities::Cities, list_path, tenth, dist_path)
    @views xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    @views xy_path   = xy_cities[list_path,:]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views is_path_from_non_prime   = cities.nprime[list_path][1:end-1]   
    # If both flags are true, *1.1, else * 1.0
    extra = dist_path .* 0.1 .* is_path_from_non_prime .* tenth
    return extra
end

"""
    calc_score(cities::Cities, list_path::Vector, tenth::Vector)

Calculate the santa score given a Cities object the path (or a subpath) and the correponding array with a one/true every tenth step.

Attention: tenth might need to be shifted for a subpath. Assumes that list_path has cities 1 indexed so not as in a normal submission file
"""
function calc_score(cities::Cities, list_path::Vector, tenth::Vector)
    @views xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    @views xy_path   = xy_cities[list_path,:]
    @views @inbounds dist_path = sqrt.(sum((xy_path[1:end-1,:] .- xy_path[2:end,:]).^2; dims=2))[:,1]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views is_path_from_non_prime   = cities.nprime[list_path[1:end-1]]
    # If both flags are true, *1.1, else * 1.0
    return sum(dist_path .* (1.0 .+ 0.1 .* is_path_from_non_prime .* tenth))
end

function calc_score(precomputed::Vector{Vector{Float64}}, mod_from::Int, from::Int, to::Int)
    return precomputed[mod_from][from] - precomputed[mod_from][to]
end

"""
    improved_by!(cfroms::Vector{Int}, cities::Cities, tour::Vector{Int}, precomputed::Vector{Vector{Float64}},
        rev_precomputed::Vector{Vector{Float64}}, swap::Vector{Int}, euc_dict::Dict{Tuple{Int,Int},Float64}, base_score::Float64,
        ltour::Int,ltour1::Int,lswap::Int,hlswap::Int
        )

This function is highly optimized therefore a lot of parameters are used.
Return the improved of a swap compared to the base score. Return -1 if no improvement was achieved (early break)
"""
function improved_by!(cfroms::Vector{Int}, cities_xy::Array{Float64,2}, cities_nprime::Vector{Float64}, tour::Vector{Int}, precomputed::Vector{Vector{Float64}},
        rev_precomputed::Vector{Vector{Float64}}, swap::Vector{Int}, euc_dict::Dict{Tuple{Int,Int},Float64}, base_score::Float64,
        ltour::Int,ltour1::Int,lswap::Int,hlswap::Int
    )
    score = 0.0
    cfrom = swap[1]

    # calculate the long subpaths first as it's faster to calculate
    for i in 1:hlswap-1
        cfroms[i] = cfrom 
        i2 = 2*i
        i2p1 = i2+1
            
        cfrom += 1

        # tours between swap edges
        if swap[i2] < swap[i2p1]
            # + 10000000 because it can be negative
            t1 = (swap[i2]-cfrom+10000000) % 10 + 1
            cfrom += swap[i2p1]-swap[i2]
            score += calc_score(precomputed, t1, swap[i2], swap[i2p1])
        else
            t1 = ((ltour1-swap[i2])-cfrom+10000000) % 10 + 1
            score += calc_score(rev_precomputed, t1,  ltour1-swap[i2], ltour1-swap[i2p1])
            cfrom += swap[i2]-swap[i2p1]
        end
    end
    if score > base_score
        return -1.0
    end
    cfroms[hlswap] = cfrom 

    # calculate the single swap edges if the path can still be better
    for i in 1:hlswap
        i2 = 2*i
        i2m1 = i2-1
        # swap edge
        dist = get(euc_dict, (swap[i2m1],swap[i2]), 0.0)
        if dist == 0.0
            dist = euclidean(cities_xy[tour[swap[i2m1]],:],cities_xy[tour[swap[i2]],:])
            euc_dict[(swap[i2m1],swap[i2])] = dist
        end
        if cfroms[i] % 10 == 0 && cities_nprime[tour[swap[i2m1]]] == 1.0
            dist *= 1.1
        end
            
        score += dist
        if score > base_score
            return -1.0
        end
    end
    return base_score-score
end

"""
	rev_cumsum(l::Vector)

Compute the reverse cumulative sum of the Vector
[1,1,3,5,2] -> [12,11,10,7,2]
"""
function rev_cumsum(l::Vector)
    q = copy(l)
	@inbounds for i = length(l)-1:-1:1
        q[i] += q[i+1]
	end
	return q
end

include("mip.jl")
include("kopt.jl")

end