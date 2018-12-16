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
    mod = ((mod_from-1) % 10)+1
    if mod == 11
        mod = 1
    end 
    return precomputed[mod][from] - precomputed[mod][to]
end

"""
    improved_score_by(cities::Cities, tour::Vector, precomputed::Vector{Vector{Float64}},
        rev_precomputed::Vector{Vector{Float64}}, swap::Vector, euc_dict::Dict{Tuple{Int,Int},Float64}, base_score::Float64)

Return -1 if not improved (maybe on early break) otherwise return by how much it improved compared to base_score
"""
function improved_score_by(cities::Cities, tour::Vector, precomputed::Vector{Vector{Float64}},
    rev_precomputed::Vector{Vector{Float64}}, swap::Vector, euc_dict::Dict{Tuple{Int,Int},Float64}, base_score::Float64)
    @views cities_xy = cities.xy
    # global euc_time, subpath_time

    cfrom = swap[1]
    score = 0.0
    ltour = length(tour)
    lswap = length(swap)
    # half length same as Int(floor(lswap/2)) but faster
    hlswap = fld(lswap,2)
    cfroms = Vector{Int}(undef, hlswap)
    cfroms[1] = cfrom

    # println("swap: ", swap)
    for i in 1:hlswap-1
        # s_time = time_ns()
        i2 = 2*i
        i2m1 = i2-1
        i2p1 = i2+1

        # for the swap edge each time
        cfrom += 1

        # tours between swap edges
        if swap[i2] < swap[i2p1]
            t1 = (swap[i2]-cfrom+10000000) % 10 + 1
            new_scoring = calc_score(precomputed, t1, swap[i2], swap[i2p1])
            cfrom += swap[i2p1]-swap[i2]
            score += new_scoring
        else
            t1 = ((1+ltour-swap[i2])-cfrom+10000000) % 10 + 1
            new_scoring = calc_score(rev_precomputed, t1,
                                    1+ltour-swap[i2], 1+ltour-swap[i2p1])
            score += new_scoring
            cfrom += swap[i2]-swap[i2p1]
        end
        cfroms[i+1] = cfrom
        # subpath_time += time_ns()-s_time

        if score > base_score
            return -1
        end
    end

    for i=1:hlswap
        # e_time = time_ns()
        i2 = 2*i
        i2m1 = i2-1
        i2p1 = i2+1
        
        # swap edge
        dist = 0.0
        tpl = (swap[i2m1],swap[i2])
        if haskey(euc_dict,tpl)
            dist = euc_dict[tpl]
        else
            dist = euclidean(cities_xy[tour[swap[i2m1]],:],cities_xy[tour[swap[i2]],:])
            euc_dict[tpl] = dist
        end
        if cfroms[i] % 10 == 0 && cities.nprime[tour[swap[i2m1]]] == 1.0
            dist *= 1.1
        end
        
        score += dist
        # euc_time += time_ns()-e_time

        if score > base_score
            return -1
        end
    end
    return base_score - score
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