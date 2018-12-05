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

"""
    calc_score_with_extra(cities::Cities, list_path, tenth)

Calculate the santa score given a Cities object the path (or a subpath) and the correponding array with a one/true every tenth step.

Attention: tenth might need to be shifted for a subpath. Assumes that list_path has cities 1 indexed so not as in a normal submission file
"""
function calc_score(cities, list_path, tenth)
    @views xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    @views xy_path   = xy_cities[list_path,:]
    @views @inbounds dist_path = sqrt.(sum((xy_path[1:end-1,:] .- xy_path[2:end,:]).^2; dims=2))[:,1]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views is_path_from_non_prime   = cities.nprime[list_path][1:end-1]   
    # If both flags are true, *1.1, else * 1.0
    return sum(dist_path .* (1.0 .+ 0.1 .* is_path_from_non_prime .* tenth))
end

include("mip.jl")

end