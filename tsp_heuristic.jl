
using CSV, DataFrames, Distances, DelimitedFiles

vcattime = 0.0
starttime = 0.0
sumtime = 0.0

mutable struct Score
    value :: Float64
    dist_path :: Vector{Float64}
    scaled_dist_path ::  Vector{Float64}
end

function Base.copy(a::Score) 
    return Score(a.value, a.dist_path, a.scaled_dist_path)
end

mutable struct Cities
    xy :: Array{Float64,2}
    nprimes :: Vector{Float64}
end

function get_score(cities, subm_path)
    global distdict
    all_ids = cities[:CityId]
    all_x = cities[:X]
    all_y = cities[:Y]

    score = 0.0
    p1 = Vector{Int}()
    p2 = Vector{Int}()
    for i in 1:length(subm_path)-1
        c_idx = subm_path[i]+1
        n_idx = subm_path[i+1]+1
        p1 = [all_x[c_idx],all_y[c_idx]]
        p2 = [all_x[n_idx],all_y[n_idx]]
        stepSize = euclidean(p1,p2)
        if i % 10 == 0 && !cities[:primes][subm_path[i]+1]
            stepSize *= 1.1
        end
#         println(stepSize)
        score += stepSize
    end
    return score
end

function calc_score(cities, list_path, tenth)
    xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    xy_path   = xy_cities[list_path,:]
    dist_path = sqrt.(sum((xy_path[1:end-1,:] .- xy_path[2:end,:]).^2; dims=2))[:,1]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    is_path_from_non_prime   = cities.nprimes[list_path][1:end-1]   
    # If both flags are true, *1.1, else * 1.0
    result = dist_path .* (1.0 .+ 0.1 .* is_path_from_non_prime .* tenth)
    return Score(sum(result), dist_path, result)
end

function calc_score_reverse!(dist_path, result, xy_cities, nprimes, list_path, tenth, normal_dist_path, reversed_dist_path, low, high)
    len_path  = length(list_path)
    # Calc Distance
    xy_path   = xy_cities[list_path,:]
    
    dist_path[1:low-2] = normal_dist_path[1:low-2]
    dist_path[low-1] = euclidean(xy_path[low-1,:],xy_path[low,:])
    @views dist_path[low:high-1] = reversed_dist_path[len_path-high+1:len_path-low]
    dist_path[high] = euclidean(xy_path[high,:],xy_path[high+1,:])
    dist_path[high+1:end] = normal_dist_path[high+1:end]
    #=
    dist_path = vcat(normal_dist_path[1:low-2], # before switch
                     [euclidean(xy_path[low-1,:],xy_path[low,:])], # new edge at the beginning
                     normal_dist_path[low:high-1][end:-1:1], # reverse this part
                     [euclidean(xy_path[high,:],xy_path[high+1,:])], # new edge in the end
                     normal_dist_path[high+1:end],
                )[:,1]
    =#
    
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views @inbounds is_path_from_non_prime  = nprimes[list_path[1:end-1]]
    # If both flags are true, *1.1, else * 1.0
    result[1:end] = dist_path 
    result .*= (1.0 .+ 0.1 .* is_path_from_non_prime .* tenth)
#     new_result = @. dist_path * (1.0 + 0.1 * is_path_from_non_prime * tenth)
    sum_result = sum(result)
    return Score(sum_result, dist_path, result)
end

function two_opt_reasonable!(old, new, xy_path, dist_path, low, high; margin=100)
    @inbounds old = dist_path[low-1]+dist_path[high]
#     @assert euclidean(xy_path[low-1,:],xy_path[low,:])+euclidean(xy_path[high,:],xy_path[high+1,:]) == old
    @inbounds new = euclidean(xy_path[low-1,:],xy_path[high,:])+euclidean(xy_path[low,:],xy_path[high+1,:])
    return new-margin < old
end

function two_opt(cities,subm_path)
    n = length(subm_path)
    path = copy(subm_path)
    println("Path: ", path[1:10])
    println("n: ", n)
    switchLow = 2
	switchHigh = n - 1
	reverses = 0
	improved = true
    tenth = [(i % 10) == 0 for i in 1:n-1]
    oldCost = calc_score(cities, path, tenth)
    reversed_dist_path = reverse(oldCost.dist_path)
    time_for_10000 = time()
    overwrite_dist_path = zeros(n-1)
    overwrite_result = zeros(n-1)
    checked_reasonable = 0
    oldReasonable, newReasonable = 0.0, 0.0
    xy_cities = cities.xy
    xy_path   = xy_cities[path,:]
    nprimes   = cities.nprimes
    
    t = time()
    print_every = 100
	while improved
		improved = false
		# we can't change the first 
		for i in switchLow:(switchHigh-1)
            if i % print_every == 0
                println("i: ", i)
                println("checked_reasonable: ", checked_reasonable)
                println("Time for one i: ", (time()-t)/print_every)
                flush(stdout)
                t = time()
            end
            
			for j in switchHigh:-1:(i+1)
                reasonable =  two_opt_reasonable!(oldReasonable, newReasonable,
                                                      xy_path, oldCost.dist_path, i, j)
                if reasonable
                    checked_reasonable += 1
                    n_path = copy(path)
                    reverse!(n_path, i, j)
                    altCost = calc_score_reverse!(overwrite_dist_path, overwrite_result,
                                                 xy_cities, nprimes, n_path, tenth, oldCost.dist_path, 
                                                 reversed_dist_path, i, j)

                   

                    if altCost.value < oldCost.value
                        altCost_dif =  calc_score(cities, n_path, tenth)
                        println(altCost.value - altCost_dif.value)
                        println("i: ", i)
                        println("j: ", j)
                        println("altCost.value: ", altCost.value)
                        println("altCost_dif.value: ", altCost_dif.value)
                        @assert altCost.value - altCost_dif.value == 0
                        
                        println("Improved by: ", oldCost.value-altCost.value)
                        path = n_path
                        oldCost = altCost_dif
                        reversed_dist_path = reverse(oldCost.dist_path)
                        xy_path   = xy_cities[path,:]
                        reverses += 1
                        improved = true
                        df = DataFrame(Path=path.-1)
                        CSV.write("tsp_improved_new.csv", df);
                    end
                end
			end
            
		end
	end
	println("Reverses 1: ", reverses)
    score = calc_score(cities, path, tenth)
	return path, score
    
end

function main()
    global vcattime, starttime, sumtime
    vcattime = 0.0
    starttime = 0.0
    sumtime = 0.0
    cities_csv = CSV.read("cities_p.csv");
    subm_df = CSV.read("tsp_improved.csv"); # 53
    subm_path = collect(skipmissing(subm_df[:Path]))
    subm_path .+= 1
    xy_cities   = zeros(size(cities_csv)[1],2)
    xy_cities[:,1] = cities_csv[:X]
    xy_cities[:,2] = cities_csv[:Y]
    cities = Cities(xy_cities, cities_csv[:nprimes])
    
    tenth = [(i % 10) == 0 for i in 1:length(subm_path)-1]
#     subm_path = vcat(subm_path[1:200], subm_path[end-200:end])
    @time new_score = calc_score(cities, subm_path, tenth)
    println("Score new: ", new_score.value)
    
    path, cost = two_opt(cities,subm_path)
    println("New cost: ", cost.value)
    return path
end
new_path = main();
println("vcattime: ", vcattime)
println("starttime: ", starttime)
println("sumtime: ", sumtime)
# println(new_path)
