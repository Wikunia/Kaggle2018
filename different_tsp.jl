
using CSV, DataFrames, Distances, DelimitedFiles, JuMP, Cbc

mutable struct Cities
    xy :: Array{Float64,2}
    nprimes :: Vector{Float64}
end

function calc_score(cities, list_path, tenth)
    @views xy_cities = cities.xy
    len_path     = length(list_path)
    # Calc Distance
    @views xy_path   = xy_cities[list_path,:]
    @views @inbounds dist_path = sqrt.(sum((xy_path[1:end-1,:] .- xy_path[2:end,:]).^2; dims=2))[:,1]
    
    # List of Primes 0 to (len_path-1)
    # Flag array, is path's from-city number non-prime?
    @views is_path_from_non_prime   = cities.nprimes[list_path][1:end-1]   
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

function get_cycle(x_val, N)
    # find cycle
    cycle_idx = Int[]
    push!(cycle_idx, 1)
    while true && length(cycle_idx) <= N+1
        v, idx = findmax(x_val[cycle_idx[end],1:N])
        v < 0.9 && println("v: ", v)
        if idx == cycle_idx[1]
            break
        else
            push!(cycle_idx,idx)
        end
    end
    return cycle_idx
end

function solved(m,x,N)
    x_val = getvalue(x)
    
    cycle_idx = get_cycle(x_val, N)
#     println("cycle_idx: ", cycle_idx)
#     println("Length: ", length(cycle_idx))
    if length(cycle_idx) < N
        @constraint(m, sum(x[cycle_idx,cycle_idx]) <= length(cycle_idx)-1)
        return false
    end
    return true
end 

function prime_solved(start, m, x, N, cities, tenth, dists, score, extra, subm_path)
    x_val = getvalue(x)
    cycle_idx = get_cycle(x_val, N)
#     println("New tour: ", cycle_idx)
    new_path = subm_path[start.+cycle_idx]
    new_score, new_extra, min_extra = calc_score(cities, new_path, tenth)
    println("Current New score: ", new_score)
    println("Current New score without extra: ", new_score-new_extra)
    if new_score < score
        return true
    else
        len_cyc = length(cycle_idx)
        @constraint(m, sum(x[cycle_idx[k],cycle_idx[k+1]] for k in 1:len_cyc-1) <= len_cyc-2)
        return false
    end
end

function run_mip(start, cities, tenth, subm_path, N)
    println("Start: ", start)
    old_subm_path = copy(subm_path)
    score, extra, min_extra = calc_score(cities, subm_path[start+1:start+N], tenth)
    println("Current distance with extra: ", score)
    println("Current distance w/o extra: ", score-extra)
    println("Current extra: ", extra)
    println("Search form normal dist < ", score-min_extra)
    
    pre_time = time()
    # generate distance matrix for the N cities
    c_pos = cities.xy[subm_path[start+1:start+N],:]
    dists = zeros(N,N)
    for i=1:N
        for j=i+1:N
            dists[i,j] = euclidean(c_pos[i,:],c_pos[j,:])
            dists[j,i] = dists[i,j]
        end
    end
    # it should be circular so N -> 1 has zero costs
    dists[N,1] = 0
    
    solver = CbcSolver()
    m = Model(solver=solver)
    @variable(m, x[1:N,1:N], Bin)
    @objective(m, Min, sum(x[i,j]*dists[i,j] for i=1:N,j=1:N))
    @constraint(m, notself[i=1:N], x[i,i] == 0)
    @constraint(m, oneout[i=1:N], sum(x[i,1:N]) == 1)
    @constraint(m, onein[j=1:N], sum(x[1:N,j]) == 1)
    for f=1:N, t=1:N
        @constraint(m, x[f,t]+x[t,f] <= 1)
    end
    # N has to be connected to 1
    
    @constraint(m, x[N,1] == 1)

    @constraint(m, sum(x[k,k+1] for k in 1:N-1) <= N-2)
    @constraint(m, sum(x[i,i+1]*dists[i,i+1] for i=1:N-1) <= score-min_extra)
    
#     @constraint(m, sum(x[i,j]*dists[i,j] for i=1:N,j=1:N) >= score-extra+0.001)
    
    
    t = time()
    status = solve(m)
    println("Status: ", status)
    
    if status == :Infeasible
        return subm_path, false
    end
    
    while !solved(m,x,N)
        status = solve(m)
    end
    
    if status == :Infeasible
        return subm_path, false
    end
    
    println("First solution with MIP: ", getobjectivevalue(m))
    
    improved = false
    counter = 0
    removed_cycles = 0
    while !prime_solved(start, m, x, N, cities, tenth, dists, score, extra, subm_path) && counter < 10
        status = solve(m)
        status == :Infeasible && break
        while !solved(m,x,N)
            status = solve(m)
            status == :Infeasible && break
            removed_cycles += 1
            if removed_cycles >= 70
                status == :Infeasible
                break
            end
        end
        status == :Infeasible && break
        counter += 1 
    end
    if counter < 10 && status != :Infeasible
        println("FOUND BETTER...") 
        println("Counter: , ", counter)
        x_val = getvalue(x)
        cycle_idx = get_cycle(x_val, N)
        println("New tour: ", cycle_idx)
        new_path = subm_path[start.+cycle_idx]
        new_score, new_extra, _ = calc_score(cities, new_path, tenth)
        @assert score-new_score > 0
        println("Improved by: ", score-new_score)
        global_tenth = [(i % 10) == 0 for i in 1:length(subm_path)-1]
        println("oldScore: ", calc_score(cities, subm_path, global_tenth))
        subm_path[start+1:start+N] = new_path 
        println("newScore: ", calc_score(cities, subm_path, global_tenth))
        improved = true
    else
       @assert subm_path == old_subm_path     
    end
    println("Time for solving: ", time()-t)

    if status != :Infeasible
        println("Obj: ", getobjectivevalue(m))
    end

    flush(stdout)
    return subm_path, improved # subm_path, improved
end           

function main(from, to)
    cities_csv = CSV.read("cities_p.csv");
    subm_df = CSV.read("tsp_improved_diff_mip.csv");
    xy_cities   = zeros(size(cities_csv)[1],2)
    xy_cities[:,1] = cities_csv[:X]
    xy_cities[:,2] = cities_csv[:Y]
    cities = Cities(xy_cities, cities_csv[:nprimes])

    subm_path = collect(skipmissing(subm_df[:Path]));
    subm_path .+= 1
    N = 102
    for i=from:75:to
        tenth = [(s % 10) == 0 for s in i+1:i+N-1]
        subm_path, improved = run_mip(i, cities, tenth, subm_path, N);
        if improved
           df = DataFrame(Path=subm_path.-1)
           CSV.write("tsp_improved_diff_mip_new.csv", df);
        end
    end
end

main(700,700)
