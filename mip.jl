
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

function get_cycle(m, x, x_vals, N, to_vec, cycle_vars; first=true)
    cycle_idx = Int[]
    if !first
        # get next cycle
        for y=1:N
            if !haskey(cycle_vars,y)
                push!(cycle_idx, y)
                to_vec,bcycle = get_vec(N,to_vec,x_vals;cycle_det=true,i=y)
                break
            end
        end
    else
        push!(cycle_idx, 1)
    end
    first = false
    while true
        idx = to_vec[cycle_idx[end]]
        cycle_vars[idx] = true
        push!(cycle_idx, idx)
        if idx == cycle_idx[1]
            break
        end
    end
    sumx = 0
    last = cycle_idx[1]
    for i=2:length(cycle_idx)
        if last < cycle_idx[i]
            sumx += x[last,cycle_idx[i]]
        else 
            sumx += x[cycle_idx[i],last]
        end
        last = cycle_idx[i]
    end
    return cycle_idx, sumx
end

function get_vec(N,to_vec,x_vals;i=1,used=Dict{Int64,Bool}(),cycle_det=false,bcycle=false)
    idx = 0
    m_idx = 0
    for j = i+1:N
        if isapprox(x_vals[i,j]-1, 0, atol=1e-6)
            m_idx = 0
            if !haskey(used,j)
                idx = j
                used[j] = true
                to_vec[i] = j
                x_vals[i,j] = 0
                to_vec, bcycle = get_vec(N,to_vec,x_vals;i=idx,used=used,cycle_det=cycle_det,bcycle=bcycle)
                break
            end
        end
    end
    if idx == 0
        for j = 1:i-1
            if isapprox(x_vals[j,i]-1, 0, atol=1e-6)
                m_idx = 0
                if !haskey(used,j)
                    idx = j
                    used[j] = true
                    x_vals[i,j] = 0
                    to_vec[i] = j
                    to_vec, bcycle = get_vec(N,to_vec,x_vals;i=idx,used=used,cycle_det=cycle_det,bcycle=bcycle)
                    break
                end
            end
        end
    end
    
    if idx == 0 && cycle_det
        for y=1:N
            if to_vec[y] == 0
                return to_vec, 1
            end
        end
        return to_vec, bcycle
    end

    if idx == 0 
        for y=1:N
            if to_vec[y] == 0
                to_vec, bcycle = get_vec(N,to_vec,x_vals;i=y,used=used,cycle_det=cycle_det,bcycle=bcycle)
                break
            end
        end
    end
    return to_vec, bcycle
end

function get_x_vals(x, N)
    x_sparse = getvalue(x)

    x_vals = zeros(N,N)
    for i=1:N, j=i+1:N
        x_vals[i,j] = x_sparse[i,j]
    end
    return x_vals
end

"""
    solved(m, x, N)

Checks whether there is a subtour if yes returns false and the subtour will be removed
in the next solve
"""
function solved(m, x, N)
    x_vals = get_x_vals(x,N)

    # find cycle
    to_vec = zeros(Int,N)
    used = Dict{Int64,Bool}()
    objective = 0
    cycle_vars = Dict{Int64,Bool}()
    to_vec, bcycle = get_vec(N,to_vec,x_vals;cycle_det=true)
    if bcycle == 1
        first = true
        while length(cycle_vars) != N
            cycle_idx, sumx = get_cycle(m, x, x_vals, N, to_vec, cycle_vars; first=first)
            first = false
            if length(cycle_idx)-1 < N
                @constraint(m, sumx <= length(cycle_idx)-2)
            end
        end
        return false
    end
    return true 
end

function prime_solved(start, m, x, N, cities, tenth, dists, score, extra, subm_path)
    to_vec = zeros(Int,N)
    x_vals = get_x_vals(x,N)
    cycle_vars = Dict{Int64,Bool}()
    to_vec, bcycle = get_vec(N,to_vec,x_vals)
    cycle_idx, _ = get_cycle(m, x, x_vals, N, to_vec, cycle_vars)
    cycle_idx = cycle_idx[1:end-1]
#     println("New tour: ", cycle_idx)
#     println("New tour length: ", length(cycle_idx))
    new_path = subm_path[start.+cycle_idx]
    new_score, new_extra, min_extra = calc_score(cities, new_path, tenth)
    println("Current New score: ", new_score)
    println("Current New score without extra: ", new_score-new_extra)
    if new_score < score
        return true
    else
        # remove the current path as it's not an improvemeny
        len_cyc = length(cycle_idx)
        sumx = 0
        for k in 1:len_cyc-1
            if cycle_idx[k] < cycle_idx[k+1]
                sumx += x[cycle_idx[k],cycle_idx[k+1]]
            else
                sumx += x[cycle_idx[k+1],cycle_idx[k]]
            end
        end
        @constraint(m, sumx <= len_cyc-2)
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
    dists[1,N] = 0
    
    # Normal MIP approach with N^2/2 edges only 1 -> 2 not 2 -> 1 
    solver = CbcSolver()
    m = Model(solver=solver)
    @variable(m, x[f=1:N,t=f+1:N], Bin)
    @objective(m, Min, sum(x[i,j]*dists[i,j] for i=1:N,j=i+1:N))
    for i=1:N
        @constraint(m, sum(x[j,i] for j=1:i-1)+sum(x[i,j] for j=i+1:N) == 2)
    end
    
    # the prime score is at least min_extra higher so we want something below 
    @constraint(m, sum(x[i,j]*dists[i,j] for i=1:N,j=i+1:N) <= score-min_extra)
    
    # N has to be connected to 1
    @constraint(m, x[1,N] == 1)

    t = time()
    status = solve(m)
    println("Status: ", status)
    
    # This shouldn't happen but also shouldn't kill the program
    if status == :Infeasible
        return subm_path, false
    end
    println("MIP Obj: ", getobjectivevalue(m))
    
    # remove cycles
    while !solved(m,x,N)
        status = solve(m)
    end
    
    # this only if we later restrict the number of cycles removed for run time maybe
    if status == :Infeasible
        return subm_path, false
    end
        
    counter = 0
    removed_cycles = 0
    max_paths = 20
    # check if the current is a better prime solution if not remove the path 
    # and try again
    while !prime_solved(start, m, x, N, cities, tenth, dists, score, extra, subm_path) && counter < max_paths
        status = solve(m)
        status == :Infeasible && break
        # again remove cycles
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
    
    improved = false
    improved_by = 0.0
    if counter < max_paths && status != :Infeasible
        println("FOUND BETTER...") 
        println("Counter: , ", counter)
        to_vec = zeros(Int,N)
        x_vals = get_x_vals(x,N)
        cycle_vars = Dict{Int64,Bool}()
        to_vec, bcycle = get_vec(N,to_vec,x_vals)
        cycle_idx, _ = get_cycle(m, x, x_vals, N, to_vec, cycle_vars)
        println("New tour: ", cycle_idx)
        new_path = subm_path[start.+cycle_idx][1:end-1]
        new_score, new_extra, _ = calc_score(cities, new_path, tenth)
        @assert score-new_score > 0
        improved_by = score-new_score
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
    return subm_path, improved, improved_by # subm_path, improved
end           

function main(from, to)
    t = time()
    cities_csv = CSV.read("cities_p.csv");
    subm_df = CSV.read("tsp_improved_mip_fast.csv");
    xy_cities   = zeros(size(cities_csv)[1],2)
    xy_cities[:,1] = cities_csv[:X]
    xy_cities[:,2] = cities_csv[:Y]
    cities = Cities(xy_cities, cities_csv[:nprimes])

    subm_path = collect(skipmissing(subm_df[:Path]));
    subm_path .+= 1
    N = 202
    total_improved_by = 0.0
    for i=from:45:to
        tenth = [(s % 10) == 0 for s in i+1:i+N-1]
        subm_path, improved, improved_by = run_mip(i, cities, tenth, subm_path, N);
        total_improved_by += improved_by
        println("Improvement per hour: ", total_improved_by*(3600/(time()-t)))
        println("Total run time: ", time()-t)
        println("Total improvement: ", total_improved_by)
        if improved
           df = DataFrame(Path=subm_path.-1)
           CSV.write("tsp_temp.csv", df);
        end
    end
end

main(0,200)
