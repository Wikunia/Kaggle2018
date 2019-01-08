
using Distributed, Distances, JuMP, Cbc

"""
    get_cycle!(cycle_vars, x_vals, to_vec, m, x, N; first=true)

Get a (sub)tour based on the current solution of the model.
Return the cycle and the number of edges (used to elimate a subtour)
"""
function get_cycle!(cycle_vars, x_vals, to_vec, m, x, N; first=true)
    cycle_idx = Int[]
    if !first
        # get start of next cycle
        for y=1:N
            if !haskey(cycle_vars,y)
                push!(cycle_idx, y)
                bcycle = get_vec!(to_vec,x_vals,cycle_vars,N;cycle_det=true,i=y)
                break
            end
        end
    else
        push!(cycle_idx, 1)
    end

    first = false
    # get the rest of the cycle done if back at the starting node
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


"""
    get_vec!(to_vec,x_vals,used,N;i=1,cycle_det=false,bcycle=false)

Get a representation of where is an edge from one node to another. Just adds starting from i.
[2,3,4,1,0,0] => Edges 1->2->3->4->1 (subtour is done therefore the rest is 0) 
At the end if cycle_det = false there are no zeros in the to_vec vector but a full representation of all subtours is obtained recursively
Attention: This changes to_vec, x_vals and used 
Return if a subtour exists
"""
function get_vec!(to_vec,x_vals,used,N;i=1,cycle_det=false,bcycle=false)
    idx = 0
    m_idx = 0
    # starting from node i
    for j = i+1:N
        # if there is an edge
        if isapprox(x_vals[i,j], 1, atol=1e-6)
            m_idx = 0
            # where the end point isn't in a tour yet
            if !haskey(used,j)
                # update the index remove the edge from x_vals and get to the next step -> recursion
                idx = j
                used[j] = true
                to_vec[i] = j
                x_vals[i,j] = 0
                bcycle = get_vec!(to_vec,x_vals,used,N;i=idx,cycle_det=cycle_det,bcycle=bcycle)
                break
            end
        end
    end
    # if there is no edge from i<->j (j > i) we check  j < i
    if idx == 0
        for j = 1:i-1
            if isapprox(x_vals[j,i], 1, atol=1e-6)
                m_idx = 0
                if !haskey(used,j)
                    idx = j
                    used[j] = true
                    x_vals[i,j] = 0
                    to_vec[i] = j
                    bcycle = get_vec!(to_vec,x_vals,used,N;i=idx,cycle_det=cycle_det,bcycle=bcycle)
                    break
                end
            end
        end
    end
    
    # if we want to detect if there is a cycle in the end
    if idx == 0 && cycle_det
        # check if there is a zero entry then we didn't get the full tour
        for y=1:N
            if to_vec[y] == 0
                return true
            end
        end
        return bcycle
    end

    # fill 
    if idx == 0 
        for y=1:N
            if to_vec[y] == 0
                bcycle = get_vec!(to_vec,x_vals,used,N;i=y,cycle_det=cycle_det,bcycle=bcycle)
                break
            end
        end
    end
    return bcycle
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
    solved(m::JuMP.Model, x::JuMP.JuMPDict, N::Int)

Checks whether there is a subtour if yes return false and the subtours will be removed
in the next solve
"""
function solved(m::JuMP.Model, x::JuMP.JuMPDict, N::Int)
    x_vals = get_x_vals(x,N)

    # find cycle
    to_vec = zeros(Int,N)
    used = Dict{Int64,Bool}()
    objective = 0
    cycle_vars = Dict{Int64,Bool}()
    bcycle = get_vec!(to_vec,x_vals,cycle_vars,N;cycle_det=true)
    # find all cycles and remove them
    if bcycle
        first = true
        while length(cycle_vars) != N
            cycle_idx, sumx = get_cycle!(cycle_vars, x_vals, to_vec, m, x, N; first=first)
            first = false
            # remove the subtour
            @constraint(m, sumx <= length(cycle_idx)-2)
        end
        return false
    end
    return true 
end

"""
    prime_improved(start, m, x, N, cities, tenth, dists, score, extra, subm_path, best_score, best_path; verbose=false)

Computes a new path and scores it with the santa score. If it's better return true else false.
@TODO might be reasonable to further check for improvements even if already an improvement was found
"""
function check_prime_improved(start, m, x, N, cities, tenth, dists, score, extra, subm_path,
            best_score, best_path; verbose=false)
    to_vec = zeros(Int,N)
    x_vals = get_x_vals(x,N)
    cycle_vars = Dict{Int64,Bool}()
    bcycle = get_vec!(to_vec,x_vals,cycle_vars,N)
    # we know that it's actually the full tour as we removed all subtours before
    cycle_idx, _ = get_cycle!(cycle_vars, x_vals, to_vec, m, x, N)
    cycle_idx = cycle_idx[1:end-1]
    new_path = subm_path[start.+cycle_idx]
    new_score, new_extra, min_extra = calc_score_with_extra(cities, new_path, tenth)
    if verbose
        println("Current New score: ", new_score)
        println("Current New score without extra: ", new_score-new_extra)
    end
   
    # remove the current path even if it might be an improvement we should look further
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
        
    if new_score < best_score
        println("FOUND better but looking for further improvement")
        println("Current improvement: ", score-new_score)
        return true, new_score, copy(new_path)
    else
        return false, best_score, best_path
    end
end

function run_mip(start, cities, tenth, subm_path, N; max_time=100, verbose=false)
    verbose && println("Start: ", start)
    println("Starting with: ", start, " N: ", N, " max time: ", max_time)
    old_subm_path = copy(subm_path)
    score, extra, min_extra = calc_score_with_extra(cities, subm_path[start+1:start+N], tenth)
    if verbose
        println("Current distance with extra: ", score)
        println("Current distance w/o extra: ", score-extra)
        println("Current extra: ", extra)
        println("Search form normal dist < ", score-min_extra)
    end
    
    mip_time = time()
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
    verbose && println("Status: ", status)
    
    # This shouldn't happen but also shouldn't kill the program
    if status == :Infeasible
        return subm_path, false
    end
    verbose && println("MIP Obj: ", getobjectivevalue(m))
    
    # remove cycles
    while !solved(m,x,N)
        status = solve(m)
    end
    
    # this only if we later restrict the number of cycles removed for run time maybe
    if status == :Infeasible
        return subm_path, false
    end

    best_score = score
    best_path = copy(subm_path[start+1:start+N])
    global_improved = false
        
    counter = 0
    # check if the current is a better prime solution if not remove the path 
    # and try again
    improved, best_score, best_path = check_prime_improved(start, m, x, N, cities, tenth, dists, score,
                                                    extra, subm_path, best_score, best_path; verbose=verbose)
    if improved 
        global_improved = true
    end
    while time()-mip_time < max_time
        status = solve(m)
        status == :Infeasible && break
        # again remove cycles
        while !solved(m,x,N)
            status = solve(m)
            status == :Infeasible && break
            if time()-mip_time >= max_time
                status = :Infeasible
                break
            end
        end
        status == :Infeasible && break
        counter += 1 
        improved, best_score, best_path = check_prime_improved(start, m, x, N, cities, tenth, dists, score,
                                                    extra, subm_path, best_score, best_path; verbose=verbose)
        if improved 
            global_improved = true
        end
    end
    
    
    improved_by = 0.0
    if global_improved
        println("FOUND BETTER...") 
        verbose && println("Counter: , ", counter)
        @assert score-best_score > 0
        improved_by = score-best_score
        println("Improved by: ", score-best_score)
        if verbose
            global_tenth = [(i % 10) == 0 for i in 1:length(subm_path)-1]
            println("oldScore: ", calc_score(cities, subm_path, global_tenth))
            subm_path[start+1:start+N] = best_path 
            println("newScore: ", calc_score(cities, subm_path, global_tenth))
        else
            subm_path[start+1:start+N] = best_path 
        end
    end
    verbose && println("Time for solving: ", time()-t)

    if verbose && status != :Infeasible
        println("Obj: ", getobjectivevalue(m))
    end
    verbose && println("Counter: , ", counter)

    return subm_path, global_improved, improved_by # subm_path, improved
end           

"""
    main_mip(path_fn_in, path_fn_out; from=0, to=197769, stride=45, N=202, max_mip_time=100)

Run x MIPs from `from` to `to` with a given stride and a length of N using a time limit, a given path of input and output path file 
"""
function main_mip(path_fn_in, path_fn_out; from=0, to=197769, stride=45, N=202, max_mip_time=100)
    t = time()
    if path_fn_in[end-3:end] != ".csv"
        path_fn_in *= ".csv"
    end
    if path_fn_out[end-3:end] != ".csv"
        path_fn_out *= ".csv"
    end
    cities, subm_path = read_cities_and_subm("cities_p.csv", "submissions/"*path_fn_in)

    total_improved_by = 0.0
    for i=from:stride:to-N-1
        tenth = [(s % 10) == 0 for s in i+1:i+N-1]
        subm_path, improved, improved_by = run_mip(i, cities, tenth, subm_path, N; max_time=max_mip_time, verbose=true);
        total_improved_by += improved_by
        println("Improvement per hour: ", total_improved_by*(3600/(time()-t)))
        println("Total run time: ", time()-t)
        println("Total improvement: ", total_improved_by)
        if improved
           df = DataFrame(Path=subm_path.-1)
           CSV.write("submissions/"*path_fn_out, df);
        end
    end
end

"""
    main_mip_parallel(path_fn_in, path_fn_out; from=0, to=197769, N=202, max_mip_time=100)

Run x MIPs from `from` to `to` and length of N using a time limit, a given path of input and output path file.
This function is run in parallel based on the number of processors available (you can add them by using addprocs(Z))
"""
function main_mip_parallel(path_fn_in, path_fn_out; from=0, to=197769, N=202, max_mip_time=100)
    next_start = from-N
    function get_next_start()
        next_start += N
        next_start > to-N-1 && return -1, []
        tenth = [(s % 10) == 0 for s in next_start+1:next_start+N-1]
        return next_start, tenth
    end
    println("From: ", from, " to ", to)
    np = nprocs()

    if length(path_fn_in) <= 4 || path_fn_in[end-3:end] != ".csv"
        path_fn_in *= ".csv"
    end
    if length(path_fn_out) <= 4 || path_fn_out[end-3:end] != ".csv"
        path_fn_out *= ".csv"
    end
    cities, subm_path = read_cities_and_subm("cities_p.csv", "submissions/"*path_fn_in)

    still_running = true
    global_start_time = time()
    total_improved_by = 0.0
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        start, tenth = get_next_start()
                        
                        if start == -1
                            still_running = false
                            break
                        end
                        println("Processor ",p, " starts with: ", start)
                        flush(stdout)

                        new_subm_path, improved, improved_by = remotecall_fetch(run_mip, p, start, cities, tenth, subm_path, N; max_time=max_mip_time, verbose=true)
                        total_improved_by += improved_by
                        println("Processor ",p, " finished with: ", start)
                        println("Improvement per hour: ", total_improved_by*(3600/(time()-global_start_time)))
                        println("Total run time: ", time()-global_start_time)
                        println("Total improvement: ", total_improved_by)
                        flush(stdout)
                        if improved
                            subm_path[start+1:start+N] = new_subm_path[start+1:start+N]
                            df = DataFrame(Path=subm_path.-1)
                            CSV.write("submissions/"*path_fn_out, df);
                        end
                        !still_running && break
                    end
                end
            end
        end
    end
end

