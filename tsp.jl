
using CSV, DataFrames, Distances, DelimitedFiles, Plots, JuMP, Cbc

cities = CSV.read("cities_p.csv");

primes = findall(cities[:primes] .== true);

subm_df = CSV.read("new_submission_91000.csv");

subm_path = collect(skipmissing(subm_df[:Path]));

function get_score(cities, subm_path)
    all_ids = cities[:CityId]
    all_x = cities[:X]
    all_y = cities[:Y]

    incs = 0
    score = 0.0
    pimp = 0.0
    for i in 1:length(subm_path)-1
        c_idx = subm_path[i]+1
        n_idx = subm_path[i+1]+1
        p1 = [all_x[c_idx],all_y[c_idx]]
        p2 = [all_x[n_idx],all_y[n_idx]]
        stepSize = euclidean(p1,p2)
        if i % 10 == 0 && !cities[:primes][subm_path[i]+1]
            pimp += 0.1*stepSize
            stepSize *= 1.1
            incs += 1
        end
        # print(stepSize)
        score += stepSize
    end
    return score
end

function solved(m,x,N)
    x_val = getvalue(x)
    
    solved = true
    for s=1:N-1
        for i=1:N,j=1:N
            if i != j
                if x_val[i,j,s] >  sum(x_val[j,1:N,s+1])
                    @constraint(m, x[i,j,s] <= sum(x[j,1:N,s+1]));
                    solved = false
                end
            end
        end
    end
    return solved
end 

function run_test(start, subm_path)
    if start % 2000 == 0
        println("Start: ", start)
    end
    N = 13
    c_pos = [zeros(2) for _ in 1:N]

    for i = 1:N
        c_pos[i] = [cities[:X][subm_path[start+i]+1],cities[:Y][subm_path[start+i]+1]]
    end
    dists = zeros(N,N)
    for i=1:N
        for j=i+1:N
            dists[i,j] = euclidean(c_pos[i],c_pos[j])
            dists[j,i] = dists[i,j]
        end
    end
    dists[N,1] = 0

    extras = ones(N,N)
    for i=1:N
        if !cities[:primes][subm_path[i]+1]
            for k = 1:convert(Int,floor(N/10))
               extras[i,10*k-(start % 10)] = 1.1 
            end
        end
    end
    
    m = Model(solver=CbcSolver(seconds=10))
    @variable(m, x[1:N,1:N,1:N], Bin)
    @objective(m, Min, sum(x[i,j,s]*dists[i,j]*extras[i,s] for i=1:N,j=1:N,s=1:N));
    @constraint(m, notself[i=1:N], sum(x[i,i,1:N]) == 0);
    @constraint(m, eachsteponce[s=1:N], sum(x[1:N,1:N,s]) == 1);
    @constraint(m, oneout[i=1:N], sum(x[i,1:N,1:N]) == 1);
    @constraint(m, onein[j=1:N], sum(x[1:N,j,1:N]) == 1);
    
    for s=1:N-1
        for i=1:N,j=1:N
            if i != j
    #             println("If ",i," -> ", j , " in ", s)
    #             println("Then next step from: ", j)
                @constraint(m, x[i,j,s] <= sum(x[j,1:N,s+1]));
            end
        end
    end
    @constraint(m, sum(x[1,1:N,1]) == 1);
    @constraint(m, sum(x[1:N,N,N-1]) == 1);

    #=
    for f=1:N, t=1:N
        @constraint(m, x[f,t]+x[t,f] <= 1)
    end
    =#

    status = solve(m)

    if status != :Optimal
        println("Time limit")
        return subm_path, false
    end
    x_val = getvalue(x)
    
    
    # find cycle
    cycle_idx = []
    push!(cycle_idx, 1)
    while true
        v, idx = findmax(x_val[cycle_idx[end],1:N,1:N])
        if idx[1] == cycle_idx[1]
            break
        else
            push!(cycle_idx,idx[1])
        end
    end

    improved = false
#     println("cycle_idx: ", cycle_idx)
    if cycle_idx != collect(1:N)
        base_score = get_score(cities, subm_path)
        new_subm_path = copy(subm_path)
        new_subm_path[start+1:start+N] = new_subm_path[start.+cycle_idx]
        new_score = get_score(cities, new_subm_path)
        println("Old: ", base_score)
        println("New: ", new_score)
        println("Improved by: ",base_score-new_score)
        if base_score-new_score > 0
            subm_path = new_subm_path 
            improved = true
        end
    end
        
#     println("Obj: ", getobjectivevalue(m))
    # println("Cus to Fac: ",getvalue(cf))

    return subm_path, improved
end            

function main()
    subm_df = CSV.read("tsp_improved.csv");

    subm_path = collect(skipmissing(subm_df[:Path]));
    for i=5:10:size(cities)[1]
        subm_path, improved = run_test(i,subm_path);
        if improved
           df = DataFrame(Path=subm_path)
           CSV.write("tsp_improved.csv", df);
        end
    end
end
main()
