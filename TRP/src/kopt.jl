
using Distances, NearestNeighbors, Combinatorics, IterTools

function kopt(cities, tour, num_opts, num_neighbors, path_fn_out)
	function init()
		print("Finding all nearest neighbors...")

		# KDTree for finding nearest neighbor
		# Need to swap dimensions of cities_xy for KDTree
		cities_xy_kd = zeros(2,size(cities_xy, 1))
		for i in 1:size(cities_xy, 1)
			cities_xy_kd[1:2, i] = cities_xy[tour[i],1:2]
		end
		kdtree = KDTree(cities_xy_kd)

		# Filter and store KDTree results for easy access
		nearest_neighbors = zeros(Int, length(tour)-1, num_neighbors)
		for i in 1:length(tour)-1
			neighbors, _ = knn(kdtree, cities_xy[tour[i], 1:2], num_neighbors+7, true)
			filter!(x->(x < i-3 || x > i+3), neighbors)
			nearest_neighbors[i,1:num_neighbors] = neighbors[1:num_neighbors]
		end
		println(" Done.")

		println("Precompute path lengths...")
		# precomputation of path length 10 forward and ten reverse
		tour_santa_score = Vector{Vector{Float64}}(undef, 10)
		rev_tour_santa_score = Vector{Vector{Float64}}(undef, 10)
		# forward
		gtenth = [(i % 10) == 0 for i in 1:length(tour)-1]
		forward_dists, forward_extras = calc_score_dists_and_extras(cities, tour, gtenth)
		tour_santa_score[1] = rev_cumsum(forward_dists .+ forward_extras)
		println("tour_santa_score[1][1]: ", tour_santa_score[1][1])
		for k = 1:9
			temp_tenth = [(i % 10) == k for i in 1:length(tour)-1]
			tour_santa_score[k+1] = rev_cumsum(forward_dists .+ calc_score_extras(cities, tour, temp_tenth, forward_dists))
		end


		# backward
		rev_tour = reverse(tour)
		backward_dists, backward_extras = calc_score_dists_and_extras(cities, rev_tour, gtenth)
		rev_tour_santa_score[1] = rev_cumsum(backward_dists .+ backward_extras)
		for k = 1:9
			temp_tenth = [(i % 10) == k for i in 1:length(tour)-1]
			rev_tour_santa_score[k+1] = rev_cumsum(backward_dists .+ calc_score_extras(cities, rev_tour, temp_tenth, backward_dists))
		end
		println(" Done.")

		euc_dict = Dict{Tuple{Int,Int},Float64}()
		return nearest_neighbors, tour_santa_score, rev_tour_santa_score, euc_dict
	end
	cities_xy = cities.xy
    cities_nprime = cities.nprime

	nearest_neighbors, tour_santa_score, rev_tour_santa_score, euc_dict = init()
	
	# Set parameters
	current_score = score_tour(cities_xy, cities_nprime, tour)
	num_subtours = num_opts-1
	subtours = zeros(Int, 2, num_subtours, 2)
	subtours_pos = collect(permutations(collect(1:num_subtours)))
	d_subtours_pos = permutations_repetition([1,2], num_subtours)
	nsubtours_pos = size(subtours_pos, 1)
	nd_subtours_pos = size(d_subtours_pos, 1)
	nneighbor_perms = nsubtours_pos*nd_subtours_pos
	total_promising_swaps = 0
	total_pos_heuristic = 0
	total_gain = 0.0
	total_time = time()

	neighbor_permutations = Vector{Vector{Int}}(undef, nneighbor_perms)
	permutation = Vector{Int}(undef, 2*num_subtours)

	# Iterate through each city in tour
	for i in 40000:length(tour)-1
		t = time()
		best_swap_score = 0.0

		# Find nearest neighbors
		neighbor_pos = nearest_neighbors[i,1:num_neighbors]

		# Find all promising k-opt swaps
		swaps = Vector{Vector{Int}}(undef,1)

		# Generate all possible k-1 neighbor combinations
		neighbor_combinations = collect(subsets(neighbor_pos, num_subtours))
		# For each combination
		for j in 1:size(neighbor_combinations, 1)

			# Sort deletion points (first city in each edge)
			del_points = sort([i;neighbor_combinations[j]])

			# Filter out deletion points that are direclty connected to each other
			if del_valid(del_points)
				
				# Set subtours
				for k in 1:num_subtours
					subtours[1,k,1], subtours[1,k,2] = del_points[k]+1,del_points[k+1] 
					subtours[2,k,1:2] = reverse(subtours[1,k,1:2])
				end

				# Generate all possible k-opt swaps by rearranging subtours
				kl = 1
				@inbounds for k in 1:nsubtours_pos
					for l in 1:nd_subtours_pos
						for m in 1:num_subtours
							permutation[2*m-1:2*m] = subtours[d_subtours_pos[l,m], subtours_pos[k][m], 1:2]
						end
						neighbor_permutations[kl] = permutation[:]
						kl += 1
					end
				end		

                # Calculate score differential for k-opt swaps, only saving promising ones
				swap = [del_points[1];neighbor_permutations[1];del_points[num_opts]+1]
				# tenth = [(s % 10) == 0 for s in swap[1]:swap[length(swap)]-1]
				# base_score, base_extra, base_min_extra = c5optalc_score_with_extra(cities, get_swap_subtour(tour, swap, num_subtours), tenth)
				base_score = calc_score(cities, tour, tour_santa_score, rev_tour_santa_score, swap, euc_dict)
				# println("base score: ", base_score)
				# println("base score1: ", base_score1)
				# println("abs(diff): ", abs(base_score-base_score1))
				# @assert isapprox(base_score,base_score1)
				# base_euc_score = score_swap(cities_xy, tour, swap)
				for k in 2:size(neighbor_permutations,1)
					@inbounds swap[2:end-1] = neighbor_permutations[k]
					# println("swap: ", swap)
					# println("swap: ", swap)
					total_pos_heuristic += 1
					# swap_subtour = get_swap_subtour(tour, swap, num_subtours)
					# old_score = calc_score(cities, swap_subtour, tenth)
					new_score = calc_score(cities, tour, tour_santa_score, rev_tour_santa_score, swap, euc_dict)
					# err
					# println("old score: ", old_score)
					# println("new score: ", new_score)
					# println("abs(diff): ", abs(old_score-new_score))
					# @assert isapprox(old_score,new_score)
					score_dif = base_score-new_score
					# println("score_dif: ", score_dif)
					# println("======================================================")
					if score_dif > best_swap_score
						best_swap_score = score_dif
						swaps[1] = copy(swap)
					end
				end
			end
		end

		# Fully score each promising k-opt swap
		best_gain = 0.0
		best_tour = []
		if best_swap_score > 0
			total_promising_swaps += 1
			for j in 1:size(swaps, 1)
				swap = swaps[j]

				# Generate full tour from swap
				new_tour = tour[1:swap[1]]
				for k in 1:num_subtours
					if swap[2*k] < swap[2*k+1]
						append!(new_tour, tour[swap[2*k]:swap[2*k+1]])
					else
						append!(new_tour, reverse(tour[swap[2*k+1]:swap[2*k]]))
					end
				end
				append!(new_tour, tour[swap[length(swap)]:length(tour)])

        	    gain = current_score - score_tour(cities_xy, cities_nprime, new_tour)
				if gain > best_gain
					best_gain = gain
					best_tour = new_tour
				end
			end

			# Set current tour to best k-opt swap
			if best_gain > 0
				total_gain += best_gain
				current_score -= best_gain
				tour = best_tour[:]
				df = DataFrame(Path=best_tour .-= 1)
				CSV.write("submissions/"*path_fn_out, df)
				nearest_neighbors, tour_santa_score, rev_tour_santa_score, euc_dict = init()
			end
		end

		t = time()-t
		println("i: $i, total pos heuristic: $total_pos_heuristic, total promising swaps: $total_promising_swaps, total gain: $total_gain, time for one i: $t")
		if i % 10 == 0
			tt = time() - total_time
			println("Total time for $total_pos_heuristic checks: ", tt, " that are  ", total_pos_heuristic/tt, " checks per second.")
		end
		flush(stdout)
	end
end

function get_swap_subtour(tour, swap, num_subtours)
	len_subtour = swap[length(swap)]-swap[1]+1
	new_tour = zeros(Int, len_subtour)
	new_tour[1] = tour[swap[1]]
	c = 2
	for k in 1:num_subtours
		swap2k, swap2k1 = swap[2*k], swap[2*k+1]
        if swap2k < swap2k1
			new_tour[c:c+swap2k1-swap2k] = tour[swap2k:swap2k1]
			c += swap2k1-swap2k+1
        else
			new_tour[c:c+swap2k-swap2k1] = reverse(tour[swap2k1:swap2k])
			c += swap2k-swap2k1+1
        end
	end
	new_tour[end] = tour[swap[length(swap)]]
    return new_tour
end

function length_of_swap(swap)
	return swap[length(swap)]-swap[1]+1
end

function score_swap(cities_xy, tour, swap)
	score = 0.0
	for i in 1:Int(size(swap, 1)/2)
		p1 = [cities_xy[tour[swap[2*i-1]],1], cities_xy[tour[swap[2*i-1]],2]]
		p2 = [cities_xy[tour[swap[2*i]],1], cities_xy[tour[swap[2*i]],2]]
		score += euclidean(p1,p2)
	end
	return score
end

function score_tour(cities_xy, cities_nprime, scored_tour)
	score = 0.0
	for i in 1:length(scored_tour)-1
		p1 = [cities_xy[scored_tour[i],1], cities_xy[scored_tour[i],2]]
		p2 = [cities_xy[scored_tour[i+1],1], cities_xy[scored_tour[i+1],2]]
		d = euclidean(p1,p2)
		if i%10==0&&cities_nprime[scored_tour[i]]==1.0
			d*=1.1
		end
		score+=d
	end
	return score
end

function del_valid(del_points)
	for i in 1:length(del_points)-1
		if del_points[i]+1 == del_points[i+1]
			return false
		end
	end
	return true
end

function permutations_repetition(v, n)
  l = length(v)
  m = l^n
  a = zeros(Int, m, n)
  for i in 1:m
    k = i-1
    j = n
    while (j > 0)
      (d,r) = divrem(k,l)
      a[i,j] = v[r+1]
      k = d
      j -= 1
    end
  end
  return a
end

function main_kopt(path_fn_in, path_fn_out, k; neighbors = 100)
    # Load city coordinates and tour
    cities, tour = read_cities_and_subm("cities_p.csv", "submissions/"*path_fn_in)

    # k-opt only considering neighbors nearest neighbors
    kopt(cities, tour, k, neighbors, path_fn_out)
end

