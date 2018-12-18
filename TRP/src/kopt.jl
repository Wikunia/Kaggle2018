
using Distances, NearestNeighbors, Combinatorics, IterTools

function kopt(cities, tour, num_opts, num_neighbors, path_fn_out, list_of_i)
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
		push!(tour_santa_score[1], 0.0)
		println("tour_santa_score[1][1]: ", tour_santa_score[1][1])
		for k = 1:9
			temp_tenth = [(i % 10) == k for i in 1:length(tour)-1]
			tour_santa_score[k+1] = rev_cumsum(forward_dists .+ calc_score_extras(cities, tour, temp_tenth, forward_dists))
			push!(tour_santa_score[k+1], 0.0)
		end


		# backward
		rev_tour = reverse(tour)
		backward_dists, backward_extras = calc_score_dists_and_extras(cities, rev_tour, gtenth)
		rev_tour_santa_score[1] = rev_cumsum(backward_dists .+ backward_extras)
		push!(rev_tour_santa_score[1], 0.0)
		for k = 1:9
			temp_tenth = [(i % 10) == k for i in 1:length(tour)-1]
			rev_tour_santa_score[k+1] = rev_cumsum(backward_dists .+ calc_score_extras(cities, rev_tour, temp_tenth, backward_dists))
			push!(rev_tour_santa_score[k+1], 0.0)
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
	total_scorings = 0
	total_gain = 0.0
	total_time = time()

	swap = Vector{Int}(undef, 2+2*num_subtours)
	hlswap = fld(2+2*num_subtours, 2)
	cfroms = Vector{Int}(undef, hlswap)
	ltour = length(tour)
    ltour1 = 1+ltour
    lswap = length(swap)
    hlswap = fld(lswap, 2)

	# Iterate through each city in tour
	@inbounds for i in list_of_i
		t = time()
		best_swap_score = 0.0

		# Find nearest neighbors
		neighbor_pos = nearest_neighbors[i,1:num_neighbors]

		# Find all promising k-opt swaps
		best_swap = Vector{Int}(undef, 2+2*num_subtours)

		# Generate all possible k-1 neighbor combinations
		neighbor_combinations = collect(subsets(neighbor_pos, num_subtours))
		# For each combination
		tc = 0.0
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
				
				swap[1] = del_points[1]
				swap[end] = del_points[num_opts]+1
				for m in 1:num_subtours
					swap[2*m] = subtours[d_subtours_pos[1,], subtours_pos[1][m], 1]
					swap[2*m+1] = subtours[d_subtours_pos[1,], subtours_pos[1][m], 2]
				end
				# current base score
				base_score = calc_score(tour_santa_score, 1, swap[1], swap[end])
				# Generate all possible k-opt swaps by rearranging subtours
				# Calculate score differential for k-opt swaps, only saving the best swap
				total_scorings += nsubtours_pos*nd_subtours_pos-1
				for k in 1:nsubtours_pos
					# don't calculate the base score again
					if k == 1
						lrange = 2:nd_subtours_pos
					else
						lrange = 1:nd_subtours_pos
					end
					for l in lrange
						# calculate the current swap don't changing del points (first and last in swap)
						# this is faster than setting swap[2*m:2*m+1] =  subtours[..,..,1:2]
						for m in 1:num_subtours
							swap[2*m] = subtours[d_subtours_pos[l,m], subtours_pos[k][m], 1]
							swap[2*m+1] = subtours[d_subtours_pos[l,m], subtours_pos[k][m], 2]
						end
						score_dif = improved_by!(cfroms, cities_xy, cities_nprime, tour, tour_santa_score, rev_tour_santa_score, swap, euc_dict, base_score,
												ltour,ltour1,lswap,hlswap)
						if score_dif > best_swap_score
							best_swap_score = score_dif
							best_swap = copy(swap)
						end
					end
				end
			end
		end

		# Fully score each promising k-opt swap
		best_gain = 0.0
		if best_swap_score > 0
			total_promising_swaps += 1

			# Generate full tour from swap
			best_tour = tour[:]
			pos = swap[1]+1
			for k in 1:num_subtours
				if best_swap[2*k] < best_swap[2*k+1]
					best_tour[pos:pos+best_swap[2*k+1]-best_swap[2*k]] = tour[best_swap[2*k]:best_swap[2*k+1]]
					pos += best_swap[2*k+1]-best_swap[2*k]+1
				else
					best_tour[pos:pos+best_swap[2*k]-best_swap[2*k+1]] = reverse(tour[best_swap[2*k+1]:best_swap[2*k]])
					pos += best_swap[2*k]-best_swap[2*k+1]+1
				end
			end

			# compute with better accuracy
			best_gain = current_score - score_tour(cities_xy, cities_nprime, best_tour)
			

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
		println("i: $i, total pos heuristic: $total_scorings, total promising swaps: $total_promising_swaps, total gain: $total_gain, time for one i: $t")
		if i % 10 == 0
			tt = time() - total_time
			println("Total time for $total_scorings checks: ", tt, " that are  ", total_scorings/tt, " checks per second.")
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

function main_kopt(path_fn_in, path_fn_out, k; list_of_i=undef, neighbors = 100)
	# Load city coordinates and tour
    cities, tour = read_cities_and_subm("cities_p.csv", "submissions/"*path_fn_in)

	if list_of_i == undef
		list_of_i = 1:197769
	end

    # k-opt only considering neighbors nearest neighbors
    kopt(cities, tour, k, neighbors, path_fn_out, list_of_i)
end

