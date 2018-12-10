
using Distances, NearestNeighbors, Combinatorics, IterTools

function kopt(cities, tour, num_opts, num_neighbors, path_fn_out)
	cities_xy = cities.xy
    cities_nprime = cities.nprime

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
		neighbors, _ = knn(kdtree, cities_xy[tour[i], 1:2], num_neighbors+3, true)
		filter!(x->x≠i&&x≠i+1&&x≠i-1, neighbors)
		nearest_neighbors[i,1:num_neighbors] = neighbors[1:num_neighbors]
	end
	println(" Done.")

	# Set parameters
	current_score = score_tour(cities_xy, cities_nprime, tour)
	num_subtours = num_opts-1
	subtours = zeros(Int, 2, num_subtours, 2)
	subtours_pos = collect(permutations(collect(1:num_subtours)))
	d_subtours_pos = permutations_repetition([1,2], num_subtours)
	total_promising_swaps = 0
	total_gain = 0.0

	# Iterate through each city in tour
	for i in 105900:length(tour)-1
		t = time()

		# Find nearest neighbors
		neighbor_pos = nearest_neighbors[i,1:10]

		# Find all promising k-opt swaps
		swaps = Vector{Vector{Int}}()

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
				neighbor_permutations = Vector{Vector{Int}}()
				for k in 1:size(subtours_pos, 1)
					for l in 1:size(d_subtours_pos, 1)
						permutation = []
						for m in 1:num_subtours
							append!(permutation, subtours[d_subtours_pos[l,m], subtours_pos[k][m], 1:2])
						end
						push!(neighbor_permutations, permutation)
					end
				end

                # Calculate score differential for k-opt swaps, only saving promising ones
				swap = [del_points[1];neighbor_permutations[1];del_points[num_opts]+1]
				tenth = [(s % 10) == 0 for s in swap[1]:swap[length(swap)]-1]
				base_score = calc_score(cities, get_swap_subtour(tour, swap, num_subtours), tenth)
				for k in 2:size(neighbor_permutations,1)
					swap = [del_points[1];neighbor_permutations[k];del_points[num_opts]+1]
					score_dif = base_score - calc_score(cities, get_swap_subtour(tour, swap, num_subtours), tenth)
                    if score_dif > 0
						push!(swaps, swap)
					end
				end
			end
		end

		# Fully score each promising k-opt swap
		total_promising_swaps += size(swaps, 1)
		best_gain = 0.0
		best_tour = []
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
			for i in 1:length(best_tour)
				best_tour[i]-=1
			end
			df = DataFrame(Path=best_tour)
        	CSV.write("submissions/"*path_fn_out, df)
		end

		t = time()-t
		println("i: $i, total promising swaps: $total_promising_swaps, total gain: $total_gain, time for one i: $t")
	end
end

function get_swap_subtour(tour, swap, num_subtours)
    new_tour = [tour[swap[1]]]
    for k in 1:num_subtours
        if swap[2*k] < swap[2*k+1]
            append!(new_tour, tour[swap[2*k]:swap[2*k+1]])
        else
            append!(new_tour, reverse(tour[swap[2*k+1]:swap[2*k]]))
        end
	end
	append!(new_tour, tour[swap[length(swap)]])
    return new_tour
end

function length_of_swap(num_subtours, swap)
	tour_length = 0
	for k in 1:num_subtours
		if swap[2*k] < swap[2*k+1]
			tour_length += swap[2*k+1]-swap[2*k]+1
		else
			tour_length += swap[2*k]-swap[2*k+1]+1
		end
	end	
	return tour_length
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

