# K-Opt

## Usage

`kopt.jl --k [k (as in k-opts)] --n [number of nearest neighbors to consider] --c [cities csv file (must include n_prime column)] --i [input tour csv file] --o [output tour csv file]`

## How does it work?

For each edge in the tour:
1. generates all possible k-1 neighbor combinations so that it ends up with k edges.
2. from those edges generates subtours (tours going from the end of one edge to the beginning of the next)
3. looks at all possible ways of rearranging subtours (including reversing)
4. scores each possibility without considering primeness
5. for each possibility that scored higher than 0, calculate the score of the entire tour. replaces current tour with the highest scoring replacement tour for that edge (if it scores more than 0)
