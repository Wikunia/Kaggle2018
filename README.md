# Kaggle Santa 2018 - TRP

## MIP
to run MIP. Start Julia
```
using TRP
```

you probably first need to add it by: `] dev TRP` .

Then you can run:

```
 TRP.main_mip("combined","mip", 0,200; N=100)
```

This reads the submission file "combined.csv" and if there is an improvement it will be stored in `submissions/mip.csv`.
Starting MIP from 0 up to 200 with a path length of 100. Current stride is 45 but will be customizable.

### Parallel
You can run it by the command `julia run_par_mip.jl -i combined -o mip -p X` where `X` is the numer of processors you want to use. Then you need to install ArgParse with `] add ArgParse`.
`submissions/combined.csv` should be the submission file you want to improve on and `submissions/mip.csv` will be created or overwritten with the new current best.
This will be done whenever a new better solution was found. You can run `julia run_par_mip.jl --help` to get more information about the commands.


## K-Opt

### Usage

`run_kopt.jl -k [k (as in k-opts)] -n [number of nearest neighbors to consider] -i [input tour csv file] -o [output tour csv file]`

### How does it work?

For each edge in the tour:
1. generates all possible k-1 neighbor combinations so that it ends up with k edges.
2. from those edges generates subtours (tours going from the end of one edge to the beginning of the next)
3. looks at all possible ways of rearranging subtours (including reversing)
4. scores each possibility without considering primeness
5. for each possibility that scored higher than 0, calculate the score of the entire tour. replaces current tour with the highest scoring replacement tour for that edge (if it scores more than 0)