using Distributed, ArgParse, CSV

function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--old"
            help = "This is the file in gets compared to. This file was already run with MIP."
            required = true
        "--in", "-i"
            help = "Input csv file which is stored in submissions. MIP improved by kopt which should now be improved by MIP again."
            required = true
        "--out", "-o"
            help = "Output csv file will be stored in submissions"
            required = true
        "--start", "-s"
            help = "Starting point default: 0 (zero index based)"
            default = 0
            arg_type = Int64
        "--end", "-e"
            help = "End point set if only a partial run should be made"
            default = 197769
            arg_type = Int64
        "--length", "-N"
            help = "Length of a path solved with MIP default=202"
            arg_type = Int
            default = 202
        "--max_time", "-t"
            help = "The maximum amount of time in seconds available to find a better solution for one MIP default 100 seconds"
            arg_type = Int
            default = 100
    end

    return parse_args(s)
end

if isinteractive() == false
    args = get_args()
    using TRP
    N = args["length"]
    stride = convert(Int, floor(N/5))
    path_fn_in = args["in"]
    path_fn_out = args["out"]
    path_fn_old = args["old"]
    if length(path_fn_old) <= 4 || path_fn_old[end-3:end] != ".csv"
        path_fn_old *= ".csv"
    end
    if length(path_fn_in) <= 4 || path_fn_in[end-3:end] != ".csv"
        path_fn_in *= ".csv"
    end
    if length(path_fn_out) <= 4 || path_fn_out[end-3:end] != ".csv"
        path_fn_out *= ".csv"
    end
    # copy in to out and then just work on out
    cp("submissions/"*path_fn_in,"submissions/"*path_fn_out; force=true)

    cities, subm_path = TRP.read_cities_and_subm("cities_p.csv", "submissions/"*path_fn_in)
    subm_old_df = CSV.read("submissions/"*path_fn_old);
    subm_old_path = collect(skipmissing(subm_old_df[:Path]));
    subm_old_path .+= 1

    i = args["start"]
    while i < args["end"]-N-1
        global i
        # find current in old and compare the next N
        old_idx = findfirst(x->x.==subm_path[i], subm_old_path)
        # if not identical run one MIP
        if subm_path[i:i+N-1] != subm_old_path[old_idx:old_idx+N-1]
            println("Change in position: i: ", i)
            TRP.main_mip(args["out"], args["out"]; from=max(1,i-10),
                to=i+N+1, N=N, max_mip_time=args["max_time"])
            i += N
        end
        i += 1
    end
   
end
