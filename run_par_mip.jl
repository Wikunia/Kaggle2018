using Distributed, ArgParse

function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--in", "-i"
            help = "Input csv file which is stored in submissions"
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
        "--processors", "-p"
            help = "Number of parallel processes to use"
            arg_type = Int
            required = true
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
    addprocs(args["processors"])
    println("Added "*string(args["processors"])*" processors")
    using TRP
    N = args["length"]
    stride = convert(Int, floor(N/5))
    path_fn_in = args["in"]
    path_fn_out = args["out"]
    if length(path_fn_in) <= 4 || path_fn_in[end-3:end] != ".csv"
        path_fn_in *= ".csv"
    end
    if length(path_fn_out) <= 4 || path_fn_out[end-3:end] != ".csv"
        path_fn_out *= ".csv"
    end
    # copy in to out and then just work on out
    cp("submissions/"*path_fn_in,"submissions/"*path_fn_out; force=true)

    for s = 0:stride:N-stride
        println("START with ", args["start"]+s)
        flush(stdout)
        TRP.main_mip_parallel(args["out"], args["out"]; from=args["start"]+s,
                            to=args["end"], N=N, max_mip_time=args["max_time"])
        println("FINISHED with ", args["start"]+s)
    end
end
