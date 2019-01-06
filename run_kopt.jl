using ArgParse, TRP

function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Input tour csv file"
            required = true
        "--output", "-o"
            help = "Output tour csv file"
            required = true
        "--k", "-k"
            help = "k as in k-opts"
            arg_type = Int
            required = true
        "--neighbors", "-n"
            help = "number of neighbors to consider each time"
            default = 100
            arg_type = Int64
        "--start", "-s"
            help = "start for i"
            default = 1
            arg_type = Int64
        "--end", "-e"
            help = "end for i"
            default = 197769
            arg_type = Int64
    end

    return parse_args(s)
end

if isinteractive() == false
    args = get_args()

    path_fn_in = args["input"]
    path_fn_out = args["output"]
    if length(path_fn_in) <= 4 || path_fn_in[end-3:end] != ".csv"
        path_fn_in *= ".csv"
    end
    if length(path_fn_out) <= 4 || path_fn_out[end-3:end] != ".csv"
        path_fn_out *= ".csv"
    end

    TRP.main_kopt(path_fn_in, path_fn_out, args["k"]; list_of_i = args["start"]:args["end"], neighbors = args["neighbors"])
end
