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
    end

    return parse_args(s)
end

if isinteractive() == false
    args = get_args()
    addprocs(args["processors"])
    println("Added "*string(args["processors"])*" processors")
    using TRP
    TRP.main_mip_parallel(args["in"], args["out"]; from=args["start"], to=args["end"], N=args["length"])
end
