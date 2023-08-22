using DrWatson
@quickactivate "MixedTSL"
using ArgParse
using Random
using SPM

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--index"
            help = "index of the simulation"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

seeds = [457214195, 770082713, 680828959, 818386247, 683336764, 833841407, 753500714, 598448832, 731658602, 207305674, 897762020, 522672142, 560656253, 123879386, 501593288, 343426746, 428053117, 823585039, 747421761, 231412467, 304035054, 373743647, 504021295, 383347445, 979776442, 887799784, 30906912, 35275305, 489084581, 874055179]