# ApproximateBisection

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ApproximateBisection

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ApproximateBisection"
```
which auto-activate the project and enable local path handling from DrWatson.

2. To replicate the entire results in the accompanying manuscript, run the following command in a Unix terminal.
   ```
   ./run_all_simulations.sh
   ```
   Once the simulations are completed, the results can be collected using
   ```
   julia scripts/collect-results.jl
   ```
   Then, the tables can be obtained by running
   ```
   Rscript scripts/calculate_tables.R
   ```
Note that the simulation scripts require a SLURM workload manager to run the individual simulations. If this is not available, the scripts can be easily adapted to be run without SLURM.
