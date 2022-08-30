using JumpProcesses, Test, Random, StableRNGs, BenchmarkTools
const JP = JumpProcesses

rng = StableRNG(12345)
dims_list = [(d, d, d) for d in [50, 100, 150]]
num_samples = 10^2
function benchmark_rand_nbr(grid, num_samples, rng)
    nb = zero(typeof(JP.rand_nbr(rng, grid, 1)))
    for _ in 1:num_samples
        for site in 1:JP.num_sites(grid)
            nb = JP.rand_nbr(rng, grid, site)
        end
    end
    nb
end
function benchmark_nth_nbr(grid, num_samples, rng)
    nb = zero(typeof(JP.rand_nbr(rng, grid, 1)))
    for _ in 1:num_samples
        for site in 1:JP.num_sites(grid)
            n = rand(1:outdegree(grid, site))
            nb = JP.nth_nbr(grid, site, n)
        end
    end
    nb
end
function benchmark_outdegree(grid, num_samples, rng)
    sum = zero(typeof(JP.outdegree(grid, 1)))
    for _ in 1:num_samples
        for site in 1:JP.num_sites(grid)
            sum += JP.outdegree(grid, site)
        end
    end
    sum
end

function run_all_benchmarks(dims_list, num_samples, rng)
    println("Using $num_samples samples.")
    for dims in dims_list
        println("  Dimension $dims")
        for grid in [JP.CartesianGridRej(dims), JP.CartesianGridRejOG(dims)]
            println("  Grid $(typeof(grid)) with $(JP.num_sites(grid)) sites.")
            b1= @benchmarkable benchmark_rand_nbr($grid, $num_samples, $rng) samples = 5 seconds = 10
            b2 = @benchmarkable benchmark_nth_nbr($grid, $num_samples, $rng) samples = 5 seconds = 10
            b3 = @benchmarkable benchmark_outdegree($grid, $num_samples, $rng) samples = 5 seconds = 10
            println("    rand_nbr: $(run(b1))")
            println("    nth_nbr: $(run(b2))")
            println("    outdegree_nbr: $(run(b3))\n")
        end
    end
end
run_all_benchmarks(dims_list, num_samples, rng)