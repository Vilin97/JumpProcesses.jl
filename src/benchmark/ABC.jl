using JumpProcesses, DiffEqBase
# using BenchmarkTools
using Test, Graphs

dim = 3
linear_size = 100
dims = Tuple(repeat([linear_size], dim))
num_nodes = prod(dims)
starting_site = trunc(Int, (linear_size^dim + 1) / 2)
u0 = [500, 500, 0]
end_time = 10.0
diffusivity = 1.0

domain_size = 1.0 #μ-meter
mesh_size = domain_size / linear_size
# ABC model A + B <--> C
num_species = 3
reactstoch = [[1 => 1, 2 => 1], [3 => 1]]
netstoch = [[1 => -1, 2 => -1, 3 => 1], [1 => 1, 2 => 1, 3 => -1]]
rates = [0.1 / mesh_size, 1.0]
majumps = MassActionJump(rates, reactstoch, netstoch)

# spatial system setup
hopping_rate = diffusivity * (linear_size / domain_size)^2

# Starting state setup
starting_state = zeros(Int, length(u0), num_nodes)
starting_state[:, starting_site] .= u0

tspan = (0.0, end_time)
prob = DiscreteProblem(starting_state, tspan, rates)
hopping_constants = [hopping_rate for i in starting_state]

# testing
grids = [CartesianGridRej(dims), JumpProcesses.CartesianGridRejOG(dims), Graphs.grid(dims)]
jump_problems = JumpProblem[JumpProblem(prob, NSM(), majumps,
                                        hopping_constants = hopping_constants,
                                        spatial_system = grid,
                                        save_positions = (false, false)) for grid in grids]
