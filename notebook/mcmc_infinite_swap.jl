using StatsBase, Printf, DelimitedFiles, Combinatorics;
using MDToolbox, JLD2;
#PyPlot.plt.style.use("seaborn-colorblind");
#ENV["COLUMNS"] = 110; #display width for MDToolbox

V(x; k=1.0) = sum(- (1.0 / 2.0) .* k .* x.^2 .+ (1.0 ./ 4.0) .* k .* x.^4)
V([0.0])

proposal_fun(x) = x .+ (rand(Float64, size(x)) .* 2.0 .- 1.0);
proposal_fun([0.0])

pai(x, T) = exp.(-V(x) ./ T)

function calc_rho!(rho, x, T, perm)
    nreplica = length(x)
    nperm = length(rho)
    mu = 0.0
    for n = 1:nperm
        rho[n] = 0.0
        rho[n] = pai(x[1], T[perm[n][1]])
        for i = 2:nreplica
            rho[n] *= pai(x[i], T[perm[n][i]])
        end
        mu += rho[n]
    end
    rho .= rho./mu
    return rho
end

function select_n(rho)
    rho_cumsum = cumsum(rho)
    r = rand()
    n = sum(r .> rho_cumsum ) + 1
    return n
end

function flush_weight(io::IOStream, m, rho, perm)
    nperm = length(rho)
    nreplica = length(perm[1])
    weight = zeros(Float64, nreplica)

    for n = 1:nperm
        #id_replica = perm[n][m]
        id_replica = findall(iszero, perm[n] .- m)[1]
        weight[id_replica] += rho[n]
    end

    for i = 1:nreplica
        @printf(io, "%f ", weight[i])
    end
    @printf(io, "\n")
end

nreplica = 4
temperature_replica = [0.01, 0.10, 0.30, 0.40];
nstep = 5000;
perm=collect(permutations(1:nreplica))
nperm=factorial(nreplica)

x_replica = []
for i = 1:nreplica
    x = [0.0]
    push!(x_replica, x)
end

io_replica = []
for i = 1:nreplica
    filename = "./ISdata/replica$(i).dat"
    io = open(filename, "w")
    push!(io_replica, io)
end

io_weight=[]
for m = 1:nreplica
    filename = "./ISdata/weight$(m).dat"
    io = open(filename, "w")
    push!(io_weight, io)
end

rho=ones(Float64, nperm)
calc_rho!(rho, x_replica, temperature_replica, perm)

icount = 0
for istep = 1:nstep
    n = select_n(rho)
    for i = 1:nreplica
        x_replica[i] = propagate_mcmc(V, proposal_fun, x_replica[i], temperature_replica[perm[n][i]], nstep=1, io=io_replica[i]);
    end
    calc_rho!(rho, x_replica, temperature_replica, perm)

    for m = 1:nreplica
        flush_weight(io_weight[m], m, rho, perm)
    end
end

for i = 1:nreplica
    close(io_replica[i])
    close(io_weight[i])
end

traj_replica = []
temp_replica = []

for i = 1:nreplica
    filename = "./ISdata/replica$(i).dat"
    data = readdlm(filename);
    push!(temp_replica, data[:, 1])
    push!(traj_replica, data[:, 2])
end

#fig, ax = subplots(figsize=(8, 6))
#for i = 1:nreplica
#    ax.plot(traj_replica[i], linewidth=0.5)
#end
#xlabel("step",fontsize=20)
#ylabel("x(step)",fontsize=20)

#ax.legend(["replica 1", "replica 2", "replica 3", "replica 4"])
#tight_layout()

filename = "./ISdata/weight1.dat"
weight_replica = readdlm(filename);

global traj = traj_replica[1]
global weight = weight_replica[:, 1]

for i = 2:nreplica
    global traj = [traj; traj_replica[i]]
    global weight = [weight; weight_replica[:, i]]
end
weight .= weight ./ sum(weight)

x_grid = range(-1.5, 1.5, length=301);
pmf_theory = V.(x_grid) ./ temperature_replica[1]
pmf_theory .= pmf_theory .- minimum(pmf_theory);

#pmf_observed, _ = getpmf(traj_replica[1], weight=weight_replica[:, 1], grid_x = collect(x_grid), bandwidth=0.05);
pmf_observed, _ = getpmf(traj, weight=weight, grid_x = collect(x_grid), bandwidth=0.05);

#fig, ax = subplots(figsize=(8, 6))
#ax.plot(x_grid, pmf_theory, linewidth=3)
#xlabel("x",fontsize=20)
#ylabel("PMF (KBT)",fontsize=20)

#ax.plot(x_grid, pmf_observed, linewidth=3)

#ax.legend(["theory", "observed"])

#ax.xaxis.set_tick_params(which="major",labelsize=15)
#ax.yaxis.set_tick_params(which="major",labelsize=15)
#ax.grid(linestyle="--", linewidth=0.5)
#tight_layout()
#savefig("mcmc_infinite_swap.png", dpi=350)
delta_pmf = pmf_observed[51] - pmf_observed[251]
println(delta_pmf)

@save "IS_run.jld2" pmf_observed x_grid temperature_replica nstep
