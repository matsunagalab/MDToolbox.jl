using StatsBase, Printf, DelimitedFiles;
using MDToolbox, JLD2;
#PyPlot.plt.style.use("seaborn-colorblind");
#ENV["COLUMNS"] = 110; #display width for MDToolbox

V(x; k=1.0) = sum(- (1.0 / 2.0) .* k .* x.^2 .+ (1.0 ./ 4.0) .* k .* x.^4)
V([0.0], k=1.0)

proposal_fun(x) = x .+ (rand(Float64, size(x)) .* 2.0 .- 1.0);
proposal_fun([0.0])

function exchange_temperature!(m2i, i2m, icount, x_replica, pot_fun::Function, temperature_replica)
    nreplica = length(x_replica)

    m_array = 1:nreplica
    t_array = temperature_replica[m_array]
    b_array = 1.0 ./ t_array
    i_array = m2i[m_array]
    v_array = map(pot_fun, x_replica[i_array])

    if mod(icount, 2) == 0
        m_lower = 1:2:(nreplica-1)
        m_higher = 2:2:nreplica
    else
        m_lower = 2:2:(nreplica-1)
        m_higher = 3:2:nreplica
    end

    iaccepted = 0
    for ipair = 1:length(m_higher)
        m1 = m_lower[ipair]
        m2 = m_higher[ipair]
        delta = (b_array[m2] - b_array[m1]) * (v_array[m1] - v_array[m2])
        if exp(-delta) > rand()
            m2i[m_array[m1]], m2i[m_array[m2]] = m2i[m_array[m2]], m2i[m_array[m1]]
            i2m[i_array[m1]], i2m[i_array[m2]] = i2m[i_array[m2]], i2m[i_array[m1]]
            iaccepted += 1
        end
    end

    return iaccepted / length(m_higher)
end

nreplica = 4
temperature_replica = [0.01, 0.10, 0.30, 0.40];
nstep = 100;
nexchange = 50
m2i = collect(1:nreplica)
i2m = collect(1:nreplica)

x_replica = []
for i = 1:nreplica
    x = [0.0]
    push!(x_replica, x)
end

io_replica = []
for i = 1:nreplica
    filename = "./MCMCdata/replica$(i).dat"
    io = open(filename, "w")
    push!(io_replica, io)
end

icount = 0
for iexchange = 1:nexchange
    global icount;
    for i = 1:nreplica
        x_replica[i] = propagate_mcmc(y -> V(y, k=1.0), proposal_fun, x_replica[i], temperature_replica[i2m[i]], nstep=nstep, io=io_replica[i]);
    end
    # do exchange
    exchange_temperature!(m2i, i2m, icount, x_replica, y -> V(y, k=1.0), temperature_replica)
    icount += 1
end

for i = 1:nreplica
    close(io_replica[i])
end

traj_replica = []
temp_replica = []

for i = 1:nreplica
    filename = "./MCMCdata/replica$(i).dat"
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

# sort trajectories according to temperature
function sort_traj(traj_replica, temp_replica)
    traj_sorted = deepcopy(traj_replica)
    temp_sorted = deepcopy(temp_replica)
    nframe = size(traj_replica[1], 1)
    for iframe = 1:nframe
        temp_snapshot = map(x -> x[iframe], temp_replica)
        p = sortperm(temp_snapshot)
        for m = 1:nreplica
            traj_sorted[m][iframe, :] .= traj_replica[p[m]][iframe, :]
            temp_sorted[m][iframe, :] .= temp_replica[p[m]][iframe, :]
        end
    end
    return traj_sorted, temp_sorted
end

traj_sorted, temp_sorted = sort_traj(traj_replica, temp_replica)

temp_sorted[1]

#fig, ax = subplots(figsize=(8, 6))
#for i = 1:nreplica
#    ax.plot(traj_sorted[i], linewidth=0.5)
#end
#xlabel("step",fontsize=20)
#ylabel("x(step)",fontsize=20)

#ax.legend(["temperature 1", "temperature 2", "temperature 3", "temperature 4"])
#tight_layout()

x_grid = range(-1.5, 1.5, length=301);
pmf_theory = V.(x_grid, k=1) ./ temperature_replica[1]
pmf_theory .= pmf_theory .- minimum(pmf_theory);

pmf_observed, _ = getpmf(traj_sorted[1], grid_x = collect(x_grid), bandwidth=0.05);

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
#savefig("mcmc_replica_exchange.png", dpi=350)

delta_pmf = pmf_observed[51] - pmf_observed[251]
println(delta_pmf)

@save "run.jld2" pmf_observed x_grid traj_sorted temp_sorted temperature_replica nstep nexchange
