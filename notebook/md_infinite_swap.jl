using StatsBase, Printf, DelimitedFiles, Combinatorics;
using MDToolbox, JLD2;
#PyPlot.plt.style.use("seaborn-colorblind");
#ENV["COLUMNS"] = 110; #display width for MDToolbox

V(x; k=1.0) = sum(- (1.0 / 2.0) .* k .* x.^2 .+ (1.0 ./ 4.0) .* k .* x.^4)
V([0.0])

grad(x; k= 1.0) =  - k .* x .+ k .* x.^3

ISgrad(x, beta_replica, beta_sigma; k=1.0) = (-k .* x .+ k .* x.^3) ./ beta_replica .* beta_sigma

proposal_fun(x) = x .+ (rand(Float64, size(x)) .* 2.0 .- 1.0);
proposal_fun([0.0])

rho(x, T) = exp.(-V(x) ./ T)

function calc_omega!(omega, x, T, perm)
    nreplica = length(x)
    nperm=length(perm)
    qu = ones(Float64,nperm)
    qu_sum = 0.0
    for n = 1:nperm
        for i = 1:nreplica
            qu[n] *= rho(x[i], T[perm[n][i]])
        end
        qu[n] = 1/nperm * qu[n]
        qu_sum += qu[n]
    end
    for n = 1:nperm
        omega[n] = qu[n] / qu_sum
    end
end

function calc_betasigma!(beta_sigma,beta_replica,perm,omega)
    beta_sigma .= 0.0
    nreplica=length(beta_replica)
    nperm=length(perm)
    for i=1:nreplica
        for n=1:nperm
            beta_sigma[i] += beta_replica[perm[n][i]] * omega[n]
        end
    end
end

function flush_weight(io::IOStream, m, omega, perm)
    nperm = length(perm)
    nreplica = length(perm[1])
    weight = zeros(Float64, nreplica)

    for n = 1:nperm
        #id_replica = perm[n][m]
        id_replica = findall(iszero, perm[n] .- m)[1]
        weight[id_replica] += omega[n]
    end

    for i = 1:nreplica
        @printf(io, "%f ", weight[i])
    end
    @printf(io, "\n")
end

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
nexchange = 50;
perm=collect(permutations(1:nreplica))
nperm=factorial(nreplica)
omega=zeros(Float64,nperm)
beta_sigma=zeros(Float64,nreplica)
beta_replica=1 ./ temperature_replica
m2i = collect(1:nreplica)
i2m = collect(1:nreplica)


x_replica = []
for i = 1:nreplica
    x = [0.0]
    push!(x_replica, x)
end

io_replica = []
for i = 1:nreplica
    filename = "./ISMDdata/replica$(i).dat"
    io = open(filename, "w")
    push!(io_replica, io)
end



io_weight=[]
for m = 1:nreplica
    filename = "./ISMDdata/weight$(m).dat"
    io = open(filename, "w")
    push!(io_weight, io)
end



global icount = 0
global acceptance_ratio = 0.0
for iexchange = 1:nexchange
    for istep = 1:nstep
        for i = 1:nreplica
            calc_omega!(omega, x_replica, temperature_replica, perm)
            calc_betasigma!(beta_sigma, beta_replica, perm, omega)
       #    x_replica[i] = propagate_md(y -> grad(y,k=2.0), x_replica[i], temperature_replica[perm[n][i]], nstep=1, io=io_replica[i]);
            x_replica[i] = propagate_md(y -> ISgrad(y,beta_replica[i2m[i]],beta_sigma[i2m[i]],k=1.0), x_replica[i], temperature_replica[i2m[i]], nstep=1, io=io_replica[i]);
        end

        for  m=1:nreplica
          flush_weight(io_weight[m], m, omega, perm)
        end
    end
    global acceptance_ratio += exchange_temperature!(m2i, i2m, icount, x_replica, y -> V(y, k=1.0), temperature_replica)
    global icount += 1
end
acceptance_ratio = acceptance_ratio / nexchange

for i = 1:nreplica
    close(io_replica[i])
    close(io_weight[i])
end

traj_replica = []
temp_replica = []

for i = 1:nreplica
    filename = "./ISMDdata/replica$(i).dat"
    data = readdlm(filename);
    push!(temp_replica, data[:, 1])
    push!(traj_replica, data[:, 2])
end

#fig, ax = subplots(figsize=(8, 6))
#for i = 1:nreplica#
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


#fig, ax = subplots(figsize=(8, 6))
#for i = 1:nreplica#
#    ax.plot(traj_sorted[i], linewidth=0.5)
#end
#xlabel("step",fontsize=20)
#ylabel("x(step)",fontsize=20)

#ax.legend(["temperature 1", "temperature 2", "temperature 3", "temperature 4"])
#tight_layout()

#fig, ax = subplots(figsize=(8, 6))

#    ax.plot(traj_sorted[1], linewidth=0.5)

#xlabel("step",fontsize=20)
#ylabel("x(step)",fontsize=20)

#ax.legend(["temperature 3"])
#tight_layout()

filename = "./ISMDdata/weight1.dat"
weight_replica = readdlm(filename);
weight_replica


global traj = traj_sorted[1]
global weight = weight_replica[:, 1]

for i = 2:nreplica
    global traj = [traj; traj_sorted[i]]
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
#savefig("exmcmc_infinite_swap.png", dpi=350)


delta_pmf = pmf_observed[51] - pmf_observed[251]
println(delta_pmf)

@save "ISMD_run.jld2" pmf_observed x_grid traj_sorted temp_sorted temperature_replica nstep nexchange
