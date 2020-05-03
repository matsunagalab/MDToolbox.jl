"""
write to a trajectory file
"""
function flush_snapshot(io::IOStream, x::Vector{Float64}, temperature::Float64)
    @printf(io, "%f", temperature)
    for i in length(x)
        @printf(io, " %f", x[i])
    end
    @printf(io, "\n")
end


"""
propagate_mcmc

Markov Chain Monte Carlo test code
"""
function propagate_mcmc(pot_fun::Function, proposal_fun::Function, x_current::Vector{Float64}, temperature::Float64; nstep=1, io=nothing, outstep=1)::Vector{Float64}
    x = copy(x_current)
    x_proposal = similar(x)
    beta = 1 / temperature
    for istep = 1:nstep
        x_proposal .= proposal_fun(x)
        delta_pot = pot_fun(x_proposal) - pot_fun(x)
        if exp( -beta * delta_pot) > rand()
            x .= x_proposal
        end
        if isa(io, IOStream) && (mod(istep, outstep) == 0)
            flush_snapshot(io, x, temperature)
        end
    end
    return x
end

"""
propagate_md

Molecular dynamics test code
"""
function propagate_md(pot_fun::Function, proposal_fun::Function, x_current::Vector{Float64}; temperature::Float64, nstep=1, delta_t=1.0::Float64)
    return nothing
end
