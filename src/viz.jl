"""
viewstruct

based on Bio3Dview.jl
"""

#using .Bio3DView

function viewstruc(ta::TrjArray; kwargs...)
    io = IOBuffer()
    writepdb(io, ta)
    return Bio3DView.view("data-type='pdb'", String(take!(io)); style=Bio3DView.defaultstyle("pdb"), kwargs...)
end
