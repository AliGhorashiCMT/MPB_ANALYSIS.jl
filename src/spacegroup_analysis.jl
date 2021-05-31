"""
$(TYPEDSIGNATURES)

Return the fraction of simulations from a specific spacegroup that have fragile topologies. 
"""
function findfragiles(sgnum::Integer, dir::AbstractString= "./symeigs/output/output")
    validfiles = String[]
    numfragiles = 0 
    for file in readdir(dir)
        contains(file, "dim2-sg$(sgnum)") || continue
        contains(file, "symeigs") || continue
        push!(validfiles, file)
    end
    for file in validfiles
        println(file[1:end-12])
        #
        try 
            (sum([x[2] for x in label_topologies(file[1:end-12], true, dir, printisbandstruct=false, verbose=false)] .== FRAGILE) >= 1) && (numfragiles +=1)
        catch
        end
    end
    return numfragiles
end

"""
$(TYPEDSIGNATURES)

Plot a random lattice from the bash file that defines it. Overloads PyPlot's plot method. 
"""
function plot_lattice(filepath::String; kwargs...)
    Rs, flat, isoval, _ = lattice_from_mpbparams(filepath)
    pygui(true)
    PyPlot.plot(flat, Rs; isoval=isoval, kwargs...) # via Crystalline.jl overload of PyPlot.jl's `plot`
end
