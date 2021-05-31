"""
$(TYPEDSIGNATURES)

Return the fraction of simulations from a specific spacegroup that have fragile topologies. 
"""
function findfragiles(sgnum::Integer, dir::AbstractString= "./symeigs/output/output")
    validfiles = String[]
    numfragiles = 0 
    numtrivial = 0 
    numnontrivial = 0 
    for file in readdir(dir)
        contains(file, "dim2-sg$(sgnum)") || continue
        contains(file, "symeigs") || continue
        push!(validfiles, file)
    end
    totalnum = length(validfiles)
    for file in validfiles
        println(file[1:end-12])
        labeledbands =try
             label_topologies(file[1:end-12], true, dir, printisbandstruct=false, verbose=false)
        catch 
            continue
        end
        isempty(labeledbands) && continue
        (sum([x[2] for x in labeledbands ] .== FRAGILE) >= 1) && (numfragiles +=1)
        (sum([x[2] for x in labeledbands ] .== NONTRIVIAL) >= 1) && (numnontrivial +=1)
        (sum([x[2] for x in labeledbands ] .== TRIVIAL) >= 1) && (numtrivial +=1)
    end
    return numfragiles/totalnum, numnontrivial/totalnum, numtrivial/totalnum
end

"""
$(TYPEDSIGNATURES)

Plot a random lattice from the bash file that defines it. Overloads PyPlot's plot method. 
"""
function plot_lattice(filepath::AbstractString; kwargs...)
    Rs, flat, isoval, _ = lattice_from_mpbparams(filepath)
    pygui(true)
    PyPlot.plot(flat, Rs; isoval=isoval, kwargs...) # via Crystalline.jl overload of PyPlot.jl's `plot`
end
