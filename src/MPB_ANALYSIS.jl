module MPB_ANALYSIS
#Dependencies
using Pkg, Plots, MPBUtils, Crystalline, DocStringExtensions, SymmetryBases, PyPlot, DelimitedFiles#, Brillouin

using Crystalline: formatirreplabel, symvec2string, label

using MPBUtils: parse_dim, parse_sgnum

include("plot_mpb.jl")
export plot_mpb, plot_bothmodes

include("plot_topologybands.jl")
export plot_topologybands, plot_manytopologybands

include("analyze_symvecs.jl")
export analyze_symvecs

include("label_topology.jl")
export label_topologies

include("dispersionpaths.jl")

include("spacegroup_analysis.jl")
export findfragiles

include("Hoti.jl")
export findpolarizations, find_hoti

export dispersionpaths

function __init__()
    println("Running Init")
    #ENV["PYTHON"]="/usr/local/bin/python3"
	#Pkg.build("PyCall")
end

end # module
