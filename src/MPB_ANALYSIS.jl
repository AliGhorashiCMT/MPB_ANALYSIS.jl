module MPB_ANALYSIS
#Dependencies
using Pkg, MPBUtils, Crystalline, DocStringExtensions, SymmetryBases, PyPlot, PyCall, DelimitedFiles#, Brillouin
using Crystalline: formatirreplabel, symvec2string, label
using MPBUtils: parse_dim, parse_sgnum

include("plot_mpb.jl")
export plot_mpb, plot_bothmodes, plot_dispersion

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
export findpolarizations, find_hoti, find_hoti_bands

export dispersionpaths

function __init__()
    println("Running Init")
    # The Python environment should be linked to pynormaliz for symmetry calculations. 
    #ENV["PYTHON"]="/usr/local/bin/python3"
	#Pkg.build("PyCall")
end

end # module
