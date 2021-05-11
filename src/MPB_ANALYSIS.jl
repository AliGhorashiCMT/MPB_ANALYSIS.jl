module MPB_ANALYSIS

#Dependencies
using Pkg, Plots, MPBUtils, Crystalline, DocStringExtensions, SymmetryBases

using Crystalline: formatirreplabel, symvec2string, label

include("plot_mpb.jl")
export plot_mpb, plot_bothmodes

include("plot_topologybands.jl")
export plot_topologybands, plot_manytopologybands

include("analyze_symvecs.jl")
export analyze_symvecs

greet() = print("Hello World!")

function __init__()
    println("Running Init")
    #ENV["PYTHON"]="/usr/local/bin/python3"
	#Pkg.build("PyCall")
end

end # module
