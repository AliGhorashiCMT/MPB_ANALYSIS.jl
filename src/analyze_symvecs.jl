function analyze_symvecs(sgnum::Integer, D::Integer=2, id::Integer=1, res::Integer=32, runtype::String="te"; symeigdir::String="./")
    BRS = bandreps(sgnum, D)
    symeig_filename = mpb_calcname(D, sgnum, id, res, runtype)
    println(symeig_filename)
    symeigsd, lgirsd = read_symdata(symeig_filename, sgnum=sgnum, D=D, dir=symeigdir)
    #println(symeigsd)
    #println(symeigsd["Γ"])
    symvec1 = convert.(Integer, round.(MPBUtils.symeigs2irreps(symeigsd["Γ"], get_lgirreps(sgnum, D)["Γ"], 1:4), digits=3))
    symvec2 = convert.(Integer, round.(MPBUtils.symeigs2irreps(symeigsd["K"], get_lgirreps(sgnum, D)["K"], 1:4), digits=3))
    symvec3 = convert.(Integer, round.(MPBUtils.symeigs2irreps(symeigsd["M"], get_lgirreps(sgnum, D)["M"], 1:4), digits=3))
    symvec = [symvec1..., symvec2..., symvec3...]
    #calc_detailed_topology(symvec, sgnum, D)
    calc_topology(symvec, BRS)
end

function analyze_symvecs(calcname::String; dir="./", has_tr::Bool = true, checkfragile::Bool = true)
    sgnum = MPBUtils.parse_sgnum(calcname) 
    D = MPBUtils.parse_dim(calcname)
    D != 2 && error("Data must be in 2 dimensions")

    symbasis, bandrepset = compatibility_basis(sgnum, 2, timereversal=has_tr)

    bandirsd, lgirsd = extract_individual_multiplicities(calcname, timereversal=has_tr, latestarts = Dict{String, Int}(), dir = dir, atol=2e-2)

    length(lgirsd) ≠ length(symbasis.klabs) && error("missing k-point data")
    length(lgirsd) ≠ length(bandrepset.klabs) && error("missing k-point data")

    bandrepsetmat = matrix(bandrepset, true)
    F = smith(bandrepsetmat)

    bandrepsetmat[end, :] ≈ dim.(bandrepset) || error("Error in band filling")

    nontopologicalbasis = nontopological_basis(F, bandrepset)
    bands, nds = collect_separable(bandirsd, lgirsd, latestarts = Dict{String, Int}(),)
    isempty(bands) && error("   ... found no isolable band candidates ...")
    permd = Dict(klab => Vector{Int}(undef, length(lgirsd[klab])) for klab in symbasis.klabs)
    for klab in symbasis.klabs
        lgirs = lgirsd[klab]
        for (i, lgir) in enumerate(lgirs)
            irlab = formatirreplabel(label(lgir))
            j = findfirst(==(irlab), symbasis.irlabs)
            j === nothing && error("Could not find irrep label $irlab")
            permd[klab][i] = j
        end
    end

    μs = length.(bands)
    ns = [Vector{Int}(undef, length(first(symbasis))) for _ in 1:length(bands)]
    for (b, (nd, μ)) in enumerate(zip(nds, μs))
        for (klab, nᵏ) in nd
            permᵏ = permd[klab]
            ns[b][permᵏ] .= nᵏ
        end
        ns[b][end] = μ
    end

    band′ = 0:0
    n′ = similar(first(ns))
    is_bs = true # at every new iteration, `is_bs` effectively means `was_prev_iter_bs`
    println(bandrepsetmat)
    for (index, n) in enumerate(ns)
        println(n)
        println("n is: ", isbandstruct(n, F))
    end

    #topo = calc_detailed_topology(n′, nontopologicalbasis, bandrepsetmat)#checkfragile ? calc_detailed_topology(n′, nontopo_M, B) : calc_topology(n′, F)
end

function fix_gammairreps()

end

#Note that the symmetry eigenvalues read from read_symdata give us 20 arrays- one for each band. Then there are as many components in each 
#array as there are operations in the little group at that k point. 