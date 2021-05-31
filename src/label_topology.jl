"""
$(TYPEDSIGNATURES)

"""
function label_topologies(calcname::AbstractString, has_tr::Bool=true, dir="./"; printisbandstruct::Bool=false)
    sgnum = MPBUtils.parse_sgnum(calcname)
    D = MPBUtils.parse_dim(calcname)
    sb, brs = compatibility_basis(sgnum, D)
    mode = contains(calcname, "te") ? "te" : "tm"
    brsmat= matrix(brs, true)
    F = smith(brsmat)
    nontopologicalbasis = nontopological_basis(F, brs)
    bandirsd, lgirsd =  mode == "tm" ? extract_individual_multiplicities(calcname, timereversal=has_tr, dir = dir, atol=2e-2) : extract_individual_multiplicities(calcname, timereversal=has_tr, latestarts = Dict{String, Int}(), dir = dir,atol=2e-2)
    if has_tr
        mode == "tm" && pushfirst!(bandirsd["Γ"], 1:1=>[1, zeros(length(realify(get_lgirreps(sgnum, D)["Γ"]))-1)...])
    else 
        mode == "tm" && pushfirst!(bandirsd["Γ"], 1:1=>[1, zeros(length(get_lgirreps(sgnum, D)["Γ"])-1)...])
    end
    bands, nds = collect_separable(bandirsd, lgirsd, latestarts = Dict{String, Int}())

    isempty(bands) && error("   ... found no isolatable band candidates ...")
    permd = Dict(klab => Vector{Int}(undef, length(lgirsd[klab])) for klab in sb.klabs)
    for klab in sb.klabs
        lgirs = lgirsd[klab]
        for (i, lgir) in enumerate(lgirs)
            irlab = formatirreplabel(label(lgir))
            j = findfirst(==(irlab), sb.irlabs)
            j === nothing && error("Could not find irrep label $irlab")
            permd[klab][i] = j
        end
    end
    μs = length.(bands)
    ns = [Vector{Integer}(undef, length(first(sb))) for _ in 1:length(bands)]
    for (b, (nd, μ)) in enumerate(zip(nds, μs))
        for (klab, nᵏ) in nd
            permᵏ = permd[klab]
            ns[b][permᵏ] .= nᵏ
        end
        ns[b][end] = μ
    end
    minband = 1
    println("bands: ", bands)
    println("ns: ", ns)
    returnbands = []
    for (indx, band) in enumerate(bands)
        minimum(band) == minband || continue
        a, b = find_mintoposet(bands, ns, indx, F )
        println(a, "   ", b)
        printisbandstruct && println(isbandstruct(b, F))
        #println(calc_detailed_topology(b, nontopologicalbasis, brs)) 
        !isnothing(b) && push!(returnbands, [a, calc_detailed_topology(b, nontopologicalbasis, brs)])
        minband = maximum(a) + 1
    end
    return returnbands
end

"""
$(TYPEDSIGNATURES)

Returns the minimum set of compatible bands starting form the index `idx`
"""
function find_mintoposet(bands::Vector{<:UnitRange{<:Integer}}, ns::Vector{<:Vector{<:Integer}}, idx::Integer, F::Smith{Int64,Array{Int64,2}})
    nprime = ns[idx]
    bandsprime = minimum(bands[idx]):maximum(bands[idx])
    for (band, n) in zip(bands[idx+1:end], ns[idx+1:end])
        isbandstruct(nprime, F) && break
        nprime = nprime + n
        bandsprime = minimum(bands[idx]):maximum(band)
    end
    isbandstruct(nprime, F) ? (bandsprime, nprime) : (bandsprime, nothing)
end

#