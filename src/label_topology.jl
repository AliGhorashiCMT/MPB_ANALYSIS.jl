function decompose_bandset(calcname::AbstractString, runtype::AbstractString, bandrange::UnitRange{Int64}; has_tr::Bool=true, dir="./")
    sgnum = parse_sgnum(calcname)
    dim = parse_dim(calcname)
    sb, brs = compatibility_basis(sgnum, dim)
    B = matrix(brs)
    F = smith(B)

    bandirsd, lgirsd =  runtype == "tm" ? extract_individual_multiplicities(calcname, timereversal=has_tr, dir = dir, atol=2e-2) : extract_individual_multiplicities(calcname, timereversal=has_tr, latestarts = Dict{String, Int}(), dir = dir,atol=2e-2)
    runtype == "tm" && pushfirst!(bandirsd["Γ"], 1:1=>[1, zeros(length(realify(get_lgirreps(sgnum, dim)["Γ"]))-1)...])
    totaln=zeros(0)
    for (b, l) in zip(bandirsd, lgirsd)
        natkpoint=zeros(length(l[2]))
        for (bands, n) in b[2]
            minimum(bands) >= minimum(bandrange) || continue
            maximum(bands) <= maximum(bandrange) || continue    
            natkpoint += n
        end
        totaln = vcat(totaln, natkpoint)
    end
    totaln = vcat(totaln, length(bandrange))
    !isbandstruct(Int.(totaln), F) && error("You did not provide a valid band range")
    return decompose(Int.(totaln), brs)
end

"""
$(TYPEDSIGNATURES)

Returns the topological, trivial, and fragile bands
"""
function label_topologies(calcname::AbstractString, has_tr::Bool=true, dir="./"; verbose::Bool=false, printisbandstruct::Bool=false)
    sgnum = parse_sgnum(calcname)
    D = parse_dim(calcname)
    sb, brs = compatibility_basis(sgnum, D)
    mode = contains(calcname, "te") ? "te" : "tm"
    brsmat= matrix(brs, true)
    F = smith(brsmat)
    bandirsd, lgirsd = mode == "tm" ? extract_individual_multiplicities(calcname, timereversal=has_tr, dir = dir, atol=2e-2) : extract_individual_multiplicities(calcname, timereversal=has_tr, latestarts = Dict{String, Int}(), dir = dir,atol=2e-2)
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
    verbose && println("bands: ", bands)
    verbose && println("ns: ", ns)
    returnbands = []
    for (indx, band) in enumerate(bands)
        minimum(band) == minband || continue
        a, b = find_mintoposet(bands, ns, indx, F )
        verbose && println(a, "   ", b)
        printisbandstruct && println(isbandstruct(b, F))
        !isnothing(b) && push!(returnbands, [a, calc_detailed_topology(b, brs)])
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

