function plot_topologybands(filename::AbstractString, highsymmetrylabels::Array{<:AbstractString})
    Bands = Vector{Vector{Float64}}()
    try
        for line in readlines(filename)
            push!(Bands, parsebands(line))
        end
    catch 
        for line in readlines(filename*".out")
            push!(Bands, parsebands(line) )
        end
    end
    reshapedBands = Array{Float64, 2}(undef, (length(Bands), length(Bands[1])))
    for (i, b) in enumerate(Bands)
        reshapedBands[i, :] = b ##Note that each row now corresponds to a point in reciprocal space, as desired
    end
    plot(reshapedBands, legend=false, xticks=false, size=(1000, 500), linestyle=:dash, color="pink", linewidth=5)
    println(length(Bands))
    highsymmetrygap = Int((length(Bands)-1)/(length(highsymmetrylabels)-1))
    highsymmetryxcoords = collect(1:highsymmetrygap:length(Bands))
    for (highsymmetrylabel, highsymmetryxcoord) in zip(highsymmetrylabels, highsymmetryxcoords)
        irlabels = String.(split(replace(replace(highsymmetrylabel, "+" => " "), "," =>" ")        ))
        for (energy, irreps) in zip(Bands[highsymmetryxcoord], irlabels)
            text(highsymmetryxcoord, energy, irreps, 10)
        end
    end
end

"""
$(TYPEDSIGNATURES)
Returns the band eigenvalues from a string corresponding to a line from an MPB output
"""
function parsebands(line::AbstractString)
    return parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end]
end

"""
$(TYPEDSIGNATURES)
Returns the kvector from a string corresponding to a line from an MPB output. 
"""
function parsekvec(line::AbstractString, dim::Integer=2)
    return parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[2:2+dim-1] #Note that dimensionality of kvector is dependent on dimensionality of simulation
end

"""
$(TYPEDSIGNATURES)
Plots MPB bands and overlays both the irreps and the topological nature of the bands by denoting the color.
"""
function plot_topologybands(sgnum::Integer, id::Integer, runtype::AbstractString; dim::Integer=2, 
    res::Integer=32, dispersiondir::String="./", symeigdir::AbstractString="./", verbose::Bool=false, kwargs...)
    whichtopologiesVEC = Symbol[] #A vector of symbols denoting the colors of each bands (corresponding to their topological classfication)

    calcname = mpb_calcname(dim, sgnum, id, res, runtype)
    dispersion_filename = calcname*"-dispersion.out" 
    
    athighsymmetry = Vector{Integer}() #A vector of indices that will correspond to the points in k space corresponding to high symmetry k points
    highsymmetrykvecs = bandreps(sgnum, dim).kvs #High symmetry kvectors
    highsymmetryklabs = bandreps(sgnum, dim).klabs #And their corresponding labels
    
    xticks = Vector{Integer}()
    xticklabels = Vector{String}()
    
    Bands = Vector{Vector{Float64}}()
    Frag, Nontop, Top, Unknown = :Red, :Blue, :Black, :Pink 
    println("Coloring scheme is: ", "Fragile: $(string(Frag)),  Nontopological: $(string(Nontop)),  Topological: $(string(Top)), and Unknown: $(string(Unknown))\n\n")

    whichtopologies = label_topologies(calcname, true, symeigdir, verbose=verbose)
    verbose && println(whichtopologies...)
    
    for (index, line) in enumerate(readlines(dispersiondir*dispersion_filename))
        push!(Bands, parsebands(line))
        sum([isapprox(parsekvec(line, dim), highsymmetrykvec.cnst, atol=1e-3) for highsymmetrykvec in highsymmetrykvecs]) == 1 && push!(athighsymmetry, index )
    end
    
    for (i, _) in enumerate(first(Bands))
        ind = findfirst(x-> i in first(x), whichtopologies)
        isnothing(ind) && (push!(whichtopologiesVEC, Unknown); continue)
        (whichtopologies[ind][2] == TRIVIAL) && (push!(whichtopologiesVEC, Nontop); continue)
        (whichtopologies[ind][2] == NONTRIVIAL) && (push!(whichtopologiesVEC, Top); continue)
        (whichtopologies[ind][2] == FRAGILE) && (push!(whichtopologiesVEC, Frag); continue)
    end
    reshapedBands = zeros(length(Bands), length(first(Bands)))
    for (i, b) in enumerate(Bands)
        reshapedBands[i, :] = b #Note that each row now corresponds to a point in reciprocal space, as desired
    end
    for (c, b) in zip(whichtopologiesVEC, eachcol(reshapedBands))
        plot(b, color=c)
    end
    
    for (index, line) in enumerate(readlines(dispersiondir*dispersion_filename))
        lab = index in athighsymmetry ? highsymmetryklabs[findfirst(x -> isapprox(x.cnst, parsekvec(line, dim), atol=1e-3), highsymmetrykvecs )] : nothing
        !isnothing(lab) && (push!(xticks, index); push!(xticklabels, lab)) 
    end
    PyPlot.xticks(xticks .- 1, xticklabels)
    has_tr = true
    bandirsd, lgirsd =  runtype == "tm" ? extract_individual_multiplicities(calcname, timereversal=has_tr, dir = symeigdir, atol=2e-2) : extract_individual_multiplicities(calcname, timereversal=has_tr, latestarts = Dict{String, Int}(), dir = symeigdir,atol=2e-2)
    runtype == "tm" && pushfirst!(bandirsd["??"], 1:1=>[1, zeros(length(realify(get_lgirreps(sgnum, dim)["??"]))-1)...])
    irlabs = Dict(klab => formatirreplabel.(label.(lgirs)) for (klab, lgirs) in lgirsd)
    labeldict = Dict(klab => [bands => symvec2string(n, irlabs[klab]; braces=false) for (bands, n) in bandirs] for (klab, bandirs) in bandirsd)
    for (tick, label) in zip(xticks, xticklabels)
        groupings = labeldict[label]
        for (i, _) in enumerate(first(Bands))
            energy = reshapedBands[tick, i]
            whichgrouping = findfirst(x -> i in x.first, groupings)
            try
                bandlabel = groupings[whichgrouping].second
                text(tick, energy, bandlabel)
            catch
            end
        end
    end    
    title("Spacegroup $(sgnum) Type: $(runtype) ID: $(id)")
end

function plot_topologybands(filename::AbstractString, highsymmetrylabels::AbstractString)
    splitted_array_string = split.(replace(replace(replace(highsymmetrylabels, "["=>""), "]"=>""), ","=>" "))
    plot_topologybands(filename, splitted_array_string)
end

function plot_manytopologybands(directory::AbstractString)
    plot_array = Vector{Tuple{String, Plots.Plot}}()
    for filename in readdir(directory)
        push!(plot_array, (filename, plot_topologybands(filename)))
    end
    return plot_array
end