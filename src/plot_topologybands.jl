function plot_topologybands(filename::String, highsymmetrylabels::Array{<:AbstractString})
    Bands = Vector{Vector{Float64}}()
    try
        for line in readlines(filename)
            push!(Bands, parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end])
        end
    catch e
        for line in readlines(filename*".out")
            push!(Bands, parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end])
        end
    end
    reshapedBands = Array{Float64, 2}(undef, (length(Bands), length(Bands[1])))
    for i in 1:length(Bands)
        reshapedBands[i, :] = Bands[i] ##Note that each row now corresponds to a point in reciprocal space, as desired
    end
    plot(reshapedBands, legend=false, xticks=false, size=(1000, 500), linestyle=:dash, color="pink", linewidth=5)
    println(length(Bands))
    highsymmetrygap = Int((length(Bands)-1)/(length(highsymmetrylabels)-1))
    highsymmetryxcoords = collect(1:highsymmetrygap:length(Bands))
    for (highsymmetrylabel, highsymmetryxcoord) in zip(highsymmetrylabels, highsymmetryxcoords)
        irlabels = String.(split(replace(replace(highsymmetrylabel, "+" => " "), "," =>" ")        ))
        for (energy, irreps) in zip(Bands[highsymmetryxcoord], irlabels)
            display(annotate!(highsymmetryxcoord, energy, text(irreps, 10)))
        end
    end
end

"""
$(TYPEDSIGNATURES)
"""
function plot_topologybands(sgnum::Integer, id::Integer, runtype::String; dim::Integer=2, res::Integer=32, dispersiondir::String="./", symeigdir::String="./", kwargs...)
    dispersion_filename = mpb_calcname(dim, sgnum, id, res, runtype)*"-dispersion.out"
    symeig_filename = mpb_calcname(dim, sgnum, id, res, runtype)*"-symeigs.out"
    println(dispersion_filename)
    athighsymmetry = Vector{Integer}()
    highsymmetrykvecs = bandreps(sgnum, dim).kvs
    highsymmetryklabs = bandreps(sgnum, dim).klabs
    xticks = Vector{Integer}()
    xticklabels = Vector{String}()
    Bands = Vector{Vector{Float64}}()
    for (index, line) in enumerate(readlines(dispersiondir*dispersion_filename))
        push!(Bands, parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end])
        sum([isapprox(parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[2:2+dim-1] , highsymmetrykvec.cnst, atol=1e-3) for highsymmetrykvec in highsymmetrykvecs]) == 1 && push!(athighsymmetry, index )
    end
    println(athighsymmetry)
    reshapedBands = Array{Float64, 2}(undef, (length(Bands), length(Bands[1])))
    for i in 1:length(Bands)
        reshapedBands[i, :] = Bands[i] ##Note that each row now corresponds to a point in reciprocal space, as desired
    end
    display(plot(reshapedBands, linewidth=5, legend=false,  size=(1000, 500), xtickfontsize=20; kwargs...))
    for (index, line) in enumerate(readlines(dispersiondir*dispersion_filename))
        lab = index in athighsymmetry ? highsymmetryklabs[findfirst(x -> isapprox(x.cnst, parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[2:2+dim-1], atol=1e-3), highsymmetrykvecs )] : nothing
        !isnothing(lab) && (push!(xticks, index); push!(xticklabels, lab)) #isplay(annotate!(index, 0, text(lab)))
    end
    display(xticks!(xticks, xticklabels))
    bandirsd, lgirsd = extract_individual_multiplicities(mpb_calcname(dim, sgnum, id, res, runtype),
    timereversal=true, dir = symeigdir, atol=1e-1, latestarts= Dict{String,Int}())
    irlabs = Dict(klab => formatirreplabel.(label.(lgirs)) for (klab, lgirs) in lgirsd)
    labeldict = Dict(klab => [bands => symvec2string(n, irlabs[klab]; braces=false) for (bands, n) in bandirs]         for (klab, bandirs) in bandirsd)
    for (tick, label) in zip(xticks, xticklabels)
        groupings = labeldict[label]
        for x in groupings
            println(x.first)
        end
        println("\n")
        for i in 1:length(Bands[1])
            energy = reshapedBands[tick, i]
            whichgrouping = findfirst(x -> i in x.first, groupings)
            #println(whichgrouping, i)
            try
                bandlabel = groupings[whichgrouping].second
                display(annotate!(tick, energy, text(bandlabel, 15)))
            catch
            end
        end
    end    
end

function plot_topologybands(filename::String, highsymmetrylabels::String)
    splitted_array_string = split.(replace(replace(replace(highsymmetrylabels, "["=>""), "]"=>""), ","=>" "))
    plot_topologybands(filename, splitted_array_string)
end

function plot_manytopologybands(directory::String)
    plot_array = Vector{Tuple{String, Plots.Plot}}()
    for filename in readdir(directory)
        push!(plot_array, (filename, plot_topologybands(filename) ))
    end
    return plot_array
end