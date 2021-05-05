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

function plot_topologybands(sgnum::Integer, id::Integer, runtype::String; dim::Integer=2, res::Integer=32, dispersiondir::String="./", symeigdir::String="./")
    dispersion_filename = mpb_calcname(dim, sgnum, id, res, runtype)*"-dispersion.out"
    symeig_filename = mpb_calcname(dim, sgnum, id, res, runtype)*"-symeigs.out"
    for line in readlines(dispersiondir*dispersion_filename)
        println(parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end])
    end
    #Extract symmetry data 
    
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