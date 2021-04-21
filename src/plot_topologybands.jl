function plot_topologybands(filename::String)
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
end

function plot_manytopologybands(directory::String)
    plot_array = Vector{Tuple{String, Plots.Plot}}()
    for filename in readdir(directory)
        push!(plot_array, (filename, plot_topologybands(filename) ))
    end
    return plot_array
end