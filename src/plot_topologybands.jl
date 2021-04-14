function plot_topologybands(filename::String)
    Bands = Vector{Vector{Float64}}()
    reshapedBands = Array{Float64, 2}(undef, (length(Bands), length(Bands[1])))
    try
        for line in readlines(filename)
            push!(Bands, parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end], "\n")
        end
    catch e
        for line in readlines(filename*".out")
            push!(Bands, parse.(Ref(Float64), string.(split(replace(line, "," => " "))))[6:end], "\n")
        end
    end
    for i in 1:length(Bands)
        reshapedBands[i, :] = Bands[i] ##Note that each row now corresponds to a point in reciprocal space, as desired
    end
    plot(reshapedBands)
end
