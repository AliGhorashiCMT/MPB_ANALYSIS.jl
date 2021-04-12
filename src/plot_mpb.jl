function plot_mpb(filename::String, mode::String, nbands::Int)
    Bands = Vector{Vector{Float64}}()
    for line in readlines(filename)
        if contains(line, mode*"freqs")
            try
                push!(Bands, parse.(Float64, string.(chop.(string.(split(line))))[end-nbands+1:end]))
            catch 
            end
        end
    end
    numkpoints = length(Bands)
    numbands = length(Bands[1])
    reshapedBands = Array{Float64, 2}(undef, (numkpoints, numbands))
    for k in 1:numkpoints
       reshapedBands[k, :] = Bands[k]
    end
    #return Bands 
    return reshapedBands
end

function plot_bothmodes(filename::String, nbands::Int)
    plot(plot_mpb(filename, "te", nbands), linewidth=5, color="red", legend=false, xticks=false)
    plot!(plot_mpb(filename, "tm", nbands), linewidth=5, color="blue", legend=false, xticks=false, ylabel="Frequency")
end