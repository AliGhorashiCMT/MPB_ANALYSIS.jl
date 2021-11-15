function plot_mpb(filename::AbstractString, mode::AbstractString, nbands::Integer)
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
    return reshapedBands
end

function plot_dispersion(filename::AbstractString)
    bands = readdlm(filename, ',')[:, 6:end]
    bands_plot = Plots.plot(bands, legend=false, linewidth=5, )
    title!(filename, titlefontsize=20)
    display(bands_plot)
end

function plot_bothmodes(filename::AbstractString, nbands::Integer)
    plot(plot_mpb(filename, "te", nbands), linewidth=5, color="red", legend=true, xticks=false, )
    plot!(plot_mpb(filename, "tm", nbands), linewidth=5, color="blue", legend=false, xticks=false, ylabel="Frequency", label = "TM FREQS")
end