function plot_mpb(filename::AbstractString, mode::AbstractString, plotorno::Bool=true; kwargs...)
    Bands = Vector{Vector{Float64}}()
    for line in readlines(filename)
        contains(line, "k index") && continue
        contains(line, mode*"freqs") || continue
        line = replace(line, "," => "")
        kpoint_freqs = split(line)[7:end] # Don't look at kpoint indices
        kpoint_freqs = string.(kpoint_freqs) 
        kpoint_freqs = parse.(Float64, kpoint_freqs) 
        push!(Bands, kpoint_freqs)
    end
    numkpoints = length(Bands)
    numbands = length(first(Bands))
    reshapedBands = Array{Float64, 2}(undef, (numkpoints, numbands))
    for k in 1:numkpoints
       reshapedBands[k, :] = Bands[k]
    end
    plotorno && (plot(reshapedBands; kwargs...); ylabel("Frequency"); xlabel("Wavevector"); xticks([]); title("Dispersion"))
    return reshapedBands
end

function plot_dispersion(filename::AbstractString)
    bands = readdlm(filename, ',')[:, 6:end]
    bands_plot = plot(bands, legend=false, linewidth=5)
    title!(filename, titlefontsize=20)
end

function plot_bothmodes(filename::AbstractString; kwargs...)
    plot_mpb(filename, "te"; color="red", label="TE Modes", kwargs...)
    plot_mpb(filename, "tm"; color="green", label="TM Modes",  kwargs...)
end