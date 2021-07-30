function findpolarizations(calcname::AbstractString, bands::UnitRange{Int64}=1:1, dir="./")
    #=readlines(dir*calcname*"-symeigs.out")
    symeigs = readdlm(dir*calcname*"-symeigs.out", ',')
    [round.(real.(parse.(ComplexF64, symeig[3:end]))) for symeig in eachrow(symeigs)]
    #readlines(dir*calcname*"-dispersion.out")
    #readdlm(dir*calcname*"-dispersion.out", ',')
    #kvs = KVec.([a[2:3] for a in eachrow(readdlm(dir*calcname*"-dispersion.out", ','))])
    =#
    println("Supplied space group is: ", parse_sgnum(calcname))
    mults = extract_multiplicities(dir*calcname, bands)
    if parse_sgnum(calcname) == 2
        Y₁ = round(Integer, mults["Y"][1])
        Γ₁ = round(Integer, mults["Γ"][1])
        A₁ = round(Integer, mults["A"][1])
        B₁ = round(Integer, mults["B"][1])
    # P₁₂ = mod(X₁ - Γ₁ - Γ₂, 2)//2
        return [mod(Y₁-Γ₁ + A₁-Γ₁, 2), mod(B₁-Γ₁ + A₁-Γ₁, 2)]
    elseif parse_sgnum(calcname) == 10
        X₁ = round(Integer, mults["X"][1])
        Γ₁ = round(Integer, mults["Γ"][1])
        Γ₂ = round(Integer, mults["Γ"][2])
        P₁₂ = mod(X₁ - Γ₁ - Γ₂, 2)
        return [P₁₂, P₁₂]    
    elseif parse_sgnum(calcname) == 16
        return 0, 0
    else 
        throw(DomainError(parse_sgnum(calcname), "Not a topological spacegroup"))
    end
end

function findcharges(calcname::AbstractString, bands::UnitRange{Int64}=1:1, dir="./")
    mults = extract_multiplicities(dir*calcname, bands)
    sgnum = parse_sgnum(calcname) == 2
    if sgnum == 2 
        Y₁ = round(Integer, mults["Y"][1])
        Γ₁ = round(Integer, mults["Γ"][1])
        A₁ = round(Integer, mults["A"][1])
        B₁ = round(Integer, mults["B"][1])
        mod((B₁-Γ₁) - (Y₁-Γ₁) + (A₁-Γ₁), 2)
    elseif sgnum == 10
        X₁ = round(Integer, mults["X"][1])
        Γ₁ = round(Integer, mults["Γ"][1])
        Γ₂ = round(Integer, mults["Γ"][2])
        mod((X₁-Γ₁-Γ₂) + 2(M₁-Γ₁) + 3(M₃-Γ₃Γ₄), 4)

    elseif sgnum == 16

    else
        throw(DomainError(sgnum, "Not a topological spacegroup"))
    end
end

"""
Returns generic information about spacegroups 2, 10, and 16
"""
function topologicalsgs()
    pgnums = (2, 6, 10)
    for pgnum in pgnums
        lgirsd = get_lgirreps(pgnum, 2)
        pgnum_kvecs = filter(isspecial, kvec.(group.(first.([lgir for (klab, lgir) in lgirsd ]))))
        println("Plane group ", pgnum)
        println("   ", pgnum_kvecs)
    end
end

function bulk_polarization_pg2(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString}) # plane group 2
    Y₁, Γ₁, A₁, B₁ = get_mult.(Ref(n), Ref(irlabs), ["Y₁", "Γ₁", "A₁", "B₁"])
    return [mod(Y₁-Γ₁ + A₁-Γ₁, 2)//2, mod(B₁-Γ₁ + A₁-Γ₁, 2)//2]
end

function bulk_polarization_pg10(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString}) # plane group 10
    X₁, Γ₁, Γ₂ = get_mult.(Ref(n), Ref(irlabs), ["X₁", "Γ₁", "Γ₂"])
    P₁₂ = mod(X₁ - Γ₁ - Γ₂, 2)//2
    return [P₁₂, P₁₂]
end

function bulk_polarization_pg16(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString})
    return [0//1, 0//1]
end

function bulk_polarization_pg16()
    return [0//1, 0//1]
end

function bulk_polarization(n::Vector{<:Integer}, sb::SymBasis)
    irlabs = sb.irlabs
    pgnum  = sb.sgnum
    pgnum == 2  && return bulk_polarization_pg2(n, irlabs)
    pgnum == 10 && return bulk_polarization_pg10(n, irlabs)
    pgnum == 16 && return bulk_polarization_pg16(n, irlabs)
    throw(DomainError(pgnum, "unsupported plane group number"))
end

function corner_charge_pg2(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString})
    Y₁, Γ₁, A₁, B₁ = get_mult.(Ref(n), Ref(irlabs), ["Y₁", "Γ₁", "A₁", "B₁"])
    return mod((B₁-Γ₁) - (Y₁-Γ₁) + (A₁-Γ₁), 2) // 4
end

function corner_charge_pg10(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString})
    X₁, Γ₁, Γ₂, Γ₃Γ₄, M₁, M₃ = get_mult.(Ref(n), Ref(irlabs), ["X₁", "Γ₁", "Γ₂", "Γ₃Γ₄", "M₁", "M₃M₄"])
    return mod((X₁-Γ₁-Γ₂) + 2(M₁-Γ₁) + 3(M₃-Γ₃Γ₄), 4) // 4
end

function corner_charge_pg16(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString})
    M₁, K₁, Γ₁, Γ₂, Γ₃Γ₅ = get_mult.(Ref(n), Ref(irlabs), ["M₁", "K₁", "Γ₁", "Γ₂", "Γ₃Γ₅"])
    return mod(6*(M₁-Γ₁-Γ₃Γ₅) + 4*(K₁-Γ₁-Γ₂), 24) // 24
end

function corner_charge(n::Vector{<:Integer}, sb::SymBasis)
    irlabs = sb.irlabs
    pgnum  = sb.sgnum
    pgnum == 2  && return corner_charge_pg2(n, irlabs)
    pgnum == 10 && return corner_charge_pg10(n, irlabs)
    pgnum == 16 && return corner_charge_pg16(n, irlabs)
    throw(DomainError(pgnum, "unsupported plane group number"))
end

function get_mult(n::Vector{<:Integer}, irlabs::Vector{<:AbstractString}, lab::AbstractString)
    idx = findfirst(==(lab), irlabs)
    idx === nothing && throw(DomainError(lab, "non-featured irrep requested"))
    return n[idx]
end

"""
$(TYPEDSIGNATURES)
"""
function find_hoti(calcname::AbstractString, has_tr::Bool=true, dir="./"; verbose::Bool=false, printisbandstruct::Bool=false)
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
    #=
    println(typeof(sb))
    println(length(ns))
    println(bands)
    =#
    polarizations = bulk_polarization.(ns, Ref(sb))
    corner_charges = corner_charge.(ns, Ref(sb))
    return bands, [float.(pol) for pol in polarizations], [float.(charge) for charge in corner_charges]
end

function find_hoti_bands(calcname::AbstractString)
    b, p, q = find_hoti(calcname)
    relevant_bands1 = (iszero).(p)
    relevant_bands2 = (!iszero).(q)
    relevant_bands = relevant_bands1 .& relevant_bands2
    return b[relevant_bands], p[relevant_bands], q[relevant_bands]
end

function extract_gaps(calcname::AbstractString, whichgap::Tuple{<:Integer, <:Integer}=(1, 2))
    #By default we look at the fundamental gap, though this can be modified by changing whichgap
    band1, band2 = whichgap
    bandsathighk = [parse.(Float64, split(line, ","))[6:end] for line in readlines(calcname)]
    #println(bandsathighk)
    min_band2 = minimum([bands[band2] for bands in bandsathighk])
    #println(min_band2)
    max_band1 = maximum([bands[band1] for bands in bandsathighk])

    return min_band2-max_band1
end