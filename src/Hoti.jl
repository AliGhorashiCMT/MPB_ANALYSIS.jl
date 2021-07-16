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

function bulk_polarization_pg2(n, irlabs) # plane group 2
    Y₁, Γ₁, A₁, B₁ = get_mult.(Ref(n), Ref(irlabs), ["Y₁", "Γ₁", "A₁", "B₁"])
    return [mod(Y₁-Γ₁ + A₁-Γ₁, 2)//2, mod(B₁-Γ₁ + A₁-Γ₁, 2)//2]
end
function bulk_polarization_pg10(n, irlabs) # plane group 10
    X₁, Γ₁, Γ₂ = get_mult.(Ref(n), Ref(irlabs), ["X₁", "Γ₁", "Γ₂"])
    P₁₂ = mod(X₁ - Γ₁ - Γ₂, 2)//2
    return [P₁₂, P₁₂]
end
function bulk_polarization_pg16(n, irlabs)
    return [0//1, 0//1]
end
function bulk_polarization(n, sb)
    irlabs = sb.irlabs
    pgnum  = sb.sgnum
    pgnum == 2  && return bulk_polarization_pg2(n, irlabs)
    pgnum == 10 && return bulk_polarization_pg10(n, irlabs)
    pgnum == 16 && return bulk_polarization_pg16(n, irlabs)
    throw(DomainError(pgnum, "unsupported plane group number"))
end