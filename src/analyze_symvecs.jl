function analyze_symvecs(sgnum::Integer, D::Integer=2, id::Integer=1, res::Integer=32, runtype::String="te"; symeigdir::String="./")
    BRS = bandreps(sgnum, D)
    symeig_filename = mpb_calcname(D, sgnum, id, res, runtype)
    println(symeig_filename)
    symeigsd, lgirsd = read_symdata(symeig_filename, sgnum=sgnum, D=D, dir=symeigdir)
    #println(symeigsd)
    #println(symeigsd["Γ"])
    symvec1 = convert.(Integer, round.(MPBUtils.symeigs2irreps(symeigsd["Γ"], get_lgirreps(sgnum, D)["Γ"], 1:4), digits=3))
    symvec2 = convert.(Integer, round.(MPBUtils.symeigs2irreps(symeigsd["K"], get_lgirreps(sgnum, D)["K"], 1:4), digits=3))
    symvec3 = convert.(Integer, round.(MPBUtils.symeigs2irreps(symeigsd["M"], get_lgirreps(sgnum, D)["M"], 1:4), digits=3))
    symvec = [symvec1..., symvec2..., symvec3...]
    #calc_detailed_topology(symvec, sgnum, D)
    calc_topology(symvec, BRS)
end