function analyze_symvecs(sgnum::Integer, D::Integer=2, id::Integer=1, res::Integer=32, runtype::String="te"; symeigdir::String="./")
    BRS = bandreps(sgnum, D)
    symeig_filename = mpb_calcname(D, sgnum, id, res, runtype)
    println(symeig_filename)
    symeigsd, lgirsd = read_symdata(symeig_filename, sgnum=sgnum, D=D, dir=symeigdir)
    #println(symeigsd)
    #println(symeigsd["Γ"])
    println(MPBUtils.symeigs2irreps(symeigsd["Γ"], get_lgirreps(sgnum, D)["Γ"], 1:2))
end