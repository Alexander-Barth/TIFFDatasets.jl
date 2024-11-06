using TIFFDatasets
using Test

@testset "TIFFDatasets.jl" begin
    include("test_cf.jl")
    include("test_nodata.jl")
    include("test_attrib.jl")
end

#=
fname = download("https://data-assimilation.net/upload/Alex/TIFF/S2_1-12-19_48MYU_0.tif")
@time NCDatasets.write("tmp.nc",TIFFDataset(fname))
=#
