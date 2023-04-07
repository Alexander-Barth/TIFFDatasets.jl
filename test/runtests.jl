using TIFFDatasets
using Test
using Downloads: download
using NCDatasets

function test_show(ds,expect)
    io = IOBuffer()
    show(io,ds)
    out = String(take!(io))
    @test occursin(expect,out)
end

@testset "TIFFDatasets.jl" begin
    fname = download("https://data-assimilation.net/upload/Alex/TIFF/S2_1-12-19_48MYU_0.tif")

    @test isfile(fname)

    ds = TIFFDataset(fname,catbands = false)
    test_show(ds,"Variables")
    test_show(ds["band1"],"grid_mapping")
    test_show(ds["x"],"projection")

    @test length(keys(ds)) == 16
    # reference values from:
    # gdal_translate -of netCDF input.tif output.nc
    # Note: gdal_translate flips the y-axis, which we don't do here
    @test ds["x"][1] == 706745
    @test ds["x"][2] == 706755
    @test ds["y"][end] == 9338405
    @test ds["y"][end-1] == 9338415
    @test ds["band1"][1,end] == 0.09404379f0

    test_show(ds["band1"],"grid_mapping")
    test_show(ds["band2"],"grid_mapping")

    @test_throws KeyError ds["does_not_exists"]
    @test ds["crs"][:] isa Integer

    # convert to NetCDF
    ncfile = tempname()
    @time NCDatasets.write(ncfile,ds)
    ncds = NCDataset(ncfile)

    geotiff_band1 = ds["band1"][:,:]
    nc_band1 = ncds["band1"][:,:]
    @test geotiff_band1 == nc_band1
    @test collect(ncds["band1"].attrib) == collect(ds["band1"].attrib)

    # concatenate bands

    ds = TIFFDataset(fname,catbands = true)
    @test haskey(ds.dim,"bands")
    @test "band" in keys(ds)

    data = ds["band"][:,:,:]
    @test data[1,1,1] == ds["band"][1,1,1]
    @test data[1,1,1:2] == ds["band"][1,1,1:2]
    @test data[:,:,1:2] == ds["band"][:,:,1:2]

    test_show(ds["band"],"grid_mapping")
    test_show(ds,"Variables")

    @test ds["band"][1,1,1] isa Float32

    @test ndims(ds["x"]) == 1
    @test ndims(ds["y"]) == 1

    @test ds["lon"][1,1] isa Float64
    @test -180 <= ds["lon"][1,1] <= 360
    @test -90 <= ds["lat"][1,1] <= 90
    @test eltype(ds["lon"]) == Float64

    test_show(ds["crs"],"crs_wkt")
    @test ds["crs"].attrib["crs_wkt"] isa String

end

#=
fname = download("https://data-assimilation.net/upload/Alex/TIFF/S2_1-12-19_48MYU_0.tif")
@time NCDatasets.write("tmp.nc",TIFFDataset(fname))
=#
