using TIFFDatasets
using Test
using Downloads: download

fname = download("https://data-assimilation.net/upload/Alex/TIFF/sample_nodata.tiff")

@test isfile(fname)
ds = TIFFDataset(fname)
@test ismissing(ds["band1"][1,end])
close(ds)

ds = TIFFDataset(fname,maskingvalue=NaN32)
@test isnan(ds["band1"][1,end])
close(ds)
