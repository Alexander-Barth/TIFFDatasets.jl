using TIFFDatasets
using Test
using Downloads: download

fname = download("https://data-assimilation.net/upload/Alex/TIFF/MCDWD_L3_F2_NRT.A2024306.h17v05.061.tif")

ds1 = TIFFDataset(fname)
ds = TIFFDataset([fname,fname],aggdim="cols")
@test size(ds["band1"]) == (9600, 4800)
@test isequal(ds["band1"][1:4800,:], ds1["band1"][:,:])

ds = TIFFDataset([fname,fname],aggdim="time",isnewdim=true)
@test size(ds["band1"]) == (4800, 4800, 2)
