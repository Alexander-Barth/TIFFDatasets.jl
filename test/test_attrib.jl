using TIFFDatasets
using Test
using Downloads: download

fname = download("https://data-assimilation.net/upload/Alex/TIFF/MCDWD_L3_F2_NRT.A2024306.h17v05.061.tif")

ds = TIFFDataset(fname)

@test haskey(ds.attrib,"RANGEBEGINNINGDATE")
@test haskey(ds["band1"].attrib,"valid_range")
