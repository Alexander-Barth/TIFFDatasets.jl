# # MODIS/Aqua+Terra Global Flood Product
#
# In this example, we will download and visualize the
# [MODIS/Aqua+Terra Global Flood Product](https://cmr.earthdata.nasa.gov/search/concepts/C2019424090-LANCEMODIS.html)
# for the Valencia region (Spain) on the 2d November, 2024.
#
# The data has been also be downloaded from [NASA World View](https://worldview.earthdata.nasa.gov). Unfortunately, only the last 7 days are visible on this sites but the download
# links below are still available. On the NASA World View page, one can click on "Add Layers", cloose "Flood (3-Day Window)" and confirm with "Add Layer".
# To download the data for a different day, click on "Data" and download the "MODIS/Aqua+Terra Global Flood Product L3 NRT 250m 3-day GeoTIFF" product, set your area of interest and go to "Download Via EarthData Search".

#
using Dates
using TIFFDatasets
using CairoMakie
using CommonDataModel: @select
using ColorSchemes
using Statistics

# The data was originally available at [nrt3.modaps.eosdis.nasa.gov](https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/allData/61/MCDWD_L3_F2_NRT/):

urls = [
    "https://data-assimilation.net/upload/Alex/TIFF/MCDWD_L3_F2_NRT.A2024306.h17v05.061.tif",
    "https://data-assimilation.net/upload/Alex/TIFF/MCDWD_L3_F2_NRT.A2024306.h18v05.061.tif",
];

# Download the data if necessary
fnames = basename.(urls)

if !all(isfile.(fnames))
    download.(urls,fnames)
end

# Inspect the first geoTIFF file
ds = TIFFDataset(fnames[1],"r")

# The size of the data array `band1`:
size(ds["band1"])

# Concatenate two geoTIFF files along the `cols` dimension (columns)
ds = TIFFDataset(fnames,aggdim = "cols");

# The band1 array is now twice as large:
size(ds["band1"])

# Make a virtual subset based on coordinate values using longitude (x) and
# latitude (y)
ds2 = @select(ds,-1 <= x <= 0.3 && 38.1 <= y <= 39.53);

# Load all the data
band1 = ds2["band1"][:,:];
lon = ds2["lon"][:,:];
lat = ds2["lat"][:,:];

# Extract the attribute "production date time" for the metadata
datetime = DateTime(ds.attrib["PRODUCTIONDATETIME"])

# Make a plot
title = "surface water " * Dates.format(datetime,"yyyy-mm-dd")

fig = Figure(size=(800,500));
ga = Axis(fig[1, 1]; title = title,
          xlabel = "longitude",
          ylabel = "latitude",
          aspect = AxisAspect(1/cosd(mean(lat))))
xlims!(ga,extrema(lon))
ylims!(ga,extrema(lat))
s = surface!(ga,lon,lat,band1,
             shading = NoShading,
             colormap = cgrad(ColorScheme([colorant"sandybrown",
                                           colorant"lightblue",
                                           colorant"darkred",
                                           colorant"red"]),
                              4, categorical = true),
             colorrange = (-0.5, 3.5))
Colorbar(fig[1,2],s,ticks = 0:3)
fig

# Load the `"description"` attribute
println(ds.attrib["description"])
