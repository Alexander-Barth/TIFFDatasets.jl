[![Build Status](https://github.com/Alexander-Barth/TIFFDatasets.jl/workflows/CI/badge.svg)](https://github.com/Alexander-Barth/TIFFDatasets.jl/actions)
[![documentation dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexander-barth.github.io/TIFFDatasets.jl/dev/)

# TIFFDatasets

TIFFDatasets pretends that a GeoTIFF file is a NetCDF file that is accessible
with the same API as [NCDatasets](https://github.com/Alexander-Barth/NCDatasets.jl):

```julia
fname = download("https://data-assimilation.net/upload/Alex/TIFF/S2_1-12-19_48MYU_0.tif")
ds = TIFFDataset(fname)
```

Output:

```
Dataset: /tmp/jl_4LSVOjATCR
Group: /

Dimensions
   cols = 256
   rows = 256

Variables
  lon   (256 × 256)
    Datatype:    Float64
    Dimensions:  cols × rows
    Attributes:
     standard_name        = longitude
     units                = degrees_east

  lat   (256 × 256)
    Datatype:    Float64
    Dimensions:  cols × rows
    Attributes:
     standard_name        = latitude
     units                = degrees_north

  x   (256)
    Datatype:    Float64
    Dimensions:  cols
    Attributes:
     standard_name        = projection_x_coordinate
     units                = m

  y   (256)
    Datatype:    Float64
    Dimensions:  rows
    Attributes:
     standard_name        = projection_y_coordinate
     units                = m

  crs
    Attributes:
     grid_mapping_name    = transverse_mercator
     longitude_of_central_meridian = 105.0
     false_easting        = 500000.0
     false_northing       = 1.0e7
     latitude_of_projection_origin = 0.0
     scale_factor_at_central_meridian = 0.9996
     longitude_of_prime_meridian = 0.0
     semi_major_axis      = 6.378137e6
     inverse_flattening   = 298.257223563
     crs_wkt              = PROJCS["WGS 84 / UTM zone 48S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32748"]]
     GeoTransform         = 706740.0 10.0 0.0 9.34096e6 0.0 -10.0

  band1   (256 × 256)
    Datatype:    Float32
    Dimensions:  cols × rows
    Attributes:
     grid_mapping         = crs

[...]
```

The dataset `ds` will also have the virtual variables `x`, `y`, `lon` and `lat`
(unless there is no projection defined in the GeoTiff file)
representing the projected coordinates and the corresponding longitude and latitude.
The projection information as the so-called [Well Known Text](https://en.wikipedia.org/wiki/Well-known_text_representation_of_coordinate_reference_systems) is available as attribute `crs_wkt` of the virtural variable `crs` following the [CF Conventions](https://cfconventions.org/cf-conventions/cf-conventions.html#use-of-the-crs-well-known-text-format).


For example, the first band and the corresponding lon/lat coordinates can be loaded with:

```julia
data = ds["band1"][:,:];
lon = ds["lon"][:,:];
lat = ds["lat"][:,:];
```
