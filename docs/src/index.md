
## TIFFDatasets


Documentation for [TIFFDatasets.jl](https://github.com/Alexander-Barth/TIFFDatasets.jl), a Julia package for loading geoTIFF files
using [ArchGDAL](https://github.com/yeesian/ArchGDAL.jl).

TIFFDatasets.jl implements for the TIFF format the interface defined
in [CommonDataModel.jl](https://github.com/JuliaGeo/CommonDataModel.jl).
The functions defined by CommonDataModel.jl are also available for TIFF data, including:
* virtually concatenating multiple files along a given dimension
* create a virtual subset (`view`) by indices or by values of coordinate variables (`CommonDataModel.select`, `CommonDataModel.@select`)
* group, map and reduce (with `mean`, standard deviation `std`, ...) a variable (`CommonDataModel.groupby`, `CommonDataModel.@groupby`) and rolling reductions like running means `CommonDataModel.rolling`).

However, TIFFDatasets cannot modify or create geoTIFF files.

## Installation

Inside the Julia shell, you can download and install TIFFDatasets using the following commands:

```julia
using Pkg
Pkg.add("TIFFDatasets")
```

Or by typing `]add TIFFDatasets` using the package manager REPL-mode.

## API reference

```@autodocs
Modules = [TIFFDatasets]
```
