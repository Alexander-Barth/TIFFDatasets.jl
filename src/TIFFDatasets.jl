module TIFFDatasets

import ArchGDAL
import Base: size, keys, getindex
import CommonDataModel:
    AbstractDataset,
    AbstractVariable,
    attrib,
    attribnames,
    cfvariable,
    dataset,
    dim,
    dimnames,
    maskingvalue,
    MFDataset,
    name,
    path,
    variable,
    varnames

using DataStructures: OrderedDict
import GDAL
import Proj
import JSON3
import DiskArrays

const SymbolOrString = Union{AbstractString, Symbol}

struct TIFFDataset{T,N,Ttrans,Tmaskingvalue} <: AbstractDataset
    fname::String
    dataset::ArchGDAL.IDataset
    dimnames::NTuple{3,Symbol}
    varname::Symbol
    nraster::Int
    width::Int
    height::Int
    trans::Ttrans
    crs_wkt::String
    geotransform::Vector{Float64}
    dim::OrderedDict{String,Int64}
    attrib::OrderedDict{String,Any}
    maskingvalue::Tmaskingvalue
    metadata::OrderedDict{SubString{String}, SubString{String}}
    metadata_domain::OrderedDict{String, OrderedDict{SubString{String}, SubString{String}}}
end

const Dataset = TIFFDataset
## DiskArray
struct Variable{T,N} <: AbstractVariable{T,N} #AbstractArray{T,N}
    parent::TIFFDataset{T,N}
    index::Int64
    attrib::OrderedDict{String,Any}
end

struct Coord{T,N,Td,Nd} <: AbstractVariable{T,N}
    parent::TIFFDataset{Td,Nd}
    name::Symbol
    attrib::OrderedDict{String,Any}
end

struct CRS{Td,Nd} <: AbstractVariable{Int32,0}
    parent::TIFFDataset{Td,Nd}
    attrib::OrderedDict{String,Any}
end



attribnames(ds::Union{Dataset,Variable,CRS,Coord}) = keys(ds.attrib)
attrib(ds::Union{Dataset,Variable,CRS,Coord},name::SymbolOrString) =
    ds.attrib[String(name)]

dimnames(ds::Dataset) = String.(keys(ds.dim))
dim(ds::Union{TIFFDataset,Variable,CRS,Coord},name::SymbolOrString) =
    ds.dim[String(name)]

path(ds::Dataset) = ds.fname
maskingvalue(ds::Dataset) = ds.maskingvalue

aligned_grid(ds::Dataset) = (ds.geotransform[3] == 0) && (ds.geotransform[5] == 0)

dataset(v::Union{Variable,CRS,Coord}) = v.parent

function cf_crs_attributes!(crs_wkt,attrib; geotransform = nothing)
    projparam(name; default = 0.) = GDAL.osrgetprojparm(gdal_proj.ptr,name,default,C_NULL)

    gdal_proj = ArchGDAL.importWKT(crs_wkt)
    gdal_projection = GDAL.osrgetattrvalue(gdal_proj.ptr,"PROJECTION",0)
    grid_mapping_name = nothing

    #pp = Proj.proj_create_from_wkt(crs_wkt)
    #pj = JSON3.read(Proj.proj_as_projjson(pp))
    #@show pj[:coordinate_system][:axis][1][:unit]

    if !isnothing(gdal_projection)
        grid_mapping_name = lowercase(gdal_projection)
        attrib["grid_mapping_name"] = grid_mapping_name
    end

    if grid_mapping_name == "transverse_mercator"
        # mapping from
        # https://github.com/OSGeo/gdal/blob/9b313f502304894f6c9be71ec39b26c2baeaaa87/frmts/netcdf/netcdfdataset.h#L592

        attrib["longitude_of_central_meridian"] = projparam("central_meridian")
        attrib["false_easting"] = projparam("false_easting")
        attrib["false_northing"] = projparam("false_northing")
        attrib["latitude_of_projection_origin"] = projparam("latitude_of_origin")
        attrib["scale_factor_at_central_meridian"] = projparam("scale_factor")
    else
        # You are welcome to contribute here! :-)
        # Have a look at GDAL's function NCDFWriteSRSVariable
        # https://github.com/OSGeo/gdal/blob/9b313f502304894f6c9be71ec39b26c2baeaaa87/frmts/netcdf/netcdfdataset.cpp#L5281
        #
        # the following files in general
        # frmts/netcdf/netcdfdataset.cpp
        # frmts/netcdf/netcdfdataset.h
    end

    attrib["longitude_of_prime_meridian"] = GDAL.osrgetprimemeridian(gdal_proj.ptr,C_NULL)
    attrib["semi_major_axis"] = GDAL.osrgetsemimajor(gdal_proj.ptr,C_NULL)
    attrib["inverse_flattening"] = GDAL.osrgetinvflattening(gdal_proj.ptr,C_NULL)
    attrib["crs_wkt"] = crs_wkt

    if !isnothing(geotransform)
        # not in NetCDF CF, but created by
        # gdal_translate -of netCDF input.tif output.nc
        attrib["GeoTransform"] = join(string.(geotransform),' ')
    end

    return attrib
end


# return the metadata as a ordered dictionary of key-value pairs
function _metadata(dataset,domain="")
    key_values = (split(item,'=',limit=2) for item in ArchGDAL.metadata(dataset; domain))
    return OrderedDict(key_values)
end

"""
    ds = TIFFDataset(fname::AbstractString,mode = "r"; varname = "band",
                     projection = "EPSG:4326",
                     catbands = false,
                     dimnames = ("cols","rows","bands"))

Create a data set structure from the GeoTIFF file name `fname`
using [ArchGDAL](https://github.com/yeesian/ArchGDAL.jl). The data variable will be called
`varname` (followed by the band number).
The dataset `ds` will also have the virtual variables `x`, `y` and `lon`, `lat`
representing the projected coordinates (as defined in the GeoTIFF file) and the
corresponding longitude/latitude (using the projection as defined
by the `projection` parameter). The names of the dimensions
can be adapted using the parameter `dimnames`. The type of projection can be
checked with the attributes (`grid_mapping_name`, `longitude_of_central_meridian`,
`false_easting`, `false_northing`, `latitude_of_projection_origin`, `scale_factor_at_central_meridian,
longitude_of_prime_meridian`, `semi_major_axis`, `inverse_flattening`, ...
depending on the projection used) of the virtual variable `crs`.
The attribute `crs_wkt` describes the coordinate reference system as
so-called ["Well-Known Text"](https://en.wikipedia.org/wiki/Well-known_text_representation_of_coordinate_reference_systems)


In addition to the CF conventions, there is also the attribute `GeoTransform` [1]
with six numbers describing a linear mapping between indices and the projected coordinates.
For the point `(i,j)` the projected coordinates are computed as:

```julia
shift = 0.5 # center
x_geo = GeoTransform[1] + (i-shift) * GeoTransform[2] + (j-shift) * GeoTransform[3]
y_geo = GeoTransform[4] + (i-shift) * GeoTransform[5] + (j-shift) * GeoTransform[6]
```

Setting `catbands` to true, will concatenate all bands into a single
variable called `varname`. If `catbands` is true, it is the user's
responsibility to check that all bands have the same no-data value,
offset and scale factor.
Check with the command line tool `gdalinfo filename` or use ArchGDAL to inspect the file.
If this is not the case, `catbands = true` should not be used.


Example:

```julia
using TIFFDatasets
fname = download("https://data-assimilation.net/upload/Alex/TIFF/S2_1-12-19_48MYU_0.tif")
ds = TIFFDataset(fname)
band1 = ds["band1"][:,:];
longitude = ds["lon"][:,:];
latitude = ds["lat"][:,:];
# x,y native projection
x = ds["x"][:,:];
y = ds["y"][:,:];
# coordinate reference system as "Well-Known Text"
crs_wkt = ds["crs"].attrib["crs_wkt"];
```

For performance reasons, one should avoid loading a large number of individual
elements (or using function that do so, like `sum`, `mean` ...)
and rather load the whole array. For example:

```
julia> @time sum(ds["band1"]) # slow
  0.062227 seconds (721.65 k allocations: 22.045 MiB)
6403.1743f0

julia> @time sum(ds["band1"][:,:]) # much faster
  0.000381 seconds (765 allocations: 302.781 KiB)
6403.5527f0
```

The dataset can be manipulated using the functions from
[CommonDataModel.jl](https://github.com/JuliaGeo/CommonDataModel.jl).
Writing to a file is currently not supported (pull-requests are welcome!).

Reference:

* [1] [GDAL GeoTransform](https://gdal.org/tutorials/geotransforms_tut.html#transformation-from-image-coordinate-space-to-georeferenced-coordinate-space)

"""
function TIFFDataset(fname::AbstractString, mode = "r"; varname = "band",
                     projection = "EPSG:4326",
                     catbands = false,
                     dimnames = ("cols","rows","bands"),
                     attrib = OrderedDict{String,Any}(
                         "Conventions" => "CF-1.8",
                     ),
                     maskingvalue = missing,
                     )

    if mode != "r"
        error("""only mode "r" is currently supported""")
    end

    dataset = ArchGDAL.read(fname)
    width = ArchGDAL.width(dataset)
    height = ArchGDAL.height(dataset)
    nraster = ArchGDAL.nraster(dataset)
    proj = ArchGDAL.getproj(dataset)

    if proj == ""
        @debug "no projection defined in file '$fname'. If you have GDAL installed, you can check with: gdalinfo '$fname'"
        trans = nothing
        crs_wkt = ""
    else
        crs_wkt = ArchGDAL.importWKT(proj)
        trans = Proj.Transformation(ArchGDAL.toPROJ4(ArchGDAL.importWKT(proj)), projection)
    end

    geotransform = ArchGDAL.getgeotransform(dataset)
    T = ArchGDAL.pixeltype(dataset)

    dim = OrderedDict{String,Int64}(
        dimnames[1] => width,
        dimnames[2] => height,
    )

    N = (catbands ? 3 : 2)

    if catbands
        dim[dimnames[3]] = nraster
    end


    metadata = _metadata(dataset)
    domains = filter(!=(""), ArchGDAL.metadatadomainlist(dataset))
    domain_metadata = OrderedDict((domain => _metadata(dataset,domain)) for domain in domains)

    for (k,v) in metadata
        # ignore attributes specific to a variable
        if !(k in ("add_offset","scale_factor","_FillValue","long_name","valid_range"))
            attrib[k] = v
        end
    end

    return TIFFDataset{T,N,typeof(trans),typeof(maskingvalue)}(
        fname,
        dataset,
        Symbol.(dimnames),
        Symbol(varname),
        nraster,
        width,
        height,
        trans,
        proj,
        geotransform,
        dim,
        attrib,
        maskingvalue,
        metadata,
        domain_metadata
    )
end

geo_referenced(ds) = !isnothing(ds.trans)

function varnames(ds::Dataset{T,N}) where {T,N}
    return [(geo_referenced(ds) ? ("lon","lat","x","y","crs") : ())...,
            (if N == 3
                 (string(ds.varname),)
             else
                 (string(ds.varname,i) for i = 1:ds.nraster)
             end)...
                 ]
end

Base.keys(ds::Dataset) = varnames(ds)

function cf_variable_attrib!(ds::TIFFDataset{T},band,attrib) where T
    fillvalue = ArchGDAL.getnodatavalue(band)
    if !isnothing(fillvalue)
        attrib["_FillValue"] = fillvalue
    end

    add_offset = ArchGDAL.getoffset(band)
    if add_offset != 0
        attrib["add_offset"] = add_offset
    end

    scale_factor = ArchGDAL.getscale(band)
    if scale_factor != 1
        attrib["scale_factor"] = scale_factor
    end

    long_name = get(ds.metadata,"long_name",nothing)
    if !isnothing(long_name)
       attrib["long_name"] = long_name
    end

    valid_range = get(ds.metadata,"valid_range",nothing)
    if !isnothing(valid_range)
        vr = valid_range
        try
            vr = parse.(T,split(vr,','))
        catch
            @warn "cannot parse valid range ", valid_range
        end
        attrib["valid_range"] = vr
    end

end

function variable(ds::Dataset{T,N},varname::SymbolOrString) where {T,N}
    vn = Symbol(varname)
    attrib = OrderedDict{String,Any}()
    geo_ref = geo_referenced(ds)

    if (vn == ds.varname) && (N == 3)
        if geo_ref
            attrib["grid_mapping"] = "crs"
        end
        band = ArchGDAL.getband(ds.dataset,1)
        cf_variable_attrib!(ds,band,attrib)
        return Variable{T,N}(ds,0,attrib)
    elseif (startswith(string(varname),string(ds.varname))) && (N == 2)
        if geo_ref
            attrib["grid_mapping"] = "crs"
        end
        index = parse(Int,replace(string(varname),string(ds.varname)=>""))
        band = ArchGDAL.getband(ds.dataset,index)
        cf_variable_attrib!(ds,band,attrib)
        return Variable{T,N}(ds,index,attrib)
    elseif (vn == :crs) && geo_ref
        cf_crs_attributes!(ds.crs_wkt,attrib; geotransform = ds.geotransform)
        return CRS(ds,attrib)
    elseif vn in (:lon,:lat) && geo_ref
        if vn == :lon
            attrib["standard_name"] = "longitude"
            attrib["units"] = "degrees_east"
        else
            attrib["standard_name"] = "latitude"
            attrib["units"] = "degrees_north"
        end
        return Coord{Float64,2,T,N}(ds,vn,attrib)
    elseif vn in (:x,:y) && geo_ref
        if vn == :x
            attrib["standard_name"] = "projection_x_coordinate"
            #attrib["units"] = "m" # always the case? -> no
        elseif vn == :y
            attrib["standard_name"] = "projection_y_coordinate"
            #attrib["units"] = "m" # always the case? -> no
        end

        Nc = (aligned_grid(ds) ? 1 : 2)
        return Coord{Float64,Nc,T,N}(ds,vn,attrib)
    else
        throw(KeyError("$varname not found"))
    end
end


Base.getindex(ds::Dataset,varname::Union{AbstractString, Symbol}) = cfvariable(ds,varname)

Base.size(v::Variable{T,N}) where {T,N} = (v.parent.width,v.parent.height,v.parent.nraster)[1:N]


function DiskArrays.readblock!(v::Variable{T, 3},
    aout,
    indexes::Vararg{OrdinalRange, 3}) where {T}

    band_index = indexes[3]
    band_out = 1:length(band_index) 

    for i in eachindex(band_index, band_out)
        k = band_index[i]
        k_out = band_out[i]
        band = ArchGDAL.getband(v.parent.dataset, k)

        aout[:,:,k_out] .= band[indexes[1], indexes[2]]
    end

    return nothing
end

function DiskArrays.readblock!(v::Variable{T, 2},
    aout,
    indexes::Vararg{OrdinalRange, 2}) where {T}

    aout .= ArchGDAL.getband(v.parent.dataset, v.index)[indexes...]
    return nothing
end


dimnames(v::Variable) = String.(v.parent.dimnames[1:ndims(v)])

function name(v::Variable{T,N}) where {T,N}
    if N == 3
        return v.parent.varname
    else
        return string(v.parent.varname,v.index)
    end
end


name(v::Coord) = string(v.name)


@inline function xy_geo(ds::Dataset,i,j)
    gt = ds.geotransform
    shift = 0.5 # center
    x_geo = gt[1] + (i-shift) * gt[2] + (j-shift) * gt[3]
    y_geo = gt[4] + (i-shift) * gt[5] + (j-shift) * gt[6]
    return (x_geo,y_geo)
end

function Base.size(v::Coord)
    if ndims(v) == 2
        return (v.parent.width,v.parent.height)
    elseif v.name == :x
        return (v.parent.width,)
    else
        return (v.parent.height,)
    end
end

function dimnames(v::Coord)
    N = ndims(v)
    if N == 2
        return String.(v.parent.dimnames[1:N])
    elseif v.name == :x
        return (String(v.parent.dimnames[1]),)
    else
        return  (String(v.parent.dimnames[2]),)
    end
end

function _get_coord(v::Coord{T,2},i::Integer,j::Integer) where T
    ds = v.parent
    if v.name == :x
        return xy_geo(ds,i,j)[1]
    elseif v.name == :y
        return xy_geo(ds,i,j)[2]
    else
        trans = v.parent.trans
        lat,lon = trans(xy_geo(ds,i,j))
        # can there be a cleverer way?
        if v.name == :lon
            return lon
        else
            return lat
        end
    end
end

function _get_coord(v::Coord{T,1},i::Integer) where T
    ds = v.parent
    if v.name == :x
        return xy_geo(ds,i,1)[1]
    elseif v.name == :y
        return xy_geo(ds,1,i)[2]
    else
        throw(KeyError("unexpected $(v.name)"))
    end
end

function DiskArrays.readblock!(v::Coord{T, 1},
    aout,
    indexes::Vararg{OrdinalRange, 1}) where {T}
    
    index = only(indexes)
    index_out = 1:length(index)

    for t in eachindex(index,index_out)
        i = index[t]
        i_out = index_out[t]
        aout[i_out] = _get_coord(v,i)
    end

    return nothing
end



function DiskArrays.readblock!(v::Coord{T, 2},
    aout,
    indexes::Vararg{OrdinalRange, 2}) where {T}

    indexes_out = ( 1:length(indexes[1]), 1:length(indexes[2]) )

    for t1 in eachindex(indexes[1],indexes_out[1])
        i = indexes[1][t1]
        i_out = indexes_out[1][t1]
        for t2 in eachindex(indexes[2],indexes_out[2])
            j = indexes[2][t2]
            j_out = indexes_out[2][t2]
            aout[i_out,j_out] = _get_coord(v,i,j)
        end
    end

    return nothing
end

Base.size(v::CRS) = ()

function DiskArrays.readblock!(v::CRS,
    aout,
    indexes::Vararg{OrdinalRange, 0})

    aout[] = 0
    return nothing
end

dimnames(v::CRS) = ()
name(v::CRS) = "crs"


"""
    TIFFDataset(fnames::AbstractVector{<:AbstractString};
                aggdim=name, isnewdim=false)

Concatenate TIFF files `fnames` along the dimension specified with `aggdim`.
The parameter `isnewdim` must be true for new dimensions not initially present
in the file.
See CommonDataModel.MFDataset for more information of the arguments.

"""
TIFFDataset(fnames::AbstractVector{<:AbstractString}, args...; kwargs...) =
   MFDataset(TIFFDataset,fnames, args...; kwargs...)


function TIFFDataset(f::Function, args...; kwargs...)
    ds = TIFFDataset(args...; kwargs...)
    try
        f(ds)
    finally
        @debug "closing TIFFDataset" args
        close(ds)
    end
end

DiskArrays.eachchunk(v::Variable{T,3}) where {T} = 
    DiskArrays.GridChunks(v, (size(v)[1:2]...,1))
DiskArrays.haschunks(v::Variable{T,3}) where {T} =  DiskArrays.Chunked()
DiskArrays.haschunks(v::Variable) =  DiskArrays.Unchunked()
DiskArrays.haschunks(v::Coord) = DiskArrays.Unchunked()
DiskArrays.haschunks(v::CRS) = DiskArrays.Unchunked()

export TIFFDataset

end
