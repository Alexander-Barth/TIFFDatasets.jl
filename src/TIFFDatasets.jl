module TIFFDatasets

import ArchGDAL
import Base: size, keys, getindex
import CommonDataModel: AbstractDataset, AbstractVariable,
    path, name,
    attribnames, attrib,
    dimnames, dim,
    varnames, variable, cfvariable

using DataStructures: OrderedDict
import GDAL
import Proj
import JSON3

struct TIFFDataset{T,N,Ttrans} <: AbstractDataset
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
end

const Dataset = TIFFDataset

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
attrib(ds::Union{Dataset,Variable,CRS,Coord},name::AbstractString) = ds.attrib[name]

dimnames(ds::Dataset) = keys(ds.dim)
dim(ds::Union{TIFFDataset,Variable,CRS,Coord},name::AbstractString) = ds.dim[name]

path(ds::Dataset) = ds.fname

aligned_grid(ds::Dataset) = (ds.geotransform[3] == 0) && (ds.geotransform[5] == 0)

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


"""
    ds = TIFFDataset(fname::AbstractString; varname = "band",
                     projection = "EPSG:4326",
                     catbands = false,
                     dimnames = ("cols","rows","bands"))

Create a data set structure from the GeoTIFF file name `fname`
using [ArchGDAL](https://github.com/yeesian/ArchGDAL.jl). The data variable will be called `varname` (followed by the band number).
The dataset `ds` will also have the virtual variables `x`, `y` and `lon`, `lat`
representing the projected coordinates (as defined in the GeoTIFF file) and the
corresponding longitude/latitude (using the projection as defined
by the `projection` parameter). The names of the dimensions
can be adapted using the parameter `dimnames`.

Setting `catbands` to true, will concatenate all bands into a single
variable called `varname`. If `catbands` is true, it is the user's
responsibility to check that all bands have the same no-data value,
offset and scale factor.
Check with the command line tool `gdalinfo filename` or use ArchGDAL to inspect the file.
If this is not the case, `catbands = true` should not be used.
"""
function TIFFDataset(fname::AbstractString; varname = "band",
                     projection = "EPSG:4326",
                     catbands = false,
                     dimnames = ("cols","rows","bands"),
                     attrib = OrderedDict{String,Any}(
                         "Conventions" => "CF-1.8",
                     ),
                     )
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

    return TIFFDataset{T,N,typeof(trans)}(
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

function cf_variable_attrib!(band,attrib)
    fillvalue = ArchGDAL.getnodatavalue(band)
    if !isnothing(fillvalue)
        attrib["FillValue"] = fillvalue
    end

    add_offset = ArchGDAL.getoffset(band)
    if add_offset != 0
        attrib["add_offset"] = add_offset
    end

    scale_factor = ArchGDAL.getscale(band)
    if scale_factor != 1
        attrib["scale_factor"] = scale_factor
    end
end

function variable(ds::Dataset{T,N},varname::AbstractString) where {T,N}
    vn = Symbol(varname)
    attrib = OrderedDict{String,Any}()
    geo_ref = geo_referenced(ds)

    if (vn == ds.varname) && (N == 3)
        if geo_ref
            attrib["grid_mapping"] = "crs"
        end
        band = ArchGDAL.getband(ds.dataset,1)
        cf_variable_attrib!(band,attrib)
        return Variable{T,N}(ds,0,attrib)
    elseif (startswith(string(varname),string(ds.varname))) && (N == 2)
        if geo_ref
            attrib["grid_mapping"] = "crs"
        end
        index = parse(Int,replace(string(varname),string(ds.varname)=>""))
        band = ArchGDAL.getband(ds.dataset,index)
        cf_variable_attrib!(band,attrib)
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


@inline function Base.getindex(v::Variable{T,2},ij...) where T
    return ArchGDAL.getband(v.parent.dataset, v.index)[ij...]
end

@inline function Base.getindex(v::Variable{T,3},i,j,k::Integer) where T
    return ArchGDAL.getband(v.parent.dataset, k)[i,j]
end

@inline function Base.getindex(v::Variable{T,3},i::Union{AbstractRange,Colon},j::Union{AbstractRange,Colon},k::AbstractRange) where T
    return cat((v[i,j,k_] for k_ = k)...,dims=3)
end

@inline function Base.getindex(v::Variable{T,3},i,j,k::Colon) where T
    return v[i,j,begin:end]
end

function Base.getindex(v::Variable{T,3},indices...) where T
    # fall-back read-all, can be more efficient
    return v[:,:,:][indices...]
end

dimnames(v::Variable) = v.parent.dimnames[1:ndims(v)]

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
        return v.parent.dimnames[1:N]
    elseif v.name == :x
        return (v.parent.dimnames[1],)
    else
        return  (v.parent.dimnames[2],)
    end
end

function Base.getindex(v::Coord{T,2},i::Integer,j::Integer) where T
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

function Base.getindex(v::Coord{T,1},i::Integer) where T
    ds = v.parent
    if v.name == :x
        return xy_geo(ds,i,1)[1]
    elseif v.name == :y
        return xy_geo(ds,1,i)[2]
    else
        throw(KeyError("unexpected $(v.name)"))
    end
end

Base.size(v::CRS) = ()
Base.getindex(v::CRS,indices...) = 0

dimnames(v::CRS) = ()
name(v::CRS) = "crs"

export TIFFDataset

end
