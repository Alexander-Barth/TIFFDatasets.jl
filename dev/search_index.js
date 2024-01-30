var documenterSearchIndex = {"docs":
[{"location":"#TIFFDatasets","page":"Home","title":"TIFFDatasets","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [TIFFDatasets]","category":"page"},{"location":"#TIFFDatasets.TIFFDataset-Tuple{AbstractString}","page":"Home","title":"TIFFDatasets.TIFFDataset","text":"ds = TIFFDataset(fname::AbstractString; varname = \"band\",\n                 projection = \"EPSG:4326\",\n                 catbands = false,\n                 dimnames = (\"cols\",\"rows\",\"bands\"))\n\nCreate a data set structure from the GeoTIFF file name fname using ArchGDAL. The data variable will be called varname (followed by the band number). The dataset ds will also have the virtual variables x, y and lon, lat representing the projected coordinates (as defined in the GeoTIFF file) and the corresponding longitude/latitude (using the projection as defined by the projection parameter). The names of the dimensions can be adapted using the parameter dimnames. The type of projection can be checked with the attributes (grid_mapping_name, longitude_of_central_meridian, false_easting, false_northing, latitude_of_projection_origin, scale_factor_at_central_meridian, longitude_of_prime_meridian, semi_major_axis, inverse_flattening, ... depending on the projection used) of the virtual variable crs. The attribute crs_wkt describes the coordinate reference system as so-called \"Well-Known Text\"\n\nIn addition to the CF conventions, there is also the attribute GeoTransform [1] with six numbers describing a linear mapping between indices and the projected coordinates. For the point (i,j) the projected coordinates are computed as:\n\nshift = 0.5 # center\nx_geo = GeoTransform[1] + (i-shift) * GeoTransform[2] + (j-shift) * GeoTransform[3]\ny_geo = GeoTransform[4] + (i-shift) * GeoTransform[5] + (j-shift) * GeoTransform[6]\n\nSetting catbands to true, will concatenate all bands into a single variable called varname. If catbands is true, it is the user's responsibility to check that all bands have the same no-data value, offset and scale factor. Check with the command line tool gdalinfo filename or use ArchGDAL to inspect the file. If this is not the case, catbands = true should not be used.\n\nExample:\n\nusing TIFFDatasets\nfname = download(\"https://data-assimilation.net/upload/Alex/TIFF/S2_1-12-19_48MYU_0.tif\")\nds = TIFFDataset(fname)\nband1 = ds[\"band1\"][:,:];\nlongitude = ds[\"lon\"][:,:];\nlatitude = ds[\"lat\"][:,:];\n# x,y native projection\nx = ds[\"x\"][:,:];\ny = ds[\"y\"][:,:];\n# coordinate reference system as \"Well-Known Text\"\ncrs_wkt = ds[\"crs\"].attrib[\"crs_wkt\"];\n\nFor performance reasons, one should avoid loading a large number of individual elements (or using function that do so, like sum, mean ...) and rather load the whole array. For example:\n\njulia> @time sum(ds[\"band1\"]) # slow\n  0.062227 seconds (721.65 k allocations: 22.045 MiB)\n6403.1743f0\n\njulia> @time sum(ds[\"band1\"][:,:]) # much faster\n  0.000381 seconds (765 allocations: 302.781 KiB)\n6403.5527f0\n\nThe dataset can be manipulated using the functions from CommonDataModel.jl. Writing to a file is currently not supported (pull-requests are welcome!).\n\nReference:\n\n[1] GDAL GeoTransform\n\n\n\n\n\n","category":"method"}]
}