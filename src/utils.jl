using Dates,StatsBase,DataFrames

include("constants.jl")

function match_dates(df::DataFrame,dates::Vector{Date};datecol::Symbol=:datetime)
    # TODO: check to make sure that dates range is within df[:,"datetime"] range
    # produces inconsistent results if not ... not sure why...
    mask = sum(map(x -> x.==df[:,datecol],dates),dims=1)[1]
    indices = findall(!iszero,mask)

    df[indices,:]

end

function upsample(img, dshape; method=mean)
    nshape = size(img)
    ydim,xdim = dshape

    window = [Int(round(v / dshape[i])) for (i,v) in enumerate(nshape)]
    wradius = Int.(round.(window / 2))

    outimg = Array{Float64,2}(undef,dshape)

    for i in 1:ydim, j in 1:xdim
        imin = (i-1) * window[1] 
        imin = imin < 1 ? 1 : imin
        imax = imin + window[1] 
        imax = imax > nshape[1] ? nshape[1] : imax
        jmin = (j-1) * window[2]
        jmin = jmin < 1 ? 1 : jmin
        jmax = jmin + window[2]
        jmax = jmax > nshape[2] ? nshape[2] : jmax

        outimg[i,j] = method(vec(img[imin:imax,jmin:jmax]))
    end

    return outimg
end

function dd_to_meters(lon::Float64, lat::Float64; scale::Float64=0.1)
    """
    Function to convert decimal degrees to meters based on the approximation
    given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system
    Args:
    inPt (list or array): A Y,X point provided in geographic coordinates
    in that order.
    Keywords:
    scale (int): Resolution of the raster value to covert into meters,
    must be in decimal degrees.
    Returns:
    list: List of Y,X resolution values converted from meters to decimal degrees
    """

    radlat = deg2rad(lat) # convert degree latitude to radians

    ba::Float64 = 0.99664719 # constant of b/a

    ss = atan(ba*tan(radlat)) # calculate the reduced latitude

    # factor to convert meters to decimal degrees for X axis
    xfct = (pi / 180.0) * EARTHRADIUS * cos(ss)

    # factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82 * cos(2*radlat)+1.175*cos(4*radlat)-
            0.0023*cos(6*radlat))

    # get meter resolution
    y_meters = scale * yfct
    x_meters = scale * xfct

    # return list of converted resolution values
    return hcat([x_meters], [y_meters])
end

function findnearest(lons,lats,x,y)
    xval,xidx = findmin(abs.(lons .- x))
    yval,yidx = findmin(abs.(lats .- y))

    return xidx,yidx
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function calc_area(lon,lat)
    dd_res = mean([abs(lon[1]-lon[2]),abs(lat[1]-lat[2])])
    lons,lats = meshgrid(lon,lat)
    sides = dd_to_meters.(lons,lats;scale=dd_res)

    return prod.(sides)
end