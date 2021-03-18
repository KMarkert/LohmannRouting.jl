
using DSP

include("constants.jl")
include("utils.jl")
include("watersheds.jl")


function calc_xmask(lon,lat)
    dd_res = mean([abs(lon[1]-lon[2]),abs(lat[1]-lat[2])])
    lons,lats = meshgrid(lon,lat)
    sides = dd_to_meters.(lons,lats;scale=dd_res)

    xdim,ydim = size(sides)

    xmask = zeros(Float64, xdim,ydim)

    for i in 1:xdim, j in 1:ydim
        xmask[i,j] = sqrt( +((sides[i,j].^2)...)) / 2
    end

    return xmask
end

function make_irf(xmask::Float64; diffusion::Float64=800.0,velocity::Float64=1.5, max_runoff_days::Int64=2)
    """
    Function for Calculate the impulse response function for a grid cells using equation 15 from Lohmann, et al. (1996) Tellus article
    Arguments:
        xmask -- longest path for any location within the grid cell to reach
                 the stream channel [m]
        diff -- diffusion value
        velo -- overland flow velocity value
    Returns:
        irf -- normalized impulse response function for a grid cell
    Keywords:
        diff -- the diffusivity parameter (default 800) [m^2/s]
        velo -- the flow velocity parameter (default 1.5) [m/s]
    """
    # get number of hours per impulse (typically 48 hours)
    nt = max_runoff_days * HOURSPERDAY
    step = collect(1:1:nt) #number of time steps

    # convert time steps to seconds
    time = @. step * SECSPERHOUR

        
    exponent = @. -1 * (velocity * time - xmask)^2 / (4.0 * diffusion * time) 
    #impulse response function
    green = @. xmask / (2.0 * time * sqrt(pi * time * diffusion)) * exp(exponent)

    irf = green / sum(green) #normalize the output IRF

    return irf

end

function make_river_irf(irfs; max_runoff_days::Int64=2, max_basin_days::Int64=50)

    t_cell = Int(max_runoff_days * HOURSPERDAY)
    t_uh = Int(max_basin_days * HOURSPERDAY)

    n_cells = length(irfs)

    river_irf = zeros(Float64,t_uh,n_cells)

    river_irf[1:t_cell, 1] = irfs[1]

    for i in 2:n_cells
        irf_temp = conv(river_irf[:, i-1], irfs[i])

        river_irf[:,i] = irf_temp[1:t_uh] / sum(irf_temp[1:t_uh])
    end

    return river_irf

end

function make_grid_uh(river_irfs,uh_box; max_runoff_days::Int64=2, max_basin_days::Int64=50)

    t_cell = Int(max_runoff_days * HOURSPERDAY)
    t_uh = Int(max_basin_days * HOURSPERDAY)

    uh_shape = size(river_irfs)

    n_cells = uh_shape[end]

    # prealloc
    uh = zeros(Float64,uh_shape...)
    irf_temp = zeros(Float64, (t_uh + t_cell))

    # set init values
    irf_temp[1:length(uh_box)] = uh_box[:]
    uh[:,1] = irf_temp[1:t_uh]

    for i in 2:n_cells
        irf_temp = conv(uh_box, river_irfs[:,i-1])
        uh[:,i] = irf_temp[1:t_uh] / sum(irf_temp[1:t_uh])
    end

    return uh
end

function rout(runoff,unit_hydrograph,area)

    n_cells = size(unit_hydrograph)[end]
    days = size(runoff)[end]
    buffer = size(unit_hydrograph)[1]

    basinuh = zeros(Float64,n_cells,days) # ...basin wide UH
    griduh = zeros(Float64,days+buffer,buffer)
    println(size(griduh))

    for i in 1:n_cells

        factor = (area[i] / (SECSPERDAY * MMPERMETER))

        for u in 1:buffer
            # calculate the percent contributions from pixel to outlet at each time period
            griduh[u:days+u-1,u] = (runoff[i,:] * factor) * unit_hydrograph[u]
        end

        basinuh[i,:] = sum(griduh[1:days,:],dims=2)

    end

    finalq = dropdims(sum(basinuh,dims=1),dims=1)

    return finalq

end