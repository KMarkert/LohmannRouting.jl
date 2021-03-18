using WhereTheWaterFlows


function upstream_catchments(dem::AbstractArray, pourpt::CartesianIndex{2}; snap_pourpt::Bool=true,search_radius=10)
    
    area, slen, dir, nout, nin, pits, c, bnds = waterflows(dem)

    xidx,yidx = Tuple(pourpt)
        
    if snap_pourpt
        window = area[yidx-search_radius:yidx+search_radius,xidx-search_radius:xidx+search_radius]
        val, outletindex = findmax(window)
        outletindex = (outletindex - CartesianIndex(search_radius+1,search_radius+1)) + CartesianIndex(yidx,xidx)
    else
        outletindex = CartesianIndex(yidx,xidx)
    end

    catchments = catchment(dir,outletindex)


    return catchments

end

function flow_sequence(fdir::Array{Int8,2}, start::CartesianIndex)

    c = falses(size(fdir))
    idxes= CartesianIndex{2}[]

    _flowseq!(c,idxes,fdir,start,1)

    return idxes
end


function _flowseq!(c,idxes, dir, ij, iter)
    push!(idxes,ij)
    c[ij] = true
    # proc upstream points
    for IJ in WhereTheWaterFlows.iterate_D9(ij, c)
        ij==IJ && continue
        if WhereTheWaterFlows.flowsinto(IJ, dir[IJ], ij)
            _flowseq!(c, idxes, dir, IJ, iter+1)
        end
    end
end