# LohmannRouting.jl

```julia
using Dates, CSV, DataFrames, NetCDF, StatsBase, Plots, LohmannRouting

fluxes = NetCDF.open("../ut-vic/data/bcc-csm1-1-m_fluxes.1980-01-01.nc");

time = NetCDF.readvar(fluxes,"time");
lon = NetCDF.readvar(fluxes,"lon");
lat = NetCDF.readvar(fluxes,"lat");
runoff = NetCDF.readvar(fluxes,"OUT_RUNOFF");
baseflow = NetCDF.readvar(fluxes,"OUT_BASEFLOW");
total_runoff = runoff + baseflow;

dshape = size(total_runoff)[1:2]

ds = NetCDF.open("../ut-vic/data/jordan_valley_terrain.nc");

elv = NetCDF.readvar(ds,"elevation");

elv_grid = upsample(elv,dshape,method=mean)

# format [x.xxxxx, y.yyyy]
pourpt = [-111.5485,39.986];
pourptidx = CartesianIndex(findnearest(lon,lat,pourpt[1],pourpt[2]))

mon_uh = CSV.read("../ut-vic/data/usgs_gauge/USGS_10167499_UH.csv",DataFrame);

xmask = calc_xmask(lon,lat)
area = calc_area(lon,lat)

flowidxs = flow_sequence(dir,pourptidx)

irfs = make_irf.(xmask[flowidxs]; diffusion=800.,velocity=1.5)

river_irfs = make_river_irf(irfs; max_basin_days=5)

unit_hydrograph = make_grid_uh(river_irfs,mon_uh[:,:unit_hydrograph];max_basin_days=5)

runoff_in = total_runoff[flowidxs]
area_in = area[flowidxs]
           
q = rout(runoff_in,unit_hydrograph,area_in)

q_time = @. Date(1,1,1) + Dates.Day.(time-1) - Dates.Day(1)

sim_df = DataFrame(time = q_time, discharge = q)

obs_df = CSV.read("../ut-vic/data/usgs_gauge/USGS_10167499_streamflow.csv",DataFrame)
obs_df[!,:cms] = obs_df[:,:discharge] * 0.0283168;

sim_matched = match_dates(sim_df,obs_df[:,:datetime];datecol=:time)

plot(obs_df[:,:datetime],obs_df[:,:cms],label="Observed",ylabel="Discharge [cms]",lw=2,color="#2e2e2e",dpi=200)
plot!(df[:,:datetime],sim_matched[:,:discharge],label="Simulated",lw=2,color="#1f77b4")
# savefig()
```