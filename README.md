# LohmannRouting.jl
Implementation of the Lohmann streamflow routing model in Julia

The purpose of this model is to take runoff and baseflow simulated by a gridded hydrology model and estimate streamflow at a given point along the stream network. 

This is a work in progress but is planned to be developed into a full-fledged Julia package. Contributions are welcome! So if you have comments/suggestions or want to contribute code, please submit an issue or pull request.

## Installation
LohmannRouting.jl is not currently a Julia registered package.

Currently only developer version is available. It can be installed using:

```
$ julia
julia> ]
(@v1.5) pkg> dev https://github.com/kmarkert/LohmannRouting.jl
```

### References

* Lohmann et al. (1996) https://doi.org/10.1034/j.1600-0870.1996.t01-3-00009.x
* Lohmann et al. (1998) https://doi.org/10.1080/02626669809492107
