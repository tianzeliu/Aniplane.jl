## Test the free surface condition of Aniplane.jl by comparing the waveforms produced using models with a free surface and the ones with an "air layer" on the top

using Base
using Plots

include("../src/ModelOperations.jl")
include("../src/WavePropagation.jl")

## Input parameters

path_air = "./layered_AirOverSimpleCrust.mod"
path_sur = "./layered_SimpleCrust.mod"

rayp = 0.06
baz = 0.
btime = -10.
etime = 20.
dt = 0.01

## Read the models

hlist_air, rholist_air, voilist_air, _ = readanilyrmod(path_air)
hlist_sur, rholist_sur, voilist_sur, _ = readanilyrmod(path_sur)

npts = round(Int, (etime-btime)/dt)

## Compute the waveforms using the model with a free surface and plot

dispmat = getaniresp(3, 0, rayp, baz, hlist_sur, rholist_sur, voilist_sur, npts, dt, btime)
rdisp_sur = dispmat[1, :]
zdisp_sur = dispmat[3, :]

timeax = range(btime, etime, length=npts)

plot(timeax, [trace_r trace_z], show=true)