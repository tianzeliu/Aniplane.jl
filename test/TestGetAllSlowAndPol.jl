using Base
include("../src/ModelOperations.jl")
include("../src/WavePropagation.jl")

## Inputs
path_in = "./layered_HexagonalSandwich_hori_ENU.dat";

rayp = 0.04
baz = 0.0

## Read the velocity model
list_h, list_rho, list_voi, flag_coord = readanilyrmod(path_in)
# dispmatrix(list_voi[2])

## Get the slownesses and polarizations
list_vslow, list_dsst = getallslowandpol(rayp, baz, list_rho, list_voi)
dispmatrix(list_vslow[1])
dispmatrix(list_vslow[2])
dispmatrix(list_vslow[3])