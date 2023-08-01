using Base
include("../src/ModelOperations.jl")
include("../src/WavePropagation.jl")

## Inputs
path_in = "./layered_HexagonalSandwich_hori_ENU.dat";

rayp = 0.04
baz = 45.

tdur = 80
t_b = -20
dt = 0.01

## Read the velocity model
list_h, list_rho, list_voi, flag_coord = readanilyrmod(path_in)

## Rotate from E-N-U to N-E-D
numlyr = length(list_h)
list_tsr = [zeros(3, 3, 3, 3) for ind in 1:numlyr]
for ind in 1:numlyr
    tsr = voigt2tensor(list_voi[ind])
    if flag_coord == 1
        tsr = enu2ned_tensor(1, tsr)
    end
    list_tsr[ind] = tsr
end

## Get the slownesses and polarizations
list_vslow, list_dsst = getallslowandpol(rayp, baz, list_rho, list_tsr)

# Get the free surface matrix
mat_dsst_surf = list_dsst[1]
mat_ref_surf = getfreesurfmat(mat_dsst_surf)
print(mat_ref_surf)