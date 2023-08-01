using Base
include("../src/ModelOperations.jl")
include("../src/WavePropagation.jl")

## Inputs
path_in = "./hexagonal_fast_horizontal_ENU.dat";

rayp = 0.04
rho = 3.3

vect_sym = zeros(3, 1);

## Read the elastic tensor
mat_voi, _ = readvoigtmatrix(path_in)
# print(mat_voi)

tsr, flag = voigt2tensor(mat_voi)

if flag == -1
    return
end

## Get the eigen values and vectors
mat_slow_out, mat_pol_out, mat_dsst, mat_dsst_inv = geteigenval(tsr, rayp, rho, vect_sym)

## Display the results
dispmatrix(mat_slow_out)
println()
dispmatrix(mat_pol_out)