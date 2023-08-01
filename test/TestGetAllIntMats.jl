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

freq = 0.1
omega = freq*2pi

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

## Get the slowness and displacement-stress matrices of Layer 
# tsr = list_tsr[3]
# rho = list_rho[3]
# mat_slow_out, _, mat_dsst, _ = geteigenval(tsr, rayp, rho, [0; 0; 0])

# dispmatrix(mat_dsst[1:3, :])

# Get the slownesses and polarizations
list_vslow, list_dsst = getallslowandpol(rayp, baz, list_rho, list_tsr)

# print(list_vslow[1])
# println()
# print(list_vslow[2])
# println()
# dispmatrix(list_dsst[2][1:3, :])
# println()
# print(list_vslow[3])
# println()

# ispmatrix(list_dsst[2])

# # Get the interface matrices of the first interfacer
# mat_dsst_top = list_dsst[2]
# mat_dsst_bot = list_dsst[3]

# # dispmatrix(mat_dsst_top[1:3, :])
# # println()
# # dispmatrix(mat_dsst_bot[1:3, :])

# mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up = getreftransmats(mat_dsst_top, mat_dsst_bot)

# dispmatrix(mat_ref_down)
# println()

# Get all the interface matrices
list_ref_down, list_tran_down, list_ref_up, list_tran_up = getallintmats(list_vslow, list_dsst)
# dispmatrix(list_ref_down[1])
# println()
# dispmatrix(list_tran_down[1])
# println()
# dispmatrix(list_ref_up[1])
# println()
# dispmatrix(list_tran_up[1])
# println()

# Get the free-surface matrix
mat_dsst_surf = list_dsst[1]
mat_ref_surf = getfreesurfmat(mat_dsst_surf)

# Multiply by the layer matrices
vect_vslow = list_vslow[2]
h = list_h[2]
mat_lyr_down, mat_lyr_up = getlayermat(vect_vslow, h, omega)

mat_ref_down_stk = list_ref_down[2]
mat_tran_down_stk = list_tran_down[2]
mat_ref_up_stk = list_ref_up[2]
mat_tran_up_stk = list_tran_up[2]
multiplylyrmats!(mat_lyr_down, mat_lyr_up, mat_ref_down_stk, mat_tran_down_stk, mat_tran_up_stk)

# Go up the stack by one layer
mat_ref_down = list_ref_down[1]
mat_tran_down = list_tran_down[1]
mat_ref_up = list_ref_up[1] 
mat_tran_up = list_tran_up[1]
propup1lyr!(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)

# print(mat_tran_up_stk[1, 1])
# println()
# print(mat_tran_up_stk[1, 2])
# println()
# print(mat_tran_up_stk[1, 3])
# println()
# print(mat_tran_up_stk[2, 1])
# println()
# print(mat_tran_up_stk[2, 2])
# println()
# print(mat_tran_up_stk[2, 3])
# println()
# print(mat_tran_up_stk[3, 1])
# println()
# print(mat_tran_up_stk[3, 2])
# println()
# print(mat_tran_up_stk[3, 3])

# Multiply layer matrices again
vect_vslow = list_vslow[1]
h = list_h[1]
mat_lyr_down, mat_lyr_up = getlayermat(vect_vslow, h, omega)
multiplylyrmats!(mat_lyr_down, mat_lyr_up, mat_ref_down_stk, mat_tran_down_stk, mat_tran_up_stk)

# print(mat_tran_up_stk[1, 1])
# println()
# print(mat_tran_up_stk[1, 2])
# println()
# print(mat_tran_up_stk[1, 3])
# println()
# print(mat_tran_up_stk[2, 1])
# println()
# print(mat_tran_up_stk[2, 2])
# println()
# print(mat_tran_up_stk[2, 3])
# println()
# print(mat_tran_up_stk[3, 1])
# println()
# print(mat_tran_up_stk[3, 2])
# println()
# print(mat_tran_up_stk[3, 3])

# Apply the free surface condition
mat_tran_up_stk = 1.0I(3)/(1.0I(3)-mat_ref_down_stk*mat_ref_surf)*mat_tran_up_stk

#dispmatrix(mat_ref_surf)
# print(mat_tran_up_stk[1, 1])
# println()
# print(mat_tran_up_stk[1, 2])
# println()
# print(mat_tran_up_stk[1, 3])
# println()
# print(mat_tran_up_stk[2, 1])
# println()
# print(mat_tran_up_stk[2, 2])
# println()
# print(mat_tran_up_stk[2, 3])
# println()
# print(mat_tran_up_stk[3, 1])
# println()
# print(mat_tran_up_stk[3, 2])
# println()
# print(mat_tran_up_stk[3, 3])

# Compute the free-surface displacement
vect_up = mat_tran_up_stk[:, 3]

# Compute the wave vectors of the down-going waves
vect_down = mat_ref_surf*vect_up

# Compute the displacement vector at the free surface
mat_disp = mat_dsst_surf[1:3, :];
vect_disp = mat_disp*[vect_down; vect_up];

dispmatrix(mat_dsst_surf[1:3, :])
println()
print(typeof([vect_down; vect_up]))
println()
print(vect_disp)