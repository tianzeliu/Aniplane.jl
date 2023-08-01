using Base
using Base
using Plots

include("../src/ModelOperations.jl")
include("../src/WavePropagation.jl")

## Inputs
path_in = "./layered_HexagonalSandwich_hori_ENU.mod"

rayp = 0.04
baz = 45.

tdur = 80
t_b = -20
dt = 0.1

## Read the velocity model
list_h, list_rho, list_voi, flag_coord = readanilyrmod(path_in)

## Rotate from E-N-U to N-E-D
numlyr = length(list_h)
list_tsr = [zeros(3, 3, 3, 3) for ind in 1:numlyr]
for ind in 1:numlyr
    # dispmatrix(list_voi[ind])
    # println()
    tsr = voigt2tensor(list_voi[ind])
    if flag_coord == 1
        tsr = enu2ned_tensor(1, tsr)
    end
    list_tsr[ind] = tsr
end

## Compute the response function
npts = Int(tdur/dt)+1
mat_disp_out = getaniresp(3, 0, rayp, baz, list_h, list_rho, list_tsr, npts, dt, t_b)

## Plot the results
trace_t = mat_disp_out[2, :]
trace_r = mat_disp_out[1, :]

t_e = t_b+tdur
grid_time = Array(range(t_b, t_e, length=npts))

print("Plotting...\n")

plotly()
plot(grid_time, [trace_t trace_r], show=true)
# print(minimum(trace_t))
# print(maximum(trace_t))
readline()
