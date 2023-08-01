### Test RotateElasticTensor.m
## The function is tested by comparing its output with the tensor computed element-wise
## The element-wise rotation is found to have the wrong sign. Fixed.
## Tianze Liu, 2023-02-13

# push!(LOAD_PATH, "../src")
using Base
using LinearAlgebra
include("../src/ModelOperations.jl")

## Inputs
path_in = "./hexagonal_fast_horizontal_ENU.dat"

axid = 2
theta = 30

## Prepare the elastic tensor
mat_voi_in, flag_coord, flat_stat = readvoigtmatrix(path_in)

tsr_in, flag = voigt2tensor(mat_voi_in)

if flag == -1
    print("The format of the input tensor is wrong! Exit.")
    return
end

## Comupute the rotated tensor using rotatetensor
print("Performing the rotation using rotatetensor..\n")
tsr1_out, flag = rotatetensor(axid, theta, tsr_in)

if flag == -1
    return
end

## Compute the rotated tensor element-wise
# Define the rotation matrix
print("Performing the rotation element-wise...\n")

if axid == 1
    mat_rot = [[1 0 0]; [0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)]]
elseif axid == 2
    mat_rot = [[cosd(theta) 0 sind(theta)]; [0 1 0]; [-sind(theta) 0 cosd(theta)]]
elseif axid == 3
    mat_rot = [[cosd(theta) -sind(theta) 0]; [sind(theta) cosd(theta) 0]; [0 0 1]]
end

dispmatrix(mat_rot)

# Compute each element
tsr2_out = zeros(3, 3, 3, 3)
for ind1 in 1:3
    for ind2 in 1:3
        for ind3 in 1:3
            for ind4 in 1:3
                for ind5 in 1:3
                    for ind6 in 1:3
                        for ind7 in 1:3
                            for ind8 in 1:3
                                tsr2_out[ind1, ind2, ind3, ind4] = tsr2_out[ind1, ind2, ind3, ind4]+tsr_in[ind5, ind6, ind7, ind8]*mat_rot[ind1, ind5]*mat_rot[ind2, ind6]*mat_rot[ind3, ind7]*mat_rot[ind4, ind8]
                            end
                        end
                    end
                end
            end
        end
    end
end

# print(typeof(tsr1_out))
# print(typeof(tsr2_out))

## Compare the Voigt matrices before and after rotation
mat1_voi_out, _ = tensor2voigt(tsr1_out)
mat2_voi_out, _ = tensor2voigt(tsr2_out)

dispmatrix(mat_voi_in)
println()
dispmatrix(mat1_voi_out)
println()
dispmatrix(mat2_voi_out)