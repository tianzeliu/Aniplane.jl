### Functions for handling 1D anisotropic velocity models
### Translated from the corresponding MATLAB functions
### Tianze Liu
### 2023-02-12

using Base
using LinearAlgebra
using Printf

## Read a Voigt matrix from a file
function readvoigtmatrix(path_in)
    # Open the file
    lines = open(path_in) do file
        readlines(file)
    end

    # Read the coordinate-system flag and the matrix
    mat_voi = zeros(6, 6)
    lineind = 0
    flag_coord = 0
    for line in lines
        if !startswith(line, "#") && !isequal(line, "")
            lineind = lineind+1
            if lineind == 1
                words = split(line)
                flag_coord = parse(Int, words[1])
                # print(flag_coord)
            else
                words = split(line)
                for colind in 1:6
                    word = words[colind]
                    mat_voi[lineind-1, colind] = parse(Float64, word)
                end
            end
        end
    end

    return mat_voi, flag_coord
end

## Read a layered anisotropic model
##
## Input and output layer thicknesses are in km
## Input and output density are in g/cm^3
## Input and output elastic moduli are in GPa

function readanilyrmod(path_in)

    io = open(path_in)
    # Read the headers
    skipchars(isspace, io, linecomment='#')
    line = readline(io)
    words = split(line)
    flag_coord = parse(Int, words[1])
    skipchars(isspace, io, linecomment='#')
    line = readline(io)
    words = split(line)
    numlyr = parse(Int, words[1])

    # Read the parameters of each layer
    list_h = zeros(numlyr, 1)
    list_rho = zeros(numlyr, 1)
    list_voi = [zeros(6, 6) for ind in 1:numlyr]
    for ind_lyr in 1:numlyr
        skipchars(isspace, io, linecomment='#')
        line = readline(io)
        words = split(line)
        h = parse(Float64, words[1])
        list_h[ind_lyr] = h

        skipchars(isspace, io, linecomment='#')
        line = readline(io)
        words = split(line)
        rho = parse(Float64, words[1])
        list_rho[ind_lyr] = rho

        skipchars(isspace, io, linecomment='#')
        mat_voi = zeros(6, 6)
        for ind_row in 1:6
            line = readline(io)
            words = split(line)
            for ind_col in 1:6
                mat_voi[ind_row, ind_col] = parse(Float64, words[ind_col])
            end
        end
        list_voi[ind_lyr] = mat_voi
    end
    
    close(io)

    return list_h, list_rho, list_voi, flag_coord
end

## Write a layered anisotropic model
##
## Input and output layer thicknesses are in km
## Input and output density are in g/cm^3
## Input and output elastic moduli are in GPa
##
## flag_ela = 1: Elastic moduli are in Voigt format, 2: tensor format

function writeanilyrmod(path_out, flag_coord, flag_ela, list_h, list_rho, list_ela, header)
    numlyr = length(list_h)

    # Open the output file
    io = open(path_out, "w")

    # Write the header
    println(io, "# ", header)
    println(io)

    # Write the coordinate system
    println(io, string(flag_coord), " # 1 = E-N-U, 2 = N-E-D")
    println(io)

    # Write the number of layers
    println(io, string(numlyr), " # Number of layers")
    println(io)

    # Write each layer
    for ind_lyr in 1:numlyr
        h = list_h[ind_lyr]
        rho = list_rho[ind_lyr]
        ela = list_ela[ind_lyr]

        if flag_ela == 2
            mat_voi = tensor2voigt(ela)
        else
            mat_voi = ela
        end

        println(io, "# Layer ", string(ind_lyr))
        println(io)
        println(io, string(h))
        println(io)
        println(io, string(rho))
        println(io)
        
        for ind = 1:6
            row = mat_voi[ind, :]
            f = Ref(Printf.Format("%3.3f"))
            println(io, join(Printf.format.(f, row), ' '))
        end

        println(io)
    end
    
    # Close the file
    close(io)
end

## Display a realmatrix

function dispmatrix(mat_voi)
    for ind1 in 1:size(mat_voi, 1)
        for ind2 in 1:size(mat_voi, 2)
            @printf("%.3e\t", mat_voi[ind1, ind2])
        end
        println()
    end
end

## Convert the indices of a Voigt matrix to the corresponding indcies of the elastic tensor
function voigtind2tensorind(ind_voi)
    if ind_voi == 1
        ind1_tsr = 1
        ind2_tsr = 1
    elseif ind_voi == 2
        ind1_tsr = 2
        ind2_tsr = 2
    elseif ind_voi == 3
        ind1_tsr = 3
        ind2_tsr = 3
    elseif ind_voi == 4
        ind1_tsr = 2
        ind2_tsr = 3
    elseif ind_voi == 5
        ind1_tsr = 1
        ind2_tsr = 3
    elseif ind_voi == 6
        ind1_tsr = 1
        ind2_tsr = 2     
    end

    return ind1_tsr, ind2_tsr
end

## Convert a Voigt matrix to the corresponding 4th order elastic tensor
function voigt2tensor(mat_voi)
    # Check if the Voigt matrix is symmetric
    if mat_voi != mat_voi'
        print("The Voigt matrix is not symmetric! Exit.")
        tsr = []
        return tsr
    end

    tsr = zeros(3, 3, 3, 3)
    num = 0
    for ind1_voi in 1:6
        ind1_tsr, ind2_tsr = voigtind2tensorind(ind1_voi)
        
        for ind2_voi in ind1_voi:6
            ind3_tsr, ind4_tsr = voigtind2tensorind(ind2_voi)
            c = mat_voi[ind1_voi, ind2_voi]

            if ind1_tsr == ind2_tsr && ind3_tsr == ind4_tsr && ind2_tsr == ind3_tsr 
                tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr] = c
                num = num+1
            elseif ind1_tsr == ind2_tsr && ind3_tsr == ind4_tsr && ind2_tsr != ind3_tsr
                tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind1_tsr, ind2_tsr] = c
                num = num+2
            elseif ind1_tsr == ind2_tsr && ind3_tsr != ind4_tsr
                tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind2_tsr, ind1_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind1_tsr, ind2_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind2_tsr, ind1_tsr] = c
                num = num+4
            elseif ind1_tsr != ind2_tsr && ind3_tsr == ind4_tsr
                tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind2_tsr, ind1_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind1_tsr, ind2_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind2_tsr, ind1_tsr] = c
                num = num+4
            elseif ind1_tsr != ind2_tsr && ind3_tsr != ind4_tsr && ind1_tsr == ind3_tsr && ind2_tsr == ind4_tsr 
                tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind1_tsr, ind2_tsr, ind4_tsr, ind3_tsr] = c
                tsr[ind2_tsr, ind1_tsr, ind4_tsr, ind3_tsr] = c
                tsr[ind2_tsr, ind1_tsr, ind3_tsr, ind4_tsr] = c
                num = num+4
            else
                tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind2_tsr, ind1_tsr, ind3_tsr, ind4_tsr] = c
                tsr[ind2_tsr, ind1_tsr, ind4_tsr, ind3_tsr] = c
                tsr[ind4_tsr, ind3_tsr, ind2_tsr, ind1_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind2_tsr, ind1_tsr] = c
                tsr[ind3_tsr, ind4_tsr, ind1_tsr, ind2_tsr] = c
                tsr[ind4_tsr, ind3_tsr, ind1_tsr, ind2_tsr] = c
                tsr[ind1_tsr, ind2_tsr, ind4_tsr, ind3_tsr] = c
                num = num+8
            end
        end
    end

    return tsr
end

## Convert an elastic tensor to a matrix in Voigt format

function tensor2voigt(tsr)

    if size(tsr) != (3, 3, 3, 3)
        mat_voi = []

        print("The elastic tensor must be (3, 3, 3, 3)! Exit.\n")
        return mat_voi
    end

    mat_voi = zeros(6, 6)
    for ind1_voi in 1:6
        ind1_tsr, ind2_tsr = voigtind2tensorind(ind1_voi)

        for ind2_voi in 1:6
            ind3_tsr, ind4_tsr = voigtind2tensorind(ind2_voi)
            mat_voi[ind1_voi, ind2_voi] = tsr[ind1_tsr, ind2_tsr, ind3_tsr, ind4_tsr]
        end
    end    

    return mat_voi
end

## Apply a rotation matrix to a tensor
function applytensorrotation(tsr_in, mat_rot)
    tsr_out = zeros(3, 3, 3, 3)
    for ind1 in 1:3
        for ind2 in 1:3
            for ind3 in 1:3
                for ind4 in 1:3
                    for ind5 in 1:3
                        for ind6 in 1:3
                            for ind7 in 1:3
                                for ind8 in 1:3
                                    tsr_out[ind1, ind2, ind3, ind4] = tsr_out[ind1, ind2, ind3, ind4]+tsr_in[ind5, ind6, ind7, ind8]*mat_rot[ind1, ind5]*mat_rot[ind2, ind6]*mat_rot[ind3, ind7]*mat_rot[ind4, ind8]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return tsr_out
end

## Rotate a elastic tensor around an axis
## The axis orientation and the rotation directions are defined following the right-hand rule! (see CoordinateRotation.pdf for details.)

function rotatetensor(axid, theta, tsr_in)

    # Define the rotation matrix
    if axid == 1
        mat_rot = [[1 0 0]; [0 cosd(theta) -sind(theta)]; [0 sind(theta) cosd(theta)]]
    elseif axid == 2
        mat_rot = [[cosd(theta) 0 sind(theta)]; [0 1 0]; [-sind(theta) 0 cosd(theta)]]
    elseif axid == 3
        mat_rot = [[cosd(theta) -sind(theta) 0]; [sind(theta) cosd(theta) 0]; [0 0 1]]
    else
        print("The axis ID has to be 1, 2, or 3! Exit.")
        tsr_out = []
        return tsr_out
    end

    # # Define the rotation tensor
    # tsr_rot = zeros(3, 3, 3, 3, 3, 3, 3, 3)
    # for ind1 in 1:3
    #     for ind2 in 1:3
    #         for ind3 in 1:3
    #             for ind4 in 1:3
    #                 for ind5 in 1:3
    #                     for ind6 in 1:3
    #                         for ind7 in 1:3
    #                             for ind8 in 1:3
    #                                 tsr_rot[ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8] = mat_rot[ind1, ind5]*mat_rot[ind2, ind6]*mat_rot[ind3, ind7]*mat_rot[ind4, ind8]
    #                             end
    #                         end
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end


    # Perform the rotation through tensor multiplication
    tsr_out = zeros(3, 3, 3, 3)

    tsr_out = applytensorrotation(tsr_in, mat_rot)

    return tsr_out
end

## Rotate the coordinate sysntem of a tensor by a certain degree
## The axis orientation and the rotation directions are defined following the right-hand rule! (see CoordinateRotation.pdf for details.)
function rotatecoord_matrix(axid, theta, mat_in)
    if axid == 1
        mat_rot = [[1 0 0]; [0 cosd(theta) sind(theta)]; [0 -sind(theta) cosd(theta)]]
    elseif axid == 2
        mat_rot = [[cosd(theta) 0 -sind(theta)]; [0 1 0]; [sind(theta) 0 cosd(theta)]]
    elseif axid == 3
        mat_rot = [[cosd(theta) sind(theta) 0]; [-sind(theta) cosd(theta) 0]; [0 0 1]]
    else
        print("The axis ID has to be 1, 2, or 3! Exit.")
        mat_out = []
        return mat_out
    end

    mat_out = mat_rot*mat_in

    return mat_out
end


## Rotate the coordinate sysntem of a tensor by a certain degree
## The axis orientation and the rotation directions are defined following the right-hand rule! (see CoordinateRotation.pdf for details.)

function rotatecoord_tensor(axid, theta, tsr_in)

    # Define the rotation matrix
    if axid == 1
        mat_rot = [[1 0 0]; [0 cosd(theta) sind(theta)]; [0 -sind(theta) cosd(theta)]]
    elseif axid == 2
        mat_rot = [[cosd(theta) 0 -sind(theta)]; [0 1 0]; [sind(theta) 0 cosd(theta)]]
    elseif axid == 3
        mat_rot = [[cosd(theta) sind(theta) 0]; [-sind(theta) cosd(theta) 0]; [0 0 1]]
    else
        print("The axis ID has to be 1, 2, or 3! Exit.")
        tsr_out = []
        return tsr_out
    end

    # # Define the rotation tensor
    # tsr_rot = zeros(3, 3, 3, 3, 3, 3, 3, 3)
    # for ind1 = 1:3
    #     for ind2 = 1:3
    #         for ind3 = 1:3
    #             for ind4 = 1:3
    #                 for ind5 = 1:3
    #                     for ind6 = 1:3
    #                         for ind7 = 1:3
    #                             for ind8 = 1:3
    #                                 tsr_rot(ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8) = mat_rot(ind1, ind5)*mat_rot(ind2, ind6)*mat_rot(ind3, ind7)*mat_rot(ind4, ind8);
    #                             end
    #                         end
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end

    # Perform the rotation
    tsr_out = zeros(3, 3, 3, 3)
    tsr_out = applytensorrotation(tsr_in, mat_rot)

    return tsr_out
end

## Coordinate transformation for a tensor between the E-N-U and N-E-D coordinate systems
##
## flag_dir = 1: E-N-U to N-E-D, 2: N-E-D to E-N-U

function enu2ned_tensor(flag_dir, tsr_in)

    if flag_dir == 1
        tsr_out = rotatecoord_tensor(1, 180, tsr_in)
        tsr_out = rotatecoord_tensor(3, -90, tsr_out)
    elseif flag_dir == 2
        tsr_out = rotatecoord_tensor(3, 90, tsr_in)
        tsr_out = rotatecoord_tensor(1, -180, tsr_out)
    else
        tsr_out = []
        print("The direction flag must be 1 or 2! Exit.")
        return tsr_out
    end

    return tsr_out
end

## Coordinate transformation for a matrix between the E-N-U and N-E-D coordinate systems
##
## flag_dir = 1: E-N-U to N-E-D, 2: N-E-D to E-N-U

function enu2ned_matrix(flag_dir, mat_in)

    if flag_dir == 1
        mat_out = rotatecoord_matrix(1, 180, mat_in)
        mat_out = rotatecoord_matrix(3, -90, mat_out)
    elseif flag_dir == 2
        mat_out = rotatecoord_matrix(3, 90, mat_in)
        mat_out = rotatecoord_matrix(1, -180, mat_out)
    else
        mat_out = []
        print("The direction flag must be 1 or 2! Exit.")
        return mat_out
    end

    return mat_out
end

## Build a Voigt matrix with a hexagonal symmetry with the isotropic velocities, anisotropy percentage, and axis orientations as inputs
# flag_ax = 1: Fast axis, 2: Slow axis

function aniperc2voigt(flag_ax, flag_coord, vp_iso, vs_iso, perc_vp, perc_vs, phi, theta, ita, rho)
    # Convert the units to SI
    vp_iso = vp_iso*1000
    vs_iso = vs_iso*1000
    rho = rho*1000

    # Compute the voigt matrix with a vertical symmetry axis in N-E-D system
    if flag_ax == 1
        c = (vp_iso)^2*(1+perc_vp/100/2)^2*rho
        a = (vp_iso)^2*(1-perc_vp/100/2)^2*rho
        l = (vs_iso)^2*(1+perc_vs/100/2)^2*rho
        n = (vs_iso)^2*(1-perc_vs/100/2)^2*rho
        h = a-2*n
        f = ita*sqrt(a-l)*sqrt(c-l)-l
    elseif flag_ax == 2
        c = (vp_iso)^2*(1-perc_vp/100/2)^2*rho
        a = (vp_iso)^2*(1+perc_vp/100/2)^2*rho
        l = (vs_iso)^2*(1-perc_vs/100/2)^2*rho
        n = (vs_iso)^2*(1+perc_vs/100/2)^2*rho
        h = a-2*n
        f = ita*sqrt(a-l)*sqrt(c-l)-l
    end
    mat_voi_init = zeros(6, 6)
    mat_voi_init[1, :] = [a h f 0 0 0]
    mat_voi_init[2, :] = [h a f 0 0 0]
    mat_voi_init[3, :] = [f f c 0 0 0]
    mat_voi_init[4, :] = [0 0 0 l 0 0]
    mat_voi_init[5, :] = [0 0 0 0 l 0]
    mat_voi_init[6, :] = [0 0 0 0 0 n]

    tsr_init = voigt2tensor(mat_voi_init)

    # Rotate the symmetry axis to the desired direction
    tsr_out = rotatetensor(2, 90-theta, tsr_init)
    tsr_out = rotatetensor(3, phi, tsr_out)

    # If the output coordinate system is E-N-U, rotate the coordinate system 
    if flag_coord == 1
        tsr_out = enu2ned_tensor(2, tsr_out)
    end
    mat_voi_out = tensor2voigt(tsr_out)
    mat_voi_out = mat_voi_out/1e9

    return mat_voi_out
end