### Functions for computing wave propagation in 1D anisotropic velocity models
### Translated from the corresponding MATLAB functions
### Tianze Liu
### 2023-02-13

using Base
using LinearAlgebra
using Printf
using FFTW

## Sort the three wave modes by velocity in descending order

function sortwaves(vect_eval, mat_evec_in, hslow, vect_sym)
    # Threshold for determining if the two S slownesses are equal
    eps = 1e-5

    # Extract the polarization vectors
    mat_pol_in = mat_evec_in[4:6, :]

    # Determine if a symmety axis is given

    if norm(vect_sym) == 0
        flag_sym = -1
    else
        flag_sym = 1
    end

    #dispmatrix(mat_evec_in)

    # Normalize the eigen vectors and the polarization vectors by the norms of the polarization vectors
    mat_pol_out = mat_pol_in
    mat_evec_out = mat_evec_in
    for ind = 1:6
        vect = mat_pol_in[:, ind]
        r = sqrt((vect'*conj(vect)));
        mat_pol_out[:, ind] = mat_pol_in[:, ind]/r
        mat_evec_out[:, ind] = mat_evec_in[:, ind]/r
    end

    # Sort the waves by decreasing real vertical slowness to put the down-going waves before the up-going ones
    list_q = vect_eval
    list_qreal = zeros(6)
    list_qreal[:] = real(list_q)

    list_ind = sortperm(list_qreal, rev=true)
    sort!(list_qreal, rev=true)

    list_q = list_q[list_ind]
    mat_pol_out = mat_pol_out[:, list_ind] # For now, I don't sort the results by imaginary vertical slowness as Peter did, 2022/10/02
    mat_evec_out = mat_evec_out[:, list_ind]

    # Sort the up- and down-going waves into qP, qSP, and qSR
    list_qreal_down = list_qreal[1:3]
    list_q_down = list_q[1:3]
    mat_pol_down_out = mat_pol_out[:, 1:3]
    mat_evec_down_out = mat_evec_out[:, 1:3]
    list_ind = sortperm(abs.(list_qreal_down))
    list_q_down = list_q_down[list_ind]
    # print(list_q_down)
    mat_pol_out[:, 1:3] = mat_pol_down_out[:, list_ind]
    mat_evec_out[:, 1:3] = mat_evec_down_out[:, list_ind]

    list_qreal_up = list_qreal[4:6]
    list_q_up = list_q[4:6]
    mat_pol_up_out = mat_pol_out[:, 4:6]
    mat_evec_up_out = mat_evec_out[:, 4:6]
    list_ind = sortperm(abs.(list_qreal_up))
    list_q_up = list_q_up[list_ind]

    mat_pol_out[:, 4:6] = mat_pol_up_out[:, list_ind]
    mat_evec_out[:, 4:6] = mat_evec_up_out[:, list_ind]

    list_q = [list_q_down; list_q_up]

    # print(list_q)
    # println()
    # dispmatrix(mat_evec_out)

    # print(list_q)
    # dispmatrix(mat_evec_out)

    # In the case of qS1 = qS2, set SH along 2 direction and SV along the direction defined by the cross product between the slowness vector and 2 direction. Then put SV before SH
    qs1 = list_q[2]
    qs2 = list_q[3]
    if abs((qs1-qs2)/qs1) < eps
        qs = qs1
        vect_pol_sh = [0; 1; 0]
        vect_slow = [1; 0; real(qs)]
        vect_pol_sv = cross(vect_slow, vect_pol_sh)
        vect_pol_sv = vect_pol_sv/norm(vect_pol_sv)
         
        vect_evec_sv = [vect_pol_sv*qs; vect_pol_sv]
        vect_evec_sh = [vect_pol_sh*qs; vect_pol_sh]

        mat_pol_out[:, 2] = vect_pol_sv
        mat_pol_out[:, 3] = vect_pol_sh

        mat_evec_out[:, 2] = vect_evec_sv
        mat_evec_out[:, 3] = vect_evec_sh       
#         vect_pol_s1 = mat_pol_out(:, 2);
#         vect_pol_s2 = mat_pol_out(:, 3);
# 
#         vect_evec_s1 = mat_evec_out(:, 2);
#         vect_evec_s2 = mat_evec_out(:, 3);
# 
#         if abs(vect_pol_s1(3)) > abs(vect_pol_s2(3))
#             vect_pol_sv = vect_pol_s1;
#             vect_pol_sh = vect_pol_s2;
# 
#             vect_evec_sv = vect_evec_s1;
#             vect_evec_sh = vect_evec_s2;
#         else
#             vect_pol_sv = vect_pol_s2;
#             vect_pol_sh = vect_pol_s1;
# 
#             vect_evec_sv = vect_evec_s2;
#             vect_evec_sh = vect_evec_s1;            
#         end
# 
#         mat_pol_out(:, 2) = vect_pol_sv;
#         mat_pol_out(:, 3) = vect_pol_sh;
# 
#         mat_evec_out(:, 2) = vect_evec_sv;
#         mat_evec_out(:, 3) = vect_evec_sh;
    end

    qs1 = list_q[5]
    qs2 = list_q[6]
    if abs((qs1-qs2)/qs1) < eps
        qs = qs1;
        vect_pol_sh = [0; 1; 0];
        vect_slow = [1; 0; real(qs)];
        vect_pol_sv = cross(vect_slow, vect_pol_sh)
        vect_pol_sv = vect_pol_sv/norm(vect_pol_sv)
         
        vect_evec_sv = [vect_pol_sv*qs; vect_pol_sv]
        vect_evec_sh = [vect_pol_sh*qs; vect_pol_sh]

        mat_pol_out[:, 5] = vect_pol_sv
        mat_pol_out[:, 6] = vect_pol_sh

        mat_evec_out[:, 5] = vect_evec_sv
        mat_evec_out[:, 6] = vect_evec_sh 
        
#         vect_pol_s1 = mat_pol_out(:, 5);
#         vect_pol_s2 = mat_pol_out(:, 6);
# 
#         vect_evec_s1 = mat_evec_out(:, 5);
#         vect_evec_s2 = mat_evec_out(:, 6);
# 
#         if abs(vect_pol_s1(3)) > abs(vect_pol_s2(3))
#             vect_pol_sv = vect_pol_s1;
#             vect_pol_sh = vect_pol_s2;
# 
#             vect_evec_sv = vect_evec_s1;
#             vect_evec_sh = vect_evec_s2;
#         else
#             vect_pol_sv = vect_pol_s2;
#             vect_pol_sh = vect_pol_s1;
# 
#             vect_evec_sv = vect_evec_s2;
#             vect_evec_sh = vect_evec_s1;            
#         end
# 
#         mat_pol_out(:, 5) = vect_pol_sv;
#         mat_pol_out(:, 6) = vect_pol_sh;
# 
#         mat_evec_out(:, 5) = vect_evec_sv;
#         mat_evec_out(:, 6) = vect_evec_sh;
    end

    # Assemble the slowness vectors
    mat_slow_out = zeros(Complex, 3, 6)
    mat_slow_out[1, 1:6] .= hslow
    for ind = 1:6
        q = list_q[ind]
        mat_slow_out[3, ind] = hslow*q
    end

    # If the symmetry axis is specified, sort the S waves by increasing
    # projection on the qSR direction
    if flag_sym == 1
        vect_pol1 = mat_pol_out[:, 2]
        vect_slow1 = mat_slow_out[:, 2]
            
        vect_pol2 = mat_pol_out[:, 3]
        vect_slow2 = mat_slow_out[:, 3]

        vect_pol3 = mat_pol_out[:, 5]
        vect_slow3 = mat_slow_out[:, 5]

        vect_pol4 = mat_pol_out[:, 6]
        vect_slow4 = mat_slow_out[:, 6]

        vect_qsr1 = cross(vect_slow1, vect_sym)
        vect_qsr2 = cross(vect_slow2, vect_sym)
        vect_qsr3 = cross(vect_slow3, vect_sym)
        vect_qsr4 = cross(vect_slow4, vect_sym)

        if norm(vect_qsr1)*norm(vect_qsr2)*norm(vect_qsr3)*norm(vect_qsr4) > 0
            vect_qsr1 = vect_qsr1/norm(vect_qsr1)
            vect_qsr2 = vect_qsr2/norm(vect_qsr2)
            vect_qsr3 = vect_qsr3/norm(vect_qsr3)
            vect_qsr4 = vect_qsr4/norm(vect_qsr4)

            dotprod1 = abs(dot(vect_qsr1, vect_pol1))
            dotprod2 = abs(dot(vect_qsr2, vect_pol2))
            dotprod3 = abs(dot(vect_qsr3, vect_pol3))
            dotprod4 = abs(dot(vect_qsr4, vect_pol4))
            
            if dotprod1 > dotprod2
                vect_slow_tmp = mat_slow_out[:, 3]
                mat_slow_out[:, 3] = mat_slow_out[:, 2]
                mat_slow_out[:, 2] = vect_slow_tmp

                vect_pol_tmp = mat_pol_out[:, 3]
                mat_pol_out[:, 3] = mat_pol_out[:, 2]
                mat_pol_out[:, 2] = vect_pol_tmp

                evec_tmp = mat_evec_out[:, 3]
                mat_evec_out[:, 3] = mat_evec_out[:, 2]
                mat_evec_out[:, 2] = evec_tmp
            end

            if dotprod3 > dotprod4
                vect_slow_tmp = mat_slow_out[:, 6]
                mat_slow_out[:, 6] = mat_slow_out[:, 5]
                mat_slow_out[:, 5] = vect_slow_tmp

                vect_pol_tmp = mat_pol_out[:, 6]
                mat_pol_out[:, 6] = mat_pol_out[:, 5]
                mat_pol_out[:, 5] = vect_pol_tmp

                evec_tmp = mat_evec_out[:, 6]
                mat_evec_out[:, 6] = mat_evec_out[:, 5]
                mat_evec_out[:, 5] = evec_tmp
            end
        else
            print("At least one slowness vector is parallel the symmetry axis. Using the default order for the qS waves.")
        end
    end

    return mat_slow_out, mat_pol_out, mat_evec_out
end

## Compute the eigen values and eigen vectors of for a general anisotropic solid and a given horizontal slowness
## tsr: Elastic tensor in GPa
## hslow: Horizontal slowness in s/km
## rho: Density in g/cm3
## vect_sym: Hexagonal symmetry axis
##
## Based on Peter Shearer's anieigen.f, which in turn is based on equations
## in Keith and Crampin (GJRAS 49, 181-208, 1977), Crampin (Wave Motion 3, 343-391,1981), and
## Garmany (GJRAS 75, 565-569, 1983).
##
## The elastic tensor is rotated so that Axis 1 is aligned with the
## direction of the horizontal slowness!

function geteigenval(tsr, hslow, rho, vect_sym)
    # Convert to SI
    tsr = tsr*1e9
    hslow = hslow/1000
    rho = rho*1000

    # Assemble the small matrices in Crampin 1981
    mat_r = zeros(3, 3)
    mat_v = zeros(3, 3)
    mat_t = zeros(3, 3)
    hvel = 1/hslow

    for ind1 in 1:3
        for ind2 in 1:3
            mat_r[ind1, ind2] = tsr[ind1, 3, ind2, 3]
            mat_v[ind1, ind2] = tsr[ind1, 3, ind2, 1]
            mat_t[ind1, ind2] = tsr[ind1, 1, ind2, 1]
        end
    end

    mat_t_hat = mat_t-rho*hvel^2*1.0I(3)
    mat_s = mat_v+mat_v'

    if det(mat_r) == 0
        print("Matrix R is singlular! Exit.\n")
        mat_slow_out = []
        mat_pol_out = []
        mat_dsst = []
        mat_dsst_inv = []
        return mat_slow_out, mat_pol_out, mat_dsst, mat_dsst_inv
    end

    mat_scrm1 = mat_r\mat_s
    mat_scrm2 = mat_r\mat_t_hat

    # Set up the big matrix for finding eigenvectors
    mat_big = [-mat_scrm1 -mat_scrm2; 1.0I(3) zeros(3, 3)];

    # Get the eigen values and vectors
    fac = eigen(mat_big)
    vect_eval = fac.values
    mat_evec = fac.vectors

    # print(vect_eval)
    # println()
    # dispmatrix(mat_evec)

    # Sort the eigen values and eigen vectors
    mat_slow_out, mat_pol_out, mat_evec_out = sortwaves(vect_eval, mat_evec, hslow, vect_sym)
    # dispmatrix(mat_slow_out*1000)
    # println()
    # dispmatrix(mat_evec_out)

    # Compute the displacement-stress matrix following Eqn. 3.11 of Crampin (1981)
    mat_e = [zeros(3, 3) 1.0I(3); mat_r mat_v]
    mat_dsst_scr = mat_e*mat_evec_out
    mat_dsst = mat_dsst_scr # For debugging purpose, 2022/10/12

#     #Remove the i*hslow coefficient from the stress terms in the Crampin
#     #definition 
#     mat_dsst_jga = mat_dsst_scr;
#     mat_dsst_jga(4:6, :) = mat_dsst_jga(4:6, :)/1i/hslow;
#     mat_dsst = mat_dsst_jga;
# 
#     #Normalize the stress-displacement matrix as in Eqn. 20 of Garmany 1983
#     for ind = 1:6
#         vect_dsst = mat_dsst(:, ind);
#         r = sqrt(2*vect_dsst(1:3)'*vect_dsst(4:6));
#         mat_dsst(:, ind) = mat_dsst(:, ind)/r;
#     end

    # Compute the inverse of the displacement-stress matrix following Eqn.23 of Garmany 1983
    mat_dsst_inv = zeros(Complex, 6, 6)
    mat_dsst_inv[1:3, 1:3] = mat_dsst[4:6, 1:3]'
    mat_dsst_inv[1:3, 4:6] = mat_dsst[1:3, 1:3]'
    mat_dsst_inv[4:6, 1:3] = mat_dsst[4:6, 4:6]'
    mat_dsst_inv[4:6, 4:6] = mat_dsst[1:3, 4:6]'

    # Convert back to s/km
    mat_slow_out = mat_slow_out*1000

    return mat_slow_out, mat_pol_out, mat_dsst, mat_dsst_inv
end

## Get the reflection coefficients for the free surface
##
## The order of the displacement-stress matrices is:
## [P-down, S1-down, S2-down, P-up, S1-up, S2-up]
##
## The order of the reflection matrix is:
## R{P,up,P}  R{S1,up,P}  R{S2,up,P} 
## R{P,up,S1} R{S1,up,S1} R{S2,up,S1}
## R{P,up,S2} R{S1,up,S2} R{S2,up,S2}

function getfreesurfmat(mat_dsst)
    # Extract the stress components of the matrix
    mat_str_down = -mat_dsst[4:6, 1:3]
    mat_str_up = mat_dsst[4:6, 4:6]

    # Compute the reflection matrix
    if abs(det(mat_str_down)) < 0.000001
        print("The scattering matrix is nearly singular! Exit.\n")
        mat_ref= []
        return mat_ref
    end

    mat_ref = mat_str_down\mat_str_up ### Am I getting this right? Eqn. 5.77 in Kennett (2009) is counter-intuitive. 2023-02-22

    return mat_ref
end

## Special case of getfreesurfmat that handles P-SV waves

function getfreesurfmat_psv(mat_dsst)
    # Extract the stress components of the matrix
    mat_str_down = -mat_dsst[3:4, 1:2]
    mat_str_up = mat_dsst[3:4, 3:4]

    # Compute the reflection matrix
    if abs(det(mat_str_down)) < 0.000001
        print("The scattering matrix is nearly singular! Exit.\n")
        mat_ref= []
        return mat_ref
    end

    mat_ref = mat_str_down\mat_str_up

    return mat_ref
end

## Special case of getfreesurfmat that handles SH waves

function getfreesurfmat_sh(mat_dsst)
    # Extract the stress components of the matrix
    str_down = -mat_dsst[2, 1]
    str_up = mat_dsst[2, 2]

    # Compute the reflection matrix
    if abs(det(str_down)) < 0.000001
        print("The scattering matrix is nearly singular! Exit.\n")
        mat_ref= []
        return mat_ref
    end

    ref = str_up/str_down

    return ref
end


## Get the reflection and transmission matrices at an interface
##
## Based on Peter Shearer's ANIFACE in aniplanesub.f with a reversed order for the reflection and transmission efficients in the
## reflection-transmission matrix
## 
## The order of the displacement-stress matrices is:
## [P-down, S1-down, S2-down, P-up, S1-up, S2-up]
##
## The order of the reflection-transmission matrix is:
## R{P,down,P}  R{S1,down,P}  R{S2,down,P}  T{P,up,P}  T{S1,up,P}  T{S2,up,P} 
## R{P,down,S1} R{S1,down,S1} R{S2,down,S1} T{P,up,S1} T{S1,up,S1} T{S2,up,S1}
## R{P,down,S2} R{S1,down,S2} R{S2,down,S2} T{P,up,S2} T{S1,up,S2} T{S2,up,S2}
## T{P,down,P}  T{S1,down,P}  T{S2,down,P}  R{P,up,P}  R{S1,up,P}  R{S2,up,P}
## T{P,down,S1} T{S1,down,S1} T{S2,down,S1} R{P,up,S1} R{S1,up,S1} R{S2,up,S1}
## T{P,down,S2} T{S1,down,S2} T{S2,down,S2} R{P,up,S2} R{S1,up,S2} R{S2,up,S2}
## whose order for reflection and transmission coefficients is opposite from
## ANIFACE!

function getreftransmats(mat_dsst_top, mat_dsst_bot)
    # Split the displacement-stress matrices into the incident and  scattered matrices
    mat_dsst_inc = [mat_dsst_top[:, 1:3] -mat_dsst_bot[:, 4:6]]
    mat_dsst_sca = [-mat_dsst_top[:, 4:6] mat_dsst_bot[:, 1:3]]

    # For debug use, 2023-02-15
    # mat_dsst_sca[3, 1] = -mat_dsst_sca[3, 1]
    # mat_dsst_sca[1, 2] = -mat_dsst_sca[1, 2]
    # mat_dsst_sca[2, 2] = -mat_dsst_sca[2, 2]
    # mat_dsst_sca[1, 3] = -mat_dsst_sca[1, 3]
    # mat_dsst_sca[2, 3] = -mat_dsst_sca[2, 3]

    # mat_dsst_sca[3, 4] = -mat_dsst_sca[3, 4]

    #dispmatrix(mat_dsst_sca[1:3, :])
    #dispmatrix(inv(mat_dsst_sca))

    # Compute the reflection-transmission matrix
    if abs(det(mat_dsst_sca)) < 0.000001
        print("The scattering matrix is nearly singular! Exit.")
        mat_ref_down = []
        mat_tran_down = []
        mat_ref_up = []
        mat_tran_up = []
        return mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up
    end
    mat_rt = mat_dsst_sca\mat_dsst_inc;
    #dispmatrix(mat_rt)

    # Separate the submatrices
    mat_ref_down = mat_rt[1:3, 1:3]
    mat_tran_down = mat_rt[4:6, 1:3]
    mat_ref_up = mat_rt[4:6, 4:6]
    mat_tran_up = mat_rt[1:3, 4:6]

    return mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up
end

## Get the reflection and transmission matrices for the P-SV system at an interface

function getreftransmats_psv(mat_dsst_top, mat_dsst_bot)
    # Split the displacement-stress matrices into the incident and  scattered matrices
    mat_dsst_inc = [mat_dsst_top[:, 1:2] -mat_dsst_bot[:, 3:4]]
    mat_dsst_sca = [-mat_dsst_top[:, 3:4] mat_dsst_bot[:, 1:2]]

    # For debug use, 2023-02-15
    # mat_dsst_sca[3, 1] = -mat_dsst_sca[3, 1]
    # mat_dsst_sca[1, 2] = -mat_dsst_sca[1, 2]
    # mat_dsst_sca[2, 2] = -mat_dsst_sca[2, 2]
    # mat_dsst_sca[1, 3] = -mat_dsst_sca[1, 3]
    # mat_dsst_sca[2, 3] = -mat_dsst_sca[2, 3]

    # mat_dsst_sca[3, 4] = -mat_dsst_sca[3, 4]

    #dispmatrix(mat_dsst_sca[1:3, :])
    #dispmatrix(inv(mat_dsst_sca))

    # Compute the reflection-transmission matrix
    if abs(det(mat_dsst_sca)) < 0.000001
        print("The scattering matrix is nearly singular! Exit.")
        mat_ref_down = []
        mat_tran_down = []
        mat_ref_up = []
        mat_tran_up = []
        return mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up
    end
    mat_rt = mat_dsst_sca\mat_dsst_inc;
    #dispmatrix(mat_rt)

    # Separate the submatrices
    mat_ref_down = mat_rt[1:2, 1:2]
    mat_tran_down = mat_rt[3:4, 1:2]
    mat_ref_up = mat_rt[3:4, 3:4]
    mat_tran_up = mat_rt[1:2, 3:4]

    return mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up
end

## Get the reflection and transmission matrices for the SH system at an interface

function getreftransmats_sh(mat_dsst_top, mat_dsst_bot)
    # Split the displacement-stress matrices into the incident and  scattered matrices
    mat_dsst_inc = [mat_dsst_top[:, 1] -mat_dsst_bot[:, 2]]
    mat_dsst_sca = [-mat_dsst_top[:, 2] mat_dsst_bot[:, 1]]

    # For debug use, 2023-02-15
    # mat_dsst_sca[3, 1] = -mat_dsst_sca[3, 1]
    # mat_dsst_sca[1, 2] = -mat_dsst_sca[1, 2]
    # mat_dsst_sca[2, 2] = -mat_dsst_sca[2, 2]
    # mat_dsst_sca[1, 3] = -mat_dsst_sca[1, 3]
    # mat_dsst_sca[2, 3] = -mat_dsst_sca[2, 3]

    # mat_dsst_sca[3, 4] = -mat_dsst_sca[3, 4]

    #dispmatrix(mat_dsst_sca[1:3, :])
    #dispmatrix(inv(mat_dsst_sca))

    # Compute the reflection-transmission matrix
    if abs(det(mat_dsst_sca)) < 0.000001
        print("The scattering matrix is nearly singular! Exit.")
        mat_ref_down = []
        mat_tran_down = []
        mat_ref_up = []
        mat_tran_up = []
        return mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up
    end
    mat_rt = mat_dsst_sca\mat_dsst_inc;
    #dispmatrix(mat_rt)

    # Separate the submatrices
    ref_down = mat_rt[1, 1]
    tran_down = mat_rt[2, 1]
    ref_up = mat_rt[2, 2]
    tran_up = mat_rt[1, 2]

    return ref_down, tran_down, ref_up, tran_up
end

## Compute the layer matrix

function getlayermat(vect_vslow, h, omega)
    vect_vslow_down = vect_vslow[1:3]
    # print(vect_vslow_down*h)
    vect_lyr_down = exp.(-1.0im*omega*vect_vslow_down*h)
    # print(omega)
    # println()
    # print(vect_lyr_down)
    # println()
    mat_lyr_down = diagm(vect_lyr_down)

    vect_vslow_up = vect_vslow[4:6]
    # print(vect_vslow_up*h)
    vect_lyr_up = exp.(1.0im*omega*vect_vslow_up*h)
    mat_lyr_up = diagm(vect_lyr_up)
    
    return mat_lyr_down, mat_lyr_up
end

## Compute the layer matrix for the P-SV system

function getlayermat_psv(vect_vslow, h, omega)
    vect_vslow_down = vect_vslow[1:2]
    # print(vect_vslow_down*h)
    vect_lyr_down = exp.(-1.0im*omega*vect_vslow_down*h)
    # print(omega)
    # println()
    # print(vect_lyr_down)
    # println()
    mat_lyr_down = diagm(vect_lyr_down)

    vect_vslow_up = vect_vslow[3:4]
    # print(vect_vslow_up*h)
    vect_lyr_up = exp.(1.0im*omega*vect_vslow_up*h)
    mat_lyr_up = diagm(vect_lyr_up)
    
    return mat_lyr_down, mat_lyr_up
end

## Compute the layer matrix for the SH system

function getlayermat_sh(vect_vslow, h, omega)
    vslow_down = vect_vslow[1]
    # print(vect_vslow_down*h)
    lyr_down = exp.(-1.0im*omega*vslow_down*h)
    # print(omega)
    # println()
    # print(vect_lyr_down)
    # println()

    vslow_up = vect_vslow[2]
    # print(vect_vslow_up*h)
    lyr_up = exp.(1.0im*omega*vslow_up*h)
    
    return lyr_down, lyr_up
end


## Get the slowness and polarization matrices of all layers
## The code runs in parallel using @threads!

function getallslowandpol(hslow, baz, list_rho, list_tsr, vect_sym)
    numlyr = length(list_rho)
    list_vslow = [zeros(6) for ind in 1:numlyr]
    list_dsst = [zeros(6, 6) for ind in 1:numlyr]

    # Compute the matrices of each layer in parallel
    #Threads.@threads for ind_lyr in 1:numlyr
    for ind_lyr in 1:numlyr
        # Extract the parameters for the layer
        rho = list_rho[ind_lyr]
        tsr = list_tsr[ind_lyr]

        # Rotate the tensor to the local coordinate system
        tsr = rotatecoord_tensor(3, baz+180, tsr)

        # Compute the slownesses and polarizations
        mat_slow_out, _, mat_dsst, _ = geteigenval(tsr, hslow, rho, vect_sym)

        # Remove the sign ambiguity of the displacement-stress matrices
        if ind_lyr > 1
            for ind_vec in 1:6
                mat_dsst_old = list_dsst[ind_lyr-1]

                vect_old = mat_dsst_old[:, ind_vec]
                vect_new = mat_dsst[:, ind_vec]
                if real(vect_old'*vect_new) < 0
                    mat_dsst[:, ind_vec] = -mat_dsst[:, ind_vec]
                end
            end
        end

        # print(mat_slow_out[3, :])
        # println()

        list_vslow[ind_lyr] = mat_slow_out[3, :]
        list_dsst[ind_lyr] = mat_dsst
    end
    return list_vslow, list_dsst
end

## Special case for getallslowandpol that hands VTI mediums
## P-SV is in the 1-3 plane, and SH is in 2 direction!

function getallslowandpol_vti(hslow, baz, list_rho, list_tsr)
    numlyr = length(list_rho)

    list_vslow_psv = [zeros(Complex, 4) for ind in 1:numlyr]
    list_vslow_sh = [zeros(Complex, 2) for ind in 1:numlyr]

    list_dsst_psv = [zeros(Complex, 4, 4) for ind in 1:numlyr]
    list_dsst_sh = [zeros(Complex, 2, 2) for ind in 1:numlyr]

    # Compute the matrices of each layer in parallel
    #Threads.@threads for ind_lyr in 1:numlyr
    for ind_lyr in 1:numlyr
        # Extract the parameters for the layer
        rho = list_rho[ind_lyr]
        tsr = list_tsr[ind_lyr]

        # Rotate the tensor to the local coordinate system
        tsr = rotatecoord_tensor(3, baz+180, tsr)

        # Compute the slownesses and polarizations
        # print(ind_lyr)
        # println()
        mat_slow_out, _, mat_dsst, _ = geteigenval(tsr, hslow, rho, [0.; 0.; 1.;])

        # Remove the sign ambiguity of the displacement-stress matrices
        # if ind_lyr > 1
        #     for ind_vec in 1:4
        #         mat_dsst_old = list_dsst_psv[ind_lyr-1]

        #         vect_old = mat_dsst_old[:, ind_vec]
        #         vect_new = mat_dsst[:, ind_vec]
        #         if real(vect_old'*vect_new) < 0
        #             mat_dsst[:, ind_vec] = -mat_dsst[:, ind_vec]
        #         end
        #     end

        #     for ind_vec in 1:2
        #         mat_dsst_old = list_dsst_sh[ind_lyr-1]

        #         vect_old = mat_dsst_old[:, ind_vec]
        #         vect_new = mat_dsst[:, ind_vec]
        #         if real(vect_old'*vect_new) < 0
        #             mat_dsst[:, ind_vec] = -mat_dsst[:, ind_vec]
        #         end
        #     end
        # end

        list_vslow_psv[ind_lyr] = mat_slow_out[3, [1, 2, 4, 5]]
        list_vslow_sh[ind_lyr] = mat_slow_out[3, [3, 6]]

        list_dsst_psv[ind_lyr] = mat_dsst[[1, 3, 4, 6], [1, 2, 4, 5]]
        list_dsst_sh[ind_lyr] = mat_dsst[[2, 5], [3, 6]]
    end

    return list_vslow_psv, list_vslow_sh, list_dsst_psv, list_dsst_sh
end

## Get the SH displacement-stress matrices of all layers 
## We use Eqn. 3.24, 3.25, and 3.33 in Kennet (2009)
## Input vs is in km/s, rho is in g/cm^3

function getalldsst_sh(hslow, list_vs, list_rho)
    numlyr = length(list_rho)

    list_vslow = [zeros(Complex, 2) for ind in 1:numlyr]
    list_dsst = [zeros(Complex, 2, 2) for ind in 1:numlyr]
    hslow = hslow/1000 # Convert to SI
    for ind in 1:numlyr
        vs = list_vs[ind]*1000 # Convert to SI
        rho = list_rho[ind]*1000 # Convert to SI
        vslow = sqrt(1/vs^2-hslow^2)

        mat_dsst = [1/vs 1/vs; im*rho*vs*vslow -im*rho*vs*vslow]
        mat_dsst = mat_dsst/sqrt(2*rho*vslow)

        list_vslow[ind] = [vslow; -vslow]
        list_dsst[ind] = mat_dsst
    end

    list_vslow = list_vslow*1000 # Convert back to s/km

    return list_vslow, list_dsst
end


## Compute all the interface matrices
## The code runs in parallel using @threads!

function getallintmats(list_dsst)
    numlyr = length(list_dsst)

    list_ref_down = [zeros(Complex, 3, 3) for ind in 1:numlyr-1]
    list_tran_down = [zeros(Complex, 3, 3) for ind in 1:numlyr-1]
    list_ref_up = [zeros(Complex, 3, 3) for ind in 1:numlyr-1]
    list_tran_up = [zeros(Complex, 3, 3) for ind in 1:numlyr-1]

    #Threads.@threads for ind_lyr in 2:numlyr
    for ind_lyr in 2:numlyr
        mat_dsst_top = list_dsst[ind_lyr-1]
        mat_dsst_bot = list_dsst[ind_lyr]
        
        mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up = getreftransmats(mat_dsst_top, mat_dsst_bot)

        list_ref_down[ind_lyr-1] = mat_ref_down
        list_tran_down[ind_lyr-1] = mat_tran_down
        list_ref_up[ind_lyr-1] = mat_ref_up
        list_tran_up[ind_lyr-1] = mat_tran_up
    end

    return list_ref_down, list_tran_down, list_ref_up, list_tran_up
end

## Compute all the interface matrices for the P-SV system

function getallintmats_psv(list_dsst)
    numlyr = length(list_dsst)

    list_ref_down = [zeros(Complex, 2, 2) for ind in 1:numlyr-1]
    list_tran_down = [zeros(Complex, 2, 2) for ind in 1:numlyr-1]
    list_ref_up = [zeros(Complex, 2, 2) for ind in 1:numlyr-1]
    list_tran_up = [zeros(Complex, 2, 2) for ind in 1:numlyr-1]

    #Threads.@threads for ind_lyr in 2:numlyr
    for ind_lyr in 2:numlyr
        mat_dsst_top = list_dsst[ind_lyr-1]
        mat_dsst_bot = list_dsst[ind_lyr]
        
        mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up = getreftransmats_psv(mat_dsst_top, mat_dsst_bot)

        list_ref_down[ind_lyr-1] = mat_ref_down
        list_tran_down[ind_lyr-1] = mat_tran_down
        list_ref_up[ind_lyr-1] = mat_ref_up
        list_tran_up[ind_lyr-1] = mat_tran_up
    end

    return list_ref_down, list_tran_down, list_ref_up, list_tran_up
end

## Compute all the interface matrices for the SH system

function getallintmats_sh(list_dsst)

    numlyr = length(list_dsst)

    list_ref_down = zeros(Complex, numlyr-1)
    list_tran_down = zeros(Complex, numlyr-1)
    list_ref_up = zeros(Complex, numlyr-1)
    list_tran_up = zeros(Complex, numlyr-1)

    #Threads.@threads for ind_lyr in 2:numlyr
    for ind_lyr in 2:numlyr
        mat_dsst_top = list_dsst[ind_lyr-1]
        mat_dsst_bot = list_dsst[ind_lyr]
        
        ref_down, tran_down, ref_up, tran_up = getreftransmats_sh(mat_dsst_top, mat_dsst_bot)

        list_ref_down[ind_lyr-1] = ref_down
        list_tran_down[ind_lyr-1] = tran_down
        list_ref_up[ind_lyr-1] = ref_up
        list_tran_up[ind_lyr-1] = tran_up
    end

    return list_ref_down, list_tran_down, list_ref_up, list_tran_up
end


## Multiply the layer matrices

function multiplylyrmats(mat_lyr_down, mat_lyr_up, mat_ref_down, mat_tran_down, mat_tran_up)
    mat_ref_down = mat_lyr_up*mat_ref_down*mat_lyr_down
    mat_tran_down = mat_tran_down*mat_lyr_down
    mat_tran_up = mat_lyr_up*mat_tran_up

    return mat_ref_down, mat_tran_down, mat_tran_up
end

function multiplylyrmats!(mat_lyr_down, mat_lyr_up, mat_ref_down, mat_tran_down, mat_tran_up)
    mat_ref_down[:] = mat_lyr_up*mat_ref_down*mat_lyr_down
    mat_tran_down[:] = mat_tran_down*mat_lyr_down
    mat_tran_up[:] = mat_lyr_up*mat_tran_up
end

## Multiply the layer matrices for the P-SV system

function multiplylyrmats_psv!(mat_lyr_down, mat_lyr_up, mat_ref_down, mat_tran_down, mat_tran_up)
    mat_ref_down[:] = mat_lyr_up*mat_ref_down*mat_lyr_down
    mat_tran_down[:] = mat_tran_down*mat_lyr_down
    mat_tran_up[:] = mat_lyr_up*mat_tran_up
end

## Multiply the layer matrices for the SH system

function multiplylyrmats_sh(lyr_down, lyr_up, ref_down, tran_down, tran_up)
    ref_down = lyr_up*ref_down*lyr_down
    tran_down = tran_down*lyr_down
    tran_up = lyr_up*tran_up

    return ref_down, tran_down, tran_up
end

## Propagate the wavefield up one layer

function propup1lyr(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
    mat_ref_down_old = mat_ref_down_stk
    mat_tran_down_old = mat_tran_down_stk
    mat_ref_up_old = mat_ref_up_stk
    mat_tran_up_old = mat_tran_up_stk

    mat_reverb_down = 1.0I(3)-mat_ref_up*mat_ref_down_old
    mat_reverb_up = 1.0I(3)-mat_ref_down_old*mat_ref_up

    # print(mat_ref_down)
    # println()
    # print(mat_tran_up)
    # println()
    # print(mat_tran_down)
    # println()

    # Kennett (2009) Eqn. 6.3
    #mat_ref_down_stk = mat_ref_down+mat_tran_up*mat_ref_down_old/mat_reverb_down*mat_tran_down
    #mat_tran_down_stk = mat_tran_down_old/mat_reverb_down*mat_tran_down
    mat_ref_down_stk = mat_ref_down+mat_tran_up*mat_ref_down_old*inv(mat_reverb_down)*mat_tran_down
    mat_tran_down_stk = mat_tran_down_old*inv(mat_reverb_down)*mat_tran_down

    # Kennett (2009) Eqn. 6.4
    # mat_ref_up_stk = mat_ref_up_old+mat_tran_down_old*mat_ref_up/mat_reverb_up*mat_tran_up_old
    # mat_tran_up_stk = mat_tran_up/mat_reverb_up*mat_tran_up_old
    mat_ref_up_stk = mat_ref_up_old+mat_tran_down_old*mat_ref_up*inv(mat_reverb_up)*mat_tran_up_old
    mat_tran_up_stk = mat_tran_up*inv(mat_reverb_up)*mat_tran_up_old

    return mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk
end

function propup1lyr!(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
    mat_ref_down_old = mat_ref_down_stk
    mat_tran_down_old = mat_tran_down_stk
    mat_ref_up_old = mat_ref_up_stk
    mat_tran_up_old = mat_tran_up_stk

    mat_reverb_down = 1.0I(3)-mat_ref_up*mat_ref_down_old
    mat_reverb_up = 1.0I(3)-mat_ref_down_old*mat_ref_up

    # print(mat_ref_down)
    # println()
    # print(mat_tran_up)
    # println()
    # print(mat_tran_down)
    # println()

    # Kennett (2009) Eqn. 6.3
    mat_ref_down_stk[:] = mat_ref_down+mat_tran_up*mat_ref_down_old*inv(mat_reverb_down)*mat_tran_down
    mat_tran_down_stk[:] = mat_tran_down_old*inv(mat_reverb_down)*mat_tran_down

    # Kennett (2009) Eqn. 6.4
    mat_ref_up_stk[:] = mat_ref_up_old+mat_tran_down_old*mat_ref_up*inv(mat_reverb_up)*mat_tran_up_old
    mat_tran_up_stk[:] = mat_tran_up*inv(mat_reverb_up)*mat_tran_up_old
end

## Propagate the P-SV wavefield up one layer

function propup1lyr_psv!(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
    mat_ref_down_old = mat_ref_down_stk
    mat_tran_down_old = mat_tran_down_stk
    mat_ref_up_old = mat_ref_up_stk
    mat_tran_up_old = mat_tran_up_stk

    mat_reverb_down = 1.0I(2)-mat_ref_up*mat_ref_down_old
    mat_reverb_up = 1.0I(2)-mat_ref_down_old*mat_ref_up

    # print(mat_ref_down)
    # println()
    # print(mat_tran_up)
    # println()
    # print(mat_tran_down)
    # println()

    # Kennett (2009) Eqn. 6.3
    mat_ref_down_stk[:] = mat_ref_down+mat_tran_up*mat_ref_down_old*inv(mat_reverb_down)*mat_tran_down
    mat_tran_down_stk[:] = mat_tran_down_old*inv(mat_reverb_down)*mat_tran_down

    # Kennett (2009) Eqn. 6.4
    mat_ref_up_stk[:] = mat_ref_up_old+mat_tran_down_old*mat_ref_up*inv(mat_reverb_up)*mat_tran_up_old
    mat_tran_up_stk[:] = mat_tran_up*inv(mat_reverb_up)*mat_tran_up_old
end

## Propagate the SH wavefield up one layer

function propup1lyr_sh(ref_down, tran_down, ref_up, tran_up, ref_down_stk, tran_down_stk, ref_up_stk, tran_up_stk)
    ref_down_old = ref_down_stk
    tran_down_old = tran_down_stk
    ref_up_old = ref_up_stk
    tran_up_old = tran_up_stk

    reverb_down = 1.0-ref_up*ref_down_old
    reverb_up = 1.0-ref_down_old*ref_up

    # print(mat_ref_down)
    # println()
    # print(mat_tran_up)
    # println()
    # print(mat_tran_down)
    # println()

    # Kennett (2009) Eqn. 6.3
    ref_down_stk = ref_down+tran_up*ref_down_old/reverb_down*tran_down
    tran_down_stk = tran_down_old/reverb_down*tran_down

    # Kennett (2009) Eqn. 6.4
    ref_up_stk = ref_up_old+tran_down_old*ref_up/reverb_up*tran_up_old
    tran_up_stk = tran_up/reverb_up*tran_up_old

    return ref_down_stk, tran_down_stk, ref_up_stk, tran_up_stk
end

## Propagate the wavefield up the whole stack
## Loop over all frequencies using @threads!

function propupstack!(list_tran_freq, list_ref_freq, vect_corr, vect_omega, list_h, list_vslow, list_ref_down, list_tran_down, list_ref_up, list_tran_up, flag_inc)
    numfreq = length(vect_omega)
    numlyr = length(list_h)

    # Loop over all frequencies using @threads
    #Threads.@threads for ind_freq = 1:numfreq
    for ind_freq = 1:numfreq

        omega = vect_omega[ind_freq]

        phs_corr = 1.0

        mat_ref_down_stk = zeros(Complex, 3, 3)
        mat_ref_up_stk = zeros(Complex, 3, 3)
        mat_tran_down_stk = zeros(Complex, 3, 3)
        mat_tran_up_stk = zeros(Complex, 3, 3)

        mat_ref_down_stk[:] = list_ref_down[numlyr-1]
        mat_ref_up_stk[:] = list_ref_up[numlyr-1]
        mat_tran_down_stk[:] = list_tran_down[numlyr-1]
        mat_tran_up_stk[:] = list_tran_up[numlyr-1]

        # if ind_freq < 10
        #     # dispmatrix(mat_ref_down_stk)
        #     # println()
        #     # dispmatrix(mat_ref_up_stk)
        #     # println()
        #     # dispmatrix(mat_tran_down_stk)
        #     # println()
        #     print(list_tran_up[numlyr-1])
        #     println()
        # end

        for ind_lyr in numlyr-1:-1:1
            h = list_h[ind_lyr]
            # print(h)
            # println()
            # print(ind_lyr)

            # Compute the layer matrix
            vect_vslow = list_vslow[ind_lyr]
    
            mat_lyr_down, mat_lyr_up = getlayermat(vect_vslow, h, omega)
        
            # Compute the phase correction matrix
            # When the incident wave is P, use the cumulative P travel time
            # When the incident wave is S, use the cumulative S travel time
            if flag_inc == 1
                vslow_corr = vect_vslow[1]
                phs_corr = phs_corr*exp(1.0im*omega*h*vslow_corr)
            else
                vslow_corr = vect_vslow[2]
                phs_corr = phs_corr*exp(1.0im*omega*h*vslow_corr)                
            end

            # Multiply the layer matrices
            #mat_ref_down_stk, mat_tran_down_stk, mat_tran_up_stk = multiplylyrmats(mat_lyr_down, mat_lyr_up, mat_ref_down_stk, mat_tran_down_stk, mat_tran_up_stk)
            # print(mat_lyr_down)
            # println()
            # print(mat_lyr_up)
            # println()

            # if ind_freq == 100
            #     print(maximum(abs.(mat_tran_up_stk)))
            #     println()
            # end

            multiplylyrmats!(mat_lyr_down, mat_lyr_up, mat_ref_down_stk, mat_tran_down_stk, mat_tran_up_stk)
            
            # print(mat_ref_down_stk)
            # println()
            # print(mat_tran_down_stk)
            # println()
            # print(mat_tran_up_stk)

            # Go up the stack by one layer
            if ind_lyr > 1
                mat_ref_down = list_ref_down[ind_lyr-1]
                mat_tran_down = list_tran_down[ind_lyr-1]
                mat_ref_up = list_ref_up[ind_lyr-1]
                mat_tran_up = list_tran_up[ind_lyr-1]

                #mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk = propup1lyr(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
                propup1lyr!(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
            end
        end

        # if ind_freq < 100
        #     dispmatrix(abs.(mat_tran_up_stk))
        #     println()
        # end
        list_tran_freq[ind_freq] = mat_tran_up_stk
        list_ref_freq[ind_freq] = mat_ref_down_stk
        vect_corr[ind_freq] = phs_corr
    end

    # print("AAAA\n")
    # for ind = 1:10
    #     dispmatrix(abs.(list_tran_freq[ind]))
    #     println()
    # end
end

## Propagate the P-SV wavefield up the whole stack

function propupstack_psv!(list_tran_freq, list_ref_freq, vect_corr, vect_omega, list_h, list_vslow, list_ref_down, list_tran_down, list_ref_up, list_tran_up, flag_inc)
    numfreq = length(vect_omega)
    numlyr = length(list_h)

    # Loop over all frequencies using @threads
    #Threads.@threads for ind_freq = 1:numfreq
    for ind_freq = 1:numfreq

        omega = vect_omega[ind_freq]

        phs_corr = 1.0

        mat_ref_down_stk = zeros(Complex, 2, 2)
        mat_ref_up_stk = zeros(Complex, 2, 2)
        mat_tran_down_stk = zeros(Complex, 2, 2)
        mat_tran_up_stk = zeros(Complex, 2, 2)

        mat_ref_down_stk[:] = list_ref_down[numlyr-1]
        mat_ref_up_stk[:] = list_ref_up[numlyr-1]
        mat_tran_down_stk[:] = list_tran_down[numlyr-1]
        mat_tran_up_stk[:] = list_tran_up[numlyr-1]


        for ind_lyr in numlyr-1:-1:1
            h = list_h[ind_lyr]
            # print(h)
            # println()
            # print(ind_lyr)

            # Compute the layer matrix
            vect_vslow = list_vslow[ind_lyr]
    
            mat_lyr_down, mat_lyr_up = getlayermat_psv(vect_vslow, h, omega)
        
            # Compute the phase correction matrix
            # When the incident wave is P, use the cumulative P travel time
            # When the incident wave is S, use the cumulative S travel time
            if flag_inc == 1
                vslow_corr = vect_vslow[1]
                phs_corr = phs_corr*exp(1.0im*omega*h*vslow_corr)
            else
                vslow_corr = vect_vslow[2]
                phs_corr = phs_corr*exp(1.0im*omega*h*vslow_corr)                
            end

            # Multiply the layer matrices
            multiplylyrmats_psv!(mat_lyr_down, mat_lyr_up, mat_ref_down_stk, mat_tran_down_stk, mat_tran_up_stk)
            
            # print(mat_ref_down_stk)
            # println()
            # print(mat_tran_down_stk)
            # println()
            # print(mat_tran_up_stk)

            # Go up the stack by one layer
            if ind_lyr > 1
                mat_ref_down = list_ref_down[ind_lyr-1]
                mat_tran_down = list_tran_down[ind_lyr-1]
                mat_ref_up = list_ref_up[ind_lyr-1]
                mat_tran_up = list_tran_up[ind_lyr-1]

                #mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk = propup1lyr(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
                propup1lyr_psv!(mat_ref_down, mat_tran_down, mat_ref_up, mat_tran_up, mat_ref_down_stk, mat_tran_down_stk, mat_ref_up_stk, mat_tran_up_stk)
            end
        end

        # if ind_freq < 100
        #     dispmatrix(abs.(mat_tran_up_stk))
        #     println()
        # end
        list_tran_freq[ind_freq] = mat_tran_up_stk
        list_ref_freq[ind_freq] = mat_ref_down_stk
        vect_corr[ind_freq] = phs_corr
    end

    # print("AAAA\n")
    # for ind = 1:10
    #     dispmatrix(abs.(list_tran_freq[ind]))
    #     println()
    # end
end

## Propagate the SH wavefield up the whole stack

function propupstack_sh!(list_tran_freq, list_ref_freq, vect_corr, vect_omega, list_h, list_vslow, list_ref_down, list_tran_down, list_ref_up, list_tran_up)
    numfreq = length(vect_omega)
    numlyr = length(list_h)

    # Loop over all frequencies using @threads
    #Threads.@threads for ind_freq = 1:numfreq
    for ind_freq = 1:numfreq

        omega = vect_omega[ind_freq]

        phs_corr = 1.0

        ref_down_stk = list_ref_down[numlyr-1]
        ref_up_stk = list_ref_up[numlyr-1]
        tran_down_stk = list_tran_down[numlyr-1]
        tran_up_stk = list_tran_up[numlyr-1]


        for ind_lyr in numlyr-1:-1:1
            h = list_h[ind_lyr]
            # print(h)
            # println()
            # print(ind_lyr)

            # Compute the layer matrix
            vect_vslow = list_vslow[ind_lyr]
    
            lyr_down, lyr_up = getlayermat_sh(vect_vslow, h, omega)
        
            # Compute the phase correction matrix
            # When the incident wave is P, use the cumulative P travel time
            # When the incident wave is S, use the cumulative S travel time
            vslow_corr = vect_vslow[1]
            phs_corr = phs_corr*exp(1.0im*omega*h*vslow_corr)

            # Multiply the layer matrices
            # print(tran_up_stk)
            # println()
            ref_down_stk, tran_down_stk, tran_up_stk = multiplylyrmats_sh(lyr_down, lyr_up, ref_down_stk, tran_down_stk, tran_up_stk)
            
            # print(mat_ref_down_stk)
            # println()
            # print(mat_tran_down_stk)
            # println()
            # print(mat_tran_up_stk)

            # Go up the stack by one layer
            if ind_lyr > 1
                ref_down = list_ref_down[ind_lyr-1]
                tran_down = list_tran_down[ind_lyr-1]
                ref_up = list_ref_up[ind_lyr-1]
                tran_up = list_tran_up[ind_lyr-1]

                ref_down_stk, tran_down_stk, ref_up_stk, tran_up_stk = propup1lyr_sh(ref_down, tran_down, ref_up, tran_up, ref_down_stk, tran_down_stk, ref_up_stk, tran_up_stk)
            end
        end

        # if ind_freq < 100
        #     dispmatrix(abs.(mat_tran_up_stk))
        #     println()
        # end
        list_tran_freq[ind_freq] = tran_up_stk
        list_ref_freq[ind_freq] = ref_down_stk
        vect_corr[ind_freq] = phs_corr
    end

    # print("AAAA\n")
    # for ind = 1:10
    #     dispmatrix(abs.(list_tran_freq[ind]))
    #     println()
    # end
end

## Apply the free-surface condition
## Loop over all frequencies using @threads!

function applyfreesurf(mat_ref_surf, list_tran_freq, list_ref_freq)
    numfreq = length(list_tran_freq)
    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        mat_tran_up = list_tran_freq[ind]
        mat_ref_down = list_ref_freq[ind]

        mat_tran_up = 1.0I(3)/(1.0I(3)-mat_ref_down*mat_ref_surf)*mat_tran_up
        list_tran_freq[ind] = mat_tran_up
    end

    return list_tran_freq
end

function applyfreesurf!(mat_ref_surf, list_tran_freq, list_ref_freq)
    numfreq = length(list_tran_freq)
    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        mat_tran_up = list_tran_freq[ind]
        mat_ref_down = list_ref_freq[ind]

        mat_tran_up[:] = 1.0I(3)/(1.0I(3)-mat_ref_down*mat_ref_surf)*mat_tran_up
        list_tran_freq[ind] = mat_tran_up
    end
end

## Apply the free-surface condition for the P-SV system

function applyfreesurf_psv!(mat_ref_surf, list_tran_freq, list_ref_freq)
    numfreq = length(list_tran_freq)
    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        mat_tran_up = list_tran_freq[ind]
        mat_ref_down = list_ref_freq[ind]

        mat_tran_up[:] = 1.0I(2)/(1.0I(2)-mat_ref_down*mat_ref_surf)*mat_tran_up
        list_tran_freq[ind] = mat_tran_up
    end
end

## Apply the free-surface condition for the SH system

function applyfreesurf_sh!(ref_surf, list_tran_freq, list_ref_freq)
    numfreq = length(list_tran_freq)
    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        tran_up = list_tran_freq[ind]
        ref_down = list_ref_freq[ind]

        tran_up = 1.0/(1.0-ref_down*ref_surf)*tran_up
        list_tran_freq[ind] = tran_up
    end
end

## Get the displacements at the free surface
## Loop over all frequencies using @threads!

function getfreesurfdisp(flag_inc, mat_ref_surf, mat_dsst_surf, list_tran_freq)
    numfreq = length(list_tran_freq)
    mat_disp_freq = zeros(Complex, 3, numfreq)

    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        # Compute the wave vectors of the up-going waves
        mat_tran_up = list_tran_freq[ind]
        if flag_inc == 1
            vect_up = mat_tran_up[:, 1]
        elseif flag_inc == 2
            vect_up = mat_tran_up[:, 2]
        elseif flag_inc == 3
            vect_up = mat_tran_up[:, 3]
        end

        # Compute the wave vectors of the down-going waves
        vect_down = mat_ref_surf*vect_up

        # Compute the displacement vector at the free surface
        mat_disp = mat_dsst_surf[1:3, :];
        vect_disp = mat_disp*[vect_down; vect_up];
        mat_disp_freq[:, ind] = vect_disp
    end

    return mat_disp_freq
end

## Get the P-SV displacements at the free surface

function getfreesurfdisp_psv(flag_inc, mat_ref_surf, mat_dsst_surf, list_tran_freq)
    numfreq = length(list_tran_freq)
    mat_disp_freq = zeros(Complex, 2, numfreq)

    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        # Compute the wave vectors of the up-going waves
        mat_tran_up = list_tran_freq[ind]
        if flag_inc == 1
            vect_up = mat_tran_up[:, 1]
        elseif flag_inc == 2
            vect_up = mat_tran_up[:, 2]
        end

        # Compute the wave vectors of the down-going waves
        vect_down = mat_ref_surf*vect_up ### Not sure if this is correct! 2023-02-20

        # Compute the displacement vector at the free surface
        mat_disp = mat_dsst_surf[1:2, :];
        vect_disp = mat_disp*[vect_down; vect_up];
        mat_disp_freq[:, ind] = vect_disp
    end

    return mat_disp_freq
end

## Get the SH displacements at the free surface

function getfreesurfdisp_sh(ref_surf, mat_dsst_surf, list_tran_freq)
    numfreq = length(list_tran_freq)
    vect_disp_freq = zeros(Complex, numfreq)

    #Threads.@threads for ind in 1:numfreq
    for ind in 1:numfreq
        # Compute the wave vectors of the up-going waves
        tran_up = list_tran_freq[ind]
        wave_up = tran_up

        # Compute the wave vectors of the down-going waves
        wave_down = ref_surf*wave_up 

        # Compute the displacement vector at the free surface
        vect_disp = mat_dsst_surf[1, :]
        disp = vect_disp'*[wave_down; wave_up]
        vect_disp_freq[ind] = disp
    end

    return vect_disp_freq
end

## Apply a phase correction to shift the first arrival to time 0

function applyphasecorr(mat_disp, vect_phs)
    mat_disp[1, :] = mat_disp[1, :].*vect_phs
    mat_disp[2, :] = mat_disp[2, :].*vect_phs
    mat_disp[3, :] = mat_disp[3, :].*vect_phs

    return mat_disp
end

function applyphasecorr!(mat_disp, vect_phs)
    mat_disp[1, :] = mat_disp[1, :].*vect_phs
    mat_disp[2, :] = mat_disp[2, :].*vect_phs
    mat_disp[3, :] = mat_disp[3, :].*vect_phs
end

## Apply a phase correction to shift the first arrival of the P-SV system to time 0

function applyphasecorr_psv!(mat_disp, vect_phs)
    mat_disp[1, :] = mat_disp[1, :].*vect_phs
    mat_disp[2, :] = mat_disp[2, :].*vect_phs
end

## Get the response of a stack of anisotropic layers to an incident wave of P (1), SV (2), or SH (3)
## The bottom layer of the model must be isotropic!
## The input elastic tensor needs to be in N-E-D coordinate sysgtem!
## The output traces are in the order of N, E, and D!
## 
## flag_inc: Incident wave type, 1 = P, 2 = SV, 3 = SH
## flag_verb: 0 = No standard output, 1 = Standard output
## rayp: Horizontal slowness of the incident wave (s/km)
## baz: Back azimuth of the incident wave (clockwise from the north)
## list_h: List of layer thicknesses
## list_rho: List of layer densities
## list_voi: List of Voigt matrices
## npts: Number of time points
## dt: Sampling interval of the output waveform
## t_b: Starting time of the time series

function getaniresp(flag_inc, flag_verb, hslow, baz, list_h, list_rho, list_ela, npts, dt, t_b)
    # Determine if the input elastic parameters are in matrix or tensor format
    numlyr = length(list_h)

    if size(list_ela[1]) == (6, 6)
        list_tsr = [zeros(3, 3, 3, 3) for ind in 1:numlyr]
        for ind = 1:numlyr
            mat_voi = list_ela[ind]
            tsr = voigt2tensor(mat_voi)
            list_tsr[ind] = tsr
        end
    else
        list_tsr = list_ela
    end

    # Compute the slownesses and polarizatin vectors for all layers

    if flag_verb > 0
        print("Computing the slownesses and polarizations for each layer...\n")
    end

    list_vslow, list_dsst = getallslowandpol(rayp, baz, list_rho, list_tsr, zeros(3))
        
    # print(list_vslow[1])
    # println()
    # print(list_vslow[2])
    # println()
    # print(list_vslow[3])
    # println()
    # dispmatrix(list_dsst[1])
    # println()
    # dispmatrix(list_dsst[2])
    # println()
    # dispmatrix(list_dsst[3])

    if flag_verb > 0
        print("Finished\n")
    end

    # Compute the free-surface matrix
    if flag_verb > 0
        print("Computing the free-surface matrix...\n")
    end

    mat_dsst_surf = list_dsst[1]
    mat_ref_surf = getfreesurfmat(mat_dsst_surf)

    if flag_verb > 0
        print("Finished\n")
    end

    # Compute the reflection and transmission matrices at each interface
    if flag_verb > 0
        print("Computing the reflection and transmission matrices at all interfaces...\n")
    end

    list_ref_down, list_tran_down, list_ref_up, list_tran_up = getallintmats(list_dsst)

    # dispmatrix(abs.(list_ref_down[1]))
    # println()
    # dispmatrix(abs.(list_tran_down[1]))
    # println()
    # dispmatrix(abs.(list_ref_down[2]))
    # println()
    # dispmatrix(abs.(list_tran_down[2]))

    if flag_verb > 0
        print("Finished\n")
    end
    
    # Propagate the wavefield up the stack
    if flag_verb > 0
        print("Propagating the wavefield up the stack...\n")
    end

    # Define the frequency vector
    fs = 1/dt
    grid_freq = fftfreq(npts, fs)
    grid_freq_pos = grid_freq[1:div(npts, 2)+1] # Bug with the length of FFT is fixed, 2023-04-16
    grid_omega_pos = 2.0pi*grid_freq_pos
    numfreq_pos = length(grid_freq_pos)

    list_tran_freq = [zeros(Complex, 3, 3) for ind in 1:numfreq_pos]
    list_ref_freq = [zeros(Complex, 3, 3) for ind in 1:numfreq_pos]
    vect_corr = zeros(Complex, numfreq_pos)

    propupstack!(list_tran_freq, list_ref_freq, vect_corr, grid_omega_pos, list_h, list_vslow, list_ref_down, list_tran_down, list_ref_up, list_tran_up, flag_inc)

    if flag_verb > 0
        print("Finished.\n")
    end

    # for ind = 1:100
    #     print(abs(list_tran_freq[ind][1, 1]))
    #     println()
    # end

    # Apply the free surface condition
    if flag_verb > 0
        print("Applying the free-surface condition...\n")
    end

    #list_tran_freq = applyfreesurf(mat_ref_surf, list_tran_freq, list_ref_freq)
    applyfreesurf!(mat_ref_surf, list_tran_freq, list_ref_freq)

    # for ind = 1:numfreq_pos
    #     print(maximum(abs.(list_tran_freq[ind])))
    #     println()
    # end

    if flag_verb > 0
        print("Finished.\n")
    end

    # Compute the displacements at the free surface
    if flag_verb > 0
        print("Computing the displacements at the free surface...\n")
    end

    mat_disp_posfreq = getfreesurfdisp(flag_inc, mat_ref_surf, mat_dsst_surf, list_tran_freq)

    if flag_verb > 0
        print("Finished.\n")
    end

    # print(size(mat_disp_posfreq))
    # println()
    # print(maximum(abs.(mat_disp_posfreq[1, :])))
    # println()
    # print(maximum(abs.(mat_disp_posfreq[2, :])))
    # println()
    # print(maximum(abs.(mat_disp_posfreq[3, :])))

    # Apply a phase shift to place the first arrival at time 0
    if flag_verb > 0
        print("Applying a phase shift to place the first arrival at time 0...\n")
    end

    vect_corr = vect_corr.*exp.(1.0im*grid_omega_pos*t_b)
    #mat_disp_posfreq = applyphasecorr(mat_disp_posfreq, vect_phs)
    applyphasecorr!(mat_disp_posfreq, vect_corr)

    if flag_verb > 0
        print("Finished.\n")
    end

    # Apply inverse fourier transform to get the response in time domain
    if flag_verb > 0
        print("Performing inverse fourier transform...\n")
    end

    println(size(Complex.(mat_disp_posfreq)))
    println(npts)
    mat_disp = irfft(Complex.(mat_disp_posfreq), npts, 2)

    if flag_verb > 0
        print("Finished.\n")
    end

    # Rotate the displacements back to N-E-D
    if flag_verb > 0
        print("Rotating the displacements back to N-E-D...\n")
    end
    
    mat_disp_out = rotatecoord_matrix(3, -baz-180, mat_disp);

    if flag_verb > 0
        print("Finished.\n")
    end

    return mat_disp_out
end


## A special case for getaniresp that handles VTI mediums, in which P-SV and SH waves are decoupled

function getaniresp_vti(flag_inc, flag_verb, hslow, baz, list_h, list_rho, list_ela, npts, dt, t_b)
    # Determine if the input elastic parameters are in matrix or tensor format
    numlyr = length(list_h)

    if size(list_ela[1]) == (6, 6)
        list_tsr = [zeros(3, 3, 3, 3) for ind in 1:numlyr]
        for ind = 1:numlyr
            mat_voi = list_ela[ind]
            tsr = voigt2tensor(mat_voi)
            list_tsr[ind] = tsr
        end
    else
        list_tsr = list_ela
    end

    # Compute the slownesses and polarizatin vectors for all layers
    

    if flag_verb > 0
        if flag_inc == 1
            println("The incident wave is P.")
        elseif flag_inc == 2
            println("The incident wave is SV.")
        else
            println("The incident wave is SH.")
        end

        println("Computing the slownesses and polarizations for each layer...")
    end

    list_vslow_psv, list_vslow_sh, list_dsst_psv, list_dsst_sh = getallslowandpol_vti(rayp, baz, list_rho, list_tsr)

    # dispmatrix(list_dsst_sh[1])
    # println()
    # dispmatrix(list_dsst_sh[2])
    # println()
        
    if flag_verb > 0
        println("Finished.")
    end

    # Compute the free-surface matrix
    if flag_verb > 0
        println("Computing the free-surface matrix...")
    end

    if flag_inc == 1 || flag_inc == 2
        mat_dsst_surf = list_dsst_psv[1]
        mat_ref_surf = getfreesurfmat_psv(mat_dsst_surf)
    else
        mat_dsst_surf = list_dsst_sh[1]
        ref_surf = getfreesurfmat_sh(mat_dsst_surf)        
    end


    if flag_verb > 0
        println("Finished")
    end

    # Compute the reflection and transmission matrices at each interface
    if flag_verb > 0
        println("Computing the reflection and transmission matrices at all interfaces...")
    end

    if flag_inc == 1 || flag_inc == 2
        list_ref_down, list_tran_down, list_ref_up, list_tran_up = getallintmats_psv(list_dsst_psv)
    else
        list_ref_down, list_tran_down, list_ref_up, list_tran_up = getallintmats_sh(list_dsst_sh)
    end

    # print(list_ref_down[1])
    # println()
    # print(list_tran_down[1])
    # println()
    # print(list_ref_down[2])
    # println()
    # print(list_tran_down[2])
    # println()


    if flag_verb > 0
        println("Finished")
    end
    
    # Propagate the wavefield up the stack
    if flag_verb > 0
        println("Propagating the wavefield up the stack...")
    end

    # Define the frequency vector
    fs = 1/dt
    grid_freq = fftfreq(npts, fs)
    grid_freq_pos = grid_freq[grid_freq .>= 0.0]
    grid_omega_pos = 2.0pi*grid_freq_pos
    numfreq_pos = length(grid_freq_pos)

    vect_corr = zeros(Complex, numfreq_pos)

    if flag_inc == 1 || flag_inc == 2
        list_tran_freq = [zeros(Complex, 2, 2) for ind in 1:numfreq_pos]
        list_ref_freq = [zeros(Complex, 2, 2) for ind in 1:numfreq_pos]

        propupstack_psv!(list_tran_freq, list_ref_freq, vect_corr, grid_omega_pos, list_h, list_vslow_psv, list_ref_down, list_tran_down, list_ref_up, list_tran_up)
    else
        list_tran_freq = zeros(Complex, numfreq_pos)
        list_ref_freq = zeros(Complex, numfreq_pos)

        propupstack_sh!(list_tran_freq, list_ref_freq, vect_corr, grid_omega_pos, list_h, list_vslow_sh, list_ref_down, list_tran_down, list_ref_up, list_tran_up)
    end

    # print(abs.(list_tran_freq))

    if flag_verb > 0
        println("Finished.")
    end

    # Apply the free surface condition
    if flag_verb > 0
        println("Applying the free-surface condition...")
    end

    if flag_inc == 1 || flag_inc == 2
        applyfreesurf_psv!(mat_ref_surf, list_tran_freq, list_ref_freq)
    else
        # print(ref_surf)
        applyfreesurf_sh!(ref_surf, list_tran_freq, list_ref_freq)
    end

    # print(abs.(list_tran_freq))

    if flag_verb > 0
        println("Finished.")
    end

    # Compute the displacements at the free surface
    if flag_verb > 0
        println("Computing the displacements at the free surface...")
    end

    if flag_inc == 1 || flag_inc == 2
        mat_disp_posfreq = getfreesurfdisp_psv(flag_inc, mat_ref_surf, mat_dsst_surf, list_tran_freq)
    else
        vect_disp_posfreq = getfreesurfdisp_sh(ref_surf, mat_dsst_surf, list_tran_freq)
    end

    # print(abs.(vect_disp_posfreq))

    if flag_verb > 0
        println("Finished.")
    end

    # print(size(mat_disp_posfreq))
    # println()
    # print(maximum(abs.(mat_disp_posfreq[1, :])))
    # println()
    # print(maximum(abs.(mat_disp_posfreq[2, :])))
    # println()
    # print(maximum(abs.(mat_disp_posfreq[3, :])))

    # Apply a phase shift to place the first arrival at time 0
    if flag_verb > 0
        println("Applying a phase shift to place the first arrival at time 0...")
    end

    vect_corr = vect_corr.*exp.(1.0im*grid_omega_pos*t_b)
    #mat_disp_posfreq = applyphasecorr(mat_disp_posfreq, vect_phs)

    if flag_inc == 1 || flag_inc == 2
        applyphasecorr_psv!(mat_disp_posfreq, vect_corr)
    else
        vect_disp_posfreq = vect_disp_posfreq.*vect_corr
    end

    if flag_verb > 0
        println("Finished.")
    end

    # Apply inverse fourier transform to get the response in time domain
    if flag_verb > 0
        println("Performing inverse fourier transform...")
    end

    if flag_inc == 1 || flag_inc == 2
        mat_disp = irfft(Complex.(mat_disp_posfreq), npts, 2)
    else
        vect_disp = irfft(Complex.(vect_disp_posfreq), npts)
    end

    if flag_verb > 0
        println("Finished.")
    end

    # Rotate the displacements back to N-E-D
    if flag_verb > 0
        println("Rotating the displacements back to N-E-D...")
    end
    
    if flag_inc == 1 || flag_inc == 2
        mat_disp = [mat_disp[1, :]'; zeros(npts)'; mat_disp[2, :]';]
        mat_disp_out = rotatecoord_matrix(3, -baz-180, mat_disp)
    else
        mat_disp_out = zeros(3, npts)
        mat_disp_out[1, :] = vect_disp.*sind(baz)
        mat_disp_out[2, :] = -vect_disp.*cosd(baz)
        mat_disp_out[3, :] = zeros(npts)
    end

    if flag_verb > 0
        println("Finished.")
    end

    return mat_disp_out
end

## Get the response of the SH system

function getshresp(flag_verb, hslow, list_h, list_vs, list_rho, npts, dt, t_b)
    numlyr = length(list_h)

    # Get all the displacement-stress matrices
    list_vslow, list_dsst = getalldsst_sh(rayp, list_vs, list_rho)
        
    if flag_verb > 0
        println("Finished.")
    end

    # Compute the free-surface matrix
    if flag_verb > 0
        println("Computing the free-surface matrix...")
    end

    mat_dsst_surf = list_dsst[1]
    ref_surf = getfreesurfmat_sh(mat_dsst_surf)        

    if flag_verb > 0
        println("Finished")
    end

    # Compute the reflection and transmission matrices at each interface
    if flag_verb > 0
        println("Computing the reflection and transmission matrices at all interfaces...")
    end

    list_ref_down, list_tran_down, list_ref_up, list_tran_up = getallintmats_sh(list_dsst)

    if flag_verb > 0
        println("Finished")
    end
    
    # Propagate the wavefield up the stack
    if flag_verb > 0
        println("Propagating the wavefield up the stack...")
    end

    # Define the frequency vector
    fs = 1/dt
    grid_freq = fftfreq(npts, fs)
    grid_freq_pos = grid_freq[grid_freq .>= 0.0]
    grid_omega_pos = 2.0pi*grid_freq_pos
    numfreq_pos = length(grid_freq_pos)

    vect_corr = zeros(Complex, numfreq_pos)

    list_tran_freq = zeros(Complex, numfreq_pos)
    list_ref_freq = zeros(Complex, numfreq_pos)

    propupstack_sh!(list_tran_freq, list_ref_freq, vect_corr, grid_omega_pos, list_h, list_vslow, list_ref_down, list_tran_down, list_ref_up, list_tran_up)

    # print(abs.(list_tran_freq))

    if flag_verb > 0
        println("Finished.")
    end

    # Apply the free surface condition
    if flag_verb > 0
        println("Applying the free-surface condition...")
    end

    # print(ref_surf)
    applyfreesurf_sh!(ref_surf, list_tran_freq, list_ref_freq)

    # print(abs.(list_tran_freq))

    if flag_verb > 0
        println("Finished.")
    end

    # Compute the displacements at the free surface
    if flag_verb > 0
        println("Computing the displacements at the free surface...")
    end

    vect_disp_posfreq = getfreesurfdisp_sh(ref_surf, mat_dsst_surf, list_tran_freq)

    if flag_verb > 0
        println("Finished.")
    end

    # Apply a phase shift to place the first arrival at time 0
    if flag_verb > 0
        println("Applying a phase shift to place the first arrival at time 0...")
    end

    vect_corr = vect_corr.*exp.(1.0im*grid_omega_pos*t_b)
    vect_disp_posfreq = vect_disp_posfreq.*vect_corr

    if flag_verb > 0
        println("Finished.")
    end

    # Apply inverse fourier transform to get the response in time domain
    if flag_verb > 0
        println("Performing inverse fourier transform...")
    end

    vect_disp = irfft(Complex.(vect_disp_posfreq), npts)

    if flag_verb > 0
        println("Finished.")
    end

    return vect_disp
end