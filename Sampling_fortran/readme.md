!###############################################################################################################
!## Here we have presented a code for particle sampling from a 3D-thermal freeze-out surface (From boost-
!## invariant **MUSIC** output, filename="suface.dat"). 
!##
!## To Run : gfortran particle_at_freezeout_step3.f90 followed by ./a.out 
!##
!###############################################################################################################
!## Our approach is as follows :
!##       1. Read the information of all surface cells (make corresponding arrays):
!##                tau_f,x_f,y_f,eta_f,d3sigma_tau,d3sigma_x,d3sigma_y,d3sigma_eta,u_tau,u_x,u_y,uetat,ef,Tf
!##                ,mu_B,ef_Plus_Pf_by_Tf,pi_tautau,pi_taux,pi_tauy,pi_tauetat,pi_xx,pi_xy,pi_xetat,pi_yy,pi_yetat
!##                ,pi_etatetat
!##       2. Read the mass, sign_B(boson/fermion) from "my_file_mass.txt", spin and name of the particle species 
!##          we want to sample.
!##
!##       3. Calculate the corresponding degeneracy and the pre-factor.
!##       4. Calculate the p_T, phi differential spectra using the formula :
!##
!##           diff_spec(p_T,phi)= Sum_{over all surface cells} (f+delta_f)*psigma_dsigma
!##
!##           where, f=1./(exp(1./Tf*(E - mu_B)) + sign_B) and delta_f= (f*(1.sign_B*f)*pre_factor_shear*w_factor)
!##
!##       5. Integrate over p_T, phi to obtain the total number of particles (of a specified species).
!##       6. Use the Poisson distribution to sample the total number of particles (of a specified species) in
!##          single event.
!##       7. Sample y, p_T, and phi of each hadron using the acceptance-rejection sampling method as follows :
!##              
!##                   Probability = diff_spec(p_T,phi)/max_diff_spec(p_T_max,phi_max)  < 1
!##              
!##                   where, max_diff_spec(p_T_max,phi_max) is the maximum value of diff_spec occured at the 
!##                   p_T_max,phi_max.
!##
!##         Note: As we consider only the boost-invariant dynamics, we randomly sampled the y (rapidity).
!##
!##       8.  Sample the surface cell information for each particle as follows :
!##
!##                  A. Randomly take a surface cell, say 'i'th.
!##                  B. Use the following acceptance-rejection method to check whether the 'i'th surface cell can
!##                     emit the particle :
!##
!##                     Probability = diff_spec(i,p_T,phi)/max_diff_spec(i,p_T_max,phi_max)  < 1
!##
!##                   where, max_diff_spec(i,p_T_max,phi_max) is the maximum value of diff_spec occured at the 
!##                   p_T_max,phi_max for 'i'th surface cell.
!##       9.   Finally, print the x_f,y_f, eta_f,tau_f, y, p_T*cos(phi) and p_T*sin(phi) of each particle.
!##         
!##    Note: A. The sampling can be much faster if the process can be parallelized.
!##           B. The sampling method used above is a possible approach (by the author). Please follow the Ref. 
!##              "Pseudorapidity distribution and decorrelation of anisotropic flow within CLVisc
!##               hydrodynamics" by Long-Gang Pang, Hannah Petersen, and Xin-Nian Wang [arXiv: 1802.04449v2] for
!##               more appropriate approach where sampling of particles are done for each surface cell. 
!##
