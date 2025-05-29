!	input and open files 

	namelist/input_which/ input_EQDSK, EQDSK_file

	namelist/input_triangularity/ tri_type, n_tri, theta_points_max,  &
								  a_elps, b_elps, k_ellipt, delta_up,  &
								  delta_down

	namelist/input_constants/  rmajor, rcenter, zcenter, x_size , z_size,  &
							grid_type,  &
							eq_type, data_dim, numerical,   &
							Broot, mass, me_ov_mi, pe_ov_p, eV

	namelist/input_flow/ mach_phi_max, alpha_mphi, mphi_min, &
						 mach_theta_max, mach_theta_edge, t_mth, w_mth, enforce_static,  &
						 omega_0

	namelist/input_magnetic/  b_phi_zero, F_opt, Fc_o_Fv, kappa, eta_P, mu_RFP,  &
							  mu_mag, psi_e, inter_switch,  Z_eff_0, delta_Z_eff, Z_eff_exp

	namelist/input_p_d_profile/ gamma, gamma_i, gamma_e, dcenter, de_o_dc, alpha, alpha_e, alpha_i,  &
													alpha_rho, beta_center, pe_o_pc, beta_e_center, pe_o_pc_e,  &
													beta_i_center, pe_o_pc_i, p_opt, p_e_opt, p_i_opt
  
	namelist/input_numerical/  numerical_opt, numerical_n, numerical_p_iso,   &
							   numerical_p_e, numerical_p_i, numerical_F,  &
							   numerical_omega, omega_option, numerical_mtheta,  &
							   numerical_psiprim, numerical_psi_diff

	namelist/input_solver/ n_min, n, min_it, max_it, accelerate, accelerate_0, fix_orp, fix_orp_0,  &
						   bc_type, max_orp, bc_switch, restart,  &
						   jump_option, lambda0_fix, dir_switch, delta_Bern_fact, psi_diff_option

	namelist/input_gravity/ gravity_type, G_gravity, M_gravity, Kepler_Omega

	namelist/input_output/	write_all, write_all_bin, write_TF_roots, q_vnorm, vnorm_label


!------------------------------------------------------------------------

!	open(unit=5,file='inputfile.dat',status='old')
