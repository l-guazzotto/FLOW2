! modified 9/3/2013

module solver

	use constant
	use magnetic
	use p_d_profile
	use flow
	use triangularity
	use exp_data
	use interpolating_functions
	use pseudo_IMSL, only : dbsval, dbsder

	implicit none

	!  private :: dofpsi, dddpsi
	!  private :: pofpsi, dpdpsi
	!  private :: bzero, dbzerodpsi
	private :: dddpsi
	!  private :: dpdpsi
	!  private :: dbzerodpsi
	private :: mach_theta, dmach_thetadpsi
	private :: mach_phi, dmach_phidpsi
	private :: sofpsi, dsdpsi, hofpsi, dhdpsi, iofpsi, didpsi
!	private :: omegaofpsi, domegadpsi, phiofpsi, dphidpsi
	private :: domegadpsi, phiofpsi, dphidpsi
	private :: bpol

	real (kind=dkind), private :: bhut_max,bhut_min
	real (kind=dkind) :: psi_pres

	! ------------------------------------------------------------------
	! stuff to speed up functions calculation

	real(kind=dkind), private :: psic_flag = 1.d6
	real(kind=dkind), private :: psi_flag = 1.d9
	real(kind=dkind), private :: psi_flag_dep = 1.d9
	real(kind=dkind), private :: psi_flag_dep2 = 1.d9
	real(kind=dkind), private :: psi_flag_ham = 1.d9

	real(kind=dkind), private :: d_loc, dp_loc, p_loc, pp_loc, b0_loc,  &
												b0p_loc, mth_loc, mthp_loc, mph_loc,  &
												mphp_loc
	real(kind=dkind), private :: s_loc, sp_loc, phi_loc, phip_loc, i_loc,  &
													ip_loc, omega_loc, omegap_loc, h_loc, hp_loc
	real(kind=dkind), private :: dc_loc, dcp_loc, pc_loc, pcp_loc,  &
													b0c_loc, b0cp_loc
	real(kind=dkind) :: Bernmax, fBernmax, delta_Bern, psi_Bern, big_Psi_Bern

	! ------------------------------------------------------------------
	! The following variables are defined to allow the bernoulli
	! function and such to be used and its parameters set by another
	! function in this module.

	real (kind=dkind), private :: g_Phi, g_r, g_dpsidx, g_dpsidz
	real (kind=dkind), private :: g_dpsidx_L, g_dpsidz_L, g_dpsidx_R, g_dpsidz_R
	real (kind=dkind), private :: g_I, g_Omega, g_S, g_H
	integer, private :: g_indi, g_indj, g_nx, g_nz
	real (kind=dkind), private :: g_dx, g_dz, g_s_e, g_s_i, g_mtheta
	real (kind=dkind) :: g_b_polc, g_b_torc, g_bfield
	real (kind=dkind), private :: g_D, g_Lambda
	real(kind=dkind) :: g_gPsi2, g_H_e, g_H_i, g_H_diff, g_psi_diff, g_H_i_prime,  &
									g_S_i_prime, g_Fstar, g_Dstar_term
	real(kind=dkind), allocatable, dimension(:) :: x_int_H, fun_int_H, w_int_H	! fun_int must also include the Jacobian
	real(kind=dkind), allocatable, dimension(:) :: x_int_phi, fun_int_phi, w_int_phi	! fun_int must also include the Jacobian

	integer, private :: m_bound

	logical, private :: inside,outside
	real(kind=dkind) :: bpol_max, bpol_min
	logical, private :: reduce=.false.
	integer :: i_bmin, j_bmin

	! stuff for the grid
	real(kind=dkind) :: g_ratio
	integer :: n_temp
	real(kind=dkind), dimension(:,:), allocatable, public :: fmax_2D
	real(kind=dkind) :: psi_max
	logical :: write_TF_roots

	real(kind=dkind) :: grad_ratio_0 = 20.d0	!20.d0
	real(kind=dkind) :: xT1, xT2, xT3, zT1, zT2, zT3, two_A_tri
	logical :: D_TF_num = .false.	!whether to use TF numerical density in SF calculation
	integer :: id_Bern = 73
	character(len=8) :: Bern_format =  '(I5.5)' ! format descriptor

	integer :: psi_diff_option = 2 ! option for solving the GS equarion for delta_psi/big_Psi
	integer :: ipsic, jpsic ! location of psi_max

!fmt = '(I5.5)'

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : bc_setup_option

	implicit none

	integer :: n

	if(bc_type>100) then
		call initialize_bc_tri_equations(n)
		return
	endif

	if(bc_setup_option>0) then
		call initialize_bc_equations(n)
	elseif(bc_setup_option<0) then
		call initialize_bc_minimum(n)
	endif

end subroutine initialize_bc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc_equations(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Aug 15 2013: adding psi_diff_bound

	use constant, only : dx, dz, x_coord, z_coord, dx_a, dz_a, bc_switch
	use pseudo_IMSL, only : DBSNAK, DBSINT

	integer :: n

	integer, dimension(n*n,4) :: ind
	real(kind=dkind), dimension(n*n,6) :: coord
	! third number is distance ratio
	real(kind=dkind), dimension(n*n) :: psi_diff_temp

	real(kind=dkind) :: ex,ez, dist0,dist
	real(kind=dkind) :: x,z, RP, zP, R0, xb, zb, th
	real(kind=dkind) ::  xs, zs ! points on the surface
	real(kind=dkind), dimension(1:2) :: xvec, fvec 
	! solution points and function values from IMSL
	real(kind=dkind), dimension(1:2) :: xguess
	real(kind=dkind) :: fnorm
	external  eq_bound, eq_bound2
	integer :: ntrial = 5000
	real(kind=dkind) :: tolx=1.d-12
	real(kind=dkind) :: tolf = 1.d-9
	integer :: zone
	real(kind=dkind) :: dummy(1:7)
	real(kind=dkind) :: r
	real(kind=dkind) :: omega_temp
	real(kind=dkind) :: dref, theta_one, theta_two, theta_three, theta_four, thetaloc, Rloc,  &
									delta_psi

	integer :: i,j

	logical :: bound, inner, truebound, zero_dist
	logical :: newton_check

	if(bc_type==17) then
	! This needs to happen here to avoid cross references

		! big_Psi = psi on the boundary + toroidal flow term assigned through omega_0 and d(Psi)

		! We need to determine the limiting angles corresponding to the corners of the domain first.

		theta_one = atan(b_elps/a_elps) ! top right/out
		theta_two = pi - theta_one ! top left/in
		theta_three = pi + theta_one ! bottom left/in
		theta_four = pi + theta_two ! bottom right/out

		big_Psib_dim =psib_dim
		ibreak_bPb = big_Psib_dim + big_Psib_ord

		allocate(big_Psib_cscoef(big_Psib_dim))

		big_Psib_data = psib_data

		dref = d_TF_ofpsi(0.d0)

		do i = 1,big_Psib_dim

			thetaloc = big_Psib_data(i,1)

			if(thetaloc<=theta_one) then

				Rloc = rmajor + a_elps

			elseif(thetaloc<=theta_two) then

				Rloc = rmajor + b_elps/tan(thetaloc)

			elseif(thetaloc<=theta_three) then

				Rloc = rmajor - a_elps

			elseif(thetaloc<=theta_four) then

				Rloc = rmajor - b_elps/tan(thetaloc)

			else

				Rloc = rmajor + a_elps

			endif

			delta_psi = mass*Rloc**2*omega_0/eV
			big_Psib_data(i,2) = big_Psib_data(i,2)  + delta_psi

		enddo

		!this sets up the knots
		call DBSNAK(big_Psib_dim, big_Psib_data(1:big_Psib_dim,1),  &
						big_Psib_ord,big_Psib_data(1:ibreak_bPb,3))

		call DBSINT(big_Psib_dim, big_Psib_data(1:big_Psib_dim,1),  &
			 big_Psib_data(1:big_Psib_dim,2), psib_ord,  &
			 big_Psib_data(1:ibreak_bPb,3),  &
			 big_Psib_cscoef(1:big_Psib_dim))

			 bc_type = 7

	endif

	! For now, explicitly shorting free-boundary as that's what worked in FLOW. 
	! But also letting the bc_type==17 condition do stuff above as a test
	! TODO: check whether this works, and if so remove all other bc_type==7 checks below as redundant
	if(bc_type==7) then 
		return
	endif

	! Note: this is coded opposite to FLOW; instead of shorting the function for some conditions,
	! we only allow specific conditions through to the rest of the function.
	if((bc_type==3).or.(bc_type==23).or.(((bc_type==4).or.(bc_type==5).or.(bc_type==7).or.  &
		(bc_type==8).or.(bc_type==14).or.(bc_type==24).or.(bc_type==34).or.(bc_type==44).or.  &
		(bc_type==54).or.(bc_type==64))  &
		.and.(n>=bc_switch))) then

		continue

	else

		return

	endif

!	if( (bc_type<3).or.( ((bc_type==4).or.(bc_type==5).or.(bc_type==7).or.  &
!		(bc_type==8).or.(bc_type==14).or.(bc_type==24)).and.(n<bc_switch)).or.(bc_type/=23) ) return

	if(allocated(ind_bound)) deallocate(ind_bound)
	if(allocated(coord_bound)) deallocate(coord_bound)
	if(allocated(psi_diff_bound)) deallocate(psi_diff_bound)

	m_bound = 0

	coord = 0.d0
	
	omega_temp = omega_0

	do j=1, n
       do i=1, n

			zero_dist = .false.

			call check_position(i,j,bound,truebound,n)

			if(truebound) then
			! real boundary point

				m_bound = m_bound+1

				call radius(i,j,n,n,ex,ez,r,dx,dz)

				if(ex**2+ez**2<r**2) then
				! the point is very close to the boundary

					R_P = x_coord(i)
					z_P = z_coord(j)

					xguess(1) = R_P
					xguess(2) = z_P

					xvec(1) = xguess(1)
					xvec(2) = xguess(2)

					call newt(ntrial,xvec,2,newton_check,tolx,tolf)

					if((bc_type==7).or.(bc_type==8)) then

						xb = xvec(1)
						zb = xvec(2)
						th = atan2(zb,xb-rmajor)
						if(th<0.d0) th = th + 2.d0*pi

						coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
										ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

					endif

					dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

					if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

					else

						! proceed as usual

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

						ind(m_bound,1) = i
						ind(m_bound,2) = j

						inner = .false.

						coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
						coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

						call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
									ind(m_bound,3:4) )

						call check_points(ind(m_bound,3),ind(m_bound,4),inner)

						if(inner) then

							coord(m_bound,3) = 1.d0

						else

							coord(m_bound,3) = -1.d0

						endif

						coord(m_bound,4) = 0.d0 !distance from point to surface
						coord(m_bound,5) = dist !distance from surface to interpolated point

						if((bc_type==23).or.(bc_type==24)) then

							if(numerical_psi_diff) then

								ex = xvec(1) - rmajor
								ez = xvec(2)

								if (ex==0.d0) then
									theta = pi/2.d0 * dsign(1.d0,ez)
								else
									theta = datan2(ez,ex)
								endif

								if(theta<0.d0) then
									theta = theta + 2.d0*pi
								endif

								psi_diff_temp(m_bound) = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
																		ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )
								
							else

								psi_diff_temp(m_bound) = xvec(1)**2*omega_temp*mass/eV

							endif

						endif

						if(bc_type==5) then

							call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
											dummy(3),zone,dummy(4),dummy(5),dummy(6))
							! NOTE: this works only because the inner boundary is a circle,
							! in the general case one should use the coordinates of the 
							! point on the boundary

						endif

					endif

				else
				! the point REALLY is on the boundary (numerically),
				! this is extremely unlikely to occur and dealt with easily

					zero_dist = .true.

					continue

				endif

			elseif(bound) then
				! (external) boundary point

				m_bound = m_bound+1

				R_P = x_coord(i)
				z_P = z_coord(j)

				xguess(1) = R_P
				xguess(2) = z_P

				xvec(1) = xguess(1)
				xvec(2) = xguess(2)

				call newt(ntrial,xvec,2,newton_check,tolx,tolf)

				if((bc_type==7).or.(bc_type==8)) then

					xb = xvec(1)
					zb = xvec(2)
					th = atan2(zb,xb-rmajor)
					if(th<0.d0) th = th + 2.d0*pi

					coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
									ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

				endif

				dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

				if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

				else

					if(dist0<0.2d0*sqrt(dx_a(i)**2+dz_a(j)**2)) then

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

					else

						dist = dist0

					endif

					ind(m_bound,1) = i
					ind(m_bound,2) = j

					inner = .false.

					coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
					coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

					call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
								ind(m_bound,3:4) )

					call check_points(ind(m_bound,3),ind(m_bound,4),inner)

					if(inner) then

						coord(m_bound,3) = 1.d0

					else

						coord(m_bound,3) = -1.d0

					endif

					coord(m_bound,4) = dist0 !distance from point to surface
					coord(m_bound,5) = dist !distance from surface to interpolated point

					if((bc_type==23).or.(bc_type==24)) then

						if(numerical_psi_diff) then

							ex = xvec(1) - rmajor
							ez = xvec(2)

							if (ex==0.d0) then
								theta = pi/2.d0 * dsign(1.d0,ez)
							else
								theta = datan2(ez,ex)
							endif

							if(theta<0.d0) then
								theta = theta + 2.d0*pi
							endif

							psi_diff_temp(m_bound) = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
																	ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )
							
						else

							psi_diff_temp(m_bound) = xvec(1)**2*omega_temp*mass/eV

						endif
					endif

					if(bc_type==5) then

						call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
										dummy(3),zone,dummy(4),dummy(5),dummy(6))
						! NOTE: this works only because the inner boundary is a circle,
						! in the general case one should use the coordinates of the 
						! point on the boundary

					endif

					continue

				endif

				continue

			endif

			if(zero_dist) then
			! take care of the boundary point
			! (what follows in coord(m,1:2) does not matter, it's just to avoid "funny" numbers)

				ind(m_bound,1) = i
				ind(m_bound,2) = j

				coord(m_bound,1) = rmajor
				coord(m_bound,2) = 0.d0
				coord(m_bound,3) = 2.d0

			endif

       end do
    end do

	allocate(ind_bound(m_bound,4))
	if((bc_type==7).or.(bc_type==8)) then
		allocate(coord_bound(m_bound,6))
	else
		allocate(coord_bound(m_bound,5))
	endif
	if((bc_type==23).or.(bc_type==24)) then
		allocate(psi_diff_bound(m_bound))
	endif

	do i=1,m_bound

		do j=1,4
			ind_bound(i,j) = ind(i,j)
		enddo

		coord_bound(i,1) = 2.d0/dx_a(ind(i,3))*(coord(i,1) - x_coord(ind(i,3))) - 1.d0
		coord_bound(i,2) = 2.d0/dz_a(ind(i,4))*(coord(i,2) - z_coord(ind(i,4))) - 1.d0

		coord_bound(i,3) = coord(i,3) ! "in" or "out"

		coord_bound(i,4) = coord(i,4) ! distance from point to surface

		coord_bound(i,5) = coord(i,5) ! distance from surface to interpolated point

		if((bc_type==7).or.(bc_type==8)) 	coord_bound(i,6) = coord(i,6) ! psi on the boundary

		if((bc_type==23).or.(bc_type==24)) psi_diff_bound(i) = psi_diff_temp(i)

		continue

	enddo

	continue

end subroutine initialize_bc_equations

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc_minimum(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, x_coord, z_coord, dx_a, dz_a, bc_switch

	integer :: n

	integer, dimension(n*n,4) :: ind
	real(kind=dkind), dimension(n*n,6) :: coord
	! third number is distance ratio

	real(kind=dkind) :: ex,ez, dist0,dist
	real(kind=dkind) :: x,z, RP, zP, R0, xb, zb, th
	real(kind=dkind) ::  xs, zs ! points on the surface
	real(kind=dkind), dimension(1:2) :: xvec, fvec
	! solution points and function values from IMSL
	real(kind=dkind), dimension(1:2) :: xguess
	real(kind=dkind) :: fnorm
	external  eq_bound, eq_bound2
	integer :: ntrial = 5000
	real(kind=dkind) :: tolx=1.d-12
	real(kind=dkind) :: tolf = 1.d-9
	integer :: zone
	real(kind=dkind) :: dummy(1:7)
	real(kind=dkind) :: r
	real(kind=dkind), dimension(1:3) :: theta_start
	real(kind=dkind) :: thetamin

	integer :: i,j

	logical :: bound, inner, truebound, zero_dist
	logical :: newton_check

	if((bc_type==3).or.(bc_type==23).or.(((bc_type==4).or.(bc_type==5).or.(bc_type==7).or.  &
		(bc_type==8).or.(bc_type==14).or.(bc_type==24).or.(bc_type==34).or.(bc_type==44).or.(bc_type==54).or.(bc_type==64))  &
		.and.(n>=bc_switch))) then

		continue

	else

		return

	endif

	if(allocated(ind_bound)) deallocate(ind_bound)
	if(allocated(coord_bound)) deallocate(coord_bound)

	m_bound = 0

	coord = 0.d0

	do j=1, n
       do i=1, n

			zero_dist = .false.

			call check_position(i,j,bound,truebound,n)

			if(truebound) then
			! real boundary point

				m_bound = m_bound+1

				call radius(i,j,n,n,ex,ez,r,dx,dz)

				if(ex**2+ez**2<r**2) then
				! the point is very close to the boundary

					R_P = x_coord(i)
					z_P = z_coord(j)

					call radius_1_3(R_P,z_P,ex,ez,theta_start(2),r,zone,dummy(1),dummy(2),dummy(3))

					theta_start(1) = theta_start(2) - pi/5.d0
					theta_start(3) = theta_start(2) + pi/5.d0

					dist = brent(theta_start(1), theta_start(2), theta_start(3), dist_fun, tolf, thetamin)

					call radius_theta(thetamin,r,xvec(1),xvec(2))

					if((bc_type==7).or.(bc_type==8)) then

						xb = xvec(1)
						zb = xvec(2)
						th = atan2(zb,xb-rmajor)
						if(th<0.d0) th = th + 2.d0*pi

						coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
										ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

					endif

					dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

					if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

					else

						! proceed as usual

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

						ind(m_bound,1) = i
						ind(m_bound,2) = j

						inner = .false.

						coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
						coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

						call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
									ind(m_bound,3:4) )

						call check_points(ind(m_bound,3),ind(m_bound,4),inner)

						if(inner) then

							coord(m_bound,3) = 1.d0

						else

							coord(m_bound,3) = -1.d0

						endif

						coord(m_bound,4) = 0.d0 !distance from point to surface
						coord(m_bound,5) = dist !distance from surface to interpolated point


						if(bc_type==5) then

							call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
											dummy(3),zone,dummy(4),dummy(5),dummy(6))
							! NOTE: this works only because the inner boundary is a circle,
							! in the general case one should use the coordinates of the
							! point on the boundary

						endif

					endif

				else
				! the point REALLY is on the boundary (numerically),
				! this is extremely unlikely to occur and dealt with easily

					zero_dist = .true.

					continue

				endif

			elseif(bound) then
				! (external) boundary point

				m_bound = m_bound+1

				R_P = x_coord(i)
				z_P = z_coord(j)

				call radius_1_3(R_P,z_P,ex,ez,theta_start(2),r,zone,dummy(1),dummy(2),dummy(3))

				theta_start(1) = theta_start(2) - pi/5.d0
				theta_start(3) = theta_start(2) + pi/5.d0

				dist = brent(theta_start(1), theta_start(2), theta_start(3), dist_fun, tolf, thetamin)

				call radius_theta(thetamin,r,xvec(1),xvec(2))

				if((bc_type==7).or.(bc_type==8)) then

					xb = xvec(1)
					zb = xvec(2)
					th = atan2(zb,xb-rmajor)
					if(th<0.d0) th = th + 2.d0*pi

					coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
									ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

				endif

				dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

				if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

				else

					if(dist0<0.2d0*sqrt(dx_a(i)**2+dz_a(j)**2)) then

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

					else

						dist = dist0

					endif

					ind(m_bound,1) = i
					ind(m_bound,2) = j

					inner = .false.

					coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
					coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

					call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
								ind(m_bound,3:4) )

					call check_points(ind(m_bound,3),ind(m_bound,4),inner)

					if(inner) then

						coord(m_bound,3) = 1.d0

					else

						coord(m_bound,3) = -1.d0

					endif

					coord(m_bound,4) = dist0 !distance from point to surface
					coord(m_bound,5) = dist !distance from surface to interpolated point

					if(bc_type==5) then

						call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
										dummy(3),zone,dummy(4),dummy(5),dummy(6))
						! NOTE: this works only because the inner boundary is a circle,
						! in the general case one should use the coordinates of the
						! point on the boundary

					endif

					continue

				endif

				continue

			endif

			if(zero_dist) then
			! take care of the boundary point
			! (what follows in coord(m,1:2) does not matter, it's just to avoid "funny" numbers)

				ind(m_bound,1) = i
				ind(m_bound,2) = j

				coord(m_bound,1) = rmajor
				coord(m_bound,2) = 0.d0
				coord(m_bound,3) = 2.d0

			endif

       end do
    end do

	allocate(ind_bound(m_bound,4))
	if((bc_type==7).or.(bc_type==8)) then
		allocate(coord_bound(m_bound,6))
	else
		allocate(coord_bound(m_bound,5))
	endif

	do i=1,m_bound

		do j=1,4
			ind_bound(i,j) = ind(i,j)
		enddo

		coord_bound(i,1) = 2.d0/dx_a(ind(i,3))*(coord(i,1) - x_coord(ind(i,3))) - 1.d0

		coord_bound(i,2) = 2.d0/dz_a(ind(i,4))*(coord(i,2) - z_coord(ind(i,4))) - 1.d0

		coord_bound(i,3) = coord(i,3) ! "in" or "out"

		coord_bound(i,4) = coord(i,4) ! distance from point to surface

		coord_bound(i,5) = coord(i,5) ! distance from surface to interpolated point

		if((bc_type==7).or.(bc_type==8)) 	coord_bound(i,6) = coord(i,6) ! psi on the boundary

		continue

	enddo

	continue

end subroutine initialize_bc_minimum

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc_tri_equations(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this initializes bcs when using triangles instead of squares for interpolating Q
! (should be more stable)
! Sept 7 2013:
! copied over from the (to date unfinished) FLOW routine
! to avoid messy nested ifs, everything is saved anyway (in bound_temp),
! just not transfered to global arrays if not needed
! once and for all, "P" is the grid point where Bc's need to be assigned, "Q" is the interpolated point.
! "Q" is in triangle, that would be "Q2" in FLOW_MARS

! the array variables:

! coord_bound:
! 1 -> RQ (physical space coordinate!)
! 2 -> ZQ (physical space coordinate!)
! 3 -> +-1 (inner/outer) or 2 (point on boundary, it never happens)
! 4 -> distance P - boundary
! 5 -> distance Q - boundary
! 6 -> triangle area (for convenience, with sign)

! ind_bound:
! 1 -> iP
! 2 -> jP

! bound_tri(m_bound,1:3,1:2): i,j coordinates of points 1-3 in the triangle enclosing Q2
! triangle does not use point P (for stability)


	use constant, only : dx, dz, x_coord, z_coord, dx_a, dz_a, bc_switch

	integer :: n

	integer, dimension(n*n,4) :: ind
	real(kind=dkind), dimension(n*n,6) :: coord
	! third number is distance ratio
	real(kind=dkind), dimension(n*n) :: psi_diff_temp
	integer, dimension(n*n,1:3,1:2) :: bound_tri_temp

	real(kind=dkind) :: ex,ez, dist0,dist
	real(kind=dkind) :: x,z, RP, zP, R0, xb, zb, th
	real(kind=dkind) ::  xs, zs ! points on the surface
	real(kind=dkind), dimension(1:2) :: xvec, fvec 
	! solution points and function values from IMSL
	real(kind=dkind), dimension(1:2) :: xguess
	real(kind=dkind) :: fnorm
	external  eq_bound, eq_bound2
	integer :: ntrial = 5000
	real(kind=dkind) :: tolx=1.d-12
	real(kind=dkind) :: tolf = 1.d-9
	integer :: zone
	real(kind=dkind) :: dummy(1:7)
	real(kind=dkind) :: r
	real(kind=dkind) :: omega_temp

	integer :: i,j,k

	logical :: bound, inner, truebound, zero_dist
	logical :: newton_check

	if((bc_type==103).or.(bc_type==113).or.(((bc_type==104).or.(bc_type==114)).and.(n>=bc_switch))) then

		continue

	else

		return

	endif

	if(allocated(ind_bound)) deallocate(ind_bound)
	if(allocated(coord_bound)) deallocate(coord_bound)
	if(allocated(psi_diff_bound)) deallocate(psi_diff_bound)
	if(allocated(bound_tri)) deallocate(bound_tri)

	m_bound = 0

	coord = 0.d0
	
	omega_temp = omega_0

	do j=1, n
       do i=1, n

			zero_dist = .false.

			call check_position(i,j,bound,truebound,n)

			if(truebound) then
			! real boundary point

!				print*, 'initialize_bc: found truebound point'
!				print*, 'i = ', i, ' j = ', j

				m_bound = m_bound+1

				call radius(i,j,n,n,ex,ez,r,dx,dz)

				! the point is very close to the boundary

				R_P = x_coord(i)
				z_P = z_coord(j)

				xguess(1) = R_P
				xguess(2) = z_P

				xvec(1) = xguess(1)
				xvec(2) = xguess(2)

				call newt(ntrial,xvec,2,newton_check,tolx,tolf)

				if(numerical_psi_diff) then

					ex = xvec(1) - rmajor
					ez = xvec(2)

					if (ex==0.d0) then
						theta = pi/2.d0 * dsign(1.d0,ez)
					else
						theta = datan2(ez,ex)
					endif

					if(theta<0.d0) then
						theta = theta + 2.d0*pi
					endif

					psi_diff_temp(m_bound) = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
															ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

				else

					psi_diff_temp(m_bound) = xvec(1)**2*omega_temp*mass/eV

				endif

				! regardless of the distance, just treat the point as a true bondary point
				zero_dist = .true.

			elseif(bound) then
				! (external) boundary point

				m_bound = m_bound+1

				R_P = x_coord(i)
				z_P = z_coord(j)

				xguess(1) = R_P
				xguess(2) = z_P

				xvec(1) = xguess(1)
				xvec(2) = xguess(2)

				call newt(ntrial,xvec,2,newton_check,tolx,tolf)

!!$				if((bc_type==7).or.(bc_type==8)) then
!!$
!!$					xb = xvec(1)
!!$					zb = xvec(2)
!!$					th = atan2(zb,xb-rmajor)
!!$					if(th<0.d0) th = th + 2.d0*pi
!!$
!!$					coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
!!$									ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )
!!$
!!$				endif

				dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

				if(dist0==0.d0) then

					! somehow we would not want this to happen,
					! but we'll treat the point as "real boundary"

					zero_dist = .true.

				else

!					print*, m_bound, i, j
					call get_Q_tri(dist0,dist,i,j,bound_tri_temp(m_bound,1,:),  &
							bound_tri_temp(m_bound,2,:),bound_tri_temp(m_bound,3,:),xvec,n,n,  &
							coord(m_bound,6))

					ind(m_bound,1) = i
					ind(m_bound,2) = j

					inner = .false.

					coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
					coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

					call check_points_tri(bound_tri_temp(m_bound,1,:),  &
							bound_tri_temp(m_bound,2,:),bound_tri_temp(m_bound,3,:),inner)

					if(inner) then

						coord(m_bound,3) = 1.d0

					else

						coord(m_bound,3) = -1.d0

					endif

					coord(m_bound,4) = dist0 !distance from point to surface
					coord(m_bound,5) = dist !distance from surface to interpolated point


					if(numerical_psi_diff) then

						ex = xvec(1) - rmajor
						ez = xvec(2)

						if (ex==0.d0) then
							theta = pi/2.d0 * dsign(1.d0,ez)
						else
							theta = datan2(ez,ex)
						endif

						if(theta<0.d0) then
							theta = theta + 2.d0*pi
						endif

						psi_diff_temp(m_bound) = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
																ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

					else

						psi_diff_temp(m_bound) = xvec(1)**2*omega_temp*mass/eV

					endif

!!$					if(bc_type==5) then
!!$
!!$						call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
!!$										dummy(3),zone,dummy(4),dummy(5),dummy(6))
!!$						! NOTE: this works only because the inner boundary is a circle,
!!$						! in the general case one should use the coordinates of the 
!!$						! point on the boundary
!!$
!!$					endif

					continue

				endif

				continue

			endif

			if(zero_dist) then
			! take care of the boundary point
			! (what follows in coord(m,1:2) does not matter, it's just to avoid "funny" numbers)

				ind(m_bound,1) = i
				ind(m_bound,2) = j

				coord(m_bound,1) = rmajor
				coord(m_bound,2) = 0.d0
				coord(m_bound,3) = 2.d0

			endif

       end do
    end do

	allocate(ind_bound(m_bound,2))
	allocate(coord_bound(m_bound,6))
	allocate(psi_diff_bound(m_bound))
	allocate(bound_tri(m_bound,1:3,1:2))

	psi_diff_bound = 0.d0

	do i=1,m_bound

		do j=1,2
			ind_bound(i,j) = ind(i,j)
		enddo

		do j = 1, 6

			! differently from other options, we store the physical coordinates of Q in coord_bound
			coord_bound(i,j) = coord(i,j)

		enddo

		psi_diff_bound(i) = psi_diff_temp(i)

		do j = 1, 3
		do k = 1, 2
			bound_tri(i,j,k) = bound_tri_temp(i,j,k)
		enddo
		enddo

		continue

	enddo

	continue

end subroutine initialize_bc_tri_equations


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dist_fun(theta) result(dist)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: theta, dist
	real(kind=dkind) :: r, Rloc, Zloc

	call radius_theta(theta,r,Rloc,Zloc)

	dist = sqrt((Rloc-R_P)**2+(Zloc-z_P)**2)

end function dist_fun

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho0(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho

	if (bc_type==1) then

		call bc_psi_rho1(psi,rho,nx,nz)

	elseif (bc_type==3) then

		call bc_psi_rho3(psi,rho,nx,nz)

	elseif (bc_type==4) then

		if(nx<bc_switch) call bc_psi_rho1(psi,rho,nx,nz)
		if(nx>=bc_switch) call bc_psi_rho3(psi,rho,nx,nz)

	elseif ((bc_type==7).or.(bc_type==8).or.(bc_type==17)) then

		call bc_psi_rho7(psi,rho,nx,nz)

	else

		print*, 'unknown option for bc_type:'
		print*, 'bc_type =     ',bc_type

	endif

	if(tri_type==10) call bc_psi_rho1_bis(psi,nx,nz)

	continue

end subroutine bc_psi_rho0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho1(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : x_coord, z_coord

	integer :: nx,nz
!	integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k
	real(kind=dkind) :: x,z,dummy(1:3),alpha
	integer :: zone


    den = dofpsi(0.0d0)

    do j=1, nz
       do i=1, nx

			  if(sort_grid(i,j)<=0) then
				 rho(i,j) = den
				 psi(i,j) = 0.0d0
			  end if

       end do
    end do

	continue

end subroutine bc_psi_rho1


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho3(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4)

		psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

		rho(iQ,jQ) = den

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho3

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho7(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP
	real(kind=dkind) :: th
	integer :: itemax = 10000
	integer :: imin, imax, istep

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	! This isn't super computationally efficient since we should already know the boundary,
	! but this only gets called rarely, and it works just fine
	
	! Commented out if statement below (NEW May 2025)
!	if(nx<bc_switch) then
	! no interpolation

		do j=1, nz
		   do i=1, nx

			  if(sort_grid(i,j)<=0) then

				th = atan2(z_coord(j),x_coord(i)-rmajor)
				if(th<0.d0) th = th + 2.d0*pi

				 rho(i,j) = den
				 psi(i,j) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
							ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )
				 !note that this is not accurate, but it should do for the coarse grids

				 continue

			  end if

		   end do
		end do

		continue
		return

	! Shouldn't need any of the below code, so commented out by Ian (NEW May 2025)
	!	endif

	
	! this would be "else"
	! first pass

	! do i=1,m_bound

	! 	psi_val = coord_bound(i,6)

	! 	xQ = coord_bound(i,1)
	! 	zQ = coord_bound(i,2)

	! 	iQ = ind_bound(i,1)
	! 	jQ = ind_bound(i,2)
	! 	iloc = ind_bound(i,3)
	! 	jloc = ind_bound(i,4)

	! 	fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
	! 	fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
	! 	fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
	! 	fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

	! 			!--------------!
	! 	if( (iloc==0).or.(jloc==0) ) then
	! 		psiloc(1) = psi_val
	! 	else
	! 		psiloc(1) = psi(iloc,jloc)
	! 	endif
	! 			!--------------!
	! 	if( (iloc==nx).or.(jloc==0) ) then
	! 		psiloc(2) = psi_val
	! 	else
	! 		psiloc(2) = psi(iloc+1,jloc)
	! 	endif
	! 			!--------------!
	! 	if( (iloc==nx).or.(jloc==nz) ) then
	! 		psiloc(3) = psi_val
	! 	else
	! 		psiloc(3) = psi(iloc+1,jloc+1)
	! 	endif
	! 			!--------------!
	! 	if( (iloc==0).or.(jloc==nz) ) then
	! 		psiloc(4) = psi_val
	! 	else
	! 		psiloc(4) = psi(iloc,jloc+1)
	! 	endif
	! 			!--------------!

	! 	psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
	! 				psiloc(3)*fi(3) + psiloc(4)*fi(4)

	! 	psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

	! 	rho(iQ,jQ) = den

	! 	if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/2.5d0 ) then
	! 		psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/2.5d0)	! 0.d0 !
	! 	endif

	! 	continue

	! enddo

	! ! second pass
	! error = 0.d0

	! k = 0

	! imin = 1
	! imax = m_bound
	! istep = 1

	! do

	! 	k = k+1

	! 	do i=imin, imax, istep

	! 		if(coord_bound(i,3)==1.d0) cycle !(this is an "internal" point)

	! 		psi_val = coord_bound(i,6)

	! 		xQ = coord_bound(i,1)
	! 		zQ = coord_bound(i,2)

	! 		iQ = ind_bound(i,1)
	! 		jQ = ind_bound(i,2)
	! 		iloc = ind_bound(i,3)
	! 		jloc = ind_bound(i,4)

	! 		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
	! 		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
	! 		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
	! 		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

	! 		psiold = psi(iQ,jQ)

	! 				!--------------!
	! 		if( (iloc==0).or.(jloc==0) ) then
	! 			psiloc(1) = psi_val
	! 		else
	! 			psiloc(1) = psi(iloc,jloc)
	! 		endif
	! 				!--------------!
	! 		if( (iloc==nx).or.(jloc==0) ) then
	! 			psiloc(2) = psi_val
	! 		else
	! 			psiloc(2) = psi(iloc+1,jloc)
	! 		endif
	! 				!--------------!
	! 		if( (iloc==nx).or.(jloc==nz) ) then
	! 			psiloc(3) = psi_val
	! 		else
	! 			psiloc(3) = psi(iloc+1,jloc+1)
	! 		endif
	! 				!--------------!
	! 		if( (iloc==0).or.(jloc==nz) ) then
	! 			psiloc(4) = psi_val
	! 		else
	! 			psiloc(4) = psi(iloc,jloc+1)
	! 		endif
	! 				!--------------!

	! 		psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
	! 					psiloc(3)*fi(3) + psiloc(4)*fi(4)

	! 		psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

	! 		if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/2.5d0 ) then
	! 				psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/2.5d0)	! 0.d0 !
	! 		endif

	! 		error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

	! 		continue

	! 	enddo

	! 	if(error<1.d-13) exit
	! 	if( (k>itemax).and.(itemax<1000001) ) then
	! 		print*, 'bc not converged, error:', error
	! 		exit
	! 	endif

	! 	if (imin == 1) then

	! 		imin = m_bound
	! 		imax = 1
	! 		istep = -1

	! 	else

	! 		imin = 1
	! 		imax = m_bound
	! 		istep = 1

	! 	endif

	! 	error = 0.d0

	! enddo

	continue

	return

end subroutine bc_psi_rho7

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho1_bis(psi,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi
    real (kind=dkind) :: ex,ez
    integer :: i,j

	psi_ext = psic*(1.d0-0.5d0*x_size/a_elps)
	psi_ext = -.1d0

    do j=1, nz
       do i=1, nx

		  call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)

		  if((ex*ex + ez*ez) >= radius_ext**2) then

             psi(i,j) = psi_ext

          end if

       end do
    end do

end subroutine bc_psi_rho1_bis

!-------------------------------------------------------------------
!							two-fluid BC
!-------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Feb 18 2013 added option that my be transported into FLOW
! i_opt controls which variable is updated:
! i_opt = 0 -> update verything
! i_opt = 1 -> update psi
! i_opt = 2 -> update big_Psi
! i_opt = 3 -> update density
! Other note: density CANNOT be set to D(psi) out of the plasma, since that would cause large gradients at the edge!
! (this must also be transported into FLOW)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt

	if (bc_type==1) then

		call bc_TF_1(psi,big_Psi,psi_diff,n_den,nx,nz)

	elseif (bc_type==3) then

		call bc_TF_3(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==4) then

		if(nx<bc_switch) call bc_TF_1(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_3(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==7) then

		call bc_TF_7(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==11) then

		call bc_TF_1b(psi,big_Psi,psi_diff,n_den,nx,nz)

	elseif (bc_type==14) then

		if(nx<bc_switch) call bc_TF_1b(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_3(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==21) then

		call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)

	elseif (bc_type==24) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_3(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==34) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_1_5(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==44) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_1_5_zero(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==54) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_1_5_den(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==64) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_1_5_den_zero(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

!!$	elseif (bc_type==74) then
!!$
!!$		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
!!$		if(nx>=bc_switch) call bc_TF_1_5_big_Psi_zero(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==103) then

		if(nx>=bc_switch) call bc_TF_tri_1(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==104) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_tri_1(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==113) then

		if(nx>=bc_switch) call bc_TF_tri_1(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	elseif (bc_type==114) then

		if(nx<bc_switch) call bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
		if(nx>=bc_switch) call bc_TF_tri_1(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)

	else

		print*, 'unknown option for bc_type:'
		print*, 'bc_type =     ',bc_type
		pause

	endif

	continue

end subroutine bc_TF_0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_1(psi,big_Psi,psi_diff,n_den,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : x_coord, z_coord

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi,  &
    				psi_diff, n_den
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: den
    integer :: i,j

    den = d_TF_ofpsi(0.0d0)

    do j=1, nz
       do i=1, nx

			  if(sort_grid(i,j)<=0) then
				 n_den(i,j) = den
				 psi(i,j) = 0.0d0
				 big_Psi(i,j) = 0.0d0
				 psi_diff(i,j) = 0.d0
			  end if

       end do
    end do

	continue

end subroutine bc_TF_1


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_1b(psi,big_Psi,psi_diff,n_den,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : x_coord, z_coord

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi,  &
    				psi_diff, n_den
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: den
    real(kind=dkind) :: omega0, x2
    integer :: i,j

    den = d_TF_ofpsi(0.0d0)
    omega0 = omegaofpsi(0.d0)

    do j=1, nz
       do i=1, nx

				x2 = x_coord(i)**2

			  if(sort_grid(i,j)<=0) then
				 n_den(i,j) = den
				 psi(i,j) = 0.0d0
				 psi_diff(i,j) = mass/eV * x2 * omega0
				 big_Psi(i,j) = psi(i,j) + psi_diff(i,j)
			  end if

       end do
    end do

	continue

end subroutine bc_TF_1b


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_3(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! to minimize writing, this routine treates two different cases.
! OLD APPROACH:
! set density to D(0) and psi_diff to zero on the boundary
! NEW APPROACH:
! set density to D(0) and psi_diff to an assigned value (theta dependent) on the boundary

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psi_diffP, n_denP, psi_diffQ, n_denQ
	integer :: zone
	integer :: itemax = 10000
	real(kind=dkind) :: diff_max, jump_psi, jump_diff
	real(kind=dkind) :: grad_ratio

	grad_ratio = grad_ratio_0 * 2.d0 / nx

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0
	psi_diff_val = 0.d0

	diff_max = abs(big_Psic - psic)
	jump_psi = (psic-psi_val) * grad_ratio
	jump_diff = (diff_max-psi_diff_val) * grad_ratio

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
			n_den(iQ,jQ) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		! these are needed for interpolation of ANY variable
		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(i_opt<=1) then
		! update psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
			endif

			if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
				psi(iQ,jQ) = psi_val
			endif

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

			if((bc_type==23).or.(bc_type==24)) then
				psi_diff_val = psi_diff_bound(i)
			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = 0.d0
			else
				psiloc(1) = psi_diff(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = 0.d0
			else
				psiloc(2) = psi_diff(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = 0.d0
			else
				psiloc(3) = psi_diff(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = 0.d0
			else
				psiloc(4) = psi_diff(iloc,jloc+1)
			endif
					!--------------!

			! we set psi_diff = 0 or assigned value
			psi_diffP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)
			psi_diffQ = psi_diffP + (psi_diff_val - psi_diffP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi_diffQ -psi_diffP)>jump_diff) then
!				psi_diffQ = psi_diffP - jump_diff
				psi_diffQ = psi_diffP - sign(jump_diff,psi_diffP-psi_diffQ)
			endif

			if(abs(psi_diffQ-psi_diffP)<abs(psi_diffP-psi_diff_val)) then
				psi_diffQ = psi_diff_val
			endif

			big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diffQ

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = n_den(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = n_den(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = n_den(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = n_den(iloc,jloc+1)
			endif
					!--------------!

			n_denP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)
			n_denQ = n_denP + (den - n_denP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
			if(n_denQ<den*1.d-2) n_denQ = den*1.d-2

			n_den(iQ,jQ) = n_denQ

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
			endif

			if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
				psi(iQ,jQ) = psi_val
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

				if((bc_type==23).or.(bc_type==24)) then
					psi_diff_val = psi_diff_bound(i)
				endif

				psiold = big_Psi(iQ,jQ)

						!--------------!
				if( (iloc==0).or.(jloc==0) ) then
					psiloc(1) = 0.d0
				else
					psiloc(1) = psi_diff(iloc,jloc)
				endif
						!--------------!
				if( (iloc==nx).or.(jloc==0) ) then
					psiloc(2) = 0.d0
				else
					psiloc(2) = psi_diff(iloc+1,jloc)
				endif
						!--------------!
				if( (iloc==nx).or.(jloc==nz) ) then
					psiloc(3) = 0.d0
				else
					psiloc(3) = psi_diff(iloc+1,jloc+1)
				endif
						!--------------!
				if( (iloc==0).or.(jloc==nz) ) then
					psiloc(4) = 0.d0
				else
					psiloc(4) = psi_diff(iloc,jloc+1)
				endif
						!--------------!

				psi_diffP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
							psiloc(3)*fi(3) + psiloc(4)*fi(4)

				! we set psi_diff = 0 or assigned value
				psi_diffP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
							psiloc(3)*fi(3) + psiloc(4)*fi(4)
				psi_diffQ = psi_diffP + (psi_diff_val - psi_diffP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

				if(abs(psi_diffQ -psi_diffP)>jump_diff) then
	!				psi_diffQ = psi_diffP - jump_diff
					psi_diffQ = psi_diffP - sign(jump_diff,psi_diffP-psi_diffQ)
				endif

				if(abs(psi_diffQ-psi_diffP)<abs(psi_diffP-psi_diff_val)) then
					psi_diffQ = psi_diff_val
				endif

				big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diffQ

				error = max( error,abs( (big_Psi(iQ,jQ)-psiold)/psiold ) )

			endif

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_3


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_7(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j
	real(kind=dkind) :: x,ex,ez,psi_val
	real(kind=dkind) :: th

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	! no interpolation is needed with these boundary conditions

	do j=1, nz
	   do i=1, nx

		  if(sort_grid(i,j)<=0) then

			th = atan2(z_coord(j),x_coord(i)-rmajor)
			if(th<0.d0) th = th + 2.d0*pi

			if(i_opt<=1) then
			! update psi

				 psi(i,j) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
							ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

			endif


			!-------------------------------------
			if((i_opt==0).or.(i_opt==2)) then
			! update big_Psi

				 big_Psi(i,j) = dbsval(th, big_Psib_ord, big_Psib_data(1:ibreak_bPb,3),  &
							ibreak_bPb-big_Psib_ord, big_Psib_cscoef(1:ibreak_bPb-big_Psib_ord) )

			endif

			!-------------------------------------
			if((i_opt==0).or.(i_opt==3)) then
			! update n_den

				n_den(i,j) = den

				! Don't know what was going on below, but commenting out for now
				! if(i==1) then

				! 	n_den(i,j) = n_den(i+1,j)

				! elseif(i==nx) then

				! 	n_den(i,j) = n_den(i-1,j)

				! elseif(j==1) then

				! 	n_den(i,j) = n_den(i,j+1)

				! elseif(j==nz) then

				! 	n_den(i,j) = n_den(i,j-1)

				! endif

			endif

			if(i_opt/=3) then
				psi_diff(i,j) = big_Psi(i,j) - psi(i,j)
			endif

			continue

		  end if

	   end do
	end do

	continue
	return

end subroutine bc_TF_7



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_13(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! NEW APPROACH:
! set density and psi_diff gradients to zero on the boundary
! (note that of course this will not set the gradient of big_Psi to zero)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP, psi_diffP, n_denP
	integer :: zone
	integer :: itemax = 10000

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
			n_den(iQ,jQ) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		! these are needed for interpolation of ANY variable
		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(i_opt<=1) then
		! update psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = 0.d0
			else
				psiloc(1) = psi_diff(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = 0.d0
			else
				psiloc(2) = psi_diff(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = 0.d0
			else
				psiloc(3) = psi_diff(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = 0.d0
			else
				psiloc(4) = psi_diff(iloc,jloc+1)
			endif
					!--------------!

			psi_diffP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			! we set grad psi_diff = 0, so psi_diffQ = psi_diffP
			big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diffP

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = n_den(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = n_den(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = n_den(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = n_den(iloc,jloc+1)
			endif
					!--------------!

			n_denP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			! we set grad n_den = 0, so n_denQ = n_denP
			n_den(iQ,jQ) = n_denP

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			!-------------------------------------
			if(i_opt==0) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

						!--------------!
				if( (iloc==0).or.(jloc==0) ) then
					psiloc(1) = 0.d0
				else
					psiloc(1) = psi_diff(iloc,jloc)
				endif
						!--------------!
				if( (iloc==nx).or.(jloc==0) ) then
					psiloc(2) = 0.d0
				else
					psiloc(2) = psi_diff(iloc+1,jloc)
				endif
						!--------------!
				if( (iloc==nx).or.(jloc==nz) ) then
					psiloc(3) = 0.d0
				else
					psiloc(3) = psi_diff(iloc+1,jloc+1)
				endif
						!--------------!
				if( (iloc==0).or.(jloc==nz) ) then
					psiloc(4) = 0.d0
				else
					psiloc(4) = psi_diff(iloc,jloc+1)
				endif
						!--------------!

				psi_diffP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
							psiloc(3)*fi(3) + psiloc(4)*fi(4)

				! we set grad psi_diff = 0, so psi_diffQ = psi_diffP
				big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diffP

			endif

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_13

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_21(psi,big_Psi,psi_diff,n_den,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this routine assigns delta_psi in each external point based on input (no interpolations)

	use constant, only : x_coord, z_coord

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi,  &
    				psi_diff, n_den
    real (kind=dkind) :: ex, ez, theta
    real (kind=dkind) :: den
    real(kind=dkind) :: omega0, psi_diff_loc, x
    integer :: i,j

    den = d_TF_ofpsi(0.0d0)
    omega0 = omega_0

	do j=1, nz
		do i=1, nx

			x = x_coord(i)

			if(sort_grid(i,j)<=0) then

				if(numerical_psi_diff) then

					ex = x_coord(i)-rmajor
					ez = z_coord(j)

					if (ex==0.d0) then
						theta = pi/2.d0 * dsign(1.d0,ez)
					else
						theta = datan2(ez,ex)
					endif

					if(theta<0.d0) then
						theta = theta + 2.d0*pi
					endif

					psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
										ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

				else

					psi_diff_loc = mass/eV * x**2 * omega0

				endif

				n_den(i,j) = den
				psi(i,j) = 0.0d0
				psi_diff(i,j) = psi_diff_loc
				big_Psi(i,j) = psi(i,j) + psi_diff(i,j)

			endif

		enddo
	enddo

	continue

end subroutine bc_TF_21

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_1_5(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! HYBRID APPROACH:
! set density to D(0) and psi_diff to an assigned value (theta dependent) on the boundary,
! BUT assigning values as in routine 21 for psi_diff

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psi_diffP, n_denP, psi_diffQ,  &
								n_denQ, psi_diff_loc
	integer :: zone
	integer :: itemax = 10000
	real(kind=dkind) :: diff_max, jump_psi, jump_diff
	real(kind=dkind) :: grad_ratio

	grad_ratio = grad_ratio_0 * 2.d0 / nx

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0
	psi_diff_val = 0.d0

	diff_max = abs(big_Psic - psic)
	jump_psi = (psic-psi_val) * grad_ratio
	jump_diff = (diff_max-psi_diff_val) * grad_ratio

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
			n_den(iQ,jQ) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		! these are needed for interpolation of ANY variable
		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(i_opt<=1) then
		! update psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

			x = x_coord(iQ)

			if(numerical_psi_diff) then

				ex = x_coord(iQ)-rmajor
				ez = z_coord(jQ)

				if (ex==0.d0) then
					theta = pi/2.d0 * dsign(1.d0,ez)
				else
					theta = datan2(ez,ex)
				endif

				if(theta<0.d0) then
					theta = theta + 2.d0*pi
				endif

				psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
									ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

			else

				psi_diff_loc = mass/eV * x**2 * omega_0

			endif

			if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
				psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
			endif

			psi_diff(iQ,jQ) = psi_diff_loc
			big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diff(iQ,jQ)

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = n_den(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = n_den(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = n_den(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = n_den(iloc,jloc+1)
			endif
					!--------------!

			n_denP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)
			n_denQ = n_denP + (den - n_denP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
			if(n_denQ<den*1.d-2) n_denQ = den*1.d-2

			n_den(iQ,jQ) = n_denQ

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

				x = x_coord(iQ)

				if(numerical_psi_diff) then

					ex = x_coord(iQ)-rmajor
					ez = z_coord(jQ)

					if (ex==0.d0) then
						theta = pi/2.d0 * dsign(1.d0,ez)
					else
						theta = datan2(ez,ex)
					endif

					if(theta<0.d0) then
						theta = theta + 2.d0*pi
					endif

					psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
										ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

				else

					psi_diff_loc = mass/eV * x**2 * omega_0

				endif

				if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
					psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
				endif

				psi_diff(iQ,jQ) = psi_diff_loc
				big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diff(iQ,jQ)

			endif

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_1_5

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_1_5_zero(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! HYBRID APPROACH:
! set density to D(0) and psi_diff to an assigned value (theta dependent) on the boundary,
! BUT assigning values as in routine 21 for psi_diff
! contrary to the other bc_TF_1_5 routine, big_Psi = psi_diff

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psi_diffP, n_denP, psi_diffQ,  &
								n_denQ, psi_diff_loc
	integer :: zone
	integer :: itemax = 10000
	real(kind=dkind) :: diff_max, jump_psi, jump_diff
	real(kind=dkind) :: grad_ratio

	grad_ratio = grad_ratio_0 * 2.d0 / nx

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0
	psi_diff_val = 0.d0

	diff_max = abs(big_Psic - psic)
	jump_psi = (psic-psi_val) * grad_ratio
	jump_diff = (diff_max-psi_diff_val) * grad_ratio

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
			n_den(iQ,jQ) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		! these are needed for interpolation of ANY variable
		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(i_opt<=1) then
		! update psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

			x = x_coord(iQ)

			if(numerical_psi_diff) then

				ex = x_coord(iQ)-rmajor
				ez = z_coord(jQ)

				if (ex==0.d0) then
					theta = pi/2.d0 * dsign(1.d0,ez)
				else
					theta = datan2(ez,ex)
				endif

				if(theta<0.d0) then
					theta = theta + 2.d0*pi
				endif

				psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
									ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

			else

				psi_diff_loc = mass/eV * x**2 * omega_0

			endif

			if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
				psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
			endif

			psi_diff(iQ,jQ) = psi_diff_loc
			big_Psi(iQ,jQ) = psi_diff(iQ,jQ)

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = n_den(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = n_den(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = n_den(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = n_den(iloc,jloc+1)
			endif
					!--------------!

			n_denP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)
			n_denQ = n_denP + (den - n_denP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
			if(n_denQ<den*1.d-2) n_denQ = den*1.d-2

			n_den(iQ,jQ) = n_denQ

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

				x = x_coord(iQ)

				if(numerical_psi_diff) then

					ex = x_coord(iQ)-rmajor
					ez = z_coord(jQ)

					if (ex==0.d0) then
						theta = pi/2.d0 * dsign(1.d0,ez)
					else
						theta = datan2(ez,ex)
					endif

					if(theta<0.d0) then
						theta = theta + 2.d0*pi
					endif

					psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
										ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

				else

					psi_diff_loc = mass/eV * x**2 * omega_0

				endif

				if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
					psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
				endif

				psi_diff(iQ,jQ) = psi_diff_loc
				big_Psi(iQ,jQ) = psi_diff(iQ,jQ)

			endif

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_1_5_zero

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_1_5_den(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! HYBRID APPROACH:
! set density to D(0) and psi_diff to an assigned value (theta dependent) on the boundary,
! BUT do not interpolate density

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psi_diffP, n_denP, psi_diffQ,  &
								n_denQ, psi_diff_loc
	integer :: zone
	integer :: itemax = 10000
	real(kind=dkind) :: diff_max, jump_psi, jump_diff
	real(kind=dkind) :: grad_ratio

	grad_ratio = grad_ratio_0 * 2.d0 / nx

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0
	psi_diff_val = 0.d0

	diff_max = abs(big_Psic - psic)
	jump_psi = (psic-psi_val) * grad_ratio
	jump_diff = (diff_max-psi_diff_val) * grad_ratio

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
			n_den(iQ,jQ) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		! these are needed for interpolation of ANY variable
		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(i_opt<=1) then
		! update psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

			x = x_coord(iQ)

			if(numerical_psi_diff) then

				ex = x_coord(iQ)-rmajor
				ez = z_coord(jQ)

				if (ex==0.d0) then
					theta = pi/2.d0 * dsign(1.d0,ez)
				else
					theta = datan2(ez,ex)
				endif

				if(theta<0.d0) then
					theta = theta + 2.d0*pi
				endif

				psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
									ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

			else

				psi_diff_loc = mass/eV * x**2 * omega_0

			endif

			if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
				psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
			endif

			psi_diff(iQ,jQ) = psi_diff_loc
			big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diff(iQ,jQ)

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

			n_den(iQ,jQ) = den

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

				x = x_coord(iQ)

				if(numerical_psi_diff) then

					ex = x_coord(iQ)-rmajor
					ez = z_coord(jQ)

					if (ex==0.d0) then
						theta = pi/2.d0 * dsign(1.d0,ez)
					else
						theta = datan2(ez,ex)
					endif

					if(theta<0.d0) then
						theta = theta + 2.d0*pi
					endif

					psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
										ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

				else

					psi_diff_loc = mass/eV * x**2 * omega_0

				endif

				if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
					psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
				endif

				psi_diff(iQ,jQ) = psi_diff_loc
				big_Psi(iQ,jQ) = psi(iQ,jQ) + psi_diff(iQ,jQ)

			endif

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_1_5_den


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_1_5_den_zero(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! HYBRID APPROACH:
! set density to D(0) and psi_diff to an assigned value (theta dependent) on the boundary,
! BUT assigning values as in routine 21 for psi_diff
! contrary to the other bc_TF_1_5 routine, big_Psi = psi_diff
! density is not interpolated

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psi_diffP, n_denP, psi_diffQ,  &
								n_denQ, psi_diff_loc
	integer :: zone
	integer :: itemax = 10000
	real(kind=dkind) :: diff_max, jump_psi, jump_diff
	real(kind=dkind) :: grad_ratio

	grad_ratio = grad_ratio_0 * 2.d0 / nx

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0
	psi_diff_val = 0.d0

	diff_max = abs(big_Psic - psic)
	jump_psi = (psic-psi_val) * grad_ratio
	jump_diff = (diff_max-psi_diff_val) * grad_ratio

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
			n_den(iQ,jQ) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		! these are needed for interpolation of ANY variable
		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(i_opt<=1) then
		! update psi

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

			x = x_coord(iQ)

			if(numerical_psi_diff) then

				ex = x_coord(iQ)-rmajor
				ez = z_coord(jQ)

				if (ex==0.d0) then
					theta = pi/2.d0 * dsign(1.d0,ez)
				else
					theta = datan2(ez,ex)
				endif

				if(theta<0.d0) then
					theta = theta + 2.d0*pi
				endif

				psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
									ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

			else

				psi_diff_loc = mass/eV * x**2 * omega_0

			endif

			if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
				psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
			endif

			psi_diff(iQ,jQ) = psi_diff_loc
			big_Psi(iQ,jQ) = psi_diff(iQ,jQ)

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

			n_den(iQ,jQ) = den

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
					psi(iQ,jQ) = psi_val
				endif
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

				x = x_coord(iQ)

				if(numerical_psi_diff) then

					ex = x_coord(iQ)-rmajor
					ez = z_coord(jQ)

					if (ex==0.d0) then
						theta = pi/2.d0 * dsign(1.d0,ez)
					else
						theta = datan2(ez,ex)
					endif

					if(theta<0.d0) then
						theta = theta + 2.d0*pi
					endif

					psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
										ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )

				else

					psi_diff_loc = mass/eV * x**2 * omega_0

				endif

				if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
					psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
				endif

				psi_diff(iQ,jQ) = psi_diff_loc
				big_Psi(iQ,jQ) = psi_diff(iQ,jQ)

			endif

			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_1_5_den_zero

!!$
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$subroutine bc_TF_1_5_big_Psi_zero(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$! HYBRID APPROACH:
!!$! set density to D(0) and psi_diff to an assigned value (theta dependent) on the boundary
!!$! neither big_Psi nor density are interpolated
!!$
!!$    integer, intent(in) :: nx,nz
!!$    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
!!$	integer :: i_opt
!!$    real (kind=dkind), dimension(1:4) :: fi,psiloc
!!$    real (kind=dkind) :: xQ, zQ
!!$    real (kind=dkind) :: den, den_in, den_out
!!$    integer :: i,j,k, iloc, jloc, iQ, jQ
!!$	real(kind=dkind) :: error,psiold
!!$	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psi_diffP, n_denP, psi_diffQ,  &
!!$								n_denQ, psi_diff_loc
!!$	integer :: zone
!!$	integer :: itemax = 10000
!!$	real(kind=dkind) :: diff_max, jump_psi, jump_diff
!!$	real(kind=dkind) :: grad_ratio
!!$
!!$	grad_ratio = grad_ratio_0 * 2.d0 / nx
!!$
!!$    den = d_TF_ofpsi(0.0d0)
!!$	psi_val = 0.d0
!!$	psi_diff_val = 0.d0
!!$
!!$	diff_max = abs(big_Psic - psic)
!!$	jump_psi = (psic-psi_val) * grad_ratio
!!$	jump_diff = (diff_max-psi_diff_val) * grad_ratio
!!$
!!$	! first pass
!!$
!!$	do i=1,m_bound
!!$
!!$		if(coord_bound(i,3)==2.d0) then
!!$		! the point is really on the boundary
!!$
!!$			iQ = ind_bound(i,1)
!!$			jQ = ind_bound(i,2)
!!$			psi(iQ,jQ) = psi_val
!!$
!!$			big_Psi(iQ,jQ) = psi(iQ,jQ) ! for want of better ideas
!!$			n_den(iQ,jQ) = den ! for want of better ideas
!!$
!!$			! then of course we can skip the rest
!!$			cycle
!!$
!!$		endif
!!$
!!$		xQ = coord_bound(i,1)
!!$		zQ = coord_bound(i,2)
!!$
!!$		iQ = ind_bound(i,1)
!!$		jQ = ind_bound(i,2)
!!$		iloc = ind_bound(i,3)
!!$		jloc = ind_bound(i,4)
!!$
!!$		! these are needed for interpolation of ANY variable
!!$		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
!!$		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
!!$		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
!!$		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)
!!$
!!$		if(i_opt<=1) then
!!$		! update psi
!!$
!!$					!--------------!
!!$			if( (iloc==0).or.(jloc==0) ) then
!!$				psiloc(1) = psi_val
!!$			else
!!$				psiloc(1) = psi(iloc,jloc)
!!$			endif
!!$					!--------------!
!!$			if( (iloc==nx).or.(jloc==0) ) then
!!$				psiloc(2) = psi_val
!!$			else
!!$				psiloc(2) = psi(iloc+1,jloc)
!!$			endif
!!$					!--------------!
!!$			if( (iloc==nx).or.(jloc==nz) ) then
!!$				psiloc(3) = psi_val
!!$			else
!!$				psiloc(3) = psi(iloc+1,jloc+1)
!!$			endif
!!$					!--------------!
!!$			if( (iloc==0).or.(jloc==nz) ) then
!!$				psiloc(4) = psi_val
!!$			else
!!$				psiloc(4) = psi(iloc,jloc+1)
!!$			endif
!!$					!--------------!
!!$
!!$			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
!!$						psiloc(3)*fi(3) + psiloc(4)*fi(4)
!!$
!!$			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
!!$
!!$			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!!$!				psi(iQ,jQ) = psiP - jump_psi
!!$				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
!!$				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
!!$					psi(iQ,jQ) = psi_val
!!$				endif
!!$			endif
!!$
!!$		endif
!!$
!!$		!-------------------------------------
!!$		if((i_opt==0).or.(i_opt==2)) then
!!$		! update big_Psi
!!$
!!$			x = x_coord(iQ)
!!$
!!$			if(numerical_psi_diff) then
!!$
!!$				ex = x_coord(iQ)-rmajor
!!$				ez = z_coord(jQ)
!!$
!!$				if (ex==0.d0) then
!!$					theta = pi/2.d0 * dsign(1.d0,ez)
!!$				else
!!$					theta = datan2(ez,ex)
!!$				endif
!!$
!!$				if(theta<0.d0) then
!!$					theta = theta + 2.d0*pi
!!$				endif
!!$
!!$				psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
!!$									ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )
!!$
!!$			else
!!$
!!$				psi_diff_loc = mass/eV * x**2 * omega_0
!!$
!!$			endif
!!$
!!$			if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
!!$				psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
!!$			endif
!!$
!!$			psi_diff(iQ,jQ) = psi_diff_loc
!!$			big_Psi(iQ,jQ) = psi_diff(iQ,jQ)
!!$
!!$		endif
!!$
!!$		!-------------------------------------
!!$		if((i_opt==0).or.(i_opt==3)) then
!!$		! update n_den
!!$
!!$			n_den(iQ,jQ) = den
!!$
!!$		endif
!!$
!!$		continue
!!$
!!$	enddo
!!$
!!$	if(i_opt>=2) then
!!$		return
!!$		! psi is not updated, thus there is no need to cycle
!!$	endif
!!$
!!$	! second pass
!!$	error = 0.d0
!!$	k = 0
!!$
!!$	do
!!$
!!$		k = k+1
!!$
!!$		do i=1,m_bound
!!$
!!$			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)
!!$
!!$			xQ = coord_bound(i,1)
!!$			zQ = coord_bound(i,2)
!!$
!!$			iQ = ind_bound(i,1)
!!$			jQ = ind_bound(i,2)
!!$			iloc = ind_bound(i,3)
!!$			jloc = ind_bound(i,4)
!!$
!!$			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
!!$			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
!!$			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
!!$			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)
!!$
!!$			psiold = psi(iQ,jQ)
!!$
!!$					!--------------!
!!$			if( (iloc==0).or.(jloc==0) ) then
!!$				psiloc(1) = psi_val
!!$			else
!!$				psiloc(1) = psi(iloc,jloc)
!!$			endif
!!$					!--------------!
!!$			if( (iloc==nx).or.(jloc==0) ) then
!!$				psiloc(2) = psi_val
!!$			else
!!$				psiloc(2) = psi(iloc+1,jloc)
!!$			endif
!!$					!--------------!
!!$			if( (iloc==nx).or.(jloc==nz) ) then
!!$				psiloc(3) = psi_val
!!$			else
!!$				psiloc(3) = psi(iloc+1,jloc+1)
!!$			endif
!!$					!--------------!
!!$			if( (iloc==0).or.(jloc==nz) ) then
!!$				psiloc(4) = psi_val
!!$			else
!!$				psiloc(4) = psi(iloc,jloc+1)
!!$			endif
!!$					!--------------!
!!$
!!$			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
!!$						psiloc(3)*fi(3) + psiloc(4)*fi(4)
!!$
!!$			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
!!$
!!$			if(abs(psi(iQ,jQ)-psiP)>jump_psi) then
!!$!				psi(iQ,jQ) = psiP - jump_psi
!!$				psi(iQ,jQ) = psiP - sign(jump_psi,psiP-psi(iQ,jQ))
!!$				if(abs(psi(iQ,jQ)-psiP)<abs(psiP-psi_val)) then
!!$					psi(iQ,jQ) = psi_val
!!$				endif
!!$			endif
!!$
!!$			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )
!!$
!!$			!-------------------------------------
!!$		if((i_opt==0).or.(i_opt==2)) then
!!$			! also update big_Psi
!!$			! there is no need to update n, since n does not depend on psi
!!$
!!$				x = x_coord(iQ)
!!$
!!$				if(numerical_psi_diff) then
!!$
!!$					ex = x_coord(iQ)-rmajor
!!$					ez = z_coord(jQ)
!!$
!!$					if (ex==0.d0) then
!!$						theta = pi/2.d0 * dsign(1.d0,ez)
!!$					else
!!$						theta = datan2(ez,ex)
!!$					endif
!!$
!!$					if(theta<0.d0) then
!!$						theta = theta + 2.d0*pi
!!$					endif
!!$
!!$					psi_diff_loc = dbsval(theta, psi_diff_ord, edge_psi_diff_data(1:ibreak_psi_diff,3),  &
!!$										ibreak_psi_diff-psi_diff_ord, psi_diff_cscoef(1:ibreak_psi_diff-psi_diff_ord) )
!!$
!!$				else
!!$
!!$					psi_diff_loc = mass/eV * x**2 * omega_0
!!$
!!$				endif
!!$
!!$				if(abs(psi_diff_loc -psi_diff(iloc,jloc))>jump_diff) then
!!$					psi_diff_loc = psi_diff(iloc,jloc) - sign(jump_diff,psi_diff_loc -psi_diff(iloc,jloc))
!!$				endif
!!$
!!$				psi_diff(iQ,jQ) = psi_diff_loc
!!$				big_Psi(iQ,jQ) = psi_diff(iQ,jQ)
!!$
!!$			endif
!!$
!!$			continue
!!$
!!$		enddo
!!$
!!$		if(error<1.d-13) exit
!!$		if( (k>itemax).and.(itemax<1000001) ) then
!!$			print*, 'bc not converged, error:', error
!!$			exit
!!$		endif
!!$
!!$		error = 0.d0
!!$
!!$	enddo
!!$
!!$	continue
!!$
!!$	return
!!$
!!$end subroutine bc_TF_1_5_big_Psi_zero



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_TF_tri_1(psi,big_Psi,psi_diff,n_den,nx,nz,i_opt)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this does the same as bc_TF_3, i.e. set psi to 0, density to D(0)
! and psi_diff to an assigned value (theta dependent) on the boundary,
! but using triangular elements to interpolate Q (inner point)
! by construction, P (point in which BCs are assigned) is never a triangle vertex
! if bc_type==113,114 all points are internal and there is no need to cycle

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, big_Psi, psi_diff, n_den
	integer :: i_opt
    real (kind=dkind), dimension(1:3) :: fi3,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iP, jP
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x, ex, ez, psi_val, psi_diff_val, psiP, psiQ, psi_diffP, n_denP, psi_diffQ, n_denQ
	integer :: zone
	real(kind=dkind) :: diff_max, jump_psi, jump_diff
	real(kind=dkind) :: grad_ratio
	integer :: ind_tri(1:3,1:2)
	integer :: itemax = 10000

	grad_ratio = grad_ratio_0 * 2.d0 / nx

    den = d_TF_ofpsi(0.0d0)
	psi_val = 0.d0
	psi_diff_val = 0.d0

	diff_max = abs(big_Psic - psic)
	jump_psi = (psic-psi_val) * grad_ratio
	jump_diff = (diff_max-psi_diff_val) * grad_ratio

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			! once and for all, BCs are assigned in point "P"
			iP = ind_bound(i,1)
			jP = ind_bound(i,2)
			psi(iP,jP) = psi_val

			big_Psi(iP,jP) = psi(iP,jP) + psi_diff_bound(i) ! for want of better ideas
			n_den(iP,jP) = den ! for want of better ideas

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)	! this is officially "Q" (the interpolated point)
		zQ = coord_bound(i,2)	! and these are the physical coordinates!

		iP = ind_bound(i,1)
		jP = ind_bound(i,2)

		! these are needed for interpolation of ANY variable
		do j = 1, 3
		do k = 1, 2
			ind_tri(j,k) = bound_tri(i,j,k)
		enddo
		enddo

		! the following are global variables...
		xT1 = x_coord(ind_tri(1,1))
		xT2 = x_coord(ind_tri(2,1))
		xT3 = x_coord(ind_tri(3,1))

		zT1 = z_coord(ind_tri(1,2))
		zT2 = z_coord(ind_tri(2,2))
		zT3 = z_coord(ind_tri(3,2))

		two_A_tri = 2.d0 * coord_bound(i,6)
		! ...until here

		fi3(1) = N1_tri(xQ,zQ)
		fi3(2) = N2_tri(xQ,zQ)
		fi3(3) = N3_tri(xQ,zQ)

		if(i_opt<=1) then
		! update psi

			!--------------!
			psiloc(1) = psi(ind_tri(1,1),ind_tri(1,2))
			psiloc(2) = psi(ind_tri(2,1),ind_tri(2,2))
			psiloc(3) = psi(ind_tri(3,1),ind_tri(3,2))
			!--------------!

			psiQ = psiloc(1)*fi3(1) + psiloc(2)*fi3(2) + psiloc(3)*fi3(3)

			psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iP,jP)-psiQ)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iP,jP) = psiQ - sign(jump_psi,psiQ-psi(iP,jP))
			endif

			if(abs(psi(iP,jP)-psiQ)<abs(psiQ-psi_val)) then
				psi(iP,jP) = psi_val
			endif

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
		! update big_Psi

			psi_diff_val = psi_diff_bound(i)

			!--------------!
			psiloc(1) = psi_diff(ind_tri(1,1),ind_tri(1,2))
			psiloc(2) = psi_diff(ind_tri(2,1),ind_tri(2,2))
			psiloc(3) = psi_diff(ind_tri(3,1),ind_tri(3,2))
			!--------------!

			! we set psi_diff = 0 or assigned value
			psi_diffQ = psiloc(1)*fi3(1) + psiloc(2)*fi3(2) + psiloc(3)*fi3(3)
			psi_diffP = psi_diffQ + (psi_diff_val - psi_diffQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi_diffQ -psi_diffP)>jump_diff) then
!				psi_diffQ = psi_diffP - jump_diff
				psi_diffP = psi_diffQ - sign(jump_diff,psi_diffQ-psi_diffP)
			endif

			if(abs(psi_diffQ-psi_diffP)<abs(psi_diffQ-psi_diff_val)) then
				psi_diffP = psi_diff_val
			endif

			big_Psi(iP,jP) = psi(iP,jP) + psi_diffP

		endif

		!-------------------------------------
		if((i_opt==0).or.(i_opt==3)) then
		! update n_den

			!--------------!
			psiloc(1) = n_den(ind_tri(1,1),ind_tri(1,2))
			psiloc(2) = n_den(ind_tri(2,1),ind_tri(2,2))
			psiloc(3) = n_den(ind_tri(3,1),ind_tri(3,2))
			!--------------!

			n_denQ = psiloc(1)*fi3(1) + psiloc(2)*fi3(2) + psiloc(3)*fi3(3)
			n_denP = n_denQ + (den - n_denQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
			if(n_denP<den*5.d-2) n_denP = den*5.d-2

			n_den(iP,jP) = n_denP

		endif

		continue

	enddo

	if(i_opt>=2) then
		return
		! psi is not updated, thus there is no need to cycle
	endif

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)	! this is officially "Q" (the interpolated point)
			zQ = coord_bound(i,2)	! and these are the physical coordinates!

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)

			! these are needed for interpolation of ANY variable
			do j = 1, 3
			do k = 1, 2
				ind_tri(j,k) = bound_tri(i,j,k)
			enddo
			enddo

			! the following are global variables...
			xT1 = x_coord(ind_tri(1,1))
			xT2 = x_coord(ind_tri(2,1))
			xT3 = x_coord(ind_tri(3,1))

			zT1 = z_coord(ind_tri(1,2))
			zT2 = z_coord(ind_tri(2,2))
			zT3 = z_coord(ind_tri(3,2))

			two_A_tri = 2.d0 * coord_bound(i,6)
			! ...until here

			fi3(1) = N1_tri(xQ,zQ)
			fi3(2) = N2_tri(xQ,zQ)
			fi3(3) = N3_tri(xQ,zQ)

			psiold = psi(iP,jP)

			!--------------!
			psiloc(1) = psi(ind_tri(1,1),ind_tri(1,2))
			psiloc(2) = psi(ind_tri(2,1),ind_tri(2,2))
			psiloc(3) = psi(ind_tri(3,1),ind_tri(3,2))
			!--------------!

			psiQ = psiloc(1)*fi3(1) + psiloc(2)*fi3(2) + psiloc(3)*fi3(3)

			psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(abs(psi(iP,jP)-psiQ)>jump_psi) then
!				psi(iQ,jQ) = psiP - jump_psi
				psi(iP,jP) = psiQ - sign(jump_psi,psiQ-psi(iP,jP))
			endif

			if(abs(psi(iP,jP)-psiQ)<abs(psiQ-psi_val)) then
				psi(iP,jP) = psi_val
			endif

			error = max( error,abs( (psi(iP,jP)-psiold)/psiold ) )

			!-------------------------------------
		if((i_opt==0).or.(i_opt==2)) then
			! also update big_Psi
			! there is no need to update n, since n does not depend on psi

				psi_diff_val = psi_diff_bound(i)

				psiold = big_Psi(iP,jP)

				!--------------!
				psiloc(1) = psi_diff(ind_tri(1,1),ind_tri(1,2))
				psiloc(2) = psi_diff(ind_tri(2,1),ind_tri(2,2))
				psiloc(3) = psi_diff(ind_tri(3,1),ind_tri(3,2))
				!--------------!

				! we set psi_diff = 0 or assigned value
				psi_diffQ = psiloc(1)*fi3(1) + psiloc(2)*fi3(2) + psiloc(3)*fi3(3)
				psi_diffP = psi_diffQ + (psi_diff_val - psi_diffQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

				if(abs(psi_diffQ -psi_diffP)>jump_diff) then
!					psi_diffQ = psi_diffP - jump_diff
					psi_diffP = psi_diffQ - sign(jump_diff,psi_diffQ-psi_diffP)
				endif

				if(abs(psi_diffQ-psi_diffP)<abs(psi_diffQ-psi_diff_val)) then
					psi_diffP = psi_diff_val
				endif

				big_Psi(iP,jP) = psi(iP,jP) + psi_diffP

				error = max( error,abs( (big_Psi(iP,jP)-psiold)/psiold ) )

			endif

			! we never updated n, but in reality n is extrapolated, too, so it should be updated
			if((i_opt==0).or.(i_opt==3)) then
			! update n_den

				!--------------!
				psiloc(1) = n_den(ind_tri(1,1),ind_tri(1,2))
				psiloc(2) = n_den(ind_tri(2,1),ind_tri(2,2))
				psiloc(3) = n_den(ind_tri(3,1),ind_tri(3,2))
				!--------------!

				n_denQ = psiloc(1)*fi3(1) + psiloc(2)*fi3(2) + psiloc(3)*fi3(3)
				n_denP = n_denQ + (den - n_denQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
				if(n_denP<den*5.d-2) n_denP = den*5.d-2

				n_den(iP,jP) = n_denP

			endif
			continue

		enddo

		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_TF_tri_1


!--------------------------------------FE routines for interpolations--------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N1_tri(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT2*zT3-xT3*zT2) + (zT2-zT3)*x + (xT3-xT2)*z
	answer = answer / two_A_tri

end function N1_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N2_tri(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT3*zT1-xT1*zT3) + (zT3-zT1)*x + (xT1-xT3)*z
	answer = answer / two_A_tri

end function N2_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N3_tri(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT1*zT2-xT2*zT1) + (zT1-zT2)*x + (xT2-xT1)*z
	answer = answer / two_A_tri

end function N3_tri


!-------------------------------------------------------------------
!							free functions
!-------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer,apsi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if((eq_type>=10).or.(D_TF_num)) then
		answer = mass*d_TF_ofpsi(psi)
		return
	endif

	if(psi==psi_flag) then
		answer = d_loc
		return
	elseif(psi==psic_flag) then
		answer = dc_loc
		return
	endif

    if (zone_loc==-1) then
		apsi = 0.d0
    endif
    
	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif


	if (numerical_n) then

		if ((apsi/psi_max)<=d_data(1,1)) then
			answer = d_data(1,2)
		elseif (apsi>=psi_max) then
			answer = d_data(ibreak_d-d_ord,2)
		else

			answer = dbsval(apsi/psi_max, d_ord, d_data(1:ibreak_d,3),  &
					ibreak_d-d_ord, d_cscoef(1:ibreak_d-d_ord) )

		endif

	else

	    if(psi <= psi_max*fraction ) then
		   answer = dedge
		else if(dabs(psi) > psi_max) then
	       answer = dcenter
		else
		   answer = dedge + (dcenter - dedge)*  &
				dabs( psi/psi_max-fraction )**alpha_rho / ( 1.d0-fraction )**alpha_rho
		end if

	endif


	continue

end function dofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dddpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer,apsi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if((eq_type>=10).or.(D_TF_num)) then
		answer = mass*dd_TF_dpsi(psi)
		return
	endif

	if(psi==psi_flag) then
		answer = dp_loc
		return
	elseif(psi==psic_flag) then
		answer = dcp_loc
		return
	endif

    if (zone_loc==-1) then
		apsi = 0.d0
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (numerical_n) then

			if ((apsi/psi_max)<=d_data(1,1)) then

				answer = dbsder(1,d_data(1,1), d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )

			elseif (apsi>=psi_max) then

				answer = dbsder(1,1.d0-1.d-6, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )

			else

				answer = dbsder(1,apsi/psi_max, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )

			endif

			answer=answer/psi_max

	else

	    if(psi <= psi_max*fraction ) then
		   answer = 0.0d0
	    else if(dabs(psi) > psi_max) then
		   answer = 0.0d0
	    else
		   answer = alpha_rho*(dcenter - dedge)/psi_max  &
					*dabs( psi/psi_max-fraction )**(alpha_rho-1.d0)/( 1.d0-fraction )**alpha_rho
		end if

	endif

	if(abs(psi)>0.d0) answer = answer*psi/abs(psi)

	continue

end function dddpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function pofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer,apsi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag) then
		answer = p_loc
		return
	elseif(psi==psic_flag) then
		answer = pc_loc
		return
	endif

    if (zone_loc==-1) then
		if(numerical_p_iso) then
			apsi=-dabs(psi)
			if(apsi/psi_max<p_iso_data(1,1)) then
				answer = p_iso_data(1,2)
			else
				answer = dbsval(apsi/psi_max, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
					ibreak_p-p_iso_ord, p_iso_cscoef(1:ibreak_p-p_iso_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

if (numerical_p_iso) then

			if ((apsi/psi_max)<=p_iso_data(1,1)) then
				answer = p_iso_data(1,2)
			elseif (apsi>=psi_max) then
				answer = p_iso_data(ibreak_p-p_iso_ord,2)
			else

				answer = dbsval(apsi/psi_max, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
						ibreak_p-p_iso_ord, p_iso_cscoef(1:ibreak_p-p_iso_ord) )

			endif

	else

		if (p_opt == 10) then
		! RFP profile

			answer = pedge + (pcenter-pedge) * (apsi/psic_13)**2 *  &
				 (6.d0 - 8.d0*apsi/psic_13 + 3.d0*(apsi/psic_13)**2)


		else

			if(apsi/psic_13 > 1.0d0) then
				answer = pcenter
			elseif(psi >= (psi_max*fraction) ) then
				answer = pedge + (pcenter-pedge)* (   (apsi/psic_13 - &
												fraction)/ (1.d0-fraction)   )**alpha
			else
				answer = pedge
			end if

		endif



	endif

	if(answer<=0.d0) answer = pedge*1.d-10


	continue

end function pofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dpdpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag) then
		answer = pp_loc
		return
	elseif(psi==psic_flag) then
		answer = pcp_loc
		return
	endif


    if (zone_loc==-1) then

		if(numerical_p_iso) then

			apsi=-dabs(psi)

			if(apsi/psi_max<p_iso_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psi_max, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(:) )
			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

	if (numerical_p_iso) then

			if ((apsi/psi_max)<=p_iso_data(1,1)) then

				answer = dbsder(1,p_iso_data(1,1), p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(:) )

			elseif (apsi>=psi_max) then

				answer = dbsder(1,1.d0-1.d-6, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(:) )

			else

				answer = dbsder(1,apsi/psi_max, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(:) )

			endif

		answer=answer/psi_max

	else

		if (p_opt == 10) then
		! RFP profile

			answer = 12.d0 * (pcenter-pedge) * (1.d0-apsi/psic_13)**2 * apsi/psic_13 /psic_13

		else

			if(dabs(psi) > dabs(psic_13).or.(dabs(psi/psi_max)==0.d0) ) then
				answer = 0.d0
			elseif(psi >= psi_max*fraction ) then
			   answer = alpha*(pcenter-pedge)/psic_13*  &
						( apsi/psic_13-fraction )**(alpha-1.0d0)/ (1.d0-fraction)**alpha
			else
				answer = 0.0d0
			end if

		endif

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

	continue

end function dpdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function bzero(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real(kind=dkind) :: apsi,x,apsim
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep) then
		answer = b0_loc
		return
	elseif(psi==psic_flag) then
		answer = b0c_loc
		return
	endif


    if (zone_loc==-1) then
		if(numerical_F) then
			apsi=-dabs(psi)
			if(apsi/psi_max<F_data(1,1)) then
				answer = F_data(1,2)
			else
				answer = dbsval(apsi/psi_max, F_ord, F_data(1:ibreak_B0,3),  &
					ibreak_B0-F_ord, F_cscoef(1:ibreak_B0-F_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif
    
	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	apsim=min(apsi,psi_max)

	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

	if(numerical_F) then

		if ((apsi/psi_max)<=F_data(1,1)) then
			answer = F_data(1,2)
		elseif((apsi/psi_max)>=1.d0) then
			answer = F_data(ibreak_B0-F_ord,2)
		else

			answer = dbsval(apsi/psi_max, F_ord, F_data(1:ibreak_B0,3),  &
					ibreak_B0-F_ord, F_cscoef(1:ibreak_B0-F_ord) )

		endif

	else

		if (F_opt==0) then
		! no toroidal field, for levitated dipole

			answer = 0.d0

		elseif (F_opt==1) then

			if (apsi/psi_max<fraction) then

				answer = fvacuum / rmajor

			elseif (apsi>psi_max) then

				answer = fcenter / rmajor

			else

				answer = (fvacuum + (fcenter - fvacuum)*  &
							( apsi/psi_max-fraction )**kappa/ (1.d0-fraction)**kappa)/rmajor

			endif

		elseif(F_opt==2) then

			answer = dsqrt(fvacuum**2-eta_P*2.d0*mu_mag*rmajor**2*pofpsi(psi))/rmajor

		elseif(F_opt==5) then
		! RFP profile

			if(abs(psi/psi_max)>1.d0) apsi = psi_max
			! to avoid problems with the powers

			answer = b_phi_zero * ( 1.d0 + mu_RFP * ( apsi/psic_13 - 1.d0 +  &
						(1.d0 - apsi/psic_13)**(kappa+1.d0)/(kappa+1.d0) ) )

		else

			print*, 'wrong option for B_zero:   ',F_opt
			pause
			stop

		endif

	    if(psi < psi_max*fraction ) answer = b_phi_zero

	endif

	continue

end function bzero

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dbzerodpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real(kind=dkind) :: apsi,x,y,apsim
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep) then
		answer = b0p_loc
		return
	elseif(psi==psic_flag) then
		answer = b0cp_loc
		return
	endif

    if (zone_loc==-1) then

		if(numerical_F) then

			apsi=-dabs(psi)

			if(apsi/psi_max<F_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psi_max, F_ord, F_data(1:ibreak_B0,3),  &
						ibreak_B0-F_ord, F_cscoef(:) )
			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif
        
	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	apsim=min(apsi,psi_max)

	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

	if(numerical_F) then

			if (apsi/psi_max<=F_data(1,1)) then

				answer = 0.d0

			elseif (apsi>=psi_max) then

				answer = dbsder(1,1.d0-1.d-6, F_ord, F_data(1:ibreak_B0,3),  &
							ibreak_B0-F_ord, F_cscoef(:) )

			else

				answer = dbsder(1,apsi/psi_max, F_ord, F_data(1:ibreak_B0,3),  &
								ibreak_B0-F_ord, F_cscoef(:) )

			endif

		answer=answer/psi_max

	else

		if (F_opt==0) then
		! no toroidal field, for levitated dipole

			answer = 0.d0

		elseif (F_opt==1) then

		   answer = kappa*(fcenter - fvacuum)/psi_max* &
				 ( apsi/psi_max-fraction )**(kappa-1.0d0)/rmajor/ (1.d0-fraction)**kappa

		elseif(F_opt==2) then

			answer = -mu_mag*rmajor**2*eta_P*dpdpsi(psi)/  &
				dsqrt(fvacuum**2-eta_P*2.d0*mu_mag*rmajor**2*pofpsi(psi))/rmajor

			if(apsi>0.d0) answer = answer*psi/apsi

		elseif(F_opt==5) then
		! RFP profile

			if(abs(psi/psi_max)>1.d0) apsi = psi_max
			! to avoid problems with the powers

			answer = - b_phi_zero*mu_RFP * ( - 1.d0 + (1.d0 - apsi/psic_13)**kappa ) / psic_13

		endif

		if(psi < psi_max*fraction ) answer = 0.d0

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dbzerodpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function mach_theta(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi !,M0
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep2) then
		answer = mth_loc
		return
	endif

    if (zone_loc==-1) then
		if(numerical_mtheta) then
			apsi=-dabs(psi)
			if(apsi/psi_max<mtheta_data(1,1)) then
				answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
			else
				answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1:ibreak_th-mtheta_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if((Broot==4).or.(Broot==5)) then
	! this case treted first to avoid messy "if" groups

		answer = mach_thetahat(psi,zone_loc)*bzero(psi,zone_loc)/b_phi_zero
		return

	endif

	if (numerical_mtheta) then

		if ((apsi/psi_max)<=mtheta_data(1,1)) then
			answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
		elseif((apsi/psi_max)>=1.d0) then
			answer =  mach_theta_max/mach_theta_num *  &
					mtheta_data(ibreak_th-mtheta_ord,2)
		else
			answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1:ibreak_th-mtheta_ord) )
		endif

	else

			! -------------------- t-2t shape --------------------

		if(apsi<=0.d0) then
			answer = mach_theta_edge
		elseif(apsi/psi_max>=2.d0*t_mth) then
			answer = 0.d0
		elseif(apsi/psi_max<=t_mth) then
			answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
					(2.d0/t_mth*apsi/psi_max - (apsi/psi_max/t_mth)**2)
		elseif(apsi/psi_max>t_mth) then
			answer = mach_theta_max *  &
					(2.d0*t_mth-apsi/psi_max)**2 * (2.d0*apsi/psi_max-t_mth) / t_mth**3
		endif

	endif

	continue

end function mach_theta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dmach_thetadpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi
	real (kind=dkind) :: b0,d,p,b0c,dc,pc
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
	

	if(psi==psi_flag_dep2) then
		answer = mthp_loc
		return
	endif

    if (zone_loc==-1) then

		if(numerical_mtheta) then

			apsi=-dabs(psi)

			if(apsi/psi_max<mtheta_data(1,1)) then
				answer = 0.d0
			else
				answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d0-1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(:) )
			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if((Broot==4).or.(Broot==5)) then
	! this case treted first to avoid messy "if" groups

		answer = ( mach_thetahat(psi,zone_loc)*dbzerodpsi(psi,zone_loc)  + &
							dmach_thetahatdpsi(psi,zone_loc)*bzero(psi,zone_loc) )  &
					/b_phi_zero

		return

	endif

	if (numerical_mtheta) then

		if (apsi/psi_max<=mtheta_data(1,1)) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,mtheta_data(1,1), mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(:) )

		elseif (apsi>=psi_max) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d0-1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(:) )

		else

			answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(:) )

		endif

		answer=answer/psi_max

	 else

			! -------------------- t-2t shape --------------------

		if(apsi<=0.d0) then
			answer = 0.d0
		elseif(apsi/psi_max>=2.d0*t_mth) then
			answer = 0.d0
		elseif(apsi/psi_max<=t_mth) then
			answer = (mach_theta_max-mach_theta_edge) *  &
					2.d0/t_mth/psi_max*(1.d0 - (apsi/psi_max/t_mth))
		elseif(apsi/psi_max>t_mth) then
			answer = 6.d0*mach_theta_max *  &
					(apsi/psi_max-2.d0*t_mth)*(apsi/psi_max-t_mth)/t_mth**3/psi_max
		endif

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dmach_thetadpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function mach_thetahat(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi, xx
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep2) then
		answer = mth_loc
		return
	endif

    if (zone_loc==-1) then
		if(numerical_mtheta) then
			apsi=-dabs(psi)
			if(apsi/psi_max<mtheta_data(1,1)) then
				answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
			else
				answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1:ibreak_th-mtheta_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if (numerical_mtheta) then

		if ((apsi/psi_max)<=mtheta_data(1,1)) then
			answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
		elseif((apsi/psi_max)>=1.d0) then
			answer =  mach_theta_max/mach_theta_num *  &
					mtheta_data(ibreak_th-mtheta_ord,2)
		else
			answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1:ibreak_th-mtheta_ord) )
		endif

	else

			! -------------------- t shape --------------------

		xx = apsi/psi_max

		if(xx<=t_mth) then
			answer = mach_theta_edge
		elseif(xx >= t_mth + 2.d0*w_mth) then
			answer = 0.d0
		elseif(xx > t_mth + w_mth) then !f2(x)
			answer = (mach_theta_max) *  &
					 (2.d0*xx-2.d0*t_mth-w_mth)*(t_mth+2.d0*w_mth-xx)**2/w_mth**3
		elseif(xx>t_mth) then !f1(x)
			answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
							(2.d0*t_mth+3.d0*w_mth-2.d0*xx)*(t_mth-xx)**2/w_mth**3
		endif

	endif

	continue

end function mach_thetahat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dmach_thetahatdpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi, xx
	real (kind=dkind) :: b0,d,p,b0c,dc,pc
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep2) then
		answer = mthp_loc
		return
	endif

    if (zone_loc==-1) then

		if(numerical_mtheta) then

			apsi=-dabs(psi)

			if(apsi/psi_max<mtheta_data(1,1)) then
				answer = 0.d0
			else
				answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(:) )
			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if (numerical_mtheta) then

		if (apsi/psi_max<=mtheta_data(1,1)) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,mtheta_data(1,1), mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(:) )

		elseif (apsi>=psi_max) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d0-1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(:) )

		else

			answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psi_max, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(:) )

		endif

		answer=answer/psi_max

	 else


			! -------------------- t shape --------------------
		xx = apsi/psi_max

		if(xx<=t_mth) then
			answer = 0.d0
		elseif(xx >= t_mth + 2.d0*w_mth) then
			answer = 0.d0
		elseif(xx > t_mth + w_mth) then !f2(x)
			answer = - (mach_theta_max) *  &
					 6.d0*(2.d0*w_mth-xx+t_mth)*(xx-w_mth-t_mth)/w_mth**3
		elseif(xx>t_mth) then !f1(x)
			answer = (mach_theta_max-mach_theta_edge) *  &
							6.d0*(xx-t_mth)*(t_mth+w_mth-xx)/w_mth**3
		endif

		answer = answer/psi_max

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dmach_thetahatdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function mach_phi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: cs,omeg,apsi, mth
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep2) then
		answer = mph_loc
		return
	endif

     if (zone_loc==-1) then

		if((numerical_omega).or.(numerical_mphi)) then
			apsi=-dabs(psi)
			
			if(numerical_omega) then

				if(apsi/psi_max<mphi_data(1,1)) then
					omeg = mphi_data(1,2)
				else
					omeg = dbsval(apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1:ibreak_fi-mphi_ord) )
				endif

				mth = mach_theta(psi,zone_loc)

				if(eq_type<10) then
					cs = dsqrt(gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc))
				elseif(eq_type>=10) then
					cs = dsqrt((gamma_i*p_i_ofpsi(psi,zone_loc)+gamma_e*p_e_ofpsi(psi,zone_loc))  &
						/(mass*d_TF_ofpsi(psi,zone_loc)))
				endif

				answer = rmajor*omeg/cs + mth

			elseif(numerical_mphi) then

				if(apsi/psi_max<mphi_data(1,1)) then
					answer = mphi_data(1,2)
				else

					answer = dbsval(apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1:ibreak_fi-mphi_ord) )

				endif

			endif

			return
			
		else
			apsi = 0.d0
		endif
		
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if ((numerical_omega).or.(Kepler_Omega)) then

		if(numerical_omega) then

			if((apsi/psi_max) >= 1.d0) then
				omeg = mphi_data(ibreak_fi-mphi_ord,2)
			elseif((apsi/psi_max) < mphi_data(1,1)) then
				omeg = mphi_data(1,2)
			else
				omeg = dbsval(apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
						ibreak_fi-mphi_ord, mphi_cscoef(1:ibreak_fi-mphi_ord) )

			endif

		elseif(Kepler_Omega) then

			omeg = sqrt(G_gravity*M_gravity/rmajor**3)

		endif

		mth = mach_theta(psi,zone_loc)

		if(eq_type<10) then
			cs = dsqrt(gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc))
		elseif(eq_type>=10) then
			cs = dsqrt((gamma_i*p_i_ofpsi(psi,zone_loc)+gamma_e*p_e_ofpsi(psi,zone_loc))  &
				/(mass*d_TF_ofpsi(psi,zone_loc)))
		endif

		answer = rmajor*omeg/cs + mth
		! THIS NEEDS TO BE FIXED FOR THE TFCASE! CONTROLLARE!!!

    elseif (numerical_mphi) then

		if((apsi/psi_max) > 1.d0) then
			answer = mphi_data(ibreak_fi-mphi_ord,2)
		elseif((apsi/psi_max) < mphi_data(1,1)) then
			answer = mphi_data(1,2)
		else

			answer = dbsval(apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
					ibreak_fi-mphi_ord, mphi_cscoef(1:ibreak_fi-mphi_ord) )

		endif

	else

		if(psi/psi_max <= fraction .or. apsi/psi_max > 1.d0) then

		   answer = mphi_min

		else

				answer = mphi_min + (mach_phi_max-mphi_min)*  &
						( apsi/psi_max-fraction)**alpha_mphi/ (1.d0-fraction)**alpha_mphi

		end if

!		if(psi/psi_max <= fraction .or. apsi/psi_max > 1.d0) then
!
!		   answer = 0.d0
!
!		else
!
!				answer = mach_phi_max*  &
!						( apsi/psi_max-fraction +mphi_min)**alpha_mphi/ (1.d0-fraction)**alpha_mphi
!
!		end if

	endif

	continue

end function mach_phi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dmach_phidpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: omeg,d_omeg
	real (kind=dkind) :: cs,dcs,p,d,apsi, dmth, p_e, p_i, dTF
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep2) then
		answer = mphp_loc
		return
	endif

    if (zone_loc==-1) then

		if((numerical_omega).or.(numerical_mphi)) then

			apsi=-dabs(psi)

			if (numerical_omega) then

				if(apsi/psi_max<mphi_data(1,1)) then
					answer = 0.d0
				else
					omeg = dbsval(apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
						ibreak_fi-mphi_ord, mphi_cscoef(1:ibreak_fi-mphi_ord) )

					d_omeg = dbsder(1,apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
						ibreak_fi-mphi_ord, mphi_cscoef(:) )
				endif

				d_omeg=d_omeg/psi_max

				if(eq_type<10) then
					p = pofpsi(psi,zone_loc)
					d = dofpsi(psi,zone_loc)
					cs = dsqrt(gamma*p/d)
					dcs = gamma*(-p*dddpsi(psi,zone_loc)+d*dpdpsi(psi,zone_loc))/  &
						(2.d0*d**2*cs)
				elseif(eq_type>=10) then
					p_e = p_e_ofpsi(psi,zone_loc)
					p_i = p_e_ofpsi(psi,zone_loc)
					dTF = d_TF_ofpsi(psi,zone_loc)
					cs = dsqrt((gamma_e*p_e+gamma_i*p_i)/(mass*dTF))
					dcs = ((gamma_e*dp_e_dpsi(psi,zone_loc) + gamma_i*dp_i_dpsi(psi,zone_loc))*dTF -  &
							(gamma_e*p_e + gamma_i*p_i)*dd_TF_dpsi(psi,zone_loc))/(2.d0*cs*mass*dTF**2)
				endif

				if (apsi>0.d0) dcs = dcs*psi/apsi

				dmth = dmach_thetadpsi(psi,zone_loc)

				answer = rmajor*(d_omeg/cs-omeg/cs**2*dcs) + dmth

			elseif (numerical_mphi) then

				if(apsi/psi_max<mphi_data(1,1)) then

					answer = 0.d0

				else

					answer = dbsder(1,apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
									ibreak_fi-mphi_ord, mphi_cscoef(:) )

				endif

				answer=answer/psi_max

			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if ((numerical_omega).or.(Kepler_Omega)) then

		if (numerical_omega) then

			if((apsi/psi_max) >= 1.d0) then

				omeg = mphi_data(ibreak_fi-mphi_ord,2)

				d_omeg = dbsder(1,1.d0, mphi_ord, mphi_data(1:ibreak_fi,3),  &
				ibreak_fi-mphi_ord, mphi_cscoef(:) )

			elseif((apsi/psi_max) < mphi_data(1,1)) then

				omeg = omega_data(1,2)

				d_omeg = 0.d0

			else

				omeg = dbsval(apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
				ibreak_fi-mphi_ord, mphi_cscoef(1:ibreak_fi-mphi_ord) )

				d_omeg = dbsder(1,apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
				ibreak_fi-mphi_ord, mphi_cscoef(:) )

			endif

		elseif(Kepler_Omega) then

			omeg = sqrt(G_gravity*M_gravity/rmajor**3)
			d_omeg = 0.d0

		endif

		d_omeg=d_omeg/psi_max

		if(eq_type<10) then
			p = pofpsi(psi,zone_loc)
			d = dofpsi(psi,zone_loc)
			cs = dsqrt(gamma*p/d)
			dcs = gamma*(-p*dddpsi(psi,zone_loc)+d*dpdpsi(psi,zone_loc))/  &
				(2.d0*d**2*cs)
		elseif(eq_type>=10) then
			p_e = p_e_ofpsi(psi,zone_loc)
			p_i = p_e_ofpsi(psi,zone_loc)
			dTF = d_TF_ofpsi(psi,zone_loc)
			cs = dsqrt((gamma_e*p_e+gamma_i*p_i)/(mass*dTF))
			dcs = ((gamma_e*dp_e_dpsi(psi,zone_loc) + gamma_i*dp_i_dpsi(psi,zone_loc))*dTF -  &
					(gamma_e*p_e + gamma_i*p_i)*dd_TF_dpsi(psi,zone_loc))/(2.d0*cs*mass*dTF**2)
		endif

		if (apsi>0.d0) dcs = dcs*psi/apsi

		dmth = dmach_thetadpsi(psi,zone_loc)

		answer = rmajor*(d_omeg/cs-omeg/cs**2*dcs) + dmth

	elseif (numerical_mphi) then

		if(apsi/psi_max>1.d0) then

			answer = dbsder(1,1.d0, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(:) )

		elseif(apsi/psi_max<mphi_data(1,1)) then

			answer = dbsder(1,mphi_data(1,1), mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(:) )

		else

			answer = dbsder(1,apsi/psi_max, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(:) )

		endif

		answer=answer/psi_max

	 else

		if(psi <= fraction*psi_max .or. apsi > psi_max) then

			answer = 0.0d0

		else

				answer = alpha_mphi * (mach_phi_max-mphi_min) *  &
					( apsi/psi_max-fraction )**(alpha_mphi-1.d0) / psi_max/ (1.d0-fraction)**alpha_mphi

		end if

!		if(psi <= fraction*psi_max .or. apsi > psi_max) then
!
!			answer = 0.0d0
!
!		else
!
!				answer = alpha_mphi * mach_phi_max *  &
!					( apsi/psi_max-fraction + mphi_min)**(alpha_mphi-1.d0) / psi_max/ (1.d0-fraction)**alpha_mphi
!
!		end if

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

	continue

end function dmach_phidpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function bpol(psi,i_ext,j_ext,scale,nxx,nzz) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer, intent (in) :: i_ext,j_ext
	real (kind=dkind), intent (in) :: scale
	integer, intent (in) :: nxx,nzz
	real (kind=dkind) :: answer
    real (kind=dkind), dimension(1:3,1:3), intent(in) :: psi
	real (kind=dkind) :: rloc,bx,bz
	integer i,j

!	print*, 'warning, function bpol needs to be fixed'

	  i=2
	  j=2

	  rloc = x_coord(i_ext)

	  ! Calculate B_r
	  if(j_ext == 1) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_r as independent of z near the boundary.
		 bx = (psi(i,j) - psi(i,j+1))/dz_a(1)/rloc
	  else if(j_ext == nzz) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_r as independent of z near the boundary.
		 bx = (psi(i,j-1) - psi(i,j))/dz_a(nzz)/rloc
	  else
		 ! use a centered difference
		 bx = ( dz_a(j_ext-1)**2*psi(i,j+1) +  &
					 (dz_a(j_ext)**2-dz_a(j_ext-1)**2)*psi(i,j) -  &
					 dz_a(j_ext)**2*psi(i,j-1) ) /  &
					 (dz_a(j_ext)*dz_a(j_ext-1)*(dz_a(j_ext)+dz_a(j_ext-1)))/rloc
	  end if

	  ! Calculate Bz
	  if(i_ext == 1) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_z as independent of x near the boundary.
		 bz = (psi(i+1,j) - psi(i,j))/dx_a(1)/(rloc + 0.5d0*dx_a(1))
	  else if(i_ext == nxx) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_z as independent of x near the boundary.
		 bz = (psi(i,j) - psi(i-1,j))/dx_a(nxx)/(rloc - 0.5d0*dx_a(nxx))
	  else
		 ! use a centered difference
		 bz = ( dx_a(i_ext-1)**2*psi(i+1,j) +  &
					 (dx_a(i_ext)**2-dx_a(i_ext-1)**2)*psi(i,j) -  &
					 dx_a(i_ext)**2*psi(i-1,j) ) /  &
					 (dx_a(i_ext)*dx_a(i_ext-1)*(dx_a(i_ext)+dx_a(i_ext-1)))/rloc
	  end if

	answer = (dsqrt(bx**2+bz**2))/scale
	!"scale" factors in the half/full grid size step

  continue

end function bpol

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function grav_potential(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only :   gravity_type, G_gravity, M_gravity

	real(kind=dkind) :: x, z, answer

	answer = 0.d0

	if(gravity_type==1) then
	! point mass

		answer = -G_gravity*M_gravity/sqrt(x**2+z**2)

	elseif(gravity_type==2) then
	! constant gravity (-R direction)

		answer = -G_gravity*M_gravity*x

	endif

end function grav_potential

  ! ------------------------------------------------------------------
  ! The following set of function are used to construct Hameiri's
  ! functions of Psi from Betti & Freidberg's set.
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function sofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer, p, d
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = s_loc
		return
	endif

    p = pofpsi(psi,zone_loc)
    d = dofpsi(psi,zone_loc)

    answer = p/d**gamma

    return

  end function sofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dsdpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: d
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = sp_loc
		return
	endif

    d = dofpsi(psi,zone_loc)

    answer = (dpdpsi(psi,zone_loc) - gamma*pofpsi(psi,zone_loc)*dddpsi(psi,zone_loc)/d)/d**gamma

    return

  end function dsdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function phiofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: b0
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = phi_loc
		return
	endif

	if(F_opt==0) then
	! bzero is 0, let's start with no poloidal flow

		answer = 0.d0

	else

		if((Broot==4).or.(Broot==5)) then
		! RFP case, bzero can be 0, better replace it with a constant

			b0 = b_phi_zero
			answer = dsqrt(mu_mag*gamma*pofpsi(psi)*dofpsi(psi))* &
				 mach_thetahat(psi)/b0

		else

			b0 = bzero(psi,zone_loc)

			answer = dsqrt(mu_mag*gamma*pofpsi(psi,zone_loc)*dofpsi(psi,zone_loc))* &
				 mach_theta(psi,zone_loc)/b0

		endif

	endif

    return

  end function phiofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dphidpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: d,p,b,m,db
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif


	if(psi==psi_flag_ham) then
		answer = phip_loc
		return
	endif

	if(F_opt==0) then
	! bzero is 0, let's start with no poloidal flow

		answer = 0.d0

	else

		d = dofpsi(psi,zone_loc)
		p = pofpsi(psi,zone_loc)

		if((Broot==4).or.(Broot==5)) then
		! RFP case, bzero can be 0, better replace it with a constant

			b = b_phi_zero

			m = mach_thetahat(psi,zone_loc)

			answer = dsqrt(mu_mag*gamma*p*d)/b*( &
				dmach_thetahatdpsi(psi,zone_loc) &
				+ 0.5d0*m*(dddpsi(psi,zone_loc)/d + dpdpsi(psi,zone_loc)/p) &
				)

		else

			b = bzero(psi,zone_loc)

			if(d > 0.0d0 .and. dabs(b) > 0.0d0 .and. p > 0.0d0) then
			   m = mach_theta(psi,zone_loc)

				answer = dsqrt(mu_mag*gamma*p*d)/b*( &
					dmach_thetadpsi(psi,zone_loc) &
					+ 0.5d0*m*(dddpsi(psi,zone_loc)/d + dpdpsi(psi,zone_loc)/p) &
					- dbzerodpsi(psi,zone_loc)*m/b)

			else
			   answer = 0.0d0
			end if

		endif

	endif

    return

end function dphidpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function omegaofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

!!$	if(psi==psi_flag_ham) then
!!$		answer = omega_loc
!!$		return
!!$	endif

	if(eq_type<10) then

		answer = dsqrt(gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc))* &
						(mach_phi(psi,zone_loc) - mach_theta(psi,zone_loc))/rmajor

	elseif(eq_type>=10) then

		answer = dsqrt((gamma_e*p_e_ofpsi(psi,zone_loc)+gamma_i*p_i_ofpsi(psi,zone_loc))  &
			/(mass*d_TF_ofpsi(psi,zone_loc)))* &
				(mach_phi(psi,zone_loc) - mach_theta(psi,zone_loc))/rmajor

	endif

    return

  end function omegaofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function domegadpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: p,d,p_i,p_e, dTF
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif


!!$	if(psi==psi_flag_ham) then
!!$		answer = omegap_loc
!!$		return
!!$	endif

	if(eq_type<10) then

		p = pofpsi(psi,zone_loc)

		if(p > 0.0d0) then

		   d = dofpsi(psi,zone_loc)

		   answer = (0.5d0*dsqrt(gamma/(p*d))*(dpdpsi(psi,zone_loc) -dddpsi(psi,zone_loc)*p/d)* &
				(mach_phi(psi,zone_loc) - mach_theta(psi,zone_loc)) + &
				dsqrt(gamma*p/d)*(dmach_phidpsi(psi,zone_loc) - dmach_thetadpsi(psi,zone_loc)))/ &
				rmajor

		else
		   answer = 0.0d0
		end if

	elseif(eq_type>=10) then

		p_e = p_e_ofpsi(psi,zone_loc)
		p_i = p_i_ofpsi(psi,zone_loc)

		if(p_i+p_e > 0.0d0) then

			dTF = d_TF_ofpsi(psi,zone_loc)

!			answer = (gamma*(2*d*(p_e + p_i)*(dmach_phidpsi(psi) - dmach_thetadpsi(psi)) +  &
!							(mach_phi(psi) - mach_theta(psi))*(-((p_e + p_i)*dd_TF_dpsi(psi)) +  &
!							d*(dp_e_dpsi(psi) + dp_e_dpsi(psi)))))/(2.*rmajor*d**2*sqrt((gamma*(p_e + p_i))/d))
!			answer = answer / sqrt(mass)

		answer = (2*(gamma_e*p_e + gamma_i*p_i)*(dmach_phidpsi(psi,zone_loc) -   &
			dmach_thetadpsi(psi,zone_loc))*dTF +  &
			(-(gamma_e*p_e*dd_TF_dpsi(psi,zone_loc)) - gamma_i*p_i*dd_TF_dpsi(psi,zone_loc) +  &
			(gamma_e*dp_e_dpsi(psi,zone_loc) + gamma_i*dp_i_dpsi(psi,zone_loc))*dTF)*  &
			(mach_phi(psi,zone_loc) - mach_theta(psi,zone_loc)))/  &
			(2.*mass*rmajor*sqrt((gamma_e*p_e + gamma_i*p_i)/(mass*dTF))*dTF**2)

		else
		   answer = 0.0d0
		end if

	endif

    return

  end function domegadpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function iofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = i_loc
		return
	endif

    answer = rmajor*bzero(psi,zone_loc)/dsqrt(mu_mag)

    return

  end function iofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function didpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = ip_loc
		return
	endif

 		  answer = rmajor*dbzerodpsi(psi,zone_loc)/dsqrt(mu_mag)

    return

  end function didpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function hofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: mt,mp
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = h_loc
		return
	endif

    mt = mach_theta(psi,zone_loc)
    mp = mach_phi(psi,zone_loc)

	answer = (1.0d0/(gamma-1.0d0) + mt*mp - 0.5d0*mp**2)* &
		gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc)


    if(answer <= 0.0d0)then
		   print *, "Error H(Psi) = ",answer
		   pause
		   stop
	end if

    return

  end function hofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dhdpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: d,p,mt,mp
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = hp_loc
		return
	endif

    d = dofpsi(psi,zone_loc)
    p = pofpsi(psi,zone_loc)
    mt = mach_theta(psi,zone_loc)
    mp = mach_phi(psi,zone_loc)


    if(d>0.0d0) then

	   answer = (1.0d0/(gamma-1.0d0) + mt*mp - 0.5d0*mp**2)* &
            gamma*(dpdpsi(psi,zone_loc) - dddpsi(psi,zone_loc)*p/d)/d + &
            ((mt-mp)*dmach_phidpsi(psi,zone_loc) + mp*dmach_thetadpsi(psi,zone_loc))* &
            gamma*p/d

	else
       answer = 0.0d0
    end if

    return

  end function dhdpsi

  ! ------------------------------------------------------------------
  !					end of one-fluid free functions
  ! ------------------------------------------------------------------

!-------------------------------------------------------------------
!							two-fluid free functions
!-------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function d_TF_ofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! due to differente definitions, what I am calling "D" in the TF derivation has differente dimensions with respect
! to the one-fluid "D". Thus, it is preferable to redefine D and avoid confusion.
! Due to the quasi-neutrality assumption, there is still only one "D"; thus, we may use the same arrays as in
! the one-fluid case.

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer,apsi
    integer, optional, intent (in) :: zone
    integer :: zone_loc

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(eq_type<10) then
		if(D_TF_num) then
			continue
		else
			answer = dofpsi(psi,zone_loc)/mass
			return
		endif
	endif


!!$	if(psi==psi_flag) then
!!$		answer = d_loc
!!$		return
!!$	elseif(psi==psic_flag) then
!!$		answer = dc_loc
!!$		return
!!$	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then
		if(numerical_n) then
			apsi=-dabs(psi)
			if(apsi/psi_max<d_data(1,1)) then
				answer = d_data(1,2)
			else
				answer = dbsval(apsi/psi_max, d_ord, d_data(1:ibreak_d,3),  &
					ibreak_d-d_ord, d_cscoef(1:ibreak_d-d_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif

	if (numerical_n) then

		if ((apsi/psi_max)<=d_data(1,1)) then
			answer = d_data(1,2)
		elseif (apsi>=psi_max) then
			answer = d_data(ibreak_d-d_ord,2)
		else

			answer = dbsval(apsi/psi_max, d_ord, d_data(1:ibreak_d,3),  &
					ibreak_d-d_ord, d_cscoef(1:ibreak_d-d_ord) )

		endif

	else

	    if(psi <= psi_max*fraction ) then
		   answer = d_TF_edge
		else if(dabs(psi) > psi_max) then
	       answer = d_TF_center
		else
		   answer = d_TF_edge + (d_TF_center - d_TF_edge)*  &
				dabs( psi/psi_max-fraction )**alpha_rho / ( 1.d0-fraction )**alpha_rho
		end if

	endif


	continue

end function d_TF_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dd_TF_dpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer,apsi
    integer, optional, intent (in) :: zone
    integer :: zone_loc

!!$	if(psi==psi_flag) then
!!$		answer = dp_loc
!!$		return
!!$	elseif(psi==psic_flag) then
!!$		answer = dcp_loc
!!$		return
!!$	endif

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(eq_type<10) then
		if(D_TF_num) then
			continue
		else
			answer = dddpsi(psi,zone_loc)/mass
			return
		endif
	endif


	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_n) then

			apsi=-dabs(psi)

			if(apsi/psi_max<d_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psi_max, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )
			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif
    
    if (numerical_n) then

			if ((apsi/psi_max)<=d_data(1,1)) then

				answer = dbsder(1,d_data(1,1), d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )

			elseif (apsi>=psi_max) then

				answer = dbsder(1,1.d0-1.d-6, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )

			else

				answer = dbsder(1,apsi/psi_max, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(:) )

			endif

			answer=answer/psi_max

	else

	    if(psi <= psi_max*fraction ) then
		   answer = 0.0d0
	    else if(dabs(psi) > psi_max) then
		   answer = 0.0d0
	    else
		   answer = alpha_rho*(d_TF_center - d_TF_edge)/psi_max  &
					*dabs( psi/psi_max-fraction )**(alpha_rho-1.d0)/( 1.d0-fraction )**alpha_rho
		end if

	endif

	if(abs(psi)>0.d0) answer = answer*psi/abs(psi)

	continue

end function dd_TF_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function p_e_ofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! electron pressure free function

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer,apsi
    integer, optional, intent (in) :: zone
    integer :: zone_loc

!!$	if(psi==psi_flag) then
!!$		answer = p_loc
!!$		return
!!$	elseif(psi==psic_flag) then
!!$		answer = pc_loc
!!$		return
!!$	endif

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then
		if(numerical_p_e) then
			apsi=-dabs(psi)
			if(apsi/psi_max<p_e_data(1,1)) then
				answer = p_e_data(1,2)
			else
				answer = dbsval(apsi/psi_max, p_e_ord, p_e_data(1:ibreak_p_e,3),  &
						ibreak_p_e-p_e_ord, p_e_cscoef(1:ibreak_p_e-p_e_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif
    
	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

if (numerical_p_e) then

			if ((apsi/psi_max)<=p_e_data(1,1)) then
				answer = p_e_data(1,2)
			elseif (apsi>=psi_max) then
				answer = p_e_data(ibreak_p_e-p_e_ord,2)
			else

				answer = dbsval(apsi/psi_max, p_e_ord, p_e_data(1:ibreak_p_e,3),  &
						ibreak_p_e-p_e_ord, p_e_cscoef(1:ibreak_p_e-p_e_ord) )

			endif

	else

		if (p_e_opt == 10) then
		! RFP profile

			answer = p_e_edge + (p_e_center-p_e_edge) * (apsi/psic_13)**2 *  &
				 (6.d0 - 8.d0*apsi/psic_13 + 3.d0*(apsi/psic_13)**2)

		else

			if(apsi/psic_13 > 1.0d0) then
				answer = p_e_center
			elseif(psi >= (psi_max*fraction) ) then
				answer = p_e_edge + (p_e_center-p_e_edge)* (   (apsi/psic_13 - &
												fraction)/ (1.d0-fraction)   )**alpha_e
			else
				answer = p_e_edge
			end if

		endif

	endif

	if(answer<=0.d0) answer = p_e_edge*1.d-10

	continue

end function p_e_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dp_e_dpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi
	integer, optional, intent (in) :: zone
    integer :: zone_loc

!!$	if(psi==psi_flag) then
!!$		answer = pp_loc
!!$		return
!!$	elseif(psi==psic_flag) then
!!$		answer = pcp_loc
!!$		return
!!$	endif

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_p_e) then

			apsi=-dabs(psi)

			if(apsi/psi_max<p_e_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psi_max, p_e_ord, p_e_data(1:ibreak_p_e,3),  &
								ibreak_p_e-p_e_ord, p_e_cscoef(:) )
			endif

			answer=answer/psi_max
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

	if (numerical_p_e) then

			if ((apsi/psi_max)<=p_e_data(1,1)) then

				answer = dbsder(1,p_e_data(1,1), p_e_ord, p_e_data(1:ibreak_p_e,3),  &
								ibreak_p_e-p_e_ord, p_e_cscoef(:) )

			elseif (apsi>=psi_max) then

				answer = dbsder(1,1.d0-1.d-6, p_e_ord, p_e_data(1:ibreak_p_e,3),  &
								ibreak_p_e-p_e_ord, p_e_cscoef(:) )

			else

				answer = dbsder(1,apsi/psi_max, p_e_ord, p_e_data(1:ibreak_p_e,3),  &
								ibreak_p_e-p_e_ord, p_e_cscoef(:) )

			endif

		answer=answer/psi_max

	else

		if (p_e_opt == 10) then
			! RFP profile

			answer = 12.d0 * (p_e_center-p_e_edge) * (1.d0-apsi/psic_13)**2 * apsi/psic_13 /psic_13

		else

			if(dabs(psi) > dabs(psic_13).or.(dabs(psi/psi_max)==0.d0) ) then
				answer = 0.d0
			elseif(psi >= psi_max*fraction ) then
			   answer = alpha_e*(p_e_center-p_e_edge)/psic_13*  &
						( apsi/psic_13-fraction )**(alpha_e-1.0d0)/ (1.d0-fraction)**alpha_e
			else
				answer = 0.0d0
			end if

		endif

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

	continue

end function dp_e_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function p_i_ofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! ion pressure free function
! since big_Psi is not set to 0 on the boundary, range and absolute values are treated differently from the one-
! fluid version

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
	real(kind=dkind) :: apsi

!!$	if(psi==psi_flag) then
!!$		answer = p_loc
!!$		return
!!$	elseif(psi==psic_flag) then
!!$		answer = pc_loc
!!$		return
!!$	endif

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

!!$	if((bc_type/=7).and.(bc_type/=8)) then
!!$		apsi=dabs(psi)
!!$	else
!!$		apsi = psi
!!$	endif

	apsi = psi

    if (zone_loc==-1) then
		if(numerical_p_i) then
			apsi=-dabs(psi)
			if(apsi/psi_max<p_i_data(1,1)) then
				answer = p_i_data(1,2)
			else
				answer = dbsval(apsi/psi_max, p_i_ord, p_i_data(1:ibreak_p_i,3),  &
						ibreak_p_i-p_i_ord, p_i_cscoef(1:ibreak_p_i-p_i_ord) )
			endif
			return
		else
			apsi = 0.d0
		endif
    endif

    
	if (numerical_p_i) then

			if ((apsi/psi_max)<=p_i_data(1,1)) then
				answer = p_i_data(1,2)
			elseif (apsi>=psi_max) then
				answer = p_i_data(ibreak_p_i-p_i_ord,2)
			else

				answer = dbsval(apsi/psi_max, p_i_ord, p_i_data(1:ibreak_p_i,3),  &
						ibreak_p_i-p_i_ord, p_i_cscoef(1:ibreak_p_i-p_i_ord) )

			endif

	else

		if (p_i_opt == 10) then
		! RFP profile

			answer = p_i_edge + (p_i_center-p_i_edge) * (apsi/psi_max)**2 *  &
				 (6.d0 - 8.d0*apsi/psi_max + 3.d0*(apsi/psi_max)**2)

		else

			if(apsi/psi_max > 1.0d0) then
				answer = p_i_center
			elseif(apsi >= (psi_max*fraction) ) then
				answer = p_i_edge + (p_i_center-p_i_edge)* (   (apsi/psi_max- &
												fraction)/ (1.d0-fraction)   )**alpha_i
			else
				answer = p_i_edge
			end if

		endif

	endif

	if(answer<=0.d0) answer = p_i_edge*1.d-10

	continue

end function p_i_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dp_i_dpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
	real(kind=dkind) :: apsi

!!$	if(psi==psi_flag) then
!!$		answer = pp_loc
!!$		return
!!$	elseif(psi==psic_flag) then
!!$		answer = pcp_loc
!!$		return
!!$	endif

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

!!$	if((bc_type/=7).and.(bc_type/=8)) then
!!$		apsi=dabs(psi)
!!$	else
		apsi = psi
!!$	endif

    if (zone_loc==-1) then

		if(numerical_p_i) then

			apsi=-dabs(psi)

			if(apsi/psi_max<p_i_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psi_max, p_i_ord, p_i_data(1:ibreak_p_i,3),  &
								ibreak_p_i-p_i_ord, p_i_cscoef(:) )
			endif

			answer=answer/psi_max
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

	if (numerical_p_i) then

			if ((apsi/psi_max)<=p_i_data(1,1)) then

				answer = dbsder(1,p_i_data(1,1), p_i_ord, p_i_data(1:ibreak_p_i,3),  &
								ibreak_p_i-p_i_ord, p_i_cscoef(:) )

			elseif (apsi>=psi_max) then

				answer = dbsder(1,1.d0-1.d-6, p_i_ord, p_i_data(1:ibreak_p_i,3),  &
								ibreak_p_i-p_i_ord, p_i_cscoef(:) )

			else

				answer = dbsder(1,apsi/psi_max, p_i_ord, p_i_data(1:ibreak_p_i,3),  &
								ibreak_p_i-p_i_ord, p_i_cscoef(:) )

			endif

		answer=answer/psi_max

	else

		if (p_i_opt == 10) then
			! RFP profile

				answer = 12.d0 * (p_i_center-p_i_edge) * (1.d0-apsi/psi_max)**2 * apsi/psi_max/psi_max

		else

			if(dabs(apsi) > dabs(psi_max).or.(dabs(apsi/psi_max)==0.d0) ) then
				answer = 0.d0
			elseif(apsi >= psi_max*fraction ) then
			   answer = alpha_i*(p_i_center-p_i_edge)/psi_max*  &
						( apsi/psi_max-fraction )**(alpha_i-1.0d0)/ (1.d0-fraction)**alpha_i
			else
				answer = 0.0d0
			end if

		endif

	endif

	continue

end function dp_i_dpsi

  ! ------------------------------------------------------------------
  ! The following set of function are used to construct Hameiri's form
  ! functions of Psi from the intuitive set. In two-fluid version, of course
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function s_e_ofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer, p, d
    integer, optional, intent (in) :: zone
    integer :: zone_loc

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
!!$	if(psi==psi_flag_ham) then
!!$		answer = s_loc
!!$		return
!!$	endif

    p = p_e_ofpsi(psi,zone_loc)
    d = d_TF_ofpsi(psi,zone_loc)

    answer = p/d**gamma_e / mass

    return

  end function s_e_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function ds_e_dpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: d
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
!!$	if(psi==psi_flag_ham) then
!!$		answer = sp_loc
!!$		return
!!$	endif

    d = d_TF_ofpsi(psi, zone_loc)

    answer = (dp_e_dpsi(psi, zone_loc) - gamma_e*p_e_ofpsi(psi, zone_loc)*dd_TF_dpsi(psi, zone_loc)/d)  &
    	/d**gamma_e / mass

    return

  end function ds_e_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function s_i_ofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer, p, d
    integer, optional, intent (in) :: zone
    integer :: zone_loc

    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

!!$	if(psi==psi_flag_ham) then
!!$		answer = s_loc
!!$		return
!!$	endif

    p = p_i_ofpsi(psi,zone_loc)
    d = d_TF_ofpsi(psi,zone_loc)

    answer = p/d**gamma_i / mass

    return

  end function s_i_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function ds_i_dpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: d
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
!!$	if(psi==psi_flag_ham) then
!!$		answer = sp_loc
!!$		return
!!$	endif

    d = d_TF_ofpsi(psi,zone_loc)

    answer = (dp_i_dpsi(psi,zone_loc) - gamma_i*p_i_ofpsi(psi,zone_loc)*  &
    	dd_TF_dpsi(psi,zone_loc)/d)/d**gamma_i / mass

    return

  end function ds_i_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function phi_TF_ofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real (kind=dkind) :: b0
	integer, optional, intent (in) :: zone
    integer :: zone_loc

    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = phi_loc
		return
	endif

	if(F_opt==0) then
	! bzero is 0, let's start with no poloidal flow

		answer = 0.d0

	else

		if((Broot==4).or.(Broot==5)) then
		! RFP case, bzero can be 0, better replace it with a constant

			b0 = b_phi_zero
			answer = dsqrt((gamma_e*p_e_ofpsi(psi,zone_loc)+gamma_i*p_i_ofpsi(psi,zone_loc)) &
				*d_TF_ofpsi(psi,zone_loc)/mass)* &
				 mach_thetahat(psi,zone_loc)/b0
		else

			b0 = bzero(psi,zone_loc)

		endif

		answer = dsqrt((gamma_e*p_e_ofpsi(psi,zone_loc)+gamma_i*p_i_ofpsi(psi,zone_loc))  &
			*d_TF_ofpsi(psi,zone_loc)/mass)* &
			 mach_theta(psi,zone_loc)/b0

	endif

    return

  end function phi_TF_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dphi_TF_dpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: dtf, dtfp, pe, pi, pep, pip, mt, mtp, b0, b0p
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc

    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif


	if(psi==psi_flag_ham) then
		answer = phip_loc
		return
	endif

	if(F_opt==0) then
	! bzero is 0, let's start with no poloidal flow

		answer = 0.d0

	else

		dtf = d_TF_ofpsi(psi,zone_loc)
		pe = p_e_ofpsi(psi,zone_loc)
		pi = p_i_ofpsi(psi,zone_loc)
		dtfp = dd_TF_dpsi(psi,zone_loc)
		pep = dp_e_dpsi(psi,zone_loc)
		pip = dp_i_dpsi(psi,zone_loc)

		if((Broot==4).or.(Broot==5)) then
		! RFP case, bzero can be 0, better replace it with a constant

			b0 = b_phi_zero
			b0p = 0.d0

			mt = mach_thetahat(psi,zone_loc)
			mtp = dmach_thetahatdpsi(psi,zone_loc)

!				answer = (gamma*(dtfp*b0*mt*(pe + pi) +  &
!								dtf*(2*mtp*b0*(pe + pi) +  &
!								mt*((pep + pip)*b0 - 2*b0p*(pe + pi)))))/  &
!								(2.d0*mass*b0**2*sqrt((gamma*dtf*(pe + pi))/mass))

				answer = (-2.d0*b0p*mt*(gamma_e*pe + dtf*gamma_i*pi) +  &
						b0*(2.d0*gamma_e*mtp*pe + gamma_e*mt*pep + dtfp*gamma_i*mt*pi +  &
						2.d0*dtf*gamma_i*mtp*pi + dtf*gamma_i*mt*pip))/  &
						(2.d0*b0**2*mass*sqrt((gamma_e*pe + dtf*gamma_i*pi)/mass))

		else

			b0 = bzero(psi,zone_loc)
			b0p = dbzerodpsi(psi,zone_loc)

			if(dtf > 0.0d0 .and. dabs(b0) > 0.0d0 .and. (pe+pi) > 0.0d0) then

				mt = mach_theta(psi,zone_loc)
				mtp = dmach_thetadpsi(psi,zone_loc)

!				answer = (gamma*(dtfp*b0*mt*(pe + pi) +  &
!								dtf*(2*mtp*b0*(pe + pi) +  &
!								mt*((pep + pip)*b0 - 2*b0p*(pe + pi)))))/  &
!								(2.d0*mass*b0**2*sqrt((gamma*dtf*(pe + pi))/mass))

				answer = (-2.d0*b0p*dtf*mt*(gamma_e*pe + gamma_i*pi) +  &
								b0*(dtfp*mt*(gamma_e*pe + gamma_i*pi) +  &
								dtf*(2*gamma_e*mtp*pe + gamma_e*mt*pep + 2*gamma_i*mtp*pi +  &
								gamma_i*mt*pip)))/  &
								(2.d0*b0**2*mass*sqrt((dtf*(gamma_e*pe + gamma_i*pi))/mass))

			else
			   answer = 0.0d0
			end if

		endif

	endif

    return

end function dphi_TF_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function phi_e_TF_ofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	answer = phi_TF_ofpsi(psi,zone_loc) - sqrt(mu_mag)*didpsi(psi,zone_loc)/(eV*mu_mag)

	continue

  end function phi_e_TF_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function h_e_ofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: om, p_e, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc

    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

!!$    if (zone_loc==-1) then
!!$		apsi = 0.d0
!!$    endif

	d = d_TF_ofpsi(psi, zone_loc)
	p_e = p_e_ofpsi(psi, zone_loc)

	if ((psi/psi_max)<=0.d0) then
		om = Omega_int_data(1,2)
	elseif((psi/psi_max)>=1.d0) then
		om =  Omega_int_data(n_H,2)
	else
		om = dbsval(psi/psi_max, Omega_int_ord, Omega_int_data(1:ibreak_He,3),  &
									n_H, Omega_int_cscoef(1:n_H) )
	endif

	answer = gamma_e/(gamma_e-1.d0)*p_e/d - eV*om

	return

  end function h_e_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function h_e_ofpsi_partial(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this version does not have the large (electric field) part

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: p_e, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	d = d_TF_ofpsi(psi,zone_loc)
	p_e = p_e_ofpsi(psi,zone_loc)

	answer = gamma_e/(gamma_e-1.d0)*p_e/d

	return

  end function h_e_ofpsi_partial

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function h_i_ofpsi(Psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: mt, mp, om, p_e, p_i, d
	real (kind=dkind) :: answer, apsi
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	d = d_TF_ofpsi(Psi,zone_loc)
	p_e = p_e_ofpsi(Psi,zone_loc)
	p_i = p_i_ofpsi(Psi,zone_loc)
    mp = mach_phi(Psi,zone_loc)
	mt = mach_theta(Psi,zone_loc)

    if (zone_loc==-1) then
		apsi = 0.d0
    endif

	if ((apsi/psi_max)<=Omega_int_data(1,1)) then
		om = Omega_int_data(1,2)
	elseif((apsi/psi_max)>=1.d0) then ! Not sure if this applies. 2/2021
		om =  Omega_int_data(n_H,2)
	else
		om = dbsval(Psi/psi_max, Omega_int_ord, Omega_int_data(1:ibreak_Hi,3),  &
									n_H, Omega_int_cscoef(1:n_H) )
	endif

	answer = ( (gamma_i*p_i+gamma_e*p_e)*(mt*mp-mp**2/2.d0) + gamma_i/(gamma_i-1.d0)*p_i )/d + eV*om

	return

  end function h_i_ofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function h_i_ofpsi_partial(Psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this version does not have the large (electric field) part

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: mt, mp, om, p_e, p_i, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	d = d_TF_ofpsi(Psi,zone_loc)
	p_e = p_e_ofpsi(Psi,zone_loc)
	p_i = p_i_ofpsi(Psi,zone_loc)
    mp = mach_phi(Psi,zone_loc)
	mt = mach_theta(Psi,zone_loc)

!	answer = gamma * ( (p_i+p_e)*(mt*mp-mp**2/2.d0)  + 1.d0/(gamma-1.d0)*p_i )/d
	answer = ( (gamma_i*p_i+gamma_e*p_e)*(mt*mp-mp**2/2.d0) + gamma_i/(gamma_i-1.d0)*p_i )/d

	return

  end function h_i_ofpsi_partial

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function dh_e_dpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: om, p_e, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	d = d_TF_ofpsi(psi,zone_loc)
	p_e = p_e_ofpsi(psi,zone_loc)
	om = omegaofpsi(psi,zone_loc)
	answer = gamma_e/(gamma_e-1.d0) *  &
					(d*dp_e_dpsi(psi,zone_loc)-p_e*dd_TF_dpsi(psi,zone_loc) )/d**2  &
					- eV*om	 &	! Omega part
					+ dp_i_dpsi(psi,zone_loc) / d  ! electrostatic ion confinement

	return

  end function dh_e_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function dh_e_dpsi_partial(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this version does not have the large (electric field) part

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: p_e, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	d = d_TF_ofpsi(psi,zone_loc)
	p_e = p_e_ofpsi(psi,zone_loc)
	answer = gamma_e/(gamma_e-1.d0) *  &
					(d*dp_e_dpsi(psi,zone_loc)-p_e*dd_TF_dpsi(psi,zone_loc) )/d**2

	return

  end function dh_e_dpsi_partial

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function dh_i_dpsi(Psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: Psi
	real (kind=dkind) :: mp, mt, om, p_e, p_i, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	mt = mach_theta(Psi,zone_loc)
	mp = mach_phi(Psi,zone_loc)
	d = d_TF_ofpsi(Psi,zone_loc)
	p_e = p_e_ofpsi(Psi,zone_loc)
	p_i = p_i_ofpsi(Psi,zone_loc)
	om = omegaofpsi(Psi,zone_loc)

!	answer = gamma/(gamma-1.d0) * (  &
!					d*(p_e+p_i)*(dmach_phidpsi(Psi)*(mt-mp)+mp*dmach_thetadpsi(Psi)) +  &
!					(d*dp_i_dpsi(Psi)-p_i*dd_TF_dpsi(Psi))/(gamma-1.d0) +  &
!					mp**2 * ( d*(dp_i_dpsi(Psi)+dp_e_dpsi(Psi)) - dd_TF_dpsi(Psi)*(p_e+p_i) )  &
!					) / d**2  &
!					+ eV*om

!!$	answer = (1.0d0/(gamma-1.0d0) + mt*mp - 0.5d0*mp**2)* &		! ion part
!!$					gamma*(dp_i_dpsi(psi) - dd_TF_dpsi(psi)*p_i/d)/d + &
!!$					((mt-mp)*dmach_phidpsi(psi) + mp*dmach_thetadpsi(psi))* &
!!$					gamma*p_i/d +  &
!!$					( mt*mp - 0.5d0*mp**2)* &		! electron part
!!$					gamma*(dp_e_dpsi(psi) - dd_TF_dpsi(psi)*p_e/d)/d + &
!!$					((mt-mp)*dmach_phidpsi(psi) + mp*dmach_thetadpsi(psi))* &
!!$					gamma*p_e/d  &
!!$					+ eV*om		! Omega part

	answer = (1.0d0/(gamma_i-1.0d0) + mt*mp - 0.5d0*mp**2)* &		! ion part
					gamma_i*(dp_i_dpsi(psi,zone_loc) - dd_TF_dpsi(psi,zone_loc)*p_i/d)/d + &
					((mt-mp)*dmach_phidpsi(psi,zone_loc) + mp*dmach_thetadpsi(psi,zone_loc))* &
					gamma_i*p_i/d +  &
					( mt*mp - 0.5d0*mp**2)* &		! electron part
					gamma_e*(dp_e_dpsi(psi,zone_loc) - dd_TF_dpsi(psi,zone_loc)*p_e/d)/d + &
					((mt-mp)*dmach_phidpsi(psi,zone_loc) + mp*dmach_thetadpsi(psi,zone_loc))* &
					gamma_e*p_e/d  &
					+ eV*om	 &	! Omega part
					- dp_i_dpsi(psi,zone_loc) / d  ! electrostatic ion confinement

	return

  end function dh_i_dpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function dh_i_dpsi_partial(Psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this version does not have the large (electric field) part

	real (kind=dkind), intent (in) :: Psi
	real (kind=dkind) :: om, mp, mt, p_e, p_i, d
	real (kind=dkind) :: answer
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	mt = mach_theta(Psi,zone_loc)
	mp = mach_phi(Psi,zone_loc)
	d = d_TF_ofpsi(Psi,zone_loc)
	p_e = p_e_ofpsi(Psi,zone_loc)
	p_i = p_i_ofpsi(Psi,zone_loc)
	om = omegaofpsi(Psi,zone_loc)

!	answer = gamma/(gamma-1.d0) * (  &
!					d*(p_e+p_i)*(dmach_phidpsi(Psi)*(mt-mp)+mp*dmach_thetadpsi(Psi)) +  &
!					(d*dp_i_dpsi(Psi)-p_i*dd_TF_dpsi(Psi))/(gamma-1.d0) +  &
!					mp**2 * ( d*(dp_i_dpsi(Psi)+dp_e_dpsi(Psi)) - dd_TF_dpsi(Psi)*(p_e+p_i) )  &
!					) / d**2

!	answer = (1.0d0/(gamma-1.0d0) + mt*mp - 0.5d0*mp**2)* &		! ion part
!					gamma*(dp_i_dpsi(psi) - dd_TF_dpsi(psi)*p_i/d)/d + &
!					((mt-mp)*dmach_phidpsi(psi) + mp*dmach_thetadpsi(psi))* &
!					gamma*p_i/d +  &
!					( mt*mp - 0.5d0*mp**2)* &		! electron part
!					gamma*(dp_e_dpsi(psi) - dd_TF_dpsi(psi)*p_e/d)/d + &
!					((mt-mp)*dmach_phidpsi(psi) + mp*dmach_thetadpsi(psi))* &
!					gamma*p_e/d 

	answer = (1.0d0/(gamma_i-1.0d0) + mt*mp - 0.5d0*mp**2)* &		! ion part
					gamma_i*(dp_i_dpsi(psi,zone_loc) - dd_TF_dpsi(psi,zone_loc)*p_i/d)/d + &
					((mt-mp)*dmach_phidpsi(psi,zone_loc) + mp*dmach_thetadpsi(psi,zone_loc))* &
					gamma_i*p_i/d +  &
					( mt*mp - 0.5d0*mp**2)* &		! electron part
					gamma_e*(dp_e_dpsi(psi,zone_loc) - dd_TF_dpsi(psi,zone_loc)*p_e/d)/d + &
					((mt-mp)*dmach_phidpsi(psi,zone_loc) + mp*dmach_thetadpsi(psi,zone_loc))* &
					gamma_e*p_e/d

	return

  end function dh_i_dpsi_partial

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_H_diff(x1,x2,H_diff, zone)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x1, x2, H_diff
	real(kind=dkind) :: dx_H, x_temp
	integer :: j
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(x1==x2) then
		H_diff = 0.d0
		return
	endif

	dx_H = x2-x1

	do j = 1, n_int_H

		x_temp = x1 + x_int_H(j)*dx_H
		fun_int_H(j) = omegaofpsi(x_temp,zone_loc)  &
							+ dp_i_dpsi(x_temp,zone_loc) / (eV * d_TF_ofpsi(x_temp,zone_loc))

	enddo

	call integrate(n_int_H,fun_int_H,w_int_H,H_diff)

	H_diff = H_diff*dx_H	! dx_H is of course the Jacobian

end subroutine get_H_diff

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 function Fstarofpsi(psi,big_Psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: psi, big_Psi, answer
	real(kind=dkind) :: Fterm, vterm, dpsi, x_temp
	integer :: j
	integer, optional, intent (in) :: zone
    integer :: zone_loc
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	Fterm = iofpsi(psi,zone_loc)*sqrt(mu_mag)

	dpsi = big_Psi-psi

	if(dpsi==0.d0) then

		vterm = 0.d0

	else

		do j = 1, n_int_phi

			x_temp = psi + x_int_phi(j)*dpsi
			fun_int_phi(j) = phi_TF_ofpsi(x_temp,zone_loc)

		enddo

		call integrate(n_int_phi,fun_int_phi,w_int_phi,vterm)

		vterm = vterm*dpsi

	endif

	answer = eV*mu_mag*vterm + Fterm

end function Fstarofpsi



  ! ------------------------------------------------------------------
  !					end of two-fluids free functions
  ! ------------------------------------------------------------------

  ! This function calculates the residual of the Bernoulli equation
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli(rho) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: residual

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[bernoulli]: rho = ",rho
       print *, "[bernoulli]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[bernoulli]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[bernoulli]: rho must be positive"
       stop
    end if

		residual = g_H + 0.5d0*(g_r*g_Omega)**2 &
					- gamma/(gamma-1.0d0)*g_S*rho**(gamma-1.0d0) &
					- 0.5d0/mu_mag*(g_dpsidx**2 + g_dpsidz**2)*(g_Phi/(rho*g_r))**2 &
					- 0.5d0*(g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2   &
					- g_Lambda

    return

  end function bernoulli

  ! ------------------------------------------------------------------

  ! This function calculates the residual of the Bernoulli equation
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_gauss(rho) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: residual

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[bernoulli_gauss]: rho = ",rho
       print *, "[bernoulli_gauss]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[bernoulli_gauss]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[bernoulli_gauss]: rho must be positive"
       stop
    end if

	residual = g_H + 0.5d0*(g_r*g_Omega)**2 &
				- gamma/(gamma-1.0d0)*g_S*rho**(gamma-1.0d0) &
				- 0.5d0/mu_mag*(g_dpsidx**2 + g_dpsidz**2)*(g_Phi/(rho*g_r))**2 &
				- 0.5d0*(g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2   &
				- g_Lambda  &
				- fBernmax * exp(-delta_Bern*((psi_Bern-psi_degen)/psi_max)**2)

    return

  end function bernoulli_gauss

  ! This function calculates the residual of the Bernoulli equation
  ! vphi is obtained from fluid GS to more closely mimic the MHD version
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_TF(n_den) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
	real(kind=dkind) :: residual
	real(kind=dkind) :: t_a, t_b, t_c, t_d, vphi

	if(n_den <= 0.0d0) then
       print *, "[bernoulli_TF]: n_den must be positive"
       stop
    end if

	! divide the equation by m_i to get larger numbers

	vphi = g_r/eV * g_H_i_prime - g_r*n_den**(gamma_i-1.d0)/(eV*(gamma_i-1.d0)) * g_S_i_prime *mass +  & ! CHECK MASS!
				g_Dstar_term + g_phi*g_Fstar/(n_den*g_r)

	t_a = 0.5d0 * g_gPsi2*(g_phi/g_r/n_den)**2
	t_b = 0.5d0 * vphi**2
	t_c =  gamma_i/(gamma_i-1.d0) * g_S_i * n_den**(gamma_i-1.d0) +  &
				gamma_e/(gamma_e-1.d0) * g_S_e * n_den**(gamma_e-1.d0)
	t_d = g_H / mass

!!$	residual = g_H / mass &
!!$				- gamma/(gamma-1.0d0)*g_S*n_den**(gamma-1.0d0) &
!!$				- 0.5d0*g_gPsi2*(g_Phi/(n_den*g_r))**2 &
!!$				- 0.5d0 * (ev/mass/g_R)**2*g_psi_diff**2  &
!!$				- g_Lambda / mass

	residual = t_d - t_a - t_b - t_c

    return

  end function bernoulli_TF

  ! This function calculates the residual of the Bernoulli equation
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_TF_gauss(n_den) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
	real(kind=dkind) :: residual
	real(kind=dkind) :: t_a, t_b, t_c, t_d, vphi

	if(n_den <= 0.0d0) then
       print *, "[bernoulli_TF]: n_den must be positive"
       stop
    end if

	! divide the equation by m_i to get larger numbers

	vphi = g_r/eV * g_H_i_prime - g_r*n_den**(gamma_i-1.d0)/(eV*(gamma_i-1.d0)) * g_S_i_prime *mass +  &
				g_Dstar_term + g_phi*g_Fstar/(n_den*g_r)

	t_a = 0.5d0 * g_gPsi2*(g_phi/g_r/n_den)**2
	t_b = 0.5d0 * vphi**2
	t_c =  gamma_i/(gamma_i-1.d0) * g_S_i * n_den**(gamma_i-1.d0) +  &
				gamma_e/(gamma_e-1.d0) * g_S_e * n_den**(gamma_e-1.d0)
	t_d = g_H / mass

	residual = t_d - t_a - t_b - t_c  &
				- fBernmax * exp(-delta_Bern*((psi_Bern-psi_degen)/psi_max)**2)

! Changed sign on June 21 2019
!	residual = t_d - t_a - t_b - t_c  &
!				+ fBernmax * exp(-delta_Bern*((big_Psi_Bern-big_Psi_degen)/psi_max)**2)

    return

  end function bernoulli_TF_gauss

  ! This function calculates the residual of the Bernoulli equation
  ! this version uses the Bernoulli equation as is, with delta psi for v_phi
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_TF_version1(n_den) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
	real(kind=dkind) :: residual
	real(kind=dkind) :: t_a, t_b, t_c, t_d

    if(n_den == g_Phi*g_Phi*mu_mag*mass .or. n_den == 0.0d0 ) then
       print *, "[bernoulli_TF]: n_den = ",n_den
       print *, "[bernoulli_TF]: n_den(Alfvenic) = ",g_Phi*g_Phi*mu_mag*mass
       print *, "[bernoulli_TF]: This is not yet handled"
       stop
    else if(n_den <= 0.0d0) then
       print *, "[bernoulli_TF]: n_den must be positive"
       stop
    end if

	! divide the equation by m_i to get larger numbers

	t_a = 0.5d0 * g_gPsi2*(g_phi/g_r/n_den)**2
	t_b = 0.5d0 * (ev/mass/g_r)**2 * g_psi_diff**2
	t_c =  gamma_i/(gamma_i-1.d0) * g_S_i * n_den**(gamma_i-1.d0) +  &
				gamma_e/(gamma_e-1.d0) * g_S_e * n_den**(gamma_e-1.d0)
!	t_c =  gamma/(gamma-1.d0) * g_S * n_den**(gamma-1.d0)
	t_d = g_H / mass

!!$	residual = g_H / mass &
!!$				- gamma/(gamma-1.0d0)*g_S*n_den**(gamma-1.0d0) &
!!$				- 0.5d0*g_gPsi2*(g_Phi/(n_den*g_r))**2 &
!!$				- 0.5d0 * (ev/mass/g_R)**2*g_psi_diff**2  &
!!$				- g_Lambda / mass

	residual = t_d - t_a - t_b - t_c

    return

  end function bernoulli_TF_version1

  ! This function calculates the first derivative of the Bernoulli
  ! equation with respect to rho.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dbern_drho(rho) result (db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: db

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[dbern_drho]: rho = ",rho
       print *, "[dbern_drho]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[dbern_drho]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[dbern_drho]: rho must be positive"
       stop
    end if

	db = - gamma*g_S*rho**(gamma-2.0d0) &
		  + (g_dpsidx**2 + g_dpsidz**2)*(g_Phi/(rho*g_r))**2/rho/mu_mag &
		  + (g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2/ &
		  (rho - g_Phi**2)

    return

  end function dbern_drho

  ! This function calculates the first derivative of the Bernoulli
  ! equation with respect to n_den (two-fluid version).
  ! vphi is obtained from fluid GS to more closely mimic the MHD version
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dbern_dn_den(n_den) result (db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
    real (kind=dkind) :: vphi_term, vphi_termb, db

     if(n_den <= 0.0d0) then
       print *, "[bernoulli_TF]: n_den must be positive"
       stop
    end if

	! divide the equation by m_i to get larger numbers

	vphi_term = -1.d0/n_den**3 *  &
						(g_Fstar*g_phi/g_r - g_Dstar_term*n_den + g_r/eV*g_H_i_prime*n_den - n_den**gamma_i*g_r/(eV*(gamma_i-1.d0))*g_S_i_prime*mass) *  & !CHECK MASS!!
						(g_Fstar*g_phi/g_r + n_den**gamma_i*g_r/eV*g_S_i_prime*mass) !CHECK MASS

	vphi_termb = (-2.d0*(eV*g_Fstar*g_phi + g_r**2*g_S_i_prime*mass*n_den**gamma_i)*(-(g_r**2*g_S_i_prime*mass*n_den**gamma_i) +  &
						(-1.d0 + gamma_i)*(eV*g_Fstar*g_phi + g_r*(eV*g_Dstar_term + g_H_i_prime*g_r)*n_den)))/(eV**2*(-1.d0 + gamma_i)*g_r**2*n_den**3)
	vphi_termb = vphi_termb/2.d0

	db = - gamma_i*g_S_i*n_den**(gamma_i-2.0d0) - gamma_e*g_S_e*n_den**(gamma_e-2.0d0) +  &
			g_gPsi2*(g_Phi/(n_den*g_r))**2/n_den - vphi_term

    return

  end function dbern_dn_den

  ! This function calculates the first derivative of the Bernoulli
  ! equation with respect to n_den (two-fluid version).
  ! this version uses the Bernoulli equation as is, with delta psi for v_phi
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dbern_dn_den_version1(n_den) result (db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
    real (kind=dkind) :: db

     if(n_den == g_Phi*g_Phi*mu_mag*mass .or. n_den == 0.0d0 ) then
       print *, "[bernoulli_TF]: n_den = ",n_den
       print *, "[bernoulli_TF]: n_den(Alfvenic) = ",g_Phi*g_Phi*mu_mag*mass
       print *, "[bernoulli_TF]: This is not yet handled"
       stop
    else if(n_den <= 0.0d0) then
       print *, "[bernoulli_TF]: n_den must be positive"
       stop
    end if

	! divide the equation by m_i to get larger numbers

	db = - gamma*g_S*n_den**(gamma-2.0d0) &
		  + g_gPsi2*(g_Phi/(n_den*g_r))**2/n_den

    return

  end function dbern_dn_den_version1

  ! ------------------------------------------------------------------

  ! This function calculates the residual of the Bernoulli equation
  ! as well as the first derivative of the Bernoulli equation with
  ! respect to rho for use in a Newton-Raphson solver.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine newt_rap(rho,b,db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind), intent(inout) :: b, db
    ! b is the residual of the bernoulli function
    ! and db is the derivative of this function with respect to rho.

    b = bernoulli(rho)

    db = dbern_drho(rho)

  end subroutine newt_rap


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine newt_rap_gauss(rho,b,db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind), intent(inout) :: b, db
    ! b is the residual of the bernoulli function
    ! and db is the derivative of this function with respect to rho.

    b = bernoulli_gauss(rho)

    db = dbern_drho(rho)

!	print*, rho, b, db
!	pause

  end subroutine newt_rap_gauss

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine newt_rap_TF(n_den,b,db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
    real (kind=dkind), intent(inout) :: b, db
    ! b is the residual of the bernoulli function
    ! and db is the derivative of this function with respect to n.

    b = bernoulli_TF(n_den)

    db = dbern_dn_den(n_den)

  end subroutine newt_rap_TF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine newt_rap_TF_gauss(n_den,b,db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: n_den
    real (kind=dkind), intent(inout) :: b, db
    ! b is the residual of the bernoulli function
    ! and db is the derivative of this function with respect to n.

    b = bernoulli_TF_gauss(n_den)

    db = dbern_dn_den(n_den)

  end subroutine newt_rap_TF_gauss

  ! ------------------------------------------------------------------

  ! Initialize u(1..nx , 1..nz) with an initial guess to the
  ! solution which respects the boundary conditions u = 0.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine guess_soln(u,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: u
    integer :: i,j,k
    real (kind=dkind) :: ex,ez,th,r_i,r_o,psie
	integer :: dummy_int
	real(kind=dkind), dimension(1:3) :: dummy

    ! Put some peaked initial distribution in the center
    do i=1, nx
       do j=1, nz
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
              u(i,j) = 0.0d0
          else

		  ! radius_1_3 isn't the best starting guess for free-boundary cases,
		  ! but it shouldn't matter after just a few iterations
		  if((tri_type==13).or.(bc_type==7).or.(bc_type==8)) then

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,th,rminor,  &
								dummy_int,dummy(1),dummy(2),dummy(3))
			!this is required to get rminor

			r_o = dbsval(th, r_ord, r_data(1:theta_points1+r_ord,3),  &
				theta_points1, r_cscoef(1,1:theta_points1) )

			r_i = dbsval(th, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

		  else

			call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)
			!this is required to get rminor

		  endif

		if((bc_type==7).or.(bc_type==8)) then

				 psie = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
							ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

				if(tri_type==13) then

					u(i,j) = (psic + (psie-psic)*(ex**2+ez**2)/rminor**2) *  &
							abs(r_i-sqrt(ex**2+ez**2))/max(abs(r_i-sqrt(ex**2+ez**2)),1.d-2)

				else

					u(i,j) = psic + (psie-psic)*(ex**2+ez**2)/rminor**2

				endif

			else

				u(i,j) = psic *(rminor**2 - ex*ex - ez*ez)/rminor**2

			endif

          end if
       end do
    end do

	psic = maxval(u)

	continue

  end subroutine guess_soln

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine guess_soln_TF(psi,big_Psi,psi_diff,n_den,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: psi,big_Psi,psi_diff,n_den
    integer :: i,j,k
    real (kind=dkind) :: ex,ez,th,r_i,r_o,psie
	integer :: dummy_int
	real(kind=dkind), dimension(1:3) :: dummy

	! Put some peaked initial distribution in the center
	do i=1, nx
		do j=1, nz
			! Only solve the inner region problem

			if(sort_grid(i,j)<=0) then

				psi(i,j) = 0.0d0

			else

				! radius_1_3 isn't the best starting guess for free-boundary cases,
		  		! but it shouldn't matter after just a few iterations
		  		if((tri_type==13).or.(bc_type==7).or.(bc_type==8)) then

					call radius_1_3(x_coord(i),z_coord(j),ex,ez,th,rminor,  &
									dummy_int,dummy(1),dummy(2),dummy(3))
					!this is required to get rminor

					r_o = dbsval(th, r_ord, r_data(1:theta_points1+r_ord,3),  &
						theta_points1, r_cscoef(1,1:theta_points1) )

					r_i = dbsval(th, r_ord, r_data(1:theta_points2+r_ord,6),  &
						theta_points2, r_cscoef(2,1:theta_points2) )
				
				else

					call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)
					!this is required to get rminor
				
				endif

				if((bc_type==7).or.(bc_type==8)) then

					psie = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
							ibreak_pb-psib_ord, psib_cscoef(1:ibreak_pb-psib_ord) )

					psi(i,j) = psic + (psie-psic)*(ex**2+ez**2)/rminor**2

				else

					psi(i,j) = psic *(rminor**2 - ex*ex - ez*ez)/rminor**2
				
				endif

				big_Psi(i,j) = psi(i,j)
				! No izone call for density since we should be in the inner region anyway
				n_den(i,j) = D_TF_ofpsi(psi(i,j))

			end if

		end do
	end do

	! Big psi is initialized to be same as psi, so psi_diff is just zero everywhere (for now)
	psi_diff = 0.d0

	psic = maxval(psi)
	big_Psic = maxval(big_Psi)
	psi_max = max(psic,big_Psic)

	continue

  end subroutine guess_soln_TF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine two_fluid_initial(psi,rho, n_den, psi_diff, big_Psi,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: nx, nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi, rho
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: n_den, psi_diff, big_Psi
	real(kind=dkind) :: vphi_loc, phic, omegac, x
	integer :: i, j

	do j = 1, nz
		do i = 1, nx

			! Added bc_type exclusion (NEW May 2025)
			if(sort_grid(i,j)<0.and.bc_type/=7) then
				cycle
			endif

			x = x_coord(i)
			! Would need to add i_zone calls here if 1-fluid run is ever done with free-boundaries implemented
			phic = phiofpsi(psi(i,j))
			omegac = omegaofpsi(psi(i,j))

			vphi_loc = (phic*iofpsi(psi(i,j))/rho(i,j)/x + x*omegac)/ &
							(1.0d0 - phic**2/rho(i,j))

			psi_diff(i,j) = mass/eV * x*vphi_loc
			big_Psi(i,j) = psi(i,j) + psi_diff(i,j)

			n_den(i,j) = rho(i,j)/mass

		enddo
	enddo

end subroutine two_fluid_initial

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine initialize_density(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: rho
    integer :: i,j

    do i=1, nx
       do j=1, nz
          rho(i,j) = dofpsi(psi(i,j))
       end do
    end do

	call bc_psi_rho0(psi,rho,nx,nz)

	continue

  end subroutine initialize_density

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine initialize_b(b,psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    use constant, only: x_coord, dx, dz
	implicit none
	integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi,rho
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: b
	real (kind=dkind) :: r,ex,ez
    integer :: i,j,k

	b(1:nx,1:nz) = 0.d0
	return

  end subroutine initialize_b

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine update_b(psi,rho,b,nx,nz,orp,anorm)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    use constant
	use triangularity
	implicit none
	integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi,rho
	real (kind=dkind), intent(in) :: orp
	real (kind=dkind), intent(inout) :: anorm
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: b
	real (kind=dkind) :: r,ex,ez,res
	real (kind=dkind), dimension(1:3,1:3) :: psi3x3
    integer :: i,j,iii,jjj,k

	continue
	return

  end subroutine update_b

  ! Test for symmetry in z
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine symmetric(q,nx,nz,name)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: q
    character (len=*), intent(in) :: name
    integer :: riz,liz
    integer :: ir
    real (kind=dkind) :: delta, mx

    mx = 0.0d0 ! Initialize the max variable

    do liz=1, nz
       riz = nz + 1 - liz
       if(riz <= liz) exit

       do ir=1, nx
          delta = dabs(q(ir,liz) - q(ir,riz))
          mx = dmax1(mx,delta)
       end do
    end do

    if(mx > 0.0d0) then
       print *, name," is NOT symmetric, Max | Difference | = ",mx
    else
       print *, name," is symmetric"
    end if

  end subroutine symmetric


  ! van Leer Slope Limiter
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function van_Leer_slope(l,c,r,dx) result (slope)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: l,c,r,dx
    real (kind=dkind) :: slope, up, down

    up = (r-c)/dx
    down = (c-l)/dx
    if((up==0.0d0) .or. (down==0.0)) then
       slope = 0.0d0
    else
       slope = (up*dabs(down) + down*dabs(up))/(dabs(up)+dabs(down))
    endif

    return
  end function van_Leer_slope

  ! ------------------------------------------------------------------

  ! van Leer Slope Limiter
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function van_Leer_slope_new(l,c,r,dxl,dxr) result (slope)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: l,c,r,dxl,dxr
    real (kind=dkind) :: slope, up, down

    up = (r-c)/dxr
    down = (c-l)/dxl
    if((up==0.0d0) .or. (down==0.0)) then
       slope = 0.0d0
    else
       slope = (up*dabs(down) + down*dabs(up))/(dabs(up)+dabs(down))
    endif

    return
  end function van_Leer_slope_new



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
	real (kind=dkind) term1
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
    real (kind=dkind) :: phic

	real (kind=dkind), dimension (1:3,1:3) :: psi3x3
	integer iii,jjj,k

	real(kind=skind), dimension(1:nx,1:nz) :: root_diff

    ! -----------------------------------------------------


		if(Broot<4) then
			mtm_limit = 1.d0
		else
			mtm_limit = 1.d2
		endif

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

!		  rho(i,j) = dofpsi(psi(i,j))
!		  cycle

		!  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
		!		.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
		!		.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
		!  Ian comment: don't remember why, but commented out above if's and replaced with below
		  if ( (tri_type==13).and.(sort_grid(i,j)<1) )  then
				if(gravity_type==0) then
					rho(i,j) = dofpsi(0.d0)
				else
					rho(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if


          !  "Setup Global Parameters"
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))

			g_D = dofpsi(psi(i,j))		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))



		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho

				rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(phic == 0.0d0) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)

!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(rhomin,rhomax,1000)
!!!					pause
!!!					rho(i,j) = dofpsi(psi(i,j))
!!!					cycle

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   ! Ian comment: if statement below is new; not sure if necessary or correct 
				   if (mach_theta_max<5.d-2*mach_theta_max_initial) then
						! just patch the point and hope for the best
!						patch_value(i,j) = .true.
						rho(i,j) = g_D
						goto 101
				   endif

				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))


                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)


!!$					   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

!!$		  root_diff(i,j) = (heavy-light)/(heavy + light)

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

!!$	call radial_plot(root_diff,psi,nx,nz,"root_diff_SF",nx/2)

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz))
	   g_I = iofpsi(psi(min_ix,min_iz))
	   g_D = dofpsi(psi(min_ix,min_iz))
       ! Calculate the derivatives of psi
!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )

       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(i),z_coord(j))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

!!$			   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

  end subroutine update_rho


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_gauss(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    real (kind=dkind) :: light_gauss, heavy_gauss ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
!	external rtsafe, rtbis
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit 
    real (kind=dkind) :: phic
	real (kind=skind), dimension(nx,nz) :: out

	integer iii,jjj,k

    ! -----------------------------------------------------

!	print*, 'in Bernoulli Gauss'

	mtm_limit = 1.d0

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j)>0)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo



!  if(allocated(fmax_2D)) then
!  	if(size(fmax_2D,1)==nx) then
!  		continue
!  	else
!	  	deallocate(fmax_2D)
!  	endif
!  endif

  if(allocated(fmax_2D)) then
	continue
  else
  	allocate(fmax_2D(1:nx,1:nz))
  	fmax_2D = 0.d0
  endif



!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

          !  "Setup Global Parameters"
		  psi_Bern = psi(i,j)
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))

			g_D = dofpsi(psi(i,j))		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))


		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho


			rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
				     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if((phic == 0.0d0).AND.(eq_type==1)) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  elseif(eq_type==3) then
             ! No poloidal flow allowed:
			 ! there is ALWAYS only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax
		  rhomax = rhomax*3.d1

          ! Calculate a minimum density
		   if(eq_type==1) then
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r	
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)
		   else
				rhomin = 1.d-32
		   endif

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif

!			print*, i,j,rhomin, rhomax
!			print*, '      '

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))
						rhomax = rhomax * 3.d1

                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r 
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm 
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

		! find Bernmax (maximum of Bernoulli function between the roots)
        Bernmax = rtbis(dbern_drho,light,heavy,1.0d-10)
		fBernmax = bernoulli(Bernmax)
		fmax_2D(i,j) = fBernmax

		! now solve Bernoulli again, using the modified Bernoulli
!          light_gauss = rtsafe(newt_rap_gauss,xb1(1),xb2(1),1.0d-14,10000)
!          heavy_gauss = rtsafe(newt_rap_gauss,xb1(2),xb2(2),1.0d-14,10000)
          light = rtsafe(newt_rap_gauss,rhomin,Bernmax,1.0d-14,10000)
          light_gauss = light
          heavy = rtsafe(newt_rap_gauss,Bernmax,rhomax,1.0d-14,10000)
          heavy_gauss = heavy


          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region 
             ! around the high density (slowly moving) core 
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy_gauss
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light_gauss
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy_gauss
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light_gauss
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

!			print*, i,j, rho(i,j), heavy_gauss, light_gauss

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our 
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz))		
	   g_I = iofpsi(psi(min_ix,min_iz))
	   g_D = dofpsi(psi(min_ix,min_iz))
       ! Calculate the derivatives of psi

!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )


       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(i),z_coord(j))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm 
		  endif 
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if


!	out(1:nx,1:nz) = fmax_2D(1:nx,1:nz)
!	call radial_plot(out,psi,nx,nz,"fmax",nz/2)

!	pause


	continue

    ! -----------------------------------------------------

  end subroutine update_rho_gauss

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density
! Jul 24 2013: add a patch for values that do not converge (should be necessary only at 
! the beginning of each grid)

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex, ez, dx2, dz2
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	integer :: i_zone = 0
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.5d-2
	endif

! Ian comment: double check that commented out code below isn't necessary,
! and that this operation gets performed BEFORE update_rho_TF gets called
!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    do j=1, nz
       do i=1, nx
			n_old(i,j) = n_den(i,j)
		enddo
	enddo

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

		if((i==i_issue_check).and.(j==j_issue_check)) then
			continue
		endif

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		! March 18 2021: bc_type==7 is dealt with in a different way using sort_grid
		  if( (tri_type==13).and.(sort_grid(i,j)==1)) then
!!$		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
!!$				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
!!$				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(i,j)-2)
				i_zone = max(-1,i_zone)
			endif		


		  if(patch_value(i,j)) then
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j),i_zone)
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		  endif

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j),i_zone)
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j),i_zone)
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j),i_zone)
          g_S_i = s_i_ofpsi(big_Psi(i,j),i_zone)
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j),i_zone)
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j),i_zone)
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)
		  g_H_i_prime = dh_i_dpsi(big_Psi(i,j),i_zone)
		  g_S_i_prime = dS_i_dpsi(big_Psi(i,j),i_zone)
		  g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j),i_zone)

			call get_Delstar_piece(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j),i_zone)

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min/10.d0**6,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(((nb == 1).and.(g_Phi**2>0.d0)) .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm

				   if (mach_theta_max<5.d-2*mach_theta_max_initial) then
						! just patch the point and hope for the best
						patch_value(i,j) = .true.
						n_den(i,j) = g_D
						goto 101
				   endif

				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif

                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j),i_zone)
				   g_H_e = h_e_ofpsi_partial(psi(i,j),i_zone)
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j),i_zone)
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
!						n_den_min = dmax1(n_den_min,1.0d-31)
						n_den_min = dmax1(n_den_min*10.d0**(-10),1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
!						dmtm = 10.0d0*dmtm
						dmtm = 2.0d0*dmtm
						! June 15 2019: changed this to 2
				   endif
				   do while(mtm_soln-dmtm<=0.d0)
						dmtm = dmtm*(1.0-1.d-3)
					enddo
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
!             if(psi(i,j) >= psi_degen) then
             if(big_Psi(i,j) >= big_Psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  if(write_TF_roots) then
			  root_diff(i,j) = (heavy-light)/(heavy + light)
			  heavy_2D(i,j) = heavy
			  light_2D(i,j) = light
		  endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)
    big_Psi_degen = big_Psi(min_ix,min_iz)

    ! -----------------------------------------------------

	if(write_TF_roots) then
		call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
		call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
		call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)
	endif

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz),i_zone)
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz),i_zone)
       g_S = g_S_e +  g_S_i
	  g_S_i_prime = dS_i_dpsi(big_Psi(min_ix,min_iz),i_zone)
	  g_Fstar = Fstarofpsi(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),i_zone)
	  g_psi_diff = psi_diff(min_ix,min_iz)

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz),i_zone)
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz),i_zone)
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz),i_zone)
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
			g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz),i_zone)

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz),i_zone)
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz),i_zone)
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz),i_zone)
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz),i_zone)

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
!!$	real(kind=dkind) :: n_bound_max_old
	real(kind=dkind) :: t_b, t_c, t_d
!!$	real(kind=dkind) :: t_b_old, t_c_old, t_d_old

	if(abs(g_phi)==0.d0) then
	! only one solution, easy to bracket
		n_bound_max = g_D * 5.d0
		return
	endif
	
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_c_old = mass * gamma/(gamma-1.d0) * g_S
!!$	t_d_old = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_i,gamma_e)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_max_old = ((t_d_old+t_b_old)/t_c_old)**(1.d0/(gamma-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_c (old, new) = ', t_c_old, t_c
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_max (old, new) = ', n_bound_max_old, n_bound_max
!!$	pause

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	


	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
!!$	real(kind=dkind) :: n_bound_min_old
	real(kind=dkind) :: t_a, t_b, t_d
!!$	real(kind=dkind) :: t_a_old, t_b_old, t_d_old

!!$	t_a_old = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_d_old = g_H

	if(abs(g_phi)==0.d0) then
	! only one solution, easy to bracket
		n_bound_min = 0.d0
		return
	endif

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_min_old = sqrt(t_a_old/(t_d_old + t_b_old))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_a (old, new) = ', t_a_old, t_a
!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_min (old, new) = ', n_bound_min_old, n_bound_min
!!$	pause


	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max_true(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! gammas not fixed, since this routine is never called

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))

	continue

end subroutine get_n_bound_max_true

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min_true(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_min = sqrt(t_a/(t_d-t_b))

	continue

end subroutine get_n_bound_min_true

 end subroutine update_rho_TF

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_frozen_March_18_2021(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! March 18 2021: We are updating the Bernoulli solver to use sort_grid instead of local calculations. Freeze the previous version
! for safety and possibly allow the use of either.
! The two-fluid routine updates the number density, not the mass density
! Jul 24 2013: add a patch for values that do not converge (should be necessary only at 
! the beginning of each grid)

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex, ez, dx2, dz2
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.5d-2
	endif

!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    do j=1, nz
       do i=1, nx
			n_old(i,j) = n_den(i,j)
		enddo
	enddo

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

		if((i==i_issue_check).and.(j==j_issue_check)) then
			continue
		endif

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

		  if(patch_value(i,j)) then
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		  endif

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j))
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j))
          g_S_i = s_i_ofpsi(big_Psi(i,j))
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j))
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)
		  g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
		  g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
		  g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

			call get_Delstar_piece(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j))

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min/10.d0**6,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(((nb == 1).and.(g_Phi**2>0.d0)) .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm

				   if (mach_theta_max<5.d-2*mach_theta_max_initial) then
						! just patch the point and hope for the best
						patch_value(i,j) = .true.
						n_den(i,j) = g_D
						goto 101
				   endif

				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif

                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j))
				   g_H_e = h_e_ofpsi_partial(psi(i,j))
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
!						n_den_min = dmax1(n_den_min,1.0d-31)
						n_den_min = dmax1(n_den_min*10.d0**(-10),1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
!						dmtm = 10.0d0*dmtm
						dmtm = 2.0d0*dmtm
						! June 15 2019: changed this to 2
				   endif
				   do while(mtm_soln-dmtm<=0.d0)
						dmtm = dmtm*(1.0-1.d-3)
					enddo
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
!             if(psi(i,j) >= psi_degen) then
             if(big_Psi(i,j) >= big_Psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  if(write_TF_roots) then
			  root_diff(i,j) = (heavy-light)/(heavy + light)
			  heavy_2D(i,j) = heavy
			  light_2D(i,j) = light
		  endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)
    big_Psi_degen = big_Psi(min_ix,min_iz)

    ! -----------------------------------------------------

	if(write_TF_roots) then
		call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
		call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
		call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)
	endif

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz))
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz))
       g_S = g_S_e +  g_S_i
	  g_S_i_prime = dS_i_dpsi(big_Psi(min_ix,min_iz))
	  g_Fstar = Fstarofpsi(psi(min_ix,min_iz),big_Psi(min_ix,min_iz))
	  g_psi_diff = psi_diff(min_ix,min_iz)

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
			g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
!!$	real(kind=dkind) :: n_bound_max_old
	real(kind=dkind) :: t_b, t_c, t_d
!!$	real(kind=dkind) :: t_b_old, t_c_old, t_d_old

	if(abs(g_phi)==0.d0) then
	! only one solution, easy to bracket
		n_bound_max = g_D * 5.d0
		return
	endif
	
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_c_old = mass * gamma/(gamma-1.d0) * g_S
!!$	t_d_old = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_i,gamma_e)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_max_old = ((t_d_old+t_b_old)/t_c_old)**(1.d0/(gamma-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_c (old, new) = ', t_c_old, t_c
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_max (old, new) = ', n_bound_max_old, n_bound_max
!!$	pause

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	


	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
!!$	real(kind=dkind) :: n_bound_min_old
	real(kind=dkind) :: t_a, t_b, t_d
!!$	real(kind=dkind) :: t_a_old, t_b_old, t_d_old

!!$	t_a_old = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_d_old = g_H

	if(abs(g_phi)==0.d0) then
	! only one solution, easy to bracket
		n_bound_min = 0.d0
		return
	endif

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_min_old = sqrt(t_a_old/(t_d_old + t_b_old))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_a (old, new) = ', t_a_old, t_a
!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_min (old, new) = ', n_bound_min_old, n_bound_min
!!$	pause


	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max_true(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! gammas not fixed, since this routine is never called

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))

	continue

end subroutine get_n_bound_max_true

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min_true(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_min = sqrt(t_a/(t_d-t_b))

	continue

end subroutine get_n_bound_min_true

 end subroutine update_rho_TF_frozen_March_18_2021

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_gauss(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density
! Jul 24 2013: add a patch for values that do not converge (should be necessary only at 
! the beginning of each grid)

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    real (kind=dkind) :: light_gauss, heavy_gauss ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex, ez, dx2, dz2
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz

	integer iii,jjj,k

    ! -----------------------------------------------------

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.5d-2
	endif

!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

  if(allocated(fmax_2D)) then
	continue
  else
  	allocate(fmax_2D(1:nx,1:nz))
  	fmax_2D = 0.d0
  endif

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    do j=1, nz
       do i=1, nx
			n_old(i,j) = n_den(i,j)
		enddo
	enddo

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

		  if(patch_value(i,j)) then
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		  endif

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  psi_Bern = psi(i,j)
		  big_Psi_Bern = big_Psi(i,j)
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j))
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j))
          g_S_i = s_i_ofpsi(big_Psi(i,j))
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j))
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)
		  g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
		  g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
		  g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

			call get_Delstar_piece(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j))

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm

				   if (mach_theta_max<5.d-2*mach_theta_max_initial) then
						! just patch the point and hope for the best
						patch_value(i,j) = .true.
						n_den(i,j) = g_D
						goto 101
				   endif

				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif

                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j))
				   g_H_e = h_e_ofpsi_partial(psi(i,j))
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
						n_den_min = dmax1(n_den_min,1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm
				   endif
				   do while(mtm_soln-dmtm<=0.d0)
						dmtm = dmtm*(1.0-1.d-3)
					enddo
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

		! find Bernmax (maximum of Bernoulli function between the roots)
        Bernmax = rtbis(dbern_dn_den,light,heavy,1.0d-10)
		fBernmax = bernoulli_TF(Bernmax)
		fmax_2D(i,j) = fBernmax

		! now solve Bernoulli again, using the modified Bernoulli
!          light_gauss = rtsafe(newt_rap_gauss,xb1(1),xb2(1),1.0d-14,10000)
!          heavy_gauss = rtsafe(newt_rap_gauss,xb1(2),xb2(2),1.0d-14,10000)
          light_gauss = rtsafe(newt_rap_TF_gauss,n_den_min,Bernmax,x_TF_acc,10000)
!         light_gauss = light
          heavy_gauss = rtsafe(newt_rap_TF_gauss,Bernmax,n_den_max,x_TF_acc,10000)
!          heavy_gauss = heavy

			if(light_gauss==0.d0) then
				light_gauss = light
			endif
			if(heavy_gauss==0.d0) then
				heavy_gauss = heavy
			endif

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
!             if(psi(i,j) >= psi_degen) then
             if(big_Psi(i,j) >= big_Psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy_gauss
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light_gauss
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy_gauss
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light_gauss
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  if(write_TF_roots) then
			  root_diff(i,j) = (heavy-light)/(heavy + light)
			  heavy_2D(i,j) = heavy
			  light_2D(i,j) = light
		  endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)
    big_Psi_degen = big_Psi(min_ix,min_iz)

    ! -----------------------------------------------------

	if(write_TF_roots) then
		call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
		call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
		call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)
	endif

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz))
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz))
       g_S = g_S_e +  g_S_i
	  g_S_i_prime = dS_i_dpsi(big_Psi(min_ix,min_iz))
	  g_Fstar = Fstarofpsi(psi(min_ix,min_iz),big_Psi(min_ix,min_iz))
	  g_psi_diff = psi_diff(min_ix,min_iz)

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
!!$	real(kind=dkind) :: n_bound_max_old
	real(kind=dkind) :: t_b, t_c, t_d
!!$	real(kind=dkind) :: t_b_old, t_c_old, t_d_old

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_c_old = mass * gamma/(gamma-1.d0) * g_S
!!$	t_d_old = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_i,gamma_e)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_max_old = ((t_d_old+t_b_old)/t_c_old)**(1.d0/(gamma-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_c (old, new) = ', t_c_old, t_c
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_max (old, new) = ', n_bound_max_old, n_bound_max
!!$	pause

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
!!$	real(kind=dkind) :: n_bound_min_old
	real(kind=dkind) :: t_a, t_b, t_d
!!$	real(kind=dkind) :: t_a_old, t_b_old, t_d_old

!!$	t_a_old = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_d_old = g_H

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_min_old = sqrt(t_a_old/(t_d_old + t_b_old))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_a (old, new) = ', t_a_old, t_a
!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_min (old, new) = ', n_bound_min_old, n_bound_min
!!$	pause


	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

 end subroutine update_rho_TF_gauss

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_gauss_2020_flow(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density
! Jul 24 2013: add a patch for values that do not converge (should be necessary only at 
! the beginning of each grid)
! 2 2020
! This version is called only for seek_mtm=1
! First we find the minimum of the Bernoulli function, then adjust mach_theta_max and finally solve for the density.

	use constant, only : dx, dz, dx_a, dz_a

	integer, intent(in) :: nx,nz,seek_mtm
	real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
	real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
	real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
	! mtm_acc controls the bisection loop seeking mach theta max
	real (kind=dkind), intent(in) :: mtm_acc
	integer :: i,j
	real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
	real (kind=dkind) :: light_gauss, heavy_gauss ! Temp. for Bernoulli roots
	! minimum delta n_den for bernoulli roots, the mean val. and location
	real (kind=dkind) :: mean_n_den, tmp
	real (kind=dkind), intent(inout) :: min_dn_den
	real(kind=dkind) :: min_Bern, Bern_ave
	real(kind=dkind) :: mtm_step(100), min_step(100)
	integer, intent(inout) :: min_ix, min_iz
	! We use n_den_ext to store n_den for which Bernoulli function is a maximum
	real (kind=dkind) :: n_den_ext ! n_den extremum
	integer :: repeat ! Logical for repeating the density calc.
	real (kind=dkind) :: ex, ez, dx2, dz2
	integer :: nb
	real (kind=dkind), dimension(1:10) :: xb1,xb2
	real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
	! The next set of variables are used for the bisection search for
	! locating the degenerate root for Mach Theta Max
	real (kind=dkind) :: dmtm, mtm_soln
	! Maximum iteration loop for Bisection search
	integer, parameter :: mtm_jmax = 100
	real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	real(kind=dkind), dimension(1:nx,1:nz)  :: n_den_min_2D, n_den_max_2D, Bernmax_2D

	real(kind=dkind) :: min_tol=10.d-9 ! default value: we will have to do trial and error on this value
!	real(kind=dkind) :: min_tol=10.d-12 !will have to do trial and error on this value
	real(kind=dkind) :: min_err, fprim
	integer :: min_ite
	logical :: min_search
	integer, save :: calls = 0
	integer iii,jjj,k

! -----------------------------------------------------

	calls = calls+1

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
	mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.5d-2
	endif

	! This may need some improvment
	if(allocated(fmax_2D)) then
		continue
	else
		allocate(fmax_2D(1:nx,1:nz))
		fmax_2D = 0.d0
	endif

	! We need the previous iteration of the density to solve for the new one
	do j=1, nz
	do i=1, nx
		n_old(i,j) = n_den(i,j)
	enddo
	enddo


! -----------------------------------------------------


	! Mach_theta_search part

	min_search = .true.
	min_ite = 0

	do while(min_search)

		min_ite = min_ite + 1
		mtm_step(min_ite) = mach_theta_max

		min_Bern = 10.d99
		Bern_ave = 0.d0
		min_ix = 1
		min_iz = 1

		do j=1, nz
		do i=1, nx
		! Only solve the inner region problem

			if(sort_grid(i,j)<=0) then
				cycle
			end if


			if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
			if(gravity_type==0) then
				n_den(i,j) = dofpsi(0.d0)
			else
				n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
					(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
					/sofpsi(0.d0)  &
					)**(1.d0/(gamma-1.d0))
			endif
				cycle
			end if

			if(patch_value(i,j)) then
				! We can get rid of this.
				n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
				print*, 'patching density in i = ', i, ' j = ', j
				cycle
			endif

			!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
			!-----------------------------------------

			!  "Setup Global Parameters"
			! Set the globals for solving Bernouli eqn.
			psi_Bern = psi(i,j)
			big_Psi_Bern = big_Psi(i,j)
			g_r = x_coord(i)
			g_Phi = phi_TF_ofpsi(big_Psi(i,j))
			! print *, "g_S"
			g_S_e = s_e_ofpsi(psi(i,j))
			g_S_i = s_i_ofpsi(big_Psi(i,j))
			g_S = g_S_e +  g_S_i
			! print *, "g_H"
			g_H_e = h_e_ofpsi_partial(psi(i,j))
			g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
			call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
			g_psi_diff = psi_diff(i,j)
			g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
			g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
			g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

			call get_Delstar_piece(i,j)

			g_Lambda = grav_potential(x_coord(i),z_coord(j))

			g_mtheta=mach_theta(psi(i,j))

			g_indi=i
			g_indj=j

			! Calculate the derivatives of big_Psi (derivative of psi is not needed)

			call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

			g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

			! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

			if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
			endif

			! Calculate a minimum density
			call get_n_bound_min(n_den_min)
			if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
			! "Only a sub-Alfvenic Root!"
			else
				n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
			end if
			if(n_den_min>n_den_max) then
				n_den_min = 1.001d0*g_phi**2*mass*mu_mag
			endif

			if(n_den_max<n_den_min) then
				continue
			endif

			!Save n_den_min and n_den_max, since we'll use them again
			n_den_min_2D(i,j) = n_den_min
			n_den_max_2D(i,j) = n_den_max

			n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
			Bernmax_2D(i,j) = n_den_ext
			! Next we evaluate the Bernoulli function
			tmp = bernoulli_TF(n_den_ext)
			fmax_2D(i,j) = tmp
			Bern_ave = Bern_ave + tmp

			if(tmp<min_Bern) then
				min_Bern = tmp
				min_ix = i
				min_iz = j
			endif

		enddo
		enddo

		min_step(min_ite) = min_Bern
		Bern_ave = Bern_ave/(nx*nz)
		min_err = abs(min_Bern/Bern_ave)

		if((min_err<min_tol).and.(min_Bern>=0.d0)) then
			!convergence
			min_search = .false.
		else
			!update mach_theta_max
			if(min_ite==1) then
			! distinguish the first update
				if(min_Bern>0.d0) then
					! we need to increase mach_theta_max
					mach_theta_max = mach_theta_max*(1.d0+min(1.d-2,min_err))
				else
					! we need to decrease mach_theta_max
					mach_theta_max = mach_theta_max*(1.d0-min(1.d-2,min_err))
				endif
			else
			! we can proceed with the secant part of the mach_thet_max search
				fprim = (min_step(min_ite)-min_step(min_ite-1))/  &
					(mtm_step(min_ite)-mtm_step(min_ite-1))
				if(fprim==0.d0) then
					print*, "Warning: fprim=0 in update_rho_TF_Gauss_2020_flow"
					if(min_Bern>0.d0) then
						! we need to increase mach_theta_max
						mach_theta_max = mach_theta_max*(1.d0+min(1.d-2,min_err))
					else
						! we need to decrease mach_theta_max
						mach_theta_max = mach_theta_max*(1.d0-min(1.d-2,min_err))
					endif
				else
					mach_theta_max = mach_theta_max - min_step(min_ite)/fprim
				endif
			endif

			! Last, make sure that mach_theta_max is not changing too much (it can happen in the first few iterations).
			if(max(mtm_step(min_ite)/mach_theta_max,mach_theta_max/mtm_step(min_ite))>1.1d0) then
				if(mach_theta_max>mtm_step(min_ite)) then
					mach_theta_max = mtm_step(min_ite) * (1.d0+1.d-2)
				else
					mach_theta_max = mtm_step(min_ite) * (1.d0-1.d-2)
				endif
			endif

		endif

	enddo

	!---------------------------------------------------------------------------------
	! Now we have mach_theta_max and can proceed to calculate the roots of the Bernoulli equation
	!---------------------------------------------------------------------------------

!	print *, min_ix, min_iz, x_coord(min_ix), z_coord(min_iz), min_Bern
	write(id_Bern,1818)  calls, min_ix, min_iz, x_coord(min_ix), z_coord(min_iz), min_Bern

1818   format(3(I5,3x), 3(e14.8, 3x))

	if(min_Bern>1.d-2) then
		continue
	endif

	do j=1, nz
	do i=1, nx
	! Only solve the inner region problem

		if(sort_grid(i,j)<=0) then
			cycle
		end if

		! Could get rid of this, but we can leave it for now to keep all bases covered.
		if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
			.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
			.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
		if(gravity_type==0) then
			n_den(i,j) = dofpsi(0.d0)
		else
			n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
				(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
				/sofpsi(0.d0)  &
				)**(1.d0/(gamma-1.d0))
		endif
			cycle
		end if

		if(patch_value(i,j)) then
			! This can be removed
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		endif

		!  "Setup Global Parameters"
		! Set the globals for solving Bernouli equation.
		! It is not worth saving all these values from the previous part.
		psi_Bern = psi(i,j)
		big_Psi_Bern = big_Psi(i,j)
		g_r = x_coord(i)
		g_Phi = phi_TF_ofpsi(big_Psi(i,j))
		g_D = D_TF_ofpsi(big_Psi(i,j))
		! print *, "g_S"
		g_S_e = s_e_ofpsi(psi(i,j))
		g_S_i = s_i_ofpsi(big_Psi(i,j))
		g_S = g_S_e +  g_S_i
		! print *, "g_H"
		g_H_e = h_e_ofpsi_partial(psi(i,j))
		g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
		call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		g_psi_diff = psi_diff(i,j)
		g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
		g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
		g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

		call get_Delstar_piece(i,j)

		g_Lambda = grav_potential(x_coord(i),z_coord(j))

		g_mtheta=mach_theta(psi(i,j))

		g_indi=i
		g_indj=j

		! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

		! Now solve the modified Bernoulli equation using all the information from the previous part.

		Bernmax = Bernmax_2D(i,j)
		fBernmax = fmax_2D(i,j)
		n_den_min = n_den_min_2D(i,j)
		n_den_max = n_den_max_2D(i,j)

		if((Broot == 0).or.(Broot == 5)) then
			! Put a low density (rapidly moving) region
			! around the high density (slowly moving) core
			!             if(psi(i,j) >= psi_degen) then
			if(big_Psi(i,j) >= big_Psi_degen) then
			! Choose the last bracket (highest density)
				heavy_gauss = rtsafe(newt_rap_TF_gauss,Bernmax,n_den_max,x_TF_acc,10000)
				n_den(i,j) = heavy_gauss
			else
			! Choose the first bracket (lowest density)
				light_gauss = rtsafe(newt_rap_TF_gauss,n_den_min,Bernmax,x_TF_acc,10000)
				n_den(i,j) = light_gauss
			end if
		elseif((Broot == 1).or.(Broot == 4)) then
		! Choose the heavy root
			heavy_gauss = rtsafe(newt_rap_TF_gauss,Bernmax,n_den_max,x_TF_acc,10000)
			n_den(i,j) = heavy_gauss
		else if(Broot == 2) then
		! choose the light root
			light_gauss = rtsafe(newt_rap_TF_gauss,n_den_min,Bernmax,x_TF_acc,10000)
			n_den(i,j) = light_gauss
		else
			! The code will never get here, unless this routine is called for all grids and no single-fluid
			! Bernoulli solver is called.
			print *, "I don't understand Broot =",Broot
			pause
			stop
		endif

	enddo
	enddo
 
 	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
!!$	real(kind=dkind) :: n_bound_max_old
	real(kind=dkind) :: t_b, t_c, t_d
!!$	real(kind=dkind) :: t_b_old, t_c_old, t_d_old

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_c_old = mass * gamma/(gamma-1.d0) * g_S
!!$	t_d_old = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_i,gamma_e)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_max_old = ((t_d_old+t_b_old)/t_c_old)**(1.d0/(gamma-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_c (old, new) = ', t_c_old, t_c
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_max (old, new) = ', n_bound_max_old, n_bound_max
!!$	pause

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
!!$	real(kind=dkind) :: n_bound_min_old
	real(kind=dkind) :: t_a, t_b, t_d
!!$	real(kind=dkind) :: t_a_old, t_b_old, t_d_old

!!$	t_a_old = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_d_old = g_H

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_min_old = sqrt(t_a_old/(t_d_old + t_b_old))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_a (old, new) = ', t_a_old, t_a
!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_min (old, new) = ', n_bound_min_old, n_bound_min
!!$	pause


	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

end subroutine update_rho_TF_gauss_2020_flow

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_gauss_2020_magnetic(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density
! Jul 24 2013: add a patch for values that do not converge (should be necessary only at 
! the beginning of each grid)
! 2 2020
! This version is called only for seek_mtm=1
! First we find the minimum of the Bernoulli function, then adjust mach_theta_max and finally solve for the density.
! This version uses a magnetic surface to divide the sub- and supersonic regions.

	use constant, only : dx, dz, dx_a, dz_a

	integer, intent(in) :: nx,nz,seek_mtm
	real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
	real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
	real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
	! mtm_acc controls the bisection loop seeking mach theta max
	real (kind=dkind), intent(in) :: mtm_acc
	integer :: i,j
	real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
	real (kind=dkind) :: light_gauss, heavy_gauss ! Temp. for Bernoulli roots
	! minimum delta n_den for bernoulli roots, the mean val. and location
	real (kind=dkind) :: mean_n_den, tmp
	real (kind=dkind), intent(inout) :: min_dn_den
	real(kind=dkind) :: min_Bern, Bern_ave
	real(kind=dkind) :: mtm_step(100), min_step(100)
	integer, intent(inout) :: min_ix, min_iz
	! We use n_den_ext to store n_den for which Bernoulli function is a maximum
	real (kind=dkind) :: n_den_ext ! n_den extremum
	integer :: repeat ! Logical for repeating the density calc.
	real (kind=dkind) :: ex, ez, dx2, dz2
	integer :: nb
	real (kind=dkind), dimension(1:10) :: xb1,xb2
	real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
	! The next set of variables are used for the bisection search for
	! locating the degenerate root for Mach Theta Max
	real (kind=dkind) :: dmtm, mtm_soln
	! Maximum iteration loop for Bisection search
	integer, parameter :: mtm_jmax = 100
	real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	real(kind=dkind), dimension(1:nx,1:nz)  :: n_den_min_2D, n_den_max_2D, Bernmax_2D

	real(kind=dkind) :: min_tol=10.d-9 ! default value: we will have to do trial and error on this value
!	real(kind=dkind) :: min_tol=10.d-12 !will have to do trial and error on this value
	real(kind=dkind) :: min_err, fprim
	integer :: min_ite
	logical :: min_search
	integer, save :: calls = 0
	integer iii,jjj,k

! -----------------------------------------------------

	calls = calls+1

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
	mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.5d-2
	endif

	! This may need some improvment
	if(allocated(fmax_2D)) then
		continue
	else
		allocate(fmax_2D(1:nx,1:nz))
		fmax_2D = 0.d0
	endif

	! We need the previous iteration of the density to solve for the new one
	do j=1, nz
	do i=1, nx
		n_old(i,j) = n_den(i,j)
	enddo
	enddo


! -----------------------------------------------------


	! Mach_theta_search part

	min_search = .true.
	min_ite = 0

	do while(min_search)

		min_ite = min_ite + 1
		mtm_step(min_ite) = mach_theta_max

		min_Bern = 10.d99
		Bern_ave = 0.d0
		min_ix = 1
		min_iz = 1

		do j=1, nz
		do i=1, nx
		! Only solve the inner region problem

			if(sort_grid(i,j)<=0) then
				cycle
			end if


			if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
			if(gravity_type==0) then
				n_den(i,j) = dofpsi(0.d0)
			else
				n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
					(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
					/sofpsi(0.d0)  &
					)**(1.d0/(gamma-1.d0))
			endif
				cycle
			end if

			if(patch_value(i,j)) then
				! We can get rid of this.
				n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
				print*, 'patching density in i = ', i, ' j = ', j
				cycle
			endif

			!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
			!-----------------------------------------

			!  "Setup Global Parameters"
			! Set the globals for solving Bernouli eqn.
			psi_Bern = psi(i,j)
			big_Psi_Bern = big_Psi(i,j)
			g_r = x_coord(i)
			g_Phi = phi_TF_ofpsi(big_Psi(i,j))
			! print *, "g_S"
			g_S_e = s_e_ofpsi(psi(i,j))
			g_S_i = s_i_ofpsi(big_Psi(i,j))
			g_S = g_S_e +  g_S_i
			! print *, "g_H"
			g_H_e = h_e_ofpsi_partial(psi(i,j))
			g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
			call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
			g_psi_diff = psi_diff(i,j)
			g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
			g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
			g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

			call get_Delstar_piece(i,j)

			g_Lambda = grav_potential(x_coord(i),z_coord(j))

			g_mtheta=mach_theta(psi(i,j))

			g_indi=i
			g_indj=j

			! Calculate the derivatives of big_Psi (derivative of psi is not needed)

			call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

			g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

			! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

			if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
			endif

			! Calculate a minimum density
			call get_n_bound_min(n_den_min)
			if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
			! "Only a sub-Alfvenic Root!"
			else
				n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
			end if
			if(n_den_min>n_den_max) then
				n_den_min = 1.001d0*g_phi**2*mass*mu_mag
			endif

			if(n_den_max<n_den_min) then
				continue
			endif

			!Save n_den_min and n_den_max, since we'll use them again
			n_den_min_2D(i,j) = n_den_min
			n_den_max_2D(i,j) = n_den_max

			! this is smooth, so we can replace bisection with secant!
			n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
			Bernmax_2D(i,j) = n_den_ext
			! Next we evaluate the Bernoulli function
			tmp = bernoulli_TF(n_den_ext)
			fmax_2D(i,j) = tmp
			Bern_ave = Bern_ave + tmp

			if(tmp<min_Bern) then
				min_Bern = tmp
				min_ix = i
				min_iz = j
			endif

		enddo
		enddo

		min_step(min_ite) = min_Bern
		Bern_ave = Bern_ave/(nx*nz)
		min_err = abs(min_Bern/Bern_ave)

		if((min_err<min_tol).and.(min_Bern>=0.d0)) then
			!convergence
			min_search = .false.
		else
			!update mach_theta_max
			if(min_ite==1) then
			! distinguish the first update
				if(min_Bern>0.d0) then
					! we need to increase mach_theta_max
					mach_theta_max = mach_theta_max*(1.d0+min(1.d-2,min_err))
				else
					! we need to decrease mach_theta_max
					mach_theta_max = mach_theta_max*(1.d0-min(1.d-2,min_err))
				endif
			else
			! we can proceed with the secant part of the mach_thet_max search
				fprim = (min_step(min_ite)-min_step(min_ite-1))/  &
					(mtm_step(min_ite)-mtm_step(min_ite-1))
				if(fprim==0.d0) then
					print*, "Warning: fprim=0 in update_rho_TF_Gauss_2020_magnetic"
					if(min_Bern>0.d0) then
						! we need to increase mach_theta_max
						mach_theta_max = mach_theta_max*(1.d0+min(1.d-2,min_err))
					else
						! we need to decrease mach_theta_max
						mach_theta_max = mach_theta_max*(1.d0-min(1.d-2,min_err))
					endif
				else
					mach_theta_max = mach_theta_max - min_step(min_ite)/fprim
				endif
			endif

			! Last, make sure that mach_theta_max is not changing too much (it can happen in the first few iterations).
			if(max(mtm_step(min_ite)/mach_theta_max,mach_theta_max/mtm_step(min_ite))>1.1d0) then
				if(mach_theta_max>mtm_step(min_ite)) then
					mach_theta_max = mtm_step(min_ite) * (1.d0+1.d-2)
				else
					mach_theta_max = mtm_step(min_ite) * (1.d0-1.d-2)
				endif
			endif

		endif

	enddo

	!---------------------------------------------------------------------------------
	! Now we have mach_theta_max and can proceed to calculate the roots of the Bernoulli equation
	!---------------------------------------------------------------------------------

	write(id_Bern,1818)  calls, min_ix, min_iz, x_coord(min_ix), z_coord(min_iz), min_Bern

1818   format(3(I5,3x), 3(e14.8, 3x))

	do j=1, nz
	do i=1, nx
	! Only solve the inner region problem

		if(sort_grid(i,j)<=0) then
			cycle
		end if

		! Could get rid of this, but we can leave it for now to keep all bases covered.
		if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
			.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
			.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
		if(gravity_type==0) then
			n_den(i,j) = dofpsi(0.d0)
		else
			n_den(i,j) = ( (gamma-1.d0)/gamma *  &		! TO BE FIXED OR REMOVED
				(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
				/sofpsi(0.d0)  &
				)**(1.d0/(gamma-1.d0))
		endif
			cycle
		end if

		if(patch_value(i,j)) then
			! This can be removed
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		endif

		!  "Setup Global Parameters"
		! Set the globals for solving Bernouli equation.
		! It is not worth saving all these values from the previous part.
		psi_Bern = psi(i,j)
		big_Psi_Bern = big_Psi(i,j)
		g_r = x_coord(i)
		g_Phi = phi_TF_ofpsi(big_Psi(i,j))
		g_D = D_TF_ofpsi(big_Psi(i,j))
		! print *, "g_S"
		g_S_e = s_e_ofpsi(psi(i,j))
		g_S_i = s_i_ofpsi(big_Psi(i,j))
		g_S = g_S_e +  g_S_i
		! print *, "g_H"
		g_H_e = h_e_ofpsi_partial(psi(i,j))
		g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
		call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		g_psi_diff = psi_diff(i,j)
		g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
		g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
		g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

		call get_Delstar_piece(i,j)

		g_Lambda = grav_potential(x_coord(i),z_coord(j))

		g_mtheta=mach_theta(psi(i,j))

		g_indi=i
		g_indj=j

		! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

		! Now solve the modified Bernoulli equation using all the information from the previous part.

		Bernmax = Bernmax_2D(i,j)
		fBernmax = fmax_2D(i,j)
		n_den_min = n_den_min_2D(i,j)
		n_den_max = n_den_max_2D(i,j)

		if((Broot == 0).or.(Broot == 5)) then
			! Put a low density (rapidly moving) region
			! around the high density (slowly moving) core
			!             if(psi(i,j) >= psi_degen) then
			if(psi(i,j) >= psi_degen) then
			! Choose the last bracket (highest density)
				heavy_gauss = rtsafe(newt_rap_TF_gauss,Bernmax,n_den_max,x_TF_acc,10000)
				n_den(i,j) = heavy_gauss
			else
			! Choose the first bracket (lowest density)
				light_gauss = rtsafe(newt_rap_TF_gauss,n_den_min,Bernmax,x_TF_acc,10000)
				n_den(i,j) = light_gauss
			end if
		elseif((Broot == 1).or.(Broot == 4)) then
		! Choose the heavy root
			heavy_gauss = rtsafe(newt_rap_TF_gauss,Bernmax,n_den_max,x_TF_acc,10000)
			n_den(i,j) = heavy_gauss
		else if(Broot == 2) then
		! choose the light root
			light_gauss = rtsafe(newt_rap_TF_gauss,n_den_min,Bernmax,x_TF_acc,10000)
			n_den(i,j) = light_gauss
		else
			! The code will never get here, unless this routine is called for all grids and no single-fluid
			! Bernoulli solver is called.
			print *, "I don't understand Broot =",Broot
			pause
			stop
		endif

	enddo
	enddo
 
 	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
!!$	real(kind=dkind) :: n_bound_max_old
	real(kind=dkind) :: t_b, t_c, t_d
!!$	real(kind=dkind) :: t_b_old, t_c_old, t_d_old

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_c_old = mass * gamma/(gamma-1.d0) * g_S
!!$	t_d_old = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_i,gamma_e)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_max_old = ((t_d_old+t_b_old)/t_c_old)**(1.d0/(gamma-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_c (old, new) = ', t_c_old, t_c
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_max (old, new) = ', n_bound_max_old, n_bound_max
!!$	pause

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
!!$	real(kind=dkind) :: n_bound_min_old
	real(kind=dkind) :: t_a, t_b, t_d
!!$	real(kind=dkind) :: t_a_old, t_b_old, t_d_old

!!$	t_a_old = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_d_old = g_H

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_min_old = sqrt(t_a_old/(t_d_old + t_b_old))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_a (old, new) = ', t_a_old, t_a
!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_min (old, new) = ', n_bound_min_old, n_bound_min
!!$	pause


	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

end subroutine update_rho_TF_gauss_2020_magnetic

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density
! this routine does both relaxing and patching


	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex, ez, dx2, dz2
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.d-2
	endif

!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    do j=1, nz
       do i=1, nx
			n_old(i,j) = n_den(i,j)
		enddo
	enddo

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

		  if(patch_value(i,j)) then
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j))
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		  endif

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j))
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j))
          g_S_i = s_i_ofpsi(big_Psi(i,j))
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j))
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)
		  g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
		  g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
		  g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

			call get_Delstar_piece(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j))

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm

				   if (mach_theta_max<5.d-2*mach_theta_max_initial) then
						! just patch the point and hope for the best
						patch_value(i,j) = .true.
						n_den(i,j) = g_D
						goto 101
				   endif

				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif

                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j))
				   g_H_e = h_e_ofpsi_partial(psi(i,j))
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
						n_den_min = dmax1(n_den_min,1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm
				   endif
				   do while(mtm_soln-dmtm<=0.d0)
						dmtm = dmtm*(1.0-1.d-3)
					enddo
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  if(write_TF_roots) then
			  root_diff(i,j) = (heavy-light)/(heavy + light)
			  heavy_2D(i,j) = heavy
			  light_2D(i,j) = light
		  endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)
    big_Psi_degen = big_Psi(min_ix,min_iz)

    ! -----------------------------------------------------

	if(write_TF_roots) then
		call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
		call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
		call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)
	endif

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz))
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz))
       g_S = g_S_e +  g_S_i
	  g_S_i_prime = dS_i_dpsi(big_Psi(min_ix,min_iz))
	  g_Fstar = Fstarofpsi(psi(min_ix,min_iz),big_Psi(min_ix,min_iz))
	  g_psi_diff = psi_diff(min_ix,min_iz)

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_e,gamma_i)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!	endif

	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

end subroutine update_rho_TF_relax

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_no_patching(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex, ez, dx2, dz2
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*1.d-2
	endif

!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    do j=1, nz
       do i=1, nx
			n_old(i,j) = n_den(i,j)
		enddo
	enddo

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j))
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j))
          g_S_i = s_i_ofpsi(big_Psi(i,j))
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j))
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)
		  g_H_i_prime = dh_i_dpsi(big_Psi(i,j))
		  g_S_i_prime = dS_i_dpsi(big_Psi(i,j))
		  g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j))

			call get_Delstar_piece(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j))

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
 !            cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j))
				   g_H_e = h_e_ofpsi_partial(psi(i,j))
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
						n_den_min = dmax1(n_den_min,1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm
				   endif
				   do while(mtm_soln-dmtm<=0.d0)
						dmtm = dmtm*(1.0-1.d-3)
					enddo
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
!             if(psi(i,j) >= psi_degen) then
             if(big_Psi(i,j) >= big_Psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  if(write_TF_roots) then
			  root_diff(i,j) = (heavy-light)/(heavy + light)
			  heavy_2D(i,j) = heavy
			  light_2D(i,j) = light
		  endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)
    big_Psi_degen = big_Psi(min_ix,min_iz)

    ! -----------------------------------------------------

	if(write_TF_roots) then
		call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
		call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
		call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)
	endif

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz))
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz))
       g_S = g_S_e +  g_S_i
	  g_S_i_prime = dS_i_dpsi(big_Psi(min_ix,min_iz))
	  g_Fstar = Fstarofpsi(psi(min_ix,min_iz),big_Psi(min_ix,min_iz))
	  g_psi_diff = psi_diff(min_ix,min_iz)

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz))

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_e,gamma_i)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!	endif

	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max_true(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))

	continue

end subroutine get_n_bound_max_true

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min_true(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_min = sqrt(t_a/(t_d-t_b))

	continue

end subroutine get_n_bound_min_true

 end subroutine update_rho_TF_no_patching

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_version_1(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) term1
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j))
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j))
          g_S_i = s_i_ofpsi(big_Psi(i,j))
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j))
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j))

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j))
				   g_H_e = h_e_ofpsi_partial(psi(i,j))
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
						n_den_min = dmax1(n_den_min,1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  root_diff(i,j) = (heavy-light)/(heavy + light)
		  heavy_2D(i,j) = heavy
		  light_2D(i,j) = light

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

	call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
	call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
	call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz))
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz))
       g_S = g_S_e +  g_S_i

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
		n_bound_max = g_D * 2.d0
	else
		n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))
	endif

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_min'
		continue
		n_bound_min = g_D/1.d2
	else
		n_bound_min = sqrt(t_a/(t_d-t_b))
	endif

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max_true(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))

	continue

end subroutine get_n_bound_max_true

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min_true(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_min = sqrt(t_a/(t_d-t_b))

	continue

end subroutine get_n_bound_min_true

 end subroutine update_rho_TF_version_1


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_relax_version_1(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) term1
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14

!		if(Broot<4) then
!			mtm_limit = 1.d0
!		else
			mtm_limit = 1.d2
!		endif

!!$	do j=1,nz
!!$		do i=1,nx
!!$				if((bc_type/=7).and.(bc_type/=8)) then
!!$					psi(i,j) = dabs(psi_in(i,j))
!!$				else
!!$					psi(i,j) = psi_in(i,j)
!!$				endif
!!$		enddo
!!$	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
				if(gravity_type==0) then
					n_den(i,j) = dofpsi(0.d0)
				else
					n_den(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j))
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j))
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j))
          g_S_i = s_i_ofpsi(big_Psi(i,j))
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j))
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j))

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(n_den_min,n_den_max,1000)
!!!					pause
!!!					n_den(i,j) = dofpsi(psi(i,j))
!!!					cycle

!!3_11				if(Broot==2) then
!!3_11
!!3_11					! The light root
!!3_11					xb1(1) = n_den_min
!!3_11					xb2(1) = n_den_ext
!!3_11					! The heavy root
!!3_11					xb1(2) = n_den_ext
!!3_11					xb2(2) = n_den_max

!!3_11					if(nx>=65) then

!!3_11						print*, 'problem in Bernoulli solver: '
!!3_11						print*, 'no solution found with prescribed Mach theta!'
!!3_11						print*, i, j
!!3_11						print*, '--------------------------'

!!3_11					endif

!!3_11					cycle

!!3_11				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j))
				   g_H_e = h_e_ofpsi_partial(psi(i,j))
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j))
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
						n_den_min = dmax1(n_den_min,1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light * fix_orp_rho + n_den(i,j) * (1.d0-fix_orp_rho)
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  root_diff(i,j) = (heavy-light)/(heavy + light)
		  heavy_2D(i,j) = heavy
		  light_2D(i,j) = light

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

	call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
	call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
	call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz))
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz))
       g_S = g_S_e +  g_S_i

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz))
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz))
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz))
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
		n_bound_max = g_D * 2.d0
	else
		n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))
	endif

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_min'
		continue
		n_bound_min = g_D/1.d2
	else
		n_bound_min = sqrt(t_a/(t_d-t_b))
	endif

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max_true(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))

	continue

end subroutine get_n_bound_max_true

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min_true(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_min = sqrt(t_a/(t_d-t_b))

	continue

end subroutine get_n_bound_min_true

 end subroutine update_rho_TF_relax_version_1


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_TF_free_boundary(psi,big_Psi,psi_diff,n_den,nx,nz,seek_mtm,mtm_acc,min_dn_den,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The two-fluid routine updates the number density, not the mass density
! Jul 24 2013: add a patch for values that do not converge (should be necessary only at 
! the beginning of each grid)
! This version distinguishes between a "vacuum" and a plasma region

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, big_Psi,psi_diff
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: n_den
    real (kind=dkind), dimension(1:nx,1:nz) :: n_old ! previous iteration density
	logical :: patch_value(1:nx,1:nz)
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta n_den for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_n_den, tmp
    real (kind=dkind), intent(inout) :: min_dn_den
    integer, intent(inout) :: min_ix, min_iz
    ! We use n_den_ext to store n_den for which Bernoulli function is a maximum
    real (kind=dkind) :: n_den_ext ! n_den extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex, ez, dx2, dz2
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: n_den_max,n_den_min
	real (kind=dkind) :: drho_ds
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit
	real(kind=dkind) :: x_TF_acc
	real(kind=skind), dimension(1:nx,1:nz) :: root_diff, light_2D, heavy_2D
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	integer :: i_zone !0 for plasma, -1 for "vacuum"
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92

	integer iii,jjj,k

    ! -----------------------------------------------------

	x_TF_acc = d_TF_center * 1.d-14
	dx2 = dx*dx
	dz2 = dz*dz

	patch_value = .false.

			mtm_limit = 1.d2

	if((Broot/=0).and.(mach_theta_max/=mach_theta_max_initial)) then
		mach_theta_max = mach_theta_max + (mach_theta_max_initial-mach_theta_max)*2.5d-2
	endif


!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_dn_den = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    do j=1, nz
       do i=1, nx
			n_old(i,j) = n_den(i,j)
		enddo
	enddo

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j)<=0) then
				 cycle
          end if

		if((i==i_issue_check).and.(j==j_issue_check)) then
			continue
		endif

!		  n_den(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then

				i_zone = -1

				 !cycle !February 2021 MAJOR CHANGE: we also solve the Bernoulli equation in the "vacuum" region
          else

				i_zone = 0

		  end if

		  if(patch_value(i,j)) then
			n_den(i,j) = D_TF_ofpsi(big_Psi(i,j),i_zone)
			print*, 'patching density in i = ', i, ' j = ', j
			cycle
		  endif

		! CANCELLARE!!!!!!!!!!!
		!-----------------------------------------
			g_D = D_TF_ofpsi(big_Psi(i,j),i_zone)
		!-----------------------------------------

          !  "Setup Global Parameters"
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phi_TF_ofpsi(big_Psi(i,j),i_zone)
          ! print *, "g_S"
          g_S_e = s_e_ofpsi(psi(i,j),i_zone)
          g_S_i = s_i_ofpsi(big_Psi(i,j),i_zone)
          g_S = g_S_e +  g_S_i
          ! print *, "g_H"
          g_H_e = h_e_ofpsi_partial(psi(i,j),i_zone)
          g_H_i = h_i_ofpsi_partial(big_Psi(i,j),i_zone)
          call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
          g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_psi_diff = psi_diff(i,j)
		  g_H_i_prime = dh_i_dpsi(big_Psi(i,j),i_zone)
		  g_S_i_prime = dS_i_dpsi(big_Psi(i,j),i_zone)
		  g_Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j),i_zone)

			call get_Delstar_piece(i,j)

		  g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  g_mtheta=mach_theta(psi(i,j),i_zone)

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(i,j,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

		  g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

          ! Calculate bounds for n_den

		  call get_n_bound_max(n_den_max)

		  if (n_den_max<=0.d0) then
				write(*,*) 'error in n_den_max = ',n_den_max
				write(*,*) 'i,j = ', i,j
				write(*,*) 'zone = ', i_zone
				pause
				stop
		  endif

          if(g_Phi == 0.0d0) then
             ! There is only 1 root given by n_den_max
!			 n_den(i,j) = n_den_max
			 continue
!             cycle
		  end if

          ! Otherwise, increase n_den_max a little mach_theta_max == mtm_limit
          n_den_max = 1.01d0*n_den_max

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.1d0*g_Phi**2*mass*mu_mag) then
				! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				if(n_den_min>n_den_max) then
					n_den_min = 1.001d0*g_phi**2*mass*mu_mag
				endif
				n_den_min = dmax1(n_den_min/10.d0**6,1.0d-31)

!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(n_den_max<n_den_min) then
				continue
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli_TF,n_den_min,n_den_max,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(((nb == 1).and.(g_Phi**2>0.d0)) .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(n_den_min,n_den_max,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to n_den.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(n_den).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
             ! Next we evaluate the Bernoulli function
             tmp = bernoulli_TF(n_den_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                n_den = n_den_ext
                mean_n_den = n_den_ext
                min_dn_den = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = n_den_min
                xb2(1) = n_den_ext
                ! The heavy root
                xb1(2) = n_den_ext
                xb2(2) = n_den_max
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif


                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm

				   if (mach_theta_max<5.d-2*mach_theta_max_initial) then
						! just patch the point and hope for the best
						patch_value(i,j) = .true.
						n_den(i,j) = g_D
						goto 101
				   endif

				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(n_den_min,n_den_max,1000)
						pause
						stop
				   endif

                   ! Set the globals which depend on Mach Theta
                   g_Phi = phi_TF_ofpsi(big_Psi(i,j),i_zone)
				   g_H_e = h_e_ofpsi_partial(psi(i,j),i_zone)
				   g_H_i = h_i_ofpsi_partial(big_Psi(i,j),i_zone)
				   call get_H_diff(psi(i,j),big_Psi(i,j),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
				   g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				   ! Calculate bounds for n_den

					call get_n_bound_max(n_den_max)

                   ! Calculate a minimum density
						call get_n_bound_min(n_den_min)
						if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
						! "Only a sub-Alfvenic Root!"
						else
							n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
						end if
!						n_den_min = dmax1(n_den_min,1.0d-31)
						n_den_min = dmax1(n_den_min*10.d0**(-10),1.0d-31)

!!$					   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
                   ! Next we evaluate the Bernoulli function
                   tmp = bernoulli_TF(n_den_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = n_den_min
                      xb2(1) = n_den_ext
                      ! The heavy root
                      xb1(2) = n_den_ext
                      xb2(2) = n_den_max
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
!						dmtm = 10.0d0*dmtm
						dmtm = 2.0d0*dmtm
						! June 15 2019: changed this to 2
				   endif
				   do while(mtm_soln-dmtm<=0.d0)
						dmtm = dmtm*(1.0-1.d-3)
					enddo
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap_TF,xb1(1),xb2(1),x_TF_acc,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap_TF,xb1(2),xb2(2),x_TF_acc,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region
             ! around the high density (slowly moving) core
!             if(psi(i,j) >= psi_degen) then
             if(big_Psi(i,j) >= big_Psi_degen) then
                ! Choose the last bracket (highest density)
                n_den(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                n_den(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             n_den(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             n_den(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

		  if((heavy<1.d0).or.(heavy>1.d2*g_D).or.(light<1.d0).or.(light>1.d2*g_D)) then
			continue
		endif

		  if(write_TF_roots) then
			  root_diff(i,j) = (heavy-light)/(heavy + light)
			  heavy_2D(i,j) = heavy
			  light_2D(i,j) = light
		  endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_dn_den < 0.0d0) then
             ! initialization
             mean_n_den = 0.5d0*(heavy + light)
             min_dn_den = (heavy - light)/mean_n_den
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_dn_den) then
                mean_n_den = tmp
                min_dn_den = (heavy - light)/mean_n_den
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: n_den Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)
    big_Psi_degen = big_Psi(min_ix,min_iz)

    ! -----------------------------------------------------

	if(write_TF_roots) then
		call radial_plot(root_diff,psi,nx,nz,"root_diff",nx/2)
		call radial_plot(heavy_2D,psi,nx,nz,"heavy_2D",nx/2)
		call radial_plot(light_2D,psi,nx,nz,"light_2D",nx/2)
	endif

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min dn_den/n_den =",min_dn_den
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       g_S_e = s_e_ofpsi(psi(min_ix,min_iz),i_zone)
       g_S_i = s_i_ofpsi(big_Psi(min_ix,min_iz),i_zone)
       g_S = g_S_e +  g_S_i
	  g_S_i_prime = dS_i_dpsi(big_Psi(min_ix,min_iz),i_zone)
	  g_Fstar = Fstarofpsi(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),i_zone)
	  g_psi_diff = psi_diff(min_ix,min_iz)

          ! Calculate the derivatives of big_Psi (derivative of psi is not needed)

		  call psi_derivative(min_ix,min_iz,nx,nz,big_Psi,g_dpsidx,g_dpsidz)

       g_gPsi2 =  g_dpsidx**2+g_dpsidz**2

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
			g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz),i_zone)
			g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz),i_zone)
			g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz),i_zone)
			call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
			g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
			g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz),i_zone)

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

		  ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)

!!$			   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit (TF)"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm
		  endif
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
		  g_Phi = phi_TF_ofpsi(big_Psi(min_ix,min_iz),i_zone)
		  g_H_e = h_e_ofpsi_partial(psi(min_ix,min_iz),i_zone)
		  g_H_i = h_i_ofpsi_partial(big_Psi(min_ix,min_iz),i_zone)
		  call get_H_diff(psi(min_ix,min_iz),big_Psi(min_ix,min_iz),g_H_diff,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
		  g_H = g_H_e + g_H_i + eV*g_H_diff	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)
		  g_H_i_prime = dh_i_dpsi(big_Psi(min_ix,min_iz),i_zone)

		  call get_Delstar_piece(min_ix,min_iz)

		  ! Calculate bounds for n_den

			call get_n_bound_max(n_den_max)

          ! Calculate a minimum density
				call get_n_bound_min(n_den_min)
				if(n_den_min > 1.001d0*g_Phi**2*mass*mu_mag) then
					! "Only a sub-Alfvenic Root!"
				else
					n_den_min = 1.001d0*g_Phi**2*mass*mu_mag
				end if
				n_den_min = dmax1(n_den_min,1.0d-31)


!!$		   if(Broot==5) n_den_min = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          n_den_ext = rtbis(dbern_dn_den,n_den_min,n_den_max,1.0d-10)
          ! Evaluate the Bernoulli function
          tmp = bernoulli_TF(n_den_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

	contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Delstar_piece(ii,jj)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

	!-------------------------------------- get Del^star piece-------------------------------------- 
	big_Psil = 0.d0
	big_Psir = 0.d0
	big_Psiu = 0.d0
	big_Psid = 0.d0

	big_Psi0 = big_Psi(ii,jj)
	if(ii>1) big_Psil = big_Psi(ii-1,jj)
	if(ii<nx) big_Psir = big_Psi(ii+1,jj)
	if(jj>1) big_Psid = big_Psi(ii,jj-1)
	if(jj<nz) big_Psiu = big_Psi(ii,jj+1)

	drho_ds = van_Leer_slope_new(n_old(ii-1,jj),n_old(ii,jj),n_old(ii+1,jj),dx_a(ii-1),dx_a(ii))

	! Right x
	phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
	n_denrx = n_old(ii,jj) + 0.5d0*dx*drho_ds
	phinrx = phirx/n_denrx

	! Left x
	philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
	n_denlx = n_old(ii,jj) - 0.5d0*dx*drho_ds
	phinlx = philx/n_denlx

	! -----------------------------------------------------

	drho_ds = van_Leer_slope_new(n_old(ii,jj-1),n_old(ii,jj),n_old(ii,jj+1),dz_a(jj-1),dz_a(jj))

	! Right z
	phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
	n_denrz = n_old(ii,jj) + 0.5d0*dz*drho_ds
	phinrz = phirz/n_denrz

	! Left z
	philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
	n_denlz = n_old(ii,jj) - 0.5d0*dz*drho_ds
	phinlz = philz/n_denlz

	! use n_old in definition at least for now to save further modifications
	g_Dstar_term = - g_r*mass*phi_TF_ofpsi(big_Psi0)/eV/n_old(ii,jj) * (  &
				( phinrx*(big_Psir-big_Psi0)/(g_r + 0.5d0*dx) &
				+ phinlx*(big_Psil-big_Psi0)/(g_r - 0.5d0*dx) )  &
				/(g_r*dx2) &
				+ ( phinrz*(big_Psiu-big_Psi0) &
				+ phinlz*(big_Psid-big_Psi0) )/(g_r*g_r*dz2)  &
				)

	!-------------------------------------- got Del^star piece-------------------------------------- 

	continue

end subroutine get_Delstar_piece

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW, USING SOMETHING SIMILAR TO THE MHD LIMIT

	real(kind=dkind) :: n_bound_max
!!$	real(kind=dkind) :: n_bound_max_old
	real(kind=dkind) :: t_b, t_c, t_d
!!$	real(kind=dkind) :: t_b_old, t_c_old, t_d_old

	if(abs(g_phi)==0.d0) then
	! only one solution, easy to bracket
		n_bound_max = g_D * 5.d0
		return
	endif
	
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * (gamma_i/(gamma_i-1.d0) * g_S_i + gamma_e/(gamma_e-1.d0) * g_S_e)
	t_d = g_h

!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_c_old = mass * gamma/(gamma-1.d0) * g_S
!!$	t_d_old = g_h

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_max'
!		continue
!		n_bound_max = g_D * 2.d0
!	else
		n_bound_max = ((t_d+t_b)/t_c)**(1.d0/(min(gamma_i,gamma_e)-1.d0))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_max_old = ((t_d_old+t_b_old)/t_c_old)**(1.d0/(gamma-1.d0))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_c (old, new) = ', t_c_old, t_c
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_max (old, new) = ', n_bound_max_old, n_bound_max
!!$	pause

!!$	if (n_bound_max<=0.d0) then
!!$		write(*,*) 'error in n_den_max = ',n_bound_max
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_max = max(g_D*1.d2, n_bound_max)
		n_bound_max = n_bound_max*2.d0
	endif

	


	continue

end subroutine get_n_bound_max

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! GET A VERY ROUGH LIMIT FOR NOW
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
!!$	real(kind=dkind) :: n_bound_min_old
	real(kind=dkind) :: t_a, t_b, t_d
!!$	real(kind=dkind) :: t_a_old, t_b_old, t_d_old

!!$	t_a_old = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
!!$	t_b_old = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
!!$	t_d_old = g_H

	if(abs(g_phi)==0.d0) then
	! only one solution, easy to bracket
		n_bound_min = 0.d0
		return
	endif

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

!	if(t_d-t_b<=0.d0) then
!		print*, 'problem in n_den_min'
!		continue
!		n_bound_min = g_D/1.d2
!	else
		n_bound_min = sqrt(t_a/(t_d + t_b))	! NOTE THE (WRONG) "+" SIGN
!!$		n_bound_min_old = sqrt(t_a_old/(t_d_old + t_b_old))	! NOTE THE (WRONG) "+" SIGN
!	endif

!!$	print*, 't_a (old, new) = ', t_a_old, t_a
!!$	print*, 't_b (old, new) = ', t_b_old, t_b
!!$	print*, 't_d (old, new) = ', t_d_old, t_d
!!$	print*, 'n_bound_min (old, new) = ', n_bound_min_old, n_bound_min
!!$	pause


	n_bound_min = n_bound_min/10.d0

!!$	if (n_bound_min<=0.d0) then
!!$		write(*,*) 'error in n_den_min = ',n_bound_min
!!$		continue
!!$	endif

	if((psi(g_indi, g_indj)/psic)<0.9d0) then
!		n_bound_min  = min(1.d0, n_bound_min)
		n_bound_min = n_bound_min/2.d0
	endif

	continue

end subroutine get_n_bound_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_max_true(n_bound_max)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! gammas not fixed, since this routine is never called

	real(kind=dkind) :: n_bound_max
	real(kind=dkind) :: t_b, t_c, t_d

	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_c = mass * gamma/(gamma-1.d0) * g_S
	t_d = g_h

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_max = ((t_d-t_b)/t_c)**(1.d0/(gamma-1.d0))

	continue

end subroutine get_n_bound_max_true

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_n_bound_min_true(n_bound_min)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SOMETHING LOOKS WRONG HERE! CHECK PHI??
! Calculate a minimum density

	real(kind=dkind) :: n_bound_min
	real(kind=dkind) :: t_a, t_b, t_d

	t_a = 0.5d0*mass * g_gPsi2*g_phi**2/g_r**2
	t_b = 0.5d0 * ev**2/mass/g_r**2 * g_psi_diff**2
	t_d = g_H

	if(t_d-t_b<=0.d0) then
		print*, 'problem in n_den_max'
		continue
	endif

	n_bound_min = sqrt(t_a/(t_d-t_b))

	continue

end subroutine get_n_bound_min_true

 end subroutine update_rho_TF_free_boundary

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve_wrapper(psi,rho,residual,b_phi,big_Psi,psi_diff,n_den,nx,nz,  &
														min_it,max_it,eps,flag,orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! depending on input, calls a different routine to solve the GS-Bernoulli system

	use constant, only : dx, dz, dx_a, dz_a

    implicit none

	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho,residual,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: big_Psi,psi_diff,n_den
    real (kind=dkind) :: eps
    real (kind=dkind) :: orp
	character(len=8) :: cnumber
	character(len= 20) :: filename

!	if(((jump_option==1).or.(jump_option==2).or.(jump_option==3).or.(jump_option==4)).and.(Broot==0)) then
!	     call ngs_solve_jump(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
!	endif

	if(eq_type<10) then
		if(jump_option==0) then
			call ngs_solve(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
		elseif((jump_option==-7).or.(jump_option==-8)) then
			if((nx>=inter_switch).and.(Broot==0)) then
				call ngs_solve_all_gauss(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
			else
				call ngs_solve(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
			endif
		endif
	elseif(eq_type>=10) then
		if((jump_option==0).or.(abs(Broot)>0)) then
			call  ngs_solve_TF(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,orp)
		else
			if(nx>=inter_switch) then
				write (cnumber,Bern_format) nx
				filename='Bern_min'//trim(cnumber)//'.dat'
				open(unit=id_Bern,file=filename,action='write')
				if(jump_option==-7) then
					call ngs_solve_TF_gauss_flow(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,orp)
				elseif(jump_option==-8) then
					call ngs_solve_TF_gauss_magnetic(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,orp)
				endif
				close(id_Bern)
			else
				call  ngs_solve_TF(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,orp)
			endif
		endif
	endif

	continue

end subroutine ngs_solve_wrapper


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    implicit none
	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
    real (kind=dkind), intent(in) :: eps
    ! Input Over Relaxation Parameter
    ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
    ! else we set orp = in_orp
    real (kind=dkind), intent(in) :: in_orp
    integer ::ipass,i,isw,j,jsw,k,p,h,alloc_stat
    real (kind=dkind) :: res, den, x
    ! res -> residual
    ! den -> denominator of the Newton-Gauss-Seidel update
    ! x -> Radial position of grid point (i,j)
    real (kind=dkind) :: dx2,dz2,mx
    real (kind=dkind) :: anorm, eps_anorm
	real (kind=dkind) :: anorm2, eps_anorm2	! these 2 are fo b_phi
    real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
    real (kind=dkind) :: orp, std_orp ! over relaxation parameter
    ! Phi @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz
    real (kind=dkind) :: phirx,philx,phirz,philz
    ! Density @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    ! The square of the Alfvenic Mach Number
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
    real (kind=dkind) :: rhoc,phic,omegac,deltac,thetac,tparc,dc
    ! by is the phi component of the magnetic field
    ! b_dot_v is the inner product of B and V
    real (kind=dkind) :: by, b_dot_v
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: last_mtm ! last mach theta max
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: min_drho, tmp
    integer :: min_ix, min_iz
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: mtm_soln
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), parameter :: mtm_acc = 1.0d-12
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    real (kind=dkind) :: drho_ds ! Used with a slope limiter
	! The following is used to allow the possibility of anisotropic pressure
	real (kind=dkind) :: deltarx,deltalx,deltarz,deltalz
	real (kind=dkind) :: bphirx,bphilx,bphirz,bphilz
	real (kind=dkind), dimension (1:3,1:3) :: psi_delta,psi3x3
	real (kind=dkind) :: b_field_l,b_pol_l,dpsidx,dpsidz,psinew
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0
	real(kind=dkind) :: bpol_min_temp,bpol_max_temp, psi_pmax_temp, B_pmax_temp
	real (kind=dkind), dimension (-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: last_anorm
	real(kind=dkind) :: inorm = 0.d0

    real (kind=dkind), dimension(:,:,:), allocatable :: psi123
	integer :: h_steps=1
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: tic,toc
	integer :: pippa
	real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
	real(kind=dkind) :: orp_3step = 5.d-2

	integer iii,jjj

	psi_max = psic

	if((tri_type==13).and.(nx>=inter_switch).and.(nz>=inter_switch)) then
!	if((tri_type==13).and.(nx>inter_switch).and.(nz>inter_switch)) then

		h_steps = 3

	endif

	if(allocated(psi123)) deallocate(psi123)
	allocate(psi123(1:h_steps,1:nx,1:nz),stat = alloc_stat)
	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi123"
		 pause
		 stop
	endif
	! Uncomment below if you want to track psi123 in debugger
	!if(allocated(psi123)) print *, psi123

	psi123 = 0.d0

	do h=1,h_steps
	do j=1,nz
	do i=1,nx

		psi123(h,i,j) = psi_in(i,j)

	enddo
	enddo
	enddo

    if(psi_degen==0.d0) then 
		psi_degen = 0.5d0*psic ! Initialization
		big_Psi_degen = 0.5d0*psic ! Initialization
	endif
    last_mtm = mach_theta_max
    eps_anorm = 0.0d0
	eps_anorm2 = 0.0d0

	fix_orp0 = fix_orp

!!$	 	call step_output(nx,psi,rho,residual)

    ! The under/over relaxation parameter
    ! Standard Case:
    if(in_orp <= 0.0d0) then
       ! Note that the parameters below are tuned for starting with
       ! a 17 x 17 grid.
       orp = 1.0d0
       if(nx <= 5) then
          orp = 0.5d0
          rjac = 0.75d0
       else if(nx <= 9) then
          orp = 0.5d0
          rjac = 0.9d0
       else if(nx <= 17) then
          rjac = 0.98d0
       else if(nx <= 33) then
          rjac = 0.995d0
       else if(nx <= 65) then
          rjac = 0.9972d0
       else
!!$ It would appear that 129 is about where the solution begins to converge
          rjac = 0.0d0
       end if
       std_orp = orp
    else
       print *, "Input Over Relaxation Parameter =",in_orp
       orp = in_orp
       std_orp = in_orp
    endif

   if (accelerate) then
		continue
   else
	    orp = fix_orp
   endif

    dx2 = dx*dx
    dz2 = dz*dz

	! Added by Ian (NEW May 2025)
	if(bc_type==7) then
		call update_sort_grid(psi123(1,:,:),nx,nz,inorm)
	endif

! Update rho before "relaxing psi" but do not seek mach theta max
	call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
						min_ix,min_iz)

    ! Iterate until convergence
    do k=1, max_it
!!$     "Update Psi"

!		call cpu_time(tic)

       ! Reset the norm and max of Psi
       mx = 0.0d0
       anorm = 0.0d0
	   anorm2 = 0.0d0

	   do j=2,nz-1
	   do i=2,nx-1
			if ( ((tri_type==13).and.(sort_grid(i,j)==2)).or.  &
				((tri_type/=13).and.(bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j)==1)) ) then
					psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			! Removed elseif to conform with FLOW code (NEW May 2025)
			!elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
			!	psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			endif
	   enddo
	   enddo

	   bpol_max_temp = 0.d0
	   bpol_min_temp = 1.d22
	   psi_pmax_temp = 0.d0

!	   open(77,file='check.dat')
!	   open(88,file='checkw.dat')

		do h=1,h_steps
		! vertical stabilization cycle

       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
             do i=isw+1, nx-1, 2
                ! "Update Psi(",i,",",j,")"
                ! Only solve the inner region problem

				if(tri_type==10) then

!!$					if((ex*ex + ez*ez) >= radius_ext**2) then
!!$						cycle
!!$					endif

				else

					! Added bc_type exclusion (NEW May 2025)
					if(sort_grid(i,j)<=0.and.bc_type/=7) then
					   cycle
					end if

				endif

				! set up local psi values
				! this should save considerable resources in looking up values

				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi123(h,i,j)
				if(i>1) psil = psi123(h,i-1,j)
				if(i<nx) psir = psi123(h,i+1,j)
				if(j>1) psid = psi123(h,i,j-1)
				if(j<nz) psiu = psi123(h,i,j+1)

				psi_around(-1,0) = psil
				psi_around(1,0) = psir
				psi_around(0,0) = psi0
				psi_around(0,-1) = psid
				psi_around(0,1) = psiu

				! set up functions of psi

				psi_flag = 1.d9
				psi_flag_dep = 1.d9
				psi_flag_dep2 = 1.d9
				psi_flag_ham = 1.d9

				d_loc = dofpsi(psi0)
				dp_loc = dddpsi(psi0)

				p_loc = pofpsi(psi0)
				pp_loc = dpdpsi(psi0)

				psi_flag = psi0

				b0_loc = bzero(psi0)
				b0p_loc = dbzerodpsi(psi0)

				psi_flag_dep = psi0

				mth_loc = mach_theta(psi0)
				mthp_loc = dmach_thetadpsi(psi0)

				mph_loc = mach_phi(psi0)
				mphp_loc = dmach_phidpsi(psi0)

				psi_flag_dep2 = psi0

				s_loc = sofpsi(psi0)
				sp_loc = dsdpsi(psi0)

				phi_loc = phiofpsi(psi0)
				phip_loc = dphidpsi(psi0)

				omega_loc = omegaofpsi(psi0)
				omegap_loc = domegadpsi(psi0)

				i_loc = iofpsi(psi0)
				ip_loc = didpsi(psi0)

				h_loc = hofpsi(psi0)
				hp_loc = dhdpsi(psi0)

				psi_flag_ham = psi0

				! end of functions set up

                x = x_coord(i)

                ! Calculate rho, phi and omega at the current location
                rhoc = rho(i,j)
                phic = phiofpsi(psi0)
                omegac = omegaofpsi(psi0)

				dc = dofpsi(psi0)
                ! Calculate B_phi = by
				by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
					 (1.0d0 - phic*phic/rhoc)

                ! -----------------------------------------------------

                ! Calculate Phi & Density, rho @ +/- 0.5*dx, +/- 0.5*dz
                ! and then the square of the Alfvenic Mach Number.
!!$  NOTE: The use of a slope limiter to interpolate rho is essential
!!$  for eliminating an oscillation in Bpol where the solution switches
!!$  from the sub-slow to super-slow root. -- 3/14/2002 -- T. Gardiner
                ! -----------------------------------------------------

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimx(0) = ( dx_a(i-1)**2*psir +  &
									(dx_a(i)**2-dx_a(i-1)**2)*psi0 -  &
									dx_a(i)**2*psil ) /  &
									( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )

					ma2c = phic*phic/rhoc

				endif

!                drho_ds = van_Leer_slope(rho(i-1,j),rho(i,j),rho(i+1,j),dx)
                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                ! Right x
                phirx = phiofpsi(0.5d0*(psir + psi0))
                rhorx = rho(i,j) + 0.5d0*dx_a(i)*drho_ds
                ma2rx = phirx*phirx/rhorx
				psiprimx(1) = (psir - psi0)/dx_a(i)

                ! Left x
                philx = phiofpsi(0.5d0*(psil + psi0))
                rholx = rho(i,j) - 0.5d0*dx_a(i-1)*drho_ds
                ma2lx = philx*philx/rholx
				psiprimx(-1) = (psi0 - psil)/dx_a(i-1)

                ! -----------------------------------------------------

!                drho_ds = van_Leer_slope(rho(i,j-1),rho(i,j),rho(i,j+1),dz)
                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimz(0) = ( dz_a(j-1)**2*psiu +  &
									(dz_a(j)**2-dz_a(j-1)**2)*psi0 -  &
									dz_a(j)**2*psid ) /  &
									( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

				endif

                ! Right z
                phirz = phiofpsi(0.5d0*(psi0 + psiu))
                rhorz = rho(i,j) + 0.5d0*dz_a(j)*drho_ds
                ma2rz = phirz*phirz/rhorz
				psiprimz(1) = (psiu - psi0)/dz_a(j)

                ! Left z
                philz = phiofpsi(0.5d0*(psi0 + psid))
                rholz = rho(i,j) - 0.5d0*dz_a(j)*drho_ds
                ma2lz = philz*philz/rholz
				psiprimz(-1) = (psi0 - psid)/dz_a(j-1)

				! -----------------------------------------------------

				! calculate the magnetic field

					call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
					b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)
					b_field_l = dsqrt(by**2+b_pol_l**2)

                ! -----------------------------------------------------

                ! Calculate B dot v
                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)

                ! -----------------------------------------------------

				! March 12 2021: we don't need to distinguish the two zones this way for bc_type==7 anymore.
				! We are using the zone label "i_zone" instead.
				! January 21 2022: Adapted to the same option for tri_type==13
				! Added bc_type exclusion and modified if statment (NEW May 2025)
				if( ((tri_type/=13).and.(bc_type/=7).AND.((psi0/psic)<fraction))  &
								.OR.  &
								! ((bc_type==7).AND.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))  &
								! 				.OR.  &
								! ((tri_type==13).AND.(sort_grid(i,j)==1)) ) then
					((tri_type==13).AND.(sort_grid(i,j)<1)) ) then
					! OUTER REGION

						if(grid_type==0) then

							res = (1.d0/mu_mag)*(( (1.0d0 )*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 )*(psil-psi0)/(x - 0.5d0*dx) )/(x*dx2) &
								+ ( (1.0d0 )*(psiu-psi0) &
								  + (1.0d0 )*(psid-psi0) )/(x*x*dz2) )

						else

							res = (2.d0/mu_mag)*( (dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
										- dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										(dz_a(j-1)**2*psiprimz(1)   &
										- dz_a(j)**2*psiprimz(-1) +  &
										(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

						endif

				else	!inner region


						if(grid_type==0) then

							term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
							  + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
							  /(x*dx2) &
							+ ( (1.0d0 - ma2rz)*(psiu-psi0) &
							  + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )

						else

							term0 = (2.d0/mu_mag)*( ((1.0d0 - ma2rx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
										-(1.0d0 - ma2lx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1))   &
										+(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x )   &
										/(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1))   &

										+((1.0d0 - ma2rz)*dz_a(j-1)**2*psiprimz(1)   &
										-(1.0d0 - ma2lz)*dz_a(j)**2*psiprimz(-1)   &
										+(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) )   &
										/(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

						endif

						term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag)
						term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0)
						term3= by*didpsi(psi0)/x/dsqrt(mu_mag)
						term4= rhoc*dhdpsi(psi0)
						term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

						res = term0+term1+term2+term3+term4-term5

					endif

				continue

				res = res*mu_mag

                ! Store the residual
                if( (flag == 1).or.(k == max_it).or. (modulo(k,25) == 0).or.(k==1) ) then
                   residual(i,j) = res
                end if

				if(res>1.d0) then
					continue
				endif

                ! -----------------------------------------------------

					den = -2.d0*( ((1.0d0 - ma2rx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
									(1.0d0 - ma2lx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
									(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

									((1.0d0 - ma2rz)*dz_a(j-1)**2/dz_a(j) +  &
									(1.0d0 - ma2lz)*dz_a(j)**2/dz_a(j-1) +  &
									(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
									/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
									(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

					continue

				if(h_steps==1) then
					psi123(1,i,j) = psi0 - orp*res/den
				elseif((h_steps==3).and.(h<3)) then

					psi123(h+1,i,j) = psi0 - res/den

				elseif((h_steps==3).and.(h==3)) then

						psi123(1,i,j) = (1.d0-orp_3step)*psi123(1,i,j) + 2.d0*orp_3step*psi123(2,i,j)  &
										- orp_3step*psi123(3,i,j)

				endif

				if(h==h_steps) then

					! Find the max absolute value of psi in the plasma
					if((tri_type==13).and.(sort_grid(i,j)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((bc_type==7).and.(sort_grid(i,j)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((tri_type/=13).and.(bc_type/=7)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					endif

				endif

                ! Calculate the norm of the residual error
!                if ( (tri_type==10).and.((ex*ex + ez*ez) > rminor**2)  ) then
!					continue
!				else
					anorm = anorm + dabs(res)
!				endif

             end do
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do

		if((h_steps==3).and.(h<3)) then

		   call bc_psi_rho0(psi123(h+1,:,:),rho,nx,nz)

		endif

	   enddo
	   ! end of the vertical stability cycle

!	   close(77)
!	   close(88)

       ! Set the new maximum value of psic */
       psic = mx
	   if(tri_type==13)  then
		psic_13 = psic - psi_e
	  else
		psic_13 = psic
	  endif
       ! "Finished: Psi Update"
		psi_max = psic

       ! -----------------------------------------------------

       ! Move the density calculation to a separate function

	! Update rho and seek mach theta max
		call update_rho(psi123(1,:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
							min_ix,min_iz)

	   call update_b(psi123(1,:,:),rho,b_phi,nx,nz,orp,anorm2)

!   	   call bc_psi_rho0(psi123(1,:,:),rho,nx,nz)

	    ! Commented out to fix free-boundary bug (NEW May 2025)
		!do j=1,nz
		!do i=1,nx

		!	psi_in(i,j) = psi123(1,i,j)

		!enddo
		!enddo

   	   !call bc_psi_rho0(psi_in,rho,nx,nz)

		!do j=1,nz
		!do i=1,nx

		!	psi123(1,i,j) = psi_in(i,j)

		!enddo
		!enddo

       if(in_orp <= 0.0d0) then
          ! Chebyshev acceleration
          if(nx >= 5) then
             std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
          end if

       endif

	   if (accelerate) then
		    orp = std_orp
	        if (bc_type==2) orp = dmin1(orp,max_orp)
	   else
	        orp = fix_orp
	   endif

       ! -----------------------------------------------------

       if(k == 1) then
          eps_anorm = eps*anorm
		  eps_anorm2 = eps*anorm2
		  eps_anorm2 = dmax1(1.d-6,eps_anorm2)
       else
          if((in_orp <= 0.0d0).and.accelerate) then
             ! As we get closer to convergence, we want the solution
             ! to relax.  So as anorm approaches eps_anorm we want
             ! the over relaxation parameter to go to some const. ~ 1
             ! Use x and mtm_soln as a temporary variable
             x = anorm/eps_anorm
             mtm_soln = 1.0d0 ! The limiting value for x ~ 1
             tmp = x/(x + orp/mtm_soln - 1.0d0)
             tmp = dmin1(tmp,1.0d0)
             orp = orp*tmp
             orp = dmax1(orp,mtm_soln)
          endif
       end if

	if( (k<=50).and.(k>=25) ) then
!		fix_orp = (fix_orp1-fix_orp0)/25.d0*(k-25.d0) + fix_orp0
		continue
	endif

!!$ 			call step_output(nx,psi123(1,:,:),rho,residual)

       if(k == 1 .or. modulo(k,25) == 0) then
          print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		  print *, k,": anorm2 = ",real(anorm2,skind)," eps_anorm2 = ",real(eps_anorm2,skind)
          print *, "The Over-Relaxation Parameter =",orp,std_orp

 			call step_output(nx,psi123(1,:,:),rho,residual)
			if ((k>25).and.(tri_type==13)) call update_interface(psi123(1,:,:),nx,inorm)
			! if statement below added to conform with FLOW code (NEW May 2025)
			if(bc_type==7) call update_sort_grid(psi123(1,:,:),nx,nz,inorm)
			continue
      end if

		write(111,*) k,anorm


!		call cpu_time(toc)
!
!		print*, 'iteration time = ', toc-tic
!		pause

!		call step_output(nx,psi,rho,residual)

       ! Check the convergence criteria
	   if(anorm < eps_anorm .and. k > min_it .and. anorm2 < eps_anorm2 .and. inorm<1.d-5) exit
    end do

	do j=1,nz
	do i=1,nx

		psi_in(i,j) = psi123(1,i,j)

	enddo
	enddo

	deallocate(psi123)

    print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
	print *, "anorm2 = ",anorm2," eps_anorm2 = ",eps_anorm2
    print *, "Average residual =",anorm/(nx*nz)

    if ((broot==0).or.(Broot==5)) print *, "Mach Theta Max =",mach_theta_max
    print *, "Final Solution has Psi Center = ",psic
    print *


	return

  end subroutine ngs_solve


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve_no_limiter(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    implicit none
	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
    real (kind=dkind), intent(in) :: eps
    ! Input Over Relaxation Parameter
    ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
    ! else we set orp = in_orp
    real (kind=dkind), intent(in) :: in_orp
    integer ::ipass,i,isw,j,jsw,k,p,h,alloc_stat
    real (kind=dkind) :: res, den, x
    ! res -> residual
    ! den -> denominator of the Newton-Gauss-Seidel update
    ! x -> Radial position of grid point (i,j)
    real (kind=dkind) :: dx2,dz2,mx
    real (kind=dkind) :: anorm, eps_anorm
	real (kind=dkind) :: anorm2, eps_anorm2	! these 2 are fo b_phi
    real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
    real (kind=dkind) :: orp, std_orp ! over relaxation parameter
    ! Phi @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz
    real (kind=dkind) :: phirx,philx,phirz,philz
    ! Density @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    ! The square of the Alfvenic Mach Number
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
    real (kind=dkind) :: rhoc,phic,omegac,deltac,thetac,tparc,dc
    ! by is the phi component of the magnetic field
    ! b_dot_v is the inner product of B and V
    real (kind=dkind) :: by, b_dot_v
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: last_mtm ! last mach theta max
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: min_drho, tmp
    integer :: min_ix, min_iz
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: mtm_soln
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), parameter :: mtm_acc = 1.0d-12
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    real (kind=dkind) :: drho_ds ! Used with a slope limiter
	! The following is used to allow the possibility of anisotropic pressure
	real (kind=dkind) :: deltarx,deltalx,deltarz,deltalz
	real (kind=dkind) :: bphirx,bphilx,bphirz,bphilz
	real (kind=dkind), dimension (1:3,1:3) :: psi_delta,psi3x3
	real (kind=dkind) :: b_field_l,b_pol_l,dpsidx,dpsidz,psinew
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0
	real(kind=dkind) :: bpol_min_temp,bpol_max_temp, psi_pmax_temp, B_pmax_temp
	real (kind=dkind), dimension (-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: last_anorm
	real(kind=dkind) :: inorm = 0.d0

    real (kind=dkind), dimension(:,:,:), allocatable :: psi123
	integer :: h_steps=1
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: tic,toc
	integer :: pippa
	real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
	real(kind=dkind) :: orp_3step = 5.d-2

	integer iii,jjj

	psi_max = psic

	if((tri_type==13).and.(nx>=inter_switch).and.(nz>=inter_switch)) then
!	if((tri_type==13).and.(nx>inter_switch).and.(nz>inter_switch)) then

		h_steps = 3

	endif

	if(allocated(psi123)) deallocate(psi123)
	allocate(psi123(1:h_steps,1:nx,1:nz),stat = alloc_stat)
	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi123"
		 pause
		 stop
	endif

	psi123 = 0.d0

	do h=1,h_steps
	do j=1,nz
	do i=1,nx

		psi123(h,i,j) = psi_in(i,j)

	enddo
	enddo
	enddo

    if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
    last_mtm = mach_theta_max
    eps_anorm = 0.0d0
	eps_anorm2 = 0.0d0

	fix_orp0 = fix_orp

!!$	 	call step_output(nx,psi,rho,residual)

    ! The under/over relaxation parameter
    ! Standard Case:
    if(in_orp <= 0.0d0) then
       ! Note that the parameters below are tuned for starting with
       ! a 17 x 17 grid.
       orp = 1.0d0
       if(nx <= 5) then
          orp = 0.5d0
          rjac = 0.75d0
       else if(nx <= 9) then
          orp = 0.5d0
          rjac = 0.9d0
       else if(nx <= 17) then
          rjac = 0.98d0
       else if(nx <= 33) then
          rjac = 0.995d0
       else if(nx <= 65) then
          rjac = 0.9972d0
       else
!!$ It would appear that 129 is about where the solution begins to converge
          rjac = 0.0d0
       end if
       std_orp = orp
    else
       print *, "Input Over Relaxation Parameter =",in_orp
       orp = in_orp
       std_orp = in_orp
    endif

   if (accelerate) then
		continue
   else
	    orp = fix_orp
   endif

    dx2 = dx*dx
    dz2 = dz*dz

! Update rho before "relaxing psi" but do not seek mach theta max
	call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
						min_ix,min_iz)

    ! Iterate until convergence
    do k=1, max_it
!!$     "Update Psi"

!		call cpu_time(tic)

       ! Reset the norm and max of Psi
       mx = 0.0d0
       anorm = 0.0d0
	   anorm2 = 0.0d0

	   do j=2,nz-1
	   do i=2,nx-1
			if ( ((tri_type==13).and.(sort_grid(i,j)==2)).or.  &
				((tri_type/=13).and.(bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j)==1)) ) then
					psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
				psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			endif
	   enddo
	   enddo

	   bpol_max_temp = 0.d0
	   bpol_min_temp = 1.d22
	   psi_pmax_temp = 0.d0

!	   open(77,file='check.dat')
!	   open(88,file='checkw.dat')

		do h=1,h_steps
		! vertical stabilization cycle

       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
             do i=isw+1, nx-1, 2
                ! "Update Psi(",i,",",j,")"
                ! Only solve the inner region problem

				if(tri_type==10) then

!!$					if((ex*ex + ez*ez) >= radius_ext**2) then
!!$						cycle
!!$					endif

				else

					if(sort_grid(i,j)<=0) then
					   cycle
					end if

				endif

				! set up local psi values
				! this should save considerable resources in looking up values

				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi123(h,i,j)
				if(i>1) psil = psi123(h,i-1,j)
				if(i<nx) psir = psi123(h,i+1,j)
				if(j>1) psid = psi123(h,i,j-1)
				if(j<nz) psiu = psi123(h,i,j+1)

				psi_around(-1,0) = psil
				psi_around(1,0) = psir
				psi_around(0,0) = psi0
				psi_around(0,-1) = psid
				psi_around(0,1) = psiu

				! set up functions of psi
				! NOTE: DIFFERENTIATE FOR EQ_TYPE=3 LATER ON

				psi_flag = 1.d9
				psi_flag_dep = 1.d9
				psi_flag_dep2 = 1.d9
				psi_flag_ham = 1.d9

				d_loc = dofpsi(psi0)
				dp_loc = dddpsi(psi0)

				p_loc = pofpsi(psi0)
				pp_loc = dpdpsi(psi0)

				psi_flag = psi0

				b0_loc = bzero(psi0)
				b0p_loc = dbzerodpsi(psi0)

				psi_flag_dep = psi0

				mth_loc = mach_theta(psi0)
				mthp_loc = dmach_thetadpsi(psi0)

				mph_loc = mach_phi(psi0)
				mphp_loc = dmach_phidpsi(psi0)

				psi_flag_dep2 = psi0

				s_loc = sofpsi(psi0)
				sp_loc = dsdpsi(psi0)

				phi_loc = phiofpsi(psi0)
				phip_loc = dphidpsi(psi0)

				omega_loc = omegaofpsi(psi0)
				omegap_loc = domegadpsi(psi0)

				i_loc = iofpsi(psi0)
				ip_loc = didpsi(psi0)

				h_loc = hofpsi(psi0)
				hp_loc = dhdpsi(psi0)

				psi_flag_ham = psi0

				! end of functions set up

                x = x_coord(i)

                ! Calculate rho, phi and omega at the current location
                rhoc = rho(i,j)
                phic = phiofpsi(psi0)
                omegac = omegaofpsi(psi0)

				dc = dofpsi(psi0)
                ! Calculate B_phi = by
				by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
					 (1.0d0 - phic*phic/rhoc)

                ! -----------------------------------------------------

                ! Calculate Phi & Density, rho @ +/- 0.5*dx, +/- 0.5*dz
                ! and then the square of the Alfvenic Mach Number.
!!$  NOTE: The use of a slope limiter to interpolate rho is essential
!!$  for eliminating an oscillation in Bpol where the solution switches
!!$  from the sub-slow to super-slow root. -- 3/14/2002 -- T. Gardiner
                ! -----------------------------------------------------

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimx(0) = ( dx_a(i-1)**2*psir +  &
									(dx_a(i)**2-dx_a(i-1)**2)*psi0 -  &
									dx_a(i)**2*psil ) /  &
									( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )

					ma2c = phic*phic/rhoc

				endif

!                drho_ds = van_Leer_slope(rho(i-1,j),rho(i,j),rho(i+1,j),dx)
                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                ! Right x
                phirx = phiofpsi(0.5d0*(psir + psi0))
                rhorx = (rho(i,j) + rho(i+1,j))/2.d0
                ma2rx = phirx*phirx/rhorx
				psiprimx(1) = (psir - psi0)/dx_a(i)

                ! Left x
                philx = phiofpsi(0.5d0*(psil + psi0))
!                rholx = rho(i,j) - 0.5d0*dx*drho_ds
                rholx = (rho(i,j)+rho(i-1,j))/2.d0
                ma2lx = philx*philx/rholx
				psiprimx(-1) = (psi0 - psil)/dx_a(i-1)

                ! -----------------------------------------------------

!                drho_ds = van_Leer_slope(rho(i,j-1),rho(i,j),rho(i,j+1),dz)
                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimz(0) = ( dz_a(j-1)**2*psiu +  &
									(dz_a(j)**2-dz_a(j-1)**2)*psi0 -  &
									dz_a(j)**2*psid ) /  &
									( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

				endif

                ! Right z
                phirz = phiofpsi(0.5d0*(psi0 + psiu))
!                rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                rhorz = (rho(i,j) + rho(i,j+1))/2.d0
                ma2rz = phirz*phirz/rhorz
				psiprimz(1) = (psiu - psi0)/dz_a(j)

                ! Left z
                philz = phiofpsi(0.5d0*(psi0 + psid))
!                rholz = rho(i,j) - 0.5d0*dz*drho_ds
                rholz = (rho(i,j) +rho(i,j-1))/2.d0
                ma2lz = philz*philz/rholz
				psiprimz(-1) = (psi0 - psid)/dz_a(j-1)

				! -----------------------------------------------------

				! calculate the magnetic field

					call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
!					call psi_derivative(i,j,nx,nz,psi123(h,:,:),dpsidx,dpsidz)
					b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)
					b_field_l = dsqrt(by**2+b_pol_l**2)

                ! Calculate B dot v
                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)

                ! -----------------------------------------------------

				if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
								.OR.  &
					((bc_type==7).AND.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))  &
								.OR.  &
					((tri_type==13).AND.(sort_grid(i,j)==1)) ) then
					! OUTER REGION

!!$					if(h_steps==1) then

						if(grid_type==0) then

							res = (1.d0/mu_mag)*(( (1.0d0 )*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 )*(psil-psi0)/(x - 0.5d0*dx) )/(x*dx2) &
								+ ( (1.0d0 )*(psiu-psi0) &
								  + (1.0d0 )*(psid-psi0) )/(x*x*dz2) )

						else

							res = (2.d0/mu_mag)*( (dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
										- dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										(dz_a(j-1)**2*psiprimz(1)   &
										- dz_a(j)**2*psiprimz(-1) +  &
										(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

						endif

				else	!inner region

!!$						if(h_steps==1) then

							if(grid_type==0) then

								term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
								  /(x*dx2) &
								+ ( (1.0d0 - ma2rz)*(psiu-psi0) &
								  + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )

							else

								term0 = (2.d0/mu_mag)*( ((1.0d0 - ma2rx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
											-(1.0d0 - ma2lx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1))   &
											+(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x )   &
											/(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1))   &

											+((1.0d0 - ma2rz)*dz_a(j-1)**2*psiprimz(1)   &
											-(1.0d0 - ma2lz)*dz_a(j)**2*psiprimz(-1)   &
											+(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) )   &
											/(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

							endif

						term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag)
						term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0)
						term3= by*didpsi(psi0)/x/dsqrt(mu_mag)
						term4= rhoc*dhdpsi(psi0)
						term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

						res = term0+term1+term2+term3+term4-term5

				continue

				endif

				res = res*mu_mag

                ! Store the residual
                if( (flag == 1).or.(k == max_it).or. (modulo(k,25) == 0).or.(k==1) ) then
                   residual(i,j) = res
                end if

				if(res>1.d0) then
					continue
				endif

                ! -----------------------------------------------------

						if(grid_type==0) then

					den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
				             (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
					      -( (1.0d0 - ma2rz) + &
						     (1.0d0 - ma2lz) )/(x*x*dz2)

				else

					den = -2.d0*( ((1.0d0 - ma2rx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
									(1.0d0 - ma2lx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
									(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

									((1.0d0 - ma2rz)*dz_a(j-1)**2/dz_a(j) +  &
									(1.0d0 - ma2lz)*dz_a(j)**2/dz_a(j-1) +  &
									(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
									/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
									(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

				endif

					continue

				if(h_steps==1) then
					psi123(1,i,j) = psi0 - orp*res/den
				elseif((h_steps==3).and.(h<3)) then

!!$						!INNER REGION

						psi123(h+1,i,j) = psi0 - res/den

!!$					endif


				elseif((h_steps==3).and.(h==3)) then

!!$						!INNER REGION

						psi123(1,i,j) = (1.d0-orp_3step)*psi123(1,i,j) + 2.d0*orp_3step*psi123(2,i,j)  &
										- orp_3step*psi123(3,i,j)

				endif

				if(h==h_steps) then

					! Find the max absolute value of psi in the plasma
					if((tri_type==13).and.(sort_grid(i,j)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((tri_type/=13).and.(bc_type/=7)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					endif

				endif

                ! Calculate the norm of the residual error
!                if ( (tri_type==10).and.((ex*ex + ez*ez) > rminor**2)  ) then
!					continue
!				else
					anorm = anorm + dabs(res)
!				endif

             end do
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do

		if((h_steps==3).and.(h<3)) then

		   call bc_psi_rho0(psi123(h+1,:,:),rho,nx,nz)

		endif

	   enddo
	   ! end of the vertical stability cycle

!	   close(77)
!	   close(88)

       ! Set the new maximum value of psic */
       psic = mx
	   if(tri_type==13)  then
		psic_13 = psic - psi_e
	  else
		psic_13 = psic
	  endif
       ! "Finished: Psi Update"
		psi_max = psic

       ! -----------------------------------------------------

       ! Move the density calculation to a separate function

	! Update rho and seek mach theta max
		call update_rho(psi123(1,:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
							min_ix,min_iz)

	   call update_b(psi123(1,:,:),rho,b_phi,nx,nz,orp,anorm2)

!   	   call bc_psi_rho0(psi123(1,:,:),rho,nx,nz)

		do j=1,nz
		do i=1,nx

			psi_in(i,j) = psi123(1,i,j)

		enddo
		enddo

   	   call bc_psi_rho0(psi_in,rho,nx,nz)

		do j=1,nz
		do i=1,nx

			psi123(1,i,j) = psi_in(i,j)

		enddo
		enddo

       if(in_orp <= 0.0d0) then
          ! Chebyshev acceleration
          if(nx >= 5) then
             std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
          end if

       endif

	   if (accelerate) then
		    orp = std_orp
	        if (bc_type==2) orp = dmin1(orp,max_orp)
	   else
	        orp = fix_orp
	   endif

       ! -----------------------------------------------------

       if(k == 1) then
          eps_anorm = eps*anorm
		  eps_anorm2 = eps*anorm2
		  eps_anorm2 = dmax1(1.d-6,eps_anorm2)
       else
          if((in_orp <= 0.0d0).and.accelerate) then
             ! As we get closer to convergence, we want the solution
             ! to relax.  So as anorm approaches eps_anorm we want
             ! the over relaxation parameter to go to some const. ~ 1
             ! Use x and mtm_soln as a temporary variable
             x = anorm/eps_anorm
             mtm_soln = 1.0d0 ! The limiting value for x ~ 1
             tmp = x/(x + orp/mtm_soln - 1.0d0)
             tmp = dmin1(tmp,1.0d0)
             orp = orp*tmp
             orp = dmax1(orp,mtm_soln)
          endif
       end if

	if( (k<=50).and.(k>=25) ) then
!		fix_orp = (fix_orp1-fix_orp0)/25.d0*(k-25.d0) + fix_orp0
		continue
	endif

!!$ 			call step_output(nx,psi123(1,:,:),rho,residual)

       if(k == 1 .or. modulo(k,25) == 0) then
          print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		  print *, k,": anorm2 = ",real(anorm2,skind)," eps_anorm2 = ",real(eps_anorm2,skind)
          print *, "The Over-Relaxation Parameter =",orp,std_orp

 			call step_output(nx,psi123(1,:,:),rho,residual)
			if ((k>25).and.(tri_type==13)) call update_interface(psi123(1,:,:),nx,inorm)
			continue
      end if

		write(111,*) k,anorm

		last_anorm = anorm
		reduce = .false.

!		call cpu_time(toc)
!
!		print*, 'iteration time = ', toc-tic
!		pause

!		call step_output(nx,psi,rho,residual)

       ! Check the convergence criteria
	   if(anorm < eps_anorm .and. k > min_it .and. anorm2 < eps_anorm2 .and. inorm<1.d-5) exit
    end do

	do j=1,nz
	do i=1,nx

		psi_in(i,j) = psi123(1,i,j)

	enddo
	enddo

	deallocate(psi123)

    print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
	print *, "anorm2 = ",anorm2," eps_anorm2 = ",eps_anorm2
    print *, "Average residual =",anorm/(nx*nz)

    if ((broot==0).or.(Broot==5)) print *, "Mach Theta Max =",mach_theta_max
    print *, "Final Solution has Psi Center = ",psic
    print *

	return

  end subroutine ngs_solve_no_limiter



! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! the Bernoulli equation is modified to have a continuous solution
! Hameiri's GS equation (modified with the new Bernoulli) is used everywhere
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_all_gauss(psi_in,rho,residual,  &
						b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz, dx_a, dz_a

  implicit none
  integer, intent(in) :: nx,nz,min_it,max_it,flag
  real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
  real (kind=dkind), intent(in) :: eps
  ! Input Over Relaxation Parameter
  ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
  ! else we set orp = in_orp
  real (kind=dkind), intent(in) :: in_orp
  integer ::ipass,i,isw,j,jsw,k,alloc_stat
  real (kind=dkind) :: res, den, x !, xl, xr
  real (kind=dkind) :: res_1, den_1
  real (kind=dkind) :: res_2, den_2
  ! res -> residual
  ! den -> denominator of the Newton-Gauss-Seidel update
  ! x -> Radial position of grid point (i,j)
  real (kind=dkind) :: dx2,dz2,mx
  real (kind=dkind) :: anorm, eps_anorm
  real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
  real (kind=dkind) :: orp, std_orp ! over relaxation parameter
  real (kind=dkind) :: rhoc,phic,omegac
  ! by is the phi component of the magnetic field
  real (kind=dkind) :: by
  real (kind=dkind) :: last_mtm ! last mach theta max
  ! minimum delta rho for bernoulli roots, the mean val. and location
  real (kind=dkind) :: min_drho, tmp
  integer :: min_ix, min_iz
  ! The next set of variables are used for the bisection search for
  ! locating the degenerate root for Mach Theta Max
  real (kind=dkind) :: mtm_soln
  ! mtm_acc controls the bisection loop seeking mach theta max
  real (kind=dkind), parameter :: mtm_acc = 1.0d-12
  ! VERY IMPORTANT
  ! psi_degen records the value of psi where the solution to the
  ! Bernoulli equation is approximately degenerate.  It is 
  ! approximately equal to 0.5d0*psic but generally slightly less.
  real (kind=dkind) :: gpsirx, gpsilx, gpsirz, gpsilz
  real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term6
  real(kind=dkind) :: last_anorm
  real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS
  integer :: dir_filter(1:nx,1:nz)
  real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
  real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
  real (kind=dkind) :: rhorx,rholx,rhorz,rholz
  real (kind=dkind) :: phirx,philx,phirz,philz
  real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
  real(kind=dkind) :: dpsidx, dpsidz, b_pol_l, b_dot_v, drho_ds
  real(kind=dkind) :: psi_distance ! how far to go from psi_degen in using the Ptot formulation
  real(kind=dkind) :: smooth_fact
  integer :: k_save = 100


  if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
  last_mtm = mach_theta_max
  eps_anorm = 0.0d0

  fix_orp0 = fix_orp

  ! The under/over relaxation parameter
  ! Standard Case:
  if(in_orp <= 0.0d0) then
     ! Note that the parameters below are tuned for starting with 
     ! a 17 x 17 grid.
     orp = 1.0d0
     if(nx <= 5) then
        orp = 0.5d0
        rjac = 0.75d0
     else if(nx <= 9) then
        orp = 0.5d0
        rjac = 0.9d0
     else if(nx <= 17) then
        rjac = 0.98d0
     else if(nx <= 33) then
        rjac = 0.995d0
     else if(nx <= 65) then
        rjac = 0.9972d0
     else
        ! It would appear that 129 is about where the solution begins to converge
        rjac = 0.0d0
     end if
     std_orp = orp
  else
     print *, "Input Over Relaxation Parameter =",in_orp
     orp = in_orp
     std_orp = in_orp
  endif

  if (accelerate) then
     continue
  else
     orp = fix_orp
  endif


  dx2 = dx*dx
  dz2 = dz*dz
  
  if(allocated(fmax_2D)) then
!  	if(size(fmax_2D,1)==nx) then
  		continue
 ! 	else
	  	deallocate(fmax_2D)
  !	endif
  endif

!  if(allocated(fmax_2D)) then
!	continue
 ! else
  	allocate(fmax_2D(1:nx,1:nz))
  	fmax_2D = 0.d0
  !endif


  ! Update rho before "relaxing psi" but do not seek mach theta max
  call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
       min_ix,min_iz)

        call step_output(nx,psi_in(:,:),rho,residual)


  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	  call GS_RHS_all_Gauss



  ! Iterate until convergence
  do k=1, max_it
     !     "Update Psi"

     !		call cpu_time(tic)

     ! Reset the norm and max of Psi
     mx = 0.0d0
     anorm = 0.0d0

     do j=2,nz-1
        do i=2,nx-1
           if(sort_grid(i,j)>0) psi_in(i,j) = dabs(psi_in(i,j))	!	segno
        enddo
     enddo


 !    if(nx>=inter_switch) call step_output(nx,psi_in(:,:),rho,residual)

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)



     jsw = 1
     do ipass = 1, 2
        isw=jsw;
        do j=2, nz-1
           do i=isw+1, nx-1, 2
              ! "Update Psi(",i,",",j,")"
              ! Only solve the inner region problem

              if(sort_grid(i,j)<=0) cycle

              ! set up local psi values

              psil = 0.d0
              psir = 0.d0
              psiu = 0.d0
              psid = 0.d0

              psi0 = psi_in(i,j)
              if(i>1) psil = psi_in(i-1,j)
              if(i<nx) psir = psi_in(i+1,j)
              if(j>1) psid = psi_in(i,j-1)
              if(j<nz) psiu = psi_in(i,j+1)


              x = x_coord(i)
!              xl = (x_coord(i)+x_coord(i-1))/2.d0
 !             xr = (x_coord(i+1)+x_coord(i))/2.d0

              ! Calculate rho, phi and omega at the current location
              rhoc = rho(i,j)
              phic = phiofpsi(psi0)
              omegac = omegaofpsi(psi0)

			  ! old formulation with new Bernoulli

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)



                 term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4= rhoc*dhdpsi(psi0) 
                 term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 term6 = -(gradpsi_RHS(1,i,j)*(fmax_2D(i+1,j)-fmax_2D(i-1,j))/dx/2.d0 +  &
                 		gradpsi_RHS(2,i,j)*(fmax_2D(i,j+1)-fmax_2D(i,j-1))/dz/2.d0)/gradpsi_RHS(0,i,j)**2 +  &
                 		2.d0*delta_Bern*(psi_in(i,j)-psi_degen)/psic**2*fmax_2D(i,j)

                 term6 = rhoc*term6 * exp(-delta_Bern*((psi_in(i,j)-psi_degen)/psic)**2)

!				if((nx>=inter_switch).and.(term6/=0.d0)) then
!					print*, i,j, term0, term1, term2, term3, term4, term5, term6
!					print*, '   '
!					pause
!				endif

                 res = term0+term1+term2+term3+term4-term5+term6


                 den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

		              continue



              !res = res*mu_mag


              ! Store the residual
              if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
!              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then
!					residual(i,j) = res/gradpsi_RHS(0,i,j)**2
!				else
					if(nx>=inter_switch) then
						residual(i,j) = res/den
					else
						residual(i,j) = res
					endif
!				endif
              end if

              if(res>1.d0) then
                 continue
              endif

              ! -----------------------------------------------------

              continue

              psi_in(i,j) = psi0 - orp*res/den

              ! Find the max absolute value of psi in the plasma
              mx = dmax1(dabs(psi_in(i,j)),mx)

              ! Calculate the norm of the residual error
				if(nx>=inter_switch) then
	              anorm = anorm + dabs(res/den)
				else
	              anorm = anorm + dabs(res)
				endif

           end do
           isw = 3 - isw
        end do
        jsw = 3 - jsw
     end do

     ! update RHS
     	if(nx>=inter_switch) then
			call GS_RHS_all_Gauss
		endif

     ! Set the new maximum value of psic */
     psic = mx
     psic_13 = psic
     ! "Finished: Psi Update"

     ! -----------------------------------------------------

     ! Move the density calculation to a separate function
     ! Update rho and seek mach theta max
	if((nx>=inter_switch).and.(Broot==0)) then
     call update_rho_gauss(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
    else
     call update_rho(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
	endif

     call bc_psi_rho0(psi_in,rho,nx,nz)

     if(in_orp <= 0.0d0) then
        ! Chebyshev acceleration 
        if(nx >= 5) then
           std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
        end if
     endif

     if ((accelerate).and.(nx<inter_switch)) then
        orp = std_orp
     else
		orp = fix_orp
     endif

        !!$ call step_output(nx,psi_in(:,:),rho,residual)

     ! -----------------------------------------------------

     if(k == 1) then
        eps_anorm = eps*anorm
     else
        if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
           ! As we get closer to convergence, we want the solution
           ! to relax.  So as anorm approaches eps_anorm we want 
           ! the over relaxation parameter to go to some const. ~ 1
           ! Use x and mtm_soln as a temporary variable
           x = anorm/eps_anorm
           mtm_soln = 1.0d0 ! The limiting value for x ~ 1
           tmp = x/(x + orp/mtm_soln - 1.0d0)
           tmp = dmin1(tmp,1.0d0)
           orp = orp*tmp
           orp = dmax1(orp,mtm_soln)
        endif
     end if

     if(k == 1 .or. modulo(k,k_save) == 0) then
        print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
        print *, "The Over-Relaxation Parameter =",orp,std_orp
        call step_output(nx,psi_in(:,:),rho,residual)
        continue
     end if

     write(111,*) k,anorm

     ! Check the convergence criteria
     if(anorm < eps_anorm .and. k > min_it) exit

  end do

	deallocate(fmax_2D)

  print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
  print *, "Average residual =",anorm/(nx*nz)

  print *, "Mach Theta Max =",mach_theta_max
  print *, "Final Solution has Psi Center = ",psic
  print *


  return


contains



  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS_all_Gauss
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind) :: tdx, tdz


    tdx = dx*2.d0
    tdz = dz*2.d0

    gradpsi_RHS = 1.d-10
    

    ! fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)

       enddo
    enddo


!	if((k==1).or.(modulo(k,k_save)==0)) then
!
!    ! write results for debugging
!
!		open(13, file='gradpsi_RHS.plt')
!
!		write(13,*)'TITLE="solution of GS equation with flow"'
!		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
!		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'
!
!		do j=1,nz
!
!		   z = z_coord(j)
!
!		   do i=1,nx
!
!			  x = x_coord(i)
!
!			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
!
!		   end do
!		end do
!
!		close(13)
!
!	endif

90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS_all_Gauss


end subroutine ngs_solve_all_gauss





  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz], big_Psi is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this is the two-fluid version
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_TF(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! cycle separately on psi and psi_diff
! IMPORTANT: solve directly for psi_diff, then calculate big_Psi as psi + psi_diff

	use constant, only : dx, dz, dx_a, dz_a

	implicit none

	integer, intent(in) :: nx,nz,min_it,max_it,flag
	real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,big_Psi,psi_diff,n_den,residual
	real (kind=dkind), intent(in) :: eps
	! Input Over Relaxation Parameter
	! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
	! else we set orp = in_orp
	real (kind=dkind), intent(in) :: in_orp
	integer ::ipass,i,isw,j,jsw,k,alloc_stat
	real (kind=dkind) :: res, den, x !, xl, xr
	! res -> residual
	! den -> denominator of the Newton-Gauss-Seidel update
	! x -> Radial position of grid point (i,j)
	real (kind=dkind) :: dx2,dz2,mx
	real (kind=dkind) :: anorm, eps_anorm, anorm_2, eps_anorm_2
	real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
	real (kind=dkind) :: orp, std_orp ! over relaxation parameter
	real (kind=dkind) :: n_denc, phic
	! by is the phi component of the magnetic field
	real (kind=dkind) :: Fstar
	real (kind=dkind) :: last_mtm ! last mach theta max
	! minimum delta rho for bernoulli roots, the mean val. and location
	real (kind=dkind) :: min_dn_den, tmp
	integer :: min_ix, min_iz
	! The next set of variables are used for the bisection search for
	! locating the degenerate root for Mach Theta Max
	real (kind=dkind) :: mtm_soln
	! mtm_acc controls the bisection loop seeking mach theta max
	real (kind=dkind), parameter :: mtm_acc = 1.0d-12
	! VERY IMPORTANT
	! psi_degen records the value of psi where the solution to the
	! Bernoulli equation is approximately degenerate.  It is 
	! approximately equal to 0.5d0*psic but generally slightly less.
	real (kind=dkind) :: gPsi2
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term0_5
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz, phinc
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	real(kind=dkind), dimension(-1:1,-1:1) :: big_Psi_around
	real(kind=dkind) :: psi_diff_old, psi_diff_new, psi_diff_limit
	real(kind=dkind) :: H_loc, H_diff_loc
	real(kind=dkind) :: dpsidx, dpsidz, drho_ds
!	integer :: psi_diff_option = 4
	real (kind=dkind), dimension(1:nx,1:nz) :: diff_error
	integer :: i_zone = 1 !This will ony be changed if bc_type==7
!	logical, save :: initialize_zones = .true.
	integer :: k_save = 25
	real(kind=dkind) :: inorm = 0.d0
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92

	psi_max = max(psic,big_Psic)

    if(psi_degen==0.d0) then 
		psi_degen = 0.5d0*psic ! Initialization
		big_Psi_degen = 0.5d0*psic ! Initialization
	endif
	last_mtm = mach_theta_max
	eps_anorm = 0.0d0
	eps_anorm_2 = 0.0d0

	fix_orp0 = fix_orp

	! The under/over relaxation parameter
	! Standard Case:
	if(in_orp <= 0.0d0) then
		! Note that the parameters below are tuned for starting with 
		! a 17 x 17 grid.
		orp = 1.0d0
		if(nx <= 5) then
			orp = 0.5d0
			rjac = 0.75d0
		else if(nx <= 9) then
			orp = 0.5d0
			rjac = 0.9d0
		else if(nx <= 17) then
			rjac = 0.98d0
		else if(nx <= 33) then
			rjac = 0.995d0
		else if(nx <= 65) then
			rjac = 0.9972d0
		else
		! It would appear that 129 is about where the solution begins to converge
			rjac = 0.0d0
		end if
		std_orp = orp
	else
		print *, "Input Over Relaxation Parameter =",in_orp
		orp = in_orp
		std_orp = in_orp
	endif

	if (accelerate) then
		continue
	else
		orp = fix_orp
	endif

	dx2 = dx*dx
	dz2 = dz*dz

	ipsic = nx/2
	jpsic = nz/2

!	if((initialize_zones).and.(bc_type==7)) then
	if(bc_type==7) then
		! Added inorm to update_sort_grid calls (NEW May 2025)
		if(LCFS==-1) then
			call update_sort_grid(psi(:,:),nx,nz,inorm)
		elseif(LCFS==1) then
			call update_sort_grid(big_Psi(:,:),nx,nz,inorm)
		endif
!		initialize_zones = .false.
	endif

	diff_error = 0.d0

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,0)

	! Update rho before "relaxing psi" but do not seek mach theta max
!!	if(Broot==0) then
!!		call update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	else
		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
!!	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)

	call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)

	! Iterate until convergence
	do k=1, max_it
	!     "Update psi"6

	!		call cpu_time(tic)

		if((k>ite_Broot).and.(Broot/=Broot_true)) then
			Broot = Broot_true
		endif

	! Reset the norm and max of psi and big_Psi
	mx = 0.0d0
	anorm = 0.0d0
	anorm_2 = 0.0d0
	res = 0.0d0 ! Added this line? (NEW May 2025)
	
	! Ian comment: not sure this is necessary here; double check? 
	do j=2,nz-1
	do i=2,nx-1
		if((bc_type/=7).and.(sort_grid(i,j)>0)) psi(i,j) = dabs(psi(i,j))	!	segno
	enddo
	enddo

	!----------------------psi cycle------------------------
	jsw = 1
	do ipass = 1, 2
		isw=jsw;
		do j=2, nz-1
		do i=isw+1, nx-1, 2
		! "Update psi(",i,",",j,")"
		! Only solve the inner region problem

			! Added bc_type exclusion (NEW May 2025)
			if(sort_grid(i,j)<=0.and.bc_type/=7) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			big_Psil = 0.d0
			big_Psir = 0.d0
			big_Psiu = 0.d0
			big_Psid = 0.d0

			big_Psi0 = big_Psi(i,j)
			if(i>1) big_Psil = big_Psi(i-1,j)
			if(i<nx) big_Psir = big_Psi(i+1,j)
			if(j>1) big_Psid = big_Psi(i,j-1)
			if(j<nz) big_Psiu = big_Psi(i,j+1)

			x = x_coord(i)

			! Set up zone. May want to do the same for tri_type==13 here. March 17 2021

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(i,j)-2)
				i_zone = max(-1,i_zone)
			endif

			! Calculate rho and phi at the current location
			n_denc = n_den(i,j)
			phic = phi_TF_ofpsi(big_Psi0,i_zone)

			! old formulation with new Bernoulli

			drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

			! March 17 2021
			! We are evaluating free functions in half grid points. Rather than testing for wich zone
			! each point is in, we will always assume that the are in the same zone as the central point. 

			! Right x
			phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0),i_zone)
			n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
			ma2rx = mu_mag*mass * phirx*phirx/n_denrx

			! Left x
			philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0),i_zone)
			n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
			ma2lx = mu_mag*mass * philx*philx/n_denlx

			! -----------------------------------------------------

			drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

			! Right z
			phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu),i_zone)
			n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
			ma2rz = mu_mag*mass * phirz*phirz/n_denrz

			! Left z
			philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid),i_zone)
			n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
			ma2lz = mu_mag*mass * philz*philz/n_denlz

			! -----------------------------------------------------

			big_Psi_around(-1,0) = big_Psil
			big_Psi_around(1,0) = big_Psir
			big_Psi_around(0,0) = big_Psi0
			big_Psi_around(0,-1) = big_Psid
			big_Psi_around(0,1) = big_Psiu

			Fstar = Fstarofpsi(psi0,big_Psi0,i_zone)

			call psi_derivative_new(i,j,nx,nz,big_Psi_around,dpsidx,dpsidz)
			gPsi2 =  dpsidx**2+dpsidz**2

			! we follow FLOW nomenclature for term names (everything should be easily recognizable)

			term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
							+ (psil-psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

			term0_5 = - (  &
							( ma2rx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ ma2lx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( ma2rz*(big_Psiu-big_Psi0) &
							+ ma2lz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

			term1 = eV*mu_mag * Fstar *(phic-phi_e_TF_ofpsi(psi0,i_zone))/x**2
			term3 = mu_mag*mass * phic * dphi_TF_dpsi(big_Psi0,i_zone) * gPsi2 / n_denc/x**2

			term2 = mu_mag * eV * n_denc * (omegaofpsi(big_Psi0,i_zone)-omegaofpsi(psi0,i_zone)) -  &
					mu_mag * n_denc * (dp_i_dpsi(big_Psi0,i_zone)/d_TF_ofpsi(big_Psi0,i_zone)-  &
					dp_i_dpsi(psi0,i_zone)/d_TF_ofpsi(psi0,i_zone))

			term4 = mu_mag * n_denc * (dh_e_dpsi_partial(psi0,i_zone) + dh_i_dpsi_partial(big_psi0,i_zone))

			term5 = mu_mag*mass * (n_denc**gamma_i/(gamma_i-1.0d0) * ds_i_dpsi(big_Psi0,i_zone) +  &
												n_denc**gamma_e/(gamma_e-1.0d0) * ds_e_dpsi(psi0,i_zone))
			! CHECK DIMENSIONS FOR ALL!!!

			!				if((nx>=inter_switch).and.(term6/=0.d0)) then
			!					print*, i,j, term0, term1, term2, term3, term4, term5, term6
			!					print*, '   '
			!					pause
			!				endif

			res = term0 + term0_5 + term1 + term2 + term3 + term4 - term5

			den = -( 1.0d0/(x + 0.5d0*dx) + &
						1.0d0/(x - 0.5d0*dx) )/(x*dx2) &
						-2.0d0/(x*x*dz2)

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

			continue

			!res = res*mu_mag

			! Store the residual
			if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
				if(nx>=inter_switch) then
					residual(i,j) = res/den
				else
					residual(i,j) = res
				endif
			end if

			if(res>1.d0) then
				continue
			endif

			! -----------------------------------------------------

			continue

			psi(i,j) = psi0 - orp*res/den

			! Find the max absolute value of psi in the plasma
			if((bc_type/=7).or.((bc_type==7).and.(sort_grid(i,j)==2))) then
				if(dabs(psi(i,j))>mx) then
					mx = dabs(psi(i,j))
					! Comment by Ian: this next if is just for tri_type==-4, which is deprecated for now
					if((tri_type==-4).and.(LCFS==-1)) then
						ipsic = i
						jpsic = j
					endif
				endif
			endif

			anorm = anorm + dabs(res)

		end do
		isw = 3 - isw
		end do
		jsw = 3 - jsw
		! Comment by Ian: this call doesn't happen in FLOW unless tri_type==13 and above inter_switch value
		call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,1)
		! only update psi on BC
	end do
	!----------------end of psi cycle--------------

	! Set the new maximum value of psic
	psic = mx
	psi_max = max(psic,big_Psic)
	psic_13 = psi_max
	! "Finished: psi Update"

	!----------------big_Psi cycle-----------------
	! we use psi_diff as dependent variable, and solve the algebraic equation for it
	! FIRST ATTEMPT: USE OLD ITERATION FOR DERIVATIVES

	! March 25 2021: added the option to use different equations in different regions of the plasma.
	! Also added red-black.
	
	if(((psi_diff_option==3).or.(psi_diff_option==4)).and.(bc_type==7)) then
	
       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
          do i=isw+1, nx-1, 2

				! Added bc_type exclusion (NEW May 2025)
				! TODO: This means that below if never triggers, so if this works this line should be removed entirely
				if(sort_grid(i,j)<=0.and.bc_type/=7) cycle

				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(i,j)-2)
				i_zone = max(-1,i_zone)

				! set up local psi values

				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi(i,j)
				if(i>1) psil = psi(i-1,j)
				if(i<nx) psir = psi(i+1,j)
				if(j>1) psid = psi(i,j-1)
				if(j<nz) psiu = psi(i,j+1)

				x = x_coord(i)

				psi_diff_old = psi_diff(i,j)

				n_denc = n_den(i,j)

				Fstar = Fstarofpsi(psi0,big_Psi(i,j),i_zone)

				if(i_zone==-1) then
					continue
				endif

				if(((i_zone==0).and.(psi_diff_option==3)).or.((i_zone==-1).and.(psi_diff_option==4))) then
				! use magnetic GS for psi_diff

					term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
									+ (psil-psi0)/(x - 0.5d0*dx) )  &
									/(x*dx2) &
									+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

					term2 = -eV*mu_mag * Fstar * phi_e_TF_ofpsi(psi0,i_zone) / x**2
					term3 = mu_mag * n_denc * dh_e_dpsi(psi0,i_zone)
					term4 = mu_mag * mass * n_denc**gamma_e * ds_e_dpsi(psi0,i_zone)/(gamma_e-1.0d0)

					psi_diff_new = (-term0 - term2 - term3 + term4) / (ev**2*mu_mag*n_denc/mass/x**2)

				elseif(((i_zone==0).and.(psi_diff_option==4)).or.((i_zone==-1).and.(psi_diff_option==3))) then
				! use fluid GS for psi_diff

					big_Psil = 0.d0
					big_Psir = 0.d0
					big_Psiu = 0.d0
					big_Psid = 0.d0

					big_Psi0 = big_Psi(i,j)
					if(i>1) big_Psil = big_Psi(i-1,j)
					if(i<nx) big_Psir = big_Psi(i+1,j)
					if(j>1) big_Psid = big_Psi(i,j-1)
					if(j<nz) big_Psiu = big_Psi(i,j+1)

					drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

					! Right x
					phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0),i_zone)
					n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
					phinrx = mass * phirx/n_denrx

					! Left x
					philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0),i_zone)
					n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
					phinlx = mass * philx/n_denlx

					! -----------------------------------------------------

					drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

					! Right z
					phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu),i_zone)
					n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
					phinrz = mass * phirz/n_denrz

					! Left z
					philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid),i_zone)
					n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
					phinlz = mass * philz/n_denlz

					term0_5 = phi_TF_ofpsi(big_Psi0,i_zone) * (  &
								( phinrx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
								+ phinlx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
								/(x*dx2) &
								+ ( phinrz*(big_Psiu-big_Psi0) &
								+ phinlz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
								)

					term2 = eV* Fstar * phi_TF_ofpsi(big_Psi0,i_zone) / x**2
					term3 = n_denc * dh_i_dpsi(big_Psi0,i_zone)
					term4 = mass * n_denc**gamma * ds_i_dpsi(big_Psi0,i_zone)/(gamma-1.0d0)

					psi_diff_new = (-term0_5 + term2 + term3 - term4) / (ev**2*n_denc/mass/x**2)
				
				endif

				if((i==i_issue_check).and.(j==j_issue_check)) then
					continue
				endif

				call get_H_diff(psi(i,j),big_Psi(i,j),H_diff_loc,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
				H_loc = h_e_ofpsi_partial(psi(i,j),i_zone) + h_i_ofpsi_partial(big_Psi(i,j),i_zone) + eV*H_diff_loc	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				if(H_loc<=0.d0) then
					print*, 'problem in psi_diff max'
					continue
					psi_diff_limit = 0.d0
				else
					psi_diff_limit = sqrt(0.9d0*2.d0*mass*H_loc) * x/eV
				endif

				anorm_2 = anorm_2 + dabs(psi_diff_old - psi_diff_new)

				if(abs(psi_diff_new)>psi_diff_limit) then
					psi_diff_new = sign(psi_diff_limit,psi_diff_new)
				endif

				if(orp<=1.d0) then
					psi_diff(i,j) = psi_diff_old*(1.d0-orp) + psi_diff_new*orp
				else
					psi_diff(i,j) = psi_diff_new
				endif

				diff_error(i,j) = psi_diff_new - psi_diff_old

             enddo
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do		

	else

		! since psi_diff appears in the same way in two equations, we can get it from the GS one, which
		! is simpler, as it does not have derivatives of big_Psi and n_den
       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
          do i=isw+1, nx-1, 2

				! Added bc_type exclusion (NEW May 2025)
				if(sort_grid(i,j)<=0.and.bc_type/=7) cycle


				if(bc_type==7) then
					! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
					i_zone = min(0,sort_grid(i,j)-2)
					i_zone = max(-1,i_zone)
				endif

				! set up local psi values

				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi(i,j)
				if(i>1) psil = psi(i-1,j)
				if(i<nx) psir = psi(i+1,j)
				if(j>1) psid = psi(i,j-1)
				if(j<nz) psiu = psi(i,j+1)

				x = x_coord(i)

				psi_diff_old = psi_diff(i,j)

				n_denc = n_den(i,j)

				Fstar = Fstarofpsi(psi0,big_Psi(i,j),i_zone)

				if(i_zone==-1) then
					continue
				endif

				if(psi_diff_option==1) then
				! use magnetic GS for psi_diff

					term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
									+ (psil-psi0)/(x - 0.5d0*dx) )  &
									/(x*dx2) &
									+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

					term2 = -eV*mu_mag * Fstar * phi_e_TF_ofpsi(psi0,i_zone) / x**2
					term3 = mu_mag * n_denc * dh_e_dpsi(psi0,i_zone)
					term4 = mu_mag * mass * n_denc**gamma_e * ds_e_dpsi(psi0,i_zone)/(gamma_e-1.0d0)

					psi_diff_new = (-term0 - term2 - term3 + term4) / (ev**2*mu_mag*n_denc/mass/x**2)

				elseif(psi_diff_option==2) then
				! use fluid GS for psi_diff

					big_Psil = 0.d0
					big_Psir = 0.d0
					big_Psiu = 0.d0
					big_Psid = 0.d0

					big_Psi0 = big_Psi(i,j)
					if(i>1) big_Psil = big_Psi(i-1,j)
					if(i<nx) big_Psir = big_Psi(i+1,j)
					if(j>1) big_Psid = big_Psi(i,j-1)
					if(j<nz) big_Psiu = big_Psi(i,j+1)

					drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

					! Right x
					phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0),i_zone)
					n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
					phinrx = mass * phirx/n_denrx

					! Left x
					philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0),i_zone)
					n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
					phinlx = mass * philx/n_denlx

					! -----------------------------------------------------

					drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

					! Right z
					phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu),i_zone)
					n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
					phinrz = mass * phirz/n_denrz

					! Left z
					philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid),i_zone)
					n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
					phinlz = mass * philz/n_denlz

					term0_5 = phi_TF_ofpsi(big_Psi0,i_zone) * (  &
								( phinrx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
								+ phinlx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
								/(x*dx2) &
								+ ( phinrz*(big_Psiu-big_Psi0) &
								+ phinlz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
								)

					term2 = eV* Fstar * phi_TF_ofpsi(big_Psi0,i_zone) / x**2
					term3 = n_denc * dh_i_dpsi(big_Psi0,i_zone)
					term4 = mass * n_denc**gamma * ds_i_dpsi(big_Psi0,i_zone)/(gamma-1.0d0)

					psi_diff_new = (-term0_5 + term2 + term3 - term4) / (ev**2*n_denc/mass/x**2)
				
				endif

				if((i==i_issue_check).and.(j==j_issue_check)) then
					continue
				endif

				call get_H_diff(psi(i,j),big_Psi(i,j),H_diff_loc,i_zone) ! dbsval, or integrate hofpsi(psi(i,j))
				H_loc = h_e_ofpsi_partial(psi(i,j),i_zone) + h_i_ofpsi_partial(big_Psi(i,j),i_zone) + eV*H_diff_loc	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

				if(H_loc<=0.d0) then
					print*, 'problem in psi_diff max'
					continue
					psi_diff_limit = 0.d0
				else
					psi_diff_limit = sqrt(0.9d0*2.d0*mass*H_loc) * x/eV
				endif

				anorm_2 = anorm_2 + dabs(psi_diff_old - psi_diff_new)

				if(abs(psi_diff_new)>psi_diff_limit) then
					psi_diff_new = sign(psi_diff_limit,psi_diff_new)
				endif

				if(orp<=1.d0) then
					psi_diff(i,j) = psi_diff_old*(1.d0-orp) + psi_diff_new*orp
				else
					psi_diff(i,j) = psi_diff_new
				endif

				diff_error(i,j) = psi_diff_new - psi_diff_old

             enddo
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do
	
	endif

	mx = 0.d0

	! If we use the GS equation for psi_diff, we might as well update big_Psi immediately, but let's keep things separated for clarity
	! and update it here.
	do j = 1, nz
	do i = 1, nx

			! Added bc_type exclusion (NEW May 2025)
			if(sort_grid(i,j)<=0.and.bc_type/=7) cycle

			big_Psi(i,j) = psi(i,j) + psi_diff(i,j)
			if((bc_type/=7).or.((bc_type==7).and.(sort_grid(i,j)==2))) then
				if(dabs(big_Psi(i,j))>mx) then
					mx = dabs(big_Psi(i,j))
					! Comment by Ian: same restrictive if statement, but this time with positive LCFS check 
					if((tri_type==-4).and.(LCFS==1)) then
						ipsic = i
						jpsic = j
					endif
				endif
			endif

	enddo
	enddo
	!----------------end of big_Psi cycle-----------------

	! Set the new maximum value of big_Psic
	big_Psic = mx
	! "Finished: big_Psi Update"

	! CALL BC HERE AGAIN FOR BIG_PSI
	! Comment by Ian: again, not sure if all these BC calls are necessary?
	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,2)

	! -----------------------------------------------------

	psi_max = max(psic, big_Psic)
	psic_13 = psi_max

	! Move the density calculation to a separate function
	! Update density and seek mach theta max

!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)
!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)

!!	if(Broot==0) then
!!		call update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	else
		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
!!	endif

	! Comment by Ian: another BC call that I'm not sure is necessary
	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)
	! here update just n_den

	if(in_orp <= 0.0d0) then
		! Chebyshev acceleration 
		if(nx >= 5) then
			std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
		end if
	endif

	if ((accelerate).and.(nx<inter_switch)) then
		orp = std_orp
	else
		orp = fix_orp
	endif

	!!$ call step_output(nx,psi(:,:),rho,residual)

	! -----------------------------------------------------

	if(k == 1) then
		eps_anorm = eps*anorm
		eps_anorm_2 = eps*anorm_2
	else
		if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
			! As we get closer to convergence, we want the solution
			! to relax.  So as anorm approaches eps_anorm we want 
			! the over relaxation parameter to go to some const. ~ 1
			! Use x and mtm_soln as a temporary variable
			x = anorm/eps_anorm
			mtm_soln = 1.0d0 ! The limiting value for x ~ 1
			tmp = x/(x + orp/mtm_soln - 1.0d0)
			tmp = dmin1(tmp,1.0d0)
			orp = orp*tmp
			orp = dmax1(orp,mtm_soln)
		endif
	end if

	if(k == 1 .or. modulo(k,k_save) == 0) then
		print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		print *, k,": anorm_2 = ",real(anorm_2,skind)," eps_anorm_2 = ",real(eps_anorm_2,skind)
		print *, "The Over-Relaxation Parameter =",orp,std_orp
		call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)
		! Comment by Ian: first 2 layers of if statements here are reasonable, 
		! but the other conditions seem overly restrictive or just plain weird

		! Added inorm to update_sort_grid call (NEW May 2025)
		if(bc_type==7) then
			if(LCFS==-1) then
				if((tri_type==-3)) then
					call update_sort_grid(psi(:,:),nx,nz,inorm)
				elseif((k>25).and.(tri_type==-2)) then
					call update_interface(psi(:,:),n,inorm)
				endif			
			elseif(LCFS==1) then
				if((tri_type==-3)) then
					call update_sort_grid(big_Psi(:,:),nx,nz,inorm)
				elseif((k>25).and.(tri_type==-2)) then
					call update_interface(big_Psi(:,:),nx,inorm)
				endif			
			endif
		endif
		continue
	end if

	write(111,*) k,anorm
	write(112,*) k,anorm_2

	! Check the convergence criteria
	if((anorm < eps_anorm) .and. (anorm_2 < eps_anorm_2).and.(inorm<1.d-5) .and. (k > min_it)) exit

	end do

	print *, k,": anorm = ", anorm, " eps_anorm = ",eps_anorm
	print *, k,": anorm_2 = ", anorm_2, " eps_anorm_2 = ", eps_anorm_2
	print *, "Average residual =", anorm/(nx*nz)
	print *, "Average difference (2) =", anorm_2/(nx*nz)

	print *, "Mach Theta Max =", mach_theta_max
	print *, "Final Solution has psi Center = ", psic
	print *, "and big_Psi Center = ", big_Psic
	print *

	return

end subroutine ngs_solve_TF


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz], big_Psi is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this is the two-fluid version of the Gauss-modified Bernoulli equation
! the transonic surface is a flow surface
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_TF_gauss_flow(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! cycle separately on psi and psi_diff
! IMPORTANT: solve directly for psi_diff, then calculate big_Psi as psi + psi_diff

	use constant, only : dx, dz, dx_a, dz_a

	implicit none

	integer, intent(in) :: nx,nz,min_it,max_it,flag
	real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,big_Psi,psi_diff,n_den,residual
	real (kind=dkind), intent(in) :: eps
	! Input Over Relaxation Parameter
	! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
	! else we set orp = in_orp
	real (kind=dkind), intent(in) :: in_orp
	integer ::ipass,i,isw,j,jsw,k,alloc_stat
	real (kind=dkind) :: res, den, x !, xl, xr
	! res -> residual
	! den -> denominator of the Newton-Gauss-Seidel update
	! x -> Radial position of grid point (i,j)
	real (kind=dkind) :: dx2,dz2,mx
	real (kind=dkind) :: anorm, eps_anorm, anorm_2, eps_anorm_2
	real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS, grad_bigPsi_RHS
	real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
	real (kind=dkind) :: orp, std_orp ! over relaxation parameter
	real (kind=dkind) :: n_denc, phic
	! by is the phi component of the magnetic field
	real (kind=dkind) :: Fstar
	real (kind=dkind) :: last_mtm ! last mach theta max
	! minimum delta rho for bernoulli roots, the mean val. and location
	real (kind=dkind) :: min_dn_den, tmp
	integer :: min_ix, min_iz
	! The next set of variables are used for the bisection search for
	! locating the degenerate root for Mach Theta Max
	real (kind=dkind) :: mtm_soln
	! mtm_acc controls the bisection loop seeking mach theta max
	real (kind=dkind), parameter :: mtm_acc = 1.0d-12
	! VERY IMPORTANT
	! psi_degen records the value of psi where the solution to the
	! Bernoulli equation is approximately degenerate. 
	real (kind=dkind) :: gPsi2
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term0_5, term6, term6a, term6b
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz, phinc
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	real(kind=dkind), dimension(-1:1,-1:1) :: big_Psi_around
	real(kind=dkind) :: psi_diff_old, psi_diff_new, psi_diff_limit
	real(kind=dkind) :: H_loc, H_diff_loc
	real(kind=dkind) :: dpsidx, dpsidz, drho_ds
!	integer :: psi_diff_option = 2
	real (kind=dkind), dimension(1:nx,1:nz) :: diff_error
	integer :: k_save = 25
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92

	psi_max = max(psic,big_Psic)

    if(psi_degen==0.d0) then 
		psi_degen = 0.5d0*psic ! Initialization
		big_Psi_degen = 0.5d0*psic ! Initialization
	endif
	last_mtm = mach_theta_max
	eps_anorm = 0.0d0
	eps_anorm_2 = 0.0d0

	fix_orp0 = fix_orp

	! The under/over relaxation parameter
	! Standard Case:
	if(in_orp <= 0.0d0) then
		! Note that the parameters below are tuned for starting with 
		! a 17 x 17 grid.
		orp = 1.0d0
		if(nx <= 5) then
			orp = 0.5d0
			rjac = 0.75d0
		else if(nx <= 9) then
			orp = 0.5d0
			rjac = 0.9d0
		else if(nx <= 17) then
			rjac = 0.98d0
		else if(nx <= 33) then
			rjac = 0.995d0
		else if(nx <= 65) then
			rjac = 0.9972d0
		else
		! It would appear that 129 is about where the solution begins to converge
			rjac = 0.0d0
		end if
		std_orp = orp
	else
		print *, "Input Over Relaxation Parameter =",in_orp
		orp = in_orp
		std_orp = in_orp
	endif

	if (accelerate) then
		continue
	else
		orp = fix_orp
	endif

	dx2 = dx*dx
	dz2 = dz*dz

	diff_error = 0.d0

	if(allocated(fmax_2D)) then
		deallocate(fmax_2D)
	endif

	allocate(fmax_2D(1:nx,1:nz))
	fmax_2D = 0.d0

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,0)

  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	if(nx>=inter_switch) then
		call GS_RHS_all_Gauss_TF
	endif

	if((nx>=inter_switch).and.(Broot==0)) then
		call update_rho_TF_gauss_2020_flow(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
	else
		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)

	call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)

	! Iterate until convergence
	do k=1, max_it
	!     "Update psi"6

	!		call cpu_time(tic)

		if((k>ite_Broot).and.(Broot/=Broot_true)) then
			Broot = Broot_true
		endif

	! Reset the norm and max of psi and big_Psi
	mx = 0.0d0
	anorm = 0.0d0
	anorm_2 = 0.0d0

	do j=2,nz-1
	do i=2,nx-1
		if(sort_grid(i,j)>0) psi(i,j) = dabs(psi(i,j))	!	segno
	enddo
	enddo

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	!----------------------psi cycle------------------------
	jsw = 1
	do ipass = 1, 2
		isw=jsw;
		do j=2, nz-1
		do i=isw+1, nx-1, 2
		! "Update psi(",i,",",j,")"
		! Only solve the inner region problem

			if(sort_grid(i,j)<=0) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			big_Psil = 0.d0
			big_Psir = 0.d0
			big_Psiu = 0.d0
			big_Psid = 0.d0

			big_Psi0 = big_Psi(i,j)
			if(i>1) big_Psil = big_Psi(i-1,j)
			if(i<nx) big_Psir = big_Psi(i+1,j)
			if(j>1) big_Psid = big_Psi(i,j-1)
			if(j<nz) big_Psiu = big_Psi(i,j+1)

			x = x_coord(i)

			! Calculate rho and phi at the current location
			n_denc = n_den(i,j)
			phic = phi_TF_ofpsi(big_Psi0)

			! old formulation with new Bernoulli

			drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

			! Right x
			phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
			n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
			ma2rx = mu_mag*mass * phirx*phirx/n_denrx

			! Left x
			philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
			n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
			ma2lx = mu_mag*mass * philx*philx/n_denlx

			! -----------------------------------------------------

			drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

			! Right z
			phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
			n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
			ma2rz = mu_mag*mass * phirz*phirz/n_denrz

			! Left z
			philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
			n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
			ma2lz = mu_mag*mass * philz*philz/n_denlz

			! -----------------------------------------------------

			big_Psi_around(-1,0) = big_Psil
			big_Psi_around(1,0) = big_Psir
			big_Psi_around(0,0) = big_Psi0
			big_Psi_around(0,-1) = big_Psid
			big_Psi_around(0,1) = big_Psiu

			Fstar = Fstarofpsi(psi0,big_Psi0)

			call psi_derivative_new(i,j,nx,nz,big_Psi_around,dpsidx,dpsidz)
			gPsi2 =  dpsidx**2+dpsidz**2

			! we follow FLOW nomenclature for term names (everything should be easily recognizable)

			term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
							+ (psil-psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

			term0_5 = - (  &
							( ma2rx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ ma2lx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( ma2rz*(big_Psiu-big_Psi0) &
							+ ma2lz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

			term1 = eV*mu_mag * Fstar *(phic-phi_e_TF_ofpsi(psi0))/x**2
			term3 = mu_mag*mass * phic * dphi_TF_dpsi(big_Psi0) * gPsi2 / n_denc/x**2

			term2 = mu_mag * eV * n_denc * (omegaofpsi(big_Psi0)-omegaofpsi(psi0)) -  &
					mu_mag * n_denc * (dp_i_dpsi(big_Psi0)/d_TF_ofpsi(big_Psi0)-dp_i_dpsi(psi0)/d_TF_ofpsi(psi0))

			term4 = mu_mag * n_denc * (dh_e_dpsi_partial(psi0) + dh_i_dpsi_partial(big_psi0))

			term5 = mu_mag*mass * (n_denc**gamma_i/(gamma_i-1.0d0) * ds_i_dpsi(big_Psi0) +  &
												n_denc**gamma_e/(gamma_e-1.0d0) * ds_e_dpsi(psi0))
			! CHECK DIMENSIONS FOR ALL!!!

			 term6a = -(grad_bigPsi_RHS(1,i,j)*(fmax_2D(i+1,j)-fmax_2D(i-1,j))/dx/2.d0 +  &
					grad_bigPsi_RHS(2,i,j)*(fmax_2D(i,j+1)-fmax_2D(i,j-1))/dz/2.d0)/grad_bigPsi_RHS(0,i,j)**2 +  &
					2.d0*delta_Bern*(big_Psi(i,j)-big_Psi_degen)/psi_max**2*fmax_2D(i,j)

			 term6b =mu_mag * n_denc * mass * exp(-delta_Bern*((big_Psi(i,j)-big_Psi_degen)/psi_max)**2)

			 term6 = term6b*term6a

			 if(k<50) then
				term6 = term6/(51.d0-k)**2
			endif

			!				if((nx>=inter_switch).and.(term6/=0.d0)) then
			!					print*, i,j, term0, term1, term2, term3, term4, term5, term6
			!					print*, '   '
			!					pause
			!				endif

			res = term0 + term0_5 + term1 + term2 + term3 + term4 - term5 + term6

			den = -( 1.0d0/(x + 0.5d0*dx) + &
						1.0d0/(x - 0.5d0*dx) )/(x*dx2) &
						-2.0d0/(x*x*dz2)

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

			continue

			!res = res*mu_mag

			! Store the residual
			if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
				if(nx>=inter_switch) then
					residual(i,j) = res/den
				else
					residual(i,j) = res
				endif
			end if

			if(res>1.d0) then
				continue
			endif

			! -----------------------------------------------------

			continue

			psi(i,j) = psi0 - orp*res/den

			! Find the max absolute value of psi in the plasma
			mx = dmax1(dabs(psi(i,j)),mx)

			anorm = anorm + dabs(res)

		end do
		isw = 3 - isw
		end do
		jsw = 3 - jsw
		call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,1)
		! only update psi on BC
	end do
	!----------------end of psi cycle--------------

	! Set the new maximum value of psic
	psic = mx
	psi_max = max(psic,big_Psic)
	psic_13 = psi_max
	! "Finished: psi Update"

	!----------------big_Psi cycle-----------------
	! we use psi_diff as dependent variable, and solve the algebraic equation for it
	! FIRST ATTEMPT: USE OLD ITERATION FOR DERIVATIVES

	! since psi_diff appears in the same way in two equations, we can get it from the GS one, which
	! is simpler, as it does not have derivatives of big_Psi and n_den
	do j = 1, nz
	do i = 1, nx

			if(sort_grid(i,j)<=0) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			x = x_coord(i)

			psi_diff_old = psi_diff(i,j)

			n_denc = n_den(i,j)

			Fstar = Fstarofpsi(psi0,big_Psi(i,j))

			if(psi_diff_option==1) then
			! use magnetic GS for psi_diff

				term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
								+ (psil-psi0)/(x - 0.5d0*dx) )  &
								/(x*dx2) &
								+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

				term2 = -eV*mu_mag * Fstar * phi_e_TF_ofpsi(psi0) / x**2
				term3 = mu_mag * n_denc * dh_e_dpsi(psi0)
				term4 = mu_mag * mass * n_denc**gamma_e * ds_e_dpsi(psi0)/(gamma_e-1.0d0)

				psi_diff_new = (-term0 - term2 - term3 + term4) / (ev**2*mu_mag*n_denc/mass/x**2)

			elseif(psi_diff_option==2) then
			! use fluid GS for psi_diff

				big_Psil = 0.d0
				big_Psir = 0.d0
				big_Psiu = 0.d0
				big_Psid = 0.d0

				big_Psi0 = big_Psi(i,j)
				if(i>1) big_Psil = big_Psi(i-1,j)
				if(i<nx) big_Psir = big_Psi(i+1,j)
				if(j>1) big_Psid = big_Psi(i,j-1)
				if(j<nz) big_Psiu = big_Psi(i,j+1)

				drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

				! Right x
				phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
				n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
				phinrx = mass * phirx/n_denrx

				! Left x
				philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
				n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
				phinlx = mass * philx/n_denlx

				! -----------------------------------------------------

				drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

				! Right z
				phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
				n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
				phinrz = mass * phirz/n_denrz

				! Left z
				philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
				n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
				phinlz = mass * philz/n_denlz

				term0_5 = phi_TF_ofpsi(big_Psi0) * (  &
							( phinrx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ phinlx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( phinrz*(big_Psiu-big_Psi0) &
							+ phinlz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

				term2 = eV* Fstar * phi_TF_ofpsi(big_Psi0) / x**2
				term3 = n_denc * dh_i_dpsi(big_Psi0)
				term4 = mass * n_denc**gamma_i * ds_i_dpsi(big_Psi0)/(gamma_i-1.0d0)

				psi_diff_new = (-term0_5 + term2 + term3 - term4) / (ev**2*n_denc/mass/x**2)
				
			endif

			call get_H_diff(psi(i,j),big_Psi(i,j),H_diff_loc) ! dbsval, or integrate hofpsi(psi(i,j))
			H_loc = h_e_ofpsi_partial(psi(i,j)) + h_i_ofpsi_partial(big_Psi(i,j)) + eV*H_diff_loc	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

			if(H_loc<=0.d0) then
				print*, 'problem in psi_diff max'
				continue
				psi_diff_limit = 0.d0
			else
				psi_diff_limit = sqrt(0.9d0*2.d0*mass*H_loc) * x/eV
			endif

			anorm_2 = anorm_2 + dabs(psi_diff_old - psi_diff_new)

			if(abs(psi_diff_new)>psi_diff_limit) then
				psi_diff_new = sign(psi_diff_limit,psi_diff_new)
			endif

			if(orp<=1.d0) then
				psi_diff(i,j) = psi_diff_old*(1.d0-orp) + psi_diff_new*orp
			else
				psi_diff(i,j) = psi_diff_new
			endif

			diff_error(i,j) = psi_diff_new - psi_diff_old

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

	enddo
	enddo

	mx = 0.d0

	! if we use the GS equation for psi_diff, we might as well update big_Psi immediately, but let's keep things separated for clarity
	do j = 1, nz
	do i = 1, nx

			if(sort_grid(i,j)<=0) cycle

			big_Psi(i,j) = psi(i,j) + psi_diff(i,j)
			mx = dmax1(big_Psi(i,j),mx)

	enddo
	enddo
	!----------------end of big_Psi cycle-----------------

	! Set the new maximum value of big_Psic
	big_Psic = mx
	! "Finished: big_Psi Update"

	! CALL BC HERE AGAIN FOR BIG_PSI
	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,2)

	! -----------------------------------------------------

	psi_max = max(psic, big_Psic)
	psic_13 = psi_max

     ! update RHS
     	if(nx>=inter_switch) then
			call GS_RHS_all_Gauss_TF
		endif


	! Move the density calculation to a separate function
	! Update density and seek mach theta max

!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)
!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)

!!	if(Broot==0) then
!!		call update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	else
!!		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	endif

	if((nx>=inter_switch).and.(Broot==0)) then
!!$		if(k<50) then
!!$			call update_rho_TF_gauss(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!$								min_ix,min_iz)
		! changed March 1 2020
		! added attempt to use both on March 3 2020
!!$		else
			call update_rho_TF_gauss_2020_flow(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
								min_ix,min_iz)
!!$		endif
	else
		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)
	! here update just n_den

	if(in_orp <= 0.0d0) then
		! Chebyshev acceleration 
		if(nx >= 5) then
			std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
		end if
	endif

	if ((accelerate).and.(nx<inter_switch)) then
		orp = std_orp
	else
		orp = fix_orp
	endif

	!!$ call step_output(nx,psi(:,:),rho,residual)

	! -----------------------------------------------------

	if(k == 1) then
		eps_anorm = eps*anorm
		eps_anorm_2 = eps*anorm_2
	else
		if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
			! As we get closer to convergence, we want the solution
			! to relax.  So as anorm approaches eps_anorm we want 
			! the over relaxation parameter to go to some const. ~ 1
			! Use x and mtm_soln as a temporary variable
			x = anorm/eps_anorm
			mtm_soln = 1.0d0 ! The limiting value for x ~ 1
			tmp = x/(x + orp/mtm_soln - 1.0d0)
			tmp = dmin1(tmp,1.0d0)
			orp = orp*tmp
			orp = dmax1(orp,mtm_soln)
		endif
	end if

	if(k == 1 .or. modulo(k,k_save) == 0) then
		print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		print *, k,": anorm_2 = ",real(anorm_2,skind)," eps_anorm_2 = ",real(eps_anorm_2,skind)
		print *, "The Over-Relaxation Parameter =",orp,std_orp
		call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)
		continue
	end if

	write(111,*) k,anorm
	write(112,*) k,anorm_2

	! Check the convergence criteria
	if((anorm < eps_anorm) .and. (anorm_2 < eps_anorm_2) .and. (k > min_it)) exit

	end do

	print *, k,": anorm = ", anorm, " eps_anorm = ",eps_anorm
	print *, k,": anorm_2 = ", anorm_2, " eps_anorm_2 = ", eps_anorm_2
	print *, "Average residual =", anorm/(nx*nz)
	print *, "Average difference (2) =", anorm_2/(nx*nz)

	print *, "Mach Theta Max =", mach_theta_max
	print *, "Final Solution has psi Center = ", psic
	print *, "and big_Psi Center = ", big_Psic
	print *

	return

	contains


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS_all_Gauss_TF
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind) :: tdx, tdz


    tdx = dx*2.d0
    tdz = dz*2.d0

    gradpsi_RHS = 1.d-10
    grad_bigPsi_RHS = 1.d-10

    ! fill grad psi
    do j = 1, nz
       do i = 1, nx

			!-------------------psi-----------------------
          if(sort_grid(i,j)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi(i+1,j)-psi(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi(i,j+1)-psi(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)

			!-------------------big_Psi-----------------------
          grad_bigPsi_RHS(1,i,j) = (big_Psi(i+1,j)-big_Psi(i-1,j))/tdx
          grad_bigPsi_RHS(2,i,j) = (big_Psi(i,j+1)-big_Psi(i,j-1))/tdz
          grad_bigPsi_RHS(0,i,j) = sqrt(grad_bigPsi_RHS(1,i,j)**2+grad_bigPsi_RHS(2,i,j)**2)

          if(grad_bigPsi_RHS(0,i,j)<1.d-6) grad_bigPsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)


       enddo
    enddo


!	if((k==1).or.(modulo(k,k_save)==0)) then
!
!    ! write results for debugging
!
!		open(13, file='gradpsi_RHS.plt')
!
!		write(13,*)'TITLE="solution of GS equation with flow"'
!		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
!		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'
!
!		do j=1,nz
!
!		   z = z_coord(j)
!
!		   do i=1,nx
!
!			  x = x_coord(i)
!
!			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
!
!		   end do
!		end do
!
!		close(13)
!
!	endif

90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS_all_Gauss_TF

end subroutine ngs_solve_TF_gauss_flow


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz], big_Psi is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this is the two-fluid version of the Gauss-modified Bernoulli equation
! the transonic surface is a magnetic surface
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_TF_gauss_magnetic(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! cycle separately on psi and psi_diff
! IMPORTANT: solve directly for psi_diff, then calculate big_Psi as psi + psi_diff

	use constant, only : dx, dz, dx_a, dz_a

	implicit none

	integer, intent(in) :: nx,nz,min_it,max_it,flag
	real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,big_Psi,psi_diff,n_den,residual
	real (kind=dkind), intent(in) :: eps
	! Input Over Relaxation Parameter
	! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
	! else we set orp = in_orp
	real (kind=dkind), intent(in) :: in_orp
	integer ::ipass,i,isw,j,jsw,k,alloc_stat
	real (kind=dkind) :: res, den, x !, xl, xr
	! res -> residual
	! den -> denominator of the Newton-Gauss-Seidel update
	! x -> Radial position of grid point (i,j)
	real (kind=dkind) :: dx2,dz2,mx
	real (kind=dkind) :: anorm, eps_anorm, anorm_2, eps_anorm_2
	real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS, grad_bigPsi_RHS
	real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
	real (kind=dkind) :: orp, std_orp ! over relaxation parameter
	real (kind=dkind) :: n_denc, phic
	! by is the phi component of the magnetic field
	real (kind=dkind) :: Fstar
	real (kind=dkind) :: last_mtm ! last mach theta max
	! minimum delta rho for bernoulli roots, the mean val. and location
	real (kind=dkind) :: min_dn_den, tmp
	integer :: min_ix, min_iz
	! The next set of variables are used for the bisection search for
	! locating the degenerate root for Mach Theta Max
	real (kind=dkind) :: mtm_soln
	! mtm_acc controls the bisection loop seeking mach theta max
	real (kind=dkind), parameter :: mtm_acc = 1.0d-12
	! VERY IMPORTANT
	! psi_degen records the value of psi where the solution to the
	! Bernoulli equation is approximately degenerate.
	real (kind=dkind) :: gPsi2
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term0_5, term6, term6a, term6b
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz, phinc
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	real(kind=dkind), dimension(-1:1,-1:1) :: big_Psi_around
	real(kind=dkind) :: psi_diff_old, psi_diff_new, psi_diff_limit
	real(kind=dkind) :: H_loc, H_diff_loc
	real(kind=dkind) :: dpsidx, dpsidz, drho_ds
!	integer :: psi_diff_option = 2
	real (kind=dkind), dimension(1:nx,1:nz) :: diff_error
	integer :: k_save = 25
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92

	psi_max = max(psic,big_Psic)

    if(psi_degen==0.d0) then 
		psi_degen = 0.5d0*psic ! Initialization
		big_Psi_degen = 0.5d0*psic ! Initialization
	endif
	last_mtm = mach_theta_max
	eps_anorm = 0.0d0
	eps_anorm_2 = 0.0d0

	fix_orp0 = fix_orp

	! The under/over relaxation parameter
	! Standard Case:
	if(in_orp <= 0.0d0) then
		! Note that the parameters below are tuned for starting with 
		! a 17 x 17 grid.
		orp = 1.0d0
		if(nx <= 5) then
			orp = 0.5d0
			rjac = 0.75d0
		else if(nx <= 9) then
			orp = 0.5d0
			rjac = 0.9d0
		else if(nx <= 17) then
			rjac = 0.98d0
		else if(nx <= 33) then
			rjac = 0.995d0
		else if(nx <= 65) then
			rjac = 0.9972d0
		else
		! It would appear that 129 is about where the solution begins to converge
			rjac = 0.0d0
		end if
		std_orp = orp
	else
		print *, "Input Over Relaxation Parameter =",in_orp
		orp = in_orp
		std_orp = in_orp
	endif

	if (accelerate) then
		continue
	else
		orp = fix_orp
	endif

	dx2 = dx*dx
	dz2 = dz*dz

	diff_error = 0.d0

	if(allocated(fmax_2D)) then
		deallocate(fmax_2D)
	endif

	allocate(fmax_2D(1:nx,1:nz))
	fmax_2D = 0.d0

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,0)

  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	if(nx>=inter_switch) then
		call GS_RHS_all_Gauss_TF
	endif

	if((nx>=inter_switch).and.(Broot==0)) then
		call update_rho_TF_gauss_2020_magnetic(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
	else
		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)

	call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)

	! Iterate until convergence
	do k=1, max_it
	!     "Update psi"6

	!		call cpu_time(tic)

		if((k>ite_Broot).and.(Broot/=Broot_true)) then
			Broot = Broot_true
		endif

	! Reset the norm and max of psi and big_Psi
	mx = 0.0d0
	anorm = 0.0d0
	anorm_2 = 0.0d0

	do j=2,nz-1
	do i=2,nx-1
		if(sort_grid(i,j)>0) psi(i,j) = dabs(psi(i,j))	!	segno
	enddo
	enddo

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	!----------------------psi cycle------------------------
	jsw = 1
	do ipass = 1, 2
		isw=jsw;
		do j=2, nz-1
		do i=isw+1, nx-1, 2
		! "Update psi(",i,",",j,")"
		! Only solve the inner region problem

			if(sort_grid(i,j)<=0) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			big_Psil = 0.d0
			big_Psir = 0.d0
			big_Psiu = 0.d0
			big_Psid = 0.d0

			big_Psi0 = big_Psi(i,j)
			if(i>1) big_Psil = big_Psi(i-1,j)
			if(i<nx) big_Psir = big_Psi(i+1,j)
			if(j>1) big_Psid = big_Psi(i,j-1)
			if(j<nz) big_Psiu = big_Psi(i,j+1)

			x = x_coord(i)

			! Calculate rho and phi at the current location
			n_denc = n_den(i,j)
			phic = phi_TF_ofpsi(big_Psi0)

			! old formulation with new Bernoulli

			drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

			! Right x
			phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
			n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
			ma2rx = mu_mag*mass * phirx*phirx/n_denrx

			! Left x
			philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
			n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
			ma2lx = mu_mag*mass * philx*philx/n_denlx

			! -----------------------------------------------------

			drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

			! Right z
			phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
			n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
			ma2rz = mu_mag*mass * phirz*phirz/n_denrz

			! Left z
			philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
			n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
			ma2lz = mu_mag*mass * philz*philz/n_denlz

			! -----------------------------------------------------

			big_Psi_around(-1,0) = big_Psil
			big_Psi_around(1,0) = big_Psir
			big_Psi_around(0,0) = big_Psi0
			big_Psi_around(0,-1) = big_Psid
			big_Psi_around(0,1) = big_Psiu

			Fstar = Fstarofpsi(psi0,big_Psi0)

			call psi_derivative_new(i,j,nx,nz,big_Psi_around,dpsidx,dpsidz)
			gPsi2 =  dpsidx**2+dpsidz**2

			! we follow FLOW nomenclature for term names (everything should be easily recognizable)

			term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
							+ (psil-psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

			term0_5 = - (  &
							( ma2rx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ ma2lx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( ma2rz*(big_Psiu-big_Psi0) &
							+ ma2lz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

			term1 = eV*mu_mag * Fstar *(phic-phi_e_TF_ofpsi(psi0))/x**2
			term3 = mu_mag*mass * phic * dphi_TF_dpsi(big_Psi0) * gPsi2 / n_denc/x**2

			term2 = mu_mag * eV * n_denc * (omegaofpsi(big_Psi0)-omegaofpsi(psi0)) -  &
					mu_mag * n_denc * (dp_i_dpsi(big_Psi0)/d_TF_ofpsi(big_Psi0)-dp_i_dpsi(psi0)/d_TF_ofpsi(psi0))

			term4 = mu_mag * n_denc * (dh_e_dpsi_partial(psi0) + dh_i_dpsi_partial(big_psi0))

			term5 = mu_mag*mass * (n_denc**gamma_i/(gamma_i-1.0d0) * ds_i_dpsi(big_Psi0) +  &
												n_denc**gamma_e/(gamma_e-1.0d0) * ds_e_dpsi(psi0))
			! CHECK DIMENSIONS FOR ALL!!!

			 term6a = -(grad_bigPsi_RHS(1,i,j)*(fmax_2D(i+1,j)-fmax_2D(i-1,j))/dx/2.d0 +  &
					grad_bigPsi_RHS(2,i,j)*(fmax_2D(i,j+1)-fmax_2D(i,j-1))/dz/2.d0)/grad_bigPsi_RHS(0,i,j)**2 +  &
					2.d0*delta_Bern*(big_Psi(i,j)-big_Psi_degen)/psi_max**2*fmax_2D(i,j)

			 term6b =mu_mag * n_denc * mass * exp(-delta_Bern*((big_Psi(i,j)-big_Psi_degen)/psi_max)**2)

			 term6 = term6b*term6a

			 if(k<50) then
				term6 = term6/(51.d0-k)**2
			endif

			!				if((nx>=inter_switch).and.(term6/=0.d0)) then
			!					print*, i,j, term0, term1, term2, term3, term4, term5, term6
			!					print*, '   '
			!					pause
			!				endif

			res = term0 + term0_5 + term1 + term2 + term3 + term4 - term5 + term6

			den = -( 1.0d0/(x + 0.5d0*dx) + &
						1.0d0/(x - 0.5d0*dx) )/(x*dx2) &
						-2.0d0/(x*x*dz2)

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

			continue

			!res = res*mu_mag

			! Store the residual
			if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
				if(nx>=inter_switch) then
					residual(i,j) = res/den
				else
					residual(i,j) = res
				endif
			end if

			if(res>1.d0) then
				continue
			endif

			! -----------------------------------------------------

			continue

			psi(i,j) = psi0 - orp*res/den

			! Find the max absolute value of psi in the plasma
			mx = dmax1(dabs(psi(i,j)),mx)

			anorm = anorm + dabs(res)

		end do
		isw = 3 - isw
		end do
		jsw = 3 - jsw
		call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,1)
		! only update psi on BC
	end do
	!----------------end of psi cycle--------------

	! Set the new maximum value of psic
	psic = mx
	psi_max = max(psic,big_Psic)
	psic_13 = psi_max
	! "Finished: psi Update"

	!----------------big_Psi cycle-----------------
	! we use psi_diff as dependent variable, and solve the algebraic equation for it
	! FIRST ATTEMPT: USE OLD ITERATION FOR DERIVATIVES

	! since psi_diff appears in the same way in two equations, we can get it from the GS one, which
	! is simpler, as it does not have derivatives of big_Psi and n_den
	do j = 1, nz
	do i = 1, nx

			if(sort_grid(i,j)<=0) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			x = x_coord(i)

			psi_diff_old = psi_diff(i,j)

			n_denc = n_den(i,j)

			Fstar = Fstarofpsi(psi0,big_Psi(i,j))

			if(psi_diff_option==1) then
			! use magnetic GS for psi_diff

				term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
								+ (psil-psi0)/(x - 0.5d0*dx) )  &
								/(x*dx2) &
								+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

				term2 = -eV*mu_mag * Fstar * phi_e_TF_ofpsi(psi0) / x**2
				term3 = mu_mag * n_denc * dh_e_dpsi(psi0)
				term4 = mu_mag * mass * n_denc**gamma_e * ds_e_dpsi(psi0)/(gamma_e-1.0d0)

				psi_diff_new = (-term0 - term2 - term3 + term4) / (ev**2*mu_mag*n_denc/mass/x**2)

			elseif(psi_diff_option==2) then
			! use fluid GS for psi_diff

				big_Psil = 0.d0
				big_Psir = 0.d0
				big_Psiu = 0.d0
				big_Psid = 0.d0

				big_Psi0 = big_Psi(i,j)
				if(i>1) big_Psil = big_Psi(i-1,j)
				if(i<nx) big_Psir = big_Psi(i+1,j)
				if(j>1) big_Psid = big_Psi(i,j-1)
				if(j<nz) big_Psiu = big_Psi(i,j+1)

				drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

				! Right x
				phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
				n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
				phinrx = mass * phirx/n_denrx

				! Left x
				philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
				n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
				phinlx = mass * philx/n_denlx

				! -----------------------------------------------------

				drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

				! Right z
				phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
				n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
				phinrz = mass * phirz/n_denrz

				! Left z
				philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
				n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
				phinlz = mass * philz/n_denlz

				term0_5 = phi_TF_ofpsi(big_Psi0) * (  &
							( phinrx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ phinlx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( phinrz*(big_Psiu-big_Psi0) &
							+ phinlz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

				term2 = eV* Fstar * phi_TF_ofpsi(big_Psi0) / x**2
				term3 = n_denc * dh_i_dpsi(big_Psi0)
				term4 = mass * n_denc**gamma_i * ds_i_dpsi(big_Psi0)/(gamma_i-1.0d0)

				psi_diff_new = (-term0_5 + term2 + term3 - term4) / (ev**2*n_denc/mass/x**2)
				
			endif

			call get_H_diff(psi(i,j),big_Psi(i,j),H_diff_loc) ! dbsval, or integrate hofpsi(psi(i,j))
			H_loc = h_e_ofpsi_partial(psi(i,j)) + h_i_ofpsi_partial(big_Psi(i,j)) + eV*H_diff_loc	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

			if(H_loc<=0.d0) then
				print*, 'problem in psi_diff max'
				continue
				psi_diff_limit = 0.d0
			else
				psi_diff_limit = sqrt(0.9d0*2.d0*mass*H_loc) * x/eV
			endif

			anorm_2 = anorm_2 + dabs(psi_diff_old - psi_diff_new)

			if(abs(psi_diff_new)>psi_diff_limit) then
				psi_diff_new = sign(psi_diff_limit,psi_diff_new)
			endif

			if(orp<=1.d0) then
				psi_diff(i,j) = psi_diff_old*(1.d0-orp) + psi_diff_new*orp
			else
				psi_diff(i,j) = psi_diff_new
			endif

			diff_error(i,j) = psi_diff_new - psi_diff_old

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

	enddo
	enddo

	mx = 0.d0

	! if we use the GS equation for psi_diff, we might as well update big_Psi immediately, but let's keep things separated for clarity
	do j = 1, nz
	do i = 1, nx

			if(sort_grid(i,j)<=0) cycle

			big_Psi(i,j) = psi(i,j) + psi_diff(i,j)
			mx = dmax1(big_Psi(i,j),mx)

	enddo
	enddo
	!----------------end of big_Psi cycle-----------------

	! Set the new maximum value of big_Psic
	big_Psic = mx
	! "Finished: big_Psi Update"

	! CALL BC HERE AGAIN FOR BIG_PSI
	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,2)

	! -----------------------------------------------------

	psi_max = max(psic, big_Psic)
	psic_13 = psi_max

     ! update RHS
     	if(nx>=inter_switch) then
			call GS_RHS_all_Gauss_TF
		endif


	! Move the density calculation to a separate function
	! Update density and seek mach theta max

!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)
!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)

!!	if(Broot==0) then
!!		call update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	else
!!		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	endif

	if((nx>=inter_switch).and.(Broot==0)) then
!!$		if(k<50) then
!!$			call update_rho_TF_gauss(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!$								min_ix,min_iz)
		! changed March 1 2020
		! added attempt to use both on March 3 2020
!!$		else
			call update_rho_TF_gauss_2020_magnetic(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
								min_ix,min_iz)
!!$		endif
	else
		call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)
	! here update just n_den

	if(in_orp <= 0.0d0) then
		! Chebyshev acceleration 
		if(nx >= 5) then
			std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
		end if
	endif

	if ((accelerate).and.(nx<inter_switch)) then
		orp = std_orp
	else
		orp = fix_orp
	endif

	!!$ call step_output(nx,psi(:,:),rho,residual)

	! -----------------------------------------------------

	if(k == 1) then
		eps_anorm = eps*anorm
		eps_anorm_2 = eps*anorm_2
	else
		if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
			! As we get closer to convergence, we want the solution
			! to relax.  So as anorm approaches eps_anorm we want 
			! the over relaxation parameter to go to some const. ~ 1
			! Use x and mtm_soln as a temporary variable
			x = anorm/eps_anorm
			mtm_soln = 1.0d0 ! The limiting value for x ~ 1
			tmp = x/(x + orp/mtm_soln - 1.0d0)
			tmp = dmin1(tmp,1.0d0)
			orp = orp*tmp
			orp = dmax1(orp,mtm_soln)
		endif
	end if

	if(k == 1 .or. modulo(k,k_save) == 0) then
		print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		print *, k,": anorm_2 = ",real(anorm_2,skind)," eps_anorm_2 = ",real(eps_anorm_2,skind)
		print *, "The Over-Relaxation Parameter =",orp,std_orp
		call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)
		continue
	end if

	write(111,*) k,anorm
	write(112,*) k,anorm_2

	! Check the convergence criteria
	if((anorm < eps_anorm) .and. (anorm_2 < eps_anorm_2) .and. (k > min_it)) exit

	end do

	print *, k,": anorm = ", anorm, " eps_anorm = ",eps_anorm
	print *, k,": anorm_2 = ", anorm_2, " eps_anorm_2 = ", eps_anorm_2
	print *, "Average residual =", anorm/(nx*nz)
	print *, "Average difference (2) =", anorm_2/(nx*nz)

	print *, "Mach Theta Max =", mach_theta_max
	print *, "Final Solution has psi Center = ", psic
	print *, "and big_Psi Center = ", big_Psic
	print *

	return

	contains


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS_all_Gauss_TF
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind) :: tdx, tdz


    tdx = dx*2.d0
    tdz = dz*2.d0

    gradpsi_RHS = 1.d-10
    grad_bigPsi_RHS = 1.d-10

    ! fill grad psi
    do j = 1, nz
       do i = 1, nx

			!-------------------psi-----------------------
          if(sort_grid(i,j)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi(i+1,j)-psi(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi(i,j+1)-psi(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)

			!-------------------big_Psi-----------------------
          grad_bigPsi_RHS(1,i,j) = (big_Psi(i+1,j)-big_Psi(i-1,j))/tdx
          grad_bigPsi_RHS(2,i,j) = (big_Psi(i,j+1)-big_Psi(i,j-1))/tdz
          grad_bigPsi_RHS(0,i,j) = sqrt(grad_bigPsi_RHS(1,i,j)**2+grad_bigPsi_RHS(2,i,j)**2)

          if(grad_bigPsi_RHS(0,i,j)<1.d-6) grad_bigPsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)


       enddo
    enddo


!	if((k==1).or.(modulo(k,k_save)==0)) then
!
!    ! write results for debugging
!
!		open(13, file='gradpsi_RHS.plt')
!
!		write(13,*)'TITLE="solution of GS equation with flow"'
!		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
!		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'
!
!		do j=1,nz
!
!		   z = z_coord(j)
!
!		   do i=1,nx
!
!			  x = x_coord(i)
!
!			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
!
!		   end do
!		end do
!
!		close(13)
!
!	endif

90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS_all_Gauss_TF

end subroutine ngs_solve_TF_gauss_magnetic


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  ! This version differentiates between "vacuum" and plasma in one of two ways:
  ! 1) Use psi/big_Psi>0 to divide the zones (and keep the 0.4*z_size criterion in)
  ! 2) Update the interface between vacuum and plasma every few iterations
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_TF_free_boundary(psi,big_Psi,psi_diff,n_den,residual,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! cycle separately on psi and psi_diff
! IMPORTANT: solve directly for psi_diff, then calculate big_Psi as psi + psi_diff

	use constant, only : dx, dz, dx_a, dz_a

	implicit none

	integer, intent(in) :: nx,nz,min_it,max_it,flag
	real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,big_Psi,psi_diff,n_den,residual
	real (kind=dkind), intent(in) :: eps
	! Input Over Relaxation Parameter
	! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
	! else we set orp = in_orp
	real (kind=dkind), intent(in) :: in_orp
	integer ::ipass,i,isw,j,jsw,k,alloc_stat
	real (kind=dkind) :: res, den, x !, xl, xr
	! res -> residual
	! den -> denominator of the Newton-Gauss-Seidel update
	! x -> Radial position of grid point (i,j)
	real (kind=dkind) :: dx2,dz2,mx
	real (kind=dkind) :: anorm, eps_anorm, anorm_2, eps_anorm_2
	real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
	real (kind=dkind) :: orp, std_orp ! over relaxation parameter
	real (kind=dkind) :: n_denc, phic
	! by is the phi component of the magnetic field
	real (kind=dkind) :: Fstar
	real (kind=dkind) :: last_mtm ! last mach theta max
	! minimum delta rho for bernoulli roots, the mean val. and location
	real (kind=dkind) :: min_dn_den, tmp
	integer :: min_ix, min_iz
	! The next set of variables are used for the bisection search for
	! locating the degenerate root for Mach Theta Max
	real (kind=dkind) :: mtm_soln
	! mtm_acc controls the bisection loop seeking mach theta max
	real (kind=dkind), parameter :: mtm_acc = 1.0d-12
	! VERY IMPORTANT
	! psi_degen records the value of psi where the solution to the
	! Bernoulli equation is approximately degenerate.  It is 
	! approximately equal to 0.5d0*psic but generally slightly less.
	real (kind=dkind) :: gPsi2
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term0_5
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: big_Psi0, big_Psil, big_Psir, big_Psid, big_Psiu
	real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real (kind=dkind) :: phinrx, phinlx, phinrz, phinlz, phinc
	real (kind=dkind) :: n_denrx,n_denlx,n_denrz,n_denlz
	real (kind=dkind) :: phirx,philx,phirz,philz
	real(kind=dkind), dimension(-1:1,-1:1) :: big_Psi_around
	real(kind=dkind) :: psi_diff_old, psi_diff_new, psi_diff_limit
	real(kind=dkind) :: H_loc, H_diff_loc
	real(kind=dkind) :: dpsidx, dpsidz, drho_ds
!	integer :: psi_diff_option = 2
	real (kind=dkind), dimension(1:nx,1:nz) :: diff_error
	integer :: k_save = 25
	integer :: i_issue_check = 101
	integer :: j_issue_check = 92 ! For debugging

	psi_max = max(psic,big_Psic)

    if(psi_degen==0.d0) then 
		psi_degen = 0.5d0*psic ! Initialization
		big_Psi_degen = 0.5d0*psic ! Initialization
	endif
	last_mtm = mach_theta_max
	eps_anorm = 0.0d0
	eps_anorm_2 = 0.0d0

	fix_orp0 = fix_orp

	! The under/over relaxation parameter
	! Standard Case:
	if(in_orp <= 0.0d0) then
		! Note that the parameters below are tuned for starting with 
		! a 17 x 17 grid.
		orp = 1.0d0
		if(nx <= 5) then
			orp = 0.5d0
			rjac = 0.75d0
		else if(nx <= 9) then
			orp = 0.5d0
			rjac = 0.9d0
		else if(nx <= 17) then
			rjac = 0.98d0
		else if(nx <= 33) then
			rjac = 0.995d0
		else if(nx <= 65) then
			rjac = 0.9972d0
		else
		! It would appear that 129 is about where the solution begins to converge
			rjac = 0.0d0
		end if
		std_orp = orp
	else
		print *, "Input Over Relaxation Parameter =",in_orp
		orp = in_orp
		std_orp = in_orp
	endif

	if (accelerate) then
		continue
	else
		orp = fix_orp
	endif

	dx2 = dx*dx
	dz2 = dz*dz

	diff_error = 0.d0

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,0)

	! Update rho before "relaxing psi" but do not seek mach theta max
!!	if(Broot==0) then
!!		call update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	else
		call update_rho_TF_free_boundary(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
!!	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)

	call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)

	! Iterate until convergence
	do k=1, max_it
	!     "Update psi"6

	!		call cpu_time(tic)

		if((k>ite_Broot).and.(Broot/=Broot_true)) then
			Broot = Broot_true
		endif

	! Reset the norm and max of psi and big_Psi
	mx = 0.0d0
	anorm = 0.0d0
	anorm_2 = 0.0d0

!!$	March 12 2021: this should NOT be necessary here
!!$	do j=2,nz-1
!!$	do i=2,nx-1
!!$		if(sort_grid(i,j)>0) psi(i,j) = dabs(psi(i,j))	!	segno
!!$	enddo
!!$	enddo

	!----------------------psi cycle------------------------
	jsw = 1
	do ipass = 1, 2
		isw=jsw;
		do j=2, nz-1
		do i=isw+1, nx-1, 2
		! "Update psi(",i,",",j,")"
		! Only solve the inner region problem

			if(sort_grid(i,j)<=0) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			big_Psil = 0.d0
			big_Psir = 0.d0
			big_Psiu = 0.d0
			big_Psid = 0.d0

			big_Psi0 = big_Psi(i,j)
			if(i>1) big_Psil = big_Psi(i-1,j)
			if(i<nx) big_Psir = big_Psi(i+1,j)
			if(j>1) big_Psid = big_Psi(i,j-1)
			if(j<nz) big_Psiu = big_Psi(i,j+1)

			x = x_coord(i)

			! Calculate rho and phi at the current location
			n_denc = n_den(i,j)
			phic = phi_TF_ofpsi(big_Psi0)

			! old formulation with new Bernoulli

			drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

			! Right x
			phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
			n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
			ma2rx = mu_mag*mass * phirx*phirx/n_denrx

			! Left x
			philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
			n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
			ma2lx = mu_mag*mass * philx*philx/n_denlx

			! -----------------------------------------------------

			drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

			! Right z
			phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
			n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
			ma2rz = mu_mag*mass * phirz*phirz/n_denrz

			! Left z
			philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
			n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
			ma2lz = mu_mag*mass * philz*philz/n_denlz

			! -----------------------------------------------------

			big_Psi_around(-1,0) = big_Psil
			big_Psi_around(1,0) = big_Psir
			big_Psi_around(0,0) = big_Psi0
			big_Psi_around(0,-1) = big_Psid
			big_Psi_around(0,1) = big_Psiu

			Fstar = Fstarofpsi(psi0,big_Psi0)

			call psi_derivative_new(i,j,nx,nz,big_Psi_around,dpsidx,dpsidz)
			gPsi2 =  dpsidx**2+dpsidz**2

			! we follow FLOW nomenclature for term names (everything should be easily recognizable)

			term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
							+ (psil-psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

			term0_5 = - (  &
							( ma2rx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ ma2lx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( ma2rz*(big_Psiu-big_Psi0) &
							+ ma2lz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

			term1 = eV*mu_mag * Fstar *(phic-phi_e_TF_ofpsi(psi0))/x**2
			term3 = mu_mag*mass * phic * dphi_TF_dpsi(big_Psi0) * gPsi2 / n_denc/x**2

			term2 = mu_mag * eV * n_denc * (omegaofpsi(big_Psi0)-omegaofpsi(psi0)) -  &
					mu_mag * n_denc * (dp_i_dpsi(big_Psi0)/d_TF_ofpsi(big_Psi0)-dp_i_dpsi(psi0)/d_TF_ofpsi(psi0))

			term4 = mu_mag * n_denc * (dh_e_dpsi_partial(psi0) + dh_i_dpsi_partial(big_psi0))

			term5 = mu_mag*mass * (n_denc**gamma_i/(gamma_i-1.0d0) * ds_i_dpsi(big_Psi0) +  &
												n_denc**gamma_e/(gamma_e-1.0d0) * ds_e_dpsi(psi0))
			! CHECK DIMENSIONS FOR ALL!!!

			!				if((nx>=inter_switch).and.(term6/=0.d0)) then
			!					print*, i,j, term0, term1, term2, term3, term4, term5, term6
			!					print*, '   '
			!					pause
			!				endif

			res = term0 + term0_5 + term1 + term2 + term3 + term4 - term5

			den = -( 1.0d0/(x + 0.5d0*dx) + &
						1.0d0/(x - 0.5d0*dx) )/(x*dx2) &
						-2.0d0/(x*x*dz2)

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

			continue

			!res = res*mu_mag

			! Store the residual
			if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
				if(nx>=inter_switch) then
					residual(i,j) = res/den
				else
					residual(i,j) = res
				endif
			end if

			if(res>1.d0) then
				continue
			endif

			! -----------------------------------------------------

			continue

			psi(i,j) = psi0 - orp*res/den

			! Find the max absolute value of psi in the plasma
			mx = dmax1(dabs(psi(i,j)),mx)

			anorm = anorm + dabs(res)

		end do
		isw = 3 - isw
		end do
		jsw = 3 - jsw
		call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,1)
		! only update psi on BC
	end do
	!----------------end of psi cycle--------------

	! Set the new maximum value of psic
	psic = mx
	psi_max = max(psic,big_Psic)
	psic_13 = psi_max
	! "Finished: psi Update"

	!----------------big_Psi cycle-----------------
	! we use psi_diff as dependent variable, and solve the algebraic equation for it
	! FIRST ATTEMPT: USE OLD ITERATION FOR DERIVATIVES

	! since psi_diff appears in the same way in two equations, we can get it from the GS one, which
	! is simpler, as it does not have derivatives of big_Psi and n_den
	do j = 1, nz
	do i = 1, nx

			if(sort_grid(i,j)<=0) cycle

			! set up local psi values

			psil = 0.d0
			psir = 0.d0
			psiu = 0.d0
			psid = 0.d0

			psi0 = psi(i,j)
			if(i>1) psil = psi(i-1,j)
			if(i<nx) psir = psi(i+1,j)
			if(j>1) psid = psi(i,j-1)
			if(j<nz) psiu = psi(i,j+1)

			x = x_coord(i)

			psi_diff_old = psi_diff(i,j)

			n_denc = n_den(i,j)

			Fstar = Fstarofpsi(psi0,big_Psi(i,j))

			if(psi_diff_option==1) then
			! use magnetic GS for psi_diff

				term0 = ( (psir-psi0)/(x + 0.5d0*dx) &
								+ (psil-psi0)/(x - 0.5d0*dx) )  &
								/(x*dx2) &
								+ ( (psiu-psi0) + (psid-psi0) )/(x*x*dz2) 

				term2 = -eV*mu_mag * Fstar * phi_e_TF_ofpsi(psi0) / x**2
				term3 = mu_mag * n_denc * dh_e_dpsi(psi0)
				term4 = mu_mag * mass * n_denc**gamma_e * ds_e_dpsi(psi0)/(gamma_e-1.0d0)

				psi_diff_new = (-term0 - term2 - term3 + term4) / (ev**2*mu_mag*n_denc/mass/x**2)

			elseif(psi_diff_option==2) then
			! use fluid GS for psi_diff

				big_Psil = 0.d0
				big_Psir = 0.d0
				big_Psiu = 0.d0
				big_Psid = 0.d0

				big_Psi0 = big_Psi(i,j)
				if(i>1) big_Psil = big_Psi(i-1,j)
				if(i<nx) big_Psir = big_Psi(i+1,j)
				if(j>1) big_Psid = big_Psi(i,j-1)
				if(j<nz) big_Psiu = big_Psi(i,j+1)

				drho_ds = van_Leer_slope_new(n_den(i-1,j),n_den(i,j),n_den(i+1,j),dx_a(i-1),dx_a(i))

				! Right x
				phirx = phi_TF_ofpsi(0.5d0*(big_Psir + big_Psi0))
				n_denrx = n_den(i,j) + 0.5d0*dx*drho_ds
				phinrx = mass * phirx/n_denrx

				! Left x
				philx = phi_TF_ofpsi(0.5d0*(big_Psil + big_Psi0))
				n_denlx = n_den(i,j) - 0.5d0*dx*drho_ds
				phinlx = mass * philx/n_denlx

				! -----------------------------------------------------

				drho_ds = van_Leer_slope_new(n_den(i,j-1),n_den(i,j),n_den(i,j+1),dz_a(j-1),dz_a(j))

				! Right z
				phirz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psiu))
				n_denrz = n_den(i,j) + 0.5d0*dz*drho_ds
				phinrz = mass * phirz/n_denrz

				! Left z
				philz = phi_TF_ofpsi(0.5d0*(big_Psi0 + big_Psid))
				n_denlz = n_den(i,j) - 0.5d0*dz*drho_ds
				phinlz = mass * philz/n_denlz

				term0_5 = phi_TF_ofpsi(big_Psi0) * (  &
							( phinrx*(big_Psir-big_Psi0)/(x + 0.5d0*dx) &
							+ phinlx*(big_Psil-big_Psi0)/(x - 0.5d0*dx) )  &
							/(x*dx2) &
							+ ( phinrz*(big_Psiu-big_Psi0) &
							+ phinlz*(big_Psid-big_Psi0) )/(x*x*dz2)  &
							)

				term2 = eV* Fstar * phi_TF_ofpsi(big_Psi0) / x**2
				term3 = n_denc * dh_i_dpsi(big_Psi0)
				term4 = mass * n_denc**gamma * ds_i_dpsi(big_Psi0)/(gamma-1.0d0)

				psi_diff_new = (-term0_5 + term2 + term3 - term4) / (ev**2*n_denc/mass/x**2)
				
			endif

			if((i==i_issue_check).and.(j==j_issue_check)) then
				continue
			endif

			call get_H_diff(psi(i,j),big_Psi(i,j),H_diff_loc) ! dbsval, or integrate hofpsi(psi(i,j))
			H_loc = h_e_ofpsi_partial(psi(i,j)) + h_i_ofpsi_partial(big_Psi(i,j)) + eV*H_diff_loc	! the R^2 Omega^2 term is hidden in H (i.e., H_diff)

			if(H_loc<=0.d0) then
				print*, 'problem in psi_diff max'
				continue
				psi_diff_limit = 0.d0
			else
				psi_diff_limit = sqrt(0.9d0*2.d0*mass*H_loc) * x/eV
			endif

			anorm_2 = anorm_2 + dabs(psi_diff_old - psi_diff_new)

			if(abs(psi_diff_new)>psi_diff_limit) then
				psi_diff_new = sign(psi_diff_limit,psi_diff_new)
			endif

			if(orp<=1.d0) then
				psi_diff(i,j) = psi_diff_old*(1.d0-orp) + psi_diff_new*orp
			else
				psi_diff(i,j) = psi_diff_new
			endif

			diff_error(i,j) = psi_diff_new - psi_diff_old

	enddo
	enddo

	mx = 0.d0

	! if we use the GS equation for psi_diff, we might as well update big_Psi immediately, but let's keep things separated for clarity
	do j = 1, nz
	do i = 1, nx

			if(sort_grid(i,j)<=0) cycle

			big_Psi(i,j) = psi(i,j) + psi_diff(i,j)
			mx = dmax1(big_Psi(i,j),mx)

	enddo
	enddo
	!----------------end of big_Psi cycle-----------------

	! Set the new maximum value of big_Psic
	big_Psic = mx
	! "Finished: big_Psi Update"

	! CALL BC HERE AGAIN FOR BIG_PSI
	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,2)

	! -----------------------------------------------------

	psi_max = max(psic, big_Psic)
	psic_13 = psi_max

	! Move the density calculation to a separate function
	! Update density and seek mach theta max

!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,0,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)
!	call update_rho_TF(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!						min_ix,min_iz)

!!	if(Broot==0) then
!!		call update_rho_TF_relax(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
!!							min_ix,min_iz)
!!	else
		call update_rho_TF_free_boundary(psi,big_Psi,psi_diff,n_den,nx,nz,1,mtm_acc,min_dn_den,  &
							min_ix,min_iz)
!!	endif

	call bc_TF_0(psi,big_Psi,psi_diff,n_den,nx,nz,3)
	! here update just n_den

	if(in_orp <= 0.0d0) then
		! Chebyshev acceleration 
		if(nx >= 5) then
			std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
		end if
	endif

	if ((accelerate).and.(nx<inter_switch)) then
		orp = std_orp
	else
		orp = fix_orp
	endif

	!!$ call step_output(nx,psi(:,:),rho,residual)

	! -----------------------------------------------------

	if(k == 1) then
		eps_anorm = eps*anorm
		eps_anorm_2 = eps*anorm_2
	else
		if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
			! As we get closer to convergence, we want the solution
			! to relax.  So as anorm approaches eps_anorm we want 
			! the over relaxation parameter to go to some const. ~ 1
			! Use x and mtm_soln as a temporary variable
			x = anorm/eps_anorm
			mtm_soln = 1.0d0 ! The limiting value for x ~ 1
			tmp = x/(x + orp/mtm_soln - 1.0d0)
			tmp = dmin1(tmp,1.0d0)
			orp = orp*tmp
			orp = dmax1(orp,mtm_soln)
		endif
	end if

	if(k == 1 .or. modulo(k,k_save) == 0) then
		print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		print *, k,": anorm_2 = ",real(anorm_2,skind)," eps_anorm_2 = ",real(eps_anorm_2,skind)
		print *, "The Over-Relaxation Parameter =",orp,std_orp
		call step_output_TF(nx, psi, n_den, residual, big_Psi, psi_diff, diff_error)
		continue
	end if

	write(111,*) k,anorm
	write(112,*) k,anorm_2

	! Check the convergence criteria
	if((anorm < eps_anorm) .and. (anorm_2 < eps_anorm_2) .and. (k > min_it)) exit

	end do

	print *, k,": anorm = ", anorm, " eps_anorm = ",eps_anorm
	print *, k,": anorm_2 = ", anorm_2, " eps_anorm_2 = ", eps_anorm_2
	print *, "Average residual =", anorm/(nx*nz)
	print *, "Average difference (2) =", anorm_2/(nx*nz)

	print *, "Mach Theta Max =", mach_theta_max
	print *, "Final Solution has psi Center = ", psic
	print *, "and big_Psi Center = ", big_Psic
	print *

	return

end subroutine ngs_solve_TF_free_boundary


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine step_output(n,psi,rho,residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: n
	real (kind=skind), dimension(n,n) :: out
	real (kind=dkind), dimension(n,n) :: rho, psi,residual
	integer :: i,j
	real(kind=dkind) :: X0loc, theta, mth, mphi,ex,ez,dx,dz, x

	out(1:n,1:n) = psi(1:n,1:n)
	call radial_plot(out,psi,n,n,"psi_step",n/2)

	out(1:n,1:n) = rho(1:n,1:n)
	call radial_plot(out,psi,n,n,"rho_step",n/2)

	out(1:n,1:n) = residual(1:n,1:n)
	call radial_plot(out,psi,n,n,"residual_step",n/2)


	continue

	return

  end subroutine step_output
  ! ------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine step_output_TF(n, psi, n_den, residual, big_Psi, psi_diff, diff_error)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: n
	real (kind=skind), dimension(n,n) :: out
	real (kind=dkind), dimension(n,n) :: n_den, psi, residual, big_Psi, psi_diff, diff_error
	integer :: i,j
	real(kind=dkind) :: X0loc, theta, mth, mphi,ex,ez,dx,dz, x

	out(1:n,1:n) = psi(1:n,1:n)
	call radial_plot(out,psi,n,n,"psi_step",n/2)

	out(1:n,1:n) = n_den(1:n,1:n)
	call radial_plot(out,psi,n,n,"n_den_step",n/2)

	out(1:n,1:n) = residual(1:n,1:n)
	call radial_plot(out,psi,n,n,"residual_step",n/2)

	out(1:n,1:n) = big_Psi(1:n,1:n)
	call radial_plot(out,psi,n,n,"big_Psi_step",n/2)

	out(1:n,1:n) = psi_diff(1:n,1:n)
	call radial_plot(out,psi,n,n,"psi_diff_step",n/2)

	out(1:n,1:n) = diff_error(1:n,1:n)
	call radial_plot(out,psi,n,n,"diff_error_step",n/2)

	return

  end subroutine step_output_TF
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------
  ! These routine calculates the various physical
  ! quantities such as the density, pressure, velocity, etc.
  ! Recall that (x, y, z) -> (r, phi, z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine physical_var(psi,rho,nx,nz,vx,vy,vz,bx,by,bz,p,  &
							malfven2,mslow,mcusp,mpol,beta, &
							j_phi,temp,mtor,j_par, j_x, j_z,cs,csp,hyper, &
							 j_phi_new)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a
	use pseudo_IMSL, only : DBSNAK, DBS2VL, DBS2IN

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, rho
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: &
         vx,vy,vz,bx,by,bz,p,malfven2,mslow,mcusp,mpol,beta,  &
		 temp,mtor,j_phi,j_par, j_x, j_z,cs,csp,  &
		 hyper, j_phi_new
    ! malfven_2 is the square of the alfvenic mach number
!	real (kind=dkind), dimension(1:nx,1:nz) :: j_x, j_z
    real (kind=dkind) :: phic ! Temp Var. for phiofpsi(psi(j,k))
    real (kind=dkind) :: omegac ! Temp Var. for omegaofpsi(psi(j,k))
    real (kind=dkind) :: vp2 ! Square of the poloidal velocity
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: dx2,dz2
    integer :: i,j,k
    real (kind=dkind) :: x,rloc,bfield,delstar
    real (kind=dkind) :: cf2,cs2,a2,ca2,ct2,cusp2
	real(kind=dkind), dimension(-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: term1,term2,term3,term4,term5,b_dot_v, term6
	real(kind=dkind) :: Alf, plocg, bpol
	real (kind=dkind), dimension(2) :: Gloc ! (1-M_A^2) grad psi /R^2
    real (kind=dkind) :: phirx,philx,phirz,philz
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: drho_ds
	integer :: i_zone

    dx2 = dx*dx
    dz2 = dz*dz

    ! Calculate the poloidal magnetic field
    do k = 1, nz
       do j = 1, nx

		 if(sort_grid(j,k)>=1) then

          ! Calculate the position of the (j,k) grid point
          x = x_coord(j)

			psiprimx(0) = ( dx_a(j-1)**2*psi(j+1,k) +  &
							(dx_a(j)**2-dx_a(j-1)**2)*psi(j,k) -  &
							dx_a(j)**2*psi(j-1,k) ) /  &
							( dx_a(j)*dx_a(j-1)*(dx_a(j)+dx_a(j-1)) )

			psiprimx(1) = (psi(j+1,k) - psi(j,k))/dx_a(j)

			psiprimx(-1) = (psi(j,k) - psi(j-1,k))/dx_a(j-1)

			psiprimz(0) = ( dz_a(k-1)**2*psi(j,k+1) +  &
							(dz_a(k)**2-dz_a(k-1)**2)*psi(j,k) -  &
							dz_a(k)**2*psi(j,k-1) ) /  &
							( dz_a(k)*dz_a(k-1)*(dz_a(k)+dz_a(k-1)) )

			psiprimz(1) = (psi(j,k+1) - psi(j,k))/dz_a(k)


			psiprimz(-1) = (psi(j,k) - psi(j,k-1))/dz_a(k-1)

          ! Calculate B_r
          if(k == 1) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_r as independent of z near the boundary.
             bx(j,k) = (psi(j,k+1) - psi(j,k))/dz_a(1)/x
          else if(k == nz) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_r as independent of z near the boundary.
             bx(j,k) = (psi(j,k) - psi(j,k-1))/dz_a(nz)/x
          else
             ! use a centered difference
			 bx(j,k) = ( dz_a(k-1)**2*psi(j,k+1) +  &
						 (dz_a(k)**2-dz_a(k-1)**2)*psi(j,k) -  &
						 dz_a(k)**2*psi(j,k-1) ) /  &
						 (dz_a(k)*dz_a(k-1)*(dz_a(k)+dz_a(k-1)))/x
          end if

		  bx(j,k) = - bx(j,k) ! WARNING: the derivatives are taken with the PLUS sign for convenience!

          ! Calculate Bz
          if(j == 1) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_z as independent of x near the boundary.
             bz(j,k) = (psi(j+1,k) - psi(j,k))/dx_a(1)/(x + 0.5d0*dx_a(1))
          else if(j == nx) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_z as independent of x near the boundary.
             bz(j,k) = (psi(j,k) - psi(j-1,k))/dx_a(nx)/(x - 0.5d0*dx_a(nx))
          else
             ! use a centered difference
			 bz(j,k) = ( dx_a(j-1)**2*psi(j+1,k) +  &
						 (dx_a(j)**2-dx_a(j-1)**2)*psi(j,k) -  &
						 dx_a(j)**2*psi(j-1,k) ) /  &
						 (dx_a(j)*dx_a(j-1)*(dx_a(j)+dx_a(j-1)))/x
          end if

		  ! -------------------------------------differentiate here for external points-------------------------------------
		  ! October 8 2021 Added further differentiation for bc_type==7 and external region with plasma, calculated
		  ! with the "inner region".

		if( (((tri_type/=13).AND.(bc_type/=7)).AND.((psi(j,k)/psic)<fraction))  &
						.OR.  &
!			((bc_type==7).AND.(sqrt((x_coord(j)-rmajor)**2+z_coord(k)**2)>0.4d0*z_size))  &
!						.OR.  &
			((tri_type==13).AND.(sort_grid(j,k)==1)) ) then

			! external point

			phic = 0.d0

			vx(j,k) = 0.d0
			vz(j,k) = 0.d0
			by(j,k) = dsqrt(mu_mag)*iofpsi(0.d0)/x_coord(j)

			vy(j,k) = 0.d0

			! Finally we calculate the pressure
				p(j,k) = sofpsi(0.d0)*rho(j,k)**gamma
			! Calculate the plasma beta
				beta(j,k) = 2.0d0*mu_mag*p(j,k)/(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)

			! And the plasma temperature
			temp(j,k) = mass*p(j,k)/rho(j,k)/eV

			! Now the toroidal current density

			if ((j>1).AND.(k>1).AND.(j<nx).AND.(k<nz)) then

				delstar = 0.D0

				j_phi(j,k) = x*delstar ! "+" is the correct sign
				j_phi_new(j,k) = 0.d0

			endif

			! Now let's calculate some Mach numbers
			! First some intermediate quantities

			vp2 = 0.0d0

			! Sound Speed
			a2 = gamma*p(j,k)/rho(j,k)
			! Poloidal Alfven Speed Squared
			if(phic .eq. 0.0d0) then
			ca2 = (bx(j,k)*bx(j,k) + bz(j,k)*bz(j,k))/mu_mag/rho(j,k)
			endif
			! Toroidal Alfven Speed Squared
			ct2 = by(j,k)*by(j,k)/mu_mag/rho(j,k)
			! Fast Mode Speed Squared
			cf2 = 0.5d0*(a2 + ca2 + ct2 &
			+ dsqrt((a2 + ct2 - ca2)**2 + 4.0d0*ca2*ct2) )
			! Slow Mode Speed Squared
			cs2 = a2*ca2/cf2

			! Calculate the square of the Alfvenic Mach Number
			malfven2(j,k) = phic**2/rho(j,k)
			!          malfven2(j,k) = mu_mag*phic**2/rho(j,k)

			! Now the slow mode Mach number
			mslow(j,k) = dsqrt(malfven2(j,k)*cf2/a2)

			mslow(j,k) = sqrt(vp2/a2)*sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

			! Now the cusp Mach number
			mcusp(j,k) = mslow(j,k)*dsqrt(1.0d0 + cs2/cf2)

			mpol(j,k) = sqrt((vx(j,k)**2+vz(j,k)**2)/a2) *sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

			! Now the toroidal Mach number
			mtor(j,k) = dsqrt(vy(j,k)**2/a2)

			! sound speeds
			cs(j,k) = sqrt(a2)
			csp(j,k) = sqrt(a2*(bx(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+by(j,k)**2+bz(j,k)**2))

			hyper(j,k) = 1.d0

			cycle

		endif

		  ! -------------------------------------proceed from here for inner points-------------------------------------
		  ! October 8 2021 Added further differentiation for bc_type==7 and external region with plasma

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(j,k)-2)
				i_zone = max(-1,i_zone)
			endif

          ! Initialize phic
          phic = phiofpsi(psi(j,k),i_zone)
          omegac = omegaofpsi(psi(j,k),i_zone)

          ! Now let's calculate the poloidal velocity components
          vx(j,k) = phic*bx(j,k)/rho(j,k)/dsqrt(mu_mag)
          vz(j,k) = phic*bz(j,k)/rho(j,k)/dsqrt(mu_mag)

          ! Next we calculate B_phi
			by(j,k) = dsqrt(mu_mag)*(iofpsi(psi(j,k),i_zone)/x + x*phic*omegac)/ &
				   (1.0d0 - phic**2/rho(j,k))

          ! Next we calculate V_phi
          vy(j,k) = (phic*iofpsi(psi(j,k),i_zone)/rho(j,k)/x + x*omegac)/ &
               (1.0d0 - phic**2/rho(j,k))

          ! Finally we calculate the pressure
			p(j,k) = sofpsi(psi(j,k),i_zone)*rho(j,k)**gamma
          ! Calculate the plasma beta
		  beta(j,k) = 2.0d0*mu_mag*p(j,k)/(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)

		  ! And the plasma temperature
			temp(j,k) = mass*p(j,k)/rho(j,k)/eV

          ! Now the toroidal current density

		 if ((j>1).AND.(k>1).AND.(j<nx).AND.(k<nz)) then

				b_dot_v = x*omegaofpsi(psi(j,k),i_zone)*by(j,k) + (phiofpsi(psi(j,k),i_zone)/rho(j,k))*  &
						( bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2 )/sqrt(mu_mag)

				term1= b_dot_v*dphidpsi(psi(j,k),i_zone)/dsqrt(mu_mag)
				term2= x*(phiofpsi(psi(j,k),i_zone)*by(j,k)/dsqrt(mu_mag) +  &
						rho(j,k)*x*omegaofpsi(psi(j,k),i_zone))*domegadpsi(psi(j,k),i_zone)
				term3= by(j,k)*didpsi(psi(j,k),i_zone)/x/dsqrt(mu_mag)
				term4= rho(j,k)*dhdpsi(psi(j,k),i_zone)
				term5= rho(j,k)**gamma*dsdpsi(psi(j,k),i_zone)/(gamma-1.0d0)

				delstar = term1+term2+term3+term4-term5

				term6 = (psi(j+1,k)-psi(j-1,k)) *  &
								(phiofpsi(psi(j+1,k),i_zone)**2/rho(j+1,k)-phiofpsi(psi(j-1,k),i_zone)**2/rho(j-1,k)) / (4.d0*dx2) +  &
							(psi(j,k+1)-psi(j,k-1)) *  &
								(phiofpsi(psi(j,k+1),i_zone)**2/rho(j,k+1)-phiofpsi(psi(j,k-1),i_zone)**2/rho(j,k-1)) / (4.d0*dz2)
				term6 = term6/x**2


			j_phi(j,k) = x * ( delstar - term6) / (1.d0-phic**2/rho(j,k))
			! "+" is the correct sign

			j_phi_new(j,k) = ( -( (psi(j+1,k)-psi(j,k))/(x+dx/2.d0) - (psi(j,k)-psi(j-1,k))/(x-dx/2.d0) )/dx2  &
										- (psi(j,k+1)-2.d0*psi(j,k)+psi(j,k-1))/dz2/x ) / mu_mag

		 endif

		  ! Now let's calculate some Mach numbers
          ! First some intermediate quantities

          ! Rather than using the poloidal velocities just calculated
          ! which involve differentiation, let's use the Bernoulli
          ! equation to determine the square of the poloidal velocity
          ! which is really all we need for the Mach numbers
          if(phic .eq. 0.0d0) then
             vp2 = 0.0d0
          else
			vp2 = 2.0d0*hofpsi(psi(j,k),i_zone) + (x*omegac)**2 &
				  - (phic*by(j,k)/rho(j,k))**2 &
				  - 2.0d0*gamma/(gamma-1.0d0)*p(j,k)/rho(j,k)
             ! Because of roundoff error this can be negative.
             vp2 = dmax1(0.0d0,vp2)
          endif

          ! Sound Speed
		  a2 = gamma*p(j,k)/rho(j,k)
          ! Poloidal Alfven Speed Squared
          if(phic .eq. 0.0d0) then
             ca2 = (bx(j,k)*bx(j,k) + bz(j,k)*bz(j,k))/mu_mag/rho(j,k)
          else
             ca2 = rho(j,k)*vp2/(mu_mag*phic**2)
          endif
          ! Toroidal Alfven Speed Squared
          ct2 = by(j,k)*by(j,k)/mu_mag/rho(j,k)
          ! Fast Mode Speed Squared
          cf2 = 0.5d0*(a2 + ca2 + ct2 &
               + dsqrt((a2 + ct2 - ca2)**2 + 4.0d0*ca2*ct2) )
          ! Slow Mode Speed Squared
          cs2 = a2*ca2/cf2

          ! Calculate the square of the Alfvenic Mach Number
          malfven2(j,k) = phic**2/rho(j,k)
!          malfven2(j,k) = mu_mag*phic**2/rho(j,k)

          ! Now the slow mode Mach number
          mslow(j,k) = dsqrt(malfven2(j,k)*cf2/a2)

		  mslow(j,k) = sqrt(vp2/a2)*sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

          ! Now the cusp Mach number
          mcusp(j,k) = mslow(j,k)*dsqrt(1.0d0 + cs2/cf2)

			mpol(j,k) = sqrt((vx(j,k)**2+vz(j,k)**2)/a2) *sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

		  ! Now the toroidal Mach number
          mtor(j,k) = dsqrt(vy(j,k)**2/a2)

		  ! sound speeds
		  cs(j,k) = sqrt(a2)
		  csp(j,k) = sqrt(a2*(bx(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+by(j,k)**2+bz(j,k)**2))

!!$          ! Note: cs2 is only zero when the poloidal field goes to zero.
!!$          ! Since rho*vp = Phi * Bp and Phi = Phi(psi) is finite,
!!$          ! this implies that vp is also zero. The correct limiting
!!$          ! value should be treated better later.

			! check for hyperbolic regions

		  Alf = malfven2(j,k)	! Alfv�n square
		  plocg = gamma*Sofpsi(psi(j,k),i_zone)*rho(j,k)**gamma
		  bpol = sqrt(bx(j,k)**2+bz(j,k)**2)
		  bfield=sqrt(bx(j,k)**2+by(j,k)**2+bz(j,k)**2)
		  hyper(j,k) = - ( Alf*(bfield**2+plocg*mu_mag)-plocg*mu_mag)  &
				/(Alf**2*bpol**2 - Alf*(bfield**2+plocg*mu_mag)+plocg*mu_mag)

			continue

       endif

	   end do
    end do

	! need to calculate the current at the end, because it requires the fields (in particular, b_phi)
    do k = 1, nz
       do j = 1, nx

		 if(sort_grid(j,k)>=1) then

          ! Calculate the position of the (j,k) grid point
          x = x_coord(j)

		 ! now the other components of the current
		 ! first J_R

          if( (k == 1).or.(sort_grid(j,k-1)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_x(j,k) = (by(j,k+1) - by(j,k))/dz_a(k)
          else if( (k == nz).or.(sort_grid(j,k+1)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_x(j,k) = (by(j,k) - by(j,k-1))/dz_a(k)
          else
             ! use a centered difference
			 j_x(j,k) = ( dz_a(k-1)**2*by(j,k+1) +  &
						 (dx_a(k)**2-dx_a(k-1)**2)*by(j,k) -  &
						 dz_a(k)**2*by(j,k-1) ) /  &
						 (dz_a(k)*dz_a(k-1)*(dz_a(k)+dz_a(k-1)))
          end if

		 j_x(j,k) =  - j_x(j,k) / mu_mag ! WARNING: the derivatives are taken with the PLUS sign for convenience!

		 ! then J_Z

          if( (j == 1).or.(sort_grid(j-1,k)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_z(j,k) = (by(j+1,k) - by(j,k))/dx_a(j)
          else if( (j == nx).or.(sort_grid(j+1,k)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_z(j,k) = (by(j,k) - by(j-1,k))/dx_a(j)
          else
             ! use a centered difference
			 j_z(j,k) = ( dx_a(j-1)**2*by(j+1,k) +  &
						 (dx_a(j)**2-dx_a(j-1)**2)*by(j,k) -  &
						 dx_a(j)**2*by(j-1,k) ) /  &
						 (dx_a(j)*dx_a(j-1)*(dx_a(j)+dx_a(j-1)))
          end if

		 j_z(j,k) = (j_z(j,k) + by(j,k)/x) / mu_mag

		 ! then use the various components of the current to get J_par

		 j_par(j,k) = ( j_x(j,k)*bx(j,k) + j_z(j,k)*bz(j,k) + j_phi(j,k)*by(j,k) ) /  &
						dsqrt(bx(j,k)**2+by(j,k)**2+bz(j,k)**2)

		endif

	  enddo
	enddo

	! calculate staggered-grid stuff after the rest (i,j) -> (i-1/2,j-1/2)

    do k = 2, nz-1
       do j = 2, nx-1
		  ! October 8 2021 Added further differentiation for bc_type==7 and external region with plasma

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(j,k)-2)
				i_zone = max(-1,i_zone)
			endif

			x = x_coord(j)

		if( (((tri_type/=13).AND.(bc_type/=7)).AND.((psi(j,k)/psic)<fraction))  &
							.OR.  &
!				((bc_type==7).AND.(sqrt((x_coord(j)-rmajor)**2+z_coord(k)**2)>0.4d0*z_size))  &
!							.OR.  &
				((tri_type==13).AND.(sort_grid(j,k)==1)) ) then

				cycle

			endif

			psi0 = psi(j,k)
			if(j>1) psil = psi(j-1,k)
			if(j<nx) psir = psi(j+1,k)
			if(k>1) psid = psi(j,k-1)
			if(k<nz) psiu = psi(j,k+1)

			drho_ds = van_Leer_slope_new(rho(j-1,k),rho(j,k),rho(j+1,k),dx_a(j-1),dx_a(j))

			phirx = phiofpsi(0.5d0*(psir + psi0),i_zone)
			if(jump_option==-1) then
				rhorx = (rho(j,k) + rho(j+1,k))/2.d0
			else
				rhorx = rho(j,k) + 0.5d0*dx_a(j)*drho_ds
			endif
			ma2rx = phirx*phirx/rhorx

			philx = phiofpsi(0.5d0*(psil + psi0),i_zone)
			if(jump_option==-1) then
				rholx = (rho(j,k) + rho(j-1,k))/2.d0
			else
				rholx = rho(j,k) - 0.5d0*dx_a(j-1)*drho_ds
			endif
			ma2lx = philx*philx/rholx

			! -----------------------------------------------------

			drho_ds = van_Leer_slope_new(rho(j,k-1),rho(j,k),rho(j,k+1),dz_a(k-1),dz_a(k))

			phirz = phiofpsi(0.5d0*(psi0 + psiu),i_zone)
			if(jump_option==-1) then
				rhorz = (rho(j,k)+rho(j,k+1))/2.d0
			else
				rhorz = rho(j,k) + 0.5d0*dz_a(k)*drho_ds
			endif
			ma2rz = phirz*phirz/rhorz

			philz = phiofpsi(0.5d0*(psi0 + psid),i_zone)
			if(jump_option==-1) then
				rholz = (rho(j,k)+rho(j,k-1))/2.d0
			else
				rholz = rho(j,k) - 0.5d0*dz_a(k)*drho_ds
			endif
			ma2lz = philz*philz/rholz

	  enddo
	enddo


	continue

  end subroutine physical_var


  ! ------------------------------------------------------------------
  ! These routine calculates the various physical
  ! quantities such as the density, pressure, velocity, etc.
  ! Recall that (x, y, z) -> (r, phi, z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine physical_var_TF(psi,big_Psi,psi_diff,n_den,nx,nz,vx,vy,vz,bx,by,bz,p,  &
							p_e,p_i,malfven2,mpol,beta,beta_e,beta_i,j_phi,temp,temp_e,temp_i,mtor,cs,csp,  &
							j_x, j_z, j_par)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a
	use pseudo_IMSL, only : DBSNAK, DBS2VL, DBS2IN

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi,big_Psi,psi_diff,n_den
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: &
         vx,vy,vz,bx,by,bz,p,p_e,p_i,malfven2,mpol,beta,  &
		 beta_e,beta_i,temp,temp_e,temp_i,mtor,j_phi,cs,csp, j_x, j_z, j_par
    ! malfven_2 is the square of the alfvenic mach number
    real (kind=dkind) :: phic ! Temp Var. for phiofpsi(psi(i,j))
    real (kind=dkind) :: omegac ! Temp Var. for omegaofpsi(psi(i,j))
    real (kind=dkind) :: vp2 ! Square of the poloidal velocity
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: dx2,dz2
    integer :: i,j
    real (kind=dkind) :: x,rloc,bfield,delstar
    real (kind=dkind) :: cf2,cs2,a2,ca2,ct2,cusp2
	real(kind=dkind), dimension(-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: term1,term2,term3,term4,term5,b_dot_v, term6
	real(kind=dkind) :: Alf, plocg, bpol
    real (kind=dkind) :: phirx,philx,phirz,philz
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: drho_ds, dpsidx, dpsidz, Fstar
	integer :: i_zone

    dx2 = dx*dx
    dz2 = dz*dz

	i_zone = 0 ! This will only be changed for bc_type==7

    ! Calculate the poloidal magnetic field
    do j = 1, nz
       do i = 1, nx

		 if(sort_grid(i,j)>=1) then

          ! Calculate the position of the (i,j) grid point
          x = x_coord(i)

			call psi_derivative(i,j,nx,nz,psi,dpsidx,dpsidz)

			bx(i,j) = -dpsidz/x
			bz(i,j) = dpsidx/x

		  ! -------------------------------------for the time being, differentiate for external points-------------------------------------

		  ! -------------------------------------proceed from here for inner points-------------------------------------
		  ! October 8 2021 Added further differentiation for bc_type==7 and external region with plasma

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(i,j)-2)
				i_zone = max(-1,i_zone)
			endif

			call psi_derivative(i,j,nx,nz,big_Psi,dpsidx,dpsidz)

          ! Initialize phic
          phic = phi_TF_ofpsi(big_Psi(i,j),i_zone)

          ! Now let's calculate the poloidal velocity components
          vx(i,j) = -phic*dpsidz/n_den(i,j)/x
          vz(i,j) = phic*dpsidx/n_den(i,j)/x

          ! Next we calculate B_phi
		  Fstar = Fstarofpsi(psi(i,j),big_Psi(i,j),i_zone)
			by(i,j) = Fstar/x

          ! Next we calculate V_phi
          vy(i,j) = eV/mass * psi_diff(i,j)/ x

          ! Finally we calculate the pressure (total pressure, ion and electron will follow)
			p_e(i,j) = mass*s_e_ofpsi(psi(i,j),i_zone)*n_den(i,j)**gamma_e
			p_i(i,j) = mass*s_i_ofpsi(big_Psi(i,j),i_zone)*n_den(i,j)**gamma_i
			p(i,j) = p_e(i,j)+p_i(i,j)
          ! Calculate the plasma beta
		  beta_e(i,j) = 2.0d0*mu_mag*p_e(i,j)/(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
		  beta_i(i,j) = 2.0d0*mu_mag*p_i(i,j)/(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
		  beta(i,j) = beta_i(i,j) + beta_e(i,j)

		  ! And the plasma temperature
			temp_e(i,j) = p_e(i,j)/n_den(i,j)/eV
			temp_i(i,j) = p_i(i,j)/n_den(i,j)/eV
			temp(i,j) = p(i,j)/n_den(i,j)/eV

          ! Now the toroidal current density

		 if ((i>1).AND.(j>1).AND.(i<nx).AND.(j<nz)) then

!!$				b_dot_v = x*omegaofpsi(psi(i,j))*by(i,j) + (phiofpsi(psi(i,j))/rho(i,j))*  &
!!$						( bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2 )/sqrt(mu_mag)
!!$
!!$				term1= b_dot_v*dphidpsi(psi(i,j))/dsqrt(mu_mag)
!!$				term2= x*(phiofpsi(psi(i,j))*by(i,j)/dsqrt(mu_mag) +  &
!!$						rho(i,j)*x*omegaofpsi(psi(i,j)))*domegadpsi(psi(i,j))
!!$				term3= by(i,j)*didpsi(psi(i,j))/x/dsqrt(mu_mag)
!!$				term4= rho(i,j)*dhdpsi(psi(i,j))
!!$				term5= rho(i,j)**gamma*dsdpsi(psi(i,j))/(gamma-1.0d0)
!!$
!!$				delstar = term1+term2+term3+term4-term5
!!$
!!$				term6 = (psi(i+1,j)-psi(i-1,j)) *  &
!!$								(phiofpsi(psi(i+1,j))**2/rho(i+1,j)-phiofpsi(psi(i-1,j))**2/rho(i-1,j)) / (4.d0*dx2) +  &
!!$							(psi(i,j+1)-psi(i,j-1)) *  &
!!$								(phiofpsi(psi(i,j+1))**2/rho(i,j+1)-phiofpsi(psi(i,j-1))**2/rho(i,j-1)) / (4.d0*dz2)
!!$				term6 = term6/x**2

			j_phi(i,j) = ( -( (psi(i+1,j)-psi(i,j))/(x+dx/2.d0) - (psi(i,j)-psi(i-1,j))/(x-dx/2.d0) )/dx2  &
										- (psi(i,j+1)-2.d0*psi(i,j)+psi(i,j-1))/dz2/x ) / mu_mag
			!!$x * ( delstar - term6) / (1.d0-phic**2/rho(i,j))
			! "+" is the correct sign

		 endif

		  ! Now let's calculate some Mach numbers
          ! First some intermediate quantities

          ! Rather than using the poloidal velocities just calculated
          ! which involve differentiation, let's use the Bernoulli
          ! equation to determine the square of the poloidal velocity
          ! which is really all we need for the Mach numbers
			! for the time being, do it the dumb way
			vp2 = vx(i,j)**2 + vz(i,j)**2
 

          ! Sound Speed
		  a2 = gamma*p(i,j)/(mass*n_den(i,j))
          ! Poloidal Alfven Speed Squared
		  ca2 = (bx(i,j)*bx(i,j) + bz(i,j)*bz(i,j))/mu_mag/(mass*n_den(i,j))
!!$          if(phic .eq. 0.0d0) then
!!$             ca2 = 
!!$          else
!!$             ca2 = rho(i,j)*vp2/(mu_mag*phic**2)
!!$          endif
          ! Toroidal Alfven Speed Squared
          ct2 = by(i,j)*by(i,j)/mu_mag/(mass*n_den(i,j))
          ! Fast Mode Speed Squared
          cf2 = 0.5d0*(a2 + ca2 + ct2 &
               + dsqrt((a2 + ct2 - ca2)**2 + 4.0d0*ca2*ct2) )
          ! Slow Mode Speed Squared
          cs2 = a2*ca2/cf2

!!$			print*, i, j
!!$			print*, vp2, a2, ca2, ct2, cf2, cs2

          ! Calculate the square of the Alfvenic Mach Number
          malfven2(i,j) = mu_mag*mass * phic**2/n_den(i,j)
!          malfven2(i,j) = mu_mag*phic**2/rho(i,j)

          ! Now the slow mode Mach number
!!$          mslow(i,j) = dsqrt(malfven2(i,j)*cf2/a2)

!!$		  mslow(i,j) = sqrt(vp2/a2)*sqrt((bx(i,j)**2+by(i,j)**2+bz(i,j)**2)/(bx(i,j)**2+bz(i,j)**2))

          ! Now the cusp Mach number
!!$          mcusp(i,j) = mslow(i,j)*dsqrt(1.0d0 + cs2/cf2)

			mpol(i,j) = sqrt((vx(i,j)**2+vz(i,j)**2)/a2) *sqrt((bx(i,j)**2+by(i,j)**2+bz(i,j)**2)/(bx(i,j)**2+bz(i,j)**2))

		  ! Now the toroidal Mach number
          mtor(i,j) = dsqrt(vy(i,j)**2/a2)

		  ! sound speeds
		  cs(i,j) = sqrt(a2)
		  csp(i,j) = sqrt(a2*(bx(i,j)**2+bz(i,j)**2)/(bx(i,j)**2+by(i,j)**2+bz(i,j)**2))

!!$          ! Note: cs2 is only zero when the poloidal field goes to zero.
!!$          ! Since rho*vp = Phi * Bp and Phi = Phi(psi) is finite,
!!$          ! this implies that vp is also zero. The correct limiting
!!$          ! value should be treated better later.

			! check for hyperbolic regions

			continue

       endif

	   end do
    end do


	! need to calculate the current at the end, because it requires the fields (in particular, b_phi)
    do j = 1, nz
       do i = 1, nx

		 if(sort_grid(i,j)>=1) then

          ! Calculate the position of the (j,k) grid point
          x = x_coord(i)

		 ! now the other components of the current
		 ! first J_R

          if( (j == 1).or.(sort_grid(i,j-1)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_x(i,j) = (by(i,j+1) - by(i,j))/dz_a(j)
          else if( (j == nz).or.(sort_grid(i,j+1)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_x(i,j) = (by(i,j) - by(i,j-1))/dz_a(j)
          else
             ! use a centered difference
			 j_x(i,j) = ( dz_a(j-1)**2*by(i,j+1) +  &
						 (dx_a(j)**2-dx_a(j-1)**2)*by(i,j) -  &
						 dz_a(j)**2*by(i,j-1) ) /  &
						 (dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)))
          end if

		 j_x(i,j) =  - j_x(i,j) / mu_mag ! WARNING: the derivatives are taken with the PLUS sign for convenience!

		 ! then J_Z

          if( (i == 1).or.(sort_grid(i-1,j)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_z(i,j) = (by(i+1,j) - by(i,j))/dx_a(i)
          else if( (i == nx).or.(sort_grid(i+1,j)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_z(i,j) = (by(i,j) - by(i-1,j))/dx_a(i)
          else
             ! use a centered difference
			 j_z(i,j) = ( dx_a(i-1)**2*by(i+1,j) +  &
						 (dx_a(i)**2-dx_a(i-1)**2)*by(i,j) -  &
						 dx_a(i)**2*by(i-1,j) ) /  &
						 (dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)))
          end if

		 j_z(i,j) = (j_z(i,j) + by(i,j)/x) / mu_mag

		 ! then use the various components of the current to get J_par

		 j_par(i,j) = ( j_x(i,j)*bx(i,j) + j_z(i,j)*bz(i,j) + j_phi(i,j)*by(i,j) ) /  &
						dsqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)

		endif

	  enddo
	enddo

	continue

end subroutine physical_var_TF


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_trapped(nx, nz, psi, bpol, bphi, trapped)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer, intent(in) :: nx, nz

	real(kind=dkind), dimension(1:nx,1:nz) :: psi, bpol, bphi, trapped
	real(kind=dkind) :: psiloc, Bloc, thcrit, Bmax

	integer :: i,j

	do j = 1, nz
	do i = 1, nx

		if(sort_grid(i,j)<0) then

			trapped(i,j) = -1.d0
			cycle

		endif

		psiloc = psi(i,j)

		if(psiloc<modB_coef(1,1)) then

			trapped(i,j) = -1.d0
			cycle

		endif

		if(psiloc>modB_coef(1,enq)) then

			trapped(i,j) = 0.d0
			cycle

		endif

		Bmax = dbsval(psiloc,modB_ord, modB_coef(3,:),  &
					enq, modB_coef(4,1:enq) )

		Bloc = sqrt( bpol(i,j)**2 + bphi(i,j)**2 )

		if(Bloc>Bmax) then
		!this should not happen...

			thcrit = 0.d0
			continue

		else

			thcrit = sqrt(1.d0-Bloc/Bmax)

		endif

		trapped(i,j) = thcrit

	enddo
	enddo

	deallocate(modB_coef)

	continue

	return

end subroutine get_trapped


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine check_flow
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!This subroutine verifies whether there is any flow in the input

	real(kind=dkind) :: vmax	!maximum Mach number (only matters whether it's =0 or /=0)
	real(kind=dkind) :: x,y1,y2

	integer :: i

	vmax = 0.d0
	static_equi = .true.

	do i=0,201

		x=i/201.d0*psic
!		print*, i, x
		y1 = mach_theta(x)
!		print*, y1
		y2 = mach_phi(x)
!		print*, y2
		vmax = max(vmax,abs(y1),abs(y2))

		if(vmax>0.d0) then

			static_equi = .false.
			exit

		endif

	enddo

	continue

end subroutine check_flow


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psiout
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!This subroutine writes the most important functions of Psi
	!and their derivatives in a text format
	!this is mainly for debugging purposes, and will not
	! be used in the normal operation of the code

	integer i
	real (kind=dkind) x,y
	integer :: nsave = 201

!!$	return

! csv formatting and output added by Ian
! Certainly not an intelligent way to do this, but the code here is already somewhat redundent,
! so it shouldn't matter too much
46  format(e12.6,",",e12.6)

	!-----------------------------------------------------
	open(44,file='d.plt',status='unknown',action='write')
	open(45,file='d.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='dddpsi.plt',status='unknown',action='write')
	open(45,file='dddpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dddpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='pofpsi.plt',status='unknown',action='write')
	open(45,file='pofpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=pofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='dpdpsi.plt',status='unknown',action='write')
	open(45,file='dpdpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dpdpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='mach_phi.plt',status='unknown',action='write')
	open(45,file='mach_phi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=mach_phi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='mach_theta.plt',status='unknown',action='write')
	open(45,file='mach_theta.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=mach_theta(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

!	-----------------------------------------------------
	open(44,file='dmach_phidpsi.plt',status='unknown',action='write')
	open(45,file='dmach_phidpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dmach_phidpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

!	-----------------------------------------------------
	open(44,file='dmach_thetadpsi.plt',status='unknown',action='write')
	open(45,file='dmach_thetadpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dmach_thetadpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='bzero.plt',status='unknown',action='write')
	open(45,file='bzero.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=bzero(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='dbzerodpsi.plt',status='unknown',action='write')
	open(45,file='dbzerodpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dbzerodpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)
!
	!-----------------------------------------------------
	open(44,file='omega.plt',status='unknown',action='write')
	open(45,file='omega.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=omegaofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='domegadpsi.plt',status='unknown',action='write')
	open(45,file='domegadpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=domegadpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo

	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='h.plt',status='unknown',action='write')
	open(45,file='h.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=hofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='dhdpsi.plt',status='unknown',action='write')
	open(45,file='dhdpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dhdpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='s.plt',status='unknown',action='write')
	open(45,file='s.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=sofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='dsdpsi.plt',status='unknown',action='write')
	open(45,file='dsdpsi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=dsdpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='phi.plt',status='unknown',action='write')
	open(45,file='phi.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=phiofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='F.plt',status='unknown',action='write')
	open(45,file='F.csv',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=i*psic/nsave
		y=iofpsi(x)*sqrt(mu_mag)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	!-----------------------------------------------------

	if(eq_type>=10) then

		!-----------------------------------------------------
		open(44,file='p_e.plt',status='unknown',action='write')
		open(45,file='p_e.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=p_e_ofpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='p_i.plt',status='unknown',action='write')
		open(45,file='p_i.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=p_i_ofpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='h_e.plt',status='unknown',action='write')
		open(45,file='h_e.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=h_e_ofpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='he_partial.plt',status='unknown',action='write')
		open(45,file='he_partial.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=h_e_ofpsi_partial(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='h_i.plt',status='unknown',action='write')
		open(45,file='h_i.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=h_i_ofpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='hi_partial.plt',status='unknown',action='write')
		open(45,file='hi_partial.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=h_i_ofpsi_partial(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='dh_i_dpsi.plt',status='unknown',action='write')
		open(45,file='dh_i_dpsi.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dh_i_dpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='dhi_partial_dpsi.plt',status='unknown',action='write')
		open(45,file='dhi_partial_dpsi.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dh_i_dpsi_partial(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='dh_e_dpsi.plt',status='unknown',action='write')
		open(45,file='dh_e_dpsi.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dh_i_dpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='dh_e_partial_dpsi.plt',status='unknown',action='write')
		open(45,file='dh_e_partial_dpsi.csv',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dh_e_dpsi_partial(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

		!-----------------------------------------------------
		open(44,file='phi_TF.plt',status='unknown',action='write')
		open(45,file='phi_TF.csv',status='unknown',action='write')
	
		write(44,*)'TITLE="solution of GS equation with flow"'
	    write(44,*)'Variables =" psi " "var"'		!, 	fname
	    write(44,*)'ZONE I=', nsave+1,',F=Point'
	
		do i=0, nsave
			x=i*psic/nsave
			y=phi_TF_ofpsi(x)
			write(44,*) x/psic,y
			write(45,46) x/psic,y
		enddo
		close(44)
		close(45)

	endif

	open(44,file='d.dat',status='unknown',action='write')
	open(45,file='d.csv',status='unknown',action='write')

	do i=0,101
		x=i/101.d0*psic
		y=dofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	open(44,file='p.dat',status='unknown',action='write')
	open(45,file='p.csv',status='unknown',action='write')

	do i=0,101
		x=i/101.d0*psic
		y=pofpsi(x)
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)


	open(44,file='b.dat',status='unknown',action='write')
	open(45,file='b.csv',status='unknown',action='write')

	do i=0,101
		x=i/101.d0*psic
		y=bzero(x)*rmajor
		write(44,*) x/psic,y
		write(45,46) x/psic,y
	enddo
	close(44)
	close(45)

	continue

end subroutine psiout

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_integrations
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	call initialize_H

	allocate(x_int_phi(n_int_phi), fun_int_phi(n_int_phi), w_int_phi(n_int_phi))

	call set_weights(n_int_phi,0.d0,1.d0,w_int_phi,x_int_phi)

	return

end subroutine initialize_integrations

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_H
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, dx_H
	real(kind=dkind) :: p_i_term, p_e_term, v_term, int_term
	real(kind=dkind) :: int_temp, x_L, x_temp
	real(kind=dkind) :: ggm1_e, ggm1_i
	integer :: i, j

	allocate(H_e_data(n_H+H_e_ord,3))
	allocate(H_i_data(n_H+H_i_ord,3))
	allocate(H_e_cscoef(n_H))
	allocate(H_i_cscoef(n_H))
	allocate(Omega_int_data(n_H+Omega_int_ord,3))
	allocate(Omega_int_cscoef(n_H))

	ibreak_He = n_H + H_e_ord
	ibreak_Hi = n_H + H_i_ord
	ibreak_Om_int = n_H + Omega_int_ord

	allocate(x_int_H(n_int_H), fun_int_H(n_int_H), w_int_H(n_int_H))

	call set_weights(n_int_H,0.d0,1.d0,w_int_H,x_int_H)

	ggm1_e = gamma_e/(gamma_e-1.d0)
	ggm1_i = gamma_i/(gamma_i-1.d0)

	int_term = 0.d0

	dx_H = 1.d0/(n_H-1.d0)

	do i = 0, n_H-1

		x = 1.d0/(n_H-1.d0)*i

		p_e_term = p_e_ofpsi(x)/d_TF_ofpsi(x)
		p_i_term = p_i_ofpsi(x)/d_TF_ofpsi(x)
		v_term = gamma*(p_i_term+p_e_term)*mach_phi(x)
		p_i_term = ggm1_i * p_i_term
		p_e_term = ggm1_e * p_e_term

		if(i==0) then

			int_term = 0.d0

		else

			do j = 1, n_int_H

				x_temp = x_L + x_int_H(j)*dx_H
				fun_int_H(j) = omegaofpsi(x_temp)  &
											- dp_i_dpsi(x_temp) / (eV * d_TF_ofpsi(x_temp))

			enddo

			call integrate(n_int_H,fun_int_H,w_int_H,int_temp)

			int_term = int_term + int_temp*dx_H	! dx_H is of course the Jacobian

		endif

		H_e_data(i+1,1) = x
		H_e_data(i+1,2) = mass*p_e_term - eV*int_term
		H_i_data(i+1,1) = x
		H_i_data(i+1,2) = mass*(v_term + p_i_term) + eV*int_term
		Omega_int_data(i+1,1) = x
		Omega_int_data(i+1,2) = int_term

		x_L = x

	enddo

	call interp_setup(n_H,H_e_ord,H_e_data(1:n_H,1),H_e_data(1:n_H,2),H_e_data(:,3),H_e_cscoef)
	call interp_setup(n_H,H_i_ord,H_i_data(1:n_H,1),H_i_data(1:n_H,2),H_i_data(:,3),H_i_cscoef)
	call interp_setup(n_H,Omega_int_ord,Omega_int_data(1:n_H,1),Omega_int_data(1:n_H,2),Omega_int_data(:,3),Omega_int_cscoef)

	continue

end subroutine initialize_H

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psi_derivative(i,j,nx,nz,psi,dpsidx,dpsidz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx_a, dz_a

	implicit none

	integer, intent(in) :: i,j,nx,nz
	real(kind=dkind), intent(inout) :: dpsidx,dpsidz
	real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi


	! d Psi / d x
    if(i == 1) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif left"
       dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(1)
    else if(i == nx) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif right"
       dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(nx)
    else
       ! use a centered difference
	   dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
					(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
					dx_a(i)**2*psi(i-1,j) ) /  &
					( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
    end if
    ! d Psi / d z
    if(j == 1) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif left"
       dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
    else if(j == nz) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif right"
       dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
    else
       ! use a centered difference
		dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
						(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
						dz_a(j)**2*psi(i,j-1) ) /  &
						( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )


    end if

	continue

	end subroutine psi_derivative

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine psi_derivative_new(i,j,nx,nz,psi,dpsidx,dpsidz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx_a, dz_a

	implicit none

	integer, intent(in) :: i,j,nx,nz
	real(kind=dkind), intent(inout) :: dpsidx,dpsidz
	real (kind=dkind), dimension(-1:1,-1:1), intent(in) :: psi


	! d Psi / d x
    if(i == 1) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif left"
       dpsidx = (psi(1,0) - psi(0,0))/dx_a(1)
    else if(i == nx) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif right"
       dpsidx = (psi(0,0) - psi(-1,0))/dx_a(nx)
    else
       ! use a centered difference
!       dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
	   dpsidx = ( dx_a(i-1)**2*psi(1,0) +  &
					(dx_a(i)**2-dx_a(i-1)**2)*psi(0,0) -  &
					dx_a(i)**2*psi(-1,0) ) /  &
					( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
    end if
    ! d Psi / d z
    if(j == 1) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif left"
       dpsidz = (psi(0,1) - psi(0,0))/dz_a(j)
    else if(j == nz) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif right"
       dpsidz = (psi(0,0) - psi(0,1))/dz_a(j)
    else
       ! use a centered difference
!       dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
		dpsidz = ( dz_a(j-1)**2*psi(0,1) +  &
						(dz_a(j)**2-dz_a(j-1)**2)*psi(0,0) -  &
						dz_a(j)**2*psi(0,-1) ) /  &
						( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )


    end if

	continue

end subroutine psi_derivative_new

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine getqstar(psi,jphi,p,bpol,bphi,beta,rho,vr,vphi,vz,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, dx, dz, inv_aspect_ratio
	use triangularity, only : volume, area, R_edge_min, R_edge_max,  &
										z_edge_min, z_edge_max, shape_ellipticity
	use interpolating_functions, only : interp_setup
	use pseudo_IMSL, only : dbsval_safe

	implicit none

	integer, intent(in) :: n
	real (kind=dkind), intent(in), dimension(1:n,1:n) :: jphi, p, psi, bpol,bphi, beta, rho,  &
																					   vr,vphi,vz
	real (kind=dkind) :: b_phi_zero
	real (kind=dkind) :: surf,qstar,dA,ex,ez	!,curr
	real (kind=dkind) :: nu,vol,bigR,epsilon,rminor,  &
						 ptot, Itot, beta_LDX, beta_int, bphiav,  &
						 I_Sauter, I_NCLASS, I_Sauter_phi, I_NCLASS_phi,  &
						B_den, mass, kin_en
	real(kind=dkind) :: bs_loc
	integer :: i,j
	real(kind=dkind) :: R_l, R_r, angle
	real(kind=dkind) :: J_boot_interp(8,0:enq+1+ord_loc)
	real(kind=dkind), dimension(1:2) :: length, Bpint
	real(kind=dkind) :: psi_integral(2)

	! 7/7/2008 total BS current added:
	! set up J_BS interpolation first

	!2/4/2009 added diagnostic for plasma mass

	! Sauter calculation
	if (write_all.and.((bootstrap_option==0).or.(bootstrap_option==2))) then

		J_boot_interp(1,0) = 0.d0
		J_boot_interp(2,0) = J_boot(2,enq) - ( J_boot(2,enq-1)-J_boot(2,enq) ) /  &
							( psival(enq-1)-psival(enq) ) * psival(enq)

		J_boot_interp(1,enq+1) = psic
		J_boot_interp(2,enq+1) = 0.d0

		do i = 1,enq

			J_boot_interp(1,enq+1-i) = psival(i)
			J_boot_interp(2,enq+1-i) = J_boot(2,i)

		enddo

		call interp_setup(enq+2,ord_loc,J_boot_interp(1,0:enq+1),J_boot_interp(2,0:enq+1),  &
								J_boot_interp(3,0:enq+1+ord_loc),J_boot_interp(4,0:enq+1))

	endif

	! NCLASS calculation
	if (write_all.and.((bootstrap_option==1).or.(bootstrap_option==2))) then

		J_boot_interp(5,0) = 0.d0
		J_boot_interp(6,0) = J_boot(6,enq) - ( J_boot(6,enq-1)-J_boot(6,enq) ) /  &
							( psival(enq-1)-psival(enq) ) * psival(enq)

		J_boot_interp(5,enq+1) = psic
		J_boot_interp(6,enq+1) = 0.d0

		do i = 1,enq

			J_boot_interp(5,enq+1-i) = psival(i)
			J_boot_interp(6,enq+1-i) = J_boot(6,i)

		enddo

		call interp_setup(enq+2,ord_loc,J_boot_interp(5,0:enq+1),J_boot_interp(6,0:enq+1),  &
								J_boot_interp(7,0:enq+1+ord_loc),J_boot_interp(8,0:enq+1))

	endif

	! first, get the current and the cross section area


	dA = dx*dz

	curr = 0.d0
	qstar = 0.d0
	surf = 0.d0
	betator = 0.d0
	nu = 0.d0
	vol = 0.d0
	Itot = 0.d0
	beta_LDX = 0.d0
	beta_int = 0.d0
	bphiav = 0.d0
	bpolav = 0.d0
	I_Sauter = 0.d0
	I_NCLASS = 0.d0
	mass = 0.d0
	kin_en = 0.d0

	if((Broot==4).or.(Broot==5)) then
		B_den = bzero(psic)
	else
		B_den = bzero(0.d0)
	endif

	do j=1,n
		do i=1,n

		    if(sort_grid(i,j)>=1) then

				if((psi(i,j)/psic)>=fraction) then	!inner region

				bigR = x_coord(i) + dx_a(i) -  dx_a(i-1)

				if(grid_type/=0) dA = 0.25d0*(dx_a(i) + dx_a(i-1)) *  &
										(dz_a(j) + dz_a(j-1))


!				if(bigR>rmajor+a_elps/3.d0) cycle	!trucco

				surf = surf+dA
				vol = vol + bigR*dA
				curr = curr +jphi(i,j)*dA
				betator = betator+p(i,j)*bigR*dA
				betapol = betapol + bpol(i,j)*bigR*dA
				ptot = ptot + p(i,j)*bigR*dA
				beta_LDX = beta_LDX + bpol(i,j)**2*bigR*dA
				beta_int = beta_int + beta(i,j)**2*bigR*dA
				bphiav = bphiav + bphi(i,j)*dA
				bpolav = bpolav + bpol(i,j)*dA
				mass = mass  + rho(i,j)*bigR*dA
				kin_en = kin_en  + rho(i,j)*(vr(i,j)**2+vphi(i,j)**2+vz(i,j)**2)*bigR*dA

				if (write_all.and.((bootstrap_option==0).or.(bootstrap_option==2)))  then

					bs_loc = dA * dbsval_safe(psi(i,j),ord_loc,J_boot_interp(3,0:enq+1+ord_loc),enq+2,  &
								J_boot_interp(4,0:enq+1),J_boot_interp(2,0:enq+1))  &
								/ sqrt(bphi(i,j)**2+bpol(i,j)**2)

					I_Sauter = I_Sauter +  bs_loc

					I_Sauter_phi = I_Sauter_phi +  bs_loc * bphi(i,j)/sqrt(bphi(i,j)**2+bpol(i,j)**2)

				endif

				if (write_all.and.((bootstrap_option==1).or.(bootstrap_option==2)))  then

					bs_loc = dA * dbsval_safe(psi(i,j),ord_loc,J_boot_interp(7,0:enq+1+ord_loc),enq+2,  &
								J_boot_interp(8,0:enq+1),J_boot_interp(6,0:enq+1))  &
								/ sqrt(bphi(i,j)**2+bpol(i,j)**2)

					I_NCLASS = I_NCLASS +  bs_loc

					I_NCLASS_phi = I_NCLASS_phi + bs_loc * bphi(i,j)/sqrt(bphi(i,j)**2+bpol(i,j)**2)

				endif

				endif

			endif

		enddo
	enddo

	qstar = 2.d0*surf*bzero(0.d0)/(mu_mag*rmajor*curr)
	betastar = betator/vol * 2.d0*mu_mag/bzero(0.d0)**2 *(betator/vol)
	betator = betator/vol * 2.d0*mu_mag/B_den**2
	betapol = betapol/vol
	betapol = betator*B_den**2/betapol**2
	epsilon = (dsqrt(surf/pi)/rmajor)
	nu = betator*qstar**2/ epsilon
	beta_LDX = 2.d0*ptot*mu_mag/beta_LDX
	ptot = ptot*surf/vol
!	Itot = Itot / mu_mag
	bphiav = bphiav/surf
	bpolav = bpolav/surf

	call radius_theta(0.d0,rminor,R_r,ez) !ez not used here
	call radius_theta(pi,rminor,R_l,ez)

	inv_aspect_ratio = (R_r-R_l)/rmajor/2.d0

	R_edge_max = -1.d0
	z_edge_max = -1.d0
	R_edge_min = 1.d6*rmajor
	z_edge_min = 1.d6*rmajor

	do i = 1,r_of_theta_points

		angle = 2.d0*pi*i/r_of_theta_points
		call radius_theta(angle, rminor, ex, ez)

		if(ex>R_edge_max) R_edge_max = ex
		if(ex<R_edge_min) R_edge_min = ex
		if(ez>z_edge_max) z_edge_max = ez
		if(ez<z_edge_min) z_edge_min = ez

	enddo

	shape_ellipticity = (z_edge_max-z_edge_min)/(R_edge_max-R_edge_min)

	if((Broot==0).and.(jump_option<=-5)) call J_jump(psi,bpol,n,length,Bpint,psi_integral)

	print*, 'qstar = ',qstar
!	print*, 'nu = ', nu
	print*, 'plasma current = ',curr
	print*, 'betator = ',betator
	print*, 'betapol = ',betapol
	print*, 'epsilon = ',epsilon

	open(99,file='qstar.out',status='unknown',action='write')

	write(99,*) 'qstar = ',qstar
	write(99,*) 'surf = ',surf
	write(99,*) 'volume = ',vol
	write(99,*) 'plasma current = ',curr
	write(99,*) 'nu = ', nu
	write(99,*) 'betator = ',betator
	write(99,*) 'betapol = ',betapol
	write(99,*) 'average epsilon = ',epsilon
	write(99,*) 'inverse aspect ratio = ',inv_aspect_ratio
	write(99,*) 'total pressure * area = ', ptot
!	write(99,*) 'total current = ', Itot
	write(99,*) 'beta_LDX = ', beta_LDX
	write(99,*) 'B_phi_average = ', bphiav
	write(99,*) 'B_pol_average = ', bpolav
	write(99,*) 'mass = ', mass
	write(99,*) 'kinetic energy = ', kin_en
	write(99,*) '               '
	write(99,*) 'I_BS_Sauter = ', I_Sauter
	write(99,*) 'I_BS_NCLASS = ', I_NCLASS
	write(99,*) 'I_BS_Sauter_phi = ', I_Sauter_phi
	write(99,*) 'I_BS_NCLASS_phi = ', I_NCLASS_phi
	write(99,*) 'f_BS_Sauter = ', I_Sauter_phi/curr
	write(99,*) 'f_BS_NCLASS = ', I_NCLASS_phi/curr

	if((Broot==0).and.(jump_option<=-5)) then

		write(99,*) 'Bpol line integral (in, out):', Bpint(1), Bpint(2)
		write(99,*) 'length line integral (in, out):', length(1), length(2)
		write(99,*) 'Bpol average line integral (in, out):', Bpint(1)/length(1), Bpint(2)/length(2)
		write(99,*) 'integral locations:', psi_integral(1), psi_integral(2)

	endif


	close(99)

	if((broot==4).or.(broot==5)) then

		open(99,file='B_av.out',status='unknown',action='write')
		write(99,*) bphiav
		write(99,*) bpolav
		close(99)

	endif

	volume = vol
	area = surf

	return

end subroutine getqstar


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine getqstar_2F(psi,jphi,p,p_i,p_e,bpol,bphi,beta,rho,vr,vphi,vz,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, dx, dz, inv_aspect_ratio
	use triangularity, only : volume, area, R_edge_min, R_edge_max,  &
										z_edge_min, z_edge_max, shape_ellipticity
	use interpolating_functions, only : interp_setup
	use pseudo_IMSL, only : dbsval_safe

	implicit none

	integer, intent(in) :: n
	real (kind=dkind), intent(in), dimension(1:n,1:n) :: jphi, p, psi, bpol, bphi, beta, rho,  &
																					   vr, vphi, vz, p_i, p_e
	real (kind=dkind) :: b_phi_zero
	real (kind=dkind) :: surf,qstar,dA,ex,ez	!,curr
	real (kind=dkind) :: nu,vol,bigR,epsilon,rminor,  &
						 ptot, Itot, beta_int, bphiav,  &
						B_den, mass, kin_en, betator_i, betator_e
	real(kind=dkind) :: bs_loc
	integer :: i,j
	real(kind=dkind) :: R_l, R_r, angle
	real(kind=dkind), dimension(1:2) :: length, Bpint
	real(kind=dkind) :: psi_integral(2)

	! 7/7/2008 total BS current added:
	! set up J_BS interpolation first

	!2/4/2009 added diagnostic for plasma mass

	! first, get the current and the cross section area


	dA = dx*dz

	curr = 0.d0
	qstar = 0.d0
	surf = 0.d0
	betator = 0.d0
	betator_i = 0.d0
	betator_e = 0.d0
	nu = 0.d0
	vol = 0.d0
	Itot = 0.d0
	beta_int = 0.d0
	bphiav = 0.d0
	bpolav = 0.d0
	mass = 0.d0
	kin_en = 0.d0

	if((Broot==4).or.(Broot==5)) then
		B_den = bzero(psic)
	else
		B_den = bzero(0.d0)
	endif

	do j=1,n
		do i=1,n

		    if(sort_grid(i,j)>=1) then

				if((psi(i,j)/psic)>=fraction) then	!inner region

					bigR = x_coord(i) + dx_a(i) -  dx_a(i-1)

					if(grid_type/=0) dA = 0.25d0*(dx_a(i) + dx_a(i-1)) *  &
											(dz_a(j) + dz_a(j-1))


	!				if(bigR>rmajor+a_elps/3.d0) cycle	!trucco

					surf = surf+dA
					vol = vol + bigR*dA
					curr = curr +jphi(i,j)*dA
					betator = betator + p(i,j)*bigR*dA
					betator_i = betator_i + p_i(i,j)*bigR*dA
					betator_e = betator_e + p_e(i,j)*bigR*dA
					betapol = betapol + bpol(i,j)*bigR*dA
					ptot = ptot + p(i,j)*bigR*dA
					beta_int = beta_int + beta(i,j)**2*bigR*dA
					bphiav = bphiav + bphi(i,j)*dA
					bpolav = bpolav + bpol(i,j)*dA
					mass = mass  + rho(i,j)*bigR*dA
					kin_en = kin_en  + rho(i,j)*(vr(i,j)**2+vphi(i,j)**2+vz(i,j)**2)*bigR*dA

				endif

			endif

		enddo
	enddo

	qstar = 2.d0*surf*bzero(0.d0)/(mu_mag*rmajor*curr)
	betastar = betator/vol * 2.d0*mu_mag/bzero(0.d0)**2 *(betator/vol)
	betator = betator/vol * 2.d0*mu_mag/B_den**2
	betator_i = betator_i/vol * 2.d0*mu_mag/B_den**2
	betator_e = betator_e/vol * 2.d0*mu_mag/B_den**2
	betapol = betapol/vol
	betapol = betator*B_den**2/betapol**2
	epsilon = (dsqrt(surf/pi)/rmajor)
	nu = betator*qstar**2/ epsilon
	ptot = ptot*surf/vol
	bphiav = bphiav/surf
	bpolav = bpolav/surf

	call radius_theta(0.d0,rminor,R_r,ez) !ez not used here
	call radius_theta(pi,rminor,R_l,ez)

	inv_aspect_ratio = (R_r-R_l)/rmajor/2.d0

	R_edge_max = -1.d0
	z_edge_max = -1.d0
	R_edge_min = 1.d6*rmajor
	z_edge_min = 1.d6*rmajor

	do i = 1,r_of_theta_points

		angle = 2.d0*pi*i/r_of_theta_points
		call radius_theta(angle, rminor, ex, ez)

		if(ex>R_edge_max) R_edge_max = ex
		if(ex<R_edge_min) R_edge_min = ex
		if(ez>z_edge_max) z_edge_max = ez
		if(ez<z_edge_min) z_edge_min = ez

	enddo

	shape_ellipticity = (z_edge_max-z_edge_min)/(R_edge_max-R_edge_min)

	if((Broot==0).and.(jump_option<=-5)) call J_jump(psi,bpol,n,length,Bpint,psi_integral)

	print*, 'qstar = ',qstar
!	print*, 'nu = ', nu
	print*, 'plasma current = ',curr
	print*, 'betator = ',betator
	print*, 'betapol = ',betapol
	print*, 'epsilon = ',epsilon

	open(99,file='qstar.out',status='unknown',action='write')

	write(99,*) 'qstar = ',qstar
	write(99,*) 'surf = ',surf
	write(99,*) 'volume = ',vol
	write(99,*) 'plasma current = ',curr
	write(99,*) 'nu = ', nu
	write(99,*) 'betator = ',betator
	write(99,*) 'betapol = ',betapol
	write(99,*) 'betator_i = ',betator_i
	write(99,*) 'betator_e = ',betator_e
	write(99,*) 'average epsilon = ',epsilon
	write(99,*) 'inverse aspect ratio = ',inv_aspect_ratio
	write(99,*) 'total pressure * area = ', ptot
	write(99,*) 'B_phi_average = ', bphiav
	write(99,*) 'B_pol_average = ', bpolav
	write(99,*) 'mass = ', mass
	write(99,*) 'kinetic energy = ', kin_en
	write(99,*) 'Mach theta max = ', mach_theta_max
	write(99,*) 'inital Mach theta max = ', mach_theta_max_initial
	write(99,*) '               '


	close(99)

	if((broot==4).or.(broot==5)) then

		open(99,file='B_av.out',status='unknown',action='write')
		write(99,*) bphiav
		write(99,*) bpolav
		close(99)

	endif

	volume = vol
	area = surf

	return

end subroutine getqstar_2F


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine J_jump(psi,bpol,n,length,Bpint,psival)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! calculates the integral of Bp dl and dl on the two sides of the discontinuity
! (for jump and current spike calculation)

	use constant, only : dkind, dx, dz
	use interpolating_functions, only : lin_interp_2D

	implicit none

	integer, intent(in) :: n
	real (kind=dkind), intent(in), dimension(1:n,1:n) :: psi, bpol
	real(kind=dkind), dimension(1:2) :: length, Bpint
	real(kind=dkind) :: psival(2)
	real(kind=dkind) :: xx(2), zz(2)
	real(kind=dkind) :: dl
	real(kind=dkind) :: xQ, zQ, Bploc
	real(kind=dkind), dimension(1:4) :: psi4, bp4
	real(kind=dkind) :: psi_cmax, psi_cmin
	integer :: i, j, k, h

	length = 0.d0
	Bpint = 0.d0

	psival(1) = psi_degen + 0.01 !(inner)
	psival(2) = psi_degen - 0.01 !(outer)

	psival(1) = psi_degen + sqrt(1.d0/delta_Bern)*psic*2.d0 !(inner)
	psival(2) = psi_degen - sqrt(1.d0/delta_Bern)*psic*2.d0 !(outer)

	do j=1,n
	do i=1,n

		if(sort_grid(i,j)<=0) cycle

		psi4(1) = psi(i,j)
		psi4(2) = psi(i+1,j)
		psi4(3) = psi(i+1,j+1)
		psi4(4) = psi(i,j+1)

		psi_cmax = maxval(psi4)
		psi_cmin = minval(psi4)

		do k = 1, 2

			h = 1

			if((psi_cmax-psival(k))*(psi_cmin-psival(k))<0.d0) then
			! useful cell

				! check sides for intersectino points

				if((psi4(1)-psival(k))*(psi4(2)-psival(k))<0.d0) then
					xx(h) = x_coord(i) + (psival(k)-psi4(1))/(psi4(2)-psi4(1))*dx
					zz(h) = z_coord(j)
					h = h+1
				endif
				if((psi4(2)-psival(k))*(psi4(3)-psival(k))<0.d0) then
					xx(h) = x_coord(i+1)
					zz(h) = z_coord(j) + (psival(k)-psi4(2))/(psi4(3)-psi4(2))*dz
					h = h+1
				endif
				if((psi4(3)-psival(k))*(psi4(4)-psival(k))<0.d0) then
					xx(h) = x_coord(i+1) + (psival(k)-psi4(3))/(psi4(4)-psi4(3))*(-dx)
					zz(h) = z_coord(j+1)
					h = h+1
				endif
				if((psi4(4)-psival(k))*(psi4(1)-psival(k))<0.d0) then
					xx(h) = x_coord(i)
					zz(h) = z_coord(j+1) + (psival(k)-psi4(4))/(psi4(1)-psi4(4))*(-dz)
					h = h+1
				endif

				if(h<3) then
				! something dumb happened
					continue
					print*, 'error in J_jump', i, j
				endif

				! midpoint of segment
				xQ = sum(xx)/2.d0
				zQ = sum(zz)/2.d0

				dl = sqrt((xx(1)-xx(2))**2+(zz(1)-zz(2))**2)

				bp4(1) = bpol(i,j)
				bp4(2) = bpol(i+1,j)
				bp4(3) = bpol(i+1,j+1)
				bp4(4) = bpol(i,j+1)

				call lin_interp_2D(bp4,xQ,zQ,Bploc)

				length(k) = length(k) + dl
				Bpint(k) = Bpint(k) + Bploc*dl

			endif

		enddo

	enddo
	enddo

end subroutine J_jump


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
subroutine magnetic_output(psi_all,bpol,bphi,rho,bx,bz,vx,vz,xmax,zmax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

    use pseudo_IMSL, only : DBS2IN, DBS2VL, DBSNAK, DBSINT, DBSDER

    real(kind=dkind), dimension(1:n,1:n) :: psi_all, bpol, bphi, psi, rho, bx, bz, vx, vz

    integer :: enqbis
    real(kind=dkind), dimension(:), allocatable :: Rknot
    real(kind=dkind) :: psimax, delta, rguessl, rguessr, rguessl0, rguessr0
    real(kind=dkind) :: psiloc, thetaloc, rloc,  rlocprim, bigrloc, zloc,  &
		bpolloc, bphiloc, B_max, B_min, rholoc, vxloc, vzloc, bxloc, bzloc,  &
		vpolloc, vpsiloc, v_ratio_loc
	integer :: i, j,k
    integer :: nsurfdim
    integer :: n_int
    real(kind=dkind), dimension(:), allocatable :: w_int, t_int, integrand
    real(kind=dkind), dimension(:), allocatable :: w_int_boot, t_int_boot
    real(kind=dkind) :: a_q, b_q, c_q
    real(kind=dkind) :: xmax, zmax
    real(kind=dkind), dimension(:,:), allocatable :: bscoef_vx, bscoef_vz, bscoef_bx, bscoef_bz
	real(kind=dkind), dimension(:,:), allocatable :: psi_of_q_coef

    integer :: itemax = 200
    real(kind=dkind) :: xtoll = 1.d-6
    real(kind=dkind) :: ftoll = 1.d-6
    real(kind=dkind) :: tol ! for max |B|
    integer :: error

    real(kind=dkind), dimension(:,:), allocatable :: viscosity_logs
    real(kind=dkind), dimension(:,:), allocatable :: poloidal_viscosity
    real(kind=dkind) :: rho_0, B_0	! for normalization

    ! first, allocate the arrays that would have been automatic,
    ! if not for the complaints of the LINUX compiler

    allocate(bscoef_psi(1:n,1:n))
    allocate(bscoef_bpol(1:n,1:n))
    allocate(bscoef_bphi(1:n,1:n))
    allocate(bscoef_rho(1:n,1:n))

	if(n_q_vnorm>0) then

		allocate(bscoef_bx(1:n,1:n))
		allocate(bscoef_bz(1:n,1:n))
		allocate(bscoef_vx(1:n,1:n))
		allocate(bscoef_vz(1:n,1:n))

	endif

    allocate(xknot_mag(1:n+ord_loc))
    allocate(zknot_mag(1:n+ord_loc))

    xmax_mag = xmax
    zmax_mag = zmax

    rho_0 = dofpsi(0.d0)
    B_0 = abs(bzero(0.d0))

    !before starting, clean up psi to avoid problems with the interpolation

    do j=1,n
       do i=1,n

          if(sort_grid(i,j) < 0) then

             psi(i,j) = -psic/10.d0

          else

             psi(i,j) = psi_all(i,j)

          endif

       enddo
    enddo


    !first set up the interpolations

    call DBSNAK(n,x_coord(1:n),ord_loc,xknot_mag)
    call DBSNAK(n,z_coord(1:n),ord_loc,zknot_mag)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),psi,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_psi)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),bpol,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_bpol)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),bphi,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_bphi)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),rho,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_rho)

	if(n_q_vnorm>0) then

		call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),bx,n,  &
			 ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_bx)

		call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),bz,n,  &
			 ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_bz)

		call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),vx,n,  &
			 ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_vx)

		call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),vz,n,  &
			 ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_vz)

	endif


	enq = (n-1)*3/8

    ! then set up the psi values

    allocate(psival(0:enq+1))
    allocate(qval(enq))
    allocate(Rleft(enq))
    allocate(Rright(enq))
!    allocate(L_ratio(enq))

    psimax = DBS2VL(xmax,zmax,ord_loc,ord_loc,xknot_mag, &
         zknot_mag,n,n,bscoef_psi)

    delta = psimax/enq

    psival(0) = psimax

    psival(1) = psimax * ( 0.99d0 + 2.5d-3*log((n-1.d0)/128.d0)/log(2.d0) )

    if(psival(1)>=psimax) psival(1) = psimax*(1.d0-1.d-4)

    do i = 2,enq
       psival(i) = psimax - delta*(i-1.d0)
    enddo

    psival(enq+1) = 0.d0

    ! then get the magnetic surfaces and q values (surface by surface)

    nsurfdim = nsurf+1+ord_surf

    n_int = nsurf/2

    allocate(rtab(enq,2,nsurfdim))
    allocate(bigrtab(enq,2,nsurfdim))
    allocate(ztab(enq,2,nsurfdim))
    allocate(thetatab(2,nsurfdim))

    allocate(viscosity_logs(2,nsurfdim))
    allocate(poloidal_viscosity(enq,n_int))

    ! if the fraction of trapped particles is also required, we need to start the setup from here

    if(write_all) then

       ibreak_modB = enq + modB_ord
       allocate(modB_coef(4,ibreak_modB))

       tol = sqrt(epsilon(1.d0))

    endif

    ! if the complete output is desired,
    ! here we also set up the bootstrap current calculation,
    ! starting from the effective trapped particle fraction

    if(write_all) then

       if(bootstrap_option==0) then
          allocate(J_boot(5,1:enq))
          !		else	if(bootstrap_option==1) then
       elseif(bootstrap_option==2) then
          allocate(J_boot(7,1:enq))
       endif

       allocate(el_resistivity(1:3,1:enq))
       el_resistivity= 0.d0

       allocate(eff_trap(0:enq+1))
       allocate(surf_length(enq))
       allocate(B2_hat_ave(enq))
       allocate(boot_ne(enq))
       allocate(boot_pe(enq))
       allocate(boot_Te(enq))
       allocate(boot_ni(enq))
       allocate(boot_pi(enq))
       allocate(boot_Ti(enq))
       allocate(boot_neprim(enq))
       allocate(boot_peprim(enq))
       allocate(boot_Teprim(enq))
       allocate(boot_niprim(enq))
       allocate(boot_piprim(enq))
       allocate(boot_Tiprim(enq))
       allocate(inv_asp_ratio(enq))
       allocate(boot_Bmin(enq))
       allocate(boot_Bp_at_Bmin(enq))
       allocate(boot_R_of_Bmin(enq))
       allocate(boot_Rcenter(enq))
       allocate(JparB_ave(1:enq))

       if(bootstrap_option>=1) then

          allocate(boot_tor_flux(0:enq+1))
          allocate(boot_tor_rho(0:enq+1))
          allocate(boot_fhat(0:enq+1))
          allocate(boot_grho_ov_B2(0:enq+1))
          allocate(B2_ave(0:enq+1))
          allocate(Bm2_ave(0:enq+1))

       endif

    endif

    ! this sets up the 1D grid and the integration routine

    do j = 1, nsurf+1
       thetatab(1,j) = (j-1.d0)*pi*2.d0/nsurf
    enddo

    call DBSNAK(nsurf+1, thetatab(1,1:nsurf+1),  &
         ord_surf,thetatab(2,1:nsurfdim))

    allocate(w_int(1:n_int))
    allocate(t_int(1:n_int))

    call set_weights(n_int,0.d0,2.d0*pi,w_int,t_int)

    if(write_all) then

       allocate(w_int_boot(1:n_int))
       allocate(t_int_boot(1:n_int))

       call set_weights(n_int,0.d0,1.d0,w_int_boot,t_int_boot)

    endif

    rguessr0 = max(x_coord(n)-xmax, xmax-x_coord(1), z_coord(n)-zmax, zmax-z_coord(1))
    rguessl0 = (x_coord(n)-x_coord(1)+z_coord(n)-z_coord(1))/5.d2

    open(69,file='poloidal_viscosity.plt',action='write')
	open(70,file='poloidal_viscosity.csv',action='write')

    write(69,*)'TITLE="poloidal viscosity for solution of GS equation with flow"'
    write(69,*)'Variables ="theta", "r", "R [m] ","z [m]", "viscosity"'
    write(69,*)'ZONE I=',n_int,',J=',enq,',F=Point'

    do i = 1, enq

       surf_index = i

       rguessl = rguessl0
       rguessr = rguessr0

       psiloc = psival(i)

       ! this cycle determines the (R,Z) coordinates of each magnetic surface
       do j = 1, nsurf+1

          thetaloc = thetatab(1,j)
          cloc = cos(thetaloc)
          sloc = sin(thetaloc)
          psiloc_mag = psiloc

          call secroot(findsurf,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

          if(error>0) then
             continue
          endif

          rtab(i,1,j) = rloc
          bigrtab(i,1,j) = xmax + rloc*cloc
          ztab(i,1,j) = zmax + rloc*sloc

          rguessl = rloc*0.5d0
          rguessr = rguessr0
!!$			rguessr = rloc*1.05d0

          ! 7/6/2010 poloidal viscosity terms

          bigrloc = bigrtab(i,1,j)
          zloc = ztab(i,1,j)

          rholoc = DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_rho)

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bphi)

          !			print*, rholoc, rho_0
          !			print*, bpolloc, bphiloc, B_0
          viscosity_logs(1,j) = log(rholoc/rho_0) - 1.5d0*log(sqrt(bpolloc**2+bphiloc**2)/B_0)

       enddo

       ! now interpolate the coordinates and minor radius

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            rtab(i,1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            rtab(i,2,1:nsurf+1))

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            bigrtab(i,1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            bigrtab(i,2,1:nsurf+1))

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            ztab(i,1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            ztab(i,2,1:nsurf+1))

       ! 7/6/2010 poloidal viscosity terms

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            viscosity_logs(1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            viscosity_logs(2,1:nsurf+1))

       ! end of poloidal viscosity terms

       ! next, set up the integrand (Jacobian included!)

       ! also include the poloidal viscosity calculation in the same loop
       ! it is actually simpler to immediately save the poloidal viscosity to a tecplot file
       ! (even though it is not elegant at all...)

       allocate (integrand(1:n_int))

       do k = 1, n_int

          thetaloc = t_int(k)

          rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, bigrtab(i,2,1:nsurf+1) )

          zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, ztab(i,2,1:nsurf+1) )

          rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bphi)

          integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
               rloc * bphiloc / (bigrloc * bpolloc)  &
               / (2.d0*pi)

          poloidal_viscosity(i,k) = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
               nsurf+1, viscosity_logs(2,1:nsurf+1) )

          poloidal_viscosity(i,k) = poloidal_viscosity(i,k) / sqrt(1.d0+(rlocprim/rloc)**2)

          write(69,1234) thetaloc, rloc, bigrloc, zloc, poloidal_viscosity(i,k)
		  write(70,1235) thetaloc, rloc, bigrloc, zloc, poloidal_viscosity(i,k)

       enddo

       ! the calculation of q is complete, now let's save some numbers

       call integrate(n_int,integrand,w_int,qval(i))

       deallocate (integrand)

       Rleft(i) = dbsval(pi,ord_surf, thetatab(2,:),  &
            nsurf+1, bigrtab(i,2,1:nsurf+1) )

       Rright(i) = dbsval(0.d0,ord_surf, thetatab(2,:),  &
            nsurf+1, bigrtab(i,2,1:nsurf+1) )

       ! if the trapped particle fraction is required, calculate |B|_MAX

       if(write_all) then

          modB_coef(1,enq+1-i) = psiloc
          B_max =  -brent(1.d-6,pi,2.d0*pi-1.d-6,modBval,tol,thetaloc)*abs(bzero(0.d0))
          ! for convenience, the field is normalized to B_phi_0 in modBval
          modB_coef(2,enq+1-i) = B_max
          continue

          ! also get B_min and details for the bootstrap current calculation

          B_min =  brent(1.d-6,pi,2.d0*pi-1.d-6,modBval_min,tol,thetaloc)*abs(b_phi_zero)

          bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, bigrtab(i,2,1:nsurf+1) )

          zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, ztab(i,2,1:nsurf+1) )

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          boot_Bmin(i) = B_min
          boot_Bp_at_Bmin(i) = bpolloc
          boot_R_of_Bmin(i) = bigrloc

       endif

       if(write_all) then

          call bootstrap_setup

       endif

    enddo

    close(69)
	close(70)

1234 format(5(e16.9, 2x))
1235 format(4(e16.9,','),(e16.9))


    if(write_all) then

       ! sets up the average for bootstrap fraction calculation
       call surf_ave_setup(bpol, n, -1.d0)

    endif

    ! finally, let's calculate the magnetic shear r q'/q

    allocate(Rknot(1:enq+ord_surf))
    allocate(shear(2,1:enq))

    call DBSNAK(enq,Rright,ord_surf,Rknot)

    call DBSINT(enq, Rright, qval, ord_surf,  &
         Rknot, shear(2,:))

    do k = 1,enq

       shear(1,k) = dbsder(1,Rright(k), ord_surf, Rknot,  &
            enq, shear(2,:) )

       rloc = Rright(k) - xmax

       shear(1,k) = shear(1,k) * rloc / qval(k)

    enddo

    ! last, save some tecplot files

    open(4, file='magnetic_R.plt', action='write', status='unknown')

    write(4,*)'TITLE="solution of GS equation with flow, magnetic postprocessing"'
    write(4,*)'Variables =" R [m] ","q", "shear"'
!    write(4,*)'Variables =" R [m] ","q", "shear", "L_ratio"'
    write(4,*)'ZONE I=',2*enq,',F=Point'

    do i = enq,1,-1
       write(4,88) Rleft(i), qval(i), shear(1,i)
!       write(4,88) Rleft(i), qval(i), shear(1,i), L_ratio(i)
    enddo

    do i = 1,enq
       write(4,88) Rright(i), qval(i), shear(1,i)
!       write(4,88) Rright(i), qval(i), shear(1,i), L_ratio(i)
    enddo

    close(4)

    open(4, file='magnetic_psi.plt', action='write', status='unknown')

    write(4,*)'TITLE="solution of GS equation with flow, magnetic postprocessing"'
    write(4,*)'Variables =" psi ","q", "shear"'
!    write(4,*)'Variables =" psi ","q", "shear", "L_ratio"'
    write(4,*)'ZONE I=',enq,',F=Point'

    do i = 1,enq
       write(4,88) psival(i), qval(i), shear(1,i)
!       write(4,88) psival(i), qval(i), shear(1,i), L_ratio(i)
    enddo

    close(4)

    ! insert here the trapped particle calculation setup

    if(write_all) then

       call DBSNAK(enq, modB_coef(1,1:enq), modB_ord, modB_coef(3,:))

       call DBSINT(enq, modB_coef(1,1:enq), modB_coef(2,1:enq), modB_ord,  &
            modB_coef(3,:), modB_coef(4,1:enq))

    endif

    ! insert here the q profile interpolation

    ! first, extrapolate q0 with a parabola

    a_q = (qval(1)-qval(2))/(2.d0*(psival(1)-psival(2))*(psival(1)-psimax))
    b_q = -2.d0*a_q*psimax
    c_q = qval(1) - a_q*psival(1)**2 - b_q*psival(1)

    q_c = a_q*psimax**2 + b_q*psimax + c_q

    ! then, extrapolate q_edge with a linear extrapolation

    qe = qval(enq) + (qval(enq)-qval(enq-1))/(psival(enq)-psival(enq-1))*(-psival(enq)) !since psi_edge=0

	! now set up the interpolation (both ways)

	enqbis = enq+2 ! added first and last point, axis and edge

	ibreak_q = enq+2+q_ord	!enqbis + q_ord
	allocate(q_coef(4,ibreak_q))
!	allocate(psi_of_q_coef(4,1:ibreak_q))

	q_coef(1,1) = 0.d0
	q_coef(2,1) = qe

!	psi_of_q_coef(1,1) = qe
!	psi_of_q_coef(2,1) = 0.d0

	do i=1,enq

		q_coef(1,enqbis-i) = psival(i)
		q_coef(2,enqbis-i) = qval(i)

!		psi_of_q_coef(1,enqbis-i) = qval(i)
!		psi_of_q_coef(2,enqbis-i) = psival(i)

	enddo

	q_coef(1,enqbis) = psimax
	q_coef(2,enqbis) = q_c

!	psi_of_q_coef(1,enqbis) = q_c
!	psi_of_q_coef(2,enqbis) = psimax

	call DBSNAK(enqbis, q_coef(1,1:enqbis), q_ord, q_coef(3,:))

	call DBSINT(enqbis, q_coef(1,1:enqbis), q_coef(2,1:enqbis), q_ord,  &
		q_coef(3,:), q_coef(4,1:enqbis))

!	call DBSNAK(enqbis, psi_of_q_coef(1,1:enqbis), q_ord, psi_of_q_coef(3,:))

!	call DBSINT(enqbis, psi_of_q_coef(1,1:enqbis), psi_of_q_coef(2,1:enqbis), q_ord,  &
!		psi_of_q_coef(3,:), psi_of_q_coef(4,1:enqbis))

	if(write_all) then

	call DBSNAK(enq, modB_coef(1,1:enq), modB_ord, modB_coef(3,:))

	call DBSINT(enq, modB_coef(1,1:enq), modB_coef(2,1:enq), modB_ord,  &
		modB_coef(3,:), modB_coef(4,1:enq))

	endif

	! Sept 25 2013
	! repeat the earlier piece of code with the psi values required for v_normal calculation
	! write output as quantities are calculated (R, Z, theta, v_p, v_psi)

    do i = 1, n_q_vnorm

		if((q_vnorm(i)<qval(1)).or.(q_vnorm(i)>qval(enq))) then
			cycle
		endif

		open(69,file='vpsi_'//trim(vnorm_label(i))//'.plt',action='write')

		write(69,*)'TITLE="normal velocity for solution of GS equation with flow"'
		write(69,345)'Variables ="theta", "r", "R [m] ","z [m]", "q", "psi", "vpol", "vnorm", "v_ratio"'
		write(69,*)'ZONE I=',nsurf+1,',F=Point'

		qloc_mag = q_vnorm(i)

		surf_index = i

		rguessl = rguessl0
		rguessr = rguessr0

		call secroot(findpsi_of_q,psival(enq),psival(1),psiloc,itemax,xtoll,ftoll,error)

		! this cycle determines the (R,Z) coordinates of each magnetic surface
		do j = 1, nsurf+1

			thetaloc = thetatab(1,j)
			cloc = cos(thetaloc)
			sloc = sin(thetaloc)
			psiloc_mag = psiloc

			call secroot(findsurf,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

			if(error>0) then
				continue
			endif

			bigrloc = xmax + rloc*cloc
			zloc = zmax + rloc*sloc

			rguessl = rloc*0.5d0
			rguessr = rguessr0

			bxloc = DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
									zknot_mag,n,n,bscoef_bx)

			bzloc = DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
									zknot_mag,n,n,bscoef_bz)

			vxloc = DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
									zknot_mag,n,n,bscoef_vx)

			vzloc = DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
									zknot_mag,n,n,bscoef_vz)

			vpolloc = sqrt(vxloc**2+vzloc**2)
			vpsiloc = (vxloc*bzloc-vzloc*bxloc)/sqrt(bxloc**2+bzloc**2)

			if(vpolloc>0.d0) then
				v_ratio_loc = vpsiloc/vpolloc
			else
				v_ratio_loc = 1.d9
			endif

			write(69,2345) thetaloc,rloc,bigrloc,zloc,qloc_mag,psiloc,vpolloc,vpsiloc,v_ratio_loc

       enddo


	enddo

	close(69)

2345 format(9(e16.9, 2x))
345 format(A100)

    ! now let's clean up all the data

    deallocate(Rknot)
    deallocate(shear)

    deallocate(w_int)
    deallocate(t_int)

    if(write_all) then

       continue

    else

       deallocate(psival)
       deallocate(qval)
       deallocate(w_int_boot)
       deallocate(t_int_boot)

    endif

    !	deallocate(bscoef_psi)
    !	deallocate(bscoef_bpol)
    !	deallocate(bscoef_bphi)

    continue

66  format(3(e26.17,3x))
88  format(4(e26.17,3x))
    !________________________

  contains

    subroutine bootstrap_setup

      implicit none

      real(kind=dkind), dimension(5,n_int) :: integrand_boot
      real(kind=dkind) :: lambda, btotloc, int_step
      real(kind=dkind) :: Rmin_loc, Rmax_loc

      integer :: it, jt

      ! step 0, compute the length of each magnetic surface (one at a time!)

      Rmin_loc = rmajor
      Rmax_loc = 0.d0

      do it = 1, n_int

         thetaloc = t_int(it)

         rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
              nsurf+1, rtab(i,2,1:nsurf+1) )

         rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
              nsurf+1, rtab(i,2,1:nsurf+1) )

         integrand_boot(1,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
              rloc / (2.d0*pi)

      enddo

      call integrate(n_int,integrand_boot(1,:),w_int,surf_length(i))

      ! first, compute the average in the denominator and the effective trapped particle fraction integral;
      ! then do directly the d_lambda integral

      do jt = 1, n_int

         lambda = t_int_boot(jt)

         do it = 1, n_int

            thetaloc = t_int(it)

            rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
                 nsurf+1, rtab(i,2,1:nsurf+1) )

            rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
                 nsurf+1, rtab(i,2,1:nsurf+1) )

            bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
                 nsurf+1, bigrtab(i,2,1:nsurf+1) )

            zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
                 nsurf+1, ztab(i,2,1:nsurf+1) )

            bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
                 zknot_mag,n,n,bscoef_bpol)

            bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
                 zknot_mag,n,n,bscoef_bphi)

            btotloc = min(sqrt(bpolloc**2+bphiloc**2),B_max)

            ! this integrates the denominator
            integrand_boot(1,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
                 sqrt(1.d0-lambda*btotloc/B_max) *  &
                 rloc / (2.d0*pi)

            ! this integrates (B/B_max)**2
            integrand_boot(2,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
                 (btotloc/B_max)**2 *  &
                 rloc / (2.d0*pi)

            if(bootstrap_option>=1) then
               ! also calculate other averages

               ! this integrates B**2
               integrand_boot(4,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
                    btotloc**2 *  &
                    rloc / (2.d0*pi)

               ! this integrates 1/B**2
               integrand_boot(5,it) = sqrt(1.d0+(rlocprim/rloc)**2) /  &
                    btotloc**2 *  &
                    rloc / (2.d0*pi)

            endif

            ! also check for the aspect ratio at this stage

            if(jt==1) then
               ! just do it once!

               Rmax_loc = max(Rmax_loc,bigrloc)

               Rmin_loc = min(bigrloc,Rmin_loc)

            endif

         enddo

         call integrate(n_int,integrand_boot(1,:),w_int,int_step)
         ! denominator

         call integrate(n_int,integrand_boot(2,:),w_int,B2_hat_ave(i))
         ! B**2

         integrand_boot(3,jt) = lambda / int_step

         if(bootstrap_option>=1) then

            call integrate(n_int,integrand_boot(4,:),w_int,B2_ave(i))
            ! B**2

            call integrate(n_int,integrand_boot(5,:),w_int,Bm2_ave(i))
            ! B**2

         endif

      enddo

      call integrate(n_int,integrand_boot(3,:),w_int_boot,int_step)

      eff_trap(i) = 1.d0 - 0.75d0 * B2_hat_ave(i)* int_step

      if(eff_trap(i)<0.d0) eff_trap(i) = 0.d0
      if(eff_trap(i)>1.d0) eff_trap(i) = 1.d0

      ! set the "aspect ratio"
      inv_asp_ratio(i) = (Rmax_loc-Rmin_loc)/(Rmax_loc+Rmin_loc)

      ! set the "center"
      boot_rcenter(i) = (Rmax_loc+Rmin_loc)/2.d0

      ! the minimum field, the poloidal field at the same location and
      ! also the location are saved in the magnetic output section

      ! then proceed with the various physical quantities (density, temperature, etc.)
      ! the case with more complicated dependencies will be differentiated later

      if(static_equi) then
         ! all depends on psi only

         boot_ne(i) = dofpsi(psival(i))/mass ! in m^-3
         ! ion mass to turn the plasma mass in a number density
         boot_pe(i) = pofpsi(psival(i))*pe_ov_p
         ! this is p_e = pe_ov_p * (total p)

         boot_Te(i) = boot_pe(i)/boot_ne(i)/eV ! in eV
         ! this is T_e = pe_ov_p * (total T)

         boot_ni(i) = dofpsi(psival(i))/mass ! in m^-3
         ! ion mass to turn the plasma mass in a number density
         boot_pi(i) = pofpsi(psival(i))*(1.d0-pe_ov_p)
         ! this is p_i = (total p) - p_e

         boot_Ti(i) = boot_pi(i)/boot_ni(i)/eV ! in eV
         ! this is p_i = (total p) - p_e

         boot_neprim(i) = dddpsi(psival(i))/mass	! * psic?
         boot_peprim(i) = dpdpsi(psival(i)) * pe_ov_p	! * psic?

         boot_Teprim(i) = ( boot_peprim(i)/boot_ne(i) -  &
              boot_neprim(i)*boot_pe(i)/boot_ne(i)**2 )/eV

         !		boot_neprim(i) = boot_neprim(i)/mass

         boot_niprim(i) = dddpsi(psival(i))/mass	! * psic?
         boot_piprim(i) = dpdpsi(psival(i)) * (1.d0-pe_ov_p)	! * psic?

         boot_Tiprim(i) = ( boot_piprim(i)/boot_ni(i) -  &
              boot_niprim(i)*boot_pi(i)/boot_ni(i)**2 )/eV

         !		boot_niprim(i) = boot_niprim(i)/mass

      endif

      continue

      return

    end subroutine bootstrap_setup




  end subroutine magnetic_output



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine plasma_boundary(psi,xmax,zmax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! for free-boundary equilibria, determine the plasma edge (psi=0)

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBSNAK, DBSINT, DBSDER


	real(kind=dkind), dimension(1:n,1:n) :: psi

	real(kind=dkind) :: psimax, rguessl, rguessr, rguessl0, rguessr0
	real(kind=dkind) :: psiloc, thetaloc, rloc
	integer :: i, j,k
	integer :: nsurfdim
	real(kind=dkind) :: xmax, zmax

	integer :: itemax = 200
	real(kind=dkind) :: xtoll = 1.d-6
	real(kind=dkind) :: ftoll = 1.d-6
	integer :: error

	! first, allocate the arrays that would have been automatic,
	! if not for the complaints of the LINUX compiler

	allocate(bscoef_psi(1:n,1:n))

	allocate(xknot_mag(1:n+ord_loc))
	allocate(zknot_mag(1:n+ord_loc))

	xmax_mag = xmax
	zmax_mag = zmax

	enq = (n-1)*3/8

	!first set up the interpolations

	call DBSNAK(n,x_coord(1:n),ord_loc,xknot_mag)
	call DBSNAK(n,z_coord(1:n),ord_loc,zknot_mag)

	call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),psi,n,  &
					ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_psi)


	psimax = DBS2VL(xmax,zmax,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_psi)

	! get the magnetic surface corresponding to the plasma edge

	nsurfdim = nsurf+1+ord_surf

	allocate(thetatab(2,nsurfdim))

	! this sets up the 1D grid and the integration routine

	do j = 1, nsurf+1
		thetatab(1,j) = (j-1.d0)*pi*2.d0/nsurf
	enddo

	call DBSNAK(nsurf+1, thetatab(1,1:nsurf+1),  &
		ord_surf,thetatab(2,1:nsurfdim))


	rguessr0 = max(x_coord(n)-xmax, xmax-x_coord(1), z_coord(n)-zmax, zmax-z_coord(1))
	rguessl0 = (x_coord(n)-x_coord(1)+z_coord(n)-z_coord(1))/5.d2

	open(33,file='r_edge_of_theta.dat')

	rguessl = rguessl0
	rguessr = rguessr0

	psiloc = 0.d0

	! this cycle determines the (R,Z) coordinates of the magnetic surface
	! note that the radius is calculated from the geometric axis (differently from the magnetic_output routine)
	do j = 1, nsurf+1

		thetaloc = thetatab(1,j)
		cloc = cos(thetaloc)
		sloc = sin(thetaloc)
		psiloc_mag = psiloc

		call secroot(findsurf2,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

		if(error>0) then
			continue
		endif

		write(33,*) thetaloc, rloc

		rguessl = rloc*0.5d0
		rguessr = rguessr0

	enddo

	close(33)

	deallocate(bscoef_psi)
	deallocate(xknot_mag)
	deallocate(zknot_mag)
	deallocate(thetatab)



end subroutine plasma_boundary




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function findsurf(r) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2VL

	real(kind=dkind) :: r
	real(kind=dkind) :: answer
	real(kind=dkind) :: bigRtemp, bigZtemp

	bigRtemp = xmax_mag + r*cloc
	bigZtemp = zmax_mag + r*sloc

	if(bigRtemp>x_coord(n)) bigRtemp = x_coord(n)
	if(bigRtemp<x_coord(1)) bigRtemp = x_coord(1)

	if(bigZtemp>z_coord(n)) bigZtemp = z_coord(n)
	if(bigZtemp<z_coord(1)) bigZtemp = z_coord(1)

	answer = DBS2VL(bigRtemp,bigZtemp, &
								ord_loc,ord_loc,xknot_mag,zknot_mag,n,n,bscoef_psi) &
										- psiloc_mag

	return

end function findsurf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function findsurf2(r) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2VL

	real(kind=dkind) :: r
	real(kind=dkind) :: answer
	real(kind=dkind) :: bigRtemp, bigZtemp

	bigRtemp = rmajor + r*cloc
	bigZtemp = r*sloc

	if(bigRtemp>x_coord(n)) bigRtemp = x_coord(n)
	if(bigRtemp<x_coord(1)) bigRtemp = x_coord(1)

	if(bigZtemp>z_coord(n)) bigZtemp = z_coord(n)
	if(bigZtemp<z_coord(1)) bigZtemp = z_coord(1)

	answer = DBS2VL(bigRtemp,bigZtemp, &
								ord_loc,ord_loc,xknot_mag,zknot_mag,n,n,bscoef_psi) &
										- psiloc_mag

	return

end function findsurf2



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function modBval(th) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			use pseudo_IMSL, only : DBS2VL

			real(kind=dkind) :: th
			real(kind=dkind) :: answer, Bp
			real(kind=dkind) :: bigrloc, zloc, rlocprim, bpolloc, bphiloc


			bigrloc = dbsval(th,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(surf_index,2,1:nsurf+1) )

			zloc = dbsval(th,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(surf_index,2,1:nsurf+1) )

			rlocprim = dbsder(1,th, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(surf_index,2,1:nsurf+1) )

			bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bpol)

			bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bphi)

			answer = -sqrt((bpolloc**2+bphiloc**2)/bzero(0.d0)**2)
			! define this way to "always" have an answer of order 1
			! and negative to exploit minimum search

			continue

			return

	end function modBval

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	function modBval_min(th) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			real(kind=dkind) :: th
			real(kind=dkind) :: answer

			answer = -modBval(th)

			continue

			return

	end function modBval_min

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function findpsi_of_q(psiloc) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : dbsval_safe

	real(kind=dkind) :: psiloc
	real(kind=dkind) :: answer

	answer  = dbsval_safe(psiloc,q_ord, q_coef(3,:),  &
					ibreak_q-q_ord, q_coef(4,1:ibreak_q-q_ord), q_coef(2,1:ibreak_q-q_ord) )  &
					- qloc_mag

	return

end function findpsi_of_q

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine bootstrap_setup_ave(psi,rho, p, Br, Bz, T, b_phi, J_par,  &
					base, nn, base_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: psi,rho, p, Br, Bz, T, b_phi,  &
							J_par, base

	real(kind=dkind), dimension(1:nn,1:nn) :: thing
	real(kind=dkind), dimension(1:enq) :: boot_ptot, boot_ptotprim

	real(kind=dkind) :: base_option

	integer :: i,j



	! new average routine
	! this routine computes the averages of density (same for ions and electrons)
	! and of total pressure

	call get_new_avs(boot_ptot, boot_ne, base, p(1:nn,1:nn), rho(1:nn,1:nn), nn, base_option,  &
										boot_ptotprim, boot_neprim)

	if(bootstrap_option>=1) then

		call get_toroidal_flux(psi(1:nn,1:nn),b_phi(1:nn,1:nn))
		call averages_NCLASS(nn,base, base_option)

	endif

	! arrange the results for ions and elctrons
	! start from the density

	do i = 1,enq

		boot_ni(i) = boot_ne(i)
		boot_niprim(i) = boot_neprim(i)

	enddo

	! then take care of the pressures

	do i = 1,enq

		boot_pe(i) = boot_ptot(i)*pe_ov_p
		boot_pi(i) = boot_ptot(i)-boot_pe(i)

		boot_peprim(i) = boot_ptotprim(i)*pe_ov_p
		boot_piprim(i) = boot_ptotprim(i)-boot_peprim(i)

	enddo

	! finally, take care of the temperatures (no further interpolations needed)

	do i = 1,enq

		boot_Te(i) = boot_pe(i)/boot_ne(i)*mass/eV ! in eV
		boot_Ti(i) = boot_pi(i)/boot_ni(i)*mass/eV ! in eV

		boot_Teprim(i) = ( boot_peprim(i)/boot_ne(i) -  &
								boot_neprim(i)*boot_pe(i)/boot_ne(i)**2 )*mass/eV

		boot_Tiprim(i) = ( boot_piprim(i)/boot_ni(i) -  &
								boot_niprim(i)*boot_pi(i)/boot_ni(i)**2 )*mass/eV

	enddo

	! last thing, change the units of the densities

	do i = 1,enq

		boot_ne(i) = boot_ne(i)/mass ! now in m^-3
		boot_ni(i) = boot_ni(i)/mass ! now in m^-3

		boot_neprim(i) = boot_neprim(i)/mass ! now in m^-3
		boot_niprim(i) = boot_niprim(i)/mass ! now in m^-3

	enddo

	! now all the quantities for the bootstrap current calculation are available
	! still compute <J_par B>, to have the bootstrap fraction

	do j = 1,nn
	do i = 1,nn

		thing(i,j) = J_par(i,j) * sqrt( Br(i,j)**2 + Bz(i,j)**2 + b_phi(i,j)**2 )

	enddo
	enddo

	call get_thing_ave(JparB_ave, base, thing, nn, base_option)

	continue

	return

	end subroutine bootstrap_setup_ave

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_toroidal_flux(psi,b_phi)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind), dimension(1:n,1:n), intent(in) :: psi, b_phi
	real(kind=dkind) :: dA, bigR
	real(kind=dkind) :: bphiloc
	integer :: i,j,k
	real(kind=dkind), dimension(2,0:enq+1) :: temp_array !just to avoid problems with interpolation


	dA = dx*dz

	boot_tor_flux = 0.d0

	do j=1,n
		do i=1,n

			if(sort_grid(i,j)>=1) then

				do k=1,enq+1

					if(grid_type/=0) dA = 0.25d0*(dx_a(i) + dx_a(i-1)) *  &
											(dz_a(j) + dz_a(j-1))

					bphiloc = b_phi(i,j)

					if(psi(i,j)>psival(k)) then	!inner region for kth surface

						bigR = x_coord(i) + dx_a(i) -  dx_a(i-1)
						boot_tor_flux(k) = boot_tor_flux(k) + bphiloc * dA

					endif

				enddo

			endif

		enddo
	enddo

	do k = 0,enq+1

		boot_tor_rho(k) = sqrt(boot_tor_flux(k)/boot_tor_flux(enq+1)) * a_elps
		!WARNING: a_elps = (R_max-R_min)/2 is defined elsewhere (readinput)

	enddo

	do k = 0,enq+1

		temp_array(1,enq+1-k) = psival(k)
		temp_array(2,enq+1-k) = boot_tor_rho(k)

	enddo


	allocate(bscoef_tor_rho(2,enq+2+q_ord))

	! interpolation setup has been moved to a separate function
	call interp_setup(enq+2, q_ord, &
		temp_array(1,0:enq+1),temp_array(2,0:enq+1), &
		bscoef_tor_rho(1,1:enq+2+q_ord),bscoef_tor_rho(2,1:enq+2))

	continue

	return

end subroutine get_toroidal_flux

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine averages_NCLASS(nn,base, base_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: base
	real(kind=dkind) :: base_option
	real(kind=dkind) :: cutoff = 1.d-12

	real(kind=dkind), dimension(:,:), allocatable :: integrand
	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, bigrlocprim, zlocprim,  &
								baseloc, bphiloc, der_val, bpolloc

	integer :: i, j, k


	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(base, nn, base_option)

	endif


	! then take care of the derivative part
	! the integration procedure is kept here to avoid additional array copying

	allocate (integrand(1:2,1:n_ave_int))
	! "1" is Fhat, "2" is grad rho/B^2

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_boot_base)

			bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bphi)

			bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bpol)

			!------------------normal and derivative------------------

			der_val = dbsder(1,psival(i), q_ord, bscoef_tor_rho(1,1:enq+2+q_ord),  &
					enq+2, bscoef_tor_rho(2,1:enq+2) )
			! note: this is (d rho / d psi), which is equal to 1/(d psi / d rho)

!			if(abs(der_val)<cutoff) der_val = sign(cutoff, der_val)

			integrand(1,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * der_val * baseloc**base_option  &
									*bphiloc / (2.d0*pi) !2pi? no  mu_mag!

			integrand(2,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * (der_val*bpolloc*bigrloc)**2/  &
									(bphiloc**2+bpolloc**2) *  &
									baseloc**base_option  &
									/ (2.d0*pi)


		enddo

		! then calculate the integrals

		call integrate(n_ave_int, integrand(1,:), w_ave_int, boot_fhat(i))

		boot_fhat(i) = boot_fhat(i) / surf_ave_base(i) * rmajor ! * pi !WARNING!!!

		call integrate(n_ave_int, integrand(2,:), w_ave_int, boot_grho_ov_B2(i))

		boot_grho_ov_B2(i) = boot_grho_ov_B2(i) / surf_ave_base(i)


		continue

	enddo



	deallocate(integrand)

	continue


	return


end subroutine averages_NCLASS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_bootstrap_NCLASS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: ierror
	real(kind=dkind) :: grad_T(2), Tloc(2), densloc(2,1), pprimloc(2,1)
	real(kind=dkind) :: den_cutoff
	real(kind=dkind), dimension(2) :: amu !atomic mass units
	real(kind=dkind) :: der_cutoff = 1.d-8
	real(kind=dkind) :: psiloc, der_val, der_sec_val, eps_loc
	real(kind=dkind) :: gradphi, gradphi_sq, EdotB, ngrth
	real(kind=dkind) :: PS_mom(3)
	integer :: i, k

	real(kind=dkind) :: external_force_fex_iz(3,mx_mi,mx_mz) = 0.d0
	integer, dimension(1:2) :: isotope_status = 0
	integer, dimension(1:2) :: charge_status
	integer :: k_potato = 0
	integer :: m_s = 2

	! see NCLASS documentation for the meaning of the following
	real(kind=dkind) :: p_etap, p_exjb
	real(kind=dkind) :: calm_i(3,3,mx_mi)
	real(kind=dkind) :: caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),  &
					capn_ii(3,3,mx_mi,mx_mi)

	real(kind=dkind) :: boot_p_coeff(mx_ms), boot_T_coeff(mx_ms)

	real(kind=dkind) :: dn_s(mx_ms), gfl_s(5,mx_ms), qfl_s(5,mx_ms), sqz_s(mx_ms),  &
								upar_s(3,3,mx_ms), utheta_s(3,3,mx_ms), vn_s(mx_ms), veb_s(mx_ms), &
								qeb_s(mx_ms), xi_s(mx_ms), ymu_s(3,3,mx_ms)
	real(kind=dkind) :: chip_ss(mx_ms,mx_ms), chit_ss(mx_ms,mx_ms),  &
								dp_ss(mx_ms,mx_ms), dt_ss(mx_ms,mx_ms)



	do i = 1,enq

		psiloc = psival(i)
		eps_loc = inv_asp_ratio(i)

		!------------ set geometry moments for PS current--------------
		! I have no clue of what this is (yet)

		do k=1,3
			PS_mom(k)=0.d0
		enddo

		if(p_eps.gt.0.0) then

			do k=1,3
				PS_mom(k)=k*((1.d0-sqrt(1.d0-eps_loc**2))/eps_loc)**(2.d0*k)  &
					*(1.d0 + k * sqrt(1.d0-eps_loc**2))/((1.d0-eps_loc**2)**1.5d0  &
					*(qval(i)*rmajor)**2)
			enddo

		!---------------end PS current stuff----------------

		endif

		amu(1) = 2.d0
		amu(2) = amu(1)*me_ov_mi !?

		charge_status(1) = 1
		charge_status(2) = -1

		den_cutoff = dofpsi(0.d0) / mass

		EdotB = 0.d0

		der_val = dbsder(1,psiloc, q_ord, bscoef_tor_rho(1,1:enq+2+q_ord),  &
						enq+2, bscoef_tor_rho(2,1:enq+2) )
				! note: this is (d rho / d psi), which is equal to 1/(d psi / d rho)

		if(abs(der_val)<der_cutoff) der_val = sign(der_cutoff, der_val)

		der_sec_val = dbsder(2,psiloc, q_ord, bscoef_tor_rho(1,1:enq+2+q_ord),  &
						enq+2, bscoef_tor_rho(2,1:enq+2) )
				! note: this is (d2 rho / d psi2), which is NOT equal to 1/(d2 psi / d rho2)

		if(abs(der_sec_val)<der_cutoff) der_sec_val = sign(der_cutoff, der_sec_val)


		gradphi = omegaofpsi(psiloc) / der_val

		gradphi_sq = ( domegadpsi(psiloc) - omegaofpsi(psiloc) * der_sec_val/der_val)  &
							/ der_val**2


		Tloc(1) = boot_Ti(i) / 1.d3 ! this is in keV
		Tloc(2) = boot_Te(i) / 1.d3 ! this is in keV

		grad_T(1) = boot_Tiprim(i) / der_val / 1.d3 ! d T / d rho
		grad_T(2) = boot_Teprim(i) / der_val / 1.d3 ! d T / d rho

		densloc(1,1) = boot_ni(i)
		densloc(2,1) = boot_ne(i)

		pprimloc(1,1) = boot_piprim(i) * mass/eV/1.d3 / der_val
		pprimloc(2,1) = boot_peprim(i) * mass/eV/1.d3 / der_val
		pprimloc(1,1) = boot_piprim(i) /eV/1.d3 / der_val	!?
		pprimloc(2,1) = boot_peprim(i) /eV/1.d3 / der_val	!?
		! in whatever ridiculous units these are supposed to be

		ngrth = 1.d0/qval(i)/rmajor	!what the heck is this?

		call NCLASS(2, k_potato, 2, 1, den_cutoff, 0.d0, 0.d0,  &
							B2_ave(i), Bm2_ave(i), EdotB, boot_fhat(i), PS_mom, eff_trap(i), boot_grho_ov_B2(i),  &
							gradphi, gradphi_sq, ngrth, amu(1:2), grad_T, Tloc,  &
							densloc, external_force_fex_iz, pprimloc, m_s, isotope_status, charge_status, J_boot(6,i),  &
							p_etap,p_exjb,calm_i,caln_ii,capm_ii,capn_ii,  &
							boot_p_coeff, boot_T_coeff ,dn_s,gfl_s,qfl_s,sqz_s,upar_s,  &
							utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,  &
							chit_ss,dp_ss,dt_ss, ierror)

		J_boot(6,i) = -J_boot(6,i) ! derivatives have the wrong sign!
		J_boot(7,i) = J_boot(6,i) / JparB_ave(i)

		el_resistivity(3,i) = p_etap ! p_etap is the resistivity

	enddo

! m_i = 2: ions and electrons?
! start without potatoes
! put neutral and charged D, 0 density for the neutrals
! <E dot B> = 0 (ideal MHD)

	continue

end subroutine get_bootstrap_NCLASS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_bootstrap2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: i
	real(kind=dkind) :: a2_el, a2_ion, a1
	real(kind=dkind) :: Coul_log_e, Coul_log_i ! Coulomb logarithm log(Lambda)
	real(kind=dkind) :: taue, vloc
	real(kind=dkind) :: p_nondim, pprim_nd, e_mass
	real(kind=dkind) :: charge_cgs = 4.8032d-10 ! to be fixed
	real(kind=dkind) :: nustar, nustarboh, nuee
	real(kind=dkind) :: ft, ft2, ft3
	real(kind=dkind) :: Z_eff = 1.d0 ! this would be Z_effective
	real(kind=dkind), dimension(:,:), allocatable :: nustar_array
	real(kind=dkind) :: sigma_neo, sigma_spitz
	real(kind=dkind) :: Z_i, n_i
	real(kind=dkind) :: el31, el32, el31_0, el32_0, boot_alpha, boot_alpha_0, el34


	e_mass = mass * me_ov_mi


	allocate(nustar_array(2,enq))

	do i = 1, enq

		Z_i = Z_eff !is this so?

!		Z_i = 1.d0+1.5d0*(abs(psival(i)/psic))**0.25d0

		Z_i = Z_eff_0 + delta_Z_eff*(abs(psival(i)/psic))**Z_eff_exp

		ft = eff_trap(i)
!		ft = eff_trap(i)/(1.d0+sqrt(nustar*Z_eff) + 0.25d0*nustar/Z_eff)

         call sigmaneo(sigma_neo,sigma_spitz,nustar,ft,boot_ne(i),boot_Te(i),Z_i, &
               abs(qval(i)),boot_Rcenter(i),inv_asp_ratio(i))

		el_resistivity(1,i) = 1.d0/sigma_spitz
		el_resistivity(2,i) = 1.d0/sigma_neo

		nustar_array(1,i) = nustar

        nustar = nustar/Z_i

		nustar_array(2,i) = nustar

! WARNING! NOTABLE ASSUMPTION!!
! assumes z_imp=6
         n_i=(6.d0-Z_i)/(6.d0-1.d0) * boot_ne(i)
         call bootstrap_coeffs(ft,abs(qval(i)),boot_rcenter(i),inv_asp_ratio(i),boot_Te(i),boot_ne(i),boot_Ti(i), &
               n_i,Z_i,1.d0,el31_0,el32_0,boot_alpha_0,el31,el32,el34,boot_alpha)


!     For RJBSOS(K,4), use the standard definition. However one mixes p' from chease inputs and ne', Te', Ti' from
!     experimental data. If ne*Te+ni*Ti exper. is not near p_chease, it can give strange profiles, as some cancellation do not appear.
!     Thus develop p'=ne' Te + ne Te' + ni' Ti + ni Ti', and assume ni'/ni = ne'/ne and take ni*Ti=p-pe:
!     => for (K,1-3) use developed formula as well as for nue*=0 case
!!$         RJBSOS(K,1) = - TMF(K)* &
!!$     &     (ZL31_0*(ZA1+ZALFA_0*(1.d0/ZRPEOP-1.d0)*ZA2I) + ZL32_0 * ZA2E)

		! zero nustar case, in physical units
         J_boot(1,i) = rmajor*bzero(psival(i)) * (boot_pi(i)+boot_pe(i)) * &
           (el31_0*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31_0+el32_0)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(1.d0+boot_alpha)*boot_Tiprim(i)/boot_Ti(i))

		! finite nustar case, in physical units
         J_boot(2,i) = rmajor*bzero(psival(i)) * (boot_pi(i)+boot_pe(i)) * &
           (el31*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31+el32)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(el31+el34*boot_alpha)*boot_Tiprim(i)/boot_Ti(i))

		! bootstrap fraction = bootstrap / <J_par B>
		J_boot(3,i) = J_boot(2,i) / JparB_ave(i)

		! zero nustar case CHEASE
		J_boot(4,i) = rmajor*bzero(psival(i)) * boot_ne(i)*boot_Te(i)/pe_ov_p* eV*mu_mag/abs(b_phi_zero)**3 * & !1.602d-19*4.d-07*CPI/B0EXP**2 * &
           (el31_0*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31_0+el32_0)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(1.d0+boot_alpha_0)*el31_0*boot_Tiprim(i)/boot_Ti(i))

		! finite nustar case CHEASE
		J_boot(5,i) = rmajor*bzero(psival(i)) * boot_ne(i)*boot_Te(i)/pe_ov_p* eV*mu_mag/abs(b_phi_zero)**3 *  & !1.602d-19*4.d-07*CPI/B0EXP**2 * &
           (el31*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31+el32)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(el31+el34*boot_alpha)*boot_Tiprim(i)/boot_Ti(i))	!ratio 1.1901e-005?


		continue

	enddo


	deallocate(nustar_array)

	return

	end subroutine get_bootstrap2

!-----------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bootstrap_coeffs(pft,pq,pR,peps,pte,pne,pti,pni,pzeff,pzion, &
      pl31_0,pl32_0,palfa_0,pl31,pl32,pl34,palfa)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! IMPORTANT: ADAPTED FROM THE CHEASE CODE (THANKS TO O. SAUTER)

!
!     WARNING: in MKSA
!
!     Compute Bootstrap coefficients using formulas from O. Sauter et al, Phys. Plasmas 7 (1999) 2834.
!
!     Assumes to compute on a single flux surface with:
! Inputs:
!     pft   : trapped fraction
!     pq    : safety factor
!     pR    : Geometrical center of given flux surface in [m]
!     peps  : Inverse aspect ratio of given flux surface
!     pte   : Electron temperature [eV]
!     pne   : Electron density [1/m**3]
!     pti   : Ion temperature [eV]
!     pni   : Main ion density [1/m**3]
!     pzeff : Effective charge (used for Z in electronic terms)
!     pzion : Main ion charge
! Outputs:
!     pl31_0  : L31 coefficient assuming nuestar=0
!     pl32_0  : L32 coefficient assuming nuestar=0
!     palfa_0 : Alfa coefficient assuming nuestar=0
!     pl31    : L31 coefficient
!     pl32    : L32 coefficient
!     pl34    : L34 coefficient (L34 for nuestar=0 is identical to L31_0)
!     palfa   : Alfa coefficient
!

  implicit none

  real(kind=dkind), intent(in) :: pft,pq,pR,peps,pte,pne,pti,pni,pzeff,pzion
  real(kind=dkind), intent(out) :: pl31_0,pl32_0,palfa_0,pl31,pl32,pl34,palfa

  real(kind=dkind)  :: znuestar, znuistar, zlnlam_e, zlnlam_i, zdummy
!-----------------------------------------------------------------------

!     basic parameters

  zlnlam_e = 17.d0
  IF ((pne.gt.0.d0).and.(pte.gt.0.d0)) then
     zlnlam_e = 31.3d0 - log(sqrt(pne)/pte)
     zlnlam_e = max(10.d0,zlnlam_e)
  ENDIF
  zlnlam_i = 17.d0
  IF ((pni.gt.0.d0).and.(pti.gt.0.d0)) then
     zlnlam_i = 30.d0 - log(pzion**3.d0*sqrt(pni)/abs(pti)**1.5d0)
     zlnlam_i = max(10.d0,zlnlam_i)
  ENDIF
  znuestar = 6.921d-18 * pq*pR*pne*pzeff*zlnlam_e / (pte*pte*peps**1.5d0)
  znuistar = 4.900d-18 * pq*pR*pni*pzion**4*zlnlam_i / (pti*pti*peps**1.5d0)

!     Compute coefficients for nustar=0

  call final_bootstrap_coeffs(pl31_0,pl32_0,zdummy,palfa_0,pft,pzeff)

!     finite nustar
  call final_bootstrap_coeffs(pl31,pl32,pl34,palfa,pft,pzeff,znuestar,znuistar)

	continue

  return

  end subroutine bootstrap_coeffs

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine final_bootstrap_coeffs(L31, L32, L34, ALFA, ft, Zeff, nuestar, nuistar)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! IMPORTANT: ADAPTED FROM THE CHEASE CODE (THANKS TO O. SAUTER)

    real(kind=dkind), intent(in)  :: ft
    real(kind=dkind), OPTIONAL, intent(in)  :: Zeff, nuestar, nuistar
    real(kind=dkind), intent(out) :: L31, L32, L34, ALFA

    real(kind=dkind)  :: ZZ, znuestar, znuistar, zsqnuest, zeffp1, zsqnui, znui2ft6
    real(kind=dkind)  :: zft31eff, zft32ee_eff, zft32ei_eff, zft34eff, zalfa0
!-----------------------------------------------------------------------
!
    ZZ = 2.d0

    if ( PRESENT(zeff) ) ZZ = zeff

    znuestar = 0.d0

    if ( PRESENT(nuestar) ) znuestar = nuestar

    znuistar = 0.d0

    if ( PRESENT(nuistar) ) znuistar = nuistar

    zsqnuest = sqrt(znuestar)

!  effective trapped fractions

    zsqnuest = sqrt(znuestar)

    zft31eff = ft / (1d0+(1.d0-0.1d0*ft)*zsqnuest &
           + 0.5d0*(1.d0-ft)*znuestar/ZZ)

    zft32ee_eff = ft / (1.d0 + 0.26d0*(1.d0-ft)*zsqnuest &
           + 0.18d0*(1.d0-0.37d0*ft)*znuestar/sqrt(ZZ))

    zft32ei_eff = ft / (1.d0 + (1.d0+0.6d0*ft)*zsqnuest &
           + 0.85d0*(1.d0-0.37d0*ft)*znuestar*(1.d0+ZZ))

    zft34eff = ft / (1.d0+(1.d0-0.1d0*ft)*zsqnuest &
           + 0.5d0*(1.d0-0.5d0*ft)*znuestar/ZZ)

    zalfa0 = - 1.17d0*(1.d0-ft) / (1.d0-0.22d0*ft-0.19d0*ft**2)

!coefficients

    zeffp1 = ZZ+1.d0

    L31 = zft31eff * ( (1.d0+1.4d0/zeffp1) &
          - zft31eff* (1.9d0/zeffp1 - zft31eff * (0.3d0/zeffp1 + 0.2d0/zeffp1 * zft31eff)))

    L32 = (0.05d0+0.62d0*ZZ)/ZZ/(1.d0+0.44d0*ZZ)*(zft32ee_eff-zft32ee_eff**4) &
          +  zft32ee_eff**2*(1.d0-1.2d0*zft32ee_eff+0.2d0*zft32ee_eff**2) &
                           /(1.d0+0.22d0*ZZ) &
          - (0.56d0+1.93d0*ZZ)/ZZ/(1.d0+0.44d0*ZZ)*(zft32ei_eff-zft32ei_eff**4) &
          +  zft32ei_eff**2*(1.d0-0.55d0*zft32ei_eff-0.45d0*zft32ei_eff**2) &
                           * 4.95d0/(1.d0+2.48d0*ZZ) &
          + 1.2d0 / (1.d0+0.5d0*ZZ) * (zft32ee_eff**4-zft32ei_eff**4)

    L34 = zft34eff * ( (1.d0+1.4d0/zeffp1) &
          - zft34eff* (1.9d0/zeffp1 - zft34eff * (0.3d0/zeffp1 + 0.2d0/zeffp1 * zft34eff)))

    zsqnui = sqrt(znuistar)
    znui2ft6 = znuistar**2 * ft**6

    ALFA = ((zalfa0 + 0.25d0*(1.d0-ft**2)*zsqnui) &
                             / (1.d0+0.5d0*zsqnui) + 0.315d0*znui2ft6) &
          / (1.d0 + 0.15d0*znui2ft6)

	continue

    return

  end subroutine final_bootstrap_coeffs



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine sigmaneo(signeo,sigsptz,nuestar,ft,ne,te,zeff,q,R,eps)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! IMPORTANT: ADAPTED FROM THE CHEASE CODE (THANKS TO O. SAUTER)

    ! all are scalar input variables
    ! te in [eV], ne in [m**-3], R in [m]

	implicit none

    real(kind=dkind), intent(out) :: signeo, sigsptz, nuestar
    real(kind=dkind), intent(in)  :: ft, ne, te
    real(kind=dkind), OPTIONAL, intent(in) :: zeff
    real(kind=dkind), OPTIONAL, intent(in) :: q
    real(kind=dkind), OPTIONAL, intent(in) :: R
    real(kind=dkind), OPTIONAL, intent(in) :: eps

    real(kind=dkind)  :: z_zeff, zNZ, zlnL, zft33eff

    z_zeff = 2.d0
    IF ( PRESENT(zeff) ) z_zeff = zeff

    zNZ = 0.58d0 + 0.74d0 / (0.76d0 + z_zeff)
    zlnL = 17.d0

    IF (ne.gt.0.d0 .and. te.gt.0.d0) THEN
       zlnL = 31.3d0 - log(sqrt(ne)/te)
    ENDIF

    sigsptz = 1.9012d4 * abs(te)**1.5d0 / (z_zeff * zNZ * zlnL)

    nuestar = 0.01d0

    IF (PRESENT(q) .AND. PRESENT(R) .AND. PRESENT(eps)) &
          nuestar = 6.921d-18 * q * R * ne * z_zeff * zlnL / (te*te * eps**1.5d0)

    zft33eff = ft / &
          (1.d0+(0.55d0-0.1d0*ft)*sqrt(nuestar) &
           + 0.45d0*(1.d0-ft)*nuestar/z_zeff**1.5d0)

    signeo = sigsptz * (1.d0 - zft33eff*(1.d0+0.36d0/z_zeff &
                                 - zft33eff*(0.59d0/z_zeff - 0.23d0/z_zeff*zft33eff)))

  end subroutine sigmaneo


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_JparB(JparB_ave, bpol, bphi, J_par, nn)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: bpol, bphi, J_par
	real(kind=dkind), dimension(1:enq) :: JparB_ave
	real(kind=dkind), dimension(1:nn,1:nn) :: thing

	integer :: i,j

	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(bpol, nn, -1.d0)

	endif

	! then create the integrand thing, and give it to the averaging routine

	do j = 1,nn
	do i = 1,nn

		thing(i,j) = J_par(i,j) * sqrt( bpol(i,j)**2 + bphi(i,j)**2 )

	enddo
	enddo

	call surf_ave(bpol, thing, nn, -1.d0, JparB_ave)

	continue

	return

end subroutine get_JparB

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_thing_ave(thing_ave, base, thing, nn, base_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBS2DR

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: base, thing
	real(kind=dkind), dimension(1:enq) :: thing_ave
	real(kind=dkind), dimension(1:nn,1:nn) :: thing_to_pass
	real(kind=dkind) :: base_option

	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(base, nn, base_option)

	endif

	! then call the averaging routine

	call surf_ave(base, thing, nn, base_option, thing_ave)

	return

end subroutine get_thing_ave

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine surf_ave_setup(ave_base, nn, ave_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: ave_base
	real(kind=dkind), dimension(:), allocatable :: integrand
	real(kind=dkind) :: ave_option !integral of ave_base^ave_option

	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, baseloc

	integer :: i, k

	! first set up the various arrays (note that the weighting stuff will stay)

	allocate(w_ave_int(1:n_ave_int))
	allocate(t_ave_int(1:n_ave_int))

	call set_weights(n_ave_int,0.d0,2.d0*pi,w_ave_int,t_ave_int)

	allocate (integrand(1:n_ave_int))

	allocate(surf_ave_base(1:enq))

	allocate(bscoef_boot_base(1:nn,1:nn))

	! then set up the interpolation of the incoming quantity

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),ave_base,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_boot_base)

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,n,n,bscoef_boot_base)

			integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * baseloc**ave_option  &
									/ (2.d0*pi)

		enddo

		! then calculate the integral

		call integrate(n_ave_int, integrand, w_ave_int, surf_ave_base(i))

		continue

	enddo



	deallocate(integrand)

	continue

	return

end subroutine surf_ave_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine surf_ave(ave_base, ave_thing, nn, ave_option, ave_result)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: ave_base, ave_thing
	real(kind=dkind), dimension(:), allocatable :: integrand
	real(kind=dkind), dimension(1:enq) :: ave_result
	real(kind=dkind) :: ave_option !integral of ave_thing * ave_base^ave_option

	real(kind=dkind), dimension(1:nn,1:nn) :: bscoef_thing
	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, baseloc, thingloc

	integer :: i, k

	! first set up the various arrays (note that the weighting stuff is given)

	allocate (integrand(1:n_ave_int))

	! then set up the interpolation of the incoming quantities

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),ave_base,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_boot_base)

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),ave_thing,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_thing)

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,n,n,bscoef_boot_base)

			thingloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,n,n,bscoef_thing)

			integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * thingloc * baseloc**ave_option  &
									/ (2.d0*pi)

		enddo

		! then calculate the integral

		call integrate(n_ave_int, integrand, w_ave_int, ave_result(i))

		ave_result(i) = ave_result(i) / surf_ave_base(i)

		continue

	enddo



	deallocate(integrand)

	continue

	return

end subroutine surf_ave




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_new_avs(p_ave, dens_ave, base, pres, dens, nn, base_option,  &
										p_prim_ave, dens_prim_ave)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! note: only electron pressure and density get here; ion's are proportional to them

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBS2DR

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: base, pres, dens
	real(kind=dkind), dimension(1:enq) :: p_ave, dens_ave
	real(kind=dkind), dimension(1:nn,1:nn) :: thing_to_pass
	real(kind=dkind), dimension(1:enq) :: p_prim_ave, dens_prim_ave
	real(kind=dkind) :: base_option

	real(kind=dkind), dimension(:,:), allocatable :: integrand
	real(kind=dkind), dimension(1:nn,1:nn) :: bscoef_pres, bscoef_dens
	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, bigrlocprim, zlocprim,  &
								baseloc, thingloc, pres_U, pres_L, dens_U, dens_L,  &
								psi_U, psi_L, cos_n, sin_n, der_dens, der_pres

	real(kind=dkind), dimension(1:2) :: U_point, L_point
	! points on the normal to the magnetic surface

	real(kind=dkind) :: dist

	integer :: i, j, k


	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(base, nn, base_option)

	endif

	! then call the averaging routine for pressure and density

	call surf_ave(base, pres, nn, base_option, p_ave)
	call surf_ave(base, dens, nn, base_option, dens_ave)

	! then take care of the derivative part
	! the integration procedure is kept here to avoid additional array copying

	dist = sqrt(dx**2+dz**2)/10.d0

	allocate (integrand(1:2,1:n_ave_int))
	! "1" is pressure, "2" is density

	! then set up the interpolation of the incoming quantities

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),pres,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_pres)

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),dens,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_dens)

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_boot_base)

			!------------------normal and derivative------------------
			cos_n = zlocprim/sqrt(zlocprim**2+bigrlocprim**2)
			sin_n = -bigrlocprim/sqrt(zlocprim**2+bigrlocprim**2)

			U_point(1) = bigrloc + dist*cos_n
			U_point(2) = zloc + dist*sin_n

			L_point(1) = bigrloc - dist*cos_n
			L_point(2) = zloc - dist*sin_n

			psi_L = DBS2VL(L_point(1),L_point(2),ord_loc,ord_loc,xknot_mag, &
											zknot_mag,nn,nn,bscoef_psi)

			psi_U = DBS2VL(U_point(1),U_point(2),ord_loc,ord_loc,xknot_mag, &
											zknot_mag,nn,nn,bscoef_psi)

			pres_L = DBS2VL(L_point(1),L_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_pres)

			pres_U = DBS2VL(U_point(1),U_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_pres)

			dens_L = DBS2VL(L_point(1),L_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_dens)

			dens_U = DBS2VL(U_point(1),U_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_dens)

			der_pres = (pres_U-pres_L)/(psi_U-psi_L)

			der_dens = (dens_U-dens_L)/(psi_U-psi_L)


!			thingloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
!											zknot_mag,n,n,bscoef_thing)

			integrand(1,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * der_pres * baseloc**base_option  &
									/ (2.d0*pi)

			integrand(2,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * der_dens * baseloc**base_option  &
									/ (2.d0*pi)

		enddo

		! then calculate the integrals

		call integrate(n_ave_int, integrand(1,:), w_ave_int, p_prim_ave(i))

		p_prim_ave(i) = p_prim_ave(i) / surf_ave_base(i)

		call integrate(n_ave_int, integrand(2,:), w_ave_int, dens_prim_ave(i))

		dens_prim_ave(i) = dens_prim_ave(i) / surf_ave_base(i)

		continue

	enddo



	deallocate(integrand)

	continue


	return

end subroutine get_new_avs


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine bootstrap_cleanup
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


	deallocate(psival)
	deallocate(qval)

	deallocate(J_boot)
	deallocate(eff_trap)
	deallocate(surf_length)
	if(allocated(B2_ave)) deallocate(B2_ave)
	if(allocated(Bm2_ave)) deallocate(Bm2_ave)
	deallocate(B2_hat_ave)
	deallocate(boot_ne)
	deallocate(boot_pe)
	deallocate(boot_Te)
	deallocate(boot_neprim)
	deallocate(boot_peprim)
	deallocate(boot_Teprim)
	deallocate(boot_ni)
	deallocate(boot_pi)
	deallocate(boot_Ti)
	deallocate(boot_niprim)
	deallocate(boot_piprim)
	deallocate(boot_Tiprim)
	if(allocated(boot_tor_flux)) deallocate(boot_tor_flux)
	if(allocated(boot_tor_rho)) deallocate(boot_tor_rho)
	if(allocated(boot_fhat)) deallocate(boot_fhat)
	if(allocated(boot_grho_ov_B2)) deallocate(boot_grho_ov_B2)
	deallocate(inv_asp_ratio)
	deallocate(boot_Bmin)
	deallocate(boot_Bp_at_Bmin)
	deallocate(boot_R_of_Bmin)
	deallocate(boot_Rcenter)
	deallocate(el_resistivity)

	deallocate(Rleft)
	deallocate(Rright)

	deallocate(rtab)
	deallocate(bigrtab)
	deallocate(ztab)
	deallocate(thetatab)

	deallocate(xknot_mag)
	deallocate(zknot_mag)

	deallocate(surf_ave_base)
	deallocate(JparB_ave)

	deallocate(bscoef_boot_base)

end subroutine bootstrap_cleanup



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine grid_delta(aa,f,df)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	real(kind=dkind)  aa, f, df

	if(grid_type==1) then

		f = g_ratio * (1.d0-aa**((n_temp-1)/2)) - 1.d0 + aa**(n_temp-1)

		df = (n_temp-1.d0) * aa**(n_temp-2) - g_ratio*((n_temp-1.d0)/2.d0) * aa**((n_temp-3)/2)

	elseif(grid_type==2) then

		f = n_temp*g_ratio/(1.d0-g_ratio) * (1.d0-aa)*aa**(n_temp-1) - (1.d0-aa**n_temp)

!		df = n_temp*aa**(n_temp-1) * (1.d0 + g_ratio/(1.d0-g_ratio) * (&
!				(n_temp-3.d0)/aa**3.d0 - (n_temp-2.d0)/aa**2 ) )

		df = n_temp*aa**(n_temp-1) * (1.d0 + g_ratio/(1.d0-g_ratio) * (&
				(n_temp-1.d0)/g_ratio - n_temp ) )

	endif

	continue

	return

	end subroutine grid_delta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine grid_delta_2(aa,f,df)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	real(kind=dkind)  aa, f, df

	f = aa**n_temp - 1.d0 - g_ratio * (aa-1.d0)

	df = n_temp * aa**(n_temp-1)-g_ratio

	continue

	return

	end subroutine grid_delta_2



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_position(i,j,bound,truebound,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		use constant, only : dx, dz

		integer :: i,j,n
		logical :: bound, truebound
		real(kind=dkind) :: ex,ez


		truebound = .false.

! i,j
		if(sort_grid(i,j)==0) then
		! the point is actually on the boundary

			truebound = .true.
			return

		endif

		if((sort_grid(i,j)==1).or.(sort_grid(i,j)==2)) then

			bound = .false.
			return

		endif

! i,j+1
		if(sort_grid(i,j+1)==1) then

			bound = .true.
			return

		endif

! i,j-1
		if(sort_grid(i,j-1)==1) then

			bound = .true.
			return

		endif

! i+1,j
		if(sort_grid(i+1,j)==1) then

			bound = .true.
			return

		endif

! i-1,j
		if(sort_grid(i-1,j)==1) then

			bound = .true.
			return

		endif

		bound = .false.

		continue

		return

	end subroutine check_position

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine point1(i,j,xQ,zQ,Qind)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		use constant, only : x_coord, z_coord, dx_a, dz_a

		integer :: i,j
		real(kind=dkind) :: xQ,zQ
		integer, dimension(1:2) :: Qind
		integer :: ii,jj

		ii = 0
		jj = 0


		do

!			if( (abs((rmajor - 0.5d0*x_size + (ii-1)*dx)-xQ)<=dx*(1+1.d-9)).and.  &
!				(abs((rmajor - 0.5d0*x_size + ii*dx)-xQ)<=dx*(1+1.d-9))		)then
			if( (abs(x_coord(ii)-xQ)<=dx_a(ii)*(1.d0+1.d-9)).and.  &
				(abs(x_coord(ii+1)-xQ)<=dx_a(ii)*(1.d0+1.d-9))		)then

				Qind(1) = ii
				exit

			endif

			ii = ii + 1

			if((ii>i+2).and.(grid_type==0)) pause 'error in point1, i'

		enddo


		do

!			if( (abs((-0.5d0*z_size + (jj-1)*dz)-zQ)<=dz*(1+1.d-9)).and.  &
!				(abs((-0.5d0*z_size + (jj-1)*dz)-zQ)<=dz*(1+1.d-9))	) then
			if( (abs(z_coord(jj)-zQ)<=dz_a(jj)*(1.d0+1.d-9)).and.  &
				(abs(z_coord(jj+1)-zQ)<=dz_a(jj)*(1.d0+1.d-9))	) then

				Qind(2) = jj
				exit

			endif

			jj = jj + 1

			if((jj>j+2).and.(grid_type==0)) pause 'error in point1, j'

		enddo

		continue

		return

	end subroutine point1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_point1(i,j,inner)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		integer :: i,j
		logical :: inner

		inner = .true.

! i,j
		if(sort_grid(i,j)<=0) then

!			print*, 'error, point 1 is external!'

			inner = .false.
			return

		endif

! i,j+1
		if(sort_grid(i,j+1)<=0) then

			inner = .false.
			return

		endif

! i+1,j+1
		if(sort_grid(i+1,j+1)<=0) then

			inner = .false.
			return

		endif

! i+1,j
		if(sort_grid(i+1,j)<=0) then

			inner = .false.
			return

		endif

		continue

		return

	end subroutine check_point1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_points(i,j,inner)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		integer :: i,j
		logical :: inner
		integer :: how_many

		how_many = 0

! i,j
		if(sort_grid(i,j)<=0) then

			how_many = how_many + 1

		endif

! i,j+1
		if(sort_grid(i,j+1)<=0) then

			how_many = how_many + 1

			if(how_many>1) then

				inner = .false.
				return

			endif

		endif

! i+1,j+1
		if(sort_grid(i+1,j+1)<=0) then

			how_many = how_many + 1

			if(how_many>1) then

				inner = .false.
				return

			endif

		endif

! i+1,j
		if(sort_grid(i+1,j)<=0) then

			how_many = how_many + 1

			if(how_many>1) then

				inner = .false.
				return

			endif

		endif

		inner = .true.

		continue

		return

	end subroutine check_points

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_points_tri(ind1,ind2,ind3,inner)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! checking if there are any external points in the triangle
! for triangular elements:
! since by construction P is NOT a triangle point,
! any ONE point will classify the triangle as external

		integer, dimension(1:2) :: ind1,ind2,ind3
		logical :: inner

! i1, j1
		if(sort_grid(ind1(1),ind1(2))<=0) then

			inner = .false.
			return

		endif

! i2, j2
		if(sort_grid(ind2(1),ind2(2))<=0) then

			inner = .false.
			return

		endif

! i3, j3
		if(sort_grid(ind3(1),ind3(2))<=0) then

			inner = .false.
			return

		endif

		inner = .true.

		continue

		return

	end subroutine check_points_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function tri_area(p1,p2,p3) result(area)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this is the actual area

	real(kind=dkind), dimension(1:2) :: p1, p2, p3
	real(kind=dkind) :: area

	area = abs( p1(1)*p2(2) + p2(1)*p3(2) + p3(1)*p1(2)  &
					- p2(1)*p1(2)- p3(1)*p2(2)- p1(1)*p3(2) )/2.d0

	continue

	return

end function tri_area

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function tri_area_sign(p1,p2,p3) result(area)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! result can be positive or negative, but this is what is needed

	real(kind=dkind), dimension(1:2) :: p1, p2, p3
	real(kind=dkind) :: area

	area = ( p1(1)*p2(2) + p2(1)*p3(2) + p3(1)*p1(2)  &
					- p2(1)*p1(2)- p3(1)*p2(2)- p1(1)*p3(2) )/2.d0

	continue

	return

end function tri_area_sign

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_Q_tri(dist0,dist,iP,jP,tri1,tri2,tri3,xvec,nx,nz,A_return)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! determines Q so that Q is contained in a triangle
! which does NOT contain P
! if bc_option==113,114, all points in the triangle are internal

	integer :: iP, jP
	integer :: nx, nz
	real(kind=dkind) :: dist, dist0, A_return
	real(kind=dkind) :: xvec(2)
	integer, dimension(1:2) :: tri1, tri2, tri3
	integer :: points_temp(1:3,1:2)
	integer, dimension(1:12,1:3) :: i_tri, j_tri
	integer :: i_closest, i_check
	real(kind=dkind) :: d_min, A1, A2, A3, Atot, delta_dist, dist_loc
	real(kind=dkind) :: RQ_temp, ZQ_temp
	real(kind=dkind), dimension(1:2) :: p1, p2, p3, CC
	logical :: incomplete
	integer :: iC, jC
	integer :: i, j, k, iii, jjj

	! --------------- first set triangles (can be done only once at the beginning)---------------
	! there are	TWELVE triangles with any point as one of the vertices
	! givenpoint is always the first for simplicity

	i_tri(:,1) = 0
	j_tri(:,1)  = 0

	! triangle 1

	i_tri(1,2) = 1
	j_tri(1,2) = 0

	i_tri(1,3) = 1
	j_tri(1,3) = 1

	! triangle 2

	i_tri(2,2) = 0
	j_tri(2,2) = 1

	i_tri(2,3) = 1
	j_tri(2,3) = 1

	! triangle 3

	i_tri(3,2) = 1
	j_tri(3,2) = 0

	i_tri(3,3) = 0
	j_tri(3,3) = 1

	! triangle 4

	i_tri(4,2) = -1
	j_tri(4,2) = 1

	i_tri(4,3) = 0
	j_tri(4,3) = 1

	! triangle 5

	i_tri(5,2) = -1
	j_tri(5,2) = 0

	i_tri(5,3) = -1
	j_tri(5,3) = 1

	! triangle 6

	i_tri(6,2) = 0
	j_tri(6,2) = 1

	i_tri(6,3) = -1
	j_tri(6,3) = 0

	! triangle 7

	i_tri(7,2) = -1
	j_tri(7,2) = 0

	i_tri(7,3) = -1
	j_tri(7,3) = -1

	! triangle 8

	i_tri(8,2) = 0
	j_tri(8,2) = -1

	i_tri(8,3) = -1
	j_tri(8,3) = -1

	! triangle 9

	i_tri(9,2) = 0
	j_tri(9,2) = -1

	i_tri(9,3) = -1
	j_tri(9,3) = 0

	! triangle 10

	i_tri(10,2) = 0
	j_tri(10,2) = -1

	i_tri(10,3) = 1
	j_tri(10,3) = -1

	! triangle 11

	i_tri(11,2) = 1
	j_tri(11,2) = 0

	i_tri(11,3) = 1
	j_tri(11,3) = -1

	! triangle 12

	i_tri(12,2) = 1
	j_tri(12,2) = 0

	i_tri(12,3) = 0
	j_tri(12,3) = -1

	! --------------- done with triangles ---------------

	if(dist0<0.2d0*sqrt(dx_a(iP)**2+dz_a(jP)**2)) then

		dist = 0.5d0*sqrt(dx_a(iP)**2+dz_a(jP)**2)

	else

		dist = dist0

	endif

	delta_dist = 1.d-1*sqrt(dx_a(iP)**2+dz_a(jP)**2)

	incomplete = .true.

	do while(incomplete)

		RQ_temp = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
		ZQ_temp = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

		d_min = 1.d9 * (dx+dz)

		do j = jP-5, jP+5

			if((j<1).or.(j>nz)) cycle

			do i = iP-5, iP+5

				if((i<1).or.(i>nx)) cycle

				dist_loc = sqrt((x_coord(i)-RQ_temp)**2+(z_coord(j)-ZQ_temp)**2)

				if(dist_loc<d_min) then

					d_min = dist_loc
					points_temp(1,1) = i
					points_temp(1,2) = j

				endif

			enddo
		enddo

		iC = points_temp(1,1)
		jC = points_temp(1,2)

!		CC(1) = x_coord(iC)
!		CC(2) = z_coord(jC)

		CC(1) = RQ_temp
		CC(2) = ZQ_temp

		! now we found the closest point, let's build a triangle with it
		! that contains Q (4 possible cases)

		! for each triangle, check if point Q2 is inside; if so, check if triangle is internal

		tri_loop: do jjj = 1, 12

			do k = 1, 3

				if((iC+i_tri(jjj,k)<1).or.(iC+i_tri(jjj,k)>nx).or.(jC+j_tri(jjj,k)<1).or.(jC+j_tri(jjj,k)>nz)) then
					! using points outside the grid, skip this case
					cycle tri_loop
				endif

			enddo

			! even before checking if point Q is in triangle, check the triangle (we don't want P to be one of the corners)
			do iii = 1, 3

				if((iC+i_tri(jjj,iii)==iP).and.(jC+j_tri(jjj,iii)==jP)) then
					! point P is one of the corners, reject triangle
					cycle tri_loop
				endif

			enddo

			if((bc_type==113).or.(bc_type==114)) then
			! also reject the triangle if any of the points is external

				do iii = 1, 3
	
					if(sort_grid(iC+i_tri(jjj,iii),jC+j_tri(jjj,iii))<=0) then
						! one of the points is external, reject triangle
						cycle tri_loop
					endif
	
				enddo

			endif

			p1(1) = x_coord(iC+i_tri(jjj,1))
			p2(1) = x_coord(iC+i_tri(jjj,2))
			p3(1) = x_coord(iC+i_tri(jjj,3))

			p1(2) = z_coord(jC+j_tri(jjj,1))
			p2(2) = z_coord(jC+j_tri(jjj,2))
			p3(2) = z_coord(jC+j_tri(jjj,3))

			Atot = tri_area(p1,p2,p3)
			A1 = tri_area(CC,p2,p3)
			A2 = tri_area(p1,CC,p3)
			A3 = tri_area(p1,p2,CC)

			if((A1+A2+A3)/Atot>(1.d0+1.d-12)) cycle
			! point Q is not in triangle

			! if we made it this far, we found our triangle!

			incomplete = .false.

			tri1(1) = iC + i_tri(jjj,1)
			tri1(2) = jC + j_tri(jjj,1)
			tri2(1) = iC + i_tri(jjj,2)
			tri2(2) = jC + j_tri(jjj,2)
			tri3(1) = iC + i_tri(jjj,3)
			tri3(2) = jC + j_tri(jjj,3)

			A_return = tri_area_sign(p1,p2,p3)

			exit tri_loop

		enddo tri_loop

		if(incomplete) then
			dist = dist + delta_dist
		endif

	enddo

	continue

end subroutine get_Q_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!subroutine get_Q_tri_first_attempt(dist0,dist,iP,jP,tri1,tri2,tri3,xvec)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! determines Q so that Q is contained in a triangle
! which does NOT contain P
!
!	integer :: iP, jP
!	real(kind=dkind) :: dist, dist0
!	rel(kind=dkind) :: xvec(2)
!	integer, dimension(1:2) :: tri1, tri2, tri3
!	integer :: points_temp(1:3,1:2)
!	integer :: i_closest, i_farthest, i_check
!	real(kind=dkind) :: d_min, d_max, p_dist(1:3)
!	integer :: i, j, k
!
!	if(dist0<0.2d0*sqrt(dx_a(iP)**2+dz_a(jP)**2)) then
!
!		dist = 0.5d0*sqrt(dx_a(iP)**2+dz_a(jP)**2)
!
!	else
!
!		dist = dist0
!
!	endif
!
!	delta_dist = 1.d-1*sqrt(dx**2+dz**2)
!
!	incomplete = .true.
!
!	RQ_temp = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
!	ZQ_temp = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P
!
!	do while(incomplete)
!
!		i_closest = 1
!		i_farthest = 3
!
!		p_dist(1) = 1.d9 * (dx+dz)
!		p_dist(2) = 2.d9 * (dx+dz)
!		p_dist(3) = 3.d9 * (dx+dz)
!
!		d_min = p_dist(1)
!		d_max = p_dist(3)
!
!		do j = jP-5, jP+5
!
!			if((j<1).or.(j>nz)) cycle
!
!			do i = iP-5, iP+5
!
!				if((i<1).or.(i>nx)) cycle
!
!				dist_loc = sqrt((x_coord(i)-RQ_temp)**2+(z_coord(i)-ZQ_temp)**2)
!
!				if(dist_loc<d_max) then
!
!					p_dist(i_farthest) = dist_loc
!					points_temp(i_farthest,1) = i
!					points_temp(i_farthest,2) = j
!
!					d_min = min(p_dist(1),p_dist(2),p_dist(3))
!					d_max = max(p_dist(1),p_dist(2),p_dist(3))
!
!					do k = 1, 3
!
!						if(p_dist(k)==d_min) i_closest = k
!						if(p_dist(k)==d_max) i_farthest = k
!
!					enddo
!
!				endif
!
!			enddo
!		enddo
!
!		incomplete = .false.
!
!		do k = 1, 3
!
!			if((points(k,1)==iP).and.(points(k,2)==jP)) then
!				incomplete = .true.
!				dist = dist + delta_dist
!			endif
!
!		enddo
!
!	enddo
!
!	i_check = max(abs(points_temp(1,1)-points_temp(2,1)),abs(points_temp(1,1)-points_temp(3,1))) *  &
!				max(abs(points_temp(1,2)-points_temp(2,2)),abs(points_temp(1,2)-points_temp(3,2)))
!
!	if(i_check/=1) then
!	 not done, yet: the three points do not make a triangle,
!	 look for a triangle in a different way starting from the closest point (4 possible cases)
!
!		incomplete = .true.
!
!		do while(incomplete)
!
!			tri1(1) = points_temp(i_closest,1)
!			tri1(2) = points_temp(i_closest,2)
!
!			 case #1
!
!		enddo
!
!	endif
!
!	tri1(1) = points_temp(1,1)
!	tri1(2) = points_temp(1,2)
!	tri2(1) = points_temp(2,1)
!	tri2(2) = points_temp(2,2)
!	tri3(1) = points_temp(3,1)
!	tri3(2) = points_temp(3,2)
!
!	continue
!
!end subroutine get_Q_tri_first_attempt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_interface(psi,n,inorm)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use exp_data, only : r_data_psi, xknot_psi, zknot_psi, psi_bscoef
	use triangularity, only : theta_temp, r_ord, r_data
	use pseudo_IMSL, only : DBSNAK, DBSINT

	implicit none

	integer :: i, k, nth, alloc_stat, error
	real(kind=dkind) :: rg1, rg2, smallr1, smallr2, r_orp
	integer :: n
	real(kind=dkind), dimension(1:n,1:n) :: psi
	real(kind=dkind) :: inorm
	real(kind=dkind) :: fnew
	logical :: look_for_r
!	external psi_sol

	! first update fraction


	if (n<inter_switch) return

	r_orp = 0.25d0

	inorm = 0.d0

	call psi_interp_setup(psi,n)

	nth = theta_points2

	allocate (r_data_psi(nth+r_ord,3),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in r_data_psi"
		 pause
		 stop
	endif


	do i=1,nth

		theta_temp = (i-1.d0)/(nth-1.d0)*2.d0*pi

		k = 0

		smallr1 = dbsval(theta_temp, r_ord, r_data(1:theta_points2+r_ord,6),  &
			theta_points2, r_cscoef(2,1:theta_points2) )

		smallr2 = dbsval(theta_temp, r_ord, r_data(1:theta_points1+r_ord,3),  &
			theta_points1, r_cscoef(1,1:theta_points1) )

!!$		if(i==1) then
!!$			rg1 = 0.d0
!!$			rg2 = smallr1 * 1.05d0
!!$		else
!!$			rg1 = 0.d0
!!$			rg2 = r_data_psi(i-1,2)*1.09d0
!!$			if(rg2>=smallr2) rg2 = 0.99d0*smallr2
!!$		endif

		rg1 = 0.d0
		rg2 = smallr1 * .5d0

		r_data_psi(i,1) = theta_temp

		look_for_r = .true.

		do while(look_for_r)
			
			call secroot(psi_sol,rg1,rg2,r_data_psi(i,2),  &
					100,1.d-9,1.d-9, error)
			
			if(error==0) then

				look_for_r = .false.

			elseif(error==1) then
		! secroot failed because the root was not bracketed, try again
!!$			rg2 = smallr2 * 0.99**k
!!$			k = k+1
				rg2 = rg2 * 1.01d0

				!if((rg2<0.5d0*smallr1).or.(rg2>0.995d0*smallr2)) then
				if(rg2>0.995d0*smallr2) then
					r_data_psi(i,2) = smallr2
					print*, 'problem in update_interface, theta = ', theta_temp
					! no need to say this should NOT happen
					look_for_r = .false.
				endif

			endif

		enddo

	enddo

	! relaxation process

	do i = 1,nth

		smallr1 = dbsval(r_data_psi(i,1), r_ord, r_data(1:theta_points2+r_ord,6),  &
			theta_points2, r_cscoef(2,1:theta_points2) )

		inorm = max(inorm,dabs( (r_orp*r_data_psi(i,2) + (1.d0-r_orp)*smallr1)/smallr1 )-1.d0 )

		r_data_psi(i,2) = r_orp*r_data_psi(i,2) + (1.d0-r_orp)*smallr1

	enddo

	! end of relaxation process

	do i=1,theta_points2 + r_ord
	! just in case

		r_data(i,4) = 0.d0
		r_data(i,5) = 0.d0
		r_data(i,6) = 0.d0

	enddo

	r_cscoef(2,:) = 0.d0

	do i=1,theta_points2

		r_data(i,4) = r_data_psi(i,1)
		r_data(i,5) = r_data_psi(i,2)

	enddo

	call DBSNAK(theta_points2, r_data(1:theta_points2,4),  &
					r_ord,r_data(1:theta_points2+r_ord,6))

	call DBSINT(theta_points2, r_data(1:theta_points2,4),  &
		 r_data(1:theta_points2,5), r_ord,  &
		 r_data(1:theta_points2+r_ord,6),  &
		 r_cscoef(2,1:theta_points2))

	deallocate(r_data_psi)
	deallocate(psi_bscoef)
	deallocate(xknot_psi)
	deallocate(zknot_psi)

	print*, 'interface error:', inorm
	print*,'      '


	continue

	return

end subroutine update_interface

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psi_interp_setup(psi,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : x_coord, z_coord
	use exp_data, only :  s_ord, nx_FLOW, nz_FLOW,  &
							xknot_psi, zknot_psi, psi_bscoef

	use pseudo_IMSL, only : DBS2IN, DBS2GD, DBSNAK

	implicit none

	integer :: alloc_stat
	integer :: i, j
	integer :: n
	real(kind=dkind), dimension(1:n,1:n) :: psi

	nx_FLOW = n
	nz_FLOW = n

!-----------------------------------------
! interpolation setup


	allocate (psi_bscoef(1:nx_FLOW,1:nz_FLOW),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi_bscoef"
		 pause
		 stop
	endif


	allocate (xknot_psi(1:nx_FLOW+s_ord),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in xknot_psi"
		 pause
		 stop
	endif

	allocate (zknot_psi(1:nz_FLOW+s_ord),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in zknot_psi"
		 pause
		 stop
	endif

	call DBSNAK(nx_FLOW,x_coord,s_ord,xknot_psi)
	call DBSNAK(nz_FLOW,z_coord,s_ord,zknot_psi)

	! (this 2 define the nodes)

! end of interpolation setup
!-----------------------------------------
!-----------------------------------------
! set psi


	call DBS2IN(nx_FLOW,x_coord,nz_FLOW,z_coord,psi,nx_FLOW,  &
				s_ord,s_ord,xknot_psi,zknot_psi,psi_bscoef(:,:) )


!-----------------------------------------



	continue

	return

end subroutine psi_interp_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_sort_grid(psi,nx,nz,inorm) ! Added inorm argument (NEW May 2025)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The code will get here to update which points are in the main plasma and which ones are in the open field line region.
! For now, only the bc_type==7 option is considered ( March 12 2021).

	use pseudo_IMSL, only : dbs2vl

	implicit none

!	real(kind=dkind), parameter :: small = 1.d-2

	integer :: nx,nz
	real(kind=dkind), intent(in) :: psi(1:nx,1:nz)
	real(kind=dkind) :: inorm
	integer :: i,j,p,k
	real(kind=dkind) :: ex, ez, rminor, dx, dz, dummy, r_in
	real(kind=dkind) :: small = 1.d-2
	integer, allocatable :: sort_grid_old(:,:), changed_index(:,:)
	integer :: index_counter, internal_counter
	integer :: n_check = 10
	real(kind=dkind) :: xloc, zloc, psiloc

!	sort_grid = -1
!	January 21 2022: commented the previous line.
!	STILL NEED TO CHECK THAT THIS DOES NOT BREAK TRI_TYPE=-1, -2

! February 4 2022: logical variables should not be needed anymore.
!!$	if((tri_type==-2).and.(initialize_r_for_tri_type_m2)) then
!!$		call initialize_r_free_boundary(psi(:,:),nx)
!!$		initialize_r_for_tri_type_m2 = .false.
!!$		tri_type_m2_ready = .true.
!!$	endif

! September 5 2022: Coding in tri_type -4. Tri_type -2 and -3 do not work well and will be ignored for now.
! Currently under modification by Ian (NEW May 2025)
!!$	if((tri_type==-2).and.(tri_type_m2_ready)) then
	if((tri_type==-2).or.(tri_type==-3)) then
		! Proceed like in old tri_type==13, interpolate psi to find the boundary
		call update_interface(psi(:,:),nx,inorm)
	endif

	if(tri_type==-4) then
		allocate(sort_grid_old(1:nx,1:nz))
		sort_grid_old(:,:) = sort_grid(1:nx,1:nz)
		allocate(changed_index(2,nx*nz))
	endif

	do j=1,nz
	do i=1,nx

		if(tri_type==13) then
			cycle
		endif

		if((tri_type==-1).or.(tri_type==-4)) then
			call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)
		elseif(tri_type==-2) then
			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
		elseif(tri_type==-3) then
			call radius_both(x_coord(i),z_coord(j),ex,ez,r_in,rminor) ! rminor is r_out
		endif

		if((ex*ex + ez*ez) > rminor**2) then

			sort_grid(i,j) = -1

		elseif( (ex*ex + ez*ez) >= (rminor - small*sqrt(dx_a(i)**2+dz_a(j)**2) )**2 ) then

			sort_grid(i,j) = 0

		else

			if(tri_type==-1) then

				 if(((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0)) then  !&
						!.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
						! Line above modified by Ian

					sort_grid(i,j) = 1 ! External zone

				else

					sort_grid(i,j) = 2 ! Internal zone

				endif

			elseif(tri_type==-2) then

				if(p==2) then
					sort_grid(i,j) = 2
					! this sets inner zone (plasma) to "2"
				else
					sort_grid(i,j) = 1
				endif

			elseif(tri_type==-3) then

				 if((((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
						.or. ((bc_type==7).and.(psi(i,j)>0.d0).and.((ex*ex + ez*ez) > rminor**2))) then

					sort_grid(i,j) = 1 ! External zone

				else

					sort_grid(i,j) = 2 ! Internal zone

				endif

			elseif(tri_type==-4) then

				 if(((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  then

					sort_grid(i,j) = 1 ! External zone

				else

					sort_grid(i,j) = 2 ! Internal zone

				endif

			endif

		endif

	enddo
	enddo

	! March 29 2022: added check to fix radius inaccuracies in the corners

	if((tri_type==-2).or.(tri_type==-3)) then

		do i=1,nx

			if(sort_grid(i,1)>0) then
				sort_grid(i,1)=0
			endif

			if(sort_grid(i,nz)>0) then
				sort_grid(i,nz)=0
			endif

		enddo

	endif

	if(tri_type==-4) then

		index_counter = 0
		internal_counter = 0

		do i = 1,nz
		do j = 1,nx

			if (sort_grid(i,j)/=sort_grid_old(i,j)) then
			! Mark all points for debugging purposes, but only check the ones marked as internal

				index_counter = index_counter+1
				changed_index(1,index_counter) = i
				changed_index(2,index_counter) = j

				if(sort_grid(i,j)==2) then
					internal_counter = internal_counter + 1
				endif

			endif

		enddo
		enddo

		if(internal_counter>0) then
		! Check that all external points are classified as such

			! First set up the psi interpolation
			call psi_interp_setup(psi,nx)

			do k = 1,index_counter

				i = changed_index(1,k)
				j = changed_index(2,k)

				if(sort_grid(i,j)==2) then
					! Check the point.

					do p = 1, n_check

						xloc = x_coord(i) - (x_coord(i)-x_coord(ipsic))/(n_check+1.d0)*p
						zloc = z_coord(j) - (z_coord(j)-z_coord(jpsic))/(n_check+1.d0)*p

						psiloc =  DBS2VL(xloc,zloc,s_ord,s_ord,xknot_psi,zknot_psi, &
								nx_FLOW,nz_FLOW,psi_bscoef)

						if(psiloc<=0.d0) then
						! Mark the point as external

							sort_grid(i,j) = 1
							exit

						endif

					enddo

				endif

			enddo

			deallocate(psi_bscoef)
			deallocate(xknot_psi)
			deallocate(zknot_psi)

		endif

		deallocate(sort_grid_old)
		deallocate(changed_index)

	endif

!	if(tri_type==11) then
!
!		do j=1,nz
!		do i=1,nx
!
!			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
!			sort_grid(i,j,1) = p
!
!		enddo
!		enddo
!
!	elseif(tri_type==13) then
!
!		do j=1,nz
!		do i=1,nx
!
!			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
!			if((p==2).and.(sort_grid(i,j,0)==1)) sort_grid(i,j) = 2
!			! this sets inner zone (plasma) to "2"
!
!		enddo
!		enddo
!
!	endif

	open(33,file='grid.plt')
	open(34,file='grid.csv')

	write(33,*)'TITLE="grid"'
	write(33,*)'Variables =" x ","y", "boh"'
!	write(33,*)'Variables =" R [m] ","z [m]", "boh"' ! visit does not like this
	write(33,*)'ZONE I=',nx,',J=',nz,',F=Point'

	do j = 1,nz
	do i = 1,nx

		write(33,*) x_coord(i), z_coord(j),  sort_grid(i,j)
		write(34,*) x_coord(i), ",", z_coord(j), ",",  sort_grid(i,j)

	enddo
	enddo

	close(33)
	close(34)

	continue

end subroutine update_sort_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_sort_grid_old(psi,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The code will get here to update which points are in the main plasma and which ones are in the open field line region.
! For now, only the bc_type==7 option is considered ( March 12 2021).

	implicit none

!	real(kind=dkind), parameter :: small = 1.d-2

	integer :: nx,nz
	real(kind=dkind), intent(in) :: psi(1:nx,1:nz)
	integer :: i,j,p
	real(kind=dkind) :: ex, ez, rminor, dx, dz, dummy
	real(kind=dkind) :: small = 1.d-2

	sort_grid = -1

	do j=1,nz
	do i=1,nx

		call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)

		if((ex*ex + ez*ez) > rminor**2) then

			sort_grid(i,j) = -1

		elseif( (ex*ex + ez*ez) >= (rminor - small*sqrt(dx_a(i)**2+dz_a(j)**2) )**2 ) then

			sort_grid(i,j) = 0

		 elseif((((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.45d0*z_size))) then
!				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then

			sort_grid(i,j) = 1 ! External zone

		else

			sort_grid(i,j) = 2 ! Internal zone

		endif

	enddo
	enddo


	open(33,file='grid.plt')

	write(33,*)'TITLE="grid"'
	write(33,*)'Variables =" R [m] ","z [m]", "boh"'
	write(33,*)'ZONE I=',nx,',J=',nz,',F=Point'

	do j = 1,nz
	do i = 1,nx

		write(33,*) x_coord(i), z_coord(j),  sort_grid(i,j)

	enddo
	enddo

	close(33)

	continue

end subroutine update_sort_grid_old

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_minor_radius
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this calculates the minor radius for bootstrap calculation with NCLASS.
! the routine assumes |delta|<1; if needed, that will be patched later

	real(kind=dkind) :: rloc, xloc, zloc, thetaloc
	real(kind=dkind) :: Rmin, Rmax
	integer :: i
	integer :: n_angles = 200

	Rmin = 2.d0*rcenter
	Rmax = rcenter - 2.d2*x_size

	do i = 1,n_angles

		thetaloc = i*2.d0*pi/n_angles

		call radius_theta(thetaloc,rloc,xloc,zloc)

		Rmin = min(Rmin, xloc)
		Rmax = max(Rmax, xloc)

	enddo

	a_elps = (Rmax-Rmin)/2.d0

	continue


	end subroutine get_minor_radius


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius(i,j,nx,nz,ex,ez,r,dxx,dzz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine computes the minor radius r of the plasma
	! as a function of the angle AND of the triangularity option

		use triangularity
		use constant

		integer, intent(in) :: i,j,nx,nz
		real (kind=dkind) :: dxx,dzz,bigR,Z, diff
		real (kind=dkind), intent(out) :: ex,ez,r
		integer k
		real (kind=dkind) :: r_in,r_out,angle

!		ex = ((i-1)*dx - 0.5d0*x_size)
!		ez = ((j-1)*dz - 0.5d0*z_size)

		if( (i==0).or.(j==0).or.(i>nx).or.(j>nz) ) then
		! somehow the code got here with an index exceeding the grid dimension
		! assign this to be an external point

			ex = 1.d0
			ez = 1.d0
			r = 1.d-2
			return

		endif

		ex = x_coord(i)-rmajor
		ez = z_coord(j)

		if (ex==0.d0) then
			angle = pi/2.d0 * dsign(1.d0,ez)
		else
			angle = datan2(ez,ex)
		endif

		r = 0.d0
		bigR = 0.d0
		Z = 0.d0

		if(tri_type==-3) then

			! plasma - vacuum case

			if(angle<0.d0) angle = angle+2.d0*pi

			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )
			
			r = r_in

		elseif(tri_type==0) then

!			ex = ((i-1.0d0)/(nx - 1.0d0) - 0.5d0)*x_size/a_elps
!			ez = ((j-1.0d0)/(nz - 1.0d0) - 0.5d0)*z_size/b_elps
			ex = (x_coord(i)-rmajor)/a_elps
			ez = z_coord(j)/b_elps

			r = 1.d0

		elseif(tri_type==1) then

			if((dsin(angle)>=0.d0)) then
				do k=0,n_tri
					r = r + rcoeff_u(k)*dcos(k*angle)
				enddo
			else
          				do k=0,n_tri
					r = r + rcoeff_d(k)*dcos(k*angle)
				enddo
			endif

		elseif(tri_type==2) then

			if((dsin(angle)>=0.d0)) then
				r = smallr0*dsqrt( ellipt**2*(dsin(angle))**2  &
						+(dcos(angle+asin_d_up*dsin(angle)))**2 )
			else
          		r = smallr0*dsqrt( ellipt**2*(dsin(angle))**2  &
						+(dcos(angle+asin_d_down*dsin(angle)))**2 )
			endif

		elseif(tri_type==3) then

			if((dsin(angle)>=0.d0)) then
				r = smallr0/dsqrt( dsin(angle)**2/ellipt**2  &
						+dcos(angle)**2-delta_u_o4*(dcos(3*angle)-dcos(angle)) )
			else
          		r = smallr0/dsqrt( dsin(angle)**2/ellipt**2  &
						+dcos(angle)**2-delta_d_o4*(dcos(3*angle)-dcos(angle)) )
			endif

		elseif(tri_type==4) then

!			if (angle>pi/2) angle=angle-2.d0*pi
			if (angle>1.7451968605014) angle=angle-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,angle,r)

		elseif((tri_type==-1).or.(tri_type==-4)) then
			ex = dabs(((i-1.0d0)/(nx - 1.0d0) - 0.5d0)*x_size/a_elps)
			ez = dabs(((j-1.0d0)/(nz - 1.0d0) - 0.5d0)*z_size/b_elps)
			ex=dmax1(ex,ez)
			ez=0.d0
			r=1.d0

        elseif(tri_type==5) then

       		if((dsin(angle)>=0.d0)) then
				do k=0,n_tri
					r = r + rcos_u(k)*dcos(k*angle)+rsin_u(k)*dsin(k*angle)
				enddo
			else
                do k=0,n_tri
					r = r + rcos_d(k)*dcos(k*angle)+rsin_d(k)*dsin(k*angle)
				enddo
			endif

	    elseif(tri_type==6) then

			if (angle>pi/2.d0) angle=angle-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,angle,r)

		elseif(tri_type==8) then

			if(angle>=r_data(theta_points,1)) then

				angle = angle - 2.d0*pi

			elseif(angle<r_data(1,1)) then

				angle = angle + 2.d0*pi

			endif

			r = dbsval(angle, r_ord, r_data(:,3),  &
						theta_points, r_cscoef(1,1:theta_points) )

		elseif((tri_type==13).or.(tri_type==-2)) then

			! plasma - vacuum case

			if(angle<0.d0) angle = angle+2.d0*pi

			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
							theta_points1, r_cscoef(1,1:theta_points1) )

			r = r_out


		endif

		continue

	end subroutine radius


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_1_3(x,z,ex,ez,thet,r,zone,alph,rprim,rsec)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! the code used to get here only for tri_type=11
	! but now this is used to get rprim for bcs for tri_type=8, too (LG 4/14/2006)
	! now, here "4" is the zone inside the inner boundary ("hole")
	! and "1" the zone outside the outer boundary ("external vacuum")

		use triangularity
		use constant

		real (kind=dkind), intent(in) :: x,z
		real (kind=dkind), intent(out) :: r
		real (kind=dkind) :: r_in,r_out,ex,ez, thet, angle
		integer :: zone
		real(kind=dkind) :: alph, rprim, rsec

		ex = x - rmajor
		ez = z

		if (ex==0.d0) then
			thet = pi/2.d0 * dsign(1.d0,ez)
			angle = thet
		else
			thet = datan2(ez,ex)
			angle = thet
		endif

		r = 0.d0
		zone = 1
		alph = 1.d0

		if(tri_type==8) then

			! r of theta

			if(angle>=r_data(theta_points,1)) then

				angle = angle - 2.d0*pi

			elseif(angle<r_data(1,1)) then

				angle = angle + 2.d0*pi

			endif

			r = dbsval(angle, r_ord, r_data(:,3),  &
				theta_points, r_cscoef(1,1:theta_points) )

			rprim = dbsder(1,angle, r_ord, r_data(:,3),  &
					theta_points, r_cscoef(1,1:theta_points) )

			rsec = dbsder(2,angle, r_ord, r_data(:,3),  &
					theta_points, r_cscoef(1,1:theta_points) )


		elseif((tri_type==13).or.(tri_type==-3)) then

			! plasma - vacuum case

			if(angle<0.d0) angle = angle + 2.d0*pi

			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
				theta_points1, r_cscoef(1,1:theta_points1) )

			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

			if( (ex**2+ez**2)>r_in**2 ) then

				zone = 1

			else

				zone = 2

			endif

			! the actual radius is always the external radius,
			! the internal one is a separation!
			r = r_out
			rprim = dbsder(1,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
					theta_points1, r_cscoef(1,1:theta_points1) )
			rsec = dbsder(2,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
					theta_points1, r_cscoef(1,1:theta_points1) )



		else

			print*, 'error in radius_1_3: tri_type=', tri_type
			pause
			stop

		endif

		if(alph/=1.d0) then
			continue
		endif

		thet = angle

		continue

	end subroutine radius_1_3



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_theta(angle,r,x,z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine computes the minor radius r of the plasma
	! as a function of the angle AND of the triangularity option
	! 1/4/2007: added x and z for eqdsk output

		use triangularity
		use constant

		real (kind=dkind), intent(out) :: r,x,z
		integer k
		real(kind=dkind) :: angle


		r = 0.d0
		theta = angle

		if(tri_type==0) then

			r = a_elps*b_elps/sqrt( (a_elps*sin(angle))**2 +  &
						(b_elps*cos(angle))**2 )

		elseif(tri_type==1) then

			if((dsin(theta)>=0.d0)) then
				do k=0,n_tri
					r = r + rcoeff_u(k)*dcos(k*theta)
				enddo
			else
          		do k=0,n_tri
					r = r + rcoeff_d(k)*dcos(k*theta)
				enddo
			endif

		elseif(tri_type==2) then

			if((dsin(theta)>=0.d0)) then
				r = smallr0*dsqrt( ellipt**2*(dsin(theta))**2  &
						+(dcos(theta+asin_d_up*dsin(theta)))**2 )
			else
          		r = smallr0*dsqrt( ellipt**2*(dsin(theta))**2  &
						+(dcos(theta+asin_d_down*dsin(theta)))**2 )
			endif

		elseif(tri_type==3) then

			if((dsin(theta)>=0.d0)) then
				r = smallr0/dsqrt( dsin(theta)**2/ellipt**2  &
						+dcos(theta)**2-delta_u_o4*(dcos(3*theta)-dcos(theta)) )
			else
          		r = smallr0/dsqrt( dsin(theta)**2/ellipt**2  &
						+dcos(theta)**2-delta_d_o4*(dcos(3*theta)-dcos(theta)) )
			endif

		elseif(tri_type==4) then

!			if (theta>pi/2) theta=theta-2.d0*pi
			if (theta>1.7451968605014) theta=theta-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,theta,r)

		elseif(tri_type==-1)then
			! square boundary, this will not be needed

        elseif(tri_type==5) then

       		if((dsin(theta)>=0.d0)) then
				do k=0,n_tri
					r = r + rcos_u(k)*dcos(k*theta)+rsin_u(k)*dsin(k*theta)
				enddo
			else
                do k=0,n_tri
					r = r + rcos_d(k)*dcos(k*theta)+rsin_d(k)*dsin(k*theta)
				enddo
			endif

	    elseif(tri_type==6) then

			if (theta>pi/2.d0) theta=theta-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,theta,r)

		elseif(tri_type==8) then

			if(theta>=r_data(theta_points,1)) then

				theta = theta - 2.d0*pi

			elseif(theta<r_data(1,1)) then

				theta = theta + 2.d0*pi

			endif

			r = dbsval(theta, r_ord, r_data(:,3),  &
				theta_points, r_cscoef(1,1:theta_points) )


		elseif((tri_type==13).or.(tri_type==-2))  then

			! plasma - vacuum case
			! this time we want the plasma radius!

			r = dbsval(theta, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

		endif

		x = rmajor + r*cos(angle)
		z = r*sin(angle)

		continue

	end subroutine radius_theta


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_prim_theta(angle,rp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine computes the first derivative
	! of the minor radius r of the plasma as a function of the angle
	! only for numerical input

		use triangularity
		use constant

		real (kind=dkind), intent(out) :: rp
		real(kind=dkind) :: angle


		rp = 0.d0
		theta = angle

		if(tri_type==8) then
			continue
		else
			return
		endif

		if(theta>=r_data(theta_points,1)) then

			theta = theta - 2.d0*pi

		elseif(theta<r_data(1,1)) then

			theta = theta + 2.d0*pi

		endif

			rp = dbsder(1,theta, r_ord, r_data(:,3),  &
				theta_points, r_cscoef(1,1:theta_points) )


		continue

	end subroutine radius_prim_theta

! Added entire subroutine back in (NEW May 2025)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_both(x,z,ex,ez,r_in,r_out)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The code will only get here for the plasma-vacuum case.
! March 25 2022: For now, only the tri_type=-3 case is needed. Other cases are coded just in case.
		
		use triangularity
		use constant
		
		real (kind=dkind), intent(in) :: x, z
		real (kind=dkind), intent(out) :: ex, ez
		real (kind=dkind) :: r_in, r_out, thet, angle
		
		ex = x - rmajor
		ez = z
		
		if (ex==0.d0) then
			thet = pi/2.d0 * dsign(1.d0,ez)
			angle = thet
		else
			thet = datan2(ez,ex)
			angle = thet
		endif
		
		if((tri_type==13).or.(tri_type==-2).or.(tri_type==-3)) then
		
			! plasma - vacuum case
		
			if(angle<0.d0) angle = angle + 2.d0*pi
		
			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
				theta_points1, r_cscoef(1,1:theta_points1) )
		
			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )
			!print *, angle, ", ", r_out
		
		else
		
			print*, 'error in radius_both: tri_type=', tri_type
			pause
			stop
		
		endif
		
		continue
		
		
	end subroutine radius_both




!-----------------------------------------stuff that used to be in root-----------------------------------------


! Given a function fx defined on the interval from x1 -> x2, subdivide
! the interval into n equally spaced segments, and search for zero
! crossings of the function.  nb is input as the maximum number of
! roots sought, and is reset to the number of bracketing pairs
! xb1(1:nb), xb2(1:nb) that are found.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine zbrak(fx,x1,x2,m,xb1,xb2,nb)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dkind
  implicit none

  integer, intent(in) :: m
  integer, intent(inout) :: nb
  real (kind=dkind) :: x1,x2,fx
  external fx
  real (kind=dkind), dimension(1:nb), intent(out) :: xb1,xb2
  integer :: i,nbb
  real (kind=dkind) :: dx,fc,fp,xp,xc

  if(x2 <= x1) then
     print *, "[zbrak]: x1 = ",x1," and x2 = ",x2
     nb = 0
     stop
  end if

  nbb=0
  xp=x1
  dx=(x2-x1)/m ! Determine the spacing appropriate to the mesh
  fp=fx(xp)

  lp: do i=1,m ! Loop over all intervals
     xc=xp+dx
     fc=fx(xc)
     ! if a sign change occurs then record values for the bounds
     if(fc*fp.lt.0.) then 
        nbb=nbb+1
        xb1(nbb)=xp
        xb2(nbb)=xc
        if(nbb.eq.nb) exit lp
     endif
     ! Shift the x position and function value
     xp=xc
     fp=fc
  end do lp

  nb=nbb
  return
end subroutine zbrak

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function rtbis(func,x1,x2,xacc) result (soln)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dkind
!  use solver, only: b_plot
  implicit none

  integer, parameter :: jmax = 400
  real (kind=dkind) :: x1,x2,xacc,func
  real (kind=dkind) :: soln
  external func
  integer :: j
  real (kind=dkind) :: dx,f,fmid,xmid

  fmid = func(x2)
  f = func(x1)
  if(f*fmid .ge. 0.0d0) then
     print *, '[rtbis]: root must be bracketed'
     soln = x1
!	 call b_plot(x1,x2,1000)
     return
  end if

  if(f .lt. 0.0d0)then
     soln = x1
     dx = x2 - x1
  else
     soln = x2
     dx = x1 - x2
  endif

  do j = 1, jmax
     dx = 0.5d0*dx
     xmid = soln + dx
     fmid = func(xmid)
     if(fmid .le. 0.0d0) soln = xmid
     if(dabs(dx) .lt. xacc .or. fmid .eq. 0.0d0) return
  end do

  print *, '[rtbis]: too many bisections'

end function rtbis

!------------------------------------------------------------------

! Using a combination of Newton-Raphson and Bisection, find the root
! of a function between x1 and x2.  The root, returned as the value 
! soln, will be refined until its accuracy is known within +/-xacc.
! "funcd" is a user supplied subroutine which returns both the 
! value of the function and its first derivative.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function rtsafe(funcd,x1,x2,xacc,maxit) result (soln)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  integer :: maxit
  real (kind=dkind) :: x1,x2,xacc
  real (kind=dkind) :: soln
  external funcd
  integer :: j
  real (kind=dkind) :: df,dxx,dxold,f,fh,fl,temp,xh,xl

  call funcd(x1,fl,df)
  call funcd(x2,fh,df)
  if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
     print *, '[rtsafe]: root is not bracketed'
!	 call b_plot(x1,x2,1000)
     soln = 0.0 ! set it to something
     ! In fact since I will be looking for a density which 
     ! I know to be positive this will work fine for an error indicator
	 return
  end if

  if(fl.eq.0.)then
     soln=x1
     return
  else if(fh.eq.0.)then
     soln=x2
     return
  else if(fl.lt.0.)then ! Orient the search so that f(xl) < 0
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
  endif

  soln=.5*(x1+x2)  ! Initialize the guess for the root, 
  dxold=abs(x2-x1) ! the "stepsize before last", 
  dxx=dxold         ! and the last step.

  call funcd(soln,f,df) ! Initialize f and df for first guess

  lp: do j=1,maxit ! Loop over the allowed iterations
     ! If the Newton step is out of range or not decreasing 
     ! fast enough, use bisection 
     if(((soln-xh)*df-f)*((soln-xl)*df-f) .ge. 0.d0 &
          .or. abs(2.d0*f) .gt. abs(dxold*df)) then
        dxold=dxx
        dxx=0.5d0*(xh-xl)
        soln=xl+dxx
        if(xl .eq. soln) return
     else
        dxold=dxx
        dxx=f/df
        temp=soln
        soln=soln-dxx
        if(temp .eq. soln) return
     endif
     ! Return if the convergence criterion is met
     if(abs(dxx) .lt. xacc) return

     ! Otherwise, evaluate f and df at our new guess
     call funcd(soln,f,df)

     ! Shrink the bracket enclosing the root
     if(f.lt.0.) then
        xl=soln
     else
        xh=soln
     endif
  end do lp

  print *, '[rtsafe]: maximum number of iterations exceeded'
!  call b_plot(x1,x2,1000)

  return

end function rtsafe

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine secroot(fx,x1,x2,x,itemax,xtoll,ftoll,error)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	real(kind=dkind) :: x1,x2,x
	integer :: i
	integer :: itemax
	real(kind=dkind) :: xtoll, ftoll
	real(kind=dkind) :: f1,f2,fx
	real(kind=dkind) :: m, x3, f3
	integer :: error

!	external fx

	error = 0

	f1 = fx(x1)
	f2 = fx(x2)

	if(f1*f2>0.d0) then

!		pause 'error in secroot, root is not bracketed'
		error = 1
		return

	endif

	if(f2<f1) then
	! (f2>0 by hypothesis)

		x3 = x2
		x2 = x1
		x1 = x3

		f3 = f2
		f2 = f1
		f1 = f3

	endif




	sloop : do i=1,itemax

		m = (f2-f1)/(x2-x1)
		x3 = x1 - f1/m
		f3 = fx(x3)

		if(f3>0.d0) then

			x2 = x3
			f2 = f3

		else

			x1 = x3
			f1 = f3

		endif

		if( (abs(f2-f1)<ftoll).or.(abs(x2-x1)<xtoll).or.  &
			 (abs(f2)<ftoll).or.(abs(f1)<ftoll) ) then
		! converged!

			exit sloop

		endif

	enddo sloop

	if(abs(f2)<abs(f1)) then

		x = x2

	else

		x = x1

	endif

end subroutine secroot

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      function derfc (x) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!!$***BEGIN PROLOGUE  DERFC
!!$***PURPOSE  Compute the complementary error function.
!!$***LIBRARY   SLATEC (FNLIB)
!!$***CATEGORY  C8A, L5A1E
!!$***TYPE      DOUBLE PRECISION (ERFC-S, DERFC-D)
!!$***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB,
!!$             SPECIAL FUNCTIONS
!!$***AUTHOR  Fullerton, W., (LANL)
!!$***DESCRIPTION
!!$
!!$ DERFC(X) calculates the double precision complementary error function
!!$ for double precision argument X.
!!$
!!$ Series for ERF        on the interval  0.          to  1.00000E+00
!!$                                        with weighted Error   1.28E-32
!!$                                         log weighted Error  31.89
!!$                               significant figures required  31.05
!!$                                    decimal places required  32.55
!!$
!!$ Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00
!!$                                        with weighted Error   2.67E-32
!!$                                         log weighted Error  31.57
!!$                               significant figures required  30.31
!!$                                    decimal places required  32.42
!!$
!!$ Series for ERFC       on the interval  0.          to  2.50000E-01
!!$                                        with weighted error   1.53E-31
!!$                                         log weighted error  30.82
!!$                               significant figures required  29.47
!!$                                    decimal places required  31.70
!!$
!!$***REFERENCES  (NONE)
!!$***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!!$***REVISION HISTORY  (YYMMDD)
!!$   770701  DATE WRITTEN
!!$   890531  Changed all specific intrinsics to generic.  (WRB)
!!$   890531  REVISION DATE from Version 3.2
!!$   891214  Prologue converted to Version 4.0 format.  (BAB)
!!$   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!!$   920618  Removed space from variable names.  (RWC, WRB)
!!$***END PROLOGUE  DERFC

	use constant, only : dkind

	implicit none

      real(kind=dkind) :: X, ERFCS(21), ERFCCS(59), ERC2CS(49), SQEPS,  &
       SQRTPI, XMAX, TXMAX, XSML, Y, answer	!DCSEVL, D1MACH,
	real (kind=dkind) :: eta
	integer ::  nterf, nterfc, nterc2
      LOGICAL FIRST
      SAVE ERFCS, ERC2CS, ERFCCS, SQRTPI, NTERF,  &
      NTERFC, NTERC2, XSML, XMAX, SQEPS, FIRST

      DATA ERFCS(  1) / -.49046121234691808039984544033376D-1     /
      DATA ERFCS(  2) / -.14226120510371364237824741899631D+0     /
      DATA ERFCS(  3) / +.10035582187599795575754676712933D-1     /
      DATA ERFCS(  4) / -.5768764699767484765082702550917D-3     /
      DATA ERFCS(  5) / +.27419931252196061034422160791471D-4     /
      DATA ERFCS(  6) / -.11043175507344507604135381295905D-5     /
      DATA ERFCS(  7) / +.38488755420345036949961311498174D-7     /
      DATA ERFCS(  8) / -.11808582533875466969631751801581D-8     /
      DATA ERFCS(  9) / +.32334215826050909646402930953354D-10    /
      DATA ERFCS( 10) / -.79910159470045487581607374708595D-12    /
      DATA ERFCS( 11) / +.17990725113961455611967245486634D-13    /
      DATA ERFCS( 12) / -.37186354878186926382316828209493D-15    /
      DATA ERFCS( 13) / +.71035990037142529711689908394666D-17    /
      DATA ERFCS( 14) / -.12612455119155225832495424853333D-18    /
      DATA ERFCS( 15) / +.20916406941769294369170500266666D-20    /
      DATA ERFCS( 16) / -.32539731029314072982364160000000D-22    /
      DATA ERFCS( 17) / +.47668672097976748332373333333333D-24    /
      DATA ERFCS( 18) / -.65980120782851343155199999999999D-26    /
      DATA ERFCS( 19) / +.86550114699637626197333333333333D-28    /
      DATA ERFCS( 20) / -.10788925177498064213333333333333D-29    /
      DATA ERFCS( 21) / +.12811883993017002666666666666666D-31    /
      DATA ERC2CS(  1) / -.6960134660230950112739150826197D-1      /
      DATA ERC2CS(  2) / -.4110133936262089348982212084666D-1      /
      DATA ERC2CS(  3) / +.3914495866689626881561143705244D-2      /
      DATA ERC2CS(  4) / -.4906395650548979161280935450774D-3      /
      DATA ERC2CS(  5) / +.7157479001377036380760894141825D-4      /
      DATA ERC2CS(  6) / -.1153071634131232833808232847912D-4      /
      DATA ERC2CS(  7) / +.1994670590201997635052314867709D-5      /
      DATA ERC2CS(  8) / -.3642666471599222873936118430711D-6      /
      DATA ERC2CS(  9) / +.6944372610005012589931277214633D-7      /
      DATA ERC2CS( 10) / -.137122090210436601953460514121D-7      /
      DATA ERC2CS( 11) / +.2788389661007137131963860348087D-8      /
      DATA ERC2CS( 12) / -.5814164724331161551864791050316D-9      /
      DATA ERC2CS( 13) / +.1238920491752753181180168817950D-9      /
      DATA ERC2CS( 14) / -.2690639145306743432390424937889D-10     /
      DATA ERC2CS( 15) / +.5942614350847910982444709683840D-11     /
      DATA ERC2CS( 16) / -.1332386735758119579287754420570D-11     /
      DATA ERC2CS( 17) / +.3028046806177132017173697243304D-12     /
      DATA ERC2CS( 18) / -.6966648814941032588795867588954D-13     /
      DATA ERC2CS( 19) / +.1620854541053922969812893227628D-13     /
      DATA ERC2CS( 20) / -.3809934465250491999876913057729D-14     /
      DATA ERC2CS( 21) / +.9040487815978831149368971012975D-15     /
      DATA ERC2CS( 22) / -.2164006195089607347809812047003D-15     /
      DATA ERC2CS( 23) / +.5222102233995854984607980244172D-16     /
      DATA ERC2CS( 24) / -.1269729602364555336372415527780D-16     /
      DATA ERC2CS( 25) / +.3109145504276197583836227412951D-17     /
      DATA ERC2CS( 26) / -.7663762920320385524009566714811D-18     /
      DATA ERC2CS( 27) / +.1900819251362745202536929733290D-18     /
      DATA ERC2CS( 28) / -.4742207279069039545225655999965D-19     /
      DATA ERC2CS( 29) / +.1189649200076528382880683078451D-19     /
      DATA ERC2CS( 30) / -.3000035590325780256845271313066D-20     /
      DATA ERC2CS( 31) / +.7602993453043246173019385277098D-21     /
      DATA ERC2CS( 32) / -.1935909447606872881569811049130D-21     /
      DATA ERC2CS( 33) / +.4951399124773337881000042386773D-22     /
      DATA ERC2CS( 34) / -.1271807481336371879608621989888D-22     /
      DATA ERC2CS( 35) / +.3280049600469513043315841652053D-23     /
      DATA ERC2CS( 36) / -.8492320176822896568924792422399D-24     /
      DATA ERC2CS( 37) / +.2206917892807560223519879987199D-24     /
      DATA ERC2CS( 38) / -.5755617245696528498312819507199D-25     /
      DATA ERC2CS( 39) / +.1506191533639234250354144051199D-25     /
      DATA ERC2CS( 40) / -.395450295901879695310428569599D-26     /
      DATA ERC2CS( 41) / +.1041529704151500979984645051733D-26     /
      DATA ERC2CS( 42) / -.275148779527876507945017890133D-27     /
      DATA ERC2CS( 43) / +.7290058205497557408997703680000D-28     /
      DATA ERC2CS( 44) / -.193693964591594780407750109866D-28     /
      DATA ERC2CS( 45) / +.5160357112051487298370054826666D-29     /
      DATA ERC2CS( 46) / -.137841932219309409938964480000D-29     /
      DATA ERC2CS( 47) / +.3691326793107069042251093333333D-30     /
      DATA ERC2CS( 48) / -.990938959062436542065322666666D-31     /
      DATA ERC2CS( 49) / +.2666491705195388413323946666666D-31     /
      DATA ERFCCS(  1) / +.715179310202924774503697709496D-1        /
      DATA ERFCCS(  2) / -.265324343376067157558893386681D-1        /
      DATA ERFCCS(  3) / +.171115397792085588332699194606D-2        /
      DATA ERFCCS(  4) / -.163751663458517884163746404749D-3        /
      DATA ERFCCS(  5) / +.198712935005520364995974806758D-4        /
      DATA ERFCCS(  6) / -.284371241276655508750175183152D-5        /
      DATA ERFCCS(  7) / +.460616130896313036969379968464D-6        /
      DATA ERFCCS(  8) / -.822775302587920842057766536366D-7        /
      DATA ERFCCS(  9) / +.159214187277090112989358340826D-7        /
      DATA ERFCCS( 10) / -.329507136225284321486631665072D-8        /
      DATA ERFCCS( 11) / +.722343976040055546581261153890D-9        /
      DATA ERFCCS( 12) / -.166485581339872959344695966886D-9        /
      DATA ERFCCS( 13) / +.401039258823766482077671768814D-10       /
      DATA ERFCCS( 14) / -.100481621442573113272170176283D-10       /
      DATA ERFCCS( 15) / +.260827591330033380859341009439D-11       /
      DATA ERFCCS( 16) / -.699111056040402486557697812476D-12       /
      DATA ERFCCS( 17) / +.192949233326170708624205749803D-12       /
      DATA ERFCCS( 18) / -.547013118875433106490125085271D-13       /
      DATA ERFCCS( 19) / +.158966330976269744839084032762D-13       /
      DATA ERFCCS( 20) / -.472689398019755483920369584290D-14       /
      DATA ERFCCS( 21) / +.143587337678498478672873997840D-14       /
      DATA ERFCCS( 22) / -.444951056181735839417250062829D-15       /
      DATA ERFCCS( 23) / +.140481088476823343737305537466D-15       /
      DATA ERFCCS( 24) / -.451381838776421089625963281623D-16       /
      DATA ERFCCS( 25) / +.147452154104513307787018713262D-16       /
      DATA ERFCCS( 26) / -.489262140694577615436841552532D-17       /
      DATA ERFCCS( 27) / +.164761214141064673895301522827D-17       /
      DATA ERFCCS( 28) / -.562681717632940809299928521323D-18       /
      DATA ERFCCS( 29) / +.194744338223207851429197867821D-18       /
      DATA ERFCCS( 30) / -.682630564294842072956664144723D-19       /
      DATA ERFCCS( 31) / +.242198888729864924018301125438D-19       /
      DATA ERFCCS( 32) / -.869341413350307042563800861857D-20       /
      DATA ERFCCS( 33) / +.315518034622808557122363401262D-20       /
      DATA ERFCCS( 34) / -.115737232404960874261239486742D-20       /
      DATA ERFCCS( 35) / +.428894716160565394623737097442D-21       /
      DATA ERFCCS( 36) / -.160503074205761685005737770964D-21       /
      DATA ERFCCS( 37) / +.606329875745380264495069923027D-22       /
      DATA ERFCCS( 38) / -.231140425169795849098840801367D-22       /
      DATA ERFCCS( 39) / +.888877854066188552554702955697D-23       /
      DATA ERFCCS( 40) / -.344726057665137652230718495566D-23       /
      DATA ERFCCS( 41) / +.134786546020696506827582774181D-23       /
      DATA ERFCCS( 42) / -.531179407112502173645873201807D-24       /
      DATA ERFCCS( 43) / +.210934105861978316828954734537D-24       /
      DATA ERFCCS( 44) / -.843836558792378911598133256738D-25       /
      DATA ERFCCS( 45) / +.339998252494520890627359576337D-25       /
      DATA ERFCCS( 46) / -.137945238807324209002238377110D-25       /
      DATA ERFCCS( 47) / +.563449031183325261513392634811D-26       /
      DATA ERFCCS( 48) / -.231649043447706544823427752700D-26       /
      DATA ERFCCS( 49) / +.958446284460181015263158381226D-27       /
      DATA ERFCCS( 50) / -.399072288033010972624224850193D-27       /
      DATA ERFCCS( 51) / +.167212922594447736017228709669D-27       /
      DATA ERFCCS( 52) / -.704599152276601385638803782587D-28       /
      DATA ERFCCS( 53) / +.297976840286420635412357989444D-28       /
      DATA ERFCCS( 54) / -.126252246646061929722422632994D-28       /
      DATA ERFCCS( 55) / +.539543870454248793985299653154D-29       /
      DATA ERFCCS( 56) / -.238099288253145918675346190062D-29       /
      DATA ERFCCS( 57) / +.109905283010276157359726683750D-29       /
      DATA ERFCCS( 58) / -.486771374164496572732518677435D-30       /
      DATA ERFCCS( 59) / +.152587726411035756763200828211D-30       /
      DATA SQRTPI / 1.77245385090551602729816748334115D0 /
      DATA FIRST /.TRUE./
!!$***FIRST EXECUTABLE STATEMENT  DERFC
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NTERF = INITDS (ERFCS, 21, ETA)
         NTERFC = INITDS (ERFCCS, 59, ETA)
         NTERC2 = INITDS (ERC2CS, 49, ETA)
!!$
         XSML = -SQRT(-LOG(SQRTPI*D1MACH(3)))
         TXMAX = SQRT(-LOG(SQRTPI*D1MACH(1)))
         XMAX = TXMAX - 0.5D0*LOG(TXMAX)/TXMAX - 0.01D0
         SQEPS = SQRT(2.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
!!$
      IF (X.GT.XSML) GO TO 20
!!$
!!$ ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
!!$
      answer = 2.0D0
      RETURN
!!$
 20   IF (X.GT.XMAX) GO TO 40
      Y = ABS(X)
      IF (Y.GT.1.0D0) GO TO 30
!!$
!!$ ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
!!$
      IF (Y.LT.SQEPS) answer = 1.0D0 - 2.0D0*X/SQRTPI
      IF (Y.GE.SQEPS) answer = 1.0D0 - X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,  &
       ERFCS, NTERF))
      RETURN
!!$
!!$ ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
!!$
 30   Y = Y*Y
      IF (Y.LE.4.D0) answer = EXP(-Y)/ABS(X) * (0.5D0 + DCSEVL (  &
       (8.D0/Y-5.D0)/3.D0, ERC2CS, NTERC2) )
      IF (Y.GT.4.D0) answer = EXP(-Y)/ABS(X) * (0.5D0 + DCSEVL (  &
       8.D0/Y-1.D0, ERFCCS, NTERFC) )
      IF (X.LT.0.D0) answer = 2.0D0 - answer
      RETURN
!!$
 40  print*, 'SLATEC', 'DERFC', 'X SO BIG ERFC UNDERFLOWS', 1, 1
      answer = 0.D0
      return
!!$
      end function derfc



 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     function initds (OS, NOS, ETA) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$***BEGIN PROLOGUE  INITDS
!!$***PURPOSE  Determine the number of terms needed in an orthogonal
!!$            polynomial series so that it meets a specified accuracy.
!!$***LIBRARY   SLATEC (FNLIB)
!!$***CATEGORY  C3A2
!!$***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
!!$***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!!$             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!!$***AUTHOR  Fullerton, W., (LANL)
!!$***DESCRIPTION
!!$
!!$  Initialize the orthogonal series, represented by the array OS, so
!!$  that INITDS is the number of terms needed to insure the error is no
!!$  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!!$  machine precision.
!!$
!!$             Input Arguments --
!!$   OS     double precision array of NOS coefficients in an orthogonal
!!$          series.
!!$   NOS    number of coefficients in OS.
!!$   ETA    single precision scalar containing requested accuracy of
!!$          series.
!!$
!!$***REFERENCES  (NONE)
!!$***ROUTINES CALLED  XERMSG
!!$***REVISION HISTORY  (YYMMDD)
!!$   770601  DATE WRITTEN
!!$   890531  Changed all specific intrinsics to generic.  (WRB)
!!$   890831  Modified array declarations.  (WRB)
!!$   891115  Modified error message.  (WRB)
!!$   891115  REVISION DATE from Version 3.2
!!$   891214  Prologue converted to Version 4.0 format.  (BAB)
!!$   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!!$***END PROLOGUE  INITDS

	use constant, only : dkind

	implicit none

    real(kind=dkind) OS(*), eta, err
	real(kind=dkind) :: answer
	integer :: nos
	integer :: i, ii

!!$***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1)  print*, 'SLATEC', 'INITDS',  &
        'Number of coefficients is less than 1', 2, 1
!!$
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
!!$
   20 IF (I .EQ. NOS) print*, 'SLATEC', 'INITDS',  &
        'Chebyshev series too short for specified accuracy', 1, 1
      answer = I
!!$
      RETURN
      END function initds




 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      function DCSEVL (X, CS, N) result(answer)
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$***BEGIN PROLOGUE  DCSEVL
!!$***PURPOSE  Evaluate a Chebyshev series.
!!$***LIBRARY   SLATEC (FNLIB)
!!$***CATEGORY  C3A2
!!$***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!!$***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!!$***AUTHOR  Fullerton, W., (LANL)
!!$***DESCRIPTION
!!$
!!$  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!!$  a method presented in the paper by Broucke referenced below.
!!$
!!$       Input Arguments --
!!$  X    value at which the series is to be evaluated.
!!$  CS   array of N terms of a Chebyshev series.  In evaluating
!!$       CS, only half the first coefficient is summed.
!!$  N    number of terms in array CS.
!!$
!!$***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!!$                 Chebyshev series, Algorithm 446, Communications of
!!$                 the A.C.M. 16, (1973) pp. 254-256.
!!$               L. Fox and I. B. Parker, Chebyshev Polynomials in
!!$                 Numerical Analysis, Oxford University Press, 1968,
!!$                 page 56.
!!$***ROUTINES CALLED  D1MACH, XERMSG
!!$***REVISION HISTORY  (YYMMDD)
!!$   770401  DATE WRITTEN
!!$   890831  Modified array declarations.  (WRB)
!!$   890831  REVISION DATE from Version 3.2
!!$   891214  Prologue converted to Version 4.0 format.  (BAB)
!!$   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!!$   900329  Prologued revised extensively and code rewritten to allow
!!$           X to be slightly outside interval (-1,+1).  (WRB)
!!$   920501  Reformatted the REFERENCES section.  (WRB)
!!$***END PROLOGUE  DCSEVL

	use constant, only : dkind

	implicit none

    real(kind=dkind) :: B0, B1, B2, CS(*), ONEPL, TWOX, X	!, D1MACH
	real(kind=dkind) :: answer
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./

	  integer :: n, i, ni

!!$***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) print*, 'SLATEC', 'DCSEVL',  'NUMBER OF TERMS .LE. 0', 2, 2
      IF (N .GT. 1000) print*, 'SLATEC', 'DCSEVL', 'NUMBER OF TERMS .GT. 1000', 3, 2
      IF (ABS(X) .GT. ONEPL) print*, 'SLATEC', 'DCSEVL', 'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1
!!$
      B1 = 0.0D0
      B0 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!!$
      answer = 0.5D0*(B0-B2)
!!$
      RETURN

      end function dcsevl



  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     function D1MACH(I)
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

      INTEGER I
	  real(kind=dkind) :: D1MACH, DMACH(5)

!!$
!!$  DOUBLE-PRECISION MACHINE CONSTANTS
!!$  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!!$  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!!$  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!!$  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!!$  D1MACH( 5) = LOG10(B)
!!$
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!!$  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
!!$  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
!!$  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
!!$  MANY MACHINES YET.
!!$  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!!$  ON THE NEXT LINE
      DATA SC/0/
!!$  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!!$  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!!$          mail netlib@research.bell-labs.com
!!$          send old1mach from blas
!!$  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!!$
!!$     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!!$      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!!$      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!!$      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!!$      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!!$      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
!!$
!!$     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!!$     32-BIT INTEGERS.
!!$      DATA SMALL(1),SMALL(2) /    8388608,           0 /
!!$      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!!$      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!!$      DATA DIVER(1),DIVER(2) /  620756992,           0 /
!!$      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
!!$
!!$     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!!$      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!!$      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!!$      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!!$      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!!$      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
!!$
!!$     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532.AND. SMALL(2) .EQ. -448790528) THEN
!!$           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532 .AND. SMALL(1) .EQ. -448790528) THEN
!!$           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935.AND. SMALL(2) .EQ. 10752) THEN
!!$               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943.AND. SMALL(2) .EQ. 704643072) THEN
!!$               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684.AND. SMALL(2) .EQ. -448790528) THEN
!!$           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074.AND. SMALL(2) .EQ. 58688) THEN
!!$           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
!!$                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
!!$    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/' appropriate for your machine.')

!!$* /* Standard C source for D1MACH -- remove the * in column 1 */
!!$*#include <stdio.h>
!!$*#include <float.h>
!!$*#include <math.h>
!!$*double d1mach_(long *i)
!!$*{
!!$*	switch(*i){
!!$*	  case 1: return DBL_MIN;
!!$*	  case 2: return DBL_MAX;
!!$*	  case 3: return DBL_EPSILON/FLT_RADIX;
!!$*	  case 4: return DBL_EPSILON;
!!$*	  case 5: return log10((double)FLT_RADIX);
!!$*	  }
!!$*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
!!$*	exit(1); return 0; /* some compilers demand return values */
!!$*}
      end function D1MACH

      SUBROUTINE I1MCRY(A, A1, B, C, D)
!!$**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END SUBROUTINE I1MCRY



! Test out the solver for the bernoulli function by finding num 
! different solutions
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine b_plot(min,max,num)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent(in) :: min,max
	integer, intent(in) :: num
	real (kind=dkind) :: rho,drho,b,db
	integer :: i

	drho = (max - min)/(num - 1)

	open (unit=13,file='b_plot.dat',status='unknown',action='write')

	open (unit=11,file='db_plot.dat',status='unknown',action='write')

	rho = min
	do i = 1,num
	call newt_rap(rho,b,db)

	write (13,*) rho,b
	write (11,*) rho,db

	rho = rho + drho
	end do

	close(13)
	close(11)

end subroutine b_plot



! Test out the solver for the bernoulli function by finding num 
! different solutions
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine b_plot_gauss(min,max,num)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent(in) :: min,max
	integer, intent(in) :: num
	real (kind=dkind) :: rho,drho,b,db,b_g,db_g
	integer :: i

	drho = (max - min)/(num - 1)

	open (unit=13,file='b_plot_gauss.dat',status='unknown',action='write')

	open (unit=11,file='db_plot_gauss.dat',status='unknown',action='write')

	rho = min
	do i = 1,num
	call newt_rap(rho,b,db)
	call newt_rap_gauss(rho,b_g,db_g)

	write (13,*) rho,b,b_g
	write (11,*) rho,db,db_g

	rho = rho + drho
	end do

	close(13)
	close(11)

end subroutine b_plot_gauss

!-----------------------------------------end of stuff that used to be in root-----------------------------------------


!-----------------------------------------stuff that used to be in inter_grid-----------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function psi_sol(r) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, rmajor
	use magnetic, only : psi_e
	use triangularity, only : theta_temp, psifun

	implicit none

	real(kind=dkind) :: x,z,r
	real(kind=dkind) :: answer

	x = rmajor + r*cos(theta_temp)
	z = r*sin(theta_temp)


	answer = psifun(x,z) - psi_e

	continue

	return

end function psi_sol



!-----------------------------------------end of stuff that used to be in inter_grid-----------------------------------------



end module solver

