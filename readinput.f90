!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine readinput(m)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use array_dimensions
	use constant
	use flow
	use magnetic
	use p_d_profile
	use exp_data
	use triangularity
	use solver

	implicit none

	integer :: m
	real(kind=dkind) :: Fc_o_Fv		!Fcenter/Fvacuum
	real(kind=dkind) :: de_o_dc		!Dedge/Dcenter
	real(kind=dkind) ::	pe_o_pc		!Pedge/Pcenter
	real(kind=dkind) ::	pe_o_pc_e		!Pedge/Pcenter
	real(kind=dkind) ::	pe_o_pc_i		!Pedge/Pcenter

	include "inp_out.f90"

	open(unit=5,file='FLOW2_inputfile.dat',status='old')

	read(5,input_which)

	read(5,input_triangularity)
	read(5,input_constants)
	read(5,input_flow)
	read(5,input_magnetic)
	read(5,input_p_d_profile)
	read(5,input_numerical)
	read(5,input_solver)
	read(5,input_gravity)
	read(5,input_output)

	close(5)

	continue

	m	= 1 + n/2
	Broot_true = Broot

	fvacuum = b_phi_zero*rmajor
	fcenter	= Fc_o_Fv*fvacuum

	dedge = dcenter * de_o_dc

	if(eq_type>=10) then
		d_TF_center = dcenter/mass
		d_TF_edge = dedge/mass
	endif

	if(F_opt==0) then
		bpol0_fix = b_phi_zero ! change later, do not use the same symbol for different things
	endif

	pcenter = beta_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag
	pedge = pcenter * pe_o_pc

	p_e_center = beta_e_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag
	p_e_edge = p_e_center * pe_o_pc_e

	p_i_center = beta_i_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag
	p_i_edge = p_i_center * pe_o_pc_i

	if(input_EQDSK) then

		! all that follows is not needed
		continue

	else

		if(numerical) then

				continue

		else

				numerical_n =		.false.
				numerical_p_iso =	.false.
				numerical_F =		.false.
				numerical_omega =	.false.
				numerical_mtheta =	.false.
				numerical_psiprim = .false.
				numerical_p_e 	=	.false.
				numerical_p_i 	=	.false.
				numerical_psi_diff = .false.

		endif

		! Commented this out (NEW May 2025)
		!! Superceded by code below
		!!if((bc_type==7).or.(bc_type==17)) then
		!!	LCFS = 1 ! big_Psi
		!!elseif((bc_type==-7).or.(bc_type==-17)) then
		!!	LCFS = -1 ! psi
		!!	bc_type = abs(bc_type)
		!!endif

!		if(((bc_type==7).or.(bc_type==8).or.(bc_type==17)).or.(((bc_type==21).or.  &
		if(((modulo(abs(bc_type),10)==7).or.(bc_type==8)).or.(((bc_type==21).or.  &
			(bc_type==23).or.(bc_type==24)).and.(numerical_psi_diff))) call read_numerical_bc

		! This sets up which surface is used to define the boundary
		if(modulo(abs(bc_type),10)==7) then
			if(bc_type>0) then
				LCFS = 1
			else
				LCFS = -1
				bc_type = abs(bc_type)
			endif
		endif

		if((bc_type/=7).and.(bc_type/=8)) then ! need to check if (bc_type==8) still does anything
			if(psi_diff_option>2) then
				psi_diff_option = 2
			endif
		endif

		continue

		! some controls on the input
		if ((eq_type==1).or.(eq_type==3).or.(eq_type==10)) then
			continue
		else
			print*,'error in eq_type :  ',eq_type
		endif

		if( (bc_type==7).and.(tri_type==13)) then
			psic_13 = psic - psi_e
		endif

		! Added tri_type==-3 to if statement (NEW May 2025)
		if( ((bc_type==17).or.(bc_type==7)).and.((tri_type==13).or.(tri_type==-1)  &
			.or.(tri_type==-2).or.(tri_type==-3).or.(tri_type==-4))) then	! WARNING!!!!!
			continue
		elseif( (bc_type>=3).and.(tri_type/=0).and.(tri_type/=8).and.(tri_type/=9)  &
					.and.(tri_type/=18).and.(tri_type/=88) ) then
			print*, 'warning, option not yet implemented for'
			print*, 'bc_type =', bc_type, 'and tri_type =', tri_type
			print*, 'bc_type changed to 1'
			bc_type = 1
		endif

		if(grid_type==3) call read_grid

	endif

	psi_max = max(psic, big_Psic)

	bc_setup_option = sign(1,bc_type)
	bc_type = abs(bc_type)

	  if (numerical) then

			if(numerical_opt==1) then

				call read_numerical

			elseif(numerical_opt>=2) then

				call read_numerical_super

			else

				pause 'error in numerical_opt'
				stop

			endif

			continue

	  endif


	if((bc_type==21).or.(bc_type==23).or.(bc_type==24).or.(bc_type==34)) then

		if(numerical_psi_diff) then

			continue

		else

			if(omega_0/=omegaofpsi(0.d0)) then
				print*, 'WARNING: for bc_type = ', bc_type
				print*, 'omega_0 is different from omegaofpsi(0)'
				pause
			endif

		endif

	endif

	  mach_theta_max_initial = mach_theta_max

	! the eqdsk file contains several data that are duplicated above, and which will be changed
	! the previous part of the reading process has been preserved to leave control over
	! the details of the FLOW input not contained in eqdsk (TO BE EDITED LATER)
	if(input_EQDSK) then
		call readeqdsk(trim(EQDSK_file))
	endif

	! last check, verifies whether there is any flow in the equilibrium

	if(enforce_static==0) then

		call check_flow

	elseif(enforce_static==1) then

		static_equi = .true.

	elseif(enforce_static==2) then

		static_equi = .false.

	endif

	! determine the number of q values required for v_norm output
	do

		if(q_vnorm(n_q_vnorm+1)==99.d0) exit

		n_q_vnorm = n_q_vnorm + 1

	enddo

	continue

end subroutine readinput

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_numerical
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use magnetic
	use exp_data
	use interpolating_functions
	use p_d_profile
	use flow
	use pseudo_IMSL, only : DBSNAK, DBSINT, dbsval
	use solver, only : p_e_ofpsi, p_i_ofpsi, D_TF_num

	implicit none

	integer i
	integer :: d_dim
	integer :: dim_default = 40

	if(numerical_n) then

		if(eq_type>=10) then
			D_TF_num = .true.
		endif

		i = 0

		open(33,file='n.dat', status='old', action='read')

		do

			i = i+1
			read(33,*,end=19) d_data(i,1), d_data(i,2)

		enddo

	19	close(33)

		d_dim = i-1
		ibreak_d = d_dim + d_ord

		allocate(d_cscoef(d_dim))

		!this sets up the knots
		call DBSNAK(d_dim, d_data(1:d_dim,1),  &
						d_ord,d_data(1:ibreak_d,3))

		call DBSINT(d_dim, d_data(1:d_dim,1),  &
			 d_data(1:d_dim,2), d_ord,  &
			 d_data(1:ibreak_d,3),  &
			 d_cscoef(1:d_dim))

	endif

!!$	if(input_EQDSK) return
!!$ March 6 2008: this return has been removed to allow for other numerical inputs LG

!----------------------------------------------------------------------------

	if(numerical_p_iso) then

	i = 0

	open(33,file='p_iso.dat', status='old', action='read')

	do 

		i = i+1
		read(33,*,end=66) p_iso_data(i,1), p_iso_data(i,2)

	enddo

66	close(33)

	data_dim = i-1
	ibreak_p = data_dim + p_iso_ord

	allocate(p_iso_cscoef(data_dim))

	!this sets up the knots
	call DBSNAK(data_dim, p_iso_data(1:data_dim,1),  &
					p_iso_ord,p_iso_data(1:ibreak_p,3))

	call DBSINT(data_dim, p_iso_data(1:data_dim,1),  &
		 p_iso_data(1:data_dim,2), p_iso_ord,  &
		 p_iso_data(1:ibreak_p,3),  &
		 p_iso_cscoef(1:data_dim))

	endif

!----------------------------------------------------------------------------

	if(numerical_p_e) then

	i = 0

	open(33,file='p_e.dat',  &
				status='old',action='read')

	do 

		i = i+1
		read(33,*,end=67) p_e_data(i,1), p_e_data(i,2)

	enddo

67	close(33)

	data_dim = i-1
	ibreak_p_e = data_dim + p_e_ord

	allocate(p_e_cscoef(data_dim))

	!this sets up the knots
	call DBSNAK(data_dim, p_e_data(1:data_dim,1),  &
					p_e_ord,p_e_data(1:ibreak_p_e,3))

	call DBSINT(data_dim, p_e_data(1:data_dim,1),  &
		 p_e_data(1:data_dim,2), p_e_ord,  &
		 p_e_data(1:ibreak_p_e,3),  &
		 p_e_cscoef(1:data_dim))

	endif

!----------------------------------------------------------------------------

	if(numerical_p_i) then

	i = 0

	open(33,file='p_i.dat',  &
				status='old',action='read')

	do 

		i = i+1
		read(33,*,end=68) p_i_data(i,1), p_i_data(i,2)

	enddo

68	close(33)

	data_dim = i-1
	ibreak_p_i = data_dim + p_i_ord

	allocate(p_i_cscoef(data_dim))

	!this sets up the knots
	call DBSNAK(data_dim, p_i_data(1:data_dim,1),  &
					p_i_ord,p_i_data(1:ibreak_p_i,3))

	call DBSINT(data_dim, p_i_data(1:data_dim,1),  &
		 p_i_data(1:data_dim,2), p_i_ord,  &
		 p_i_data(1:ibreak_p_i,3),  &
		 p_i_cscoef(1:data_dim))

	endif

!----------------------------------------------------------------------------

	if(numerical_F) then

	open(33,file='b0.dat', status='old', action='read')

	i = 0

	do
		i = i+1
		read(33,*,end=69) F_data(i,1),F_data(i,2)
		F_data(i,2) = F_data(i,2)/rmajor
	enddo

69	close(33)

	data_dim = i-1
	ibreak_B0 = data_dim + F_ord

	allocate(F_cscoef(data_dim))

	!this sets up the knots
	call DBSNAK(data_dim, F_data(1:data_dim,1),  &
					F_ord,F_data(1:ibreak_B0,3))

	call DBSINT(data_dim, F_data(1:data_dim,1),  &
		 F_data(1:data_dim,2), F_ord,  &
		 F_data(1:ibreak_B0,3),  &
		 F_cscoef(1:data_dim))

	endif

!----------------------------------------------------------------------------

	if(numerical_omega) then

		if(omega_option==1) then

			numerical_mphi = .true.
			numerical_omega = .false.

		endif

		open(33,file='mphi.dat', status='old', action='read')

		i = 0

		do

			i = i+1

			read(33,*,end=98) mphi_data(i,1),mphi_data(i,2)

		enddo

		98	close(33)

		data_dim = i-1
		ibreak_fi = data_dim + mphi_ord

		allocate(mphi_cscoef(data_dim))

		!this sets up the knots
		call DBSNAK(data_dim, mphi_data(1:data_dim,1),  &
						mphi_ord,mphi_data(1:ibreak_fi,3))

		call DBSINT(data_dim, mphi_data(1:data_dim,1),  &
			 mphi_data(1:data_dim,2), mphi_ord,  &
			 mphi_data(1:ibreak_fi,3),  &
			 mphi_cscoef(1:data_dim))

	endif

!----------------------------------------------------------------------------

	if(numerical_mtheta) then

		open(33,file='mtheta.dat', status='old', action='read')

		i = 0

		do

			i = i+1

			read(33,*,end=109) mtheta_data(i,1),mtheta_data(i,2)

		enddo

		109	close(33)

		data_dim = i-1
		ibreak_th = data_dim + mtheta_ord

		allocate(mtheta_cscoef(data_dim))

		! interpolation setup has been moved to a separate function
		call interp_setup(data_dim,mtheta_ord, &
			mtheta_data(1:data_dim,1),mtheta_data(1:data_dim,2), &
			mtheta_data(1:ibreak_th,3),mtheta_cscoef(1:data_dim))

		if((Broot==0).or.(Broot==5)) then
		! set maximum value of Mach_theta

			mach_theta_num = 0.d0

			do i=1,data_dim

				if(abs(mtheta_data(i,2))>mach_theta_num) then

					mach_theta_num = mtheta_data(i,2)

				endif

			enddo

			mach_theta_max = mach_theta_num

		else

			! any two numbers will do, as long as they are equal
			mach_theta_max = 1.d-1
			mach_theta_num = 1.d-1

		endif

	endif

!----------------------------------------------------------------------------

	if((numerical_psiprim).and.(bc_type==5)) then

		i = 0

		open(33,file='psiprim.dat', status='old', action='read')

		do

			i = i+1
			read(33,*,end=199) psiprim_data(i,1), psiprim_data(i,2)

		enddo

	199	close(33)

		psiprim_dim = i-1
		ibreak_pp = psiprim_dim + psiprim_ord

		allocate(psiprim_cscoef(psiprim_dim))

		!this sets up the knots
		call DBSNAK(psiprim_dim, psiprim_data(1:psiprim_dim,1),  &
						psiprim_ord,psiprim_data(1:ibreak_pp,3))

		call DBSINT(psiprim_dim, psiprim_data(1:psiprim_dim,1),  &
			 psiprim_data(1:psiprim_dim,2), psiprim_ord,  &
			 psiprim_data(1:ibreak_pp,3),  &
			 psiprim_cscoef(1:psiprim_dim))

	endif

!----------------------------------------------------------------------------

! for consistency, if either p_i or p_e is numerical, p_iso also needs to be numerical

	if(numerical_p_i .or. numerical_p_e) then

		if(numerical_p_iso) then

			continue

		else

			numerical_p_iso = .true.

			do i = 1, dim_default

				p_iso_data(i,1) = (i-1.d0)/(dim_default-1.d0)
				p_iso_data(i,2) = p_e_ofpsi(p_iso_data(i,1)*psic) + p_i_ofpsi(p_iso_data(i,1)*psic)

			enddo

			data_dim = i-1
			ibreak_p = data_dim + p_iso_ord

			allocate(p_iso_cscoef(data_dim))

			!this sets up the knots
			call DBSNAK(data_dim, p_iso_data(1:data_dim,1),  &
							p_iso_ord,p_iso_data(1:ibreak_p,3))

			call DBSINT(data_dim, p_iso_data(1:data_dim,1),  &
				 p_iso_data(1:data_dim,2), p_iso_ord,  &
				 p_iso_data(1:ibreak_p,3),  &
				 p_iso_cscoef(1:data_dim))

		endif

	endif

	continue

end subroutine read_numerical

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_numerical_bc
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use magnetic
	use exp_data
	use interpolating_functions
	use p_d_profile
	use pseudo_IMSL, only : DBSNAK, DBSINT, dbsval

	implicit none

	integer i

	if((bc_type==7).or.(bc_type==8).or.(bc_type==17)) then

		i = 0

		open(33,file='psioftheta.dat', status='old', action='read')

		do

			i = i+1
			read(33,*,end=299) psib_data(i,1), psib_data(i,2)

		enddo

	299	close(33)

		psib_dim = i-1
		ibreak_pb = psib_dim + psib_ord

		allocate(psib_cscoef(psib_dim))

		!this sets up the knots
		call DBSNAK(psib_dim, psib_data(1:psib_dim,1),  &
						psib_ord,psib_data(1:ibreak_pb,3))

		call DBSINT(psib_dim, psib_data(1:psib_dim,1),  &
			 psib_data(1:psib_dim,2), psib_ord,  &
			 psib_data(1:ibreak_pb,3),  &
			 psib_cscoef(1:psib_dim))

	endif

	! Set up TF part
	if(bc_type==7) then
	! big_Psi = psi on the boundary (no toroidal flow)

		big_Psib_dim =psib_dim
		ibreak_bPb = big_Psib_dim + big_Psib_ord

		allocate(big_Psib_cscoef(big_Psib_dim))

		big_Psib_data = psib_data

		!this sets up the knots
		call DBSNAK(big_Psib_dim, big_Psib_data(1:big_Psib_dim,1),  &
						big_Psib_ord,big_Psib_data(1:ibreak_bPb,3))

		call DBSINT(big_Psib_dim, big_Psib_data(1:big_Psib_dim,1),  &
			 big_Psib_data(1:big_Psib_dim,2), psib_ord,  &
			 big_Psib_data(1:ibreak_bPb,3),  &
			 big_Psib_cscoef(1:big_Psib_dim))

	! The (bc_type==17) case requires quantities that are not immediately available here.
	! For that reason, it is dealt with later.

	endif


	if(((bc_type==21).or.(bc_type==23).or.(bc_type==24)).and.(numerical_psi_diff)) then
	! need to read input vphi_edge and transform it into delta_psi

		i = 0

		open(33,file='vphi_edge.dat', status='old', action='read')
		! the file contains "theta", v_phi_edge and psi_diff;
		! v_phi is stored in column 4, since it's not used
		! NOTE: input should be checked...
		! (keep in mind that the conversion only depends on geomtry and physical constants)

		do

			i = i+1
			read(33,*,end=399) edge_psi_diff_data(i,1), edge_psi_diff_data(i,4), edge_psi_diff_data(i,2)

		enddo

	399	close(33)

		data_dim = i-1
		ibreak_psi_diff = data_dim + psi_diff_ord

		allocate(psi_diff_cscoef(data_dim))

		call interp_setup(data_dim,psi_diff_ord,edge_psi_diff_data(1:data_dim,1),  &
									edge_psi_diff_data(1:data_dim,2),edge_psi_diff_data(1:ibreak_psi_diff,3),  &
									psi_diff_cscoef(1:data_dim))

	endif

end subroutine read_numerical_bc


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_numerical_super
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this one needs to read B0, P, M_theta, M_phi
! data_dim is not assigned, but is the same for all data,
! so it's derivad from the data

	use constant
	use magnetic
	use exp_data
	use interpolating_functions
	use p_d_profile
!	use IMSL, only : dcsdec, dc2con
	use pseudo_IMSL, only : DBSNAK, DBSINT
	implicit none
	integer i

	integer ::  option = 1 !0=normal, 1=shape preserving

	real(kind=dkind), dimension(:,:), allocatable :: dummy
	integer, dimension(data_dim_lim) :: dumint
	integer :: itmax=2000

	if(input_EQDSK) return

! start from the pressure

	i = 0

	open(33,file='p_iso.dat', status='old', action='read')

	do 

		i = i+1
		read(33,*,end=99) p_iso_data(i,1), p_iso_data(i,2)

	enddo

99	close(33)

	data_dim = i-1

	ibreak_p = data_dim + p_iso_ord
	ibreak_B0 = data_dim + F_ord
	ibreak_fi = data_dim + mphi_ord
	ibreak_th = data_dim + mtheta_ord

	allocate(p_iso_cscoef(data_dim))
	allocate(F_cscoef(data_dim))
	allocate(mphi_cscoef(data_dim))
	allocate(mtheta_cscoef(data_dim))

	!this sets up the knots
	call DBSNAK(data_dim, p_iso_data(1:data_dim,1),  &
					p_iso_ord,p_iso_data(1:ibreak_p,3))

	call DBSINT(data_dim, p_iso_data(1:data_dim,1),  &
		 p_iso_data(1:data_dim,2), p_iso_ord,  &
		 p_iso_data(1:ibreak_p,3),  &
		 p_iso_cscoef(1:data_dim))

! next B0

	open(33,file='b0.dat', status='old', action='read')

	do i=1,data_dim
		read(33,*) F_data(i,1),F_data(i,2)
		F_data(i,2) = F_data(i,2)/rmajor
	enddo

	close(33)

	!this sets up the knots
	call DBSNAK(data_dim, F_data(1:data_dim,1),  &
					F_ord,F_data(1:ibreak_B0,3))

	call DBSINT(data_dim, F_data(1:data_dim,1),  &
		 F_data(1:data_dim,2), F_ord,  &
		 F_data(1:ibreak_B0,3),  &
		 F_cscoef(1:data_dim))

! next M_theta

	open(33,file='mtheta.dat', status='old', action='read')

	do i=1,data_dim
		read(33,*) mtheta_data(i,1),mtheta_data(i,2)
!		mtheta_data(i,2) = mtheta_data(i,2)*.0d0
	enddo

	close(33)

	!this sets up the knots
	call DBSNAK(data_dim, mtheta_data(1:data_dim,1),  &
					mtheta_ord,mtheta_data(1:ibreak_th,3))

	call DBSINT(data_dim, mtheta_data(1:data_dim,1),  &
		 mtheta_data(1:data_dim,2), mtheta_ord,  &
		 mtheta_data(1:ibreak_th,3),  &
		 mtheta_cscoef(1:data_dim))

! next M_phi

	open(33,file='mphi.dat', status='old', action='read')

	do i=1,data_dim
		read(33,*) mphi_data(i,1),mphi_data(i,2)
	enddo

	close(33)

	!this sets up the knots
	call DBSNAK(data_dim, mphi_data(1:data_dim,1),  &
					mphi_ord,mphi_data(1:ibreak_fi,3))

	call DBSINT(data_dim, mphi_data(1:data_dim,1),  &
		 mphi_data(1:data_dim,2), mphi_ord,  &
		 mphi_data(1:ibreak_fi,3),  &
		 mphi_cscoef(1:data_dim))

	if((numerical_omega).and.(omega_option==1) ) then

	 	numerical_mphi = .true.
		numerical_omega = .false.

	endif

!	set other flags

!	numerical_p_iso = .true.
!	numerical_F = .true.
!	numerical_mtheta = .true.

 	numerical_n 		=	.false. 			!density
!!$	numerical_omega 	=	.false. 			!toroidal rotation

	continue

	return

end subroutine read_numerical_super

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_grid
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant

	implicit none

	integer :: i, dsize
	real(kind=dkind), dimension(5000) :: xtemp

	i = 0

	open(33,file='input_xgrid.dat', status='old', action='read')

	do

		i = i+1
		read(33,*,end=111) xtemp(i)

	enddo

111	close(33)

	dsize = i-1

	allocate(dx_ext(1:dsize))

	do i=1,dsize

		dx_ext(i) = xtemp(i)

	enddo

	allocate(dz_ext(1:dsize))

	open(33,file='input_zgrid.dat', status='old', action='read')

	do i = 1,dsize

		read(33,*) dz_ext(i)

		if(i>2000) then
			continue
		endif

	enddo

	continue

end subroutine read_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine readeqdsk(filename)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use triangularity
	use magnetic
	use exp_data

	implicit none

	real(kind=dkind), dimension(:,:), allocatable :: psi_eqdsk
	real(kind=skind), dimension(:,:), allocatable :: out
	real(kind=dkind), dimension(:), allocatable :: dummy_arr
	real(kind=dkind), dimension(:), allocatable :: R_input, Z_input
	real(kind=dkind) :: dummy, psimax, rmin, rmax
	real(kind=dkind), dimension(1:20) :: dummy_list
	real(kind=dkind) :: zmin, zmax, rcenter_one,rcenter_two, z_center

	character  cdate*8, dummychar*48
	character(len=*) :: filename

	integer :: i, j
	integer :: iomap
	integer :: idummy, nxloc, nzloc
	integer :: grid_type_save	! this allows saving the incoming data

	grid_type_save = grid_type

	if(grid_type<0) grid_type = -(10 + grid_type)

	iomap = 22

	open(iomap,file=filename,form='formatted',action='read')

	! as far as FLOW is concerned, the first line is garbage

	read(iomap,9380)dummychar,cdate,idummy,nxloc,nzloc

	! the second line contains useful staff: x_size,z_size,rmajor, dummies

	read(iomap,9381) x_size,z_size,rmajor,dummy_list(1), dummy_list(2)
	! in FLOW the grid is in physical units
	rcenter = rmajor ! this would be the default in most codes

	! the third line is essentially garbage, but we can read psimax, just in case

     read(iomap,9381) dummy,dummy, psimax,dummy,b_phi_zero

	! as far as FLOW is concerned, the fourth and fifth lines are garbage

	read(iomap,9381) dummy_list(3),dummy_list(4),dummy_list(5),dummy_list(6),dummy_list(7)	!4th line
	read(iomap,9381) dummy_list(8),dummy_list(9),dummy_list(10),dummy_list(11),dummy_list(12)	!5th line

	! the sixth entry is F(psi), this one is needed by FLOW

	numerical_F = .true.

	ibreak_B0 = nxloc + F_ord

	if(allocated(F_cscoef)) deallocate(F_cscoef)
	! F_cscoef could have been created earlier in the read_numerical part:
	! if so, it needs to be replaced
	allocate(F_cscoef(nxloc))
	! F_data is already allocated, that will be changed
	F_data(:,:) = 0.d0

	! need to set up the "x" points for the interpolation

	do i=1,nxloc
		F_data(i,1) = 1.d0/(nxloc-1.d0)*(i-1.d0)
	enddo

	read(iomap,9381) (F_data(nxloc+1-i,2), i=1,nxloc)

	do i=1,nxloc
		F_data(i,2) = F_data(i,2)/rmajor
		! FLOW uses B0, not F
	enddo

	! interpolation setup has been moved to a separate routine
	call interp_setup(nxloc,F_ord, &
		F_data(1:nxloc,1),F_data(1:nxloc,2), &
		F_data(1:nxloc+F_ord,3),F_cscoef(1:nxloc))

	! the seventh entry is P(psi), this one is also needed by FLOW

	numerical_p_iso = .true.

	ibreak_p = nxloc + p_iso_ord

	if(allocated(p_iso_cscoef)) deallocate(p_iso_cscoef)
	! p_iso_cscoef could have been created earlier in the read_numerical part:
	! if so, it needs to be replaced
	allocate(p_iso_cscoef(nxloc))
	! p_iso_data is already allocated, that will be changed
	p_iso_data(:,:) = 0.d0

	! need to set up the "x" points for the interpolation

	do i=1,nxloc
		p_iso_data(i,1) = 1.d0/(nxloc-1.d0)*(i-1.d0)
	enddo

	read(iomap,9381) (p_iso_data(nxloc+1-i,2), i=1,nxloc)

	! interpolation setup has been moved to a separate routine
	call interp_setup(nxloc,p_iso_ord, &
		p_iso_data(1:nxloc,1),p_iso_data(1:nxloc,2), &
		p_iso_data(1:nxloc+p_iso_ord,3),p_iso_cscoef(1:nxloc))

	! entries 8 to 11 are not needed

	allocate(dummy_arr(1:nxloc))
	allocate(psi_eqdsk(1:nxloc,1:nzloc))

	read(iomap,9381) (dummy_arr(i), i=1,nxloc)	! 8th
	read(iomap,9381) (dummy_arr(i), i=1,nxloc)	! 9th
	read(iomap,9381) ((psi_eqdsk(i,j), i=1,nxloc), j=1,nzloc)	! 10th
	read(iomap,9381) (dummy_arr(i), i=1,nxloc)	! 11th

	deallocate(dummy_arr)

	! the twelwth entry is the shape of the boundary, this one is needed

	read(iomap,1204) theta_points, idummy

	allocate(R_input(theta_points))
	allocate(Z_input(theta_points))

	read(iomap,9381) (R_input(i),Z_input(i), i=1,theta_points)

	! the next entry is the limiter, not needed by FLOW,
	! and the last lines contain some equilibrium values, also not needed by FLOW

	! read the limiter nevertheless, to check the grid

	allocate(dummy_arr(2*idummy))

	read(iomap,9381) (dummy_arr(i),dummy_arr(i+idummy),i=1,idummy)

	continue

	close(unit=iomap,status='keep')

	rmin = rmajor
	rmax = 0.d0

	do i = 1,theta_points

		if(R_input(i)>rmax) rmax = R_input(i)
		if(R_input(i)<rmin) rmin = R_input(i)

	enddo

	rcenter_one = (rmax+rmin)/2.d0

	rmin = rmajor
	rmax = 0.d0
	zmin = rmajor
	zmax = 0.d0

	do i = 1,idummy

		if(dummy_arr(i)>rmax) rmax = dummy_arr(i)
		if(dummy_arr(i)<rmin) rmin = dummy_arr(i)

		if(dummy_arr(i+idummy)>zmax) zmax = dummy_arr(i+idummy)
		if(dummy_arr(i+idummy)<zmin) zmin = dummy_arr(i+idummy)

	enddo

	rcenter_two = (rmax+rmin)/2.d0
	z_center = (zmax+zmin)/2.d0	! just out of curiosity

	deallocate(dummy_arr)

	! the input eqdsk psi is saved three times, since the definition of the grid is ambiguous

	allocate(out(nxloc,nzloc))

	out(1:nxloc,1:nzloc) = psi_eqdsk(1:nxloc,1:nzloc)

	call set_grid(nxloc,nzloc)

	call radial_plot(out,psi_eqdsk,nxloc,nzloc,"psi_eqdsk1",nxloc/2)
	call psi_boundary_plot(nxloc,nzloc,psi_eqdsk)

	rcenter = rcenter_two

	call set_grid(nxloc,nzloc)

	out(1:nxloc,1:nzloc) = psi_eqdsk(1:nxloc,1:nzloc)
	call radial_plot(out,psi_eqdsk,nxloc,nzloc,"psi_eqdsk2",nxloc/2)

	rcenter = rcenter_one	! this is the final rcenter
	call set_grid(nxloc,nzloc)
	call radial_plot(out,psi_eqdsk,nxloc,nzloc,"psi_eqdsk3",nxloc/2)

	deallocate(out)
	deallocate(psi_eqdsk)

	rcenter = rcenter_one	! this is the final rcenter, repeated here for clarity

	call plasma_shape_conversion(R_input,Z_input)

	deallocate(R_input)
	deallocate(Z_input)

	grid_type = grid_type_save

	continue

	return

1204    format(2I5)	!format(5I5)
9380    format(A40,A8,3I4)
9381    format(1P,5E16.9)

end subroutine readeqdsk

