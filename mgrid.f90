module mgrid_mod

! Feb 20 2013
! mgrid has been turned into a module, because an explicit interface is needed
! for optional arguments in mgrid, and I don't know how to create an interface.
! some routines were movedto inter_grid (the ones that are called from within
! trans_solve, which is itself referenced from here...)

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine mgrid_wrapper(psi,rho,residual,b_phi, n_den, psi_diff, big_Psi)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this is used to add a single-fluid calculation at the beginning of the two-fluid one, to get a good initial guess

	use constant
	use solver

	implicit none

	real (kind=dkind), dimension(1:n,1:n), intent(inout) :: psi,rho,residual,b_phi, n_den, psi_diff, big_Psi
	real (kind=dkind), dimension(:,:), allocatable :: psi_temp,rho_temp,residual_temp,b_phi_temp,  &
																					n_den_temp, psi_diff_temp, big_Psi_temp
	integer :: eq_type_temp, n_res, bc_type_temp
	logical :: accelerate_temp !!$, Broot_zero
	real (kind=dkind) :: fix_orp_temp

!!$	Broot_zero = .false.

	if(eq_type>=10) then
	! start with one-fluid for guess
	! NOTE: WE MAY WANT TO RUN MORE THAN ONE GRID IN ONE-FLUID: TO BE VERIFIED!

		eq_type_temp = eq_type
		n_res = n
		bc_type_temp = bc_type

		if(Broot==0) then
			n = inter_switch
			! June 12 2019
			! Don't remember why this is here. Try to stick with n_min.
			n = n_min
		else
			n = n_min
		endif

		eq_type = 1

		accelerate_temp = accelerate
		accelerate = accelerate_0
		fix_orp_temp = fix_orp
		fix_orp = fix_orp_0

		! allocation section for temporaty files

		allocate(psi_temp(1:n,1:n))
		psi_temp(1:n,1:n) = 0.d0

		allocate(rho_temp(1:n,1:n))
		rho_temp(1:n,1:n) = 0.d0

		allocate(residual_temp(1:n,1:n))
		residual_temp(1:n,1:n) = 0.0

		allocate(b_phi_temp(1:n,1:n))
		b_phi_temp(1:n,1:n) = 0.0

		allocate(n_den_temp(1:n,1:n))
		n_den_temp(1:n,1:n) = 0.0

		allocate(psi_diff_temp(1:n,1:n))
		psi_diff_temp(1:n,1:n) = 0.0

		allocate(big_Psi_temp(1:n,1:n))
		big_Psi_temp(1:n,1:n) = 0.0

		! end of allocation section

		print*, '   '
		print*, '----------------------------------------'
		print*, 'running one-fluid equilibrium'
		print*, 'as initialization of two-fluid one'
		print*, '----------------------------------------'
!!$		if(Broot==0) then
!!$			Broot = 1
!!$			print*, 'changing Broot from 0 to 1 for initial guess'
!!$			print*, '----------------------------------------'
!!$			Broot_zero = .true.
!!$		endif
		print*, '   '
		if(bc_type>20) then
			bc_type=1
		elseif(bc_type==17) then
			bc_type = 7
		endif

!		call mgrid(psi_temp,rho_temp,residual_temp,b_phi_temp, n_den_temp, psi_diff_temp, big_Psi_temp)
		call mgrid(psi=psi_temp,rho=rho_temp,residual=residual_temp,b_phi=b_phi_temp, n_den=n_den_temp,  &
						psi_diff=psi_diff_temp, big_Psi=big_Psi_temp, nguess=1)

!!$		if(Broot_zero) Broot = 0

		call two_fluid_initial(psi_temp,rho_temp, n_den_temp, psi_diff_temp, big_Psi_temp, n, n)

		big_Psic = maxval(big_Psi_temp)

		eq_type = eq_type_temp
		n = n_res
		if(Broot==0) then
!			n = inter_switch
			! June 12 2019
			! Don't remember why this is here. Try to stick with n.
!			n = n_min
		endif
		bc_type = bc_type_temp

		accelerate = accelerate_temp
		fix_orp = fix_orp_temp

		print*, '   '
		print*, '----------------------------------------'
		print*, 'entering two-fluid equilibrium'
		print*, '----------------------------------------'
		print*, '   '

		call mgrid(psi=psi,rho=rho,residual=residual,b_phi=b_phi, n_den=n_den, psi_diff=psi_diff, big_Psi=big_Psi,  &
						psi_guess=psi_temp, n_den_guess=n_den_temp, psi_diff_guess=psi_diff_temp, big_Psi_guess=big_Psi_temp, nguess=n_min)

	else

	!	call mgrid(psi,rho,residual,b_phi, n_den, psi_diff, big_Psi)
		call mgrid(psi=psi,rho=rho,residual=residual,b_phi=b_phi, n_den=n_den, psi_diff=psi_diff, big_Psi=big_Psi, nguess=1)

	endif

end subroutine mgrid_wrapper

! This is a VERY simple multi-grid type of calculation, wherein the
! solution is found first on a coarse grid, interpolated to a finer
! grid and refined.  This is repeated until the desired resolution is
! reached. On input psi[1..n][1..n] must have a dimension n of the
! form 2^j + 1 for some integer j.  (j is actually the number of grid
! levels used in the solution, called ng below.)  */

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine mgrid(psi,rho,residual,b_phi, n_den, psi_diff, big_Psi,  &
							psi_guess, n_den_guess, psi_diff_guess, big_Psi_guess, nguess)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	!CORREGGERE CONDIZIONI AL CONTORNO!!!!!!

  use constant
  use solver
 ! use triangularity, only : ind_bound, coord_bound

  implicit none

  real (kind=dkind) :: eps = 1.0d-7
  ! small_eps is for the first solution and should be SMALL
  ! for 17 x 17 small_eps = 1.0d-9 works pretty well
  real (kind=dkind), parameter :: small_eps = 1.0d-9
  integer, parameter :: ngmax = 15 ! Maximum number of allowable grids
  ! I highly recommend puting n_min > 9 since the convergence properties
  ! for small grids requires special treatment for the over-relaxation 
  ! parameter and the Chebyshev acceleration.  Moreover, the value of 
  ! Psi_center for small grids differs by as much as 50% from the 
  ! value for grids with n ~ 17.
  real (kind=dkind), dimension(1:n,1:n), intent(inout) :: psi,rho,residual,b_phi, n_den, psi_diff, big_Psi
  integer :: nguess
  real (kind=dkind), dimension(1:nguess,1:nguess), intent(in), optional :: psi_guess, n_den_guess, psi_diff_guess, big_Psi_guess
  logical :: all_present
  real (kind=dkind) :: orp=-1.0d0
  integer :: err
  integer :: ng ! Total number of grid levels in the calculation
  integer :: i,j,nn,nc
  integer :: ii, jj
  ! epsi and opsi are even and odd grids
  real (kind=dkind), dimension(:,:), allocatable :: epsi,opsi
  real (kind=dkind), dimension(:,:), allocatable :: erho,orho
  real (kind=dkind), dimension(:,:), allocatable :: eb_phi,ob_phi
  real (kind=dkind), dimension(:,:), allocatable :: en_den,on_den
  real (kind=dkind), dimension(:,:), allocatable :: epsi_diff,opsi_diff
  real (kind=dkind), dimension(:,:), allocatable :: ebig_Psi,obig_Psi
  integer :: alloc_stat
  integer :: acc_switch = 65
  real(kind=dkind) :: inorm ! Unutilized here, but important for convergence in ngs_solve

	if(grid_type<=-10) call pack_grid

  if(restart) then
	open (17, file='FLOW_n.dat', status='old', action='read')
	read(17,*) n_min !this still assumes same number of points in each direction
	close(17)
  endif

  nn = n
  ng = 0
  do
     nc = nn/2 + 1
     ng = ng + 1
     if(nc < n_min) exit
     nn = nc
  end do

	if((present(psi_guess)).and.(present(n_den_guess)).and.(present(psi_diff_guess)).and.(present(big_Psi_guess))) then
		all_present = .true.
	endif

  if(n == nn) then
	 if(restart) then
		 pause 'useless restart: n_min = n_max! The program will abort'
		 stop
	 endif
     ! This is a silly special case
	 call set_grid(n,n)
	 if(eq_type<10) then
     ! print *, "guess the solution"
		 call guess_soln(psi,n,n)
		 ! print *, "Initialize the density"
		 call initialize_density(psi,rho,n,n)
		 call initialize_b(b_phi,psi,rho,n,n)
	 elseif(eq_type>=10) then
		if(all_present) then
			do jj = 1,n
			do ii = 1,n
				psi(ii,jj) = psi_guess(ii,jj)
				big_Psi(ii,jj) = big_Psi_guess(ii,jj)
				psi_diff(ii,jj) = psi_diff_guess(ii,jj)
				n_den(ii,jj) = n_den_guess(ii,jj)
			enddo
			enddo
		else
			call guess_soln_TF(psi,big_Psi,psi_diff,n_den,n,n)
			do jj = 1,n
			do ii = 1,n
				if(present(psi_guess)) psi(ii,jj) = psi_guess(ii,jj)
				if(present(big_Psi_guess)) then
					big_Psi(ii,jj) = big_Psi_guess(ii,jj)
					psi_diff(ii,jj) = big_Psi_guess(ii,jj) - psi(ii,jj)
				endif
				if(present(psi_diff_guess)) psi_diff(ii,jj) = psi_diff_guess(ii,jj)
				if(present(n_den_guess)) n_den(ii,jj) = n_den_guess(ii,jj)
			enddo
			enddo
		endif
	 endif
	 ! print *, "solve the problem"
     ! (Newton) Gauss-Seidel solution (Correcting PsiC)
    call ngs_solve_wrapper(psi,rho,residual,b_phi,big_Psi,psi_diff,n_den,n,n,min_it,max_it,small_eps,0,orp)
     return
  end if

  if(n .ne. 1 + (nn-1)*2**(ng-1)) then
     print *, "[mgrid]: n-1 must be a power of 2"
     print *, "[mgrid}: n =",n,"and nn =",nn
     print *, nn*2**(ng-1) - (1 + 2**(ng-1))
     return
  end if

  if(ng > ngmax) then
     print *, "[mgrid]: increase ngmax to at least ",ng
     return
  end if

  ! Find the solution on the coarsest grid: nn x nn
  ! nn = 3
  allocate(opsi(1:nn,1:nn),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
     return
  end if

  allocate(orho(1:nn,1:nn),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
     return
  end if

  allocate(ob_phi(1:nn,1:nn),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
     return
  end if

  allocate(on_den(1:nn,1:nn),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
     return
  end if

  allocate(opsi_diff(1:nn,1:nn),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
     return
  end if

  allocate(obig_Psi(1:nn,1:nn),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
     return
  end if

  print *, "Grid  1  - Size ",nn,"x",nn

  call set_grid(nn,nn)
  call initialize_bc(nn)

 if(grid_type>=1) then

	if(nn<=acc_switch) then
		eps = eps/10.d0
	else
		eps = eps/5.d0
	endif

 endif

  if(restart) then

	call read_restart_data(opsi,nn,nn,"psi")
	call read_restart_data(orho,nn,nn,"rho")
	call read_restart_data(ob_phi,nn,nn,"b_phi")

	psic = 0.d0

	do j=1,nn
	do i=1,nn

		if(opsi(i,j)>psic) psic=opsi(i,j)

	enddo
	enddo

  else

		if(eq_type<10) then
			call guess_soln(opsi,nn,nn)
			call initialize_density(opsi,orho,nn,nn)
			call initialize_b(ob_phi,opsi,orho,nn,nn)
		elseif(eq_type>=10) then
			if(all_present) then
				do jj = 1,nn
				do ii = 1,nn
					opsi(ii,jj) = psi_guess(ii,jj)
					obig_Psi(ii,jj) = big_Psi_guess(ii,jj)
					opsi_diff(ii,jj) = psi_diff_guess(ii,jj)
					on_den(ii,jj) = n_den_guess(ii,jj)
				enddo
				enddo
			else
				call guess_soln_TF(opsi,obig_Psi,opsi_diff,n_den,nn,nn)
				do jj = 1,nn
				do ii = 1,nn
					if(present(psi_guess)) opsi(ii,jj) = psi_guess(ii,jj)
					if(present(big_Psi_guess)) then
						obig_Psi(ii,jj) = big_Psi_guess(ii,jj)
						opsi_diff(ii,jj) = big_Psi_guess(ii,jj) - psi(ii,jj)
					endif
					if(present(psi_diff_guess)) opsi_diff(ii,jj) = psi_diff_guess(ii,jj)
					if(present(n_den_guess)) on_den(ii,jj) = n_den_guess(ii,jj)
				enddo
				enddo
			endif
		endif

 	  call ngs_solve_wrapper(opsi,orho,residual,ob_phi,obig_Psi,opsi_diff,on_den,nn,nn,min_it,max_it,small_eps,0,orp)

  endif

  ! The value of 1.0d-10 is chosen empirically.
  ! Smaller values do not appear to be accessible, i.e.
  ! this is the noise level for the rest of the routines...
  ! given the initialization proceedure.

  ! Nested iteration loop
  do j=2, ng-1
     ! Allocate storage for the next level finer grids
     nc = nn
     nn = 2*nn - 1

     print *, "Grid ",j," - Size ",nn,"x",nn

     if(nn <= 33) then 
        eps = 1.0d-7
     else if(nn <= 65) then 
        eps = 1.0d-5
     else if(nn <= 129) then 
        eps = 1.0d-2
        orp = 0.8d0
	 else if(nn <= 257) then 
        eps = 5.0d-2
        orp = 0.8d0
     else 
        eps = 1.0d-1
        orp = 0.8d0
     endif

	 if(grid_type>=1) then

		if(nn<=acc_switch) then
			eps = eps/10.d0
		else
			eps = eps/5.d0
		endif

	 endif

     if(modulo(j,2) == 0) then
        ! For j even we solve on the grid epsi 
        allocate(epsi(1:nn,1:nn))
        if(alloc_stat > 0) then
           print *, "Allocation Error"
           return
        end if

        allocate(erho(1:nn,1:nn),stat = alloc_stat)
        if(alloc_stat > 0) then
           print *, "Allocation Error"
           return
        end if

		allocate(eb_phi(1:nn,1:nn),stat = alloc_stat)
        if(alloc_stat > 0) then
           print *, "Allocation Error"
           return
        end if

		allocate(en_den(1:nn,1:nn),stat = alloc_stat)
		if(alloc_stat > 0) then
		 print *, "Allocation Error"
		 return
		end if

		allocate(epsi_diff(1:nn,1:nn),stat = alloc_stat)
		if(alloc_stat > 0) then
		 print *, "Allocation Error"
		 return
		end if

		allocate(ebig_Psi(1:nn,1:nn),stat = alloc_stat)
		if(alloc_stat > 0) then
		 print *, "Allocation Error"
		 return
		end if
		! set grid
		call set_grid(nn,nn)

        ! Interpolate the solution to the next level finer grid
        call interp_nonuni(opsi,nc,epsi,nn,err)
        if(err > 0) return
        call interp_nonuni(orho,nc,erho,nn,err)
        if(err > 0) return
		call interp_nonuni(ob_phi,nc,eb_phi,nn,err)
        if(err > 0) return
		call interp_nonuni(on_den,nc,en_den,nn,err)
        if(err > 0) return
		call interp_nonuni(opsi_diff,nc,epsi_diff,nn,err)
        if(err > 0) return
		call interp_nonuni(obig_Psi,nc,ebig_Psi,nn,err)
        if(err > 0) return

		! For free-boundary calculations, fix sort_grid as it was overwritten by set_grid
		if(bc_type==17) then
			if(LCFS==-1) then
				call update_sort_grid(epsi,nn,nn,inorm)
			elseif(LCFS==1) then
				call update_sort_grid(ebig_Psi,nn,nn,inorm)
			endif
		endif

        ! Free up the previous grid
        deallocate(opsi)
        deallocate(orho)
		deallocate(ob_phi)
        deallocate(obig_psi)
        deallocate(on_den)
        deallocate(opsi_diff)

		! initialize boundary conditions
		call initialize_bc(nn)


        ! Apply boundary conditions to psi and rho
		if(eq_type<10) then
	        call bc_psi_rho0(epsi,erho,nn,nn)
		elseif(eq_type>=10) then
			call bc_TF_0(epsi,ebig_Psi,epsi_diff,en_den,nn,nn,0)
		endif

        ! (Newton) Gauss-Seidel solution (Correcting PsiC)

 	  call ngs_solve_wrapper(epsi,erho,residual,eb_phi,ebig_Psi,epsi_diff,en_den,nn,nn,min_it,max_it,eps,0,orp)
     else
        ! Repeat switching epsi and opsi 
        allocate(opsi(1:nn,1:nn),stat = alloc_stat)
        if(alloc_stat > 0) then
           print *, "Allocation Error"
           return
        end if

        allocate(orho(1:nn,1:nn),stat = alloc_stat)
        if(alloc_stat > 0) then
           print *, "Allocation Error"
           return
        end if

        allocate(ob_phi(1:nn,1:nn),stat = alloc_stat)
        if(alloc_stat > 0) then
           print *, "Allocation Error"
           return
        end if

		allocate(on_den(1:nn,1:nn),stat = alloc_stat)
		if(alloc_stat > 0) then
		 print *, "Allocation Error"
		 return
		end if

		allocate(opsi_diff(1:nn,1:nn),stat = alloc_stat)
		if(alloc_stat > 0) then
		 print *, "Allocation Error"
		 return
		end if

		allocate(obig_Psi(1:nn,1:nn),stat = alloc_stat)
		if(alloc_stat > 0) then
		 print *, "Allocation Error"
		 return
		end if

		! set grid
		call set_grid(nn,nn)

        ! Interpolate the solution to the next level finer grid
        call interp_nonuni(epsi,nc,opsi,nn,err)
        if(err > 0) return
        call interp_nonuni(erho,nc,orho,nn,err)
        if(err > 0) return
		call interp_nonuni(eb_phi,nc,ob_phi,nn,err)
        if(err > 0) return
		call interp_nonuni(en_den,nc,on_den,nn,err)
        if(err > 0) return
		call interp_nonuni(epsi_diff,nc,opsi_diff,nn,err)
        if(err > 0) return
		call interp_nonuni(ebig_Psi,nc,obig_Psi,nn,err)
        if(err > 0) return

		! For free-boundary calculations, fix sort_grid as it was overwritten by set_grid
		if(bc_type==17) then
			if(LCFS==-1) then
				call update_sort_grid(opsi,nn,nn,inorm)
			elseif(LCFS==1) then
				call update_sort_grid(obig_Psi,nn,nn,inorm)
			endif
		endif

        ! Free up the previous grid
        deallocate(epsi)
        deallocate(erho)
		deallocate(eb_phi)
        deallocate(ebig_psi)
        deallocate(en_den)
        deallocate(epsi_diff)

		! initialize boundary conditions
		call initialize_bc(nn)

        ! Apply boundary conditions to psi and rho
		if(eq_type<10) then
	        call bc_psi_rho0(opsi,orho,nn,nn)
		elseif(eq_type>=10) then
			call bc_TF_0(opsi,obig_Psi,opsi_diff,on_den,nn,nn,0)
		endif

        ! (Newton) Gauss-Seidel solution (Correcting PsiC)
	    call ngs_solve_wrapper(opsi,orho,residual,ob_phi,obig_Psi,opsi_diff,on_den,nn,nn,min_it,max_it,eps,0,orp)
     end if
  end do

 print *, "Grid ",ng," - Size ",n,"x",n
! if(n > 100) eps = dmax1(eps,1.0d-2)

	 if(n <= 33) then 
		eps = 1.0d-7
	 else if(n <= 65) then 
		eps = 1.0d-5
	 else if(n <= 129) then 
		eps = 1.0d-2
		orp = 0.8d0
	 else if(n <= 257) then 
		eps = 5.0d-2
		orp = 0.8d0
	 else 
		eps = 1.0d-1
		orp = 0.8d0
	 endif

	 if(grid_type>=1) then

		if(n<=acc_switch) then
			eps = eps/10.d0
		else
			eps = eps/5.d0
		endif

	 endif

 ! set grid
 call set_grid(n,n)

 ! last step, interpolate onto the final grid
 if(modulo(ng,2) == 0) then
    ! Interpolate the solution to the next level finer grid
    call interp_nonuni(opsi,nn,psi,n,err)
    if(err > 0) return
    call interp_nonuni(orho,nn,rho,n,err)
    if(err > 0) return
	call interp_nonuni(ob_phi,nn,b_phi,n,err)
    if(err > 0) return
	call interp_nonuni(on_den,nn,n_den,n,err)
    if(err > 0) return
	call interp_nonuni(opsi_diff,nn,psi_diff,n,err)
    if(err > 0) return
	call interp_nonuni(obig_Psi,nn,big_Psi,n,err)
    if(err > 0) return

    ! Free up the previous grid
    deallocate(opsi)
    deallocate(orho)
	deallocate(ob_phi)
    deallocate(obig_psi)
    deallocate(on_den)
    deallocate(opsi_diff)

 else
    ! Interpolate the solution to the next level finer grid
    call interp_nonuni(epsi,nn,psi,n,err)
    if(err > 0) return
    call interp_nonuni(erho,nn,rho,n,err)
    if(err > 0) return
	call interp_nonuni(eb_phi,nn,b_phi,n,err)
    if(err > 0) return
	call interp_nonuni(en_den,nn,n_den,n,err)
    if(err > 0) return
	call interp_nonuni(epsi_diff,nn,psi_diff,n,err)
    if(err > 0) return
	call interp_nonuni(ebig_Psi,nn,big_Psi,n,err)
    if(err > 0) return

    ! Free up the previous grid
    deallocate(epsi)
    deallocate(erho)
	deallocate(eb_phi)
    deallocate(ebig_psi)
    deallocate(en_den)
    deallocate(epsi_diff)

 end if

 ! For free-boundary calculations, fix sort_grid as it was overwritten by set_grid
 if(bc_type==17) then
	if(LCFS==-1) then
		call update_sort_grid(psi,n,n,inorm)
	elseif(LCFS==1) then
		call update_sort_grid(big_Psi,n,n,inorm)
	endif
 endif

 ! initialize boundary conditions
 call initialize_bc(n)

    ! Apply boundary conditions to psi and rho
	if(eq_type<10) then
	    call bc_psi_rho0(psi,rho,nn,nn)
	elseif(eq_type>=10) then
		call bc_TF_0(psi,big_Psi,psi_diff,n_den,n,n,0)
	endif

 ! Initialize the residual matrix to 0
 residual(1:n,1:n) = 0.0d0

 ! (Newton) Gauss-Seidel solution (Correcting PsiC)
 call ngs_solve_wrapper(psi,rho,residual,b_phi,big_Psi,psi_diff,n_den,n,n,min_it,max_it,eps,1,orp)

	! just in case

	if(allocated(opsi)) deallocate(opsi)
    if(allocated(orho)) deallocate(orho)
	if(allocated(ob_phi)) deallocate(ob_phi)
	if(allocated(on_den)) deallocate(on_den)
	if(allocated(opsi_diff)) deallocate(opsi_diff)
	if(allocated(obig_Psi)) deallocate(obig_Psi)

    if(allocated(epsi)) deallocate(epsi)
    if(allocated(erho)) deallocate(erho)
	if(allocated(eb_phi)) deallocate(eb_phi)
	if(allocated(en_den)) deallocate(en_den)
	if(allocated(epsi_diff)) deallocate(epsi_diff)
	if(allocated(ebig_Psi)) deallocate(ebig_Psi)

end subroutine mgrid


end module mgrid_mod