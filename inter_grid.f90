! Course-to-fine prolongation by bilinear interpolation.  The coarse
! grid has nc grid points and the fine grid has nf = 2*nc-1 grid
! points.  The coarse-grid solution is input as uc and the 
! fine-grid solution is returned in uf.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine interp(uc,nc,uf,nf,err)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  implicit none

  integer, intent(in) :: nc,nf
  integer, intent(inout) :: err
  real (kind=dkind), dimension(1:nc,1:nc), intent(in) :: uc ! Coarse Grid
  real (kind=dkind), dimension(1:nf,1:nf), intent(out) :: uf ! Fine Grid
  integer ic,iif,jc,jf

  if ( nf .ne. 2*nc - 1 ) then
     print *, "[interp]: Array Dimensions Incorrect"
     err = 1 ! Error
     return
  end if

  ! Do elements that are copies
  jf = 1
  do jc = 1, nc
     iif = 1
     do ic = 1, nc
        uf(iif,jf) = uc(ic,jc)

        iif = iif + 2
     end do
     jf = jf + 2
  end do

  ! Interpolate horizontally for iif odd and jf even
  do jf = 1, nf, 2
     do iif = 2, nf-1, 2
        uf(iif,jf) = 0.5*(uf(iif-1,jf) + uf(iif+1,jf))
     end do
  end do

  ! Interpolate vertically for all iif and jf odd
  do jf = 2, nf-1, 2
     do iif = 1, nf
        uf(iif,jf) = 0.5*(uf(iif,jf-1) + uf(iif,jf+1))
     end do
  end do

  err = 0 ! Success
end subroutine interp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine interp_nonuni(uc,nc,uf,nf,err)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  implicit none

  integer, intent(in) :: nc,nf
  integer, intent(inout) :: err
  real (kind=dkind), dimension(1:nc,1:nc), intent(in) :: uc ! Coarse Grid
  real (kind=dkind), dimension(1:nf,1:nf), intent(out) :: uf ! Fine Grid
  integer ic,iif,jc,jf

  if ( nf .ne. 2*nc - 1 ) then
     print *, "[interp]: Array Dimensions Incorrect"
     err = 1 ! Error
     return
  end if

  ! Do elements that are copies
  jf = 1
  do jc = 1, nc
     iif = 1
     do ic = 1, nc
        uf(iif,jf) = uc(ic,jc)

        iif = iif + 2
     end do
     jf = jf + 2
  end do

  ! Interpolate horizontally for iif odd and jf even
  do jf = 1, nf, 2
     do iif = 2, nf-1, 2
!        uf(iif,jf) = 0.5*(uf(iif-1,jf) + uf(iif+1,jf))
		uf(iif,jf) = (dx_a(iif)*uf(iif-1,jf) + dx_a(iif-1)*uf(iif+1,jf))  &
						/(dx_a(iif-1) + dx_a(iif+1))
     end do
  end do

  ! Interpolate vertically for all iif and jf odd
  do jf = 2, nf-1, 2
     do iif = 1, nf
!        uf(iif,jf) = 0.5*(uf(iif,jf-1) + uf(iif,jf+1))
		uf(iif,jf) = (dz_a(jf)*uf(iif,jf-1) + dz_a(jf-1)*uf(iif,jf+1))  &
						/(dz_a(jf-1) + dz_a(jf+1))
     end do
  end do

  err = 0 ! Success
end subroutine interp_nonuni

!-------------------------------------------------------------------------------------------------------------------------------------
! routines formerly in mgrid have also been moved here (they really weren't related to the multigrid...)
!-------------------------------------------------------------------------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine set_grid(nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : rcenter, zcenter, x_coord, z_coord, x_size, z_size,  &
						 grid_type, dx, dz, dx_a, dz_a, rmajor, dkind,  &
						 mx0, px0, mz0, n_min, sort_grid, dx_ext, dz_ext

	use solver, only : g_ratio, n_temp, grid_delta, grid_delta_2,  &
						radius, radius_1_3, rtsafe

	use array_dimensions

	implicit none

	integer :: nx,nz
	integer :: i, j, m, p, q

	real(kind=dkind) :: x_L, x_R, alpha, beta, delta, dx0, dz0
	real(kind=dkind) :: L1, L2, L3
	real(kind=dkind) :: xg(2)
	real(kind=dkind) :: ex,ez
	real(kind=dkind) :: temp, dummy, rminor

	if(allocated(x_coord)) deallocate(x_coord)
	if(allocated(z_coord)) deallocate(z_coord)

	if(allocated(dx_a)) deallocate(dx_a)
	if(allocated(dz_a)) deallocate(dz_a)

	allocate(x_coord(0:nx+1))
	allocate(z_coord(0:nz+1))

	allocate(dx_a(0:nx+1))
	allocate(dz_a(0:nz+1))

	dx = x_size/(nx - 1.0d0)
	dz = z_size/(nz - 1.0d0)

	x_L = rcenter - 0.5d0*x_size
	x_R = rcenter + 0.5d0*x_size

	if(grid_type==0) then

		dx_a = dx
		dz_a = dz

	elseif(grid_type==1) then

		temp = 0.25d0*(rmajor+x_size)
		g_ratio = x_size/(temp-x_L)
		n_temp = nx

		xg(1) = 1.d0+1.d-8
		xg(2) = 2.d0

		alpha = rtsafe(grid_delta,xg(1),xg(2),1.d-16,100000)

		dx0 = (temp-x_L)*(1.d0-alpha)/(1.d0-alpha**((nx-1)/2))

		dx_a(1) = dx0
		dx_a(0) = dx_a(1)

		do i=2,nx+1

			dx_a(i) = dx_a(i-1)*alpha

		enddo

		dz_a = dz

		continue

	elseif(grid_type==2) then

		g_ratio = 0.4d0
		n_temp = (nx-1)/2

		xg(1) = 1.d0+1.d-8
		xg(2) = 2.d0

		alpha = rtsafe(grid_delta,xg(1),xg(2),1.d-16,100000)

		dx0 = (1-g_ratio)/(n_temp*alpha**(n_temp-1.d0))*x_size
		dx_a(1) = dx0
		dx_a(0) = dx_a(1)

		do i = 2,n_temp

			dx_a(i) = dx_a(i-1)*alpha

		enddo

		do i = n_temp+1,nx

			dx_a(i) = dx_a(i-1)

		enddo

		dz_a = dz

		continue

	elseif(grid_type==3) then

		j = 0

		do

			j = j+1
			if(dx_ext(j)==nx*1.d0) exit

		enddo

		do i = 1,nx-1

			dx_a(i) = dx_ext(j+i)*x_size

		enddo

		dx_a(0) = dx_a(1);
		dx_a(nx) = dx_a(nx-1);
		dx_a(nx+1) = dx_a(nx);


		do i = 1,nz-1

			dz_a(i) = dz_ext(j+i)*z_size

		enddo

		dz_a(0) = dz_a(1);
		dz_a(nz) = dz_a(nz-1);
		dz_a(nz+1) = dz_a(nz);

		continue

	else

		print*, 'error in grid_type: grid_type=', grid_type
		pause
		stop

	endif

	x_coord(1) = x_L
	z_coord(1) = zcenter - 0.5d0*z_size 

	x_coord(0) = x_coord(1) - dx_a(1)
	z_coord(0) = z_coord(1) - dz_a(1)

	do i=2,nx+1

		x_coord(i) = x_coord(i-1) + dx_a(i-1)

	enddo

	do i=2,nz+1

		z_coord(i) = z_coord(i-1) + dz_a(i-1)

	enddo

	!check grid

	if ( (dabs(x_coord(nx)-x_R)/x_R>1.d-12).or.  &
		(dabs(z_coord(nz)-z_coord(1)-z_size)/z_size>1.d-12) ) then

		print*, 'problem in set_grid:'
		print*, 'x_R=', x_R
		print*, 'x_coord(nx)=', x_coord(nx)
		print*, 'z_size=', z_size
		print*, 'z_coord(nz)=', z_coord(nz)
		pause
		stop

	endif

	! assign grid point location
	! 1 = internal, -1 = external, 0 = on the edge (or very very close)

	if(allocated(sort_grid)) deallocate(sort_grid)

	allocate(sort_grid(0:nx+1,0:nz+1))

	call set_sort_grid(nx,nz)
	! this has been moved to a separate routine
	! to simplify things in the separatrix case

	continue

!	open(33,file='grid.plt')
!
!	write(33,*)'TITLE="grid"'
!	write(33,*)'Variables =" R [m] ","z [m]", "boh"'
!	write(33,*)'ZONE I=',nx,',J=',nz,',F=Point'
!
!	do j = 1,nz
!	do i = 1,nx
!
!		write(33,*) x_coord(i), z_coord(j),  sort_grid(i,j,0)
!
!	enddo
!	enddo
!
!	close(33)

	return

end subroutine set_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine set_sort_grid(nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, sort_grid, dx_a, dz_a, x_coord, z_coord
	use triangularity, only : tri_type
	use solver, only : radius, radius_1_3
	implicit none

	real(kind=dkind), parameter :: small = 1.d-2

	integer :: nx,nz
	integer :: i,j,p
	real(kind=dkind) :: ex, ez, rminor, dx, dz, dummy

	sort_grid = -1

	do j=1,nz
	do i=1,nx

		call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)

		if((ex*ex + ez*ez) > rminor**2) then

			sort_grid(i,j) = -1

		elseif( (ex*ex + ez*ez) >= (rminor - small*sqrt(dx_a(i)**2+dz_a(j)**2) )**2 ) then

			sort_grid(i,j) = 0

		else

			sort_grid(i,j) = 1

		end if

	enddo
	enddo

	if(tri_type==13) then

		do j=1,nz
		do i=1,nx

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
			if((p==2).and.(sort_grid(i,j)==1)) sort_grid(i,j) = 2
			! this sets inner zone (plasma) to "2"

		enddo
		enddo

	endif

	continue

end subroutine set_sort_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine eq_bound(n,x,f)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, rmajor
	use triangularity, only : R_P, z_P, tri_type, a_elps, b_elps

	implicit none

	integer :: n
	real(kind=dkind), dimension(1:n) :: x,f

		if(tri_type==0) then

		f(1) =  x(2)/b_elps**2 * (R_P - x(1))  -  &
				(x(1)-rmajor)/a_elps**2 * (z_P - x(2))

		f(2) = ( (x(1) - rmajor)/a_elps )**2 + ( x(2)/b_elps )**2 - 1.d0

	endif

	continue

end subroutine eq_bound

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine eq_bound2(x,f,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, rmajor
	use triangularity, only : R_P, z_P, tri_type, a_elps, b_elps

	implicit none

	integer :: n
	real(kind=dkind), dimension(1:n) :: x,f

	f(1) =  x(2)/b_elps**2 * (R_P - x(1))  -  &
			(x(1)-rmajor)/a_elps**2 * (z_P - x(2))

	f(2) = ( (x(1) - rmajor)/a_elps )**2 + ( x(2)/b_elps )**2 - 1.d0

	continue

end subroutine eq_bound2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine pack_grid
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use solver, only : radius_theta

	implicit none

	real(kind=dkind) :: x_L, x_R, z_D, z_U
	real(kind=dkind) :: theta, r
	integer :: i
	integer, parameter :: ncheck = 1000
	real(kind=dkind), dimension(0:ncheck) :: xvec, zvec

	x_L = rmajor
	x_R = 0.d0

	z_D = 1.d10
	z_U = -1.d10

	do i = 0,ncheck

		theta = 2.d0*pi*i/ncheck
		call radius_theta(theta, r, xvec(i), zvec(i))

		x_L = min(x_L, xvec(i))
		x_R = max(x_R, xvec(i))

		z_D = min(z_D, zvec(i))
		z_U = max(z_U, zvec(i))

	enddo

	rcenter = (x_R + x_L)/2.d0
	zcenter = (z_U + z_D)/2.d0

	x_size = (x_R - x_L) * 1.0473d0
	z_size = (z_U - z_D) * 1.0473d0

	grid_type = -(10 + grid_type)

	continue

	return

end subroutine pack_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine usrfun(x,n,m,f,jac)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, rmajor
	use triangularity, only : R_P, z_P, tri_type, a_elps, b_elps
	use solver, only : radius_1_3

	implicit none

	integer :: n,m	! just to keep compatibility with the library
	real(kind=dkind), dimension(1:m) :: x,f
	real(kind=dkind), dimension(1:m,1:m) :: jac
	integer :: zone
	real(kind=dkind) :: theta, ex, ez, r, rprim, rsec, rloc2
	real(kind=dkind) :: alph 
	! if the point is inside->(ex,ez) = alph*(ex,ez)

	if(tri_type==0) then

		! first the function

		f(1) =  x(2)/b_elps**2 * (R_P - x(1))  -  &
				(x(1)-rmajor)/a_elps**2 * (z_P - x(2))

		f(2) = ( (x(1) - rmajor)/a_elps )**2 + ( x(2)/b_elps )**2 - 1.d0

		! then the Jacobian

		jac(1,1) = - x(2)/b_elps**2   - (z_P - x(2))/a_elps**2

		jac(1,2) = (R_P - x(1))/b_elps**2  +  &
				(x(1)-rmajor)/a_elps**2

		jac(2,1) = 2.d0*( (x(1) - rmajor)/a_elps**2 )

		jac(2,2) = 2.d0*( x(2)/b_elps**2 )

	elseif(tri_type==8) then

		call radius_1_3(x(1),x(2),ex,ez,theta,r,zone,alph,rprim,rsec)
		! NOTE: ez=x(2)!

		rloc2 = ex**2+ez**2

		! first the function

		f(1) = ( ez - ex/rloc2*r*rprim ) * (R_P - x(1)) -  &
				( ex + ez/rloc2*r*rprim ) * (z_P - x(2))

!		f(2) = (alph*ex)**2 + (alph*ez)**2 - r**2
		f(2) = ex**2 + ez**2 - r**2

		! then the Jacobian

		jac(1,1) = -(rloc2**2*z_P + ez*(ex*(x(1)-R_P) +  &
					x(2)*(x(2)-z_P)) * rprim**2 +  &
					r * ( ((rmajor-R_P)*(ex**2-ez**2) - 2*ex*ez*z_P) *  &
					rprim + ez*(ex*(x(1)-R_P)+ez*(x(2)-z_P))*rsec))  &
					/rloc2**2

		jac(1,2) = (R_P-rmajor) + ( ex*(x(1)**2 +  &
						rmajor*R_P - x(1)*(R_P+rmajor) +  &
						ez**2-ez*z_P)*rprim**2 + &
					r*(-(2.d0*ex*(rmajor-R_P)*ez +  &
					(ex**2-ez**2)*z_P) * rprim +  &
					ex*( x(1)**2 + rmajor*R_P - x(1)*(R_P+rmajor) +  &
						ez**2-ez*z_P) * rsec ))  &
						/rloc2**2

		jac(2,1) = 2.d0 * ( ex + ez/rloc2*r*rprim )

		jac(2,2) = 2.d0 * ( ez - ex/rloc2*r*rprim )

	endif

	continue

end subroutine usrfun

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine funcv(n,x,f)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!this only computes the value of the function, to be used by newt

	use constant, only : dkind, rmajor
	use triangularity, only : R_P, z_P, tri_type, a_elps, b_elps
	use solver, only : radius_1_3

	implicit none

	integer :: n
	real(kind=dkind), dimension(1:n) :: x,f
	integer :: zone
	real(kind=dkind) :: theta, ex, ez, r, rprim, rsec, rloc2
	real(kind=dkind) :: alph 
	! if the point is inside->(ex,ez) = alph*(ex,ez)

	if(tri_type==0) then

		f(1) =  x(2)/b_elps**2 * (R_P - x(1))  -  &
				(x(1)-rmajor)/a_elps**2 * (z_P - x(2))

		f(2) = ( (x(1) - rmajor)/a_elps )**2 + ( x(2)/b_elps )**2 - 1.d0

	elseif(tri_type==8) then

		call radius_1_3(x(1),x(2),ex,ez,theta,r,zone,alph,rprim,rsec)
		! NOTE: ez=x(2)!

		rloc2 = ex**2+ez**2

		f(1) = ( ez - ex/rloc2*r*rprim ) * (R_P - x(1)) -  &
				( ex + ez/rloc2*r*rprim ) * (z_P - x(2))

		f(2) = ex**2 + ez**2 - r**2

	endif

	continue

end subroutine funcv

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine fdjac(n,x,f,m,jac)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! only computes the Jacobian, f is not needed,
! but passed for library compatibility

	use constant, only : dkind, rmajor
	use triangularity, only : R_P, z_P, tri_type, a_elps, b_elps
	use solver, only : radius_1_3

	implicit none

	integer :: n,m	! just to keep compatibility with the library, n=m
	real(kind=dkind), dimension(1:m) :: x,f
	real(kind=dkind), dimension(1:m,1:m) :: jac
	integer :: zone
	real(kind=dkind) :: theta, ex, ez, r, rprim, rsec, rloc2
	real(kind=dkind) :: alph 

	if(tri_type==0) then

		jac(1,1) = - x(2)/b_elps**2   - (z_P - x(2))/a_elps**2

		jac(1,2) = (R_P - x(1))/b_elps**2  +  &
				(x(1)-rmajor)/a_elps**2

		jac(2,1) = 2.d0*( (x(1) - rmajor)/a_elps**2 )

		jac(2,2) = 2.d0*( x(2)/b_elps**2 )

	elseif(tri_type==8) then

		call radius_1_3(x(1),x(2),ex,ez,theta,r,zone,alph,rprim,rsec)
		! NOTE: ez=x(2)!

		rloc2 = ex**2+ez**2

		jac(1,1) = -(rloc2**2*z_P + ez*(ex*(x(1)-R_P) +  &
					x(2)*(x(2)-z_P)) * rprim**2 +  &
					r * ( ((rmajor-R_P)*(ex**2-ez**2) - 2*ex*ez*z_P) *  &
					rprim + ez*(ex*(x(1)-R_P)+ez*(x(2)-z_P))*rsec))  &
					/rloc2**2

		jac(1,2) = (R_P-rmajor) + ( ex*(x(1)**2 +  &
						rmajor*R_P - x(1)*(R_P+rmajor) +  &
						ez**2-ez*z_P)*rprim**2 + &
					r*(-(2.d0*ex*(rmajor-R_P)*ez +  &
					(ex**2-ez**2)*z_P) * rprim +  &
					ex*( x(1)**2 + rmajor*R_P - x(1)*(R_P+rmajor) +  &
						ez**2-ez*z_P) * rsec ))  &
						/rloc2**2

		jac(2,1) = 2.d0 * ( ex + ez/rloc2*r*rprim )

		jac(2,2) = 2.d0 * ( ez - ex/rloc2*r*rprim )

	elseif(tri_type==13) then

		call radius_1_3(x(1),x(2),ex,ez,theta,r,zone,alph,rprim,rsec)
		! NOTE: ez=x(2)!

		rloc2 = ex**2+ez**2

		jac(1,1) = -(rloc2**2*z_P + ez*(ex*(x(1)-R_P) +  &
					x(2)*(x(2)-z_P)) * rprim**2 +  &
					r * ( ((rmajor-R_P)*(ex**2-ez**2) - 2*ex*ez*z_P) *  &
					rprim + ez*(ex*(x(1)-R_P)+ez*(x(2)-z_P))*rsec))  &
					/rloc2**2

		jac(1,2) = (R_P-rmajor) + ( ex*(x(1)**2 +  &
						rmajor*R_P - x(1)*(R_P+rmajor) +  &
						ez**2-ez*z_P)*rprim**2 + &
					r*(-(2.d0*ex*(rmajor-R_P)*ez +  &
					(ex**2-ez**2)*z_P) * rprim +  &
					ex*( x(1)**2 + rmajor*R_P - x(1)*(R_P+rmajor) +  &
						ez**2-ez*z_P) * rsec ))  &
						/rloc2**2

		jac(2,1) = 2.d0 * ( ex + ez/rloc2*r*rprim )

		jac(2,2) = 2.d0 * ( ez - ex/rloc2*r*rprim )

	endif

	continue

end subroutine fdjac

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine mnewt(ntrial,x,n,tolx,tolf)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	integer, parameter :: np=15
	real(kind=dkind) :: x(np)
	integer :: indx(np)
	integer :: n, ntrial, k, i
	real(kind=dkind) :: errf, errx, fjac(np,np), fvec(np), p(np)
	real(kind=dkind) :: tolf, d, tolx

      do 13  k=1,ntrial
        call usrfun(x,n,np,fvec,fjac)
        errf=0.d0
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
		do i=1,n
			p(i) = -fvec(i)
		enddo
        call ludcmp(fjac,n,np,indx,d)
        call lubksb(fjac,n,np,indx,p)
        errx=0.d0
        do 12 i=1,n
          errx=errx+abs(p(i))
          x(i)=x(i)+p(i)
12      continue
        if(errx.le.tolx)return
13    continue

	  if(k>=ntrial) then
		print*, 'problem in mnewt'
	  endif

      return

end subroutine mnewt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ludcmp(a,n,np,indx,d)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	integer, parameter :: nmax=100
	integer :: j, np, i, n, k, imax
	integer :: indx(np)
	real(kind=dkind), parameter :: tiny=1.0d-20
	real(kind=dkind) :: a(np,np)
	real(kind=dkind) :: vv(nmax)
	real(kind=dkind) :: d
	real(kind=dkind) :: aamax, sum, dum

      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix.'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.)a(n,n)=tiny
      return

end subroutine ludcmp


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine lubksb(a,n,np,indx,b)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none
 
 	integer :: ii, ll, j, np, i, n, k, imax
	real(kind=dkind) :: a(np,np)
	real(kind=dkind) :: b(np)
	real(kind=dkind) :: aamax, sum, dum
	integer :: indx(np)


      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return

end subroutine lubksb

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine newt(maxits,x,n,check,tolx,tolf)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	INTEGER n,nn,NP,MAXITS
	LOGICAL check
	REAL(kind=dkind) :: x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
	PARAMETER  (NP=40,TOLMIN=1.d-6,STPMX=100.d0)
	COMMON /newtv/ fvec(NP),nn
	SAVE /newtv/
!	CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
	INTEGER i,its,j,indx(NP)
	REAL(kind=dkind) d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),  &
					p(NP),xold(NP),fmin
	EXTERNAL fmin

	nn=n
	f=fmin(x)
	test=0.d0

	do 11 i=1,n
	if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
	if(test.lt..01d0*TOLF)return
	sum=0.d0
	do 12 i=1,n
	sum=sum+x(i)**2
12    continue
	stpmax=STPMX*max(sqrt(sum),real(n,dkind))
	do 21 its=1,MAXITS
	call fdjac(n,x,fvec,NP,fjac)
	do 14 i=1,n
	  sum=0.d0
	  do 13 j=1,n
		sum=sum+fjac(j,i)*fvec(j)
13        continue
	  g(i)=sum
14      continue
	do 15 i=1,n
	  xold(i)=x(i)
15      continue
	fold=f
	do 16 i=1,n
	  p(i)=-fvec(i)
16      continue
	call ludcmp(fjac,n,NP,indx,d)
	call lubksb(fjac,n,NP,indx,p)
	call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
	test=0.d0
	do 17 i=1,n
	  if(abs(fvec(i)).gt.test)test=abs(fvec(i))
17      continue
	if(test.lt.TOLF)then
	  check=.false.
	  return
	endif
	if(check)then
	  test=0.d0
	  den=max(f,.5d0*n)
	  do 18 i=1,n
		temp=abs(g(i))*max(abs(x(i)),1.d0)/den
		if(temp.gt.test)test=temp
18        continue
	  if(test.lt.TOLMIN)then
		check=.true.
	  else
		check=.false.
	  endif
	  return
	endif
	test=0.d0
	do 19 i=1,n
	  temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
	  if(temp.gt.test)test=temp
	19      continue
	if(test.lt.TOLX)return
	21    continue

	print*, 'MAXITS exceeded in newt'

end subroutine newt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	INTEGER n
	LOGICAL check
	REAL(kind=dkind) f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
	PARAMETER (ALF=1.d-4,TOLX=1.d-7)
	EXTERNAL func
!	CU    USES func
	INTEGER i
	REAL(kind=dkind) a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,  &
					test,tmplam

	check=.false.
	sum=0.d0
	do 11 i=1,n
	sum=sum+p(i)*p(i)
11    continue
	sum=sqrt(sum)
	if(sum.gt.stpmax)then
	do 12 i=1,n
	  p(i)=p(i)*stpmax/sum
12      continue
	endif
	slope=0.d0
	do 13 i=1,n
	slope=slope+g(i)*p(i)
13    continue
	test=0.d0
	do 14 i=1,n
	temp=abs(p(i))/max(abs(xold(i)),1.d0)
	if(temp.gt.test)test=temp
14    continue
	alamin=TOLX/test
	alam=1.d0
	1     continue
	do 15 i=1,n
	  x(i)=xold(i)+alam*p(i)
15      continue
	f=func(x)
	if(alam.lt.alamin)then
	  do 16 i=1,n
		x(i)=xold(i)
16        continue
	  check=.true.
	  return
	else if(f.le.fold+ALF*alam*slope)then
	  return
	else
	  if(alam.eq.1.d0)then
		tmplam=-slope/(2.d0*(f-fold-slope))
	  else
		rhs1=f-fold-alam*slope
		rhs2=f2-fold2-alam2*slope
		a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
		b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
		if(a.eq.0.d0)then
		  tmplam=-slope/(2.d0*b)
		else
		  disc=b*b-3.d0*a*slope
		  if(disc.lt.0.d0) pause 'roundoff problem in lnsrch'
		  tmplam=(-b+sqrt(disc))/(3.d0*a)
		endif
		if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
	  endif
	endif
	alam2=alam
	f2=f
	fold2=fold
	alam=max(tmplam,.1d0*alam)
	goto 1

end subroutine lnsrch

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function fmin(x)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	INTEGER n,NP
	REAL(kind=dkind) ::  fmin,x(*),fvec
	PARAMETER (NP=40)
	COMMON /newtv/ fvec(NP),n
	SAVE /newtv/

!	CU    USES funcv
	INTEGER i
	REAL(kind=dkind) :: sum
	call funcv(n,x,fvec)
	sum=0.d0
	do 11 i=1,n
	sum=sum+fvec(i)**2
	11    continue
	fmin=0.5d0*sum
	return

end function fmin