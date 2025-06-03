module exp_data

use constant

integer, parameter :: data_dim_lim = 2000
integer :: ibreak_p, ibreak_fi, ibreak_th, ibreak_b0, ibreak_pp, ibreak_d,  &
			ibreak_p_e, ibreak_p_i, ibreak_pb, ibreak_He, ibreak_Hi, ibreak_Om_int, ibreak_psi_diff,  &
			ibreak_bPb
integer :: n_H = 201
integer :: n_int_H = 10	! these two can be set up from input later
integer :: n_int_phi = 10

real(kind=dkind) p_iso_data(data_dim_lim,3)	!numerical data for p(psi)
real(kind=dkind) p_e_data(data_dim_lim,3)	!numerical data for p_e(psi)
real(kind=dkind) p_i_data(data_dim_lim,3)	!numerical data for p_i(psi)
real(kind=dkind) omega_data(data_dim_lim,3)	!numerical data for omega(psi)
real(kind=dkind) F_data(data_dim_lim,3)	!numerical data for F(psi)=R0*B0(psi)
real(kind=dkind) d_data(data_dim_lim,3)	!numerical data for D(psi)
real(kind=dkind) mphi_data(data_dim_lim,3)	!numerical data for Mphi(psi)
real(kind=dkind) mtheta_data(data_dim_lim,3)	!numerical data for Mtheta(psi)
real(kind=dkind) psiprim_data(data_dim_lim,3)	!numerical data for bc_type=5 (LDX)
real(kind=dkind) psib_data(data_dim_lim,3)	!numerical data for bc_type=7, 8 (LDX, free-boundary)
real(kind=dkind) big_Psib_data(data_dim_lim,3)	!numerical data for bc_type=7, 8 (free-boundary)
real(kind=dkind), allocatable :: H_e_data(:,:) ! numerical interpolation of H_e(psi)
real(kind=dkind), allocatable :: H_i_data(:,:) ! numerical interpolation of H_i(psi)
real(kind=dkind), allocatable :: Omega_int_data(:,:) ! numerical interpolation of integral of Omega, for H_e, H_i
real(kind=dkind) :: edge_psi_diff_data(1:data_dim_lim,4)

logical ::	numerical_n , numerical_p_iso , numerical_p_e ,  &
			numerical_p_i , numerical_F , numerical_omega,  &
			numerical_psi_diff	!whether to use numerical data

logical :: numerical_psiprim

logical :: numerical_mphi, numerical_mtheta


real(kind=dkind), dimension(:,:), allocatable :: psi_bscoef

real(kind=dkind), dimension(:), allocatable :: p_iso_cscoef, F_cscoef,   &
												 mphi_cscoef, mtheta_cscoef,  &
												 psiprim_cscoef, d_cscoef,  &
												 p_e_cscoef, p_i_cscoef,  &
												 psib_cscoef, H_e_cscoef, H_i_cscoef,  &
												 Omega_int_cscoef, psi_diff_cscoef,  &
												 big_Psib_cscoef
! for IMSL spline routine
! these have been changeed from FLOW to one dimension

integer :: p_iso_ord = 4	! interpolation order
integer :: F_ord = 4		! interpolation order
integer :: mphi_ord = 4		! interpolation order
integer :: mtheta_ord = 4	! interpolation order
integer :: psiprim_ord = 4	! interpolation order
integer :: d_ord = 4		! interpolation order
integer :: p_e_ord = 4	! interpolation order
integer :: p_i_ord = 4	! interpolation order
integer :: psib_ord = 4		! interpolation order
integer :: big_Psib_ord = 4		! interpolation order
integer :: s_ord = 4	! interpolation order for psi
integer :: H_e_ord = 4	! interpolation order
integer :: H_i_ord = 4	! interpolation order
integer :: Omega_int_ord = 4	! interpolation order
integer :: psi_diff_ord = 4	! interpolation order
!integer :: p_iso_ord = 3	! interpolation order


integer :: omega_option		!1=mphi, 2=omega

real(kind=dkind), dimension(:,:), allocatable :: r_data_psi
! to update the interface for the separatrix case (tri_type==13 or tri_type==-2)

real(kind=dkind), dimension(:), allocatable :: xknot_psi, zknot_psi
integer :: nx_FLOW, nz_FLOW

end module exp_data