module flow

  use constant
  implicit none

  ! Maximum value of Mach_phi
  real (kind=dkind) :: mach_phi_max
  ! exponent in the equation for Mach phi
  real (kind=dkind) :: alpha_mphi
  ! constant related to the minimum value of Mach phi
  real (kind=dkind) :: mphi_min

  ! Maximum and edge value of Mach_theta
  real (kind=dkind) :: mach_theta_max
  real (kind=dkind) :: mach_theta_edge
  real (kind=dkind) :: mach_theta_max_initial

  !Maximum value of Mach_theta with numerical input
  real(kind=dkind) :: mach_theta_num

  ! Mach Alfv�n poloidal for super-Alfv�nic equilibria
  real (kind=dkind) :: 	mach_alf_0

  real (kind=dkind) :: M0

  real(kind=dkind) :: t_mth != 0.1d0
  real(kind=dkind) :: w_mth! = 0.35d0
  ! for "t shape" in mach_theta(psi)

  logical :: static_equi
  ! states whether any flow is included in the equilibrium calculation

  integer :: enforce_static = 0
  ! 0 -> no effect, 1 -> static_equi = .true., 2 -> static_equi = .false.

  real(kind=dkind) :: omega_0
  ! for psi_diff boudnary condition (bc_type 21-24)

end module flow
