module p_d_profile

  use constant, only: dkind
!  use magnetic
  implicit none

  ! The adiabatic exponent, $\gamma$
  real (kind=dkind) :: gamma
  real (kind=dkind) :: gamma_i	! polytropic index for ions
  real (kind=dkind) :: gamma_e	! polytropic index for electrons

  real (kind=dkind) :: dcenter, d_TF_center
  real (kind=dkind) :: dedge, d_TF_edge

  ! The Pressure P = pedge + pcenter * (psi/psic)^alpha
  ! The Density D = dedge + dcenter * (psi/psic)^alpha_rho
  ! pcenter is the Gas Pressure at the Magnetic Axis
  ! dcenter is the Gas Density at the Magnetic Axis
  real (kind=dkind) :: alpha, alpha_e, alpha_i
  real (kind=dkind) :: alpha_rho

  real (kind=dkind) :: beta_center, beta_e_center, beta_i_center

  real (kind=dkind) :: pcenter, p_e_center, p_i_center
!      pcenter =  beta_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag

  real (kind=dkind) :: pedge, p_e_edge, p_i_edge

  integer :: p_opt, p_e_opt, p_i_opt
  ! option for pressure equation


end module p_d_profile
