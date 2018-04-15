! hooks.f90 to test the rk_stepper

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  subroutine echo_error(pf, level, state)
    use iso_c_binding
    use feval

    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    ! to be implemented
    
  end subroutine echo_error


  subroutine echo_residual(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(pfdp),              pointer :: r(:)

    r => array1(level%R(level%nnodes-1))

    print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
         state%step+1, state%iter, level%index, maxval(abs(r))

  end subroutine echo_residual


end module hooks
