!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level, state)
    use feval, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(pfdp) :: yexact(level%mpibuflen)
    real(pfdp), pointer :: y_end(:)
    real(pfdp) :: maxerr,err0

    !>  compute the error at last end point
    y_end => get_array1d(level%qend)
    call exact(state%t0+state%dt, yexact)
    maxerr = maxval(abs(y_end-yexact))

    !>  compute the error at in initial condition point
    y_end => get_array1d(level%q0)
    call exact(state%t0, yexact)
    err0 = maxval(abs(y_end-yexact))
    
    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7,"  err0: ",es14.7," res: ",es18.10e4)', &
         state%step+1, state%iter,level%index, maxerr,err0,level%residual
    call flush
  end subroutine echo_error

  !>  Output the current residual in the solution
  subroutine echo_residual(pf, level, state)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    print '("resid: time: ", f8.4," rank: ",i3.3," step: ",i5.5," iter: ",i3.3," level: ",i1.1," resid: ",es18.10e4)', &
         state%t0+state%dt, pf%rank, state%step+1, state%iter, level%index, level%residual
    call flush
  end subroutine echo_residual

end module hooks
