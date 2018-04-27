!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

!> Module to do interpolation between pfasst levels
module pf_mod_interpolate
  use pf_mod_dtype
  use pf_mod_restrict
  use pf_mod_timer
  use pf_mod_hooks
  use pf_mod_utils
  implicit none
contains

  !> Subroutine to interpolate (in time and space) level_index-1 to level_index
  !! Interpolation is done by interpolating increments.  
  !! The fine function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, level_index, F_INTERP)
    type(pf_pfasst_t), intent(inout),target :: pf      !< main pfasst structure
    real(pfdp),        intent(in)    :: t0             !< time at beginning of time interval
    real(pfdp),        intent(in)    :: dt             !< time step
    integer,           intent(in)    :: level_index    !< defines which level to interpolate to
    logical,           intent(in)    :: F_INTERP !<  Flag, if true, then do interp on f not sol

    !  Local variables
    class(pf_level_t), pointer :: c_lev_ptr   !  Pointer to coarse level
    class(pf_level_t), pointer :: f_lev_ptr   !  Pointer to fine level

    integer    :: m, p
    real(pfdp), allocatable :: c_times(:)   ! coarse level node times
    real(pfdp), allocatable :: f_times(:)   ! fine level node times

    class(pf_encap_t), allocatable :: c_delta(:)    !  Coarse in time and space
    class(pf_encap_t), allocatable :: cf_delta(:)   !  Coarse in time but fine in space

    class(pf_encap_t), allocatable :: Ftmp(:) 

    f_lev_ptr => pf%levels(level_index)   ! fine level
    c_lev_ptr => pf%levels(level_index-1) ! coarse level

    call call_hooks(pf, level_index, PF_PRE_INTERP_ALL)
    call start_timer(pf, TINTERPOLATE + level_index - 1)

    !> create workspaces
    call c_lev_ptr%ulevel%factory%create_array(c_delta,  c_lev_ptr%nnodes, &
      c_lev_ptr%index, SDC_KIND_CORRECTION, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_array(cf_delta, c_lev_ptr%nnodes, &
      f_lev_ptr%index, SDC_KIND_CORRECTION, f_lev_ptr%nvars, f_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_array(Ftmp, f_lev_ptr%nnodes, &
      f_lev_ptr%index, SDC_KIND_CORRECTION, f_lev_ptr%nvars, f_lev_ptr%shape)

    !> set time at coarse and fine nodes
    allocate(c_times(c_lev_ptr%nnodes))
    allocate(f_times(f_lev_ptr%nnodes))

    c_times = t0 + dt*c_lev_ptr%nodes
    f_times = t0 + dt*f_lev_ptr%nodes

    ! needed for amr
    do m = 1, c_lev_ptr%nnodes
       call c_delta(m)%setval(0.0_pfdp)
       call cf_delta(m)%setval(0.0_pfdp)
    end do

    !>  interpolate coarse level correction in space only
    do m = 1, c_lev_ptr%nnodes
       call c_delta(m)%copy(c_lev_ptr%Q(m))
       call c_delta(m)%axpy(-1.0_pfdp, c_lev_ptr%pQ(m))
       call f_lev_ptr%ulevel%interpolate(f_lev_ptr,c_lev_ptr, cf_delta(m), c_delta(m), c_times(m))
    end do  

    !> interpolate corrections in time
    call pf_apply_mat(f_lev_ptr%Q, 1.0_pfdp, f_lev_ptr%tmat, cf_delta, .false.)

    !> either interpolate function values or recompute them
    if (F_INTERP) then         !  Interpolating F
      do p = 1,size(c_lev_ptr%F(1,:))
          do m = 1, c_lev_ptr%nnodes
             call c_delta(m)%setval(0.0_pfdp)
             call cf_delta(m)%setval(0.0_pfdp)
          end do
          ! interpolate coarse corrections  in space
          do m = 1, c_lev_ptr%nnodes
            call c_delta(m)%copy(c_lev_ptr%F(m,p))
            call c_delta(m)%axpy(-1.0_pfdp, c_lev_ptr%pF(m,p))
            call f_lev_ptr%ulevel%interpolate(f_lev_ptr, c_lev_ptr, cf_delta(m), c_delta(m), c_times(m))
         end do

         ! interpolate corrections  in time
         !call pf_apply_mat(f_lev_ptr%F(:,p), 1.0_pfdp, f_lev_ptr%tmat, cf_delta, .false.)

         call pf_apply_mat(Ftmp, 1.0_pfdp, f_lev_ptr%tmat, cf_delta, .false.)
         
         do m = 1, f_lev_ptr%nnodes
            call f_lev_ptr%F(m,p)%copy(Ftmp(m))
         end do

       end do !  Loop on npieces
    else    ! recompute function values
       call f_lev_ptr%ulevel%sweeper%evaluate_all(f_lev_ptr, f_times)
    end if  !  Feval

    !>  reset qend so that it is up to date
    call f_lev_ptr%qend%copy(f_lev_ptr%Q(f_lev_ptr%nnodes))

    !> destroy local data structures
    call c_lev_ptr%ulevel%factory%destroy_array(c_delta,  c_lev_ptr%nnodes, &
      c_lev_ptr%index, SDC_KIND_CORRECTION, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_array(cf_delta, c_lev_ptr%nnodes, &
      f_lev_ptr%index, SDC_KIND_CORRECTION, f_lev_ptr%nvars, f_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_array(Ftmp, f_lev_ptr%nnodes,     &
      f_lev_ptr%index, SDC_KIND_CORRECTION, f_lev_ptr%nvars, f_lev_ptr%shape)


    call end_timer(pf, TINTERPOLATE + f_lev_ptr%index - 1)
    call call_hooks(pf, f_lev_ptr%index, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  !>  Subroutine to update the fine initial condition from coarse increment by spatial interpolation
  subroutine interpolate_q0(pf, f_lev_ptr, c_lev_ptr)

    type(pf_pfasst_t), intent(inout) :: pf          !<  main pfasst structure
    class(pf_level_t),  intent(inout) :: f_lev_ptr  !<  fine level
    class(pf_level_t),  intent(inout) :: c_lev_ptr  !<  coarse level

    class(pf_encap_t), allocatable ::    c_delta    !<  coarse correction
    class(pf_encap_t), allocatable ::    f_delta    !<  fine correction
    class(pf_encap_t), allocatable ::    c_q0       !<  coarse initial condition
    class(pf_encap_t), allocatable ::    f_q0       !<  fine  initial condition

    call call_hooks(pf, f_lev_ptr%index, PF_PRE_INTERP_Q0)
    call start_timer(pf, TINTERPOLATE + f_lev_ptr%index - 1)

    !> create local workspace
    call c_lev_ptr%ulevel%factory%create_single(c_q0,  c_lev_ptr%index, SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_single(f_q0,  f_lev_ptr%index, SDC_KIND_SOL_NO_FEVAL, f_lev_ptr%nvars, f_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%create_single(c_delta, c_lev_ptr%index, SDC_KIND_CORRECTION,   c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_single(f_delta, f_lev_ptr%index, SDC_KIND_CORRECTION,   f_lev_ptr%nvars, f_lev_ptr%shape)

    ! needed for amr
    call f_q0%setval(0.0_pfdp)
    call c_q0%setval(0.0_pfdp)
    call c_delta%setval(0.0_pfdp)
    call f_delta%setval(0.0_pfdp)

    call c_q0%copy(c_lev_ptr%q0)
    call f_q0%copy(f_lev_ptr%q0)

    !>  restrict fine initial data to coarse
    call f_lev_ptr%ulevel%restrict(f_lev_ptr, c_lev_ptr, f_q0, c_delta, pf%state%t0)
    !>  get coarse level correction
    call c_delta%axpy(-1.0_pfdp, c_q0)

    !>  interpolate correction in space
    call f_lev_ptr%ulevel%interpolate(f_lev_ptr, c_lev_ptr, f_delta, c_delta, pf%state%t0)

    !> update fine inital condition
    call f_q0%axpy(-1.0_pfdp, f_delta)
    call f_lev_ptr%q0%copy(f_q0)

    call end_timer(pf, TINTERPOLATE + f_lev_ptr%index - 1)
    call call_hooks(pf, f_lev_ptr%index, PF_POST_INTERP_Q0)

    !> destroy local workspace
    call c_lev_ptr%ulevel%factory%destroy_single(c_q0,  f_lev_ptr%index, SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_single(f_q0,  f_lev_ptr%index, SDC_KIND_SOL_NO_FEVAL, f_lev_ptr%nvars, f_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%destroy_single(c_delta, c_lev_ptr%index, SDC_KIND_CORRECTION,   c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_single(f_delta, f_lev_ptr%index, SDC_KIND_CORRECTION,   f_lev_ptr%nvars, f_lev_ptr%shape)

  end subroutine interpolate_q0

end module pf_mod_interpolate
