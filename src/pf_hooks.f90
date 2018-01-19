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
!>  This module enables the calling of user defined routines from various places in the pfasst algorithm
module pf_mod_hooks
  use pf_mod_dtype
  implicit none

  !>  Define hook indices
  integer, parameter :: &
       PF_PRE_PREDICTOR     = 1, &
       PF_POST_PREDICTOR    = 2, &
       PF_PRE_ITERATION     = 3, &
       PF_POST_ITERATION    = 4, &
       PF_PRE_SWEEP         = 5, &
       PF_POST_SWEEP        = 6, &
       PF_PRE_STEP          = 7, &
       PF_POST_STEP         = 8, &
       PF_PRE_INTERP_ALL    = 9, &
       PF_POST_INTERP_ALL   = 10, &
       PF_PRE_INTERP_Q0     = 11, &
       PF_POST_INTERP_Q0    = 12, &
       PF_PRE_RESTRICT_ALL  = 13, &
       PF_POST_RESTRICT_ALL = 14, &
       PF_PRE_CONVERGENCE   = 15, &
       PF_POST_CONVERGENCE  = 16, &
       PF_MAX_HOOK          = 16

  integer, parameter :: &
       PF_HOOK_LOG_ONE  = 1, &
       PF_HOOK_LOG_ALL  = 7, &
       PF_HOOK_LOG_LAST = PF_MAX_HOOK

  !>  Define hook names
  character(len=20), parameter :: hook_names(PF_HOOK_LOG_LAST) = (/ &
       'pre-predictor      ',  &
       'post-predictor     ',  &
       'pre-iteration      ',  &
       'post-iteration     ',  &
       'pre-sweep          ',  &
       'post-sweep         ',  &
       'pre-step           ',  &
       'post-step          ',  &
       'pre-interp-all     ',  &
       'post-interp-all    ',  &
       'pre-interp-q0      ',  &
       'post-interp-q0     ',  &
       'pre-restrict-all   ',  &
       'post-restrict-all  ',  &
       'pre-convergence    ',  &
       'post-convergence   ' /)

contains

  !> Subroutine to add a procedure to the hook on the given level
  subroutine pf_add_hook(pf, level_ind, hook, proc)
    type(pf_pfasst_t), intent(inout) :: pf            !< main pfasst structure
    integer,           intent(in)    :: level_ind     !< which pfasst level to add hook
    integer,           intent(in)    :: hook          !< which hook to add
    procedure(pf_hook_p)             :: proc          !< precudre to call from hook

    integer :: l   !

    if (level_ind == -1) then  ! Do to all levels
       do l = 1, pf%nlevels
          pf%nhooks(l,hook) = pf%nhooks(l,hook) + 1
          pf%hooks(l,hook,pf%nhooks(l,hook))%proc => proc
       end do
    else  ! Do to just level level_ind
       pf%nhooks(level_ind,hook) = pf%nhooks(level_ind,hook) + 1
       pf%hooks(level_ind,hook,pf%nhooks(level_ind,hook))%proc => proc
    end if

  end subroutine pf_add_hook

  !> Subroutine to call hooks associated with the hook and level
  subroutine call_hooks(pf, level_ind, hook)
    use pf_mod_timer
    type(pf_pfasst_t), intent(inout), target :: pf         !< main pfasst structure
    integer,           intent(in)            :: level_ind  !< which pfasst level to call hook
    integer,           intent(in)            :: hook       !< which hook to call

    integer :: i  !<  hook loop index
    integer :: l  !<  level loop index

    call start_timer(pf, THOOKS)

    pf%state%hook = hook

    if (level_ind == -1) then  ! Do to all levels
       do l = 1, pf%nlevels
          do i = 1, pf%nhooks(l,hook)
             call pf%hooks(l,hook,i)%proc(pf, pf%levels(l), pf%state)
          end do
       end do
    else  ! Do to just level level_ind
       do i = 1, pf%nhooks(level_ind,hook)
          call pf%hooks(level_ind,hook,i)%proc(pf, pf%levels(level_ind), pf%state)
       end do
    end if

    call end_timer(pf, THOOKS)
  end subroutine call_hooks

end module pf_mod_hooks
