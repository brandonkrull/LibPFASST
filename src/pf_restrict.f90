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
! Restriction and FAS routines.

!
! Notes:
!
!   2013-04-30 - Matthew Emmett
!
!    The pf_residual subroutine is now called after each SDC sweep,
!     and it computes the '0 to node' integrals and stores them in
!     'F%I' while it is computing the full SDC residual.  Furthermore,
!     these 'F%I' integrals also contain the appropriate tau corrections.
!
!     This means that when computing FAS corrections: the fine
!     integral term is already done for us, and it is already FAS
!     corrected, so we dont't have to "bring down fas corrections"
!     from finer levels.
!
!
!   2013-04-17 - Matthew Emmett
!
!     Time restriction was switched from point injection to polynomial
!     interpolation (ie, using the 'rmat's in each level) so that we
!     can use proper nodes for each level.
!
!     To recover point injection (ie, use copy instead of axpy)
!     properly we should really do some masking trickery with the
!     restriction matrices (rmat).  XXX.
!
!     Finally, perhaps the workspaces should be preallocated (along
!     with interpolation workspaces...).  XXX
!
!>  Module to restrict solutions between pfasst levels and create the FAS tau correction
module pf_mod_restrict
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_hooks
  implicit none
contains



  !> Restrict (in time and space) fine level to coarse and set coarse level FAS correction.
  !>
  !> The coarse function values are re-evaluated after restriction.
  !> Note that even if the number of variables and nodes is the same,
  !> we should still compute the FAS correction since the function
  !> evaluations may be different.
  !
  subroutine restrict_time_space_fas(pf, t0, dt, level_index)
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp),        intent(in)    :: t0            !<  time at beginning of step
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: level_index   !< defines which level to restrict

    !  Local variables
    class(pf_level_t), pointer :: c_lev_ptr    
    class(pf_level_t), pointer :: f_lev_ptr

    integer    :: m

    real(pfdp), allocatable :: c_times(:)
    real(pfdp), allocatable :: f_times(:)
    
    f_lev_ptr => pf%levels(level_index);
    c_lev_ptr => pf%levels(level_index-1)

    call call_hooks(pf, level_index, PF_PRE_RESTRICT_ALL)
    call start_timer(pf, TRESTRICT + level_index - 1)
    !
    ! create workspaces
    !
    
    if (f_lev_ptr%restrict_workspace_allocated .eqv. .false.) then 
       
       ! allocate memory 
       call c_lev_ptr%ulevel%factory%create_array(f_lev_ptr%c_tmp_array, c_lev_ptr%nnodes, &
            c_lev_ptr%index, SDC_KIND_INTEGRAL, c_lev_ptr%nvars, c_lev_ptr%shape)
       call c_lev_ptr%ulevel%factory%create_array(f_lev_ptr%f_int_arrayr, c_lev_ptr%nnodes, &
            c_lev_ptr%index, SDC_KIND_INTEGRAL, c_lev_ptr%nvars, c_lev_ptr%shape)

       ! set the flag
       f_lev_ptr%restrict_workspace_allocated = .true. 

    end if

    allocate(c_times(c_lev_ptr%nnodes))
    allocate(f_times(f_lev_ptr%nnodes))
    !
    ! restrict q's and recompute f's
    !
    c_times = t0 + dt*c_lev_ptr%nodes
    f_times = t0 + dt*f_lev_ptr%nodes

    call restrict_sdc(f_lev_ptr, c_lev_ptr, f_lev_ptr%Q, c_lev_ptr%Q, .false.,f_times)

    !  Recompute the functions
     call c_lev_ptr%ulevel%sweeper%evaluate_all(c_lev_ptr, c_times)

    !
    ! fas correction
    !
    do m = 1, c_lev_ptr%nnodes-1
       call c_lev_ptr%tauQ(m)%setval(0.0_pfdp)
    end do
    if (pf%state%iter >= pf%taui0)  then
       ! compute '0 to node' integral on the coarse level
      call c_lev_ptr%ulevel%sweeper%integrate(c_lev_ptr, c_lev_ptr%Q, &
        c_lev_ptr%F, dt, f_lev_ptr%c_tmp_array)
       ! compute '0 to node' integral on the fine level
      call f_lev_ptr%ulevel%sweeper%integrate(f_lev_ptr, f_lev_ptr%Q, &
        f_lev_ptr%F, dt, f_lev_ptr%I)
       !  put tau in
       if (allocated(f_lev_ptr%tauQ)) then
          do m = 1, f_lev_ptr%nnodes-1
             call f_lev_ptr%I(m)%axpy(1.0_pfdp, f_lev_ptr%tauQ(m))
          end do
       end if

       ! restrict '0 to node' integral on the fine level (which was
       call restrict_sdc(f_lev_ptr, c_lev_ptr, f_lev_ptr%I, f_lev_ptr%f_int_arrayr, .true.,f_times)

      ! compute '0 to node' tau correction
       do m = 1, c_lev_ptr%nnodes-1
          call c_lev_ptr%tauQ(m)%axpy(1.0_pfdp, f_lev_ptr%f_int_arrayr(m))
          call c_lev_ptr%tauQ(m)%axpy(-1.0_pfdp, f_lev_ptr%c_tmp_array(m))
       end do
    end if

    call end_timer(pf, TRESTRICT + level_index - 1)
    call call_hooks(pf, level_index, PF_POST_RESTRICT_ALL)

    deallocate(c_times)
    deallocate(f_times)
  end subroutine restrict_time_space_fas


    !> Restrict (in time and space) f_sol_array  to c_sol_array
  !! Depending on the flag INTEGRAL, we may be restricting solutions, or integrals of F
  subroutine restrict_sdc(f_lev_ptr, c_lev_ptr, f_encap_array, c_encap_array, IS_INTEGRAL,f_time)


    class(pf_level_t),  intent(inout) :: f_lev_ptr   !<   pointer to fine level
    class(pf_level_t),  intent(inout) :: c_lev_ptr   !<   pointer to coarse level
    class(pf_encap_t),  intent(inout) :: f_encap_array(:)   !< array of fine level data to be restricted
    class(pf_encap_t),  intent(inout) :: c_encap_array(:)   !< array of coarse level data to be computed
    logical,            intent(in)    :: IS_INTEGRAL       !< flag determines if it is integral data being restricted
    real(pfdp),         intent(in) :: f_time(:)             !< time at the fine nodes

    class(pf_encap_t), allocatable :: f_encap_array_c(:)  !<  fine solution restricted in space only
    integer :: m
    integer :: f_nnodes

    f_nnodes = f_lev_ptr%nnodes
    !>  create works space
    if (IS_INTEGRAL) then   ! Restriction of integrals 

      call c_lev_ptr%ulevel%factory%create_array(f_encap_array_c, f_nnodes-1, c_lev_ptr%index, SDC_KIND_INTEGRAL, &
           c_lev_ptr%nvars, c_lev_ptr%shape)
    else    ! Restriction of solutions
       call c_lev_ptr%ulevel%factory%create_array(f_encap_array_c, f_nnodes, c_lev_ptr%index, &
            SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    end if

    !>  do the restriction
    if (IS_INTEGRAL) then   ! Restriction of integrals 
       !  spatial restriction
       do m = 1, f_nnodes-1
          call f_lev_ptr%ulevel%restrict(f_lev_ptr, c_lev_ptr, f_encap_array(m), f_encap_array_c(m), f_time(m))
       end do

       ! temporal restriction
       ! when restricting '0 to node' integral terms, skip the first entry since it is zero
       call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat(2:,2:), f_encap_array_c)

    else
       !  spatial restriction
       do m = 1, f_nnodes
          call f_lev_ptr%ulevel%restrict(f_lev_ptr, c_lev_ptr, f_encap_array(m), f_encap_array_c(m), f_time(m))
       end do! temporal restriction
       call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat, f_encap_array_c)
    end if

    !>  destroy workspace
    if (IS_INTEGRAL) then   ! Restriction of integrals 
       call c_lev_ptr%ulevel%factory%destroy_array(f_encap_array_c, f_nnodes-1, c_lev_ptr%index, &
         SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    else
       call c_lev_ptr%ulevel%factory%destroy_array(f_encap_array_c, f_nnodes, c_lev_ptr%index, &
            SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    end if

  end subroutine restrict_sdc

  
  
  !> Apply an matrix (tmat or rmat) to src to get dst.
  subroutine pf_apply_mat(dst, a, mat, src, zero)
    class(pf_encap_t), intent(inout) :: dst(:)
    real(pfdp),        intent(in)    :: a, mat(:, :)
    class(pf_encap_t), intent(in)    :: src(:)
    logical,           intent(in), optional :: zero
    
    logical :: lzero
    integer :: n, m, i, j

    lzero = .true.; if (present(zero)) lzero = zero

    n = size(mat, dim=1)
    m = size(mat, dim=2)

    ! XXX: test for nan's in matrices...

    do i = 1, n
       if (lzero) call dst(i)%setval(0.0_pfdp)
       do j = 1, m
          call dst(i)%axpy(a * mat(i, j), src(j))
       end do
    end do
  end subroutine pf_apply_mat

end module pf_mod_restrict
