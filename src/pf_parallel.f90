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

!> Module of parallel PFASST routines.

module pf_mod_parallel
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_pfasst

  implicit none
contains

  !>  This is the main interface to pfasst.
  !!  It examines the parameters and decides which subroutine to call
  !!  to execute the code correctly
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf   !<  The complete PFASST structure
    class(pf_encap_t), intent(in   )           :: q0   !<  The initial condition
    real(pfdp),        intent(in   )           :: dt   !<  The time step for each processor
    real(pfdp),        intent(in   )           :: tend !<  The final time of run
    integer,           intent(in   ), optional :: nsteps  !<  The number of time steps
    class(pf_encap_t), intent(inout), optional :: qend    !<  The computed solution at tend

    !  Local variables
    integer :: nproc  !<  Total number of processors
    integer :: nsteps_loc  !<  local number of time steps    
    real(pfdp) :: tend_loc !<  The final time of run

    ! make a local copy of nproc
    nproc = pf%comm%nproc

    !>  Set the number of time steps to do
    !!  The user can either pass in the number of time steps or
    !!  pass in the time step size and length of run
    if (present(nsteps)) then
      nsteps_loc = nsteps
      tend_loc=dble(nsteps_loc*dt)
    else
      nsteps_loc = ceiling(1.0*tend/dt)
      !  Do  sanity check on steps
      if (abs(dble(nsteps_loc)-tend/dt) > dt/100.0) then
        print *,'dt=',dt
        print *,'nsteps=',nsteps_loc
        print *,'tend=',tend
        stop "Invalid nsteps"
      end if
    end if
    pf%state%nsteps = nsteps_loc

    !  do sanity checks on Nproc
    if (mod(nsteps,nproc) > 0) stop "ERROR: nsteps must be multiple of nproc (pf_parallel.f90)."

    ! figure out what routine to call
    if (pf%nlevels .eq. 1) then
       print *,'Calling pipelined SDC with 1 level'
       if (present(qend)) then
          call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc, qend)
       else
          call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc)
       end if
    else
       if (pf%Vcycle) then
          !  Right now, we just call the old routine
          
          ! Use Parareal
          if (pf%use_rk_stepper) then
             if (present(qend)) then
                call pf_parareal_run(pf, q0, dt, tend_loc, nsteps_loc, qend)
             else
                call pf_parareal_run(pf, q0, dt, tend_loc, nsteps_loc)
             end if

          ! Use PFASST
          else
             if (present(qend)) then
                call pf_pfasst_run_old(pf, q0, dt, tend_loc, nsteps_loc, qend)
             else
                call pf_pfasst_run_old(pf, q0, dt, tend_loc, nsteps_loc)
             end if
          end if

       else
          print *,'Calling pipelined SDC with multiple levels'
          if (present(qend)) then
             call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc, qend)
          else
             call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc)
          end if
       end if
    end if
    !  What we would like to do is check for
    !  1.  nlevels==1  and nprocs ==1 -> Serial SDC
    !      Predictor is either spreadQ or nothing
    !      Then we just call a loop on sweeps
    !      Communication is copy
    !  2.  nlevels > 1  and nprocs ==1 -> Serial MLSDC
    !      Predictor is needed to populate levels (or nothing)
    !      Then we just call a loop on MLSDC sweeps
    !      Communication is copy
    !  3.  nlevels == 1  and nprocs > 1 -> Pipelined SDC
    !      Predictor is just like PFASST, but on finest (only) level (or nothing)
    !  4.  nlevels > 1  and nprocs > 1 -> PFASST
  end subroutine pf_pfasst_run
  
    !
  ! Parareal Predictor. Currently the same structure as pf_pfasst_predictor
  !
  !>  Subroutine  to initialize the solution on each processor
  !! Spreads the fine initial condition (F%q0) to all levels and all
  !! nodes.  If we're running with more than one processor, performs
  !! sweeps on the coarsest level.
  !!
  !! No time communication is performed during the predictor since all
  !! procesors can do the work themselves
  !!
  !!  The iteration count is reset to 0, and the status is reset to
  !!  ITERATING.
  subroutine pf_pfasst_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf     !< PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !< Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !< time step

    class(pf_level_t), pointer :: coarse_lev_p
    class(pf_level_t), pointer :: fine_lev_p     !<
    integer                   :: j, k            !<  Loop indices
    integer                   :: level_index     !<  Local variable for looping over levels
    real(pfdp)                :: t0k             !<  Initial time at time step k

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    fine_lev_p => pf%levels(pf%nlevels)
    call spreadq0(fine_lev_p, t0)

    !>  If we are doing a single level, then we only spreadq0 and return
    if (pf%nlevels > 1) then
       do level_index = pf%nlevels, 2, -1
          fine_lev_p => pf%levels(level_index);
          coarse_lev_p => pf%levels(level_index-1)
          call pf_residual(pf, fine_lev_p, dt)
          call restrict_time_space_fas(pf, t0, dt, level_index)
          call save(coarse_lev_p)
          call coarse_lev_p%q0%copy(coarse_lev_p%Q(1))
       end do  !  level_index = pf%nlevels, 2, -1

       if (pf%comm%nproc > 1) then
          coarse_lev_p => pf%levels(1)
          if (pf%Pipeline_G .and. (coarse_lev_p%nsweeps_pred > 1)) then
             !  This is the weird choice.  We burn in without communication, then do extra sweeps
             level_index=1
             coarse_lev_p => pf%levels(level_index)
             !  Burn in the processors without communciation
             do k = 1, pf%rank + 1
                pf%state%iter = -k

                ! Reset initial value (skip on first iteration)
                ! This is equivalent to getting value from previous processor
                if (k > 1) then
                   call coarse_lev_p%q0%copy(coarse_lev_p%qend)
                   ! spread new coarse value if we are not doing PFASST_pred
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(coarse_lev_p, t0)
                   end if
                end if
                !  Do a sweep
                call coarse_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt,1)
             end do  ! k = 1, pf%rank + 1
             ! Now we have mimicked the burn in and we must do pipe-lined sweeps
             do k = 1, coarse_lev_p%nsweeps_pred-1
                pf%state%pstatus = PF_STATUS_ITERATING
                pf%state%status = PF_STATUS_ITERATING
                pf%state%iter =-(pf%rank + 1) -k

                !  Get new initial conditions
                call pf_recv(pf, coarse_lev_p, coarse_lev_p%index*20000+pf%rank+k, .true.)

                !  Do a sweep
                call coarse_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt,1 )
                !  Send forward
                call pf_send(pf, coarse_lev_p,  coarse_lev_p%index*20000+pf%rank+1+k, .false.)

             end do ! k = 1, coarse_lev_p%nsweeps_pred-1
          else  ! (pf%Pipeline_G .and. (coarse_lev_p%nsweeps_pred > 1)) then
             ! Normal predictor burn in
             level_index=1
             coarse_lev_p => pf%levels(level_index)
             do k = 1, pf%rank + 1
                pf%state%iter = -k
                t0k = t0-(pf%rank)*dt + (k-1)*dt

                ! Get new initial value (skip on first iteration)
                if (k > 1) then
                   call coarse_lev_p%q0%copy(coarse_lev_p%qend)
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(coarse_lev_p, t0k)
                   end if
                end if

                call coarse_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt,coarse_lev_p%nsweeps_pred)
             end do
          end if ! (pf%Pipeline_G .and. (coarse_lev_p%nsweeps_pred > 1)) then

          ! Return to fine level...
          call pf_pfasst_v_cycle_post_predictor(pf, t0, dt)

       else 
          ! Single processor... sweep on coarse and return to fine level.
          level_index = 1
          coarse_lev_p => pf%levels(level_index)
          do k = 1, pf%rank + 1
             pf%state%iter = -k
             t0k = t0-(pf%rank)*dt + (k-1)*dt
             call coarse_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt,coarse_lev_p%nsweeps_pred)
          end do

          ! Return to fine level...
          call pf_pfasst_v_cycle_post_predictor(pf, t0, dt)

       end if

    else ! nlevels == 1

       level_index = 1
       coarse_lev_p => pf%levels(level_index)
       do k = 1, pf%rank+1
          pf%state%iter = -k
          t0k = t0 - (pf%rank)*dt + (k-1)*dt
          if (k > 1) then
             call coarse_lev_p%q0%copy(coarse_lev_p%qend)
             if (.not. pf%PFASST_pred) then
                call spreadq0(coarse_lev_p, t0k)
             end if
          end if

          call coarse_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt, coarse_lev_p%nsweeps_pred)
       enddo

       ! print*, 'rank ', pf%rank, ' just did ', k-1, ' iterations in the predictor'

    end if  ! (pf%nlevels > 1) then

    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

  end subroutine pf_pfasst_predictor

  !
  ! Parareal Predictor.
  !
  !>  Subroutine  to initialize the solution on each processor
  !! Spreads the fine initial condition (F%q0) to all levels and all
  !! nodes.  If we're running with more than one processor, performs
  !! sweeps on the coarsest level.
  !!
  !! No time communication is performed during the predictor since all
  !! procesors can do the work themselves
  !!
  !!  The iteration count is reset to 0, and the status is reset to
  !!  ITERATING.
  subroutine pf_parareal_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf     !< PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !< Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !< time step

    class(pf_level_t), pointer :: coarse_lev_p
    class(pf_level_t), pointer :: fine_lev_p     !<
    integer                   :: j, k            !<  Loop indices
    integer                   :: level_index     !<  Local variable for looping over levels
    real(pfdp)                :: t0k             !<  Initial time at time step k

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    fine_lev_p => pf%levels(pf%nlevels)
    call spreadq0(fine_lev_p, t0)

    !>  If we are doing a single level, then we only spreadq0 and return
    if (pf%nlevels > 1) then
       do level_index = pf%nlevels, 2, -1
          fine_lev_p => pf%levels(level_index);
          coarse_lev_p => pf%levels(level_index-1)
          call pf_residual(pf, fine_lev_p, dt)
          call restrict_time_space_fas(pf, t0, dt, level_index)
          call save(coarse_lev_p)
          call coarse_lev_p%q0%copy(coarse_lev_p%Q(1))
       end do  !  level_index = pf%nlevels, 2, -1

       if (pf%comm%nproc > 1) then
          coarse_lev_p => pf%levels(1)
          if (pf%Pipeline_G .and. (coarse_lev_p%nsteps_pred > 1)) then
             !  This is the weird choice.  We burn in without communication, then do extra sweeps
             level_index=1
             coarse_lev_p => pf%levels(level_index)
             !  Burn in the processors without communciation
             do k = 1, pf%rank + 1
                pf%state%iter = -k

                ! Reset initial value (skip on first iteration)
                ! This is equivalent to getting value from previous processor
                if (k > 1) then
                   call coarse_lev_p%q0%copy(coarse_lev_p%qend)
                   ! spread new coarse value if we are not doing PFASST_pred
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(coarse_lev_p, t0)
                   end if
                end if
                !  Do a step
                call coarse_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0, dt, coarse_lev_p%nsteps_rk)
             end do  ! k = 1, pf%rank + 1
             ! Now we have mimicked the burn in and we must do pipe-lined sweeps
             do k = 1, coarse_lev_p%nsteps_pred-1
                pf%state%pstatus = PF_STATUS_ITERATING
                pf%state%status = PF_STATUS_ITERATING
                pf%state%iter =-(pf%rank + 1) -k

                !  Get new initial conditions
                call pf_recv(pf, coarse_lev_p, coarse_lev_p%index*20000+pf%rank+k, .true.)

                !  Do a sweep
                call coarse_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0, dt, coarse_lev_p%nsteps_rk)

                !  Send forward
                call pf_send(pf, coarse_lev_p,  coarse_lev_p%index*20000+pf%rank+1+k, .false.)

             end do ! k = 1, coarse_lev_p%nsweeps_pred-1
          else  ! (pf%Pipeline_G .and. (coarse_lev_p%nsweeps_pred > 1)) then
             ! Normal predictor burn in
             level_index=1
             coarse_lev_p => pf%levels(level_index)
             do k = 1, pf%rank + 1
                pf%state%iter = -k
                t0k = t0-(pf%rank)*dt + (k-1)*dt

                ! Get new initial value (skip on first iteration)
                if (k > 1) then
                   call coarse_lev_p%q0%copy(coarse_lev_p%qend)
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(coarse_lev_p, t0k)
                   end if
                end if
                
                call coarse_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0k, dt, coarse_lev_p%nsteps_rk)

             end do
          end if ! (pf%Pipeline_G .and. (coarse_lev_p%nsweeps_pred > 1)) then

          ! Return to fine level...
          call pf_parareal_v_cycle_post_predictor(pf, t0, dt)

       else 
          ! Single processor... sweep on coarse and return to fine level.
          level_index = 1
          coarse_lev_p => pf%levels(level_index)
          do k = 1, pf%rank + 1
             pf%state%iter = -k
             t0k = t0-(pf%rank)*dt + (k-1)*dt
             
             call coarse_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0k, dt, coarse_lev_p%nsteps_rk)
          end do

          ! Return to fine level...
          call pf_parareal_v_cycle_post_predictor(pf, t0, dt)

       end if

    else ! nlevels == 1

       level_index = 1
       coarse_lev_p => pf%levels(level_index)
       do k = 1, pf%rank+1
          pf%state%iter = -k
          t0k = t0 - (pf%rank)*dt + (k-1)*dt
          if (k > 1) then
             call coarse_lev_p%q0%copy(coarse_lev_p%qend)
             if (.not. pf%PFASST_pred) then
                call spreadq0(coarse_lev_p, t0k)
             end if
          end if

          call coarse_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0k, dt, coarse_lev_p%nsteps_rk)
       enddo

       ! print*, 'rank ', pf%rank, ' just did ', k-1, ' iterations in the predictor'

    end if  ! (pf%nlevels > 1) then

    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

  end subroutine pf_parareal_predictor
  
  !>
  !> Test residuals to determine if the current processor has converged.
  !>
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor hasn't converged yet either,
  !> regardless of the residual.
  !>
  subroutine pf_check_convergence(pf, k, dt, residual, energy, converged)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual, energy
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: k
    logical,           intent(out)   :: converged   !<  True if this processor is done
    real(pfdp)     :: residual1
    converged = .false.

    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    ! Check to see if tolerances are met
    residual1 = pf%levels(pf%nlevels)%residual
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if ( (abs(1.0_pfdp - abs(residual1/residual)) < pf%rel_res_tol) .or. &
            (abs(residual1)                          < pf%abs_res_tol) ) then
          pf%state%status = PF_STATUS_CONVERGED
       end if
    end if
    residual = residual1

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)
    call pf_recv_status(pf, 8000+k)

    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING) &
         pf%state%status = PF_STATUS_ITERATING

    call pf_send_status(pf, 8000+k)
    call call_hooks(pf, 1, PF_POST_CONVERGENCE)

    ! XXX: this ain't so pretty, perhaps we should use the
    ! 'nmoved' thinger to break this cycle if everyone is
    ! done...

    if (pf%state%status == PF_STATUS_CONVERGED) then
       converged = .true.
       return
    end if

    if (0 == pf%comm%nproc) then
       pf%state%status = PF_STATUS_PREDICTOR
       converged = .true.
       return
    end if

  end subroutine pf_check_convergence

  !
  ! Run in parallel using PFASST.
  !
  subroutine pf_pfasst_run_old(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt, tend
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer                   :: j, k
    integer                   :: level_index
    real(pfdp)                :: residual, energy
    integer                   ::  ierror  !<  Warning flag for communication routines

    logical :: converged, qbroadcast
    logical :: did_post_step_hook

    call start_timer(pf, TTOTAL)


    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
    pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66

    residual = -1
    energy   = -1
    did_post_step_hook = .false.

    fine_lev_p => pf%levels(pf%nlevels)
    call fine_lev_p%q0%copy(q0)

    do k = 1, 666666666   !  Loop over blocks of time steps
       !  Check to see if we should do one more hook 
       if (pf%state%status == PF_STATUS_CONVERGED .and. .not. did_post_step_hook) then
         call call_hooks(pf, -1, PF_POST_STEP)
         did_post_step_hook = .true.
         pf%state%itcnt = pf%state%itcnt + pf%state%iter-1
         pf%state%mysteps = pf%state%mysteps + 1
       end if

       ! in block mode, jump to next block if we've reached the max iteration count
       qbroadcast = .false.
       if (pf%state%iter >= pf%niters) then

          if (.not. did_post_step_hook) then
            call call_hooks(pf, -1, PF_POST_STEP)
            pf%state%itcnt = pf%state%itcnt + pf%state%iter-1
            pf%state%mysteps = pf%state%mysteps + 1
          end if
          did_post_step_hook = .false.

          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          if (pf%state%step >= pf%state%nsteps) exit

          pf%state%status = PF_STATUS_PREDICTOR
          qbroadcast = .true.
       end if

       !  Do this when starting a new block, broadcast new initial conditions to all procs
       if (k > 1 .and. qbroadcast) then
          fine_lev_p => pf%levels(pf%nlevels)
          call pf%comm%wait(pf, pf%nlevels,ierror)             !<  make sure everyone is done
          call fine_lev_p%qend%pack(fine_lev_p%send)    !<  Pack away your last solution
          call pf_broadcast(pf, fine_lev_p%send, fine_lev_p%nvars, pf%comm%nproc-1)
          call fine_lev_p%q0%unpack(fine_lev_p%send)    !<  Everyone resets their q0
       end if

       ! predictor, if requested or we are starting new bloc
       if (pf%state%status == PF_STATUS_PREDICTOR) &
            call pf_pfasst_predictor(pf, pf%state%t0, dt)

       !
       ! perform fine sweeps
       !
       pf%state%iter  = pf%state%iter + 1

       call start_timer(pf, TITERATION)

       ! XXX: this if statement is necessary for block mode cycling...
       if (pf%state%status /= PF_STATUS_CONVERGED) then

          fine_lev_p => pf%levels(pf%nlevels)
          call fine_lev_p%ulevel%sweeper%sweep(pf,      pf%nlevels, pf%state%t0, dt, fine_lev_p%nsweeps)

       end if

       !
       ! check convergence, continue with iteration
       !

       call pf_check_convergence(pf, k, dt, residual, energy, converged)

       if (pf%state%step >= pf%state%nsteps) exit

       if (.not. converged) then
          do level_index = 2, pf%nlevels
             fine_lev_p => pf%levels(level_index)
             call pf_post(pf, fine_lev_p, fine_lev_p%index*10000+k)
          end do
          
          if (pf%state%status /= PF_STATUS_CONVERGED) then
             
             fine_lev_p => pf%levels(pf%nlevels)
             call pf_send(pf, fine_lev_p, fine_lev_p%index*10000+k, .false.)
             
             if (pf%nlevels > 1) then
                coarse_lev_p => pf%levels(pf%nlevels-1)
                call restrict_time_space_fas(pf, pf%state%t0, dt, pf%nlevels)
                call save(coarse_lev_p)
             end if
             
          end if
          
          call pf_pfasst_v_cycle(pf, k, pf%state%t0, dt)
          call call_hooks(pf, -1, PF_POST_ITERATION)
          call end_timer(pf, TITERATION)
       end if
    end do  !   Loop on k over blocks of time steps

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       fine_lev_p => pf%levels(pf%nlevels)
       call qend%copy(fine_lev_p%qend)
    end if
  end subroutine pf_pfasst_run_old

  !
  ! Run in parallel using Parareal. Currently the same structure as pf_pfasst_run_old
  !
  subroutine pf_parareal_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt, tend
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer                   :: j, k
    integer                   :: level_index
    real(pfdp)                :: residual, energy
    integer                   ::  ierror  !<  Warning flag for communication routines

    logical :: converged, qbroadcast
    logical :: did_post_step_hook

    call start_timer(pf, TTOTAL)


    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
    pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66

    residual = -1
    energy   = -1
    did_post_step_hook = .false.

    fine_lev_p => pf%levels(pf%nlevels)
    call fine_lev_p%q0%copy(q0)

    do k = 1, 666666666   !  Loop over blocks of time steps
       !  Check to see if we should do one more hook 
       if (pf%state%status == PF_STATUS_CONVERGED .and. .not. did_post_step_hook) then
         call call_hooks(pf, -1, PF_POST_STEP)
         did_post_step_hook = .true.
         pf%state%itcnt = pf%state%itcnt + pf%state%iter-1
         pf%state%mysteps = pf%state%mysteps + 1
       end if

       ! in block mode, jump to next block if we've reached the max iteration count
       qbroadcast = .false.
       if (pf%state%iter >= pf%niters) then

          if (.not. did_post_step_hook) then
            call call_hooks(pf, -1, PF_POST_STEP)
            pf%state%itcnt = pf%state%itcnt + pf%state%iter-1
            pf%state%mysteps = pf%state%mysteps + 1
          end if
          did_post_step_hook = .false.

          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          if (pf%state%step >= pf%state%nsteps) exit

          pf%state%status = PF_STATUS_PREDICTOR
          qbroadcast = .true.
       end if

       !  Do this when starting a new block, broadcast new initial conditions to all procs
       if (k > 1 .and. qbroadcast) then
          fine_lev_p => pf%levels(pf%nlevels)
          call pf%comm%wait(pf, pf%nlevels,ierror)             !<  make sure everyone is done
          call fine_lev_p%qend%pack(fine_lev_p%send)    !<  Pack away your last solution
          call pf_broadcast(pf, fine_lev_p%send, fine_lev_p%nvars, pf%comm%nproc-1)
          call fine_lev_p%q0%unpack(fine_lev_p%send)    !<  Everyone resets their q0
       end if

       ! predictor, if requested or we are starting new bloc
       if (pf%state%status == PF_STATUS_PREDICTOR) then
          call pf_parareal_predictor(pf, pf%state%t0, dt)
          print *, pf%rank, ' is done with the predictor'
       end if

       !
       ! perform fine sweeps
       !
       pf%state%iter  = pf%state%iter + 1

       call start_timer(pf, TITERATION)

       ! XXX: this if statement is necessary for block mode cycling...
       if (pf%state%status /= PF_STATUS_CONVERGED) then

          fine_lev_p => pf%levels(pf%nlevels)
          
          ! decide between sweep and RK steps
          call fine_lev_p%ulevel%stepper%do_n_steps(pf, pf%nlevels, pf%state%t0, dt, fine_lev_p%nsteps_rk)

          print *, pf%rank, ' is done with the ARK step'

       end if

       !
       ! check convergence, continue with iteration
       !

       call pf_check_convergence(pf, k, dt, residual, energy, converged)

       if (pf%state%step >= pf%state%nsteps) exit

       if (.not. converged) then
          do level_index = 2, pf%nlevels
             fine_lev_p => pf%levels(level_index)
             call pf_post(pf, fine_lev_p, fine_lev_p%index*10000+k)
          end do
          
          if (pf%state%status /= PF_STATUS_CONVERGED) then
             
             fine_lev_p => pf%levels(pf%nlevels)
             call pf_send(pf, fine_lev_p, fine_lev_p%index*10000+k, .false.)
             print *, pf%rank, ' is just sent the fine initial condition'
             call fine_lev_p%qend%eprint()
             

             if (pf%nlevels > 1) then
                coarse_lev_p => pf%levels(pf%nlevels-1)
                call restrict_time_space_fas(pf, pf%state%t0, dt, pf%nlevels)
                call save(coarse_lev_p)
             end if
             
          end if
          
          call pf_parareal_v_cycle(pf, k, pf%state%t0, dt)
          
          print *, pf%rank, ' is done with the v-cycle'

          call call_hooks(pf, -1, PF_POST_ITERATION)
          call end_timer(pf, TITERATION)
       end if
    end do  !   Loop on k over blocks of time steps

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       fine_lev_p => pf%levels(pf%nlevels)
       call qend%copy(fine_lev_p%qend)
    end if
  end subroutine pf_parareal_run

  !
  !> Run single level SDC in pipeline fashion
  !
  subroutine pf_pipeline_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt, tend
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend

    class(pf_level_t), pointer :: lev_p  !<  pointer to the one level we are operating on
    integer                   :: j, k
    real(pfdp)                :: residual
    integer                   :: nblocks !<  The number of blocks of steps to do
    integer                   :: nproc   !<  The number of processors being used
    integer                   :: ierror  !<  Warning flag for communication routines

    logical :: converged   !<  True when this processor is converged to residual

    call start_timer(pf, TTOTAL)

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
    pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66
    residual = 1

!    if (pf%nlevels > 1) stop "ERROR: nlevels  must be 1 to run pipeline mode (pf_parallel.f90)"

    !  pointer to fine level on which we will iterate
    lev_p => pf%levels(pf%nlevels)
    call lev_p%q0%copy(q0)

    nproc = pf%comm%nproc
    nblocks = nsteps/nproc
    do k = 1, nblocks   !  Loop over blocks of time steps
       ! print *,'Starting  step=',pf%state%step,'  block k=',k      
       ! Each block will consist of
       !  1.  predictor
       !  2.  A loop until max iterations, or tolerances met
       !      2a.  Recieve
       !      2b.  Do SDC sweep(s)
       !      2c.  Send
       !  3.  Move solution to next block


       !>  When starting a new block, broadcast new initial conditions to all procs
       !>  For initial block, this is done when initial conditions are set
       if (k > 1) then
          if (nproc > 1)  then
             call lev_p%qend%pack(lev_p%send)    !<  Pack away your last solution
             call pf_broadcast(pf, lev_p%send, lev_p%nvars, pf%comm%nproc-1)
             call lev_p%q0%unpack(lev_p%send)    !<  Everyone resets their q0
          else
             call lev_p%q0%copy(lev_p%qend)    !<  Just stick qend in q0
          end if

          !>  Update the step and t0 variables for new block
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          pf%state%status = PF_STATUS_PREDICTOR
          residual = 1
       end if

       !> Call the predictor
       !> Currently the predictor will do nothing but spread q0 to all the nodes
       if (pf%state%status == PF_STATUS_PREDICTOR) then
          if (pf%use_rk_stepper .eqv. .true.) then
             call pf_parareal_predictor(pf, pf%state%t0, dt)
          else
             call pf_pfasst_predictor(pf, pf%state%t0, dt)
          end if
       end if

       !>  Start the loops over SDC sweeps
       pf%state%iter = 0
       converged = .FALSE.
       pf%state%status = PF_STATUS_ITERATING

       if (pf%use_rk_stepper .eqv. .true. .and. pf%niters > 1) then
          pf%niters = 1
          print *, "Number of iterations redefined, now equal to 1"
       end if

       call start_timer(pf, TITERATION)
       do j = 1, pf%niters
          call call_hooks(pf, -1, PF_PRE_ITERATION)

          pf%state%iter = j

          !<  Get new initial condition unless this is the first step or this processor is done
          if (pf%state%pstatus .ne. PF_STATUS_CONVERGED) then
             call pf_recv_status(pf, 8000+k)
             call pf_recv(pf, lev_p, lev_p%index*10000+100*k+pf%state%iter, .true.)
          end if

          !<  Decide between sweep and RK step
          if (pf%use_rk_stepper .eqv. .true.) then
             call lev_p%ulevel%stepper%do_n_steps(pf, pf%nlevels, pf%state%t0, dt, lev_p%nsteps_rk)
          else           
             call lev_p%ulevel%sweeper%sweep(pf,      pf%nlevels, pf%state%t0, dt, lev_p%nsweeps)
          end if

          if (pf%use_rk_stepper .eqv. .true.) then
             converged = .false.
          else
             call pf_check_convergence_pipeline(pf, lev_p%residual, converged)
          end if

          if (pf%state%status .ne. PF_STATUS_CONVERGED) then
             call pf_send(pf, lev_p, lev_p%index*10000+100*k+pf%state%iter, .false.)
             if (converged) then
                pf%state%status = PF_STATUS_CONVERGED
             endif
             call pf_send_status(pf, 8000+k)
          endif

          call call_hooks(pf, -1, PF_POST_ITERATION)

          if (pf%state%status == PF_STATUS_CONVERGED) exit
       end do  !  Loop over the iteration in this bloc

       call end_timer(pf, TITERATION)

    end do

    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       call qend%copy(lev_p%qend)
    end if
  end subroutine pf_pipeline_run

  !>
  !> Test residuals to determine if the current processor has converged.
  !>
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor hasn't converged yet either,
  !> regardless of the residual.
  !>
  subroutine pf_check_convergence_pipeline(pf, residual, converged)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual
    logical,           intent(out)   :: converged   !<  True if this processor is done
    real(pfdp)     :: residual1
    ! converged = .false.

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)

    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    ! Check to see if tolerances are met
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if (abs(residual) < pf%abs_res_tol) converged = .true.
       print*, pf%rank, 'r =', residual
       print*, pf%rank, 'converged = ', converged
    end if

    call call_hooks(pf, 1, PF_POST_CONVERGENCE)

  end subroutine pf_check_convergence_pipeline

  !
  ! In PFASST, After predictor, return to fine level.
  !
  subroutine pf_pfasst_v_cycle_post_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer :: j
    integer :: level_index

    if (pf%nlevels <= 1) return

    do level_index = 2, pf%nlevels-1
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, coarse_lev_p%Finterp)
       !       call fine_lev_p%Q(1)%pack(fine_lev_p%q0)
       call fine_lev_p%q0%copy(fine_lev_p%Q(1))
                 
       call fine_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, fine_lev_p%nsweeps_pred)

    end do
    coarse_lev_p => pf%levels(pf%nlevels-1)
    fine_lev_p => pf%levels(pf%nlevels)
    call interpolate_time_space(pf, t0, dt, pf%nlevels, coarse_lev_p%Finterp)
    call fine_lev_p%q0%copy(fine_lev_p%Q(1))

  end subroutine pf_pfasst_v_cycle_post_predictor

  !
  ! In Parareal, after predictor, return to fine level.
  !
  subroutine pf_parareal_v_cycle_post_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer :: j
    integer :: level_index

    if (pf%nlevels <= 1) return

    do level_index = 2, pf%nlevels-1
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, coarse_lev_p%Finterp)
       !       call fine_lev_p%Q(1)%pack(fine_lev_p%q0)
       call fine_lev_p%q0%copy(fine_lev_p%Q(1))

       call fine_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0, dt, fine_lev_p%nsteps_rk)

    end do
    coarse_lev_p => pf%levels(pf%nlevels-1)
    fine_lev_p => pf%levels(pf%nlevels)
    call interpolate_time_space(pf, t0, dt, pf%nlevels, coarse_lev_p%Finterp)
    call fine_lev_p%q0%copy(fine_lev_p%Q(1))

  end subroutine pf_parareal_v_cycle_post_predictor


  !
  ! In PFASST, execute a V-cycle, starting and ending from the middle level.
  !
  subroutine pf_pfasst_v_cycle(pf, iteration, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer :: j
    integer :: level_index
    
    !  For a single level, just get new initial conditions and return
    if (pf%nlevels == 1) then
       fine_lev_p => pf%levels(1)
       call pf_recv(pf, fine_lev_p, fine_lev_p%index*10000+iteration, .true.)
       return
    end if

    !
    ! down (fine to coarse)
    !
    do level_index = pf%nlevels-1, 2, -1
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
      call fine_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, fine_lev_p%nsweeps)
      call pf_send(pf, fine_lev_p, level_index*10000+iteration, .false.)
      call restrict_time_space_fas(pf, t0, dt, level_index)
      call save(coarse_lev_p)
    end do

    !
    ! bottom  (coarsest level)
    !
    level_index=1
    fine_lev_p => pf%levels(level_index)
    
    if (pf%Pipeline_G) then
       do j = 1, fine_lev_p%nsweeps
          call pf_recv(pf, fine_lev_p, fine_lev_p%index*10000+iteration+j, .true.)
          call fine_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt,1)
          call pf_send(pf, fine_lev_p, fine_lev_p%index*10000+iteration+j, .false.)
       end do
    else
       call pf_recv(pf, fine_lev_p, fine_lev_p%index*10000+iteration, .true.)
       call fine_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, fine_lev_p%nsweeps)
       call pf_send(pf, fine_lev_p, level_index*10000+iteration, .false.)
    endif

    !
    ! up  (coarse to fine)
    !
    do level_index = 2, pf%nlevels
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
      call interpolate_time_space(pf, t0, dt, level_index,coarse_lev_p%Finterp)
      call pf_recv(pf, fine_lev_p, level_index*10000+iteration, .false.)

       if (pf%rank /= 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          call interpolate_q0(pf,fine_lev_p, coarse_lev_p)
       end if

       if (level_index < pf%nlevels) then
          call call_hooks(pf, level_index, PF_PRE_SWEEP)
          ! compute residual
          ! do while residual > tol and j < nswps
          ! assuming residual computed at end of sweep 
          call fine_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, fine_lev_p%nsweeps)
       end if
    end do

  end subroutine pf_pfasst_v_cycle

  !
  ! In Parareal, execute a V-cycle, starting and ending from the middle level.
  !
  subroutine pf_parareal_v_cycle(pf, iteration, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer :: j
    integer :: level_index
    
    !  For a single level, just get new initial conditions and return
    if (pf%nlevels == 1) then
       fine_lev_p => pf%levels(1)
       call pf_recv(pf, fine_lev_p, fine_lev_p%index*10000+iteration, .true.)
       return
    end if

    !
    ! down (fine to coarse)
    !
    do level_index = pf%nlevels-1, 2, -1
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
      call fine_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0, dt, fine_lev_p%nsteps_rk)
      call pf_send(pf, fine_lev_p, level_index*10000+iteration, .false.)
      call restrict_time_space_fas(pf, t0, dt, level_index)
      call save(coarse_lev_p)
    end do

    !
    ! bottom  (coarsest level)
    !
    level_index=1
    fine_lev_p => pf%levels(level_index)    
    call pf_recv(pf, fine_lev_p, fine_lev_p%index*10000+iteration, .true.)
    call fine_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0, dt, fine_lev_p%nsteps_rk)
    call pf_send(pf, fine_lev_p, level_index*10000+iteration, .false.)

    !
    ! up  (coarse to fine)
    !
    do level_index = 2, pf%nlevels
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
      call interpolate_time_space(pf, t0, dt, level_index,coarse_lev_p%Finterp)
      call pf_recv(pf, fine_lev_p, level_index*10000+iteration, .false.)
      print *, pf%rank, ' is just received the fine initial condition'
      call fine_lev_p%q0%eprint()

       if (pf%rank /= 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          call interpolate_q0(pf,fine_lev_p, coarse_lev_p)
       end if

       if (level_index < pf%nlevels) then
          call call_hooks(pf, level_index, PF_PRE_SWEEP)
          ! compute residual
          ! do while residual > tol and j < nswps
          ! assuming residual computed at end of sweep 
          call fine_lev_p%ulevel%stepper%do_n_steps(pf, level_index, t0, dt, fine_lev_p%nsteps_rk)
          
       end if
    end do

  end subroutine pf_parareal_v_cycle

  !
  !> Communication helpers
  !
  !>  Subroutine to post a receive request for a new initial condition to be received after doing some work
  subroutine pf_post(pf, level, tag)

    type(pf_pfasst_t), intent(in)    :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    integer ::  ierror 
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING) then
       call pf%comm%post(pf, level, tag,ierror)
       if (ierror /= 0) then
          print *, pf%rank, 'warning: error during post', ierror
          stop "pf_parallel:pf_post"
       endif
    end if
  end subroutine pf_post

  !>  Subroutine to send this processor's convergence status to the next processor
  subroutine pf_send_status(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer ::  istatus
    integer ::  ierror

    istatus = pf%state%status
    if (pf%rank /= pf%comm%nproc-1) then
       ! print*, pf%rank,  'is sending status'
       call pf%comm%send_status(pf, tag,istatus,ierror)
       if (ierror /= 0) then
          print *, pf%rank, 'warning: error during send_status', ierror
          stop "pf_parallel:pf_send_status"
       endif
    end if
  end subroutine pf_send_status

  !>  Subroutine to receive the convergence status from the previous processor
  subroutine pf_recv_status(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer ::  ierror, istatus
    if (pf%rank /= 0 .and. pf%state%pstatus .ne. PF_STATUS_CONVERGED) then
       ! print*, pf%rank,  'is receiving status'
       call pf%comm%recv_status(pf, tag,istatus,ierror)
       ! print*, pf%rank,  'has received status'
       if (ierror .eq. 0) then
          pf%state%pstatus = istatus
       else
          print *, pf%rank, 'warning: error during recv_status', ierror
          stop "pf_parallel:pf_recv_status"
       endif
    end if
  end subroutine pf_recv_status

  !>  Subroutine to send the solution to the next processor
  subroutine pf_send(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer ::  ierror
    call start_timer(pf, TSEND + level%index - 1)
    if (pf%rank /= pf%comm%nproc-1 &
         .and. pf%state%status == PF_STATUS_ITERATING) then
       ! print*, pf%rank,  'is sending soln'
       call pf%comm%send(pf, level, tag, blocking,ierror)
       if (ierror /= 0) then
          print *, pf%rank, 'warning: error during send', ierror
          stop "pf_parallel:pf_send"
       endif
    end if
    call end_timer(pf, TSEND + level%index - 1)
  end subroutine pf_send

  !>  Subroutine to recieve the solution from the previous processor
  subroutine pf_recv(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer ::  ierror
    call start_timer(pf, TRECEIVE + level%index - 1)
    if (pf%rank /= 0 .and.  pf%state%pstatus == PF_STATUS_ITERATING) then
       ! print*, pf%rank,  'is receiving soln'
       call pf%comm%recv(pf, level,tag, blocking,ierror)
       ! print*, pf%rank,  'has received soln'
       !  Unpack the sent data into the intial condition
       if (ierror .eq. 0) then
          call level%q0%unpack(level%recv)
       else
          print *, pf%rank, 'warning: mpi error during receive', ierror
          stop "pf_parallel:pf_recv"
       endif
    end if
    call end_timer(pf, TRECEIVE + level%index - 1)
  end subroutine pf_recv

  !>  Subroutine to broadcast the initial condition to all processors
  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root
    integer :: ierror
    call start_timer(pf, TBROADCAST)
    call pf%comm%broadcast(pf, y, nvar, root,ierror)
       if (ierror /= 0) then
          print *, pf%rank, 'warning:  error during broadcast', ierror
          stop "pf_parallel:pf_broadcast"
       endif
    call end_timer(pf, TBROADCAST)
  end subroutine pf_broadcast


  !> Save current solution and function value so that future corrections can be computed
  subroutine save(lev)
    class(pf_level_t), intent(inout) :: lev  !<  Level to save on
    
    integer :: m, p
    
    if (lev%Finterp) then
       if (allocated(lev%pFflt)) then
          do m = 1, lev%nnodes
             do p = 1,size(lev%F(1,:))
                call lev%pF(m,p)%copy(lev%F(m,p))
             end do
             call lev%pQ(m)%copy(lev%Q(m))
          end do
       end if
    else
       if (allocated(lev%pQ)) then
          do m = 1, lev%nnodes
             call lev%pQ(m)%copy(lev%Q(m))
          end do
       end if
    end if
  end subroutine save
  
end module pf_mod_parallel
