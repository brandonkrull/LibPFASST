module pf_mod_optimization
  use pfasst
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_mpi
  use pf_mod_ndarray_oc
  use feval
  use solutions
  implicit none

  real(pfdp), parameter, private :: armijoDecrease = 1e-4
  real(pfdp), parameter, private :: minStepSize    = 1e-6
  real(pfdp), parameter, private :: maxStepSize    = 1e6
  real(pfdp), parameter, private :: stepIncFactor  = 10
  real(pfdp), parameter, private :: c2             = 0.9 
  integer,    parameter, private :: maxIterZoom    = 20
 
contains

  subroutine evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objective, L2NormUSq, savedStatesFlat, &
                                L2errorState, LinfErrorState, L2exactState, LinfExactState)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(ndarray_oc), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha
    integer,                   intent(in   ) :: nsteps
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(  out) :: objective, L2NormUSq
    real(pfdp),                intent(inout) :: savedStatesFlat(:,:,:)    !nsteps, nnodes, product(shape)
    real(pfdp),                intent(inout) :: L2errorState, LinfErrorState, L2exactState, LinfExactState
    integer             :: m, step, nnodes, thisstep
    real(pfdp), pointer :: obj(:)
    real(pfdp) :: stepobjective
    logical :: local_predict
    
    nnodes = pf%levels(pf%nlevels)%nnodes
    allocate(obj(nnodes))    
    objective = 0
    
    L2errorState = 0.0
    LinfErrorState = 0.0
    L2exactState = 0.0
    LinfEXactState = 0.0
    
    !if (pf%rank .eq. 0) then
    call initial(q1, 0.0_pfdp, dt) !dt should be given by rank !initial_rd
    !else
    !  q1%yflatarray = 0.0_pfdp
    !end if
    if ((pf%rank .eq. 0) .or. predict) then
      call pf%levels(pf%nlevels)%q0%copy(q1, 1) !do this only for rank 0 or predict = true?
    end if
    
    if ((do_both .eq. 1) .and. (pf%rank .eq. pf%comm%nproc-1)) then
      q1%pflatarray = 0.0_pfdp
      call pf%levels(pf%nlevels)%qend%copy(q1, 2)
    end if
    
    pf%state%itcnt   = 0  ! reset iteration count here, not inside pf_pfasst_block_oc
                          ! how to count skipped sweeps for state?
    pf%state%skippedy = 0
    
    if (pf%rank .eq. 0) print *, ' *********  solve state **************'
    
    do step = 1, nsteps
!       local_predict = .true.
!       if(warmstart == 1 .and. step == 1 .and. predict == .false. ) local_predict = .false.
      thisstep = (step-1)*pf%comm%nproc + pf%rank
!       print *, pf%rank, "thisstep = ", thisstep
!       call pf_pfasst_block_oc(pf, dt, step, .true., 1, step-1)      
      if (do_both .eq. 1 ) then
        call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., flags=0, step=thisstep)            
      else
        call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., flags=1, step=thisstep)      
      end if
     
      ! evaluate objective on current step    
      do m = 1, nnodes
        call objective_function(pf%levels(pf%nlevels)%ulevel%sweeper, pf%levels(pf%nlevels)%Q(m), &
                                    pf%levels(pf%nlevels)%shape, m, obj(m), step)

        ! record state solution
        if (do_both .eq. 1 ) then ! don't save state values, but the adjoint in this case
          call pf%levels(pf%nlevels)%Q(m)%pack(savedStatesFlat(step,m,:), 2)
        else
          call pf%levels(pf%nlevels)%Q(m)%pack(savedStatesFlat(step,m,:), 1)
        end if  
        
        !
      end do
      stepobjective = 0
      do m = 1, pf%levels(pf%nlevels)%nnodes-1
        stepobjective = stepobjective + &
                 (obj(m)+obj(m+1))*(pf%levels(pf%nlevels)%nodes(m+1)-pf%levels(pf%nlevels)%nodes(m))*dt
      end do
      objective = objective + stepobjective
      
      call dump_exact_state(pf%levels(pf%nlevels)%ulevel%sweeper, pf, 0.0_pfdp, dt, 'yex', 'ycomp', 'ydiff', &
                          thisstep, L2errorState, LinfErrorState, L2exactState, LinfExactState)
      
      ! copy qend to q0 for next step -> this neads to be broadcast in general
!       call pf%levels(pf%nlevels)%q0%copy(pf%levels(pf%nlevels)%qend,1)
      if( step < nsteps ) then
        call pf%levels(pf%nlevels)%qend%pack(pf%levels(pf%nlevels)%send, 1)    !<  Pack away your last solution
        call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, pf%comm%nproc-1)
        call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%send, 1)    !<  Everyone resets their q0     
        
        if (do_both .eq. 1 ) then ! we could broadcast here as well instead of setting adjoint back to zero
                                   ! but that would be a different algorithm (sth. like paraexp only on the blocks)
!           q1%pflatarray = 0.0_pfdp
!           call pf%levels(pf%nlevels)%qend%copy(q1, 2)
          call pf%levels(pf%nlevels)%q0%pack(pf%levels(pf%nlevels)%send, 2)    !<  Pack away your last solution
          call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, 0)
          call pf%levels(pf%nlevels)%qend%unpack(pf%levels(pf%nlevels)%send, 2)    !<  Everyone resets their q0    
        end if
      end if
    end do
    
    

    objective = 0.5*objective !0.5 for trapezoidal rule
    call control_L2Q(pf%levels(pf%nlevels)%ulevel%sweeper, dt, pf%levels(pf%nlevels)%nodes, &
                                  pf%levels(pf%nlevels)%shape, L2NormUSq)
    objective = 0.5*objective + 0.5*alpha*L2NormUSq
    deallocate(obj)
  end subroutine evaluate_objective


  subroutine evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStatesFlat, &
                                L2errorState, LinfErrorState, L2exactState, LinfExactState)
    type(pf_pfasst_t),        intent(inout) :: pf
    type(ndarray_oc), target, intent(inout) :: q1
    real(pfdp),               intent(in   ) :: dt
    integer,                  intent(in   ) :: nsteps
    logical,                  intent(in   ) :: predict
    real(pfdp),               intent(inout) :: gradient(:,:,:,:,:), LinftyNormGrad, L2NormGradSq
    real(pfdp),               intent(inout) :: savedStatesFlat(:,:,:)
    real(pfdp),               intent(inout) :: L2errorState, LinfErrorState, L2exactState, LinfExactState

    type(ndarray_oc) :: q_tmp
    integer :: m, step, thisstep, nnodes
    integer :: dest, source, ierror, stat(MPI_STATUS_SIZE)
    real(pfdp), pointer :: tmp(:), terminal(:)
    real(pfdp) :: nodes(pf%levels(pf%nlevels)%nnodes), tend, t

    allocate(terminal(product(pf%levels(pf%nlevels)%shape)))
    allocate(tmp(product(pf%levels(pf%nlevels)%shape)))
      
    call ndarray_oc_build(q_tmp, pf%levels(pf%nlevels)%shape) 
      
    q1%pflatarray = 0.0_pfdp
    terminal = 0.0_pfdp
    
    if(pf%rank==pf%comm%nproc-1) &
      call pf%levels(pf%nlevels)%qend%copy(q1, 2)
    ! do this only for final step or predict = true?
    
    if (do_both .ne. 1 ) pf%state%itcnt   = 0 ! reset iteration count here, not inside pf_pfasst_block_oc
    
    if (pf%rank .eq. 0) print *, '*********  solve adjoint *************'
    
     do step = nsteps, 1, -1
        thisstep = (step-1)*pf%comm%nproc + pf%rank
        
        if (do_both .eq. 1 ) then ! nothing needs to be done
          do m = 1, pf%levels(pf%nlevels)%nnodes
            call pf%levels(pf%nlevels)%Q(m)%unpack(savedStatesFlat(step,m,:), 2)
          end do
!           if(pf%rank /= pf%comm%nproc-1) then
!             source = modulo(pf%rank+1, pf%comm%nproc)
!             call mpi_recv(terminal, pf%levels(pf%nlevels)%mpibuflen, MPI_REAL8, &
!                           source, 1, pf%comm%comm, stat, ierror) ! use step as tag?
! 
!           end if
!           ! apply matrix exponential to compute value at Q(1;2)
!           
!           call q_tmp%unpack(terminal, 2)
!           m = pf%levels(pf%nlevels)%nnodes          
!           call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%copy(q_tmp, 2)
!           call pf%levels(pf%nlevels)%qend%copy(pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes), 2)
!           m = 1          
!           call compute_exp_lap_dt_times_vec(pf%levels(pf%nlevels)%ulevel%sweeper, dt, q_tmp)
!           call pf%levels(pf%nlevels)%Q(1)%axpy(1.0_pfdp, q_tmp, 2)
!           call pf%levels(pf%nlevels)%q0%copy(pf%levels(pf%nlevels)%Q(1), 2)
!        
!           ! send new Q(1;2)
!           if (pf%rank /= 0) then
!             call pf%levels(pf%nlevels)%Q(1)%pack(tmp, 2)
!             dest = modulo(pf%rank-1, pf%comm%nproc)
!             call mpi_send(tmp, pf%levels(pf%nlevels)%mpibuflen, MPI_REAL8, &   !mpi_send(savedAdjoint(1,:),...
!                         dest, 1, pf%comm%comm, stat, ierror)  ! use step as tag?
!           end if
!        
!           tend = (thisstep+1)*dt 
!           ! compute remaining values
! !           if(pf%rank /= pf%comm%nproc-1) then ! can only skip on last if nproc=nsteps!
!             nnodes = pf%levels(pf%nlevels)%nnodes
!             nodes = pf%levels(pf%nlevels)%nodes
!             do m = nnodes-1, 2, -1
!               call q_tmp%unpack(terminal, 2)
!               t = thisstep*dt + dt*(nodes(m)-nodes(1))
!               call compute_exp_lap_dt_times_vec(pf%levels(pf%nlevels)%ulevel%sweeper, &
!                                             tend-t, q_tmp)
!               call pf%levels(pf%nlevels)%Q(m)%axpy(1.0_pfdp, q_tmp, 2)
!             end do
!           end if
          call call_hooks(pf, -1, PF_POST_STEP)
          
          do m = 1, pf%levels(pf%nlevels)%nnodes
            call pf%levels(pf%nlevels)%Q(m)%pack(tmp, 2)
            gradient(step,m,:,:,:) = reshape(tmp, (/pf%levels(pf%nlevels)%shape(1), &
                                           pf%levels(pf%nlevels)%shape(2), pf%levels(pf%nlevels)%shape(3)/))
          end do
          
          ! need to copy here as below, this works only for nproc=nsteps
          ! don't need broadcast but only have to send from rank 0 to rank nproc-1
!           if(step > 1) then
!             if (pf%comm%nproc > 1 ) then
!               if (pf%rank == 0) then
!                 call pf%levels(pf%nlevels)%Q(1)%pack(tmp, 2)
!                 dest = pf%comm%nproc-1
!                 call mpi_send(tmp, pf%levels(pf%nlevels)%mpibuflen, MPI_REAL8, &   !mpi_send(savedAdjoint(1,:),...
!                           dest, step, pf%comm%comm, stat, ierror)  ! use step as tag?
!               else if (pf%rank == pf%comm%nproc-1) then
!                 source = 0
!                 call mpi_recv(terminal, pf%levels(pf%nlevels)%mpibuflen, MPI_REAL8, &
!                           source, step, pf%comm%comm, stat, ierror) ! use step as tag?
!               end if
!             else
!               call pf%levels(pf%nlevels)%Q(1)%pack(terminal, 2)
!             end if
! !             call pf%levels(pf%nlevels)%q0%pack(pf%levels(pf%nlevels)%send, 2)    !<  Pack away your initial condition
! !             call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, 0)
! !             call pf%levels(pf%nlevels)%qend%unpack(pf%levels(pf%nlevels)%send, 2)    !<  Everyone resets their qend
!           
!           end if
          
        else ! if not do_mixed
          ! assign savedStates to [Q(m), 1]
          do m = 1, pf%levels(pf%nlevels)%nnodes
            call pf%levels(pf%nlevels)%Q(m)%unpack(savedStatesFlat(step,m,:), 1)
          end do
          !         call restrict_for_adjoint(pf, (nsteps-step)*dt, dt, 1, step) ! save states at all levels instead of restricting?      
          call restrict_for_adjoint(pf, thisstep*dt, dt, 1, thisstep+1) ! save states at all levels instead of restricting?      
  
          pf%q0_style = 1
          call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., flags=2, step=thisstep) 
          pf%q0_style = 0
        
          do m = 1, pf%levels(pf%nlevels)%nnodes
            call pf%levels(pf%nlevels)%Q(m)%pack(tmp, 2)
            gradient(step,m,:,:,:) = reshape(tmp, (/pf%levels(pf%nlevels)%shape(1), &
                                           pf%levels(pf%nlevels)%shape(2), pf%levels(pf%nlevels)%shape(3)/))
!            call pf%levels(pf%nlevels)%Q(m)%pack(gradient(step, m,:,:), 2)
          end do
                                      
!          call pf%levels(pf%nlevels)%qend%copy(pf%levels(pf%nlevels)%q0, 2)
          if(step > 1) then
            call pf%levels(pf%nlevels)%q0%pack(pf%levels(pf%nlevels)%send, 2)    !<  Pack away your last solution
            call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, 0)
            call pf%levels(pf%nlevels)%qend%unpack(pf%levels(pf%nlevels)%send, 2)    !<  Everyone resets their q0
          end if
        end if ! do_mixed 
     end do

    call construct_gradient(pf%levels(pf%nlevels)%ulevel%sweeper, gradient, pf%levels(pf%nlevels)%nodes, &
                                      LinftyNormGrad, L2NormGradSq)

    call ndarray_oc_destroy(q_tmp)
    deallocate(tmp)
    deallocate(terminal)
     
  end subroutine evaluate_gradient
  
  

 subroutine armijo_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStatesFlat, alpha, &
                         globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall, &
                                L2errorState, LinfErrorState, L2exactState, LinfExactState)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(ndarray_oc), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
    integer,                   intent(in   ) :: nsteps
    integer,                   intent(inout) :: itersState, itersAdjoint
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(in   ) :: searchDir(:,:,:,:,:)  
    real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq
    real(pfdp),                intent(inout) :: stepSize, LinftyNormGrad, L2NormGradSq, gradient(:,:,:,:,:), savedStatesFlat(:,:,:)
    logical,                   intent(inout) :: stepTooSmall
    real(pfdp),                intent(inout) :: L2errorState, LinfErrorState, L2exactState, LinfExactState

       
    real(pfdp) :: globObjNew, directionTimesGradient, globDirXGradNew
    integer    :: l, ierror
    !real(pfdp) :: armijoDecrease ! this should be in a structure global to the module, and set in
    !                             ! something like init_optimization, along with other parameters 
    !real(pfdp) :: minStepSize
    !armijoDecrease = 1e-4
    !minStepSize    = 1e-6

    do
      call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, stepSize)
      !restrict control
      do l = pf%nlevels-1,1,-1
	call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
      end do

      call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStatesFlat, &
                              L2errorState, LinfErrorState, L2exactState, LinfExactState)
!       itersState = itersState + pf%state%itcnt        
      itersState = itersState + pf%state%itcnt - pf%state%skippedy
           
      call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
      if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew

      if (globObjNew < globObj + armijoDecrease*stepSize*globDirXGrad) then
        ! evaluate gradient to be consistent with wolfe_powell_step
        call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStatesFlat, &
                               L2errorState, LinfErrorState, L2exactState, LinfExactState)
        itersAdjoint = itersAdjoint + pf%state%itcnt
        return 
      end if

       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
       ! no need to restrict here, coarse levels get overwritten above
       stepSize = 0.5 * stepSize
       if (stepSize < minStepSize) then
         stepTooSmall = .true.
         if (pf%rank .eq. 0) print *, 'no step found, stopping'
         !call write_control(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
         return
       end if
     end do
   end subroutine armijo_step

   
   subroutine strong_wolfe_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStatesFlat, &
                                 alpha, globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, &
                                 stepSize, stepTooSmall, &
                                 L2errorState, LinfErrorState, L2exactState, LinfExactState)
     type(pf_pfasst_t),         intent(inout) :: pf
     type(ndarray_oc), target,  intent(inout) :: q1
     real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
     integer,                   intent(in   ) :: nsteps
     integer,                   intent(inout) :: itersState, itersAdjoint
     logical,                   intent(in   ) :: predict
     real(pfdp),                intent(in   ) :: searchDir(:,:,:,:,:) 
     real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq
     real(pfdp),                intent(inout) :: stepSize, gradient(:,:,:,:,:), savedStatesFlat(:,:,:) !savedAdjoint(:,:)
     logical,                   intent(inout) :: stepTooSmall
     real(pfdp),                intent(inout) :: L2errorState, LinfErrorState, L2exactState, LinfExactState
     
     integer    :: l, ierror
     real(pfdp) :: low, high
     real(pfdp) :: directionTimesGradient, globObjNew, globDirXGradNew
     logical    :: first
     
     first = .true.
     
     low  = 0.0_pfdp
     high = stepSize
     !high = 1.0_pfdp
     !if (stepSize > high) high = stepSize
     
     do
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, high) !stepSize)
       !restrict control
       do l = pf%nlevels-1,1,-1
       call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
       end do
 !       call dump_control(pf%levels(pf%nlevels)%ctx, pf, 'u1')
 
       !call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStatesFlat, &
                               L2errorState, LinfErrorState, L2exactState, LinfExactState)
       first = .false.
       
!        itersState = itersState + pf%state%itcnt
       itersState = itersState + pf%state%itcnt - pf%state%skippedy

       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
       if(pf%rank == 0) print *, high, 'objectiveNew (L2) = ', globObjNew ! *, stepSize,
 
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -high) !-stepSize)
       
       if (globObjNew < globObj + armijoDecrease*stepSize*globDirXGrad) then
         high = 2.0*high
       else ! Armijo condition not satisfied, exit loop and reduce stepsize
         exit
       end if
     end do
       
     do
       stepSize = 0.5*(low+high)
       
       if (stepSize < minStepSize ) then
         stepTooSmall = .true.
         if (pf%rank .eq. 0) print *, 'no step found (< min), stopping'
         return
       end if
       
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, stepSize)
       do l = pf%nlevels-1,1,-1
       call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
       end do
 
       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStatesFlat, &
                               L2errorState, LinfErrorState, L2exactState, LinfExactState)
       !(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
!        itersState = itersState + pf%state%itcnt
       itersState = itersState + pf%state%itcnt - pf%state%skippedy

       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
       if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
       
       if (globObjNew >= globObj + armijoDecrease*stepSize*globDirXGrad) then
         call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
         high = stepSize
         cycle
       end if
 
       ! now Armijo is satisfied again
       !call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
       call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStatesFlat, &
                              L2errorState, LinfErrorState, L2exactState, LinfExactState)
       itersAdjoint = itersAdjoint + pf%state%itcnt
       directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
       call mpi_allreduce(directionTimesGradient, globDirXGradNew, 1, MPI_REAL8, MPI_SUM, &
                         pf%comm%comm, ierror)          
                          
       if (pf%rank .eq. 0) print *, 'globDirXGradNew = ', globDirXGradNew, '-c2*globDirXGrad', -c2*globDirXGrad
 
           
 !       if (abs(globDirXGradNew) <= -c2*globDirXGrad) then  ! strong Wolfe conditions satisfied
 !         return
 !       end if
    
       if (globDirXGradNew < c2*globDirXGrad) then
 !         if(pf%rank == 0) print *, stepSize, high, low
         low = stepSize
         if (pf%rank .eq. 0) print *, 'set new stepLow;  stepHigh = ', high, 'stepLow = ', low, 'stepSize = ', stepSize 
       elseif (globDirXGradNew > -c2*globDirXGrad) then
         high = stepSize
         if (pf%rank .eq. 0) print *, 'set new stepHigh; stepHigh = ', high, 'stepLow = ', low, 'stepSize = ', stepSize 
       else
         exit
       end if
       
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
     end do
     
   end subroutine strong_wolfe_step
   
end module pf_mod_optimization

