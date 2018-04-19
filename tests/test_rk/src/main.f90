! test of the RK schemes

program main
  use pf_mod_mpi!, only: mpi_init, mpi_finalize
  use feval
  integer ::  ierror

  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."
  call test_rk()
  call mpi_finalize(ierror)

contains

  subroutine test_rk()
    use pfasst
    use feval
    use hooks
    use pf_mod_mpi

    implicit none

    integer, parameter :: maxlevs = 2

    type(pf_pfasst_t)             :: pf
    type(pf_comm_t)               :: comm
    type(ndarray), allocatable    :: q0
    type(ndarray), pointer        :: ndarray_obj
    integer                       :: nvars(maxlevs), nnodes(maxlevs), nsteps, l, numprocs, mpi_stat
    double precision              :: dt
    class(ad_sweeper_t), pointer  :: ad_sweeper_ptr

    !
    ! initialize pfasst
    !

    nvars  = [ 1, 1 ]   ! number of dofs on the time/space levels
    nnodes = [ 4, 7 ]   ! number of sdc nodes on time/space levels
    nsteps = 0
    dt     = 0.05_pfdp

    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, maxlevs)

    pf%qtype       = SDC_GAUSS_LOBATTO
    pf%niters      = 8
    pf%abs_res_tol = 0
    pf%rel_res_tol = 0
    pf%Pipeline_G  = .false.

    ! test rk stepper
    pf%use_rk_stepper = .true.

    do l = 1, pf%nlevels
       pf%levels(l)%nsweeps   = 1
       pf%levels(l)%nsteps_rk = 1

       pf%levels(l)%nvars  = nvars(maxlevs-pf%nlevels+l)
       pf%levels(l)%nnodes = nnodes(maxlevs-pf%nlevels+l)

       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
       allocate(ad_stepper_t::pf%levels(l)%ulevel%stepper)

       !call setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%nvars)

       ! select the order of the stepper
       if (l == pf%nlevels) then
          pf%levels(l)%ulevel%stepper%order = 4
       else
          pf%levels(l)%ulevel%stepper%order = 2
       end if

       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = pf%levels(l)%nvars       
    end do


    call pf_pfasst_setup(pf)

    do l = 1, pf%nlevels
       
       if (pf%nlevels == 1) then
          pf%levels(l)%nsweeps            = 1
          pf%levels(l)%nsweeps_pred       = 1   
          pf%levels(l)%nsteps_pred        = 1
       else
          if (l > 1) then

             pf%levels(l)%nsweeps         = 1             
             pf%levels(l)%nsweeps_pred    = 1
             pf%levels(l)%nsteps_pred     = 1
          else
             pf%levels(l)%nsweeps         = 1
             if (comm%nproc .eq. 1) then 
                pf%levels(l)%nsweeps_pred = 0
                pf%levels(l)%nsteps_pred  = 0
             else
                pf%levels(l)%nsweeps_pred = 1
                pf%levels(l)%nsteps_pred  = 1
             end if
          end if
       end if
      
    end do

    !
    ! compute initial condition, add hooks, run
    !
    allocate(q0)
    call ndarray_build(q0, [ pf%levels(pf%nlevels)%nvars ])
    call initial(q0)

    call pf_add_hook(pf, pf%nlevels, PF_POST_SWEEP, echo_error)
    call pf_add_hook(pf, -1,         PF_POST_SWEEP, echo_residual)
    
    call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

    nsteps = comm%nproc
    call pf_pfasst_run(pf, q0, dt, tend=0.d0, nsteps=nsteps)

    call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

    if (pf%rank == comm%nproc-1) then
       ndarray_obj => cast_as_ndarray(pf%levels(pf%nlevels)%qend)
       print *, 'final = ', ndarray_obj%flatarray
       call exact(nsteps*dt, q0%flatarray)
       print *, 'final = ', q0%flatarray
    end if

    deallocate(q0%flatarray)
    deallocate(q0%shape)
    deallocate(q0)

    !
    ! cleanup
    !
    call pf_pfasst_destroy(pf)

  end subroutine test_rk

end program

