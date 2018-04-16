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

    integer, parameter :: maxlevs = 1

    type(pf_pfasst_t)             :: pf
    type(pf_comm_t)               :: comm
    type(ndarray), allocatable    :: q0
    type(ndarray), pointer        :: ndarray_obj
    integer                       :: nvars(maxlevs), nnodes(maxlevs), nsteps, l
    double precision              :: dt
    class(ad_sweeper_t), pointer  :: ad_sweeper_ptr

    !
    ! initialize pfasst
    !

    nvars  = [ 1 ]   ! number of dofs on the time/space levels
    nnodes = [ 4 ]   ! number of sdc nodes on time/space levels
    nsteps = 0
    dt     = 0.05_pfdp

    
    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, maxlevs)

    pf%qtype       = SDC_GAUSS_LOBATTO
    pf%niters      = 3
    pf%abs_res_tol = 1.0D-12    
    pf%rel_res_tol = 1.0D-12

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

       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = pf%levels(l)%nvars

    end do

    if (pf%nlevels > 1) then
       pf%levels(1)%nsweeps      = 2 
       pf%levels(1)%nsweeps_pred = 2
    end if

    call pf_pfasst_setup(pf)

    !
    ! compute initial condition, add hooks, run
    !
    allocate(q0)
    call ndarray_build(q0, [ pf%levels(pf%nlevels)%nvars ])
    call initial(q0)

    call pf_add_hook(pf, pf%nlevels, PF_POST_SWEEP, echo_error)
    call pf_add_hook(pf, -1,         PF_POST_SWEEP, echo_residual)
    
    nsteps = 2*comm%nproc
    call pf_pfasst_run(pf, q0, dt, tend=0.d0, nsteps=nsteps)

    ndarray_obj => cast_as_ndarray(pf%levels(pf%nlevels)%Q(nnodes(1)))
    print *, ndarray_obj%flatarray
    call exact(nsteps*dt, q0%flatarray)
    print *, q0%flatarray

    deallocate(q0%flatarray)
    deallocate(q0%shape)
    deallocate(q0)

    !
    ! cleanup
    !
    call pf_pfasst_destroy(pf)

  end subroutine test_rk

end program

