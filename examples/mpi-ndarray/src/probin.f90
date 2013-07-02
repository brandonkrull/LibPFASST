module probin

  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: two_pi = 6.2831853071795862d0

  integer, parameter :: maxlevs = 3

  integer, parameter :: PROB_AD    = 1
  integer, parameter :: PROB_HEAT  = 2
  integer, parameter :: PROB_VB    = 3
  integer, parameter :: PROB_WAVE  = 4
  integer, parameter :: PROB_KS    = 5
  integer, parameter :: PROB_SHEAR = 6

  integer, save :: problem

  character(len=64), save :: problem_type

  double precision, save :: v      ! advection velocity (PROB_AD only)
  double precision, save :: Lx     ! domain size
  double precision, save :: nu     ! viscosity
  double precision, save :: t0     ! initial time for exact solution (PROB_AD only)
  double precision, save :: sigma  ! initial condition parameter
  double precision, save :: dt     ! time step

  integer, save :: nlevs           ! number of pfasst levels
  integer, save :: nnodes(maxlevs) ! number of nodes
  integer, save :: nvars(maxlevs)  ! number of grid points
  integer, save :: nsteps          ! number of time steps
  integer, save :: niters          ! number of iterations

  character(len=64), save :: output ! directory name for output
  

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename

    integer :: un

    namelist /prbin/ &
         problem_type, output, &
         v, Lx, nu, t0, &
         nlevs, nnodes, nvars, nsteps, niters


    !
    ! defaults
    !

    problem_type = "ad"
    output       = ""

    nlevs   = 2
    nnodes  = [ 2, 3, 5 ]
    nvars   = [ 32, 64, 128 ]
    niters  = 8
    nsteps  = -1

    v       = 1.d0
    Lx      = 1.d0
    nu      = 0.02d0
    sigma   = 0.004d0
    t0      = 0.25d0
    dt      = 0.01d0


    !
    ! read
    !

    un = 66
    open(unit=un, file=filename, status='old', action='read')
    read(unit=un, nml=prbin)
    close(unit=un)


    !
    ! init
    !
    select case (problem_type)
    case ("ad")
       problem = PROB_AD
    case ("heat")
       problem = PROB_HEAT
    case ("burgers")
       problem = PROB_VB
    case ("wave")
       problem = PROB_WAVE
    case ("ks")
       problem = PROB_KS
       nu      = -1.d0
    case ("shear")
       problem = PROB_SHEAR
    end select

  end subroutine probin_init

end module probin
