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
!>  Module to define the main parameters, data types, and interfaces in pfasst
!!
module pf_mod_dtype
  use iso_c_binding
  implicit none

  !  static pfasst paramters
!  integer, parameter :: pfdp = c_long_double
  integer, parameter :: pfdp = c_double

  real(pfdp), parameter :: ZERO   = 0.0_pfdp
  real(pfdp), parameter :: ONE    = 1.0_pfdp
  real(pfdp), parameter :: TWO    = 2.0_pfdp
  real(pfdp), parameter :: THREE  = 3.0_pfdp
  real(pfdp), parameter :: HALF   = 0.5_pfdp
 
  integer, parameter :: PF_MAX_HOOKS = 32
  !> Quadrature node types
  integer, parameter :: SDC_GAUSS_LOBATTO   = 1
  integer, parameter :: SDC_GAUSS_RADAU     = 2
  integer, parameter :: SDC_CLENSHAW_CURTIS = 3
  integer, parameter :: SDC_UNIFORM         = 4
  integer, parameter :: SDC_GAUSS_LEGENDRE  = 5
  integer, parameter :: SDC_PROPER_NODES    = 2**8
  integer, parameter :: SDC_COMPOSITE_NODES = 2**9
  integer, parameter :: SDC_NO_LEFT         = 2**10

  !> Types of encaps node types
  integer, parameter :: SDC_KIND_SOL_FEVAL    = 1
  integer, parameter :: SDC_KIND_SOL_NO_FEVAL = 2
  integer, parameter :: SDC_KIND_FEVAL        = 3
  integer, parameter :: SDC_KIND_INTEGRAL     = 4
  integer, parameter :: SDC_KIND_CORRECTION   = 5

  !> States of operatrion
  integer, parameter :: PF_STATUS_ITERATING = 1
  integer, parameter :: PF_STATUS_CONVERGED = 2
  integer, parameter :: PF_STATUS_PREDICTOR = 3

  type, bind(c) :: pf_state_t
    real(pfdp) :: t0  !<  Time at beginning of this time step
    real(pfdp) :: dt  !<  Time step size
    integer :: nsteps   !< total number of time steps
    integer :: cycle    !< deprecated?
    integer :: iter     !< current iteration number
    integer :: step     !< current time step number assigned to processor
    integer :: level    !< which level is currently being operated on
    integer :: hook     !< which hook
    integer :: proc     !< which processor
    integer :: sweep    !< sweep number
    integer :: status   !< status (iterating, converged etc)
    integer :: pstatus  !< previous rank's status
    integer :: itcnt    !< iteration counter
    integer :: mysteps  !< steps I did
  end type pf_state_t

  type :: pf_hook_t
     procedure(pf_hook_p), pointer, nopass :: proc
  end type pf_hook_t

  !<  Defines the base sweeper type
  type, abstract :: pf_sweeper_t
     integer     :: npieces
   contains
     procedure(pf_sweep_p),        deferred :: sweep
     procedure(pf_initialize_p),   deferred :: initialize
     procedure(pf_evaluate_p),     deferred :: evaluate
     procedure(pf_integrate_p),    deferred :: integrate
     procedure(pf_evaluate_all_p), deferred :: evaluate_all
     procedure(pf_residual_p),     deferred :: residual
     procedure(pf_destroy_p),      deferred :: destroy
  end type pf_sweeper_t

  !<  Defines the base stepper type
  type, abstract :: pf_stepper_t
     integer     :: npieces
     integer     :: order
   contains
     procedure(pf_do_n_steps_p),           deferred :: do_n_steps
     procedure(pf_initialize_stepper_p),   deferred :: initialize
     procedure(pf_destroy_stepper_p),      deferred :: destroy
  end type pf_stepper_t

  !<  Defines the base data type 
  type, abstract :: pf_encap_t
   contains
     procedure(pf_encap_setval_p),  deferred :: setval
     procedure(pf_encap_copy_p),    deferred :: copy
     procedure(pf_encap_norm_p),    deferred :: norm
     procedure(pf_encap_pack_p),    deferred :: pack
     procedure(pf_encap_unpack_p),  deferred :: unpack
     procedure(pf_encap_axpy_p),    deferred :: axpy
     procedure(pf_encap_eprint_p),  deferred :: eprint
  end type pf_encap_t

  !<  Generic type for creation and desstruction of objects
  type, abstract :: pf_factory_t
   contains
     procedure(pf_encap_create_single_p),  deferred :: create_single
     procedure(pf_encap_create_array_p),   deferred :: create_array
     procedure(pf_encap_destroy_single_p), deferred :: destroy_single
     procedure(pf_encap_destroy_array_p),  deferred :: destroy_array
  end type pf_factory_t

  !<  The user level which is inherited  to include problem dependent stuff
  type, abstract :: pf_user_level_t
     class(pf_factory_t), allocatable :: factory
     class(pf_sweeper_t), allocatable :: sweeper
     class(pf_stepper_t), allocatable :: stepper
   contains
     procedure(pf_transfer_p), deferred :: restrict
     procedure(pf_transfer_p), deferred :: interpolate
  end type pf_user_level_t

  !>  Data type of a PFASST level
  type :: pf_level_t
     integer  :: index        = -1   !< level number (1 is the coarsest)
     integer  :: nvars        = -1   !< number of variables (dofs)
     integer  :: nnodes       = -1   !< number of sdc nodes
     integer  :: nsteps_rk    = -1   !< number of rk steps to perform
     integer  :: nsteps_pred  = -1   !< number of times do_n_steps is called in predictor
     integer  :: nsweeps      =  1   !< number of sdc sweeps to perform
     integer  :: nsweeps_pred =  1   !< number of sdc sweeps to perform (predictor)
     logical     :: Finterp = .false.   !< interpolate functions instead of solutions

     real(pfdp)  :: residual            !< holds the user defined residual
     real(pfdp)  :: residual0           !< residual at beginning of PFASST call
                                        ! used for relative residual tolerance calculation

     class(pf_user_level_t), allocatable :: ulevel  !<  user customized level info

     !>  Simple data storage at each level
     real(pfdp), allocatable :: &
          send(:),    &                 !< send buffer
          recv(:),    &                 !< recv buffer
          nodes(:),   &                 !< list of SDC nodes
          t_sdc(:),   &                 !< time at the SDC nodes
          qmat(:,:),  &                 !< spectral integration matrix (0 to node)
          qmatFE(:,:),  &               !< Forward Euler integration matrix (0 to node)
          qmatBE(:,:),  &               !< Backward Euler matrix (0 to node)
          LUmat(:,:), &                 !< LU factorization (replaces BE matrix in Q form)
          s0mat(:,:), &                 !< integration matrix (node to node)
          rmat(:,:),  &                 !< time restriction matrix
          tmat(:,:)                     !< time interpolation matrix

     integer, allocatable :: &
          nflags(:)                     !< sdc node flags

     !>  Solution variable storage
     class(pf_encap_t), allocatable :: &
          Q(:),     &           !< solution at sdc nodes
          pQ(:),    &           !< unknowns at sdc nodes, previous sweep
          R(:),     &           !< full residuals
          I(:),     &           !< 0 to node integrals
          Fflt(:),  &           !< functions values at sdc nodes (flat)
          tauQ(:),  &           !< fas correction in Q form
          pFflt(:), &           !< functions at sdc nodes, previous sweep (flat)
          q0,       &           !< initial condition 
          qend                  !< solution at end time

     !>  Function  storage
     class(pf_encap_t), pointer :: &
          F(:,:), &                     !< functions values at sdc nodes
          pF(:,:)                       !< functions at sdc nodes, previous sweep

     !>  Interpolation and restriction data structures
     logical :: interp_workspace_allocated = .false.
     class(pf_encap_t), allocatable :: &
          f_delta(:),  &   ! delta on the fine level
          cf_delta(:), &   ! delta fine in space and coarse in time
          c_delta(:)       ! delta on the coarse level
     
     logical :: interp_q0_workspace_allocated = .false.
     class(pf_encap_t), allocatable :: &
          c_delta_q0, &    !<  coarse correction
          f_delta_q0, &    !<  fine correction
          c_q0,       &    !<  coarse initial condition
          f_q0             !<  fine  initial condition

     logical :: restrict_workspace_allocated = .false.
     class(pf_encap_t), allocatable :: &
          c_tmp_array(:), &  ! coarse integral of coarse function values
          f_int_arrayr(:)    ! coarse integral of restricted fine function values

     integer, allocatable :: shape(:)   !< user defined shape array

     logical :: allocated = .false.
  end type pf_level_t

  !>  Data type to define the communicator
  type :: pf_comm_t
     integer :: nproc = -1              ! total number of processors

     integer :: comm = -1               ! communicator
     integer, pointer :: &
          recvreq(:), &                 ! receive requests (indexed by level)
          sendreq(:)                    ! send requests (indexed by level)
     integer :: statreq                 ! status send request

     ! fakie, needs modernization
     type(c_ptr), pointer :: pfs(:)     ! pfasst objects (indexed by rank)
     type(c_ptr), pointer :: pfpth(:,:)

     !> Procedure interfaces
     procedure(pf_post_p),        pointer, nopass :: post
     procedure(pf_recv_p),        pointer, nopass :: recv
     procedure(pf_recv_status_p), pointer, nopass :: recv_status
     procedure(pf_send_p),        pointer, nopass :: send
     procedure(pf_send_status_p), pointer, nopass :: send_status
     procedure(pf_wait_p),        pointer, nopass :: wait
     procedure(pf_broadcast_p),   pointer, nopass :: broadcast
  end type pf_comm_t

  !>  The main data type which includes pretty much everything
  type :: pf_pfasst_t
     !>  Parameters
     integer :: nlevels = -1            !< number of pfasst levels
     integer :: niters  = 5             !< number of PFASST iterations to do
     integer :: rank    = -1            !< rank of current processor
     integer :: qtype   = SDC_GAUSS_LOBATTO  !< type of nodes

     real(pfdp) :: abs_res_tol = 0.d0   !<  absolute convergence tolerance
     real(pfdp) :: rel_res_tol = 0.d0   !<  relative convergence tolerance

     !>  run options  (should be set before pfasst_run is called)
     logical :: Pipeline_G =  .false.   !<  decides if coarsest level sweeps are pipelined
     logical :: PFASST_pred = .false.   !<  decides if the PFASST type predictor is used
     logical :: Vcycle = .true.         !<  decides if Vcycles are done
     logical :: do_predictor = .true.   !<  decides if the predictor is done

     !> RK and Parareal options
     logical :: use_rk_stepper = .false. !< decides if RK steps are used instead of the sweeps

     integer     :: taui0 = -999999     !< iteration cutoff for tau inclusion

     !> pf objects
     type(pf_state_t), allocatable :: state   !<  Describes where in the algorithm proc is
     type(pf_level_t), allocatable :: levels(:) !< Holds the levels
     type(pf_comm_t),  pointer :: comm    !< Points to communicator

     !> hooks variables
     type(pf_hook_t), allocatable :: hooks(:,:,:)  !<  Holds the hooks
     integer,  allocatable :: nhooks(:,:)   !<  Holds the number hooks

     !> timing variables
     logical    :: echo_timings  = .false.
     integer :: timers(100)   = 0
     integer :: runtimes(100) = 0

     !> misc
     logical :: debug = .false.
     character(512) :: outdir

  end type pf_pfasst_t

  !> Interfaces for subroutines
  interface
    !> hooks subroutines
    subroutine pf_hook_p(pf, level, state)
       use iso_c_binding
       import pf_pfasst_t, pf_level_t, pf_state_t
       type(pf_pfasst_t), intent(inout) :: pf
       class(pf_level_t), intent(inout) :: level
       type(pf_state_t),  intent(in)    :: state
     end subroutine pf_hook_p
     
     !> SDC sweeper subroutines
     subroutine pf_sweep_p(this, pf, level_index, t0, dt,nsweeps)
       import pf_pfasst_t, pf_sweeper_t, pf_level_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       real(pfdp),          intent(in)    :: dt
       real(pfdp),          intent(in)    :: t0
       integer,             intent(in)    :: level_index
       integer,             intent(in)    :: nsweeps
     end subroutine pf_sweep_p

     subroutine pf_evaluate_p(this, lev, t, m)
       import pf_sweeper_t, pf_level_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: lev
       real(pfdp),          intent(in)    :: t
       integer,      intent(in)    :: m
     end subroutine pf_evaluate_p

     subroutine pf_evaluate_all_p(this, lev, t)
       import pf_sweeper_t, pf_level_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: lev
       real(pfdp),          intent(in)    :: t(:)
     end subroutine pf_evaluate_all_p

     subroutine pf_initialize_p(this, lev)
       import pf_sweeper_t, pf_level_t
       class(pf_sweeper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: lev
     end subroutine pf_initialize_p

     subroutine pf_destroy_sweeper_p(this)
       import pf_sweeper_t
       class(pf_sweeper_t), intent(inout) :: this
     end subroutine pf_destroy_sweeper_p

     subroutine pf_integrate_p(this, lev, qSDC, fSDC, dt, fintSDC)
       import pf_sweeper_t, pf_level_t, pf_encap_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       class(pf_level_t),   intent(in)    :: lev
       class(pf_encap_t),   intent(in)    :: qSDC(:), fSDC(:, :)
       real(pfdp),          intent(in)    :: dt !<  Time step size
       class(pf_encap_t),   intent(inout) :: fintSDC(:)
     end subroutine pf_integrate_p

     subroutine pf_residual_p(this, lev, dt)
       import pf_sweeper_t, pf_level_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: Lev
       real(pfdp),          intent(in)    :: dt !<  Time step size
     end subroutine pf_residual_p

     subroutine pf_destroy_p(this, lev)
       import pf_sweeper_t, pf_level_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: Lev
     end subroutine pf_destroy_p

     !> RK sweeper subroutines
     subroutine pf_do_n_steps_p(this, pf, level_index, t0, big_dt,nsteps_rk)
       import pf_pfasst_t, pf_stepper_t, pf_level_t, pfdp
       class(pf_stepper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       real(pfdp),          intent(in)    :: big_dt !<  Time step size
       real(pfdp),          intent(in)    :: t0
       integer,             intent(in)    :: level_index
       integer,             intent(in)    :: nsteps_rk
     end subroutine pf_do_n_steps_p

     subroutine pf_initialize_stepper_p(this, lev)
       import pf_stepper_t, pf_level_t
       class(pf_stepper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: lev
     end subroutine pf_initialize_stepper_p
     
     subroutine pf_destroy_stepper_p(this, lev)
       import pf_stepper_t, pf_level_t, pfdp
       class(pf_stepper_t), intent(inout) :: this
       class(pf_level_t),   intent(inout) :: Lev
     end subroutine pf_destroy_stepper_p

     !> transfer interfaces used for restriction and interpolation
     subroutine pf_transfer_p(this, levelF, levelG, qF, qG, t)
       import pf_user_level_t, pf_level_t, pf_encap_t, pfdp
       class(pf_user_level_t), intent(inout) :: this
       class(pf_level_t), intent(inout)      :: levelF, levelG
       class(pf_encap_t),   intent(inout)    :: qF, qG
       real(pfdp),          intent(in)       :: t
     end subroutine pf_transfer_p

     !> encapsulation interfaces
     subroutine pf_encap_create_single_p(this, x, level, kind, nvars, shape)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x
       integer,      intent(in   )              :: level, kind, nvars, shape(:)
     end subroutine pf_encap_create_single_p

     subroutine pf_encap_create_array_p(this, x, n, level, kind, nvars, shape)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x(:)
       integer,      intent(in   )              :: n, level, kind, nvars, shape(:)
     end subroutine pf_encap_create_array_p

     subroutine pf_encap_destroy_single_p(this, x, level, kind, nvars, shape)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x
       integer,      intent(in   )              :: level, kind, nvars, shape(:)
     end subroutine pf_encap_destroy_single_p

     subroutine pf_encap_destroy_array_p(this, x, n, level, kind, nvars, shape)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x(:)
       integer,      intent(in   )              :: n, level, kind, nvars, shape(:)
     end subroutine pf_encap_destroy_array_p

     subroutine pf_encap_setval_p(this, val, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout)        :: this
       real(pfdp),        intent(in)           :: val
       integer,    intent(in), optional :: flags
     end subroutine pf_encap_setval_p

     subroutine pf_encap_copy_p(this, src, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout)           :: this
       class(pf_encap_t), intent(in   )           :: src
       integer,    intent(in   ), optional :: flags
     end subroutine pf_encap_copy_p

     function pf_encap_norm_p(this) result (norm)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(in   ) :: this
       real(pfdp) :: norm
     end function pf_encap_norm_p

     subroutine pf_encap_pack_p(this, z)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(in   ) :: this
       real(pfdp),        intent(  out) :: z(:)
     end subroutine pf_encap_pack_p

     subroutine pf_encap_unpack_p(this, z)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout) :: this
       real(pfdp),        intent(in   ) :: z(:)
     end subroutine pf_encap_unpack_p

     subroutine pf_encap_axpy_p(this, a, x, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout)  :: this
       class(pf_encap_t), intent(in   )  :: x
       real(pfdp),  intent(in)           :: a
       integer, intent(in), optional :: flags
     end subroutine pf_encap_axpy_p

     subroutine pf_encap_eprint_p(this)
       import pf_encap_t
       class(pf_encap_t), intent(inout) :: this
     end subroutine pf_encap_eprint_p

     !> communicator interfaces
     subroutine pf_post_p(pf, level, tag,ierror)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(in)    :: pf
       class(pf_level_t), intent(inout) :: level
       integer,    intent(in)    :: tag
       integer,    intent(inout)       :: ierror
     end subroutine pf_post_p

     subroutine pf_recv_p(pf, level, tag, blocking,ierror)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       class(pf_level_t), intent(inout) :: level
       integer,    intent(in)    :: tag
       logical,           intent(in)    :: blocking
       integer,    intent(inout)       :: ierror
     end subroutine pf_recv_p

     subroutine pf_recv_status_p(pf, tag,istatus,ierror)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer,    intent(in)    :: tag
       integer,    intent(inout)       :: istatus
       integer,    intent(inout)       :: ierror
     end subroutine pf_recv_status_p

     subroutine pf_send_p(pf, level, tag, blocking,ierror)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       class(pf_level_t), intent(inout) :: level
       integer,    intent(in)    :: tag
       logical,           intent(in)    :: blocking
       integer,    intent(inout)       :: ierror
     end subroutine pf_send_p

     subroutine pf_send_status_p(pf, tag,istatus,ierror)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer,    intent(in)    :: tag
       integer,    intent(in)       :: istatus
       integer,    intent(inout)       :: ierror
     end subroutine pf_send_status_p

     subroutine pf_wait_p(pf, level,ierror)
       import pf_pfasst_t
       type(pf_pfasst_t), intent(in) :: pf
       integer,    intent(in) :: level
       integer,    intent(inout)       :: ierror
     end subroutine pf_wait_p

     subroutine pf_broadcast_p(pf, y, nvar, root,ierror)
       import pf_pfasst_t, pfdp
       type(pf_pfasst_t), intent(inout) :: pf
       integer,    intent(in)    :: nvar, root
       real(pfdp)  ,      intent(in)    :: y(nvar)
       integer,    intent(inout)       :: ierror
     end subroutine pf_broadcast_p

  end interface

end module pf_mod_dtype
