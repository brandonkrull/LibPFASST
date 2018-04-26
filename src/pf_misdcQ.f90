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

module pf_mod_misdcQ
  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_misdcQ_t
     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
     real(pfdp), allocatable :: dtsdc(:)

     logical                 :: use_LUq = .true.

     class(pf_encap_t), allocatable :: rhs   !> holds rhs for implicit solve
     class(pf_encap_t), allocatable :: S3(:)

   contains 
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: sweep        => misdcQ_sweep
     procedure :: initialize   => misdcQ_initialize
     procedure :: evaluate     => misdcQ_evaluate
     procedure :: integrate    => misdcQ_integrate
     procedure :: evaluate_all => misdcQ_evaluate_all
     procedure :: residual     => misdcQ_residual
     procedure :: destroy      => misdcQ_destroy
     procedure :: misdcQ_destroy
  end type pf_misdcQ_t
  
  interface
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       import pf_misdcQ_t, pf_encap_t, pfdp
       class(pf_misdcQ_t),  intent(inout) :: this
       class(pf_encap_t),   intent(in   ) :: y
       real(pfdp),          intent(in   ) :: t
       integer,             intent(in   ) :: level_index
       class(pf_encap_t),   intent(inout) :: f
       integer,             intent(in   ) :: piece
     end subroutine pf_f_eval_p
      subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_misdcQ_t, pf_encap_t, pfdp
       class(pf_misdcQ_t),  intent(inout) :: this
       class(pf_encap_t),   intent(inout) :: y
       real(pfdp),          intent(in   ) :: t
       real(pfdp),          intent(in   ) :: dtq
       class(pf_encap_t),   intent(in   ) :: rhs
       integer,             intent(in   ) :: level_index
       class(pf_encap_t),   intent(inout) :: f
       integer,             intent(in   ) :: piece
     end subroutine pf_f_comp_p
  end interface
  
contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine misdcQ_sweep(this, pf, level_index, t0, dt, nsweeps)
    use pf_mod_timer
    class(pf_misdcQ_t),   intent(inout) :: this
    type(pf_pfasst_t),    intent(inout), target :: pf
    real(pfdp),           intent(in)    :: dt, t0
    integer,              intent(in)    :: level_index  !>  which level this is
    integer,              intent(in)    :: nsweeps      !>  number of sweeps to do


    integer                        :: m, n
    real(pfdp)                     :: t

    class(pf_level_t), pointer    :: lev
    
    lev => pf%levels(level_index)   !<  Assign level pointer

    call start_timer(pf, TLEVEL+lev%index-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1

       call lev%I(m)%setval(0.0_pfdp)
       call this%S3(m)%setval(0.0d0)
       do n = 1, lev%nnodes
          call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
          call lev%I(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2))
          call lev%I(m)%axpy(dt*lev%qmat(m,n),    lev%F(n,3))
          call this%S3(m)%axpy(dt*this%QtilI(m,n),     lev%F(n,3))
          !  Note we have to leave off the -dt*Qtil here and put it in after f2comp
       end do
       if (allocated(lev%tauQ)) then
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%copy(lev%q0)

    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1),1)
    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2),2)
    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,3),3)

    t = t0
    do m = 1, lev%nnodes-1
       t = t + dt*this%dtsdc(m)

       call this%rhs%setval(0.0_pfdp)
       do n = 1, m
          call this%rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1))  
          call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2))  
       end do
       !  Add the tau term
       call this%rhs%axpy(1.0_pfdp, lev%I(m))
       !  Add the starting value
       call this%rhs%axpy(1.0_pfdp, lev%Q(1))

       call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,2),2)

       !  Now we need to do the final subtraction for the f3 piece
       call this%rhs%copy(Lev%Q(m+1))       
       do n = 1, m
          call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,3))  
       end do

       call this%rhs%axpy(-1.0_pfdp, this%S3(m))

       call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,3),3)
       call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1),1)
       call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,2),2)
    end do
                         
    call pf_residual(pf, lev, dt)
    call lev%qend%copy(lev%Q(lev%nnodes))

    call end_timer(pf, TLEVEL+lev%index-1)

  end subroutine misdcQ_sweep

  ! Evaluate function values
  subroutine misdcQ_evaluate(this, lev, t, m)
    use pf_mod_dtype
    class(pf_misdcQ_t), intent(inout) :: this
    real(pfdp),        intent(in)    :: t
    integer,           intent(in)    :: m
    class(pf_level_t),  intent(inout) :: lev

    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,3),3)
  end subroutine misdcQ_evaluate

  subroutine misdcQ_evaluate_all(this, lev, t)
    class(pf_misdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),       intent(in)    :: t(:)

    call pf_generic_evaluate_all(this, lev, t)
  end subroutine misdcQ_evaluate_all

     
  ! Initialize matrices
  subroutine misdcQ_initialize(this, lev)
    class(pf_misdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, n, nnodes

    this%npieces = 3

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes)) ! S-FE
    allocate(this%QdiffI(nnodes-1,nnodes)) ! S-BE 
    allocate(this%QtilE(nnodes-1,nnodes)) ! S-FE
    allocate(this%QtilI(nnodes-1,nnodes)) ! S-BE
    allocate(this%dtsdc(nnodes-1))
    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    ! Implicit matrix
    if (this%use_LUq) then 
       ! Get the LU
       call myLUq(lev%qmat,lev%LUmat,lev%nnodes,0)
       this%QtilI = lev%LUmat
    else 
       this%QtilI = lev%qmatBE
    end if
    
    !! Explicit matrix
    this%QtilE=lev%qmatFE

    this%QdiffE = lev%qmat-this%QtilE
    this%QdiffI = lev%qmat-this%QtilI

    !>  Make space for rhs
    call lev%ulevel%factory%create_single(this%rhs,         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%create_array( this%S3, lev%nnodes-1, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

  end subroutine misdcQ_initialize

  subroutine misdcQ_destroy(this, lev)
    class(pf_misdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)

    call lev%ulevel%factory%destroy_single(this%rhs,        lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%destroy_array(this%S3, lev%nnodes-1, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    
  end subroutine misdcQ_destroy

  subroutine misdcQ_residual(this, lev, dt)
    !> Compute  Residual
    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t),   intent(inout) :: lev
    real(pfdp),          intent(in   ) :: dt
    call pf_generic_residual(this, lev, dt)
  end subroutine misdcQ_residual

  ! Compute SDC integral
  subroutine misdcQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t),   intent(in)    :: lev
    class(pf_encap_t),   intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),          intent(in)    :: dt
    class(pf_encap_t),   intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine misdcQ_integrate

end module pf_mod_misdcQ
