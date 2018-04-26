! Test feval.f90 for the RK stepper

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_misdcQ
  use pf_mod_rkstepper
  implicit none

  real(pfdp), parameter :: &
       a = 1.0_pfdp, & ! "advection"
       d = 2.0_pfdp, & ! "diffusion"
       r = 3.0_pfdp    ! "reaction"

  ! Define the derived user_lever type
  type, extends(pf_user_level_t) :: ad_level_t

   contains

     procedure :: restrict    => restrict
     procedure :: interpolate => interpolate

  end type ad_level_t

  ! Define the derived sweeper type
  type, extends(pf_misdcQ_t) :: ad_sweeper_t

   contains

     procedure :: f_eval => ad_sweeper_f_eval
     procedure :: f_comp => ad_sweeper_f_comp

  end type ad_sweeper_t

  ! Define the derived stepper type
  type, extends(pf_ark_t) :: ad_stepper_t
     
   contains
       
     procedure :: f_eval     => ad_stepper_f_eval
     procedure :: f_comp     => ad_stepper_f_comp
     procedure :: f_finalize => ad_stepper_f_finalize
     
  end type ad_stepper_t

contains

  ! Helper functions

  function as_ad_sweeper(sweeper) result(ptr)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(ad_sweeper_t), pointer :: ptr
    select type(sweeper)
    type is (ad_sweeper_t)
       ptr => sweeper
    class default
       stop
    end select
  end function as_ad_sweeper

  subroutine destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev
    
    call this%misdcQ_destroy(lev)

  end subroutine destroy

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0

    call exact(0.0_pfdp, q0%flatarray)

  end subroutine initial

  subroutine exact(t, yex)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    yex(1) = dexp((a+d+r)*t)
  end subroutine exact

  ! Sweeper functions

  ! Evaluate the explicit function at y, t.
  subroutine ad_sweeper_f_eval(this, y, t, level_index, f, piece)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece    

    real(pfdp),          pointer       :: yvec(:), fvec(:)
    real(pfdp)                         :: val

    yvec => array1(y)
    fvec => array1(f)

    if (piece == 1) then
       fvec = a*yvec
    else if (piece == 2) then
       fvec = d*yvec
    else if (piece == 3) then
       fvec = r*yvec
    else 
       stop "piece not implemented"   
    end if

  end subroutine ad_sweeper_f_eval

  ! Solve for y and return f2 also.
  subroutine ad_sweeper_f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),          pointer       :: yvec(:), rhsvec(:), fvec(:)
    real(pfdp)                         :: val

    yvec   => array1(y)
    rhsvec => array1(rhs)
    fvec   => array1(f)
    
    if (piece == 2) then
       yvec = rhsvec / (1.0_pfdp - d*dtq)
       fvec = (yvec - rhsvec) / dtq 
    else if (piece == 3) then
       yvec = rhsvec / (1.0_pfdp - r*dtq)
       fvec = (yvec - rhsvec) / dtq
    else if (piece == 1) then       
       stop("Error: Piece 1 is explicit")
    else
       stop "piece not implemented"
    end if

  end subroutine ad_sweeper_f_comp

  ! RK stepper functions

  ! Evaluate the explicit function at y, t.
  subroutine ad_stepper_f_eval(this, y, t, level_index, f, piece)
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece    

    real(pfdp),          pointer       :: yvec(:), fvec(:)
    real(pfdp)                         :: val

    yvec => array1(y)
    fvec => array1(f)

    if (piece == 1) then
       fvec = a*yvec
    else 
       fvec = d*yvec
    end if

  end subroutine ad_stepper_f_eval

  ! Solve for y and return f2 also.
  subroutine ad_stepper_f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    
    real(pfdp),          pointer       :: yvec(:), rhsvec(:), fvec(:)
    real(pfdp)                         :: val

    yvec   => array1(y)
    rhsvec => array1(rhs)
    fvec   => array1(f)
    
    if (piece == 2) then
       yvec = rhsvec / (1.0_pfdp - d*dtq)
       fvec = (yvec - rhsvec) / dtq
    else 
       stop("Error: Piece 1 is explicit")
    end if

  end subroutine ad_stepper_f_comp

  subroutine ad_stepper_f_finalize(this, y, t, dtq, level_index)
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    integer,             intent(in   ) :: level_index

    real(pfdp),          pointer       :: yvec(:)
    
    yvec   => array1(y)
    
    yvec = yvec / (1 - r * dtq)
    
  end subroutine ad_stepper_f_finalize


  ! User_level functions

  subroutine interpolate(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    
    ! just copy the solution to the fine level
    call qF%copy(qG)

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t

    ! just copy the solution to the coarse level
    call qG%copy(qF)

  end subroutine restrict

end module feval

